# See COPYRIGHT file at the top of the source tree.
#
# This file is part of fgcmcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Make the final fgcmcal output products.

This task takes the final output from fgcmFitCycle and produces the following
outputs for use in the DM stack: the FGCM standard stars in a reference
catalog format; the model atmospheres in "transmission_atmosphere_fgcm"
format; and the zeropoints in "fgcm_photoCalib" format.  Optionally, the
task can transfer the 'absolute' calibration from a reference catalog
to put the fgcm standard stars in units of Jansky.  This is accomplished
by matching stars in a sample of healpix pixels, and applying the median
offset per band.
"""
import sys
import traceback
import copy

import numpy as np
import healpy as hp
import esutil
from astropy import units

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
from lsst.afw.image import TransmissionCurve
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.meas.algorithms import ReferenceObjectLoader
from lsst.pipe.tasks.photoCal import PhotoCalTask
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.meas.algorithms import IndexerRegistry
from lsst.meas.algorithms import DatasetConfig
from lsst.meas.algorithms.ingestIndexReferenceTask import addRefCatMetadata

from .utilities import computeApproxPixelAreaFields
from .utilities import lookupStaticCalibrations

import fgcm

__all__ = ['FgcmOutputProductsConfig', 'FgcmOutputProductsTask', 'FgcmOutputProductsRunner']


class FgcmOutputProductsConnections(pipeBase.PipelineTaskConnections,
                                    dimensions=("instrument",),
                                    defaultTemplates={"cycleNumber": "0"}):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibrations,
        isCalibration=True,
    )

    fgcmLookUpTable = connectionTypes.PrerequisiteInput(
        doc=("Atmosphere + instrument look-up-table for FGCM throughput and "
             "chromatic corrections."),
        name="fgcmLookUpTable",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmVisitCatalog = connectionTypes.PrerequisiteInput(
        doc="Catalog of visit information for fgcm",
        name="fgcmVisitCatalog",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStandardStars = connectionTypes.PrerequisiteInput(
        doc="Catalog of standard star data from fgcm fit",
        name="fgcmStandardStars{cycleNumber}",
        storageClass="SimpleCatalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmZeropoints = connectionTypes.PrerequisiteInput(
        doc="Catalog of zeropoints from fgcm fit",
        name="fgcmZeropoints{cycleNumber}",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmAtmosphereParameters = connectionTypes.PrerequisiteInput(
        doc="Catalog of atmosphere parameters from fgcm fit",
        name="fgcmAtmosphereParameters{cycleNumber}",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    refCat = connectionTypes.PrerequisiteInput(
        doc="Reference catalog to use for photometric calibration",
        name="cal_ref_cat",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    fgcmBuildStarsTableConfig = connectionTypes.PrerequisiteInput(
        doc="Config used to build FGCM input stars",
        name="fgcmBuildStarsTable_config",
        storageClass="Config",
    )

    fgcmPhotoCalib = connectionTypes.Output(
        doc="Per-visit photoCalib exposure catalogs produced from fgcm calibration",
        name="fgcmPhotoCalibCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit",),
        multiple=True,
    )

    fgcmTransmissionAtmosphere = connectionTypes.Output(
        doc="Per-visit atmosphere transmission files produced from fgcm calibration",
        name="transmission_atmosphere_fgcm",
        storageClass="TransmissionCurve",
        dimensions=("instrument",
                    "visit",),
        multiple=True,
    )

    fgcmOffsets = connectionTypes.Output(
        doc="Per-band offsets computed from doReferenceCalibration",
        name="fgcmReferenceCalibrationOffsets",
        storageClass="Catalog",
        dimensions=("instrument",),
        multiple=False,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if str(int(config.connections.cycleNumber)) != config.connections.cycleNumber:
            raise ValueError("cycleNumber must be of integer format")
        if config.connections.refCat != config.refObjLoader.ref_dataset_name:
            raise ValueError("connections.refCat must be the same as refObjLoader.ref_dataset_name")

        if config.doRefcatOutput:
            raise ValueError("FgcmOutputProductsTask (Gen3) does not support doRefcatOutput")

        if not config.doReferenceCalibration:
            self.prerequisiteInputs.remove("refCat")
        if not config.doAtmosphereOutput:
            self.prerequisiteInputs.remove("fgcmAtmosphereParameters")
        if not config.doZeropointOutput:
            self.prerequisiteInputs.remove("fgcmZeropoints")
        if not config.doReferenceCalibration:
            self.outputs.remove("fgcmOffsets")


class FgcmOutputProductsConfig(pipeBase.PipelineTaskConfig,
                               pipelineConnections=FgcmOutputProductsConnections):
    """Config for FgcmOutputProductsTask"""

    cycleNumber = pexConfig.Field(
        doc="Final fit cycle from FGCM fit",
        dtype=int,
        default=None,
    )

    # The following fields refer to calibrating from a reference
    # catalog, but in the future this might need to be expanded
    doReferenceCalibration = pexConfig.Field(
        doc=("Transfer 'absolute' calibration from reference catalog? "
             "This afterburner step is unnecessary if reference stars "
             "were used in the full fit in FgcmFitCycleTask."),
        dtype=bool,
        default=False,
    )
    doRefcatOutput = pexConfig.Field(
        doc="Output standard stars in reference catalog format",
        dtype=bool,
        default=True,
    )
    doAtmosphereOutput = pexConfig.Field(
        doc="Output atmospheres in transmission_atmosphere_fgcm format",
        dtype=bool,
        default=True,
    )
    doZeropointOutput = pexConfig.Field(
        doc="Output zeropoints in fgcm_photoCalib format",
        dtype=bool,
        default=True,
    )
    doComposeWcsJacobian = pexConfig.Field(
        doc="Compose Jacobian of WCS with fgcm calibration for output photoCalib?",
        dtype=bool,
        default=True,
    )
    doApplyMeanChromaticCorrection = pexConfig.Field(
        doc="Apply the mean chromatic correction to the zeropoints?",
        dtype=bool,
        default=True,
    )
    refObjLoader = pexConfig.ConfigurableField(
        target=LoadIndexedReferenceObjectsTask,
        doc="reference object loader for 'absolute' photometric calibration",
    )
    photoCal = pexConfig.ConfigurableField(
        target=PhotoCalTask,
        doc="task to perform 'absolute' calibration",
    )
    referencePixelizationNside = pexConfig.Field(
        doc="Healpix nside to pixelize catalog to compare to reference catalog",
        dtype=int,
        default=64,
    )
    referencePixelizationMinStars = pexConfig.Field(
        doc=("Minimum number of stars per healpix pixel to select for comparison"
             "to the specified reference catalog"),
        dtype=int,
        default=200,
    )
    referenceMinMatch = pexConfig.Field(
        doc="Minimum number of stars matched to reference catalog to be used in statistics",
        dtype=int,
        default=50,
    )
    referencePixelizationNPixels = pexConfig.Field(
        doc=("Number of healpix pixels to sample to do comparison. "
             "Doing too many will take a long time and not yield any more "
             "precise results because the final number is the median offset "
             "(per band) from the set of pixels."),
        dtype=int,
        default=100,
    )
    datasetConfig = pexConfig.ConfigField(
        dtype=DatasetConfig,
        doc="Configuration for writing/reading ingested catalog",
    )

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

        # In order to transfer the "absolute" calibration from a reference
        # catalog to the relatively calibrated FGCM standard stars (one number
        # per band), we use the PhotoCalTask to match stars in a sample of healpix
        # pixels.  These basic settings ensure that only well-measured, good stars
        # from the source and reference catalogs are used for the matching.

        # applyColorTerms needs to be False if doReferenceCalibration is False,
        # as is the new default after DM-16702
        self.photoCal.applyColorTerms = False
        self.photoCal.fluxField = 'instFlux'
        self.photoCal.magErrFloor = 0.003
        self.photoCal.match.referenceSelection.doSignalToNoise = True
        self.photoCal.match.referenceSelection.signalToNoise.minimum = 10.0
        self.photoCal.match.sourceSelection.doSignalToNoise = True
        self.photoCal.match.sourceSelection.signalToNoise.minimum = 10.0
        self.photoCal.match.sourceSelection.signalToNoise.fluxField = 'instFlux'
        self.photoCal.match.sourceSelection.signalToNoise.errField = 'instFluxErr'
        self.photoCal.match.sourceSelection.doFlags = True
        self.photoCal.match.sourceSelection.flags.good = []
        self.photoCal.match.sourceSelection.flags.bad = ['flag_badStar']
        self.photoCal.match.sourceSelection.doUnresolved = False
        self.datasetConfig.ref_dataset_name = 'fgcm_stars'
        self.datasetConfig.format_version = 1

    def validate(self):
        super().validate()

        # Force the connections to conform with cycleNumber
        self.connections.cycleNumber = str(self.cycleNumber)


class FgcmOutputProductsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmOutputProductsTask

    fgcmOutputProductsTask.run() takes one argument, the butler, and
    does not run on any data in the repository.
    This runner does not use any parallelization.
    """

    @staticmethod
    def getTargetList(parsedCmd):
        """
        Return a list with one element, the butler.
        """
        return [parsedCmd.butler]

    def __call__(self, butler):
        """
        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`

        Returns
        -------
        exitStatus: `list` with `pipeBase.Struct`
           exitStatus (0: success; 1: failure)
           if self.doReturnResults also
           results (`np.array` with absolute zeropoint offsets)
        """
        task = self.TaskClass(butler=butler, config=self.config, log=self.log)

        exitStatus = 0
        if self.doRaise:
            results = task.runDataRef(butler)
        else:
            try:
                results = task.runDataRef(butler)
            except Exception as e:
                exitStatus = 1
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)

        if self.doReturnResults:
            # The results here are the zeropoint offsets for each band
            return [pipeBase.Struct(exitStatus=exitStatus,
                                    results=results)]
        else:
            return [pipeBase.Struct(exitStatus=exitStatus)]

    def run(self, parsedCmd):
        """
        Run the task, with no multiprocessing

        Parameters
        ----------
        parsedCmd: `lsst.pipe.base.ArgumentParser` parsed command line
        """

        resultList = []

        if self.precall(parsedCmd):
            targetList = self.getTargetList(parsedCmd)
            # make sure that we only get 1
            resultList = self(targetList[0])

        return resultList


class FgcmOutputProductsTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    """
    Output products from FGCM global calibration.
    """

    ConfigClass = FgcmOutputProductsConfig
    RunnerClass = FgcmOutputProductsRunner
    _DefaultName = "fgcmOutputProducts"

    def __init__(self, butler=None, **kwargs):
        super().__init__(**kwargs)

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataRefDict = {}
        dataRefDict['camera'] = butlerQC.get(inputRefs.camera)
        dataRefDict['fgcmLookUpTable'] = butlerQC.get(inputRefs.fgcmLookUpTable)
        dataRefDict['fgcmVisitCatalog'] = butlerQC.get(inputRefs.fgcmVisitCatalog)
        dataRefDict['fgcmStandardStars'] = butlerQC.get(inputRefs.fgcmStandardStars)

        if self.config.doZeropointOutput:
            dataRefDict['fgcmZeropoints'] = butlerQC.get(inputRefs.fgcmZeropoints)
            photoCalibRefDict = {photoCalibRef.dataId.byName()['visit']:
                                 photoCalibRef for photoCalibRef in outputRefs.fgcmPhotoCalib}

        if self.config.doAtmosphereOutput:
            dataRefDict['fgcmAtmosphereParameters'] = butlerQC.get(inputRefs.fgcmAtmosphereParameters)
            atmRefDict = {atmRef.dataId.byName()['visit']: atmRef for
                          atmRef in outputRefs.fgcmTransmissionAtmosphere}

        if self.config.doReferenceCalibration:
            refConfig = self.config.refObjLoader
            self.refObjLoader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                               for ref in inputRefs.refCat],
                                                      refCats=butlerQC.get(inputRefs.refCat),
                                                      config=refConfig,
                                                      log=self.log)
        else:
            self.refObjLoader = None

        dataRefDict['fgcmBuildStarsTableConfig'] = butlerQC.get(inputRefs.fgcmBuildStarsTableConfig)

        fgcmBuildStarsConfig = butlerQC.get(inputRefs.fgcmBuildStarsTableConfig)
        filterMap = fgcmBuildStarsConfig.filterMap

        if self.config.doComposeWcsJacobian and not fgcmBuildStarsConfig.doApplyWcsJacobian:
            raise RuntimeError("Cannot compose the WCS jacobian if it hasn't been applied "
                               "in fgcmBuildStarsTask.")
        if not self.config.doComposeWcsJacobian and fgcmBuildStarsConfig.doApplyWcsJacobian:
            self.log.warn("Jacobian was applied in build-stars but doComposeWcsJacobian is not set.")

        struct = self.run(dataRefDict, filterMap, returnCatalogs=True)

        # Output the photoCalib exposure catalogs
        if struct.photoCalibCatalogs is not None:
            self.log.info("Outputting photoCalib catalogs.")
            for visit, expCatalog in struct.photoCalibCatalogs:
                butlerQC.put(expCatalog, photoCalibRefDict[visit])
            self.log.info("Done outputting photoCalib catalogs.")

        # Output the atmospheres
        if struct.atmospheres is not None:
            self.log.info("Outputting atmosphere transmission files.")
            for visit, atm in struct.atmospheres:
                butlerQC.put(atm, atmRefDict[visit])
            self.log.info("Done outputting atmosphere files.")

        if self.config.doReferenceCalibration:
            # Turn offset into simple catalog for persistence if necessary
            schema = afwTable.Schema()
            schema.addField('offset', type=np.float64,
                            doc="Post-process calibration offset (mag)")
            offsetCat = afwTable.BaseCatalog(schema)
            offsetCat.resize(len(struct.offsets))
            offsetCat['offset'][:] = struct.offsets

            butlerQC.put(offsetCat, outputRefs.fgcmOffsets)

        return

    @pipeBase.timeMethod
    def runDataRef(self, butler):
        """
        Make FGCM output products for use in the stack

        Parameters
        ----------
        butler:  `lsst.daf.persistence.Butler`
        cycleNumber: `int`
           Final fit cycle number, override config.

        Returns
        -------
        offsets: `lsst.pipe.base.Struct`
           A structure with array of zeropoint offsets

        Raises
        ------
        RuntimeError:
           Raised if any one of the following is true:

           - butler cannot find "fgcmBuildStars_config" or
             "fgcmBuildStarsTable_config".
           - butler cannot find "fgcmFitCycle_config".
           - "fgcmFitCycle_config" does not refer to
             `self.config.cycleNumber`.
           - butler cannot find "fgcmAtmosphereParameters" and
             `self.config.doAtmosphereOutput` is `True`.
           - butler cannot find "fgcmStandardStars" and
             `self.config.doReferenceCalibration` is `True` or
             `self.config.doRefcatOutput` is `True`.
           - butler cannot find "fgcmZeropoints" and
             `self.config.doZeropointOutput` is `True`.
        """
        if self.config.doReferenceCalibration:
            # We need the ref obj loader to get the flux field
            self.makeSubtask("refObjLoader", butler=butler)

        # Check to make sure that the fgcmBuildStars config exists, to retrieve
        # the visit and ccd dataset tags
        if not butler.datasetExists('fgcmBuildStarsTable_config') and \
                not butler.datasetExists('fgcmBuildStars_config'):
            raise RuntimeError("Cannot find fgcmBuildStarsTable_config or fgcmBuildStars_config, "
                               "which is prereq for fgcmOutputProducts")

        if butler.datasetExists('fgcmBuildStarsTable_config'):
            fgcmBuildStarsConfig = butler.get('fgcmBuildStarsTable_config')
        else:
            fgcmBuildStarsConfig = butler.get('fgcmBuildStars_config')
        visitDataRefName = fgcmBuildStarsConfig.visitDataRefName
        ccdDataRefName = fgcmBuildStarsConfig.ccdDataRefName
        filterMap = fgcmBuildStarsConfig.filterMap

        if self.config.doComposeWcsJacobian and not fgcmBuildStarsConfig.doApplyWcsJacobian:
            raise RuntimeError("Cannot compose the WCS jacobian if it hasn't been applied "
                               "in fgcmBuildStarsTask.")

        if not self.config.doComposeWcsJacobian and fgcmBuildStarsConfig.doApplyWcsJacobian:
            self.log.warn("Jacobian was applied in build-stars but doComposeWcsJacobian is not set.")

        # And make sure that the atmosphere was output properly
        if (self.config.doAtmosphereOutput
                and not butler.datasetExists('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)):
            raise RuntimeError(f"Atmosphere parameters are missing for cycle {self.config.cycleNumber}.")

        if not butler.datasetExists('fgcmStandardStars',
                                    fgcmcycle=self.config.cycleNumber):
            raise RuntimeError("Standard stars are missing for cycle %d." %
                               (self.config.cycleNumber))

        if (self.config.doZeropointOutput
                and (not butler.datasetExists('fgcmZeropoints', fgcmcycle=self.config.cycleNumber))):
            raise RuntimeError("Zeropoints are missing for cycle %d." %
                               (self.config.cycleNumber))

        dataRefDict = {}
        # This is the _actual_ camera
        dataRefDict['camera'] = butler.get('camera')
        dataRefDict['fgcmLookUpTable'] = butler.dataRef('fgcmLookUpTable')
        dataRefDict['fgcmVisitCatalog'] = butler.dataRef('fgcmVisitCatalog')
        dataRefDict['fgcmStandardStars'] = butler.dataRef('fgcmStandardStars',
                                                          fgcmcycle=self.config.cycleNumber)

        if self.config.doZeropointOutput:
            dataRefDict['fgcmZeropoints'] = butler.dataRef('fgcmZeropoints',
                                                           fgcmcycle=self.config.cycleNumber)
        if self.config.doAtmosphereOutput:
            dataRefDict['fgcmAtmosphereParameters'] = butler.dataRef('fgcmAtmosphereParameters',
                                                                     fgcmcycle=self.config.cycleNumber)

        struct = self.run(dataRefDict, filterMap, butler=butler, returnCatalogs=False)

        if struct.photoCalibs is not None:
            self.log.info("Outputting photoCalib files.")

            filterMapping = {}
            for visit, detector, filtername, photoCalib in struct.photoCalibs:
                if filtername not in filterMapping:
                    # We need to find the mapping from encoded filter to dataid filter,
                    # and this trick allows us to do that.
                    dataId = {visitDataRefName: visit,
                              ccdDataRefName: detector}
                    dataRef = butler.dataRef('raw', dataId=dataId)
                    filterMapping[filtername] = dataRef.dataId['filter']

                butler.put(photoCalib, 'fgcm_photoCalib',
                           dataId={visitDataRefName: visit,
                                   ccdDataRefName: detector,
                                   'filter': filterMapping[filtername]})

            self.log.info("Done outputting photoCalib files.")

        if struct.atmospheres is not None:
            self.log.info("Outputting atmosphere transmission files.")
            for visit, atm in struct.atmospheres:
                butler.put(atm, "transmission_atmosphere_fgcm",
                           dataId={visitDataRefName: visit})
            self.log.info("Done outputting atmosphere transmissions.")

        return pipeBase.Struct(offsets=struct.offsets)

    def run(self, dataRefDict, filterMap, returnCatalogs=True, butler=None):
        """Run the output products task.

        Parameters
        ----------
        dataRefDict : `dict`
            All dataRefs are `lsst.daf.persistence.ButlerDataRef` (gen2) or
            `lsst.daf.butler.DeferredDatasetHandle` (gen3)
            dataRef dictionary with keys:

            ``"camera"``
                Camera object (`lsst.afw.cameraGeom.Camera`)
            ``"fgcmLookUpTable"``
                dataRef for the FGCM look-up table.
            ``"fgcmVisitCatalog"``
                dataRef for visit summary catalog.
            ``"fgcmStandardStars"``
                dataRef for the output standard star catalog.
            ``"fgcmZeropoints"``
                dataRef for the zeropoint data catalog.
            ``"fgcmAtmosphereParameters"``
                dataRef for the atmosphere parameter catalog.
            ``"fgcmBuildStarsTableConfig"``
                Config for `lsst.fgcmcal.fgcmBuildStarsTableTask`.
        filterMap : `dict`
            Dictionary of mappings from filter to FGCM band.
        returnCatalogs : `bool`, optional
            Return photoCalibs as per-visit exposure catalogs.
        butler : `lsst.daf.persistence.Butler`, optional
            Gen2 butler used for reference star outputs

        Returns
        -------
        retStruct : `lsst.pipe.base.Struct`
            Output structure with keys:

            offsets : `np.ndarray`
                Final reference offsets, per band.
            atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
                Generator that returns (visit, transmissionCurve) tuples.
            photoCalibs : `generator` [(`int`, `int`, `str`, `lsst.afw.image.PhotoCalib`)]
                Generator that returns (visit, ccd, filtername, photoCalib) tuples.
                (returned if returnCatalogs is False).
            photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
                Generator that returns (visit, exposureCatalog) tuples.
                (returned if returnCatalogs is True).
        """
        stdCat = dataRefDict['fgcmStandardStars'].get()
        md = stdCat.getMetadata()
        bands = md.getArray('BANDS')

        if self.config.doReferenceCalibration:
            offsets = self._computeReferenceOffsets(stdCat, bands)
        else:
            offsets = np.zeros(len(bands))

        # This is Gen2 only, and requires the butler.
        if self.config.doRefcatOutput and butler is not None:
            self._outputStandardStars(butler, stdCat, offsets, bands, self.config.datasetConfig)

        del stdCat

        if self.config.doZeropointOutput:
            zptCat = dataRefDict['fgcmZeropoints'].get()
            visitCat = dataRefDict['fgcmVisitCatalog'].get()

            pcgen = self._outputZeropoints(dataRefDict['camera'], zptCat, visitCat, offsets, bands,
                                           filterMap, returnCatalogs=returnCatalogs)
        else:
            pcgen = None

        if self.config.doAtmosphereOutput:
            atmCat = dataRefDict['fgcmAtmosphereParameters'].get()
            atmgen = self._outputAtmospheres(dataRefDict, atmCat)
        else:
            atmgen = None

        retStruct = pipeBase.Struct(offsets=offsets,
                                    atmospheres=atmgen)
        if returnCatalogs:
            retStruct.photoCalibCatalogs = pcgen
        else:
            retStruct.photoCalibs = pcgen

        return retStruct

    def generateTractOutputProducts(self, dataRefDict, tract,
                                    visitCat, zptCat, atmCat, stdCat,
                                    fgcmBuildStarsConfig,
                                    returnCatalogs=True,
                                    butler=None):
        """
        Generate the output products for a given tract, as specified in the config.

        This method is here to have an alternate entry-point for
        FgcmCalibrateTract.

        Parameters
        ----------
        dataRefDict : `dict`
            All dataRefs are `lsst.daf.persistence.ButlerDataRef` (gen2) or
            `lsst.daf.butler.DeferredDatasetHandle` (gen3)
            dataRef dictionary with keys:

            ``"camera"``
                Camera object (`lsst.afw.cameraGeom.Camera`)
            ``"fgcmLookUpTable"``
                dataRef for the FGCM look-up table.
        tract : `int`
            Tract number
        visitCat : `lsst.afw.table.BaseCatalog`
            FGCM visitCat from `FgcmBuildStarsTask`
        zptCat : `lsst.afw.table.BaseCatalog`
            FGCM zeropoint catalog from `FgcmFitCycleTask`
        atmCat : `lsst.afw.table.BaseCatalog`
            FGCM atmosphere parameter catalog from `FgcmFitCycleTask`
        stdCat : `lsst.afw.table.SimpleCatalog`
            FGCM standard star catalog from `FgcmFitCycleTask`
        fgcmBuildStarsConfig : `lsst.fgcmcal.FgcmBuildStarsConfig`
            Configuration object from `FgcmBuildStarsTask`
        returnCatalogs : `bool`, optional
            Return photoCalibs as per-visit exposure catalogs.
        butler: `lsst.daf.persistence.Butler`, optional
            Gen2 butler used for reference star outputs

        Returns
        -------
        retStruct : `lsst.pipe.base.Struct`
            Output structure with keys:

            offsets : `np.ndarray`
                Final reference offsets, per band.
            atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
                Generator that returns (visit, transmissionCurve) tuples.
            photoCalibs : `generator` [(`int`, `int`, `str`, `lsst.afw.image.PhotoCalib`)]
                Generator that returns (visit, ccd, filtername, photoCalib) tuples.
                (returned if returnCatalogs is False).
            photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
                Generator that returns (visit, exposureCatalog) tuples.
                (returned if returnCatalogs is True).
        """
        filterMap = fgcmBuildStarsConfig.filterMap

        md = stdCat.getMetadata()
        bands = md.getArray('BANDS')

        if self.config.doComposeWcsJacobian and not fgcmBuildStarsConfig.doApplyWcsJacobian:
            raise RuntimeError("Cannot compose the WCS jacobian if it hasn't been applied "
                               "in fgcmBuildStarsTask.")

        if not self.config.doComposeWcsJacobian and fgcmBuildStarsConfig.doApplyWcsJacobian:
            self.log.warn("Jacobian was applied in build-stars but doComposeWcsJacobian is not set.")

        if self.config.doReferenceCalibration:
            offsets = self._computeReferenceOffsets(stdCat, bands)
        else:
            offsets = np.zeros(len(bands))

        if self.config.doRefcatOutput and butler is not None:
            # Create a special config that has the tract number in it
            datasetConfig = copy.copy(self.config.datasetConfig)
            datasetConfig.ref_dataset_name = '%s_%d' % (self.config.datasetConfig.ref_dataset_name,
                                                        tract)
            self._outputStandardStars(butler, stdCat, offsets, bands, datasetConfig)

        if self.config.doZeropointOutput:
            pcgen = self._outputZeropoints(dataRefDict['camera'], zptCat, visitCat, offsets, bands,
                                           filterMap, returnCatalogs=returnCatalogs)
        else:
            pcgen = None

        if self.config.doAtmosphereOutput:
            atmgen = self._outputAtmospheres(dataRefDict, atmCat)
        else:
            atmgen = None

        retStruct = pipeBase.Struct(offsets=offsets,
                                    atmospheres=atmgen)
        if returnCatalogs:
            retStruct.photoCalibCatalogs = pcgen
        else:
            retStruct.photoCalibs = pcgen

        return retStruct

    def _computeReferenceOffsets(self, stdCat, bands):
        """
        Compute offsets relative to a reference catalog.

        This method splits the star catalog into healpix pixels
        and computes the calibration transfer for a sample of
        these pixels to approximate the 'absolute' calibration
        values (on for each band) to apply to transfer the
        absolute scale.

        Parameters
        ----------
        stdCat : `lsst.afw.table.SimpleCatalog`
            FGCM standard stars
        bands : `list` [`str`]
            List of band names from FGCM output
        Returns
        -------
        offsets : `numpy.array` of floats
            Per band zeropoint offsets
        """

        # Only use stars that are observed in all the bands that were actually used
        # This will ensure that we use the same healpix pixels for the absolute
        # calibration of each band.
        minObs = stdCat['ngood'].min(axis=1)

        goodStars = (minObs >= 1)
        stdCat = stdCat[goodStars]

        self.log.info("Found %d stars with at least 1 good observation in each band" %
                      (len(stdCat)))

        # We have to make a table for each pixel with flux/fluxErr
        # This is a temporary table generated for input to the photoCal task.
        # These fluxes are not instFlux (they are top-of-the-atmosphere approximate and
        # have had chromatic corrections applied to get to the standard system
        # specified by the atmosphere/instrumental parameters), nor are they
        # in Jansky (since they don't have a proper absolute calibration: the overall
        # zeropoint is estimated from the telescope size, etc.)
        sourceMapper = afwTable.SchemaMapper(stdCat.schema)
        sourceMapper.addMinimalSchema(afwTable.SimpleTable.makeMinimalSchema())
        sourceMapper.editOutputSchema().addField('instFlux', type=np.float64,
                                                 doc="instrumental flux (counts)")
        sourceMapper.editOutputSchema().addField('instFluxErr', type=np.float64,
                                                 doc="instrumental flux error (counts)")
        badStarKey = sourceMapper.editOutputSchema().addField('flag_badStar',
                                                              type='Flag',
                                                              doc="bad flag")

        # Split up the stars
        # Note that there is an assumption here that the ra/dec coords stored
        # on-disk are in radians, and therefore that starObs['coord_ra'] /
        # starObs['coord_dec'] return radians when used as an array of numpy float64s.
        theta = np.pi/2. - stdCat['coord_dec']
        phi = stdCat['coord_ra']

        ipring = hp.ang2pix(self.config.referencePixelizationNside, theta, phi)
        h, rev = esutil.stat.histogram(ipring, rev=True)

        gdpix, = np.where(h >= self.config.referencePixelizationMinStars)

        self.log.info("Found %d pixels (nside=%d) with at least %d good stars" %
                      (gdpix.size,
                       self.config.referencePixelizationNside,
                       self.config.referencePixelizationMinStars))

        if gdpix.size < self.config.referencePixelizationNPixels:
            self.log.warn("Found fewer good pixels (%d) than preferred in configuration (%d)" %
                          (gdpix.size, self.config.referencePixelizationNPixels))
        else:
            # Sample out the pixels we want to use
            gdpix = np.random.choice(gdpix, size=self.config.referencePixelizationNPixels, replace=False)

        results = np.zeros(gdpix.size, dtype=[('hpix', 'i4'),
                                              ('nstar', 'i4', len(bands)),
                                              ('nmatch', 'i4', len(bands)),
                                              ('zp', 'f4', len(bands)),
                                              ('zpErr', 'f4', len(bands))])
        results['hpix'] = ipring[rev[rev[gdpix]]]

        # We need a boolean index to deal with catalogs...
        selected = np.zeros(len(stdCat), dtype=np.bool)

        refFluxFields = [None]*len(bands)

        for p, pix in enumerate(gdpix):
            i1a = rev[rev[pix]: rev[pix + 1]]

            # the stdCat afwTable can only be indexed with boolean arrays,
            # and not numpy index arrays (see DM-16497).  This little trick
            # converts the index array into a boolean array
            selected[:] = False
            selected[i1a] = True

            for b, band in enumerate(bands):

                struct = self._computeOffsetOneBand(sourceMapper, badStarKey, b, band, stdCat,
                                                    selected, refFluxFields)
                results['nstar'][p, b] = len(i1a)
                results['nmatch'][p, b] = len(struct.arrays.refMag)
                results['zp'][p, b] = struct.zp
                results['zpErr'][p, b] = struct.sigma

        # And compute the summary statistics
        offsets = np.zeros(len(bands))

        for b, band in enumerate(bands):
            # make configurable
            ok, = np.where(results['nmatch'][:, b] >= self.config.referenceMinMatch)
            offsets[b] = np.median(results['zp'][ok, b])
            # use median absolute deviation to estimate Normal sigma
            # see https://en.wikipedia.org/wiki/Median_absolute_deviation
            madSigma = 1.4826*np.median(np.abs(results['zp'][ok, b] - offsets[b]))
            self.log.info("Reference catalog offset for %s band: %.12f +/- %.12f" %
                          (band, offsets[b], madSigma))

        return offsets

    def _computeOffsetOneBand(self, sourceMapper, badStarKey,
                              b, band, stdCat, selected, refFluxFields):
        """
        Compute the zeropoint offset between the fgcm stdCat and the reference
        stars for one pixel in one band

        Parameters
        ----------
        sourceMapper: `lsst.afw.table.SchemaMapper`
           Mapper to go from stdCat to calibratable catalog
        badStarKey: `lsst.afw.table.Key`
           Key for the field with bad stars
        b: `int`
           Index of the band in the star catalog
        band: `str`
           Name of band for reference catalog
        stdCat: `lsst.afw.table.SimpleCatalog`
           FGCM standard stars
        selected: `numpy.array(dtype=np.bool)`
           Boolean array of which stars are in the pixel
        refFluxFields: `list`
           List of names of flux fields for reference catalog
        """

        sourceCat = afwTable.SimpleCatalog(sourceMapper.getOutputSchema())
        sourceCat.reserve(selected.sum())
        sourceCat.extend(stdCat[selected], mapper=sourceMapper)
        sourceCat['instFlux'] = 10.**(stdCat['mag_std_noabs'][selected, b]/(-2.5))
        sourceCat['instFluxErr'] = (np.log(10.)/2.5)*(stdCat['magErr_std'][selected, b]
                                                      * sourceCat['instFlux'])
        # Make sure we only use stars that have valid measurements
        # (This is perhaps redundant with requirements above that the
        # stars be observed in all bands, but it can't hurt)
        badStar = (stdCat['mag_std_noabs'][selected, b] > 90.0)
        for rec in sourceCat[badStar]:
            rec.set(badStarKey, True)

        exposure = afwImage.ExposureF()
        exposure.setFilter(afwImage.Filter(band))

        if refFluxFields[b] is None:
            # Need to find the flux field in the reference catalog
            # to work around limitations of DirectMatch in PhotoCal
            ctr = stdCat[0].getCoord()
            rad = 0.05*lsst.geom.degrees
            refDataTest = self.refObjLoader.loadSkyCircle(ctr, rad, band)
            refFluxFields[b] = refDataTest.fluxField

        # Make a copy of the config so that we can modify it
        calConfig = copy.copy(self.config.photoCal.value)
        calConfig.match.referenceSelection.signalToNoise.fluxField = refFluxFields[b]
        calConfig.match.referenceSelection.signalToNoise.errField = refFluxFields[b] + 'Err'
        calTask = self.config.photoCal.target(refObjLoader=self.refObjLoader,
                                              config=calConfig,
                                              schema=sourceCat.getSchema())

        struct = calTask.run(exposure, sourceCat)

        return struct

    def _outputStandardStars(self, butler, stdCat, offsets, bands, datasetConfig):
        """
        Output standard stars in indexed reference catalog format.
        This is not currently supported in Gen3.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
        stdCat : `lsst.afw.table.SimpleCatalog`
            FGCM standard star catalog from fgcmFitCycleTask
        offsets : `numpy.array` of floats
            Per band zeropoint offsets
        bands : `list` [`str`]
            List of band names from FGCM output
        datasetConfig : `lsst.meas.algorithms.DatasetConfig`
            Config for reference dataset
        """

        self.log.info("Outputting standard stars to %s" % (datasetConfig.ref_dataset_name))

        indexer = IndexerRegistry[self.config.datasetConfig.indexer.name](
            self.config.datasetConfig.indexer.active)

        # We determine the conversion from the native units (typically radians) to
        # degrees for the first star.  This allows us to treat coord_ra/coord_dec as
        # numpy arrays rather than Angles, which would we approximately 600x slower.
        # TODO: Fix this after DM-16524 (HtmIndexer.indexPoints should take coords
        # (as Angles) for input
        conv = stdCat[0]['coord_ra'].asDegrees()/float(stdCat[0]['coord_ra'])
        indices = np.array(indexer.indexPoints(stdCat['coord_ra']*conv,
                                               stdCat['coord_dec']*conv))

        formattedCat = self._formatCatalog(stdCat, offsets, bands)

        # Write the master schema
        dataId = indexer.makeDataId('master_schema',
                                    datasetConfig.ref_dataset_name)
        masterCat = afwTable.SimpleCatalog(formattedCat.schema)
        addRefCatMetadata(masterCat)
        butler.put(masterCat, 'ref_cat', dataId=dataId)

        # Break up the pixels using a histogram
        h, rev = esutil.stat.histogram(indices, rev=True)
        gd, = np.where(h > 0)
        selected = np.zeros(len(formattedCat), dtype=np.bool)
        for i in gd:
            i1a = rev[rev[i]: rev[i + 1]]

            # the formattedCat afwTable can only be indexed with boolean arrays,
            # and not numpy index arrays (see DM-16497).  This little trick
            # converts the index array into a boolean array
            selected[:] = False
            selected[i1a] = True

            # Write the individual pixel
            dataId = indexer.makeDataId(indices[i1a[0]],
                                        datasetConfig.ref_dataset_name)
            butler.put(formattedCat[selected], 'ref_cat', dataId=dataId)

        # And save the dataset configuration
        dataId = indexer.makeDataId(None, datasetConfig.ref_dataset_name)
        butler.put(datasetConfig, 'ref_cat_config', dataId=dataId)

        self.log.info("Done outputting standard stars.")

    def _formatCatalog(self, fgcmStarCat, offsets, bands):
        """
        Turn an FGCM-formatted star catalog, applying zeropoint offsets.

        Parameters
        ----------
        fgcmStarCat : `lsst.afw.Table.SimpleCatalog`
            SimpleCatalog as output by fgcmcal
        offsets : `list` with len(self.bands) entries
            Zeropoint offsets to apply
        bands : `list` [`str`]
            List of band names from FGCM output

        Returns
        -------
        formattedCat: `lsst.afw.table.SimpleCatalog`
           SimpleCatalog suitable for using as a reference catalog
        """

        sourceMapper = afwTable.SchemaMapper(fgcmStarCat.schema)
        minSchema = LoadIndexedReferenceObjectsTask.makeMinimalSchema(bands,
                                                                      addCentroid=False,
                                                                      addIsResolved=True,
                                                                      coordErrDim=0)
        sourceMapper.addMinimalSchema(minSchema)
        for band in bands:
            sourceMapper.editOutputSchema().addField('%s_nGood' % (band), type=np.int32)
            sourceMapper.editOutputSchema().addField('%s_nTotal' % (band), type=np.int32)
            sourceMapper.editOutputSchema().addField('%s_nPsfCandidate' % (band), type=np.int32)

        formattedCat = afwTable.SimpleCatalog(sourceMapper.getOutputSchema())
        formattedCat.reserve(len(fgcmStarCat))
        formattedCat.extend(fgcmStarCat, mapper=sourceMapper)

        # Note that we don't have to set `resolved` because the default is False

        for b, band in enumerate(bands):
            mag = fgcmStarCat['mag_std_noabs'][:, b].astype(np.float64) + offsets[b]
            # We want fluxes in nJy from calibrated AB magnitudes
            # (after applying offset).  Updated after RFC-549 and RFC-575.
            flux = (mag*units.ABmag).to_value(units.nJy)
            fluxErr = (np.log(10.)/2.5)*flux*fgcmStarCat['magErr_std'][:, b].astype(np.float64)

            formattedCat['%s_flux' % (band)][:] = flux
            formattedCat['%s_fluxErr' % (band)][:] = fluxErr
            formattedCat['%s_nGood' % (band)][:] = fgcmStarCat['ngood'][:, b]
            formattedCat['%s_nTotal' % (band)][:] = fgcmStarCat['ntotal'][:, b]
            formattedCat['%s_nPsfCandidate' % (band)][:] = fgcmStarCat['npsfcand'][:, b]

        addRefCatMetadata(formattedCat)

        return formattedCat

    def _outputZeropoints(self, camera, zptCat, visitCat, offsets, bands,
                          filterMap, returnCatalogs=True,
                          tract=None):
        """Output the zeropoints in fgcm_photoCalib format.

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`
            Camera from the butler.
        zptCat : `lsst.afw.table.BaseCatalog`
            FGCM zeropoint catalog from `FgcmFitCycleTask`.
        visitCat : `lsst.afw.table.BaseCatalog`
            FGCM visitCat from `FgcmBuildStarsTask`.
        offsets : `numpy.array`
            Float array of absolute calibration offsets, one for each filter.
        bands : `list` [`str`]
            List of band names from FGCM output.
        filterMap : `dict`
            Dictionary of mappings from filter to FGCM band.
        returnCatalogs : `bool`, optional
            Return photoCalibs as per-visit exposure catalogs.
        tract: `int`, optional
            Tract number to output.  Default is None (global calibration)

        Returns
        -------
        photoCalibs : `generator` [(`int`, `int`, `str`, `lsst.afw.image.PhotoCalib`)]
            Generator that returns (visit, ccd, filtername, photoCalib) tuples.
            (returned if returnCatalogs is False).
        photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
            Generator that returns (visit, exposureCatalog) tuples.
            (returned if returnCatalogs is True).
        """
        # Select visit/ccds where we have a calibration
        # This includes ccds where we were able to interpolate from neighboring
        # ccds.
        cannot_compute = fgcm.fgcmUtilities.zpFlagDict['CANNOT_COMPUTE_ZEROPOINT']
        selected = (((zptCat['fgcmFlag'] & cannot_compute) == 0)
                    & (zptCat['fgcmZptVar'] > 0.0))

        # Log warnings for any visit which has no valid zeropoints
        badVisits = np.unique(zptCat['visit'][~selected])
        goodVisits = np.unique(zptCat['visit'][selected])
        allBadVisits = badVisits[~np.isin(badVisits, goodVisits)]
        for allBadVisit in allBadVisits:
            self.log.warn(f'No suitable photoCalib for visit {allBadVisit}')

        # Get a mapping from filtername to the offsets
        offsetMapping = {}
        for f in filterMap:
            # Not every filter in the map will necesarily have a band.
            if filterMap[f] in bands:
                offsetMapping[f] = offsets[bands.index(filterMap[f])]

        # Get a mapping from "ccd" to the ccd index used for the scaling
        ccdMapping = {}
        for ccdIndex, detector in enumerate(camera):
            ccdMapping[detector.getId()] = ccdIndex

        # And a mapping to get the flat-field scaling values
        scalingMapping = {}
        for rec in visitCat:
            scalingMapping[rec['visit']] = rec['scaling']

        if self.config.doComposeWcsJacobian:
            approxPixelAreaFields = computeApproxPixelAreaFields(camera)

        # The zptCat is sorted by visit, which is useful
        lastVisit = -1
        zptCounter = 0
        zptVisitCatalog = None
        for rec in zptCat[selected]:

            # Retrieve overall scaling
            scaling = scalingMapping[rec['visit']][ccdMapping[rec['detector']]]

            # The postCalibrationOffset describe any zeropoint offsets
            # to apply after the fgcm calibration.  The first part comes
            # from the reference catalog match (used in testing).  The
            # second part comes from the mean chromatic correction
            # (if configured).
            postCalibrationOffset = offsetMapping[rec['filtername']]
            if self.config.doApplyMeanChromaticCorrection:
                postCalibrationOffset += rec['fgcmDeltaChrom']

            fgcmSuperStarField = self._getChebyshevBoundedField(rec['fgcmfZptSstarCheb'],
                                                                rec['fgcmfZptChebXyMax'])
            # Convert from FGCM AB to nJy
            fgcmZptField = self._getChebyshevBoundedField((rec['fgcmfZptCheb']*units.AB).to_value(units.nJy),
                                                          rec['fgcmfZptChebXyMax'],
                                                          offset=postCalibrationOffset,
                                                          scaling=scaling)

            if self.config.doComposeWcsJacobian:

                fgcmField = afwMath.ProductBoundedField([approxPixelAreaFields[rec['detector']],
                                                         fgcmSuperStarField,
                                                         fgcmZptField])
            else:
                # The photoCalib is just the product of the fgcmSuperStarField and the
                # fgcmZptField
                fgcmField = afwMath.ProductBoundedField([fgcmSuperStarField, fgcmZptField])

            # The "mean" calibration will be set to the center of the ccd for reference
            calibCenter = fgcmField.evaluate(fgcmField.getBBox().getCenter())
            calibErr = (np.log(10.0)/2.5)*calibCenter*np.sqrt(rec['fgcmZptVar'])
            photoCalib = afwImage.PhotoCalib(calibrationMean=calibCenter,
                                             calibrationErr=calibErr,
                                             calibration=fgcmField,
                                             isConstant=False)

            if not returnCatalogs:
                # Return individual photoCalibs
                yield (int(rec['visit']), int(rec['detector']), rec['filtername'], photoCalib)
            else:
                # Return full per-visit exposure catalogs
                if rec['visit'] != lastVisit:
                    # This is a new visit.  If the last visit was not -1, yield
                    # the ExposureCatalog
                    if lastVisit > -1:
                        yield (int(lastVisit), zptVisitCatalog)
                    else:
                        # We need to create a new schema
                        zptExpCatSchema = afwTable.ExposureTable.makeMinimalSchema()
                        zptExpCatSchema.addField('visit', type='I', doc='Visit number')
                        zptExpCatSchema.addField('detector_id', type='I', doc='Detector number')

                    # And start a new one
                    zptVisitCatalog = afwTable.ExposureCatalog(zptExpCatSchema)
                    zptVisitCatalog.resize(len(camera))

                    # Reset the counter
                    zptCounter = 0

                    lastVisit = int(rec['visit'])

                zptVisitCatalog[zptCounter].setPhotoCalib(photoCalib)
                zptVisitCatalog[zptCounter]['visit'] = int(rec['visit'])
                zptVisitCatalog[zptCounter]['detector_id'] = int(rec['detector'])

                zptCounter += 1

        # Final output of last exposure catalog
        if returnCatalogs:
            yield (int(lastVisit), zptVisitCatalog)

    def _getChebyshevBoundedField(self, coefficients, xyMax, offset=0.0, scaling=1.0):
        """
        Make a ChebyshevBoundedField from fgcm coefficients, with optional offset
        and scaling.

        Parameters
        ----------
        coefficients: `numpy.array`
           Flattened array of chebyshev coefficients
        xyMax: `list` of length 2
           Maximum x and y of the chebyshev bounding box
        offset: `float`, optional
           Absolute calibration offset.  Default is 0.0
        scaling: `float`, optional
           Flat scaling value from fgcmBuildStars.  Default is 1.0

        Returns
        -------
        boundedField: `lsst.afw.math.ChebyshevBoundedField`
        """

        orderPlus1 = int(np.sqrt(coefficients.size))
        pars = np.zeros((orderPlus1, orderPlus1))

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0.0, 0.0),
                               lsst.geom.Point2I(*xyMax))

        pars[:, :] = (coefficients.reshape(orderPlus1, orderPlus1)
                      * (10.**(offset/-2.5))*scaling)

        boundedField = afwMath.ChebyshevBoundedField(bbox, pars)

        return boundedField

    def _outputAtmospheres(self, dataRefDict, atmCat):
        """
        Output the atmospheres.

        Parameters
        ----------
        dataRefDict : `dict`
            All dataRefs are `lsst.daf.persistence.ButlerDataRef` (gen2) or
            `lsst.daf.butler.DeferredDatasetHandle` (gen3)
            dataRef dictionary with keys:

            ``"fgcmLookUpTable"``
                dataRef for the FGCM look-up table.
        atmCat : `lsst.afw.table.BaseCatalog`
            FGCM atmosphere parameter catalog from fgcmFitCycleTask.

        Returns
        -------
        atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
            Generator that returns (visit, transmissionCurve) tuples.
        """
        # First, we need to grab the look-up table and key info
        lutCat = dataRefDict['fgcmLookUpTable'].get()

        atmosphereTableName = lutCat[0]['tablename']
        elevation = lutCat[0]['elevation']
        atmLambda = lutCat[0]['atmLambda']
        lutCat = None

        # Make the atmosphere table if possible
        try:
            atmTable = fgcm.FgcmAtmosphereTable.initWithTableName(atmosphereTableName)
            atmTable.loadTable()
        except IOError:
            atmTable = None

        if atmTable is None:
            # Try to use MODTRAN instead
            try:
                modGen = fgcm.ModtranGenerator(elevation)
                lambdaRange = np.array([atmLambda[0], atmLambda[-1]])/10.
                lambdaStep = (atmLambda[1] - atmLambda[0])/10.
            except (ValueError, IOError) as e:
                raise RuntimeError("FGCM look-up-table generated with modtran, "
                                   "but modtran not configured to run.") from e

        zenith = np.degrees(np.arccos(1./atmCat['secZenith']))

        for i, visit in enumerate(atmCat['visit']):
            if atmTable is not None:
                # Interpolate the atmosphere table
                atmVals = atmTable.interpolateAtmosphere(pmb=atmCat[i]['pmb'],
                                                         pwv=atmCat[i]['pwv'],
                                                         o3=atmCat[i]['o3'],
                                                         tau=atmCat[i]['tau'],
                                                         alpha=atmCat[i]['alpha'],
                                                         zenith=zenith[i],
                                                         ctranslamstd=[atmCat[i]['cTrans'],
                                                                       atmCat[i]['lamStd']])
            else:
                # Run modtran
                modAtm = modGen(pmb=atmCat[i]['pmb'],
                                pwv=atmCat[i]['pwv'],
                                o3=atmCat[i]['o3'],
                                tau=atmCat[i]['tau'],
                                alpha=atmCat[i]['alpha'],
                                zenith=zenith[i],
                                lambdaRange=lambdaRange,
                                lambdaStep=lambdaStep,
                                ctranslamstd=[atmCat[i]['cTrans'],
                                              atmCat[i]['lamStd']])
                atmVals = modAtm['COMBINED']

            # Now need to create something to persist...
            curve = TransmissionCurve.makeSpatiallyConstant(throughput=atmVals,
                                                            wavelengths=atmLambda,
                                                            throughputAtMin=atmVals[0],
                                                            throughputAtMax=atmVals[-1])

            yield (int(visit), curve)
