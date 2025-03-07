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
import copy

import numpy as np
import hpgeom as hpg
from astropy import units
from astropy.table import Table
import esutil

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
from lsst.afw.image import TransmissionCurve
from lsst.meas.algorithms import ReferenceObjectLoader, LoadReferenceObjectsConfig
from lsst.pipe.tasks.photoCal import PhotoCalTask
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.skymap import BaseSkyMap

from .utilities import computeApproxPixelAreaFields
from .utilities import FGCM_ILLEGAL_VALUE

import fgcm

__all__ = ['FgcmOutputProductsConfig', 'FgcmOutputProductsTask']


class FgcmOutputProductsConnections(pipeBase.PipelineTaskConnections,
                                    dimensions=("instrument",),
                                    defaultTemplates={"cycleNumber": "0"}):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )

    skymap = connectionTypes.Input(
        doc="Skymap used for tract sharding of output catalog.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    fgcmLookUpTable = connectionTypes.PrerequisiteInput(
        doc=("Atmosphere + instrument look-up-table for FGCM throughput and "
             "chromatic corrections."),
        name="fgcmLookUpTable",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmVisitCatalog = connectionTypes.Input(
        doc="Catalog of visit information for fgcm",
        name="fgcmVisitCatalog",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStandardStars = connectionTypes.Input(
        doc="Catalog of standard star data from fgcm fit",
        name="fgcm_Cycle{cycleNumber}_StandardStars",
        storageClass="SimpleCatalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmZeropoints = connectionTypes.Input(
        doc="Catalog of zeropoints from fgcm fit",
        name="fgcm_Cycle{cycleNumber}_Zeropoints",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmAtmosphereParameters = connectionTypes.Input(
        doc="Catalog of atmosphere parameters from fgcm fit",
        name="fgcm_Cycle{cycleNumber}_AtmosphereParameters",
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

    fgcmPhotoCalib = connectionTypes.Output(
        doc=("Per-visit photometric calibrations derived from fgcm calibration. "
             "These catalogs use detector id for the id and are sorted for "
             "fast lookups of a detector."),
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

    fgcmTractStars = connectionTypes.Output(
        doc="Per-tract fgcm calibrated stars.",
        name="fgcm_standard_star",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "tract", "skymap"),
        multiple=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if str(int(config.connections.cycleNumber)) != config.connections.cycleNumber:
            raise ValueError("cycleNumber must be of integer format")

        if not config.doReferenceCalibration:
            self.prerequisiteInputs.remove("refCat")
        if not config.doAtmosphereOutput:
            self.inputs.remove("fgcmAtmosphereParameters")
        if not config.doZeropointOutput:
            self.inputs.remove("fgcmZeropoints")
        if not config.doReferenceCalibration:
            self.outputs.remove("fgcmOffsets")
        if not config.doTractStars:
            del self.skymap
            del self.fgcmTractStars

    def getSpatialBoundsConnections(self):
        return ("fgcmPhotoCalib",)


class FgcmOutputProductsConfig(pipeBase.PipelineTaskConfig,
                               pipelineConnections=FgcmOutputProductsConnections):
    """Config for FgcmOutputProductsTask"""

    physicalFilterMap = pexConfig.DictField(
        doc="Mapping from 'physicalFilter' to band.",
        keytype=str,
        itemtype=str,
        default={},
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
    doTractStars = pexConfig.Field(
        doc="Output tract-sharded standard stars?",
        dtype=bool,
        default=True,
        # default=False,
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
        self.photoCal.match.sourceSelection.doRequirePrimary = False


class FgcmOutputProductsTask(pipeBase.PipelineTask):
    """
    Output products from FGCM global calibration.
    """

    ConfigClass = FgcmOutputProductsConfig
    _DefaultName = "fgcmOutputProducts"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        handleDict = {}
        handleDict['camera'] = butlerQC.get(inputRefs.camera)
        handleDict['fgcmLookUpTable'] = butlerQC.get(inputRefs.fgcmLookUpTable)
        handleDict['fgcmVisitCatalog'] = butlerQC.get(inputRefs.fgcmVisitCatalog)
        handleDict['fgcmStandardStars'] = butlerQC.get(inputRefs.fgcmStandardStars)

        if self.config.doZeropointOutput:
            handleDict['fgcmZeropoints'] = butlerQC.get(inputRefs.fgcmZeropoints)
            photoCalibRefDict = {photoCalibRef.dataId['visit']:
                                 photoCalibRef for photoCalibRef in outputRefs.fgcmPhotoCalib}

        if self.config.doAtmosphereOutput:
            handleDict['fgcmAtmosphereParameters'] = butlerQC.get(inputRefs.fgcmAtmosphereParameters)
            atmRefDict = {atmRef.dataId['visit']: atmRef for
                          atmRef in outputRefs.fgcmTransmissionAtmosphere}

        if self.config.doReferenceCalibration:
            refConfig = LoadReferenceObjectsConfig()
            self.refObjLoader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                               for ref in inputRefs.refCat],
                                                      refCats=butlerQC.get(inputRefs.refCat),
                                                      name=self.config.connections.refCat,
                                                      log=self.log,
                                                      config=refConfig)
        else:
            self.refObjLoader = None

        if self.config.doTractStars:
            handleDict['skymap'] = butlerQC.get(inputRefs.skymap)
            tractStarRefDict = {tractStarRef.dataId["tract"]: tractStarRef for
                                tractStarRef in outputRefs.fgcmTractStars}

        struct = self.run(handleDict, self.config.physicalFilterMap)

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

        if self.config.doTractStars:
            self.log.info("Outputting standard stars per-tract.")
            for tractId, catalog in struct.tractStars:
                butlerQC.put(catalog, tractStarRefDict[tractId])

        return

    def run(self, handleDict, physicalFilterMap):
        """Run the output products task.

        Parameters
        ----------
        handleDict : `dict`
            All handles are `lsst.daf.butler.DeferredDatasetHandle`
            handle dictionary with keys:

            ``"camera"``
                Camera object (`lsst.afw.cameraGeom.Camera`)
            ``"fgcmLookUpTable"``
                handle for the FGCM look-up table.
            ``"fgcmVisitCatalog"``
                handle for visit summary catalog.
            ``"fgcmStandardStars"``
                handle for the output standard star catalog.
            ``"fgcmZeropoints"``
                handle for the zeropoint data catalog.
            ``"fgcmAtmosphereParameters"``
                handle for the atmosphere parameter catalog.
            ``"fgcmBuildStarsTableConfig"``
                Config for `lsst.fgcmcal.fgcmBuildStarsTableTask`.
            ``"skymap"``
                Skymap for sharding standard stars (optional).

        physicalFilterMap : `dict`
            Dictionary of mappings from physical filter to FGCM band.

        Returns
        -------
        retStruct : `lsst.pipe.base.Struct`
            Output structure with keys:

            offsets : `np.ndarray`
                Final reference offsets, per band.
            atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
                Generator that returns (visit, transmissionCurve) tuples.
            photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
                Generator that returns (visit, exposureCatalog) tuples.
        """
        stdCat = handleDict['fgcmStandardStars'].get()
        md = stdCat.getMetadata()
        bands = md.getArray('BANDS')

        if self.config.doReferenceCalibration:
            lutCat = handleDict['fgcmLookUpTable'].get()
            offsets = self._computeReferenceOffsets(stdCat, lutCat, physicalFilterMap, bands)
        else:
            offsets = np.zeros(len(bands))

        if self.config.doZeropointOutput:
            zptCat = handleDict['fgcmZeropoints'].get()
            visitCat = handleDict['fgcmVisitCatalog'].get()

            pcgen = self._outputZeropoints(handleDict['camera'], zptCat, visitCat, offsets, bands,
                                           physicalFilterMap)
        else:
            pcgen = None

        if self.config.doAtmosphereOutput:
            atmCat = handleDict['fgcmAtmosphereParameters'].get()
            atmgen = self._outputAtmospheres(handleDict, atmCat)
        else:
            atmgen = None

        if self.config.doTractStars:
            skymap = handleDict['skymap']
            tractstargen = self._outputTractStars(skymap, stdCat)
        else:
            tractstargen = None

        retStruct = pipeBase.Struct(offsets=offsets,
                                    atmospheres=atmgen,
                                    tractStars=tractstargen)
        retStruct.photoCalibCatalogs = pcgen

        return retStruct

    def generateTractOutputProducts(self, handleDict, tract,
                                    visitCat, zptCat, atmCat, stdCat,
                                    fgcmBuildStarsConfig):
        """
        Generate the output products for a given tract, as specified in the config.

        This method is here to have an alternate entry-point for
        FgcmCalibrateTract.

        Parameters
        ----------
        handleDict : `dict`
            All handles are `lsst.daf.butler.DeferredDatasetHandle`
            handle dictionary with keys:

            ``"camera"``
                Camera object (`lsst.afw.cameraGeom.Camera`)
            ``"fgcmLookUpTable"``
                handle for the FGCM look-up table.
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

        Returns
        -------
        retStruct : `lsst.pipe.base.Struct`
            Output structure with keys:

            offsets : `np.ndarray`
                Final reference offsets, per band.
            atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
                Generator that returns (visit, transmissionCurve) tuples.
            photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
                Generator that returns (visit, exposureCatalog) tuples.
        """
        physicalFilterMap = fgcmBuildStarsConfig.physicalFilterMap

        md = stdCat.getMetadata()
        bands = md.getArray('BANDS')

        if self.config.doComposeWcsJacobian and not fgcmBuildStarsConfig.doApplyWcsJacobian:
            raise RuntimeError("Cannot compose the WCS jacobian if it hasn't been applied "
                               "in fgcmBuildStarsTask.")

        if not self.config.doComposeWcsJacobian and fgcmBuildStarsConfig.doApplyWcsJacobian:
            self.log.warning("Jacobian was applied in build-stars but doComposeWcsJacobian is not set.")

        if self.config.doReferenceCalibration:
            lutCat = handleDict['fgcmLookUpTable'].get()
            offsets = self._computeReferenceOffsets(stdCat, lutCat, bands, physicalFilterMap)
        else:
            offsets = np.zeros(len(bands))

        if self.config.doZeropointOutput:
            pcgen = self._outputZeropoints(handleDict['camera'], zptCat, visitCat, offsets, bands,
                                           physicalFilterMap)
        else:
            pcgen = None

        if self.config.doAtmosphereOutput:
            atmgen = self._outputAtmospheres(handleDict, atmCat)
        else:
            atmgen = None

        retStruct = pipeBase.Struct(offsets=offsets,
                                    atmospheres=atmgen)
        retStruct.photoCalibCatalogs = pcgen

        return retStruct

    def _computeReferenceOffsets(self, stdCat, lutCat, physicalFilterMap, bands):
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
        lutCat : `lsst.afw.table.SimpleCatalog`
            FGCM Look-up table
        physicalFilterMap : `dict`
            Dictionary of mappings from physical filter to FGCM band.
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

        # Associate each band with the appropriate physicalFilter and make
        # filterLabels
        filterLabels = []

        lutPhysicalFilters = lutCat[0]['physicalFilters'].split(',')
        lutStdPhysicalFilters = lutCat[0]['stdPhysicalFilters'].split(',')
        physicalFilterMapBands = list(physicalFilterMap.values())
        physicalFilterMapFilters = list(physicalFilterMap.keys())
        for band in bands:
            # Find a physical filter associated from the band by doing
            # a reverse lookup on the physicalFilterMap dict
            physicalFilterMapIndex = physicalFilterMapBands.index(band)
            physicalFilter = physicalFilterMapFilters[physicalFilterMapIndex]
            # Find the appropriate fgcm standard physicalFilter
            lutPhysicalFilterIndex = lutPhysicalFilters.index(physicalFilter)
            stdPhysicalFilter = lutStdPhysicalFilters[lutPhysicalFilterIndex]
            filterLabels.append(afwImage.FilterLabel(band=band,
                                                     physical=stdPhysicalFilter))

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
        ipring = hpg.angle_to_pixel(
            self.config.referencePixelizationNside,
            stdCat['coord_ra'],
            stdCat['coord_dec'],
            degrees=False,
        )
        h, rev = fgcm.fgcmUtilities.histogram_rev_sorted(ipring)

        gdpix, = np.where(h >= self.config.referencePixelizationMinStars)

        self.log.info("Found %d pixels (nside=%d) with at least %d good stars" %
                      (gdpix.size,
                       self.config.referencePixelizationNside,
                       self.config.referencePixelizationMinStars))

        if gdpix.size < self.config.referencePixelizationNPixels:
            self.log.warning("Found fewer good pixels (%d) than preferred in configuration (%d)" %
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
        selected = np.zeros(len(stdCat), dtype=bool)

        refFluxFields = [None]*len(bands)

        for p_index, pix in enumerate(gdpix):
            i1a = rev[rev[pix]: rev[pix + 1]]

            # the stdCat afwTable can only be indexed with boolean arrays,
            # and not numpy index arrays (see DM-16497).  This little trick
            # converts the index array into a boolean array
            selected[:] = False
            selected[i1a] = True

            for b_index, filterLabel in enumerate(filterLabels):
                struct = self._computeOffsetOneBand(sourceMapper, badStarKey, b_index,
                                                    filterLabel, stdCat,
                                                    selected, refFluxFields)
                results['nstar'][p_index, b_index] = len(i1a)
                results['nmatch'][p_index, b_index] = len(struct.arrays.refMag)
                results['zp'][p_index, b_index] = struct.zp
                results['zpErr'][p_index, b_index] = struct.sigma

        # And compute the summary statistics
        offsets = np.zeros(len(bands))

        for b_index, band in enumerate(bands):
            # make configurable
            ok, = np.where(results['nmatch'][:, b_index] >= self.config.referenceMinMatch)
            offsets[b_index] = np.median(results['zp'][ok, b_index])
            # use median absolute deviation to estimate Normal sigma
            # see https://en.wikipedia.org/wiki/Median_absolute_deviation
            madSigma = 1.4826*np.median(np.abs(results['zp'][ok, b_index] - offsets[b_index]))
            self.log.info("Reference catalog offset for %s band: %.12f +/- %.12f",
                          band, offsets[b_index], madSigma)

        return offsets

    def _computeOffsetOneBand(self, sourceMapper, badStarKey,
                              b_index, filterLabel, stdCat, selected, refFluxFields):
        """
        Compute the zeropoint offset between the fgcm stdCat and the reference
        stars for one pixel in one band

        Parameters
        ----------
        sourceMapper : `lsst.afw.table.SchemaMapper`
            Mapper to go from stdCat to calibratable catalog
        badStarKey : `lsst.afw.table.Key`
            Key for the field with bad stars
        b_index : `int`
            Index of the band in the star catalog
        filterLabel : `lsst.afw.image.FilterLabel`
            filterLabel with band and physical filter
        stdCat : `lsst.afw.table.SimpleCatalog`
            FGCM standard stars
        selected : `numpy.array(dtype=bool)`
            Boolean array of which stars are in the pixel
        refFluxFields : `list`
            List of names of flux fields for reference catalog
        """

        sourceCat = afwTable.SimpleCatalog(sourceMapper.getOutputSchema())
        sourceCat.reserve(selected.sum())
        sourceCat.extend(stdCat[selected], mapper=sourceMapper)
        sourceCat['instFlux'] = 10.**(stdCat['mag_std_noabs'][selected, b_index]/(-2.5))
        sourceCat['instFluxErr'] = (np.log(10.)/2.5)*(stdCat['magErr_std'][selected, b_index]
                                                      * sourceCat['instFlux'])
        # Make sure we only use stars that have valid measurements
        # (This is perhaps redundant with requirements above that the
        # stars be observed in all bands, but it can't hurt)
        badStar = (stdCat['mag_std_noabs'][selected, b_index] > 90.0)
        for rec in sourceCat[badStar]:
            rec.set(badStarKey, True)

        exposure = afwImage.ExposureF()
        exposure.setFilter(filterLabel)

        if refFluxFields[b_index] is None:
            # Need to find the flux field in the reference catalog
            # to work around limitations of DirectMatch in PhotoCal
            ctr = stdCat[0].getCoord()
            rad = 0.05*lsst.geom.degrees
            refDataTest = self.refObjLoader.loadSkyCircle(ctr, rad, filterLabel.bandLabel)
            refFluxFields[b_index] = refDataTest.fluxField

        # Make a copy of the config so that we can modify it
        calConfig = copy.copy(self.config.photoCal.value)
        calConfig.match.referenceSelection.signalToNoise.fluxField = refFluxFields[b_index]
        calConfig.match.referenceSelection.signalToNoise.errField = refFluxFields[b_index] + 'Err'
        calTask = self.config.photoCal.target(refObjLoader=self.refObjLoader,
                                              config=calConfig,
                                              schema=sourceCat.getSchema())

        struct = calTask.run(exposure, sourceCat)

        return struct

    def _outputZeropoints(self, camera, zptCat, visitCat, offsets, bands,
                          physicalFilterMap, tract=None):
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
        physicalFilterMap : `dict`
            Dictionary of mappings from physical filter to FGCM band.
        tract: `int`, optional
            Tract number to output.  Default is None (global calibration)

        Returns
        -------
        photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
            Generator that returns (visit, exposureCatalog) tuples.
        """
        # Select visit/ccds where we have a calibration
        # This includes ccds where we were able to interpolate from neighboring
        # ccds.
        cannot_compute = fgcm.fgcmUtilities.zpFlagDict['CANNOT_COMPUTE_ZEROPOINT']
        selected = (((zptCat['fgcmFlag'] & cannot_compute) == 0)
                    & (zptCat['fgcmZptVar'] > 0.0)
                    & (zptCat['fgcmZpt'] > FGCM_ILLEGAL_VALUE))

        # Log warnings for any visit which has no valid zeropoints
        badVisits = np.unique(zptCat['visit'][~selected])
        goodVisits = np.unique(zptCat['visit'][selected])
        allBadVisits = badVisits[~np.isin(badVisits, goodVisits)]
        for allBadVisit in allBadVisits:
            self.log.warning(f'No suitable photoCalib for visit {allBadVisit}')

        # Get a mapping from filtername to the offsets
        offsetMapping = {}
        for f in physicalFilterMap:
            # Not every filter in the map will necesarily have a band.
            if physicalFilterMap[f] in bands:
                offsetMapping[f] = offsets[bands.index(physicalFilterMap[f])]

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
        zptVisitCatalog = None

        metadata = dafBase.PropertyList()
        metadata.add("COMMENT", "Catalog id is detector id, sorted.")
        metadata.add("COMMENT", "Only detectors with data have entries.")

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

            # Return full per-visit exposure catalogs
            if rec['visit'] != lastVisit:
                # This is a new visit.  If the last visit was not -1, yield
                # the ExposureCatalog
                if lastVisit > -1:
                    # ensure that the detectors are in sorted order, for fast lookups
                    zptVisitCatalog.sort()
                    yield (int(lastVisit), zptVisitCatalog)
                else:
                    # We need to create a new schema
                    zptExpCatSchema = afwTable.ExposureTable.makeMinimalSchema()
                    zptExpCatSchema.addField('visit', type='L', doc='Visit number')

                # And start a new one
                zptVisitCatalog = afwTable.ExposureCatalog(zptExpCatSchema)
                zptVisitCatalog.setMetadata(metadata)

                lastVisit = int(rec['visit'])

            catRecord = zptVisitCatalog.addNew()
            catRecord['id'] = int(rec['detector'])
            catRecord['visit'] = rec['visit']
            catRecord.setPhotoCalib(photoCalib)

        # Final output of last exposure catalog
        # ensure that the detectors are in sorted order, for fast lookups
        zptVisitCatalog.sort()
        yield (int(lastVisit), zptVisitCatalog)

    @staticmethod
    def _getChebyshevBoundedField(coefficients, xyMax, offset=0.0, scaling=1.0):
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

        bbox = lsst.geom.Box2I(minimum=lsst.geom.Point2I(0, 0),
                               maximum=lsst.geom.Point2I(*xyMax))

        pars[:, :] = (coefficients.reshape(orderPlus1, orderPlus1)
                      * (10.**(offset/-2.5))*scaling)

        boundedField = afwMath.ChebyshevBoundedField(bbox, pars)

        return boundedField

    def _outputAtmospheres(self, handleDict, atmCat):
        """
        Output the atmospheres.

        Parameters
        ----------
        handleDict : `dict`
            All data handles are `lsst.daf.butler.DeferredDatasetHandle`
            The handleDict has the follownig keys:

            ``"fgcmLookUpTable"``
                handle for the FGCM look-up table.
        atmCat : `lsst.afw.table.BaseCatalog`
            FGCM atmosphere parameter catalog from fgcmFitCycleTask.

        Returns
        -------
        atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
            Generator that returns (visit, transmissionCurve) tuples.
        """
        # First, we need to grab the look-up table and key info
        lutCat = handleDict['fgcmLookUpTable'].get()

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

    def _outputTractStars(self, skymap, stdCat):
        """Output the tract-sharded stars.

        Parameters
        ----------
        skymap : `lsst.skymap.SkyMap`
            Skymap for tract id information.
        stdCat : `lsst.afw.table.SimpleCatalog`
            FGCM standard star catalog from ``FgcmFitCycleTask``

        Returns
        -------
        tractgen : `generator` [(`int`, `astropy.table.Table`)]
            Generator that returns (tractId, Table) tuples.
        """
        md = stdCat.getMetadata()
        bands = md.getArray('BANDS')

        dtype = [
            ("fgcm_id", "i8"),
            ("isolated_star_id", "i8"),
            ("ra", "f8"),
            ("dec", "f8"),
        ]

        for band in bands:
            dtype.extend(
                (
                    (f"mag_{band}", "f4"),
                    (f"magErr_{band}", "f4"),
                    (f"ngood_{band}", "i4"),
                ),
            )

        tractIds = skymap.findTractIdArray(stdCat["coord_ra"], stdCat["coord_dec"])

        h, rev = esutil.stat.histogram(tractIds, rev=True)

        good, = np.where(h > 0)

        for index in good:
            i1a = rev[rev[index]: rev[index + 1]]
            tractId = tractIds[i1a[0]]

            table = Table(np.zeros(len(i1a), dtype=dtype))
            table["fgcm_id"] = stdCat["id"][i1a]
            table["isolated_star_id"] = stdCat["isolated_star_id"][i1a]
            table["ra"] = np.rad2deg(stdCat["coord_ra"][i1a])*units.degree
            table["dec"] = np.rad2deg(stdCat["coord_dec"][i1a])*units.degree

            for i, band in enumerate(bands):
                table[f"mag_{band}"] = stdCat["mag_std_noabs"][i1a, i]*units.ABmag
                table[f"magErr_{band}"] = stdCat["magErr_std"][i1a, i]*units.ABmag
                table[f"ngood_{band}"] = stdCat["ngood"][i1a, i]

            yield (int(tractId), table)
