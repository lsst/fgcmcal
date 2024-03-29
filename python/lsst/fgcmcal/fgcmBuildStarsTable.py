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
"""Build star observations for input to FGCM using sourceTable_visit.

This task finds all the visits and sourceTable_visits in a repository (or a
subset based on command line parameters) and extracts all the potential
calibration stars for input into fgcm.  This task additionally uses fgcm to
match star observations into unique stars, and performs as much cleaning of the
input catalog as possible.
"""

import time
import warnings

import numpy as np
import collections

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
import lsst.afw.table as afwTable
from lsst.meas.algorithms import ReferenceObjectLoader, LoadReferenceObjectsConfig

from .fgcmBuildStarsBase import FgcmBuildStarsConfigBase, FgcmBuildStarsBaseTask
from .utilities import computeApproxPixelAreaFields, computeApertureRadiusFromName

__all__ = ['FgcmBuildStarsTableConfig', 'FgcmBuildStarsTableTask']


class FgcmBuildStarsTableConnections(pipeBase.PipelineTaskConnections,
                                     dimensions=("instrument",),
                                     defaultTemplates={}):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
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

    sourceSchema = connectionTypes.InitInput(
        doc="Schema for source catalogs",
        name="src_schema",
        storageClass="SourceCatalog",
    )

    refCat = connectionTypes.PrerequisiteInput(
        doc="Reference catalog to use for photometric calibration",
        name="cal_ref_cat",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    sourceTable_visit = connectionTypes.Input(
        doc="Source table in parquet format, per visit",
        name="sourceTable_visit",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
        deferLoad=True,
        multiple=True,
    )

    visitSummary = connectionTypes.Input(
        doc=("Per-visit consolidated exposure metadata.  These catalogs use "
             "detector id for the id and must be sorted for fast lookups of a "
             "detector."),
        name="visitSummary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
        deferLoad=True,
        multiple=True,
    )

    fgcmVisitCatalog = connectionTypes.Output(
        doc="Catalog of visit information for fgcm",
        name="fgcmVisitCatalog",
        storageClass="Catalog",
        dimensions=("instrument",),
    )

    fgcmStarObservations = connectionTypes.Output(
        doc="Catalog of star observations for fgcm",
        name="fgcmStarObservations",
        storageClass="Catalog",
        dimensions=("instrument",),
    )

    fgcmStarIds = connectionTypes.Output(
        doc="Catalog of fgcm calibration star IDs",
        name="fgcmStarIds",
        storageClass="Catalog",
        dimensions=("instrument",),
    )

    fgcmStarIndices = connectionTypes.Output(
        doc="Catalog of fgcm calibration star indices",
        name="fgcmStarIndices",
        storageClass="Catalog",
        dimensions=("instrument",),
    )

    fgcmReferenceStars = connectionTypes.Output(
        doc="Catalog of fgcm-matched reference stars",
        name="fgcmReferenceStars",
        storageClass="Catalog",
        dimensions=("instrument",),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doReferenceMatches:
            self.prerequisiteInputs.remove("refCat")
            self.prerequisiteInputs.remove("fgcmLookUpTable")

        if not config.doReferenceMatches:
            self.outputs.remove("fgcmReferenceStars")

    def getSpatialBoundsConnections(self):
        return ("visitSummary",)


class FgcmBuildStarsTableConfig(FgcmBuildStarsConfigBase, pipeBase.PipelineTaskConfig,
                                pipelineConnections=FgcmBuildStarsTableConnections):
    """Config for FgcmBuildStarsTableTask"""

    referenceCCD = pexConfig.Field(
        doc="Reference CCD for checking PSF and background",
        dtype=int,
        default=40,
    )

    def setDefaults(self):
        super().setDefaults()

        # The names here correspond to the post-transformed
        # sourceTable_visit catalogs, which differ from the raw src
        # catalogs.  Therefore, all field and flag names cannot
        # be derived from the base config class.
        self.instFluxField = 'apFlux_12_0_instFlux'
        self.localBackgroundFluxField = 'localBackground_instFlux'
        self.apertureInnerInstFluxField = 'apFlux_12_0_instFlux'
        self.apertureOuterInstFluxField = 'apFlux_17_0_instFlux'
        self.psfCandidateName = 'calib_psf_candidate'

        sourceSelector = self.sourceSelector["science"]

        fluxFlagName = self.instFluxField[0: -len('instFlux')] + 'flag'

        sourceSelector.flags.bad = ['pixelFlags_edge',
                                    'pixelFlags_interpolatedCenter',
                                    'pixelFlags_saturatedCenter',
                                    'pixelFlags_crCenter',
                                    'pixelFlags_bad',
                                    'pixelFlags_interpolated',
                                    'pixelFlags_saturated',
                                    'centroid_flag',
                                    fluxFlagName]

        if self.doSubtractLocalBackground:
            localBackgroundFlagName = self.localBackgroundFluxField[0: -len('instFlux')] + 'flag'
            sourceSelector.flags.bad.append(localBackgroundFlagName)

        sourceSelector.signalToNoise.fluxField = self.instFluxField
        sourceSelector.signalToNoise.errField = self.instFluxField + 'Err'

        sourceSelector.isolated.parentName = 'parentSourceId'
        sourceSelector.isolated.nChildName = 'deblend_nChild'

        sourceSelector.requireFiniteRaDec.raColName = 'ra'
        sourceSelector.requireFiniteRaDec.decColName = 'dec'

        sourceSelector.unresolved.name = 'extendedness'

        sourceSelector.doRequirePrimary = True


class FgcmBuildStarsTableTask(FgcmBuildStarsBaseTask):
    """
    Build stars for the FGCM global calibration, using sourceTable_visit catalogs.
    """
    ConfigClass = FgcmBuildStarsTableConfig
    _DefaultName = "fgcmBuildStarsTable"

    canMultiprocess = False

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(initInputs=initInputs, **kwargs)
        if initInputs is not None:
            self.sourceSchema = initInputs["sourceSchema"].schema

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputRefDict = butlerQC.get(inputRefs)

        sourceTableHandles = inputRefDict['sourceTable_visit']

        self.log.info("Running with %d sourceTable_visit handles",
                      len(sourceTableHandles))

        sourceTableHandleDict = {sourceTableHandle.dataId['visit']: sourceTableHandle for
                                 sourceTableHandle in sourceTableHandles}

        if self.config.doReferenceMatches:
            # Get the LUT handle
            lutHandle = inputRefDict['fgcmLookUpTable']

            # Prepare the reference catalog loader
            refConfig = LoadReferenceObjectsConfig()
            refConfig.filterMap = self.config.fgcmLoadReferenceCatalog.filterMap
            refObjLoader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                          for ref in inputRefs.refCat],
                                                 refCats=butlerQC.get(inputRefs.refCat),
                                                 name=self.config.connections.refCat,
                                                 log=self.log,
                                                 config=refConfig)
            self.makeSubtask('fgcmLoadReferenceCatalog',
                             refObjLoader=refObjLoader,
                             refCatName=self.config.connections.refCat)
        else:
            lutHandle = None

        # Compute aperture radius if necessary.  This is useful to do now before
        # any heave lifting has happened (fail early).
        calibFluxApertureRadius = None
        if self.config.doSubtractLocalBackground:
            try:
                calibFluxApertureRadius = computeApertureRadiusFromName(self.config.instFluxField)
            except RuntimeError as e:
                raise RuntimeError("Could not determine aperture radius from %s. "
                                   "Cannot use doSubtractLocalBackground." %
                                   (self.config.instFluxField)) from e

        visitSummaryHandles = inputRefDict['visitSummary']
        visitSummaryHandleDict = {visitSummaryHandle.dataId['visit']: visitSummaryHandle for
                                  visitSummaryHandle in visitSummaryHandles}

        camera = inputRefDict['camera']
        groupedHandles = self._groupHandles(sourceTableHandleDict,
                                            visitSummaryHandleDict)

        visitCat = self.fgcmMakeVisitCatalog(camera, groupedHandles)

        rad = calibFluxApertureRadius
        fgcmStarObservationCat = self.fgcmMakeAllStarObservations(groupedHandles,
                                                                  visitCat,
                                                                  self.sourceSchema,
                                                                  camera,
                                                                  calibFluxApertureRadius=rad)

        butlerQC.put(visitCat, outputRefs.fgcmVisitCatalog)
        butlerQC.put(fgcmStarObservationCat, outputRefs.fgcmStarObservations)

        fgcmStarIdCat, fgcmStarIndicesCat, fgcmRefCat = self.fgcmMatchStars(visitCat,
                                                                            fgcmStarObservationCat,
                                                                            lutHandle=lutHandle)

        butlerQC.put(fgcmStarIdCat, outputRefs.fgcmStarIds)
        butlerQC.put(fgcmStarIndicesCat, outputRefs.fgcmStarIndices)
        if fgcmRefCat is not None:
            butlerQC.put(fgcmRefCat, outputRefs.fgcmReferenceStars)

    def _groupHandles(self, sourceTableHandleDict, visitSummaryHandleDict):
        """Group sourceTable and visitSummary handles.

        Parameters
        ----------
        sourceTableHandleDict : `dict` [`int`, `str`]
            Dict of source tables, keyed by visit.
        visitSummaryHandleDict : `dict` [int, `str`]
            Dict of visit summary catalogs, keyed by visit.

        Returns
        -------
        groupedHandles : `dict` [`int`, `list`]
            Dictionary with sorted visit keys, and `list`s with
            `lsst.daf.butler.DeferredDataSetHandle`.  The first
            item in the list will be the visitSummary ref, and
            the second will be the source table ref.
        """
        groupedHandles = collections.defaultdict(list)
        visits = sorted(sourceTableHandleDict.keys())

        for visit in visits:
            groupedHandles[visit] = [visitSummaryHandleDict[visit],
                                     sourceTableHandleDict[visit]]

        return groupedHandles

    def fgcmMakeAllStarObservations(self, groupedHandles, visitCat,
                                    sourceSchema,
                                    camera,
                                    calibFluxApertureRadius=None):
        startTime = time.time()

        if self.config.doSubtractLocalBackground and calibFluxApertureRadius is None:
            raise RuntimeError("Must set calibFluxApertureRadius if doSubtractLocalBackground is True.")

        # To get the correct output schema, we use the legacy code.
        # We are not actually using this mapper, except to grab the outputSchema
        sourceMapper = self._makeSourceMapper(sourceSchema)
        outputSchema = sourceMapper.getOutputSchema()

        # Construct mapping from ccd number to index
        ccdMapping = {}
        for ccdIndex, detector in enumerate(camera):
            ccdMapping[detector.getId()] = ccdIndex

        approxPixelAreaFields = computeApproxPixelAreaFields(camera)

        fullCatalog = afwTable.BaseCatalog(outputSchema)

        visitKey = outputSchema['visit'].asKey()
        ccdKey = outputSchema['ccd'].asKey()
        instMagKey = outputSchema['instMag'].asKey()
        instMagErrKey = outputSchema['instMagErr'].asKey()
        deltaMagAperKey = outputSchema['deltaMagAper'].asKey()

        # Prepare local background if desired
        if self.config.doSubtractLocalBackground:
            localBackgroundArea = np.pi*calibFluxApertureRadius**2.

        columns = None

        k = 2.5/np.log(10.)

        for counter, visit in enumerate(visitCat):
            expTime = visit['exptime']

            handle = groupedHandles[visit['visit']][-1]

            if columns is None:
                inColumns = handle.get(component='columns')
                columns = self._get_sourceTable_visit_columns(inColumns)
            df = handle.get(parameters={'columns': columns})

            goodSrc = self.sourceSelector.selectSources(df)

            # Need to add a selection based on the local background correction
            # if necessary
            if self.config.doSubtractLocalBackground:
                localBackground = localBackgroundArea*df[self.config.localBackgroundFluxField].values
                use, = np.where((goodSrc.selected)
                                & ((df[self.config.instFluxField].values - localBackground) > 0.0))
            else:
                use, = np.where(goodSrc.selected)

            tempCat = afwTable.BaseCatalog(fullCatalog.schema)
            tempCat.resize(use.size)

            tempCat['ra'][:] = np.deg2rad(df['ra'].values[use])
            tempCat['dec'][:] = np.deg2rad(df['dec'].values[use])
            tempCat['x'][:] = df['x'].values[use]
            tempCat['y'][:] = df['y'].values[use]
            # The "visit" name in the parquet table is hard-coded.
            tempCat[visitKey][:] = df['visit'].values[use]
            tempCat[ccdKey][:] = df['detector'].values[use]
            tempCat['psf_candidate'] = df[self.config.psfCandidateName].values[use]

            with warnings.catch_warnings():
                # Ignore warnings, we will filter infinites and nans below
                warnings.simplefilter("ignore")

                instMagInner = -2.5*np.log10(df[self.config.apertureInnerInstFluxField].values[use])
                instMagErrInner = k*(df[self.config.apertureInnerInstFluxField + 'Err'].values[use]
                                     / df[self.config.apertureInnerInstFluxField].values[use])
                instMagOuter = -2.5*np.log10(df[self.config.apertureOuterInstFluxField].values[use])
                instMagErrOuter = k*(df[self.config.apertureOuterInstFluxField + 'Err'].values[use]
                                     / df[self.config.apertureOuterInstFluxField].values[use])
                tempCat[deltaMagAperKey][:] = instMagInner - instMagOuter
                # Set bad values to illegal values for fgcm.
                tempCat[deltaMagAperKey][~np.isfinite(tempCat[deltaMagAperKey][:])] = 99.0

            if self.config.doSubtractLocalBackground:
                # At the moment we only adjust the flux and not the flux
                # error by the background because the error on
                # base_LocalBackground_instFlux is the rms error in the
                # background annulus, not the error on the mean in the
                # background estimate (which is much smaller, by sqrt(n)
                # pixels used to estimate the background, which we do not
                # have access to in this task).  In the default settings,
                # the annulus is sufficiently large such that these
                # additional errors are are negligibly small (much less
                # than a mmag in quadrature).

                # This is the difference between the mag with local background correction
                # and the mag without local background correction.
                tempCat['deltaMagBkg'] = (-2.5*np.log10(df[self.config.instFluxField].values[use]
                                                        - localBackground[use]) -
                                          -2.5*np.log10(df[self.config.instFluxField].values[use]))
            else:
                tempCat['deltaMagBkg'][:] = 0.0

            # Need to loop over ccds here
            for detector in camera:
                ccdId = detector.getId()
                # used index for all observations with a given ccd
                use2 = (tempCat[ccdKey] == ccdId)
                tempCat['jacobian'][use2] = approxPixelAreaFields[ccdId].evaluate(tempCat['x'][use2],
                                                                                  tempCat['y'][use2])
                scaledInstFlux = (df[self.config.instFluxField].values[use[use2]]
                                  * visit['scaling'][ccdMapping[ccdId]])
                tempCat[instMagKey][use2] = (-2.5*np.log10(scaledInstFlux) + 2.5*np.log10(expTime))

            # Compute instMagErr from instFluxErr/instFlux, any scaling
            # will cancel out.
            tempCat[instMagErrKey][:] = k*(df[self.config.instFluxField + 'Err'].values[use]
                                           / df[self.config.instFluxField].values[use])

            # Apply the jacobian if configured
            if self.config.doApplyWcsJacobian:
                tempCat[instMagKey][:] -= 2.5*np.log10(tempCat['jacobian'][:])

            fullCatalog.extend(tempCat)

            deltaOk = (np.isfinite(instMagInner) & np.isfinite(instMagErrInner)
                       & np.isfinite(instMagOuter) & np.isfinite(instMagErrOuter))

            visit['deltaAper'] = np.median(instMagInner[deltaOk] - instMagOuter[deltaOk])
            visit['sources_read'] = True

            self.log.info("  Found %d good stars in visit %d (deltaAper = %0.3f)",
                          use.size, visit['visit'], visit['deltaAper'])

        self.log.info("Found all good star observations in %.2f s" %
                      (time.time() - startTime))

        return fullCatalog

    def _get_sourceTable_visit_columns(self, inColumns):
        """
        Get the sourceTable_visit columns from the config.

        Parameters
        ----------
        inColumns : `list`
            List of columns available in the sourceTable_visit

        Returns
        -------
        columns : `list`
            List of columns to read from sourceTable_visit.
        """
        # Some names are hard-coded in the parquet table.
        columns = ['visit', 'detector',
                   'ra', 'dec', 'x', 'y', self.config.psfCandidateName,
                   self.config.instFluxField, self.config.instFluxField + 'Err',
                   self.config.apertureInnerInstFluxField, self.config.apertureInnerInstFluxField + 'Err',
                   self.config.apertureOuterInstFluxField, self.config.apertureOuterInstFluxField + 'Err']
        if self.sourceSelector.config.doFlags:
            columns.extend(self.sourceSelector.config.flags.bad)
        if self.sourceSelector.config.doUnresolved:
            columns.append(self.sourceSelector.config.unresolved.name)
        if self.sourceSelector.config.doIsolated:
            columns.append(self.sourceSelector.config.isolated.parentName)
            columns.append(self.sourceSelector.config.isolated.nChildName)
        if self.sourceSelector.config.doRequirePrimary:
            columns.append(self.sourceSelector.config.requirePrimary.primaryColName)
        if self.config.doSubtractLocalBackground:
            columns.append(self.config.localBackgroundFluxField)

        return columns
