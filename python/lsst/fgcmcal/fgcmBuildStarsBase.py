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
"""Base class for BuildStars using src tables or sourceTable_visit tables.
"""

import os
import sys
import traceback
import abc

import numpy as np

import lsst.daf.persistence as dafPersist
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.geom as geom
from lsst.daf.base import PropertyList
from lsst.daf.base.dateTime import DateTime
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .utilities import computeApertureRadiusFromDataRef
from .fgcmLoadReferenceCatalog import FgcmLoadReferenceCatalogTask

import fgcm

REFSTARS_FORMAT_VERSION = 1

__all__ = ['FgcmBuildStarsConfigBase', 'FgcmBuildStarsRunner', 'FgcmBuildStarsBaseTask']


class FgcmBuildStarsConfigBase(pexConfig.Config):
    """Base config for FgcmBuildStars tasks"""

    instFluxField = pexConfig.Field(
        doc=("Faull name of the source instFlux field to use, including 'instFlux'. "
             "The associated flag will be implicitly included in badFlags"),
        dtype=str,
        default='slot_CalibFlux_instFlux',
    )
    minPerBand = pexConfig.Field(
        doc="Minimum observations per band",
        dtype=int,
        default=2,
    )
    matchRadius = pexConfig.Field(
        doc="Match radius (arcseconds)",
        dtype=float,
        default=1.0,
    )
    isolationRadius = pexConfig.Field(
        doc="Isolation radius (arcseconds)",
        dtype=float,
        default=2.0,
    )
    densityCutNside = pexConfig.Field(
        doc="Density cut healpix nside",
        dtype=int,
        default=128,
    )
    densityCutMaxPerPixel = pexConfig.Field(
        doc="Density cut number of stars per pixel",
        dtype=int,
        default=1000,
    )
    matchNside = pexConfig.Field(
        doc="Healpix Nside for matching",
        dtype=int,
        default=4096,
    )
    coarseNside = pexConfig.Field(
        doc="Healpix coarse Nside for partitioning matches",
        dtype=int,
        default=8,
    )
    filterMap = pexConfig.DictField(
        doc="Mapping from 'filterName' to band.",
        keytype=str,
        itemtype=str,
        default={},
        deprecated=("This field is no longer used, and has been deprecated by "
                    "DM-28088.  It will be removed after v22.  Use "
                    "physicalFilterMap instead.")
    )
    # The following config will not be necessary after Gen2 retirement.
    # In the meantime, obs packages should set to 'filterDefinitions.filter_to_band'
    # which is easiest to access in the config file.
    physicalFilterMap = pexConfig.DictField(
        doc="Mapping from 'physicalFilter' to band.",
        keytype=str,
        itemtype=str,
        default={},
    )
    requiredBands = pexConfig.ListField(
        doc="Bands required for each star",
        dtype=str,
        default=(),
    )
    primaryBands = pexConfig.ListField(
        doc=("Bands for 'primary' star matches. "
             "A star must be observed in one of these bands to be considered "
             "as a calibration star."),
        dtype=str,
        default=None
    )
    visitDataRefName = pexConfig.Field(
        doc="dataRef name for the 'visit' field, usually 'visit'.",
        dtype=str,
        default="visit"
    )
    ccdDataRefName = pexConfig.Field(
        doc="dataRef name for the 'ccd' field, usually 'ccd' or 'detector'.",
        dtype=str,
        default="ccd"
    )
    doApplyWcsJacobian = pexConfig.Field(
        doc="Apply the jacobian of the WCS to the star observations prior to fit?",
        dtype=bool,
        default=True
    )
    doModelErrorsWithBackground = pexConfig.Field(
        doc="Model flux errors with background term?",
        dtype=bool,
        default=True
    )
    psfCandidateName = pexConfig.Field(
        doc="Name of field with psf candidate flag for propagation",
        dtype=str,
        default="calib_psf_candidate"
    )
    doSubtractLocalBackground = pexConfig.Field(
        doc=("Subtract the local background before performing calibration? "
             "This is only supported for circular aperture calibration fluxes."),
        dtype=bool,
        default=False
    )
    localBackgroundFluxField = pexConfig.Field(
        doc="Full name of the local background instFlux field to use.",
        dtype=str,
        default='base_LocalBackground_instFlux'
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources",
        default="science"
    )
    apertureInnerInstFluxField = pexConfig.Field(
        doc=("Full name of instFlux field that contains inner aperture "
             "flux for aperture correction proxy"),
        dtype=str,
        default='base_CircularApertureFlux_12_0_instFlux'
    )
    apertureOuterInstFluxField = pexConfig.Field(
        doc=("Full name of instFlux field that contains outer aperture "
             "flux for aperture correction proxy"),
        dtype=str,
        default='base_CircularApertureFlux_17_0_instFlux'
    )
    doReferenceMatches = pexConfig.Field(
        doc="Match reference catalog as additional constraint on calibration",
        dtype=bool,
        default=True,
    )
    fgcmLoadReferenceCatalog = pexConfig.ConfigurableField(
        target=FgcmLoadReferenceCatalogTask,
        doc="FGCM reference object loader",
    )
    nVisitsPerCheckpoint = pexConfig.Field(
        doc="Number of visits read between checkpoints",
        dtype=int,
        default=500,
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["science"]
        sourceSelector.setDefaults()

        sourceSelector.doFlags = True
        sourceSelector.doUnresolved = True
        sourceSelector.doSignalToNoise = True
        sourceSelector.doIsolated = True

        sourceSelector.signalToNoise.minimum = 10.0
        sourceSelector.signalToNoise.maximum = 1000.0

        # FGCM operates on unresolved sources, and this setting is
        # appropriate for the current base_ClassificationExtendedness
        sourceSelector.unresolved.maximum = 0.5


class FgcmBuildStarsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for FgcmBuildStars tasks

    fgcmBuildStarsTask.run() and fgcmBuildStarsTableTask.run() take a number of
    arguments, one of which is the butler (for persistence and mapper data),
    and a list of dataRefs extracted from the command line.  Note that FGCM
    runs on a large set of dataRefs, and not on single dataRef/tract/patch.
    This class transforms the process arguments generated by the ArgumentParser
    into the arguments expected by FgcmBuildStarsTask.run().  This runner does
    not use any parallelization.
    """
    @staticmethod
    def getTargetList(parsedCmd):
        """
        Return a list with one element: a tuple with the butler and
        list of dataRefs
        """
        # we want to combine the butler with any (or no!) dataRefs
        return [(parsedCmd.butler, parsedCmd.id.refList)]

    def __call__(self, args):
        """
        Parameters
        ----------
        args: `tuple` with (butler, dataRefList)

        Returns
        -------
        exitStatus: `list` with `lsst.pipe.base.Struct`
           exitStatus (0: success; 1: failure)
        """
        butler, dataRefList = args

        task = self.TaskClass(config=self.config, log=self.log)

        exitStatus = 0
        if self.doRaise:
            task.runDataRef(butler, dataRefList)
        else:
            try:
                task.runDataRef(butler, dataRefList)
            except Exception as e:
                exitStatus = 1
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)

        # The task does not return any results:
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
            resultList = self(targetList[0])

        return resultList


class FgcmBuildStarsBaseTask(pipeBase.PipelineTask, pipeBase.CmdLineTask, abc.ABC):
    """
    Base task to build stars for FGCM global calibration

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
    """
    def __init__(self, butler=None, initInputs=None, **kwargs):
        super().__init__(**kwargs)

        self.makeSubtask("sourceSelector")
        # Only log warning and fatal errors from the sourceSelector
        self.sourceSelector.log.setLevel(self.sourceSelector.log.WARN)

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def runDataRef(self, butler, dataRefs):
        """
        Cross-match and make star list for FGCM Input

        Parameters
        ----------
        butler:  `lsst.daf.persistence.Butler`
        dataRefs: `list` of `lsst.daf.persistence.ButlerDataRef`
           Source data references for the input visits.

        Raises
        ------
        RuntimeErrror: Raised if `config.doReferenceMatches` is set and
           an fgcmLookUpTable is not available, or if computeFluxApertureRadius()
           fails if the calibFlux is not a CircularAperture flux.
        """
        datasetType = dataRefs[0].butlerSubset.datasetType
        self.log.info("Running with %d %s dataRefs", len(dataRefs), datasetType)

        if self.config.doReferenceMatches:
            self.makeSubtask("fgcmLoadReferenceCatalog", butler=butler)
            # Ensure that we have a LUT
            if not butler.datasetExists('fgcmLookUpTable'):
                raise RuntimeError("Must have fgcmLookUpTable if using config.doReferenceMatches")
        # Compute aperture radius if necessary.  This is useful to do now before
        # any heavy lifting has happened (fail early).
        calibFluxApertureRadius = None
        if self.config.doSubtractLocalBackground:
            try:
                calibFluxApertureRadius = computeApertureRadiusFromDataRef(dataRefs[0],
                                                                           self.config.instFluxField)
            except RuntimeError as e:
                raise RuntimeError("Could not determine aperture radius from %s. "
                                   "Cannot use doSubtractLocalBackground." %
                                   (self.config.instFluxField)) from e

        camera = butler.get('camera')
        groupedDataRefs = self._findAndGroupDataRefsGen2(butler, camera, dataRefs)

        # Make the visit catalog if necessary
        # First check if the visit catalog is in the _current_ path
        # We cannot use Gen2 datasetExists() because that checks all parent
        # directories as well, which would make recovering from faults
        # and fgcmcal reruns impossible.
        visitCatDataRef = butler.dataRef('fgcmVisitCatalog')
        filename = visitCatDataRef.get('fgcmVisitCatalog_filename')[0]
        if os.path.exists(filename):
            # This file exists and we should continue processing
            inVisitCat = visitCatDataRef.get()
            if len(inVisitCat) != len(groupedDataRefs):
                raise RuntimeError("Existing visitCatalog found, but has an inconsistent "
                                   "number of visits.  Cannot continue.")
        else:
            inVisitCat = None

        visitCat = self.fgcmMakeVisitCatalog(camera, groupedDataRefs,
                                             visitCatDataRef=visitCatDataRef,
                                             inVisitCat=inVisitCat)

        # Persist the visitCat as a checkpoint file.
        visitCatDataRef.put(visitCat)

        starObsDataRef = butler.dataRef('fgcmStarObservations')
        filename = starObsDataRef.get('fgcmStarObservations_filename')[0]
        if os.path.exists(filename):
            inStarObsCat = starObsDataRef.get()
        else:
            inStarObsCat = None

        rad = calibFluxApertureRadius
        sourceSchemaDataRef = butler.dataRef('src_schema')
        fgcmStarObservationCat = self.fgcmMakeAllStarObservations(groupedDataRefs,
                                                                  visitCat,
                                                                  sourceSchemaDataRef,
                                                                  camera,
                                                                  calibFluxApertureRadius=rad,
                                                                  starObsDataRef=starObsDataRef,
                                                                  visitCatDataRef=visitCatDataRef,
                                                                  inStarObsCat=inStarObsCat)
        visitCatDataRef.put(visitCat)
        starObsDataRef.put(fgcmStarObservationCat)

        # Always do the matching.
        if self.config.doReferenceMatches:
            lutDataRef = butler.dataRef('fgcmLookUpTable')
        else:
            lutDataRef = None
        fgcmStarIdCat, fgcmStarIndicesCat, fgcmRefCat = self.fgcmMatchStars(visitCat,
                                                                            fgcmStarObservationCat,
                                                                            lutDataRef=lutDataRef)

        # Persist catalogs via the butler
        butler.put(fgcmStarIdCat, 'fgcmStarIds')
        butler.put(fgcmStarIndicesCat, 'fgcmStarIndices')
        if fgcmRefCat is not None:
            butler.put(fgcmRefCat, 'fgcmReferenceStars')

    @abc.abstractmethod
    def _findAndGroupDataRefsGen2(self, butler, camera, dataRefs):
        """
        Find and group dataRefs (by visit); Gen2 only.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
            Gen2 butler.
        camera : `lsst.afw.cameraGeom.Camera`
            Camera from the butler.
        dataRefs : `list` of `lsst.daf.persistence.ButlerDataRef`
            Data references for the input visits.

        Returns
        -------
        groupedDataRefs : `dict` [`int`, `list`]
            Dictionary with sorted visit keys, and `list`s of
            `lsst.daf.persistence.ButlerDataRef`
        """
        raise NotImplementedError("_findAndGroupDataRefsGen2 not implemented.")

    @abc.abstractmethod
    def fgcmMakeAllStarObservations(self, groupedDataRefs, visitCat,
                                    sourceSchemaDataRef,
                                    camera,
                                    calibFluxApertureRadius=None,
                                    visitCatDataRef=None,
                                    starObsDataRef=None,
                                    inStarObsCat=None):
        """
        Compile all good star observations from visits in visitCat.  Checkpoint files
        will be stored if both visitCatDataRef and starObsDataRef are not None.

        Parameters
        ----------
        groupedDataRefs: `dict` of `list`s
           Lists of `~lsst.daf.persistence.ButlerDataRef` or
           `~lsst.daf.butler.DeferredDatasetHandle`, grouped by visit.
        visitCat: `~afw.table.BaseCatalog`
           Catalog with visit data for FGCM
        sourceSchemaDataRef: `~lsst.daf.persistence.ButlerDataRef` or
                             `~lsst.daf.butler.DeferredDatasetHandle`
           DataRef for the schema of the src catalogs.
        camera: `~lsst.afw.cameraGeom.Camera`
        calibFluxApertureRadius: `float`, optional
           Aperture radius for calibration flux.
        visitCatDataRef: `~lsst.daf.persistence.ButlerDataRef`, optional
           Dataref to write visitCat for checkpoints
        starObsDataRef: `~lsst.daf.persistence.ButlerDataRef`, optional
           Dataref to write the star observation catalog for checkpoints.
        inStarObsCat: `~afw.table.BaseCatalog`
           Input observation catalog.  If this is incomplete, observations
           will be appended from when it was cut off.

        Returns
        -------
        fgcmStarObservations: `afw.table.BaseCatalog`
           Full catalog of good observations.

        Raises
        ------
        RuntimeError: Raised if doSubtractLocalBackground is True and
           calibFluxApertureRadius is not set.
        """
        raise NotImplementedError("fgcmMakeAllStarObservations not implemented.")

    def fgcmMakeVisitCatalog(self, camera, groupedDataRefs, bkgDataRefDict=None,
                             visitCatDataRef=None, inVisitCat=None):
        """
        Make a visit catalog with all the keys from each visit

        Parameters
        ----------
        camera: `lsst.afw.cameraGeom.Camera`
           Camera from the butler
        groupedDataRefs: `dict`
           Dictionary with visit keys, and `list`s of
           `lsst.daf.persistence.ButlerDataRef`
        bkgDataRefDict: `dict`, optional
           Dictionary of gen3 dataRefHandles for background info.
        visitCatDataRef: `lsst.daf.persistence.ButlerDataRef`, optional
           Dataref to write visitCat for checkpoints
        inVisitCat: `afw.table.BaseCatalog`, optional
           Input (possibly incomplete) visit catalog

        Returns
        -------
        visitCat: `afw.table.BaseCatalog`
        """

        self.log.info("Assembling visitCatalog from %d %ss" %
                      (len(groupedDataRefs), self.config.visitDataRefName))

        nCcd = len(camera)

        if inVisitCat is None:
            schema = self._makeFgcmVisitSchema(nCcd)

            visitCat = afwTable.BaseCatalog(schema)
            visitCat.reserve(len(groupedDataRefs))
            visitCat.resize(len(groupedDataRefs))

            visitCat['visit'] = list(groupedDataRefs.keys())
            visitCat['used'] = 0
            visitCat['sources_read'] = False
        else:
            visitCat = inVisitCat

        # No matter what, fill the catalog. This will check if it was
        # already read.
        self._fillVisitCatalog(visitCat, groupedDataRefs,
                               bkgDataRefDict=bkgDataRefDict,
                               visitCatDataRef=visitCatDataRef)

        return visitCat

    def _fillVisitCatalog(self, visitCat, groupedDataRefs, bkgDataRefDict=None,
                          visitCatDataRef=None):
        """
        Fill the visit catalog with visit metadata

        Parameters
        ----------
        visitCat: `afw.table.BaseCatalog`
           Catalog with schema from _makeFgcmVisitSchema()
        groupedDataRefs: `dict`
           Dictionary with visit keys, and `list`s of
           `lsst.daf.persistence.ButlerDataRef`
        visitCatDataRef: `lsst.daf.persistence.ButlerDataRef`, optional
           Dataref to write visitCat for checkpoints
        bkgDataRefDict: `dict`, optional
           Dictionary of gen3 dataRefHandles for background info. FIXME
        """
        bbox = geom.BoxI(geom.PointI(0, 0), geom.PointI(1, 1))

        for i, visit in enumerate(groupedDataRefs):
            # We don't use the bypasses since we need the psf info which does
            # not have a bypass
            # TODO: When DM-15500 is implemented in the Gen3 Butler, this
            # can be fixed

            # Do not read those that have already been read
            if visitCat['used'][i]:
                continue

            if (i % self.config.nVisitsPerCheckpoint) == 0:
                self.log.info("Retrieving metadata for %s %d (%d/%d)" %
                              (self.config.visitDataRefName, visit, i, len(groupedDataRefs)))
                # Save checkpoint if desired
                if visitCatDataRef is not None:
                    visitCatDataRef.put(visitCat)

            dataRef = groupedDataRefs[visit][0]
            if isinstance(dataRef, dafPersist.ButlerDataRef):
                # Gen2: calexp dataRef
                # The first dataRef in the group will be the reference ccd (if available)
                exp = dataRef.get(datasetType='calexp_sub', bbox=bbox)
                visitInfo = exp.getInfo().getVisitInfo()
                label = dataRef.get(datasetType='calexp_filterLabel')
                physicalFilter = label.physicalLabel
                psf = exp.getPsf()
                psfSigma = psf.computeShape().getDeterminantRadius()
            else:
                # Gen3: use the visitSummary dataRef
                summary = dataRef.get()

                summaryRow = summary.find(self.config.referenceCCD)
                if summaryRow is None:
                    # Take the first available ccd if reference isn't available
                    summaryRow = summary[0]

                visitInfo = summaryRow.getVisitInfo()
                physicalFilter = summaryRow['physical_filter']
                # Compute the median psf sigma if possible
                goodSigma, = np.where(summary['psfSigma'] > 0)
                if goodSigma.size > 2:
                    psfSigma = np.median(summary['psfSigma'][goodSigma])
                elif goodSigma > 0:
                    psfSigma = np.mean(summary['psfSigma'][goodSigma])
                else:
                    psfSigma = 0.0

            rec = visitCat[i]
            rec['visit'] = visit
            rec['physicalFilter'] = physicalFilter
            # TODO DM-26991: when gen2 is removed, gen3 workflow will make it
            # much easier to get the wcs's necessary to recompute the pointing
            # ra/dec at the center of the camera.
            radec = visitInfo.getBoresightRaDec()
            rec['telra'] = radec.getRa().asDegrees()
            rec['teldec'] = radec.getDec().asDegrees()
            rec['telha'] = visitInfo.getBoresightHourAngle().asDegrees()
            rec['telrot'] = visitInfo.getBoresightRotAngle().asDegrees()
            rec['mjd'] = visitInfo.getDate().get(system=DateTime.MJD)
            rec['exptime'] = visitInfo.getExposureTime()
            # convert from Pa to millibar
            # Note that I don't know if this unit will need to be per-camera config
            rec['pmb'] = visitInfo.getWeather().getAirPressure() / 100
            # Flag to signify if this is a "deep" field.  Not currently used
            rec['deepFlag'] = 0
            # Relative flat scaling (1.0 means no relative scaling)
            rec['scaling'][:] = 1.0
            # Median delta aperture, to be measured from stars
            rec['deltaAper'] = 0.0
            rec['psfSigma'] = psfSigma

            if self.config.doModelErrorsWithBackground:
                foundBkg = False
                if isinstance(dataRef, dafPersist.ButlerDataRef):
                    det = dataRef.dataId[self.config.ccdDataRefName]
                    if dataRef.datasetExists(datasetType='calexpBackground'):
                        bgList = dataRef.get(datasetType='calexpBackground')
                        foundBkg = True
                else:
                    det = dataRef.dataId['detector']
                    try:
                        bkgRef = bkgDataRefDict[(visit, det)]
                        bgList = bkgRef.get()
                        foundBkg = True
                    except KeyError:
                        pass

                if foundBkg:
                    bgStats = (bg[0].getStatsImage().getImage().array
                               for bg in bgList)
                    rec['skyBackground'] = sum(np.median(bg[np.isfinite(bg)]) for bg in bgStats)
                else:
                    self.log.warn('Sky background not found for visit %d / ccd %d' %
                                  (visit, det))
                    rec['skyBackground'] = -1.0
            else:
                rec['skyBackground'] = -1.0

            rec['used'] = 1

    def _makeSourceMapper(self, sourceSchema):
        """
        Make a schema mapper for fgcm sources

        Parameters
        ----------
        sourceSchema: `afwTable.Schema`
           Default source schema from the butler

        Returns
        -------
        sourceMapper: `afwTable.schemaMapper`
           Mapper to the FGCM source schema
        """

        # create a mapper to the preferred output
        sourceMapper = afwTable.SchemaMapper(sourceSchema)

        # map to ra/dec
        sourceMapper.addMapping(sourceSchema['coord_ra'].asKey(), 'ra')
        sourceMapper.addMapping(sourceSchema['coord_dec'].asKey(), 'dec')
        sourceMapper.addMapping(sourceSchema['slot_Centroid_x'].asKey(), 'x')
        sourceMapper.addMapping(sourceSchema['slot_Centroid_y'].asKey(), 'y')
        # Add the mapping if the field exists in the input catalog.
        # If the field does not exist, simply add it (set to False).
        # This field is not required for calibration, but is useful
        # to collate if available.
        try:
            sourceMapper.addMapping(sourceSchema[self.config.psfCandidateName].asKey(),
                                    'psf_candidate')
        except LookupError:
            sourceMapper.editOutputSchema().addField(
                "psf_candidate", type='Flag',
                doc=("Flag set if the source was a candidate for PSF determination, "
                     "as determined by the star selector."))

        # and add the fields we want
        sourceMapper.editOutputSchema().addField(
            "visit", type=np.int32, doc="Visit number")
        sourceMapper.editOutputSchema().addField(
            "ccd", type=np.int32, doc="CCD number")
        sourceMapper.editOutputSchema().addField(
            "instMag", type=np.float32, doc="Instrumental magnitude")
        sourceMapper.editOutputSchema().addField(
            "instMagErr", type=np.float32, doc="Instrumental magnitude error")
        sourceMapper.editOutputSchema().addField(
            "jacobian", type=np.float32, doc="Relative pixel scale from wcs jacobian")
        sourceMapper.editOutputSchema().addField(
            "deltaMagBkg", type=np.float32, doc="Change in magnitude due to local background offset")

        return sourceMapper

    def fgcmMatchStars(self, visitCat, obsCat, lutDataRef=None):
        """
        Use FGCM code to match observations into unique stars.

        Parameters
        ----------
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for fgcm
        obsCat: `afw.table.BaseCatalog`
           Full catalog of star observations for fgcm
        lutDataRef: `lsst.daf.persistence.ButlerDataRef` or
                    `lsst.daf.butler.DeferredDatasetHandle`, optional
           Data reference to fgcm look-up table (used if matching reference stars).

        Returns
        -------
        fgcmStarIdCat: `afw.table.BaseCatalog`
           Catalog of unique star identifiers and index keys
        fgcmStarIndicesCat: `afwTable.BaseCatalog`
           Catalog of unique star indices
        fgcmRefCat: `afw.table.BaseCatalog`
           Catalog of matched reference stars.
           Will be None if `config.doReferenceMatches` is False.
        """
        # get filter names into a numpy array...
        # This is the type that is expected by the fgcm code
        visitFilterNames = np.zeros(len(visitCat), dtype='a30')
        for i in range(len(visitCat)):
            visitFilterNames[i] = visitCat[i]['physicalFilter']

        # match to put filterNames with observations
        visitIndex = np.searchsorted(visitCat['visit'],
                                     obsCat['visit'])

        obsFilterNames = visitFilterNames[visitIndex]

        if self.config.doReferenceMatches:
            # Get the reference filter names, using the LUT
            lutCat = lutDataRef.get()

            stdFilterDict = {filterName: stdFilter for (filterName, stdFilter) in
                             zip(lutCat[0]['physicalFilters'].split(','),
                                 lutCat[0]['stdPhysicalFilters'].split(','))}
            stdLambdaDict = {stdFilter: stdLambda for (stdFilter, stdLambda) in
                             zip(lutCat[0]['stdPhysicalFilters'].split(','),
                                 lutCat[0]['lambdaStdFilter'])}

            del lutCat

            referenceFilterNames = self._getReferenceFilterNames(visitCat,
                                                                 stdFilterDict,
                                                                 stdLambdaDict)
            self.log.info("Using the following reference filters: %s" %
                          (', '.join(referenceFilterNames)))

        else:
            # This should be an empty list
            referenceFilterNames = []

        # make the fgcm starConfig dict
        starConfig = {'logger': self.log,
                      'filterToBand': self.config.physicalFilterMap,
                      'requiredBands': self.config.requiredBands,
                      'minPerBand': self.config.minPerBand,
                      'matchRadius': self.config.matchRadius,
                      'isolationRadius': self.config.isolationRadius,
                      'matchNSide': self.config.matchNside,
                      'coarseNSide': self.config.coarseNside,
                      'densNSide': self.config.densityCutNside,
                      'densMaxPerPixel': self.config.densityCutMaxPerPixel,
                      'primaryBands': self.config.primaryBands,
                      'referenceFilterNames': referenceFilterNames}

        # initialize the FgcmMakeStars object
        fgcmMakeStars = fgcm.FgcmMakeStars(starConfig)

        # make the primary stars
        # note that the ra/dec native Angle format is radians
        # We determine the conversion from the native units (typically
        # radians) to degrees for the first observation.  This allows us
        # to treate ra/dec as numpy arrays rather than Angles, which would
        # be approximately 600x slower.
        conv = obsCat[0]['ra'].asDegrees() / float(obsCat[0]['ra'])
        fgcmMakeStars.makePrimaryStars(obsCat['ra'] * conv,
                                       obsCat['dec'] * conv,
                                       filterNameArray=obsFilterNames,
                                       bandSelected=False)

        # and match all the stars
        fgcmMakeStars.makeMatchedStars(obsCat['ra'] * conv,
                                       obsCat['dec'] * conv,
                                       obsFilterNames)

        if self.config.doReferenceMatches:
            fgcmMakeStars.makeReferenceMatches(self.fgcmLoadReferenceCatalog)

        # now persist

        objSchema = self._makeFgcmObjSchema()

        # make catalog and records
        fgcmStarIdCat = afwTable.BaseCatalog(objSchema)
        fgcmStarIdCat.reserve(fgcmMakeStars.objIndexCat.size)
        for i in range(fgcmMakeStars.objIndexCat.size):
            fgcmStarIdCat.addNew()

        # fill the catalog
        fgcmStarIdCat['fgcm_id'][:] = fgcmMakeStars.objIndexCat['fgcm_id']
        fgcmStarIdCat['ra'][:] = fgcmMakeStars.objIndexCat['ra']
        fgcmStarIdCat['dec'][:] = fgcmMakeStars.objIndexCat['dec']
        fgcmStarIdCat['obsArrIndex'][:] = fgcmMakeStars.objIndexCat['obsarrindex']
        fgcmStarIdCat['nObs'][:] = fgcmMakeStars.objIndexCat['nobs']

        obsSchema = self._makeFgcmObsSchema()

        fgcmStarIndicesCat = afwTable.BaseCatalog(obsSchema)
        fgcmStarIndicesCat.reserve(fgcmMakeStars.obsIndexCat.size)
        for i in range(fgcmMakeStars.obsIndexCat.size):
            fgcmStarIndicesCat.addNew()

        fgcmStarIndicesCat['obsIndex'][:] = fgcmMakeStars.obsIndexCat['obsindex']

        if self.config.doReferenceMatches:
            refSchema = self._makeFgcmRefSchema(len(referenceFilterNames))

            fgcmRefCat = afwTable.BaseCatalog(refSchema)
            fgcmRefCat.reserve(fgcmMakeStars.referenceCat.size)

            for i in range(fgcmMakeStars.referenceCat.size):
                fgcmRefCat.addNew()

            fgcmRefCat['fgcm_id'][:] = fgcmMakeStars.referenceCat['fgcm_id']
            fgcmRefCat['refMag'][:, :] = fgcmMakeStars.referenceCat['refMag']
            fgcmRefCat['refMagErr'][:, :] = fgcmMakeStars.referenceCat['refMagErr']

            md = PropertyList()
            md.set("REFSTARS_FORMAT_VERSION", REFSTARS_FORMAT_VERSION)
            md.set("FILTERNAMES", referenceFilterNames)
            fgcmRefCat.setMetadata(md)

        else:
            fgcmRefCat = None

        return fgcmStarIdCat, fgcmStarIndicesCat, fgcmRefCat

    def _makeFgcmVisitSchema(self, nCcd):
        """
        Make a schema for an fgcmVisitCatalog

        Parameters
        ----------
        nCcd: `int`
           Number of CCDs in the camera

        Returns
        -------
        schema: `afwTable.Schema`
        """

        schema = afwTable.Schema()
        schema.addField('visit', type=np.int32, doc="Visit number")
        schema.addField('physicalFilter', type=str, size=30, doc="Physical filter")
        schema.addField('telra', type=np.float64, doc="Pointing RA (deg)")
        schema.addField('teldec', type=np.float64, doc="Pointing Dec (deg)")
        schema.addField('telha', type=np.float64, doc="Pointing Hour Angle (deg)")
        schema.addField('telrot', type=np.float64, doc="Camera rotation (deg)")
        schema.addField('mjd', type=np.float64, doc="MJD of visit")
        schema.addField('exptime', type=np.float32, doc="Exposure time")
        schema.addField('pmb', type=np.float32, doc="Pressure (millibar)")
        schema.addField('psfSigma', type=np.float32, doc="PSF sigma (reference CCD)")
        schema.addField('deltaAper', type=np.float32, doc="Delta-aperture")
        schema.addField('skyBackground', type=np.float32, doc="Sky background (ADU) (reference CCD)")
        # the following field is not used yet
        schema.addField('deepFlag', type=np.int32, doc="Deep observation")
        schema.addField('scaling', type='ArrayD', doc="Scaling applied due to flat adjustment",
                        size=nCcd)
        schema.addField('used', type=np.int32, doc="This visit has been ingested.")
        schema.addField('sources_read', type='Flag', doc="This visit had sources read.")

        return schema

    def _makeFgcmObjSchema(self):
        """
        Make a schema for the objIndexCat from fgcmMakeStars

        Returns
        -------
        schema: `afwTable.Schema`
        """

        objSchema = afwTable.Schema()
        objSchema.addField('fgcm_id', type=np.int32, doc='FGCM Unique ID')
        # Will investigate making these angles...
        objSchema.addField('ra', type=np.float64, doc='Mean object RA (deg)')
        objSchema.addField('dec', type=np.float64, doc='Mean object Dec (deg)')
        objSchema.addField('obsArrIndex', type=np.int32,
                           doc='Index in obsIndexTable for first observation')
        objSchema.addField('nObs', type=np.int32, doc='Total number of observations')

        return objSchema

    def _makeFgcmObsSchema(self):
        """
        Make a schema for the obsIndexCat from fgcmMakeStars

        Returns
        -------
        schema: `afwTable.Schema`
        """

        obsSchema = afwTable.Schema()
        obsSchema.addField('obsIndex', type=np.int32, doc='Index in observation table')

        return obsSchema

    def _makeFgcmRefSchema(self, nReferenceBands):
        """
        Make a schema for the referenceCat from fgcmMakeStars

        Parameters
        ----------
        nReferenceBands: `int`
           Number of reference bands

        Returns
        -------
        schema: `afwTable.Schema`
        """

        refSchema = afwTable.Schema()
        refSchema.addField('fgcm_id', type=np.int32, doc='FGCM Unique ID')
        refSchema.addField('refMag', type='ArrayF', doc='Reference magnitude array (AB)',
                           size=nReferenceBands)
        refSchema.addField('refMagErr', type='ArrayF', doc='Reference magnitude error array',
                           size=nReferenceBands)

        return refSchema

    def _getReferenceFilterNames(self, visitCat, stdFilterDict, stdLambdaDict):
        """
        Get the reference filter names, in wavelength order, from the visitCat and
        information from the look-up-table.

        Parameters
        ----------
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for FGCM
        stdFilterDict: `dict`
           Mapping of filterName to stdFilterName from LUT
        stdLambdaDict: `dict`
           Mapping of stdFilterName to stdLambda from LUT

        Returns
        -------
        referenceFilterNames: `list`
           Wavelength-ordered list of reference filter names
        """

        # Find the unique list of filter names in visitCat
        filterNames = np.unique(visitCat.asAstropy()['physicalFilter'])

        # Find the unique list of "standard" filters
        stdFilterNames = {stdFilterDict[filterName] for filterName in filterNames}

        # And sort these by wavelength
        referenceFilterNames = sorted(stdFilterNames, key=stdLambdaDict.get)

        return referenceFilterNames
