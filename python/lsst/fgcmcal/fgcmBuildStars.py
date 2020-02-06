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
"""Build star observations for input to FGCM.

This task finds all the visits and calexps in a repository (or a subset
based on command line parameters) and extract all the potential calibration
stars for input into fgcm.  This task additionally uses fgcm to match
star observations into unique stars, and performs as much cleaning of
the input catalog as possible.
"""

import sys
import time
import traceback

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.geom as geom
from lsst.daf.base import PropertyList
from lsst.daf.base.dateTime import DateTime
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .fgcmLoadReferenceCatalog import FgcmLoadReferenceCatalogTask
from .utilities import computeApproxPixelAreaFields, computeApertureRadius

import fgcm

REFSTARS_FORMAT_VERSION = 1

__all__ = ['FgcmBuildStarsConfig', 'FgcmBuildStarsTask', 'FgcmBuildStarsRunner']


class FgcmBuildStarsConfig(pexConfig.Config):
    """Config for FgcmBuildStarsTask"""

    instFluxField = pexConfig.Field(
        doc=("Name of the source instFlux field to use.  The associated flag field "
             "('<name>_flag') will be implicitly included in badFlags"),
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
    referenceCCD = pexConfig.Field(
        doc="Reference CCD for scanning visits",
        dtype=int,
        default=13,
    )
    checkAllCcds = pexConfig.Field(
        doc=("Check repo for all CCDs for each visit specified.  To be used when the "
             "full set of ids (visit/ccd) are not specified on the command line.  For "
             "Gen2, specifying one ccd and setting checkAllCcds=True is significantly "
             "faster than the alternatives."),
        dtype=bool,
        default=True,
    )
    visitDataRefName = pexConfig.Field(
        doc="dataRef name for the 'visit' field",
        dtype=str,
        default="visit"
    )
    ccdDataRefName = pexConfig.Field(
        doc="dataRef name for the 'ccd' field",
        dtype=str,
        default="ccd"
    )
    applyJacobian = pexConfig.Field(
        doc="Apply Jacobian correction?",
        dtype=bool,
        deprecated=("This field is no longer used, and has been deprecated by DM-20163. "
                    "It will be removed after v20."),
        default=False
    )
    jacobianName = pexConfig.Field(
        doc="Name of field with jacobian correction",
        dtype=str,
        deprecated=("This field is no longer used, and has been deprecated by DM-20163. "
                    "It will be removed after v20."),
        default="base_Jacobian_value"
    )
    doApplyWcsJacobian = pexConfig.Field(
        doc="Apply the jacobian of the WCS to the star observations prior to fit?",
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
        doc="Name of the local background instFlux field to use.",
        dtype=str,
        default='base_LocalBackground_instFlux'
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources",
        default="science"
    )
    apertureInnerInstFluxField = pexConfig.Field(
        doc="Field that contains inner aperture for aperture correction proxy",
        dtype=str,
        default='base_CircularApertureFlux_12_0_instFlux'
    )
    apertureOuterInstFluxField = pexConfig.Field(
        doc="Field that contains outer aperture for aperture correction proxy",
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

    def setDefaults(self):
        sourceSelector = self.sourceSelector["science"]
        sourceSelector.setDefaults()

        fluxFlagName = self.instFluxField[0: -len('instFlux')] + 'flag'

        sourceSelector.flags.bad = ['base_PixelFlags_flag_edge',
                                    'base_PixelFlags_flag_interpolatedCenter',
                                    'base_PixelFlags_flag_saturatedCenter',
                                    'base_PixelFlags_flag_crCenter',
                                    'base_PixelFlags_flag_bad',
                                    'base_PixelFlags_flag_interpolated',
                                    'base_PixelFlags_flag_saturated',
                                    'slot_Centroid_flag',
                                    fluxFlagName]

        if self.doSubtractLocalBackground:
            localBackgroundFlagName = self.localBackgroundFluxField[0: -len('instFlux')] + 'flag'
            sourceSelector.flags.bad.append(localBackgroundFlagName)

        sourceSelector.doFlags = True
        sourceSelector.doUnresolved = True
        sourceSelector.doSignalToNoise = True
        sourceSelector.doIsolated = True

        sourceSelector.signalToNoise.fluxField = self.instFluxField
        sourceSelector.signalToNoise.errField = self.instFluxField + 'Err'
        sourceSelector.signalToNoise.minimum = 10.0
        sourceSelector.signalToNoise.maximum = 1000.0

        # FGCM operates on unresolved sources, and this setting is
        # appropriate for the current base_ClassificationExtendedness
        sourceSelector.unresolved.maximum = 0.5


class FgcmBuildStarsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmBuildStarsTask

    fgcmBuildStarsTask.run() takes a number of arguments, one of which is the
    butler (for persistence and mapper data), and a list of dataRefs
    extracted from the command line.  Note that FGCM runs on a large set of
    dataRefs, and not on single dataRef/tract/patch.
    This class transforms the process arguments generated by the ArgumentParser
    into the arguments expected by FgcmBuildStarsTask.run().
    This runner does not use any parallelization.

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


class FgcmBuildStarsTask(pipeBase.CmdLineTask):
    """
    Build stars for the FGCM global calibration
    """

    ConfigClass = FgcmBuildStarsConfig
    RunnerClass = FgcmBuildStarsRunner
    _DefaultName = "fgcmBuildStars"

    def __init__(self, butler=None, **kwargs):
        """
        Instantiate an `FgcmBuildStarsTask`.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.makeSubtask("sourceSelector")
        # Only log warning and fatal errors from the sourceSelector
        self.sourceSelector.log.setLevel(self.sourceSelector.log.WARN)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="Data ID, e.g. --id visit=6789")

        return parser

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
           Data references for the input visits.
           If this is an empty list, all visits with src catalogs in
           the repository are used.
           Only one individual dataRef from a visit need be specified
           and the code will find the other source catalogs from
           each visit.

        Raises
        ------
        RuntimeErrror: Raised if `config.doReferenceMatches` is set and
           an fgcmLookUpTable is not available, or if computeFluxApertureRadius()
           fails if the calibFlux is not a CircularAperture flux.
        """

        if self.config.doReferenceMatches:
            # Ensure that we have a LUT
            if not butler.datasetExists('fgcmLookUpTable'):
                raise RuntimeError("Must have fgcmLookUpTable if using config.doReferenceMatches")
        # Compute aperture radius if necessary.  This is useful to do now before
        # any heavy lifting has happened (fail early).
        calibFluxApertureRadius = None
        if self.config.doSubtractLocalBackground:
            sourceSchema = butler.get('src_schema').schema
            try:
                calibFluxApertureRadius = computeApertureRadius(sourceSchema,
                                                                self.config.instFluxField)
            except (RuntimeError, LookupError):
                raise RuntimeError("Could not determine aperture radius from %s. "
                                   "Cannot use doSubtractLocalBackground." %
                                   (self.config.instFluxField))

        groupedDataRefs = self.findAndGroupDataRefs(butler, dataRefs)

        camera = butler.get('camera')

        # Make the visit catalog if necessary
        if not butler.datasetExists('fgcmVisitCatalog'):
            # we need to build visitCat
            visitCat = self.fgcmMakeVisitCatalog(camera, groupedDataRefs)
        else:
            self.log.info("Found fgcmVisitCatalog.")
            visitCat = butler.get('fgcmVisitCatalog')

        # Compile all the stars
        if not butler.datasetExists('fgcmStarObservations'):
            rad = calibFluxApertureRadius
            fgcmStarObservationCat = self.fgcmMakeAllStarObservations(groupedDataRefs,
                                                                      visitCat,
                                                                      calibFluxApertureRadius=rad)
        else:
            self.log.info("Found fgcmStarObservations")
            fgcmStarObservationCat = butler.get('fgcmStarObservations')

        if not butler.datasetExists('fgcmStarIds') or not butler.datasetExists('fgcmStarIndices'):
            fgcmStarIdCat, fgcmStarIndicesCat, fgcmRefCat = self.fgcmMatchStars(butler,
                                                                                visitCat,
                                                                                fgcmStarObservationCat)
        else:
            self.log.info("Found fgcmStarIds and fgcmStarIndices")

        # Persist catalogs via the butler
        butler.put(visitCat, 'fgcmVisitCatalog')
        butler.put(fgcmStarObservationCat, 'fgcmStarObservations')
        butler.put(fgcmStarIdCat, 'fgcmStarIds')
        butler.put(fgcmStarIndicesCat, 'fgcmStarIndices')
        if fgcmRefCat is not None:
            butler.put(fgcmRefCat, 'fgcmReferenceStars')

    def fgcmMakeVisitCatalog(self, camera, groupedDataRefs):
        """
        Make a visit catalog with all the keys from each visit

        Parameters
        ----------
        camera: `lsst.afw.cameraGeom.Camera`
           Camera from the butler
        groupedDataRefs: `dict`
           Dictionary with visit keys, and `list`s of
           `lsst.daf.persistence.ButlerDataRef`

        Returns
        -------
        visitCat: `afw.table.BaseCatalog`
        """

        nCcd = len(camera)

        schema = self._makeFgcmVisitSchema(nCcd)

        visitCat = afwTable.BaseCatalog(schema)
        visitCat.reserve(len(groupedDataRefs))

        self._fillVisitCatalog(visitCat, groupedDataRefs)

        return visitCat

    def _fillVisitCatalog(self, visitCat, groupedDataRefs):
        """
        Fill the visit catalog with visit metadata

        Parameters
        ----------
        visitCat: `afw.table.BaseCatalog`
           Catalog with schema from _makeFgcmVisitSchema()
        groupedDataRefs: `dict`
           Dictionary with visit keys, and `list`s of
           `lsst.daf.persistence.ButlerDataRef`
        """

        bbox = geom.BoxI(geom.PointI(0, 0), geom.PointI(1, 1))

        for i, visit in enumerate(sorted(groupedDataRefs)):
            # We don't use the bypasses since we need the psf info which does
            # not have a bypass
            # TODO: When DM-15500 is implemented in the Gen3 Butler, this
            # can be fixed

            # Note that the reference ccd is first in the list (if available).

            # The first dataRef in the group will be the reference ccd (if available)
            dataRef = groupedDataRefs[visit][0]

            exp = dataRef.get(datasetType='calexp_sub', bbox=bbox,
                              flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

            visitInfo = exp.getInfo().getVisitInfo()
            f = exp.getFilter()
            psf = exp.getPsf()

            rec = visitCat.addNew()
            rec['visit'] = visit
            rec['filtername'] = f.getName()
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

            rec['psfSigma'] = psf.computeShape().getDeterminantRadius()

            if dataRef.datasetExists(datasetType='calexpBackground'):
                # Get background for reference CCD
                # This approximation is good enough for now
                bgStats = (bg[0].getStatsImage().getImage().array
                           for bg in dataRef.get(datasetType='calexpBackground'))
                rec['skyBackground'] = sum(np.median(bg[np.isfinite(bg)]) for bg in bgStats)
            else:
                self.log.warn('Sky background not found for visit %d / ccd %d' %
                              (visit, dataRef.dataId[self.config.ccdDataRefName]))
                rec['skyBackground'] = -1.0

    def findAndGroupDataRefs(self, butler, dataRefs):
        """
        Find and group dataRefs (by visit).  If dataRefs is an empty list,
        this will look for all source catalogs in a given repo.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        dataRefs: `list` of `lsst.daf.persistence.ButlerDataRef`
           Data references for the input visits.
           If this is an empty list, all visits with src catalogs in
           the repository are used.

        Returns
        -------
        groupedDataRefs: `dict`
           Dictionary with visit keys, and `list`s of `lsst.daf.persistence.ButlerDataRef`
        """

        camera = butler.get('camera')

        ccdIds = []
        for detector in camera:
            ccdIds.append(detector.getId())

        # TODO: related to DM-13730, this dance of looking for source visits
        # will be unnecessary with Gen3 Butler.  This should be part of
        # DM-13730.

        groupedDataRefs = {}
        for dataRef in dataRefs:
            visit = dataRef.dataId[self.config.visitDataRefName]
            # If we don't have the dataset, just continue
            if not dataRef.datasetExists(datasetType='src'):
                continue
            # If we need to check all ccds, do it here
            if self.config.checkAllCcds:
                dataId = dataRef.dataId.copy()
                # For each ccd we must check that a valid source catalog exists.
                for ccdId in ccdIds:
                    dataId[self.config.ccdDataRefName] = ccdId
                    if butler.datasetExists('src', dataId=dataId):
                        goodDataRef = butler.dataRef('src', dataId=dataId)
                        if visit in groupedDataRefs:
                            if (goodDataRef.dataId[self.config.ccdDataRefName] not in
                               [d.dataId[self.config.ccdDataRefName] for d in groupedDataRefs[visit]]):
                                groupedDataRefs[visit].append(goodDataRef)
                        else:
                            groupedDataRefs[visit] = [goodDataRef]
            else:
                # We have already confirmed that the dataset exists, so no need
                # to check here.
                if visit in groupedDataRefs:
                    if (dataRef.dataId[self.config.ccdDataRefName] not in
                       [d.dataId[self.config.ccdDataRefName] for d in groupedDataRefs[visit]]):
                        groupedDataRefs[visit].append(dataRef)
                else:
                    groupedDataRefs[visit] = [dataRef]

        # Put them in ccd order, with the reference ccd first (if available)
        def ccdSorter(dataRef):
            ccdId = dataRef.dataId[self.config.ccdDataRefName]
            if ccdId == self.config.referenceCCD:
                return -100
            else:
                return ccdId

        # If we did not check all ccds, put them in ccd order
        if not self.config.checkAllCcds:
            for visit in groupedDataRefs:
                groupedDataRefs[visit] = sorted(groupedDataRefs[visit], key=ccdSorter)

        return groupedDataRefs

    def fgcmMakeAllStarObservations(self, groupedDataRefs, visitCat,
                                    calibFluxApertureRadius=None):
        """
        Compile all good star observations from visits in visitCat.

        Parameters
        ----------
        groupedDataRefs: `dict` of `list`s
           Lists of `lsst.daf.persistence.ButlerDataRef`, grouped by visit.
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for FGCM
        calibFluxApertureRadius: `float`, optional
           Aperture radius for calibration flux.  Default is None.

        Returns
        -------
        fgcmStarObservations: `afw.table.BaseCatalog`
           Full catalog of good observations.

        Raises
        ------
        RuntimeError: Raised if doSubtractLocalBackground is True and
           calibFluxApertureRadius is not set.
        """
        startTime = time.time()

        if self.config.doSubtractLocalBackground and calibFluxApertureRadius is None:
            raise RuntimeError("Must set calibFluxApertureRadius if doSubtractLocalBackground is True.")

        # create our source schema.  Use the first valid dataRef
        dataRef = groupedDataRefs[list(groupedDataRefs.keys())[0]][0]
        sourceSchema = dataRef.get('src_schema', immediate=True).schema

        # Construct a mapping from ccd number to index
        camera = dataRef.get('camera')
        ccdMapping = {}
        for ccdIndex, detector in enumerate(camera):
            ccdMapping[detector.getId()] = ccdIndex

        approxPixelAreaFields = computeApproxPixelAreaFields(camera)

        sourceMapper = self._makeSourceMapper(sourceSchema)

        # We also have a temporary catalog that will accumulate aperture measurements
        aperMapper = self._makeAperMapper(sourceSchema)

        outputSchema = sourceMapper.getOutputSchema()
        fullCatalog = afwTable.BaseCatalog(outputSchema)

        # FGCM will provide relative calibration for the flux in config.instFluxField

        instFluxKey = sourceSchema[self.config.instFluxField].asKey()
        instFluxErrKey = sourceSchema[self.config.instFluxField + 'Err'].asKey()
        visitKey = outputSchema['visit'].asKey()
        ccdKey = outputSchema['ccd'].asKey()
        instMagKey = outputSchema['instMag'].asKey()
        instMagErrKey = outputSchema['instMagErr'].asKey()

        # Prepare local background if desired
        if self.config.doSubtractLocalBackground:
            localBackgroundFluxKey = sourceSchema[self.config.localBackgroundFluxField].asKey()
            localBackgroundArea = np.pi*calibFluxApertureRadius**2.
        else:
            localBackground = 0.0

        aperOutputSchema = aperMapper.getOutputSchema()

        instFluxAperInKey = sourceSchema[self.config.apertureInnerInstFluxField].asKey()
        instFluxErrAperInKey = sourceSchema[self.config.apertureInnerInstFluxField + 'Err'].asKey()
        instFluxAperOutKey = sourceSchema[self.config.apertureOuterInstFluxField].asKey()
        instFluxErrAperOutKey = sourceSchema[self.config.apertureOuterInstFluxField + 'Err'].asKey()
        instMagInKey = aperOutputSchema['instMag_aper_inner'].asKey()
        instMagErrInKey = aperOutputSchema['instMagErr_aper_inner'].asKey()
        instMagOutKey = aperOutputSchema['instMag_aper_outer'].asKey()
        instMagErrOutKey = aperOutputSchema['instMagErr_aper_outer'].asKey()

        k = 2.5 / np.log(10.)

        # loop over visits
        for visit in visitCat:
            expTime = visit['exptime']

            nStarInVisit = 0

            # Reset the aperture catalog (per visit)
            aperVisitCatalog = afwTable.BaseCatalog(aperOutputSchema)

            for dataRef in groupedDataRefs[visit['visit']]:

                ccdId = dataRef.dataId[self.config.ccdDataRefName]

                sources = dataRef.get(datasetType='src', flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

                # If we are subtracting the local background, then correct here
                # before we do the s/n selection.  This ensures we do not have
                # bad stars after local background subtraction.

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

                    localBackground = localBackgroundArea*sources[localBackgroundFluxKey]
                    sources[instFluxKey] -= localBackground

                goodSrc = self.sourceSelector.selectSources(sources)

                tempCat = afwTable.BaseCatalog(fullCatalog.schema)
                tempCat.reserve(goodSrc.selected.sum())
                tempCat.extend(sources[goodSrc.selected], mapper=sourceMapper)
                tempCat[visitKey][:] = visit['visit']
                tempCat[ccdKey][:] = ccdId

                # Compute "instrumental magnitude" by scaling flux with exposure time.
                scaledInstFlux = (sources[instFluxKey][goodSrc.selected] *
                                  visit['scaling'][ccdMapping[ccdId]])
                tempCat[instMagKey][:] = (-2.5*np.log10(scaledInstFlux) + 2.5*np.log10(expTime))

                # Compute instMagErr from instFluxErr / instFlux, any scaling
                # will cancel out.

                tempCat[instMagErrKey][:] = k*(sources[instFluxErrKey][goodSrc.selected] /
                                               sources[instFluxKey][goodSrc.selected])

                # Compute the jacobian from an approximate PixelAreaBoundedField
                tempCat['jacobian'] = approxPixelAreaFields[ccdId].evaluate(tempCat['x'],
                                                                            tempCat['y'])

                # Apply the jacobian if configured
                if self.config.doApplyWcsJacobian:
                    tempCat[instMagKey][:] -= 2.5*np.log10(tempCat['jacobian'][:])

                fullCatalog.extend(tempCat)

                # And the aperture information
                # This does not need the jacobian because it is all locally relative
                tempAperCat = afwTable.BaseCatalog(aperVisitCatalog.schema)
                tempAperCat.reserve(goodSrc.selected.sum())
                tempAperCat.extend(sources[goodSrc.selected], mapper=aperMapper)

                with np.warnings.catch_warnings():
                    # Ignore warnings, we will filter infinities and
                    # nans below.
                    np.warnings.simplefilter("ignore")

                    tempAperCat[instMagInKey][:] = -2.5*np.log10(
                        sources[instFluxAperInKey][goodSrc.selected])
                    tempAperCat[instMagErrInKey][:] = (2.5/np.log(10.))*(
                        sources[instFluxErrAperInKey][goodSrc.selected] /
                        sources[instFluxAperInKey][goodSrc.selected])
                    tempAperCat[instMagOutKey][:] = -2.5*np.log10(
                        sources[instFluxAperOutKey][goodSrc.selected])
                    tempAperCat[instMagErrOutKey][:] = (2.5/np.log(10.))*(
                        sources[instFluxErrAperOutKey][goodSrc.selected] /
                        sources[instFluxAperOutKey][goodSrc.selected])

                aperVisitCatalog.extend(tempAperCat)

                nStarInVisit += len(tempCat)

            # Compute the median delta-aper
            if not aperVisitCatalog.isContiguous():
                aperVisitCatalog = aperVisitCatalog.copy(deep=True)

            instMagIn = aperVisitCatalog[instMagInKey]
            instMagErrIn = aperVisitCatalog[instMagErrInKey]
            instMagOut = aperVisitCatalog[instMagOutKey]
            instMagErrOut = aperVisitCatalog[instMagErrOutKey]

            ok = (np.isfinite(instMagIn) & np.isfinite(instMagErrIn) &
                  np.isfinite(instMagOut) & np.isfinite(instMagErrOut))

            visit['deltaAper'] = np.median(instMagIn[ok] - instMagOut[ok])

            self.log.info("  Found %d good stars in visit %d (deltaAper = %.3f)" %
                          (nStarInVisit, visit['visit'], visit['deltaAper']))

        self.log.info("Found all good star observations in %.2f s" %
                      (time.time() - startTime))

        return fullCatalog

    def fgcmMatchStars(self, butler, visitCat, obsCat):
        """
        Use FGCM code to match observations into unique stars.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for fgcm
        obsCat: `afw.table.BaseCatalog`
           Full catalog of star observations for fgcm

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

        if self.config.doReferenceMatches:
            # Make a subtask for reference loading
            self.makeSubtask("fgcmLoadReferenceCatalog", butler=butler)

        # get filter names into a numpy array...
        # This is the type that is expected by the fgcm code
        visitFilterNames = np.zeros(len(visitCat), dtype='a10')
        for i in range(len(visitCat)):
            visitFilterNames[i] = visitCat[i]['filtername']

        # match to put filterNames with observations
        visitIndex = np.searchsorted(visitCat['visit'],
                                     obsCat['visit'])

        obsFilterNames = visitFilterNames[visitIndex]

        if self.config.doReferenceMatches:
            # Get the reference filter names, using the LUT
            lutCat = butler.get('fgcmLookUpTable')

            stdFilterDict = {filterName: stdFilter for (filterName, stdFilter) in
                             zip(lutCat[0]['filterNames'].split(','),
                                 lutCat[0]['stdFilterNames'].split(','))}
            stdLambdaDict = {stdFilter: stdLambda for (stdFilter, stdLambda) in
                             zip(lutCat[0]['stdFilterNames'].split(','),
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
                      'filterToBand': self.config.filterMap,
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
        # Note that the FGCM code currently handles filternames up to 2 characters long
        schema.addField('filtername', type=str, size=10, doc="Filter name")
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

        return schema

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
        sourceMapper.addMapping(sourceSchema[self.config.psfCandidateName].asKey(),
                                'psf_candidate')

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

        return sourceMapper

    def _makeAperMapper(self, sourceSchema):
        """
        Make a schema mapper for fgcm aperture measurements

        Parameters
        ----------
        sourceSchema: `afwTable.Schema`
           Default source schema from the butler

        Returns
        -------
        aperMapper: `afwTable.schemaMapper`
           Mapper to the FGCM aperture schema
        """

        aperMapper = afwTable.SchemaMapper(sourceSchema)
        aperMapper.addMapping(sourceSchema['coord_ra'].asKey(), 'ra')
        aperMapper.addMapping(sourceSchema['coord_dec'].asKey(), 'dec')
        aperMapper.editOutputSchema().addField('instMag_aper_inner', type=np.float64,
                                               doc="Magnitude at inner aperture")
        aperMapper.editOutputSchema().addField('instMagErr_aper_inner', type=np.float64,
                                               doc="Magnitude error at inner aperture")
        aperMapper.editOutputSchema().addField('instMag_aper_outer', type=np.float64,
                                               doc="Magnitude at outer aperture")
        aperMapper.editOutputSchema().addField('instMagErr_aper_outer', type=np.float64,
                                               doc="Magnitude error at outer aperture")

        return aperMapper

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
        filterNames = np.unique(visitCat.asAstropy()['filtername'])

        # Find the unique list of "standard" filters
        stdFilterNames = {stdFilterDict[filterName] for filterName in filterNames}

        # And sort these by wavelength
        referenceFilterNames = sorted(stdFilterNames, key=stdLambdaDict.get)

        return referenceFilterNames
