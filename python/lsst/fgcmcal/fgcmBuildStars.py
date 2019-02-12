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
import lsst.afw.geom as afwGeom
from lsst.daf.base.dateTime import DateTime
import lsst.daf.persistence.butlerExceptions as butlerExceptions
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .fgcmLoadReferenceCatalog import FgcmLoadReferenceCatalogTask

import fgcm

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
    filterToBand = pexConfig.DictField(
        doc="filterName to band mapping",
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
        doc="Check all CCDs.  Necessary for testing",
        dtype=bool,
        default=False,
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
        default=False
    )
    jacobianName = pexConfig.Field(
        doc="Name of field with jacobian correction",
        dtype=str,
        default="base_Jacobian_value"
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
        default=False,
    )
    fgcmLoadReferenceCatalog = pexConfig.ConfigurableField(
        target=FgcmLoadReferenceCatalogTask,
        doc="FGCM reference object loader",
    )
    referenceBands = pexConfig.ListField(
        doc="bands, in wavelength order, that will be calibrated",
        dtype=str,
        default=(),
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

        sourceSelector.doFlags = True
        sourceSelector.doUnresolved = True
        sourceSelector.doSignalToNoise = True
        sourceSelector.doIsolated = True

        sourceSelector.signalToNoise.fluxField = self.instFluxField
        sourceSelector.signalToNoise.errField = self.instFluxField + 'Err'
        sourceSelector.signalToNoise.minimum = 10.0
        sourceSelector.signalToNoise.maximum = 1000.0

        # FGCM operates on unresolved sources, and this setting is
        # appropriate for the current base_classificationExtendedness
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
        exitStatus: `list` with `pipeBase.Struct`
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
        parsedCmd: ArgumentParser parsed command line
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
        Instantiate an FgcmBuildStarsTask.

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
        parser.add_id_argument("--id", "calexp", help="Data ID, e.g. --id visit=6789 (optional)")

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
        """

        if self.config.doReferenceMatches:
            self.makeSubtask("fgcmLoadReferenceCatalog", butler=butler)

        # Make the visit catalog if necessary
        if not butler.datasetExists('fgcmVisitCatalog'):
            # we need to build visitCat
            visitCat = self._fgcmMakeVisitCatalog(butler, dataRefs)
        else:
            self.log.info("Found fgcmVisitCatalog.")
            visitCat = butler.get('fgcmVisitCatalog')

        # Compile all the stars
        if not butler.datasetExists('fgcmStarObservations'):
            self._fgcmMakeAllStarObservations(butler, visitCat)
        else:
            self.log.info("Found fgcmStarObservations")

        if not butler.datasetExists('fgcmStarIds') or not butler.datasetExists('fgcmStarIndices'):
            self._fgcmMatchStars(butler, visitCat)
        else:
            self.log.info("Found fgcmStarIds and fgcmStarIndices")

    def _fgcmMakeVisitCatalog(self, butler, dataRefs):
        """
        Make a visit catalog with all the key data from each visit

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        dataRefs: `list` of `lsst.daf.persistence.ButlerDataRef`
           Data references for the input visits.
           If this is an empty list, all visits with src catalogs in
           the repository are used.
           Only one individual dataRef from a visit need be specified
           and the code will find the other source catalogs from
           each visit.

        Returns
        -------
        visitCat: `afw.table.BaseCatalog`
        """

        startTime = time.time()

        camera = butler.get('camera')
        nCcd = len(camera)

        # TODO: related to DM-13730, this dance of looking for source visits
        # will be unnecessary with Gen3 Butler.  This should be part of
        # DM-13730.

        if len(dataRefs) == 0:
            srcVisits, srcCcds = self._findSourceVisits(butler, nCcd)
        else:
            # get the visits from the datarefs, only for referenceCCD
            srcVisits = [d.dataId[self.config.visitDataRefName] for d in dataRefs if
                         d.dataId[self.config.ccdDataRefName] == self.config.referenceCCD]
            srcCcds = [self.config.referenceCCD] * len(srcVisits)

        # Sort the visits for searching/indexing
        srcVisits.sort()

        self.log.info("Found %d visits in %.2f s" %
                      (len(srcVisits), time.time() - startTime))

        schema = self._makeFgcmVisitSchema(nCcd)

        visitCat = afwTable.BaseCatalog(schema)
        visitCat.table.preallocate(len(srcVisits))

        startTime = time.time()

        self._fillVisitCatalog(butler, visitCat, srcVisits, srcCcds)

        self.log.info("Found all VisitInfo in %.2f s" % (time.time() - startTime))

        return visitCat

    def _findSourceVisits(self, butler, nCcd):
        """
        Find all source catalogs in the repository

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        nCcd: `int`
           Number of CCDs in the camera

        Returns
        -------
        (srcVisits, srcCcds): `tuple` of `list`s
        """

        # We did not specify any datarefs, so find all of them
        if not self.config.checkAllCcds:
            # Faster mode, scan through referenceCCD
            allVisits = butler.queryMetadata('src',
                                             format=[self.config.visitDataRefName, 'filter'],
                                             dataId={self.config.ccdDataRefName:
                                                     self.config.referenceCCD})
            srcVisits = []
            srcCcds = []
            for dataset in allVisits:
                if (butler.datasetExists('src', dataId={self.config.visitDataRefName: dataset[0],
                                                        self.config.ccdDataRefName:
                                                        self.config.referenceCCD})):
                    srcVisits.append(dataset[0])
                    srcCcds.append(self.config.referenceCCD)
        else:
            # Slower mode, check all CCDs
            allVisits = butler.queryMetadata('src',
                                             format=[self.config.visitDataRefName, 'filter'])
            srcVisits = []
            srcCcds = []

            for dataset in allVisits:
                if dataset[0] in srcVisits:
                    continue
                for ccd in range(nCcd):
                    if (butler.datasetExists('src', dataId={self.config.visitDataRefName: dataset[0],
                                                            self.config.ccdDataRefName:
                                                                ccd})):
                        srcVisits.append(dataset[0])
                        srcCcds.append(ccd)
                        # Once we find that a butler dataset exists, break out
                        break

        return (srcVisits, srcCcds)

    def _fillVisitCatalog(self, butler, visitCat, srcVisits, srcCcds):
        """
        Fill the visit catalog with visit metadata

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        visitCat: `afw.table.BaseCatalog`
           Catalog with schema from _createFgcmVisitSchema()
        srcVisits: `list'
           List of source visits
        srcCcds: `list`
           List of source CCDs
        """

        bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.PointI(1, 1))

        # now loop over visits and get the information
        for i, srcVisit in enumerate(srcVisits):
            # We don't use the bypasses since we need the psf info which does
            # not have a bypass
            # TODO: When DM-15500 is implemented in the Gen3 Butler, this
            # can be fixed

            dataId = {self.config.visitDataRefName: srcVisit,
                      self.config.ccdDataRefName: srcCcds[i]}

            exp = butler.get('calexp_sub', dataId=dataId, bbox=bbox,
                             flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
            visitInfo = exp.getInfo().getVisitInfo()
            f = exp.getFilter()
            psf = exp.getPsf()

            rec = visitCat.addNew()
            rec['visit'] = srcVisit
            rec['filtername'] = f.getName()
            radec = visitInfo.getBoresightRaDec()
            rec['telra'] = radec.getRa().asDegrees()
            rec['teldec'] = radec.getDec().asDegrees()
            rec['telha'] = visitInfo.getBoresightHourAngle().asDegrees()
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

            if butler.datasetExists('calexpBackground', dataId=dataId):
                # Get background for reference CCD
                # This approximation is good enough for now
                bgStats = (bg[0].getStatsImage().getImage().array
                           for bg in butler.get('calexpBackground',
                                                dataId=dataId))
                rec['skyBackground'] = sum(np.median(bg[np.isfinite(bg)]) for bg in bgStats)
            else:
                rec['skyBackground'] = -1.0

    def _fgcmMakeAllStarObservations(self, butler, visitCat):
        """
        Compile all good star observations from visits in visitCat

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for FGCM
        """

        startTime = time.time()

        # create our source schema
        sourceSchema = butler.get('src_schema', immediate=True).schema

        sourceMapper = self._makeSourceMapper(sourceSchema)

        # We also have a temporary catalog that will accumulate aperture measurements
        aperMapper = self._makeAperMapper(sourceSchema)

        # we need to know the ccds...
        camera = butler.get('camera')

        outputSchema = sourceMapper.getOutputSchema()
        fullCatalog = afwTable.BaseCatalog(outputSchema)

        # FGCM will provide relative calibration for the flux in config.instFluxField

        instFluxKey = sourceSchema[self.config.instFluxField].asKey()
        instFluxErrKey = sourceSchema[self.config.instFluxField + 'Err'].asKey()
        visitKey = outputSchema['visit'].asKey()
        ccdKey = outputSchema['ccd'].asKey()
        instMagKey = outputSchema['instMag'].asKey()
        instMagErrKey = outputSchema['instMagErr'].asKey()

        aperOutputSchema = aperMapper.getOutputSchema()

        instFluxAperInKey = sourceSchema[self.config.apertureInnerInstFluxField].asKey()
        instFluxErrAperInKey = sourceSchema[self.config.apertureInnerInstFluxField + 'Err'].asKey()
        instFluxAperOutKey = sourceSchema[self.config.apertureOuterInstFluxField].asKey()
        instFluxErrAperOutKey = sourceSchema[self.config.apertureOuterInstFluxField + 'Err'].asKey()
        instMagInKey = aperOutputSchema['instMag_aper_inner'].asKey()
        instMagErrInKey = aperOutputSchema['instMagErr_aper_inner'].asKey()
        instMagOutKey = aperOutputSchema['instMag_aper_outer'].asKey()
        instMagErrOutKey = aperOutputSchema['instMagErr_aper_outer'].asKey()

        # loop over visits
        for visit in visitCat:
            expTime = visit['exptime']

            nStarInVisit = 0

            # Reset the aperture catalog (per visit)
            aperVisitCatalog = afwTable.BaseCatalog(aperOutputSchema)

            # loop over CCDs
            for ccdIndex, detector in enumerate(camera):

                ccdId = detector.getId()

                try:
                    sources = butler.get('src', dataId={self.config.visitDataRefName:
                                                        visit['visit'],
                                                        self.config.ccdDataRefName: ccdId},
                                         flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                except butlerExceptions.NoResults:
                    # this is not a problem if this ccd isn't there
                    continue

                goodSrc = self.sourceSelector.selectSources(sources)

                tempCat = afwTable.BaseCatalog(fullCatalog.schema)
                tempCat.reserve(goodSrc.selected.sum())
                tempCat.extend(sources[goodSrc.selected], mapper=sourceMapper)
                tempCat[visitKey][:] = visit['visit']
                tempCat[ccdKey][:] = ccdId
                # Compute "magnitude" by scaling flux with exposure time.
                # Add an arbitrary zeropoint that needs to be investigated.
                scaledInstFlux = sources[instFluxKey][goodSrc.selected] * visit['scaling'][ccdIndex]
                tempCat[instMagKey][:] = (-2.5 * np.log10(scaledInstFlux) + 2.5 * np.log10(expTime))
                # instMagErr is computed with original (unscaled) flux
                k = 2.5 / np.log(10.)
                tempCat[instMagErrKey][:] = k * (sources[instFluxErrKey][goodSrc.selected] /
                                                 sources[instFluxKey][goodSrc.selected])

                if self.config.applyJacobian:
                    tempCat[instMagKey][:] -= 2.5 * np.log10(tempCat['jacobian'][:])

                fullCatalog.extend(tempCat)

                # And the aperture information
                tempAperCat = afwTable.BaseCatalog(aperVisitCatalog.schema)
                tempAperCat.reserve(goodSrc.selected.sum())
                tempAperCat.extend(sources[goodSrc.selected], mapper=aperMapper)

                with np.warnings.catch_warnings():
                    # Ignore warnings, we will filter infinities and
                    # nans below.
                    np.warnings.simplefilter("ignore")

                    tempAperCat[instMagInKey][:] = -2.5 * \
                        np.log10(sources[instFluxAperInKey][goodSrc.selected])
                    tempAperCat[instMagErrInKey][:] = (2.5 / np.log(10.)) * (
                        sources[instFluxErrAperInKey][goodSrc.selected] /
                        sources[instFluxAperInKey][goodSrc.selected])
                    tempAperCat[instMagOutKey][:] = -2.5 * np.log10(
                        sources[instFluxAperOutKey][goodSrc.selected])
                    tempAperCat[instMagErrOutKey][:] = (2.5 / np.log(10.)) * (
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

        # Write all the observations
        butler.put(fullCatalog, 'fgcmStarObservations')

        # And overwrite the visitCatalog with delta_aper info
        butler.put(visitCat, 'fgcmVisitCatalog')

        self.log.info("Done with all stars in %.2f s" %
                      (time.time() - startTime))

    def _fgcmMatchStars(self, butler, visitCat):
        """
        Use FGCM code to match observations into unique stars.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for FGCM

        Outputs
        -------
        Butler will output fgcmStarIds (matched star identifiers) and
        fgcmStarIndices (index values linking fgcmStarObservations to
        fgcmStarIds).
        """

        obsCat = butler.get('fgcmStarObservations')

        # get filter names into a numpy array...
        # This is the type that is expected by the fgcm code
        visitFilterNames = np.zeros(len(visitCat), dtype='a2')
        for i in range(len(visitCat)):
            visitFilterNames[i] = visitCat[i]['filtername']

        # match to put filterNames with observations
        visitIndex = np.searchsorted(visitCat['visit'],
                                     obsCat['visit'])

        obsFilterNames = visitFilterNames[visitIndex]

        # make the fgcm starConfig dict

        starConfig = {'logger': self.log,
                      'filterToBand': self.config.filterToBand,
                      'requiredBands': self.config.requiredBands,
                      'minPerBand': self.config.minPerBand,
                      'matchRadius': self.config.matchRadius,
                      'isolationRadius': self.config.isolationRadius,
                      'matchNSide': self.config.matchNside,
                      'coarseNSide': self.config.coarseNside,
                      'densNSide': self.config.densityCutNside,
                      'densMaxPerPixel': self.config.densityCutMaxPerPixel,
                      'primaryBands': self.config.primaryBands,
                      'referenceBands': self.config.referenceBands}

        # initialize the FgcmMakeStars object
        fgcmMakeStars = fgcm.FgcmMakeStars(starConfig)

        # make the primary stars
        #  note that the ra/dec native Angle format is radians
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

        butler.put(fgcmStarIdCat, 'fgcmStarIds')

        obsSchema = self._makeFgcmObsSchema()

        fgcmStarIndicesCat = afwTable.BaseCatalog(obsSchema)
        fgcmStarIndicesCat.reserve(fgcmMakeStars.obsIndexCat.size)
        for i in range(fgcmMakeStars.obsIndexCat.size):
            fgcmStarIndicesCat.addNew()

        fgcmStarIndicesCat['obsIndex'][:] = fgcmMakeStars.obsIndexCat['obsindex']

        butler.put(fgcmStarIndicesCat, 'fgcmStarIndices')

        if self.config.doReferenceMatches:
            refSchema = self._makeFgcmRefSchema(len(self.config.referenceBands))

            fgcmRefCat = afwTable.BaseCatalog(refSchema)
            fgcmRefCat.reserve(fgcmMakeStars.referenceCat.size)

            for i in range(fgcmMakeStars.referenceCat.size):
                fgcmRefCat.addNew()

            fgcmRefCat['fgcm_id'][:] = fgcmMakeStars.referenceCat['fgcm_id']
            fgcmRefCat['refMag'][:, :] = fgcmMakeStars.referenceCat['refMag']
            fgcmRefCat['refMagErr'][:, :] = fgcmMakeStars.referenceCat['refMagErr']

            butler.put(fgcmRefCat, 'fgcmReferenceStars')

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
        schema.addField('filtername', type=str, size=2, doc="Filter name")
        schema.addField('telra', type=np.float64, doc="Pointing RA (deg)")
        schema.addField('teldec', type=np.float64, doc="Pointing Dec (deg)")
        schema.addField('telha', type=np.float64, doc="Pointing Hour Angle (deg)")
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
        sourceMapper.addMapping(sourceSchema[self.config.jacobianName].asKey(),
                                'jacobian')
        sourceMapper.addMapping(sourceSchema['slot_Centroid_x'].asKey(), 'x')
        sourceMapper.addMapping(sourceSchema['slot_Centroid_y'].asKey(), 'y')

        # and add the fields we want
        sourceMapper.editOutputSchema().addField(
            "visit", type=np.int32, doc="Visit number")
        sourceMapper.editOutputSchema().addField(
            "ccd", type=np.int32, doc="CCD number")
        sourceMapper.editOutputSchema().addField(
            "instMag", type=np.float32, doc="Instrumental magnitude")
        sourceMapper.editOutputSchema().addField(
            "instMagErr", type=np.float32, doc="Instrumental magnitude error")

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
