# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from past.builtins import xrange

import sys
import traceback

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from lsst.daf.base.dateTime import DateTime
import lsst.daf.persistence.butlerExceptions as butlerExceptions
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

import time

import fgcm

__all__ = ['FgcmBuildStarsConfig', 'FgcmBuildStarsTask']


class FgcmBuildStarsConfig(pexConfig.Config):
    """Config for FgcmBuildStarsTask"""

    fluxField = pexConfig.Field(
        doc=("Name of the source flux field to use.  The associated flag field\n"
             "('<name>_flag') will be implicitly included in badFlags"),
        dtype=str,
        default='slot_CalibFlux_instFlux',
    )
    remake = pexConfig.Field(
        doc="Remake visit catalog and stars even if they are already in the butler tree.",
        dtype=bool,
        default=False,
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
    referenceBands = pexConfig.ListField(
        doc="Reference bands for primary matches",
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
    renormalizeFlats = pexConfig.Field(
        doc="Renormalize large-scale flat-field changes?",
        dtype=bool,
        default=False
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources",
        default="science"
    )
    apertureInnerFluxField = pexConfig.Field(
        doc="Field that contains inner aperture for aperture correction proxy",
        dtype=str,
        default='base_CircularApertureFlux_12_0_instFlux'
    )
    apertureOuterFluxField = pexConfig.Field(
        doc="Field that contains outer aperture for aperture correction proxy",
        dtype=str,
        default='base_CircularApertureFlux_17_0_instFlux'
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["science"]
        sourceSelector.setDefaults()

        fluxFlagName = self.fluxField[0: -len('instFlux')] + 'flag'

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

        sourceSelector.signalToNoise.fluxField = self.fluxField
        sourceSelector.signalToNoise.errField = self.fluxField + 'Err'
        sourceSelector.signalToNoise.minimum = 10.0
        sourceSelector.signalToNoise.maximum = 1000.0

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

    # TaskClass = FgcmBuildStarsTask

    # only need a single butler instance to run on
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
        args: Tuple with (butler, dataRefList)

        Returns
        -------
        None if self.doReturnResults is False
        A pipe.base.Struct containing these fields if self.doReturnResults is True:
           dataRefList: the provided data references
        """
        butler, dataRefList = args

        task = self.TaskClass(config=self.config, log=self.log)

        exitStatus = 0
        if self.doRaise:
            results = task.runDataRef(butler, dataRefList)
        else:
            try:
                results = task.runDataRef(butler, dataRefList)
            except Exception as e:
                exitStatus = 1
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)

        if self.doReturnResults:
            return [pipeBase.Struct(exitStatus=exitStatus,
                                    results=results)]
        else:
            return [pipeBase.Struct(exitStatus=exitStatus)]

    # turn off any multiprocessing

    def run(self, parsedCmd):
        """
        Run the task, with no multiprocessing

        Parameters
        ----------
        parsedCmd: ArgumentParser parsed command line
        """

        resultList = []

        if self.precall(parsedCmd):
            # profileName = parsedCmd.profile if hasattr(parsedCmd, "profile") else None
            # log = parsedCmd.log
            targetList = self.getTargetList(parsedCmd)
            # And call the runner on the first (and only) item in the list,
            #  which is a tuple of the butler and any dataRefs
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
        butler : lsst.daf.persistence.Butler
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.makeSubtask("sourceSelector")
        # Only log fatal errors from the sourceSelector
        self.sourceSelector.log.setLevel(self.sourceSelector.log.FATAL)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="Data ID, e.g. --id visit=6789 (optional)")

        return parser

    # no saving of the config for now
    # def _getConfigName(self):
    #     return None

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def runDataRef(self, butler, dataRefs):
        """
        Cross-match and make star list for FGCM Input

        Parameters
        ----------
        butler:  lsst.daf.persistence.Butler
        dataRefs: list of lsst.daf.persistence.ButlerDataRef
           Data references for the input visits
           If this is an empty list, all visits with src catalogs in
           the repository are used.
           Only one individual dataRef from a visit need be specified
           and the code will find the other source catalogs from
           each visit

        Returns
        -------
        pipe.base.Struct
            struct containing:
            * dataRefs: the provided data references consolidated
        """

        # Make the visit catalog if necessary
        if self.config.remake or not butler.datasetExists('fgcmVisitCatalog'):
            # we need to build visitCat
            visitCat = self._fgcmMakeVisitCatalog(butler, dataRefs)
        else:
            self.log.info("Found fgcmVisitCatalog.")
            visitCat = butler.get('fgcmVisitCatalog')

        # Compile all the stars
        if self.config.remake or not butler.datasetExists('fgcmStarObservations'):
            self._fgcmMakeAllStarObservations(butler, visitCat)
        else:
            self.log.info("Found fgcmStarObservations")

        if self.config.remake or (not butler.datasetExists('fgcmStarIds') or
                                  not butler.datasetExists('fgcmStarIndices')):
            self._fgcmMatchStars(butler, visitCat)
        else:
            self.log.info("Found fgcmStarIds and fgcmStarIndices")

        # The return value could be the visitCat, if anybody wants that.
        return visitCat

    def _fgcmMakeVisitCatalog(self, butler, dataRefs):
        """
        Make a visit catalog with all the key data from each visit

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
        dataRefs: list of lsst.daf.persistence.ButlerDataRef
           Data references for the input visits
           If this is an empty list, all visits with src catalogs in
           the repository are used.
           Only one individual dataRef from a visit need be specified
           and the code will find the other source catalogs from
           each visit

        Returns
        -------
        visitCat: afw.table.BaseCatalog
        """

        startTime = time.time()

        camera = butler.get('camera')
        nCcd = len(camera)

        if len(dataRefs) == 0:
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
                    for ccd in xrange(nCcd):
                        if (butler.datasetExists('src', dataId={self.config.visitDataRefName: dataset[0],
                                                                self.config.ccdDataRefName:
                                                                    ccd})):
                            srcVisits.append(dataset[0])
                            srcCcds.append(ccd)
                            # Once we find that a butler dataset exists, break out
                            break
        else:
            # get the visits from the datarefs, only for referenceCCD
            srcVisits = [d.dataId[self.config.visitDataRefName] for d in dataRefs if
                         d.dataId[self.config.ccdDataRefName] == self.config.referenceCCD]
            srcCcds = [self.config.referenceCCD] * len(srcVisits)

        # Sort the visits for searching/indexing
        srcVisits.sort()

        self.log.info("Found %d visits in %.2f s" %
                      (len(srcVisits), time.time() - startTime))

        schema = afwTable.Schema()
        schema.addField('visit', type=np.int32, doc="Visit number")
        schema.addField('filtername', type=str, size=2, doc="Filter name")
        schema.addField('telra', type=np.float64, doc="Pointing RA (deg)")
        schema.addField('teldec', type=np.float64, doc="Pointing Dec (deg)")
        schema.addField('telha', type=np.float64, doc="Pointing Hour Angle (deg)")
        schema.addField('mjd', type=np.float64, doc="MJD of visit")
        schema.addField('exptime', type=np.float32, doc="Exposure time")
        schema.addField('pmb', type=np.float32, doc="Pressure (millibar)")
        schema.addField('psfsigma', type=np.float32, doc="PSF sigma (reference CCD)")
        schema.addField('deltaaper', type=np.float32, doc="Delta-aperture")
        schema.addField('skybackground', type=np.float32, doc="Sky background (ADU) (reference CCD)")
        # the following field is not used yet
        schema.addField('deepflag', type=np.int32, doc="Deep observation")
        schema.addField('scaling', type='ArrayD', doc="Scaling applied due to flat adjustment",
                        size=len(camera))

        visitCat = afwTable.BaseCatalog(schema)
        visitCat.table.preallocate(len(srcVisits))

        startTime = time.time()
        # reading in a small bbox is faster for non-gzipped images
        bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.PointI(1, 1))

        # now loop over visits and get the information
        for i, srcVisit in enumerate(srcVisits):
            # We don't use the bypasses since we need the psf info which does
            # not have a bypass

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
            rec['deepflag'] = 0
            rec['scaling'][:] = 1.0
            rec['deltaaper'] = 0.0

            rec['psfsigma'] = psf.computeShape().getDeterminantRadius()

            if butler.datasetExists('calexpBackground', dataId=dataId):
                # Get background for reference CCD
                # This approximation is good enough for now
                bgStats = (bg[0].getStatsImage().getImage().array
                           for bg in butler.get('calexpBackground',
                                                dataId=dataId))
                rec['skybackground'] = sum(np.median(bg[np.isfinite(bg)]) for bg in bgStats)
            else:
                rec['skybackground'] = -1.0

        # Compute flat scaling if desired...
        if self.config.renormalizeFlats:
            self.log.info("Reading flats for renormalizeFlats")
            scalingValues = self._computeFlatScaling(butler, visitCat)
            visitCat['scaling'] *= scalingValues

        self.log.info("Found all VisitInfo in %.2f s" % (time.time() - startTime))

        # and now persist it
        butler.put(visitCat, 'fgcmVisitCatalog')

        return visitCat

    def _fgcmMakeAllStarObservations(self, butler, visitCat):
        """
        Compile all good star observations from visits in visitCat

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
        visitCat: afw.table.BaseCatalog
           Catalog with visit data for FGCM

        Returns
        -------
        None
        """

        startTime = time.time()

        # create our source schema
        sourceSchema = butler.get('src_schema', immediate=True).schema

        # create a mapper to the preferred output
        sourceMapper = afwTable.SchemaMapper(sourceSchema)

        # map to ra/dec
        sourceMapper.addMapping(sourceSchema.find('coord_ra').key, 'ra')
        sourceMapper.addMapping(sourceSchema.find('coord_dec').key, 'dec')
        sourceMapper.addMapping(sourceSchema.find(self.config.jacobianName).key,
                                'jacobian')

        # and add the fields we want
        sourceMapper.editOutputSchema().addField(
            "visit", type=np.int32, doc="Visit number")
        sourceMapper.editOutputSchema().addField(
            "ccd", type=np.int32, doc="CCD number")
        sourceMapper.editOutputSchema().addField(
            "mag", type=np.float32, doc="Raw magnitude")
        sourceMapper.editOutputSchema().addField(
            "magerr", type=np.float32, doc="Raw magnitude error")

        # We also have a temporary catalog that will accumulate aperture measurements

        aperMapper = afwTable.SchemaMapper(sourceSchema)
        aperMapper.addMapping(sourceSchema.find('coord_ra').key, 'ra')
        aperMapper.addMapping(sourceSchema.find('coord_dec').key, 'dec')
        aperMapper.editOutputSchema().addField('mag_aper_inner', type=np.float64,
                                               doc="Magnitude at inner aperture")
        aperMapper.editOutputSchema().addField('magerr_aper_inner', type=np.float64,
                                               doc="Magnitude error at inner aperture")
        aperMapper.editOutputSchema().addField('mag_aper_outer', type=np.float64,
                                               doc="Magnitude at outer aperture")
        aperMapper.editOutputSchema().addField('magerr_aper_outer', type=np.float64,
                                               doc="Magnitude error at outer aperture")

        # we need to know the ccds...
        camera = butler.get('camera')

        started = False

        # loop over visits
        for visit in visitCat:
            expTime = visit['exptime']

            nStarInVisit = 0

            # The temporary aperture catalog needs to be reset
            aperStarted = False

            # loop over CCDs
            for ccdIndex, detector in enumerate(camera):
                ccdId = detector.getId()

                try:
                    # Need to cast visit['visit'] to python int because butler
                    # can't use numpy ints
                    sources = butler.get('src', dataId={self.config.visitDataRefName:
                                                        int(visit['visit']),
                                                        self.config.ccdDataRefName: ccdId},
                                         flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                except butlerExceptions.NoResults:
                    # this is not a problem if this ccd isn't there
                    continue

                if not started:
                    # get the keys for quicker look-up

                    # Calibration is based on configuration fluxField
                    fluxKey = sources.schema[self.config.fluxField].asKey()
                    fluxErrKey = sources.schema[self.config.fluxField + 'Err'].asKey()

                    outputSchema = sourceMapper.getOutputSchema()
                    visitKey = outputSchema['visit'].asKey()
                    ccdKey = outputSchema['ccd'].asKey()
                    magKey = outputSchema['mag'].asKey()
                    magErrKey = outputSchema['magerr'].asKey()

                    # and the final part of the sourceMapper
                    sourceMapper.addMapping(sources.schema['slot_Centroid_x'].asKey(), 'x')
                    sourceMapper.addMapping(sources.schema['slot_Centroid_y'].asKey(), 'y')

                    # Create a stub of the full catalog
                    fullCatalog = afwTable.BaseCatalog(sourceMapper.getOutputSchema())

                    started = True

                if not aperStarted:
                    # And the aperture catalog
                    fluxAperInKey = sources.schema[self.config.apertureInnerFluxField].asKey()
                    fluxErrAperInKey = sources.schema[self.config.apertureInnerFluxField + 'Err'].asKey()
                    fluxAperOutKey = sources.schema[self.config.apertureOuterFluxField].asKey()
                    fluxErrAperOutKey = sources.schema[self.config.apertureOuterFluxField + 'Err'].asKey()

                    aperOutputSchema = aperMapper.getOutputSchema()
                    magInKey = aperOutputSchema['mag_aper_inner'].asKey()
                    magErrInKey = aperOutputSchema['magerr_aper_inner'].asKey()
                    magOutKey = aperOutputSchema['mag_aper_outer'].asKey()
                    magErrOutKey = aperOutputSchema['magerr_aper_outer'].asKey()

                    aperVisitCatalog = afwTable.BaseCatalog(aperMapper.getOutputSchema())

                    aperStarted = True

                goodSrc = self.sourceSelector.selectSources(sources)

                tempCat = afwTable.BaseCatalog(fullCatalog.schema)
                tempCat.reserve(goodSrc.selected.sum())
                tempCat.extend(sources[goodSrc.selected], mapper=sourceMapper)
                tempCat[visitKey][:] = visit['visit']
                tempCat[ccdKey][:] = ccdId
                # Compute "magnitude" by scaling flux with exposure time.
                # Add an arbitrary zeropoint that needs to be investigated.
                scaledFlux = sources[fluxKey][goodSrc.selected] * visit['scaling'][ccdIndex]
                tempCat[magKey][:] = (-2.5 * np.log10(scaledFlux) +
                                      2.5 * np.log10(expTime))
                # magErr is computed with original (unscaled) flux
                tempCat[magErrKey][:] = (2.5 / np.log(10.)) * (sources[fluxErrKey][goodSrc.selected] /
                                                               sources[fluxKey][goodSrc.selected])

                if self.config.applyJacobian:
                    tempCat[magKey][:] -= 2.5 * np.log10(tempCat['jacobian'][:])

                fullCatalog.extend(tempCat)

                # And the aperture information
                tempAperCat = afwTable.BaseCatalog(aperVisitCatalog.schema)
                tempAperCat.reserve(goodSrc.selected.sum())
                tempAperCat.extend(sources[goodSrc.selected], mapper=aperMapper)
                tempAperCat[magInKey][:] = -2.5 * np.log10(sources[fluxAperInKey][goodSrc.selected])
                tempAperCat[magErrInKey][:] = (2.5 / np.log(10.)) * (
                    sources[fluxErrAperInKey][goodSrc.selected] /
                    sources[fluxAperInKey][goodSrc.selected])
                tempAperCat[magOutKey][:] = -2.5 * np.log10(
                    sources[fluxAperOutKey][goodSrc.selected])
                tempAperCat[magErrOutKey][:] = (2.5 / np.log(10.)) * (
                    sources[fluxErrAperOutKey][goodSrc.selected] /
                    sources[fluxAperOutKey][goodSrc.selected])

                aperVisitCatalog.extend(tempAperCat)

                nStarInVisit += len(tempCat)

            # Compute the median delta-aper
            if not aperVisitCatalog.isContiguous():
                aperVisitCatalog = aperVisitCatalog.copy(deep=True)

            magIn = aperVisitCatalog[magInKey]
            magErrIn = aperVisitCatalog[magErrInKey]
            magOut = aperVisitCatalog[magOutKey]
            magErrOut = aperVisitCatalog[magErrOutKey]

            ok = (np.isfinite(magIn) & np.isfinite(magErrIn) & np.isfinite(magOut) & np.isfinite(magErrOut))

            visit['deltaaper'] = np.median(magIn[ok] - magOut[ok])

            self.log.info("  Found %d good stars in visit %d (delta_aper = %.3f)" %
                          (nStarInVisit, visit['visit'], visit['deltaaper']))

        self.log.info("Found all good star observations in %.2f s" %
                      (time.time() - startTime))

        # Write all the observations
        butler.put(fullCatalog, 'fgcmStarObservations')

        # And overwrite the visitCatalog with delta_aper info
        butler.put(visitCat, 'fgcmVisitCatalog')

        self.log.info("Done with all stars in %.2f s" %
                      (time.time() - startTime))
        return None

    def _fgcmMatchStars(self, butler, visitCat):
        """
        Use FGCM code to match observations into unique stars.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
        visitCat: afw.table.BaseCatalog
           Catalog with visit data for FGCM

        Returns
        -------
        None
        """

        obsCat = butler.get('fgcmStarObservations')

        # get filter names into a numpy array...
        visitFilterNames = np.zeros(len(visitCat), dtype='a2')
        for i in xrange(len(visitCat)):
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
                      'referenceBands': self.config.referenceBands}

        # initialize the FgcmMakeStars object
        fgcmMakeStars = fgcm.FgcmMakeStars(starConfig)

        # make the reference stars
        #  note that the ra/dec native Angle format is radians
        fgcmMakeStars.makeReferenceStars(np.rad2deg(obsCat['ra']),
                                         np.rad2deg(obsCat['dec']),
                                         filterNameArray=obsFilterNames,
                                         bandSelected=False)

        # and match all the stars
        fgcmMakeStars.makeMatchedStars(np.rad2deg(obsCat['ra']),
                                       np.rad2deg(obsCat['dec']),
                                       obsFilterNames)

        # now persist

        # afwTable for objects
        objSchema = afwTable.Schema()
        objSchema.addField('fgcm_id', type=np.int32, doc='FGCM Unique ID')
        # FIXME: should be angle?
        objSchema.addField('ra', type=np.float64, doc='Mean object RA')
        objSchema.addField('dec', type=np.float64, doc='Mean object Dec')
        objSchema.addField('obsarrindex', type=np.int32,
                           doc='Index in obsIndexTable for first observation')
        objSchema.addField('nobs', type=np.int32, doc='Total number of observations')

        # make catalog and records
        fgcmStarIdCat = afwTable.BaseCatalog(objSchema)
        fgcmStarIdCat.table.preallocate(fgcmMakeStars.objIndexCat.size)
        for i in xrange(fgcmMakeStars.objIndexCat.size):
            fgcmStarIdCat.addNew()

        # fill the catalog
        fgcmStarIdCat['fgcm_id'][:] = fgcmMakeStars.objIndexCat['fgcm_id']
        fgcmStarIdCat['ra'][:] = fgcmMakeStars.objIndexCat['ra']
        fgcmStarIdCat['dec'][:] = fgcmMakeStars.objIndexCat['dec']
        fgcmStarIdCat['obsarrindex'][:] = fgcmMakeStars.objIndexCat['obsarrindex']
        fgcmStarIdCat['nobs'][:] = fgcmMakeStars.objIndexCat['nobs']

        butler.put(fgcmStarIdCat, 'fgcmStarIds')

        # afwTable for observation indices
        obsSchema = afwTable.Schema()
        obsSchema.addField('obsindex', type=np.int32, doc='Index in observation table')

        fgcmStarIndicesCat = afwTable.BaseCatalog(obsSchema)
        fgcmStarIndicesCat.table.preallocate(fgcmMakeStars.obsIndexCat.size)
        for i in xrange(fgcmMakeStars.obsIndexCat.size):
            fgcmStarIndicesCat.addNew()

        fgcmStarIndicesCat['obsindex'][:] = fgcmMakeStars.obsIndexCat['obsindex']

        butler.put(fgcmStarIndicesCat, 'fgcmStarIndices')

        # and we're done with the stars
        return None

    def _computePsfSigma(self, butler, visitCat):
        """
        """
        startTime = time.time()

        camera = butler.get('camera')

        bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.PointI(1, 1))

        psfSigma = np.zeros((len(visitCat), len(camera)))

        for visitIndex, vis in enumerate(visitCat):
            self.log.info(' Working on %d' % (vis['visit']))
            visit = vis['visit']
            for ccdIndex, detector in enumerate(camera):
                dataId = {'visit': int(visit),
                          'ccd': detector.getId()}

                if not butler.datasetExists('calexp', dataId=dataId):
                    continue
                exp = butler.get('calexp_sub', dataId=dataId, bbox=bbox,
                                 flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                psfSigma[visitIndex, ccdIndex] = exp.getPsf().computeShape().getDeterminantRadius()

        self.log.info("Computed psfs from %d visits in %.2f s" %
                      (len(visitCat), time.time() - startTime))
        return psfSigma

    def _computeSkyBackground(self, butler, visitCat):
        """
        """
        startTime = time.time()

        camera = butler.get('camera')

        skyBackground = np.zeros((len(visitCat), len(camera)))

        for visitIndex, vis in enumerate(visitCat):
            visit = vis['visit']
            for ccdIndex, detector in enumerate(camera):
                dataId = {'visit': int(visit),
                          'ccd': detector.getId()}

                if not butler.datasetExists('calexpBackground', dataId=dataId):
                    continue

                bgStats = (bg[0].getStatsImage().getImage().array
                           for bg in butler.get('calexpBackground',
                                                dataId=dataId))
                skyBackground[visitIndex, ccdIndex] = sum(np.median(bg[np.isfinite(bg)]) for bg in bgStats)

        self.log.info("Computed background from %d visits in %.2f s" %
                      (len(visitCat), time.time() - startTime))
        return skyBackground

    def _computeFlatScaling(self, butler, visitCat):
        """
        """

        # Get the maximum ccd id.
        # Note that this (and other assumptions) will have to be
        # rethought if we don't have integer ccds...
        camera = butler.get('camera')
        nCcd = len(camera)

        # A dictionary to map Flats to visit/ccd pars
        visitCcdFlatDict = {}

        # A dictionary to store flat values
        flatValueDict = {}

        # And a string that will be unique for my internal keys
        joiner = '++'

        for vis in visitCat:
            visit = vis['visit']
            for ccdIndex, detector in enumerate(camera):
                dataId = {'visit': int(visit),
                          'ccd': detector.getId()}
                if not butler.datasetExists('src', dataId=dataId):
                    continue

                flatRef = butler.dataRef('flat', dataId=dataId)

                fName = flatRef.dataId['filter']
                cDate = flatRef.dataId['calibDate']

                flatKey = '%s%s%s%s%04d' % (cDate, joiner, fName, joiner, ccdIndex)

                if flatKey not in flatValueDict:
                    # Read in the flat and record the value
                    self.log.info("Found new flat: %s" % (flatKey))
                    flat = flatRef.get()
                    flatValueDict[flatKey] = np.median(flat.getImage().getArray())

                visitCcdKey = visit * (nCcd + 1) + ccdIndex
                visitCcdFlatDict[visitCcdKey] = flatKey

        # Group flats together in a dict
        flatFields = {}
        for key in flatValueDict:
            parts = key.split(joiner)
            flatName = parts[0] + joiner + parts[1]
            if flatName not in flatFields:
                flatFields[flatName] = (parts[1], np.zeros(nCcd))
            flatFields[flatName][1][int(parts[2])] = flatValueDict[key]

        # And group by filter (uniquely)...
        flatFilters = {flatFields[x][0] for x in flatFields}

        # Compute scaling
        # Here's the math...
        # When flat fields are applied, we have:
        #  maskedImage.scaledDivides(1.0 / flatScale, flatMaskedImage)
        # Right now, I'm assuming that flatScale == 1.0, but I don't know where
        # this is recorded.
        # And scaledDivides is "Divide lhs by c * rhs"
        #  maskedImageNew = maskedImage / ((1.0 / flatScale) * flatMaskedImage)
        # To remove the flat (on CCD scale) we need to *multiply* by median(flatField)
        # Then to apply the reference flat we need to *divide* by median(referenceFlat)
        # So the (multiplicative) scaling to switch from flatField to referenceFlat is:
        #  scaling = median(flatField) / median(referenceFlat)

        flatScaleDict = flatValueDict.copy()
        for flatFilter in flatFilters:
            # get all the flats which have this...
            flatKeys = [key for key in flatFields if flatFields[key][0] == flatFilter]

            # take the first one as the arbitrary reference
            referenceFlat = flatFields[flatKeys[0]][1]

            # Loop over all the flats and scale to reference
            for flatKey in flatKeys:
                scaleVals = np.ones_like(referenceFlat)
                u, = np.where((referenceFlat > 0) & (flatFields[flatKey][1] > 0))
                scaleVals[u] = flatFields[flatKey][1][u] / referenceFlat[u]

                for ccdIndex in range(nCcd):
                    flatCcdKey = '%s%s%04d' % (flatKey, joiner, ccdIndex)
                    if flatCcdKey in flatScaleDict:
                        flatScaleDict[flatCcdKey] = scaleVals[ccdIndex]

        # And compute scaling values for each visit/ccd pair
        scalingValues = np.ones((len(visitCat), nCcd))
        for visitIndex, vis in enumerate(visitCat):
            visit = vis['visit']
            for ccdIndex, detector in enumerate(camera):
                visitCcdKey = visit * (nCcd + 1) + ccdIndex
                try:
                    scalingValues[visitIndex, ccdIndex] = flatScaleDict[visitCcdFlatDict[visitCcdKey]]
                except KeyError:
                    pass

        return scalingValues
