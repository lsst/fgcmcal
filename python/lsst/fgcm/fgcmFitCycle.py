# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import sys
import traceback

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as afwCameraGeom

import fgcm

__all__ = ['FgcmFitCycleConfig', 'FgcmFitCycleTask']


class FgcmFitCycleConfig(pexConfig.Config):
    """Config for FgcmFitCycle"""

    bands = pexConfig.ListField(
        doc="Bands to run calibration (in wavelength order)",
        dtype=str,
        default=("NO_DATA",),
    )
    fitFlag = pexConfig.ListField(
        doc="Flag for bands to fit",
        dtype=int,
        default=(0,),
    )
    filterToBand = pexConfig.DictField(
        doc="filterName to band mapping",
        keytype=str,
        itemtype=str,
        default={},
    )
    nCore = pexConfig.Field(
        doc="Number of cores to use",
        dtype=int,
        default=4,
    )
    nStarPerRun = pexConfig.Field(
        doc="Number of stars to run in each chunk",
        dtype=int,
        default=200000,
    )
    nExpPerRun = pexConfig.Field(
        doc="Number of exposures to run in each chunk",
        dtype=int,
        default=1000,
    )
    reserveFraction = pexConfig.Field(
        doc="Fraction of stars to reserve for testing",
        dtype=float,
        default=0.1,
    )
    freezeStdAtmosphere = pexConfig.Field(
        doc="Freeze atmosphere parameters to standard (for testing)",
        dtype=bool,
        default=False,
    )
    precomputeSuperStarInitialCycle = pexConfig.Field(
        doc="Precompute superstar flat for initial cycle",
        dtype=bool,
        default=False,
    )
    cycleNumber = pexConfig.Field(
        doc="Fit Cycle Number",
        dtype=int,
        default=None,
    )
    maxIter = pexConfig.Field(
        doc="Max iterations",
        dtype=int,
        default=50,
    )
    utBoundary = pexConfig.Field(
        doc="Boundary (in UTC) from day-to-day",
        dtype=float,
        default=None,
    )
    washMjds = pexConfig.ListField(
        doc="Mirror wash MJDs",
        dtype=float,
        default=(0.0,),
    )
    epochMjds = pexConfig.ListField(
        doc="Epoch boundaries in MJD",
        dtype=float,
        default=(0.0,),
    )
    minObsPerBand = pexConfig.Field(
        doc="Minimum good observations per band",
        dtype=int,
        default=2,
    )
    latitude = pexConfig.Field(
        doc="Observatory latitude",
        dtype=float,
        default=None,
    )
    pixelScale = pexConfig.Field(
        doc="Pixel scale (arcsec/pixel) (temporary)",
        dtype=float,
        default=None,
    )
    brightObsGrayMax = pexConfig.Field(
        doc="Maximum gray extinction to be considered bright observation",
        dtype=float,
        default=0.15,
    )
    minStarPerCcd = pexConfig.Field(
        doc="Minimum number of good stars per CCD for calibration",
        dtype=int,
        default=5,
    )
    minCcdPerExp = pexConfig.Field(
        doc="Minimum number of good CCDs per exposure",
        dtype=int,
        default=5,
    )
    maxCcdGrayErr = pexConfig.Field(
        doc="Maximum error on CCD gray offset to be considered good",
        dtype=float,
        default=0.05,
    )
    minStarPerExp = pexConfig.Field(
        doc="Minimum number of good stars per exposure to be considered good",
        dtype=int,
        default=600,
    )
    minExpPerNight = pexConfig.Field(
        doc="Minimum number of good exposures to consider a good night",
        dtype=int,
        default=10,
    )
    expGrayInitialCut = pexConfig.Field(
        doc="Maximum exposure gray value for initial cut",
        dtype=float,
        default=-0.25,
    )
    expGrayPhotometricCut = pexConfig.ListField(
        doc="Negative exposure gray cut for photometric selection",
        dtype=float,
        default=(0.0,),
    )
    expGrayHighCut = pexConfig.ListField(
        doc="Positive exposure gray cut for photometric selection",
        dtype=float,
        default=(0.0,),
    )
    expGrayRecoverCut = pexConfig.Field(
        doc="Maximum exposure gray to be able to recover bad ccds",
        dtype=float,
        default=-1.0,
    )
    expVarGrayPhotometricCut = pexConfig.Field(
        doc="Maximum exposure variance to be considered possibly photometric",
        dtype=float,
        default=0.0005,
    )
    expGrayErrRecoverCut = pexConfig.Field(
        doc="Maximum exposure gray error to be able to recover bad ccds",
        dtype=float,
        default=0.05,
    )
    illegalValue = pexConfig.Field(
        doc="Sentinal value for no-values",
        dtype=float,
        default=-9999.0,
    )
    aperCorrFitNBins = pexConfig.Field(
        doc="Aperture correction number of bins",
        dtype=int,
        default=None,
    )
    sedFudgeFactors = pexConfig.ListField(
        doc="Fudge factors for computing linear SED from colors",
        dtype=float,
        default=(0,),
    )
    sigFgcmMaxErr = pexConfig.Field(
        doc="Maximum mag error for fitting sigma_FGCM",
        dtype=float,
        default=0.01,
    )
    sigFgcmMaxEGray = pexConfig.Field(
        doc="Maximum (absolute) gray value for observation in sigma_FGCM",
        dtype=float,
        default=0.05,
    )
    ccdGrayMaxStarErr = pexConfig.Field(
        doc="Maximum error on a star observation to use in ccd gray computation",
        dtype=float,
        default=0.10,
    )
    approxThroughput = pexConfig.Field(
        doc="Approximate overall throughput at start of calibration observations",
        dtype=float,
        default=1.0,
    )
    sigma0Cal = pexConfig.Field(
        doc="Systematic error floor for all observations",
        dtype=float,
        default=0.003,
    )
    sigma0Phot = pexConfig.Field(
        doc="Systematic error floor for all zeropoints",
        dtype=float,
        default=0.003,
    )
    mapLongitudeRef = pexConfig.Field(
        doc="Reference longitude for plotting maps",
        dtype=float,
        default=0.0,
    )
    mapNSide = pexConfig.Field(
        doc="Healpix nside for plotting maps",
        dtype=int,
        default=256,
    )
    varNSig = pexConfig.Field(
        doc="Number of sigma to be tested as a variable",
        dtype=float,
        default=4.0,
    )
    varMinBand = pexConfig.Field(
        doc="Minimum number of bands with variability to be flagged as variable",
        dtype=int,
        default=2,
    )
    cameraGain = pexConfig.Field(
        doc="Gain value for the typical CCD",
        dtype=float,
        default=None,
    )
    outfileBase = pexConfig.Field(
        doc="Filename start for plot output files",
        dtype=str,
        default=None,
    )
    starColorCuts = pexConfig.ListField(
        doc="Encoded star-color cuts (to be cleaned up)",
        dtype=str,
        default=("NO_DATA",),
    )

    def setDefaults(self):
        pass


class FgcmFitCycleRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmFitCycleTask

    fgcmFitCycleTask.run() takes one argument, the butler, and uses
    stars and visits previously extracted from dataRefs by
    fgcmBuildStars.
    This runner does not use any parallelization, although the
    FGCM code uses multiprocessing
    """

    @staticmethod
    def getTargetList(parsedCmd):
        """
        Return a list with one element, the butler.
        """
        return [parsedCmd.butler]

    def precall(self, parsedCmd):
        return True

    def __call__(self, butler):
        """
        Parameters
        ----------
        butler: lsst.daf.persistence.Butler

        Returns
        -------
        None if self.doReturnResults is False
        An empty list if self.doReturnResults is True
        """

        task = self.TaskClass(config=self.config, log=self.log)
        if self.doRaise:
            results = task.run(butler)
        else:
            try:
                results = task.run(butler)
            except Exception as e:
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)
        if self.doReturnResults:
            return results

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
            # make sure that we only get 1
            resultList = self(targetList[0])

        return resultList


class FgcmFitCycleTask(pipeBase.CmdLineTask):
    """
    Run Single fit cycle for FGCM global calibration
    """

    ConfigClass = FgcmFitCycleConfig
    RunnerClass = FgcmFitCycleRunner
    _DefaultName = "fgcmFitCycle"

    def __init__(self, butler=None, **kwargs):
        """
        Instantiate an fgcmFitCycle.

        Parameters
        ----------
        butler : lsst.daf.persistence.Butler
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)

        return parser

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def run(self, butler):
        """
        Run a single fit cycle for FGCM

        Parameters
        ----------
        butler:  lsst.daf.persistence.Butler

        Returns
        -------
        Empty list
        """

        self._fgcmFitCycle(butler)

        return []

    def _fgcmFitCycle(self, butler):
        """
        Run the fit cycle

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
           (used for mapper information)

        """

        # FIXME:
        #   more sensible configuration options related to bands/fitbands
        #   check that array lengths are matched

        # FIXME:
        #   Need to be able to turn off plots if desired

        #  TBD:
        #   updating configuration at the end for next cycle?
        #   where to output logging and figures?

        if not butler.datasetExists('fgcmVisitCatalog'):
            raise ValueError("Could not find fgcmVisitCatalog in repo!")
        if not butler.datasetExists('fgcmStarObservations'):
            raise ValueError("Could not find fgcmStarObservations in repo!")
        if not butler.datasetExists('fgcmStarIds'):
            raise ValueError("Could not find fgcmStarIds in repo!")
        if not butler.datasetExists('fgcmStarIndices'):
            raise ValueError("Could not find fgcmStarIndices in repo!")
        if not butler.datasetExists('fgcmLookUpTable'):
            raise ValueError("Could not find fgcmLookUpTable in repo!")

        # Need additional datasets if we are not the initial cycle
        if (self.config.cycleNumber > 0):
            if not butler.datasetExists('fgcmFitParameters',
                                        fgcmcycle=self.config.cycleNumber-1):
                raise ValueError("Could not find fgcmFitParameters for previous cycle (%d) in repo!" %
                                 (self.config.cycleNumber-1))
            if not butler.datasetExists('fgcmFlaggedStars',
                                        fgcmcycle=self.config.cycleNumber-1):
                raise ValueError("Could not find fgcmFlaggedStars for previous cycle (%d) in repo!" %
                                 (self.config.cycleNumber-1))

        # FIXME:
        #  check config variables for valid ranges

        bands = np.array(self.config.bands)
        fitFlag = np.array(self.config.fitFlag, dtype=np.bool)

        camera = butler.get('camera')

        # process the starColorCuts
        starColorCutList = []
        for ccut in self.config.starColorCuts:
            parts = ccut.split(',')
            starColorCutList.append([parts[0], parts[1], float(parts[2]), float(parts[3])])

        # create a configuration dictionary for fgcmFitCycle
        configDict = {'outfileBase': self.config.outfileBase,
                      'logger': self.log,
                      'exposureFile': None,
                      'obsFile': None,
                      'indexFile': None,
                      'lutFile': None,
                      'ccdOffsetFile': None,
                      'mirrorArea': np.pi*(camera.telescopeDiameter/2.)**2.,
                      'cameraGain': self.config.cameraGain,
                      'ccdStartIndex': camera[0].getId(),
                      'expField': 'VISIT',
                      'ccdField': 'CCD',
                      'seeingField': 'SEEING',  # FIXME
                      'deepFlag': 'DEEPFLAG',  # unused
                      'bands': bands,
                      'fitBands': bands[fitFlag],
                      'extraBands': bands[~fitFlag],
                      'filterToBand': self.config.filterToBand,
                      'logLevel': 'INFO',  # FIXME
                      'nCore': self.config.nCore,
                      'nStarPerRun': self.config.nStarPerRun,
                      'nExpPerRun': self.config.nExpPerRun,
                      'reserveFraction': self.config.reserveFraction,
                      'freezeStdAtmosphere': self.config.freezeStdAtmosphere,
                      'precomputeSuperStarInitialCycle': self.config.precomputeSuperStarInitialCycle,
                      'cycleNumber': self.config.cycleNumber,
                      'maxIter': self.config.maxIter,
                      'UTBoundary': self.config.utBoundary,
                      'washMJDs': self.config.washMjds,
                      'epochMJDs': self.config.epochMjds,
                      'minObsPerBand': self.config.minObsPerBand,
                      'latitude': self.config.latitude,
                      'brightObsGrayMax': self.config.brightObsGrayMax,
                      'minStarPerCCD': self.config.minStarPerCcd,
                      'minCCDPerExp': self.config.minCcdPerExp,
                      'maxCCDGrayErr': self.config.maxCcdGrayErr,
                      'minStarPerExp': self.config.minStarPerExp,
                      'minExpPerNight': self.config.minExpPerNight,
                      'expGrayInitialCut': self.config.expGrayInitialCut,
                      'expGrayPhotometricCut': np.array(self.config.expGrayPhotometricCut),
                      'expGrayHighCut': np.array(self.config.expGrayHighCut),
                      'expGrayRecoverCut': self.config.expGrayRecoverCut,
                      'expVarGrayPhotometricCut': self.config.expVarGrayPhotometricCut,
                      'expGrayErrRecoverCut': self.config.expGrayErrRecoverCut,
                      'illegalValue': self.config.illegalValue,
                      'starColorCuts': starColorCutList,
                      'aperCorrFitNBins': self.config.aperCorrFitNBins,
                      'sedFitBandFudgeFactors': np.array(self.config.sedFudgeFactors)[fitFlag],
                      'sedExtraBandFudgeFactors': np.array(self.config.sedFudgeFactors)[~fitFlag],
                      'sigFgcmMaxErr': self.config.sigFgcmMaxErr,
                      'sigFgcmMaxEGray': self.config.sigFgcmMaxEGray,
                      'ccdGrayMaxStarErr': self.config.ccdGrayMaxStarErr,
                      'approxThroughput': self.config.approxThroughput,
                      'sigma0Cal': self.config.sigma0Cal,
                      'sigma0Phot': self.config.sigma0Phot,
                      'mapLongitudeRef': self.config.mapLongitudeRef,
                      'mapNSide': self.config.mapNSide,
                      'varNSig': self.config.varNSig,
                      'varMinBand': self.config.varMinBand,
                      'useRetrievedPWV': False,
                      'useNightlyRetrievedPWV': False,
                      'pwvRetrievalSmoothBlock': 25,
                      'useRetrievedTauInit': False,
                      'tauRetrievalMinCCDPerNight': 500,
                      'printOnly': False,
                      'outputStars': False,
                      'clobber': True,
                      'useSedLUT': False,
                      'resetParameters': True}

        # set up the look-up-table
        lutCat = butler.get('fgcmLookUpTable')

        # first we need the lutIndexVals
        lutFilterNames = np.array(lutCat[0]['filternames'].split(','))
        lutStdFilterNames = np.array(lutCat[0]['stdfilternames'].split(','))

        # FIXME: check that lutBands equal listed bands!

        lutIndexVals = np.zeros(1, dtype=[('FILTERNAMES', lutFilterNames.dtype.str,
                                           lutFilterNames.size),
                                          ('STDFILTERNAMES', lutStdFilterNames.dtype.str,
                                           lutStdFilterNames.size),
                                          ('PMB', 'f8', lutCat[0]['pmb'].size),
                                          ('PMBFACTOR', 'f8', lutCat[0]['pmbfactor'].size),
                                          ('PMBELEVATION', 'f8'),
                                          ('LAMBDANORM', 'f8'),
                                          ('PWV', 'f8', lutCat[0]['pwv'].size),
                                          ('O3', 'f8', lutCat[0]['o3'].size),
                                          ('TAU', 'f8', lutCat[0]['tau'].size),
                                          ('ALPHA', 'f8', lutCat[0]['alpha'].size),
                                          ('ZENITH', 'f8', lutCat[0]['zenith'].size),
                                          ('NCCD', 'i4')])

        lutIndexVals['FILTERNAMES'][:] = lutFilterNames
        lutIndexVals['STDFILTERNAMES'][:] = lutStdFilterNames
        lutIndexVals['PMB'][:] = lutCat[0]['pmb']
        lutIndexVals['PMBFACTOR'][:] = lutCat[0]['pmbfactor']
        lutIndexVals['PMBELEVATION'] = lutCat[0]['pmbelevation']
        lutIndexVals['LAMBDANORM'] = lutCat[0]['lambdanorm']
        lutIndexVals['PWV'][:] = lutCat[0]['pwv']
        lutIndexVals['O3'][:] = lutCat[0]['o3']
        lutIndexVals['TAU'][:] = lutCat[0]['tau']
        lutIndexVals['ALPHA'][:] = lutCat[0]['alpha']
        lutIndexVals['ZENITH'][:] = lutCat[0]['zenith']
        lutIndexVals['NCCD'] = lutCat[0]['nccd']

        # now we need the Standard Values
        lutStd = np.zeros(1, dtype=[('PMBSTD', 'f8'),
                                    ('PWVSTD', 'f8'),
                                    ('O3STD', 'f8'),
                                    ('TAUSTD', 'f8'),
                                    ('ALPHASTD', 'f8'),
                                    ('ZENITHSTD', 'f8'),
                                    ('LAMBDARANGE', 'f8', 2),
                                    ('LAMBDASTEP', 'f8'),
                                    ('LAMBDASTD', 'f8', lutFilterNames.size),
                                    ('LAMBDASTDFILTER', 'f8', lutStdFilterNames.size),
                                    ('I0STD', 'f8', lutFilterNames.size),
                                    ('I1STD', 'f8', lutFilterNames.size),
                                    ('I10STD', 'f8', lutFilterNames.size),
                                    ('LAMBDAB', 'f8', lutFilterNames.size),
                                    ('ATMLAMBDA', 'f8', lutCat[0]['atmlambda'].size),
                                    ('ATMSTDTRANS', 'f8', lutCat[0]['atmstdtrans'].size)])
        lutStd['PMBSTD'] = lutCat[0]['pmbstd']
        lutStd['PWVSTD'] = lutCat[0]['pwvstd']
        lutStd['O3STD'] = lutCat[0]['o3std']
        lutStd['TAUSTD'] = lutCat[0]['taustd']
        lutStd['ALPHASTD'] = lutCat[0]['alphastd']
        lutStd['ZENITHSTD'] = lutCat[0]['zenithstd']
        lutStd['LAMBDARANGE'][:] = lutCat[0]['lambdarange'][:]
        lutStd['LAMBDASTEP'] = lutCat[0]['lambdastep']
        lutStd['LAMBDASTD'][:] = lutCat[0]['lambdastd']
        lutStd['LAMBDASTDFILTER'][:] = lutCat[0]['lambdastdfilter']
        lutStd['I0STD'][:] = lutCat[0]['i0std']
        lutStd['I1STD'][:] = lutCat[0]['i1std']
        lutStd['I10STD'][:] = lutCat[0]['i10std']
        lutStd['LAMBDAB'][:] = lutCat[0]['lambdab']
        lutStd['ATMLAMBDA'][:] = lutCat[0]['atmlambda'][:]
        lutStd['ATMSTDTRANS'][:] = lutCat[0]['atmstdtrans'][:]

        lutTypes = []
        for row in lutCat:
            lutTypes.append(row['luttype'])

        # And the flattened look-up-table
        lutFlat = np.zeros(lutCat[0]['lut'].size, dtype=[('I0', 'f4'),
                                                         ('I1', 'f4')])

        try:
            lutFlat['I0'][:] = lutCat[lutTypes.index('I0')]['lut'][:]
            lutFlat['I1'][:] = lutCat[lutTypes.index('I1')]['lut'][:]
        except:
            # need to raise exception
            pass

        lutDerivFlat = np.zeros(lutCat[0]['lut'].size, dtype=[('D_PWV', 'f4'),
                                                              ('D_O3', 'f4'),
                                                              ('D_LNTAU', 'f4'),
                                                              ('D_ALPHA', 'f4'),
                                                              ('D_SECZENITH', 'f4'),
                                                              ('D_PWV_I1', 'f4'),
                                                              ('D_O3_I1', 'f4'),
                                                              ('D_LNTAU_I1', 'f4'),
                                                              ('D_ALPHA_I1', 'f4'),
                                                              ('D_SECZENITH_I1', 'f4')])

        try:
            for name in lutDerivFlat.dtype.names:
                lutDerivFlat[name][:] = lutCat[lutTypes.index(name)]['lut'][:]
        except:
            # raise a helpful exception
            pass

        # and clear out the memory from the big object
        lutCat = None

        fgcmLut = fgcm.FgcmLUT(lutIndexVals, lutFlat, lutDerivFlat, lutStd)

        # and clear out the memory of the big created objects
        lutFlat = None
        lutDerivFlat = None

        # next we need the exposure/visit information
        visitCat = butler.get('fgcmVisitCatalog')

        fgcmExpInfo = np.zeros(len(visitCat), dtype=[('VISIT', 'i8'),
                                                     ('MJD', 'f8'),
                                                     ('EXPTIME', 'f8'),
                                                     ('SEEING', 'f8'),
                                                     ('DEEPFLAG', 'i2'),
                                                     ('TELHA', 'f8'),
                                                     ('TELRA', 'f8'),
                                                     ('TELDEC', 'f8'),
                                                     ('PMB', 'f8'),
                                                     ('FILTERNAME', 'a2')])
        fgcmExpInfo['VISIT'][:] = visitCat['visit']
        fgcmExpInfo['MJD'][:] = visitCat['mjd']
        fgcmExpInfo['EXPTIME'][:] = visitCat['exptime']
        fgcmExpInfo['SEEING'][:] = visitCat['fwhm']
        fgcmExpInfo['DEEPFLAG'][:] = visitCat['deepflag']
        fgcmExpInfo['TELHA'][:] = visitCat['telha']
        fgcmExpInfo['TELRA'][:] = visitCat['telra']
        fgcmExpInfo['TELDEC'][:] = visitCat['teldec']
        fgcmExpInfo['PMB'][:] = visitCat['pmb']
        # Note that we have to go through asAstropy() to get a string
        #  array out of an afwTable
        fgcmExpInfo['FILTERNAME'][:] = visitCat.asAstropy()['filtername']

        # and we need to know the ccd offsets from the camera geometry
        ccdOffsets = np.zeros(lutIndexVals['NCCD'], dtype=[('CCDNUM', 'i4'),
                                                           ('DELTA_RA', 'f8'),
                                                           ('DELTA_DEC', 'f8'),
                                                           ('RA_SIZE', 'f8'),
                                                           ('DEC_SIZE', 'f8')])

        camera = butler.get('camera')

        extent = afwGeom.Extent2D(self.config.pixelScale, self.config.pixelScale)

        for i, detector in enumerate(camera):
            # new version, using proper rotations
            #  but I worry this only works with HSC, as there's a unit inconsistency

            camPoint = detector.getCenter(afwCameraGeom.PIXELS)
            bbox = detector.getBBox()
            orient = detector.getOrientation()

            ccdOffsets['CCDNUM'][i] = detector.getId()

            xform = orient.makePixelFpTransform(extent)
            pointXform = xform.applyForward(camPoint.getPoint())
            # this requires a pixelScale
            ccdOffsets['DELTA_RA'][i] = pointXform.getX() * self.config.pixelScale / 3600.0
            ccdOffsets['DELTA_DEC'][i] = pointXform.getY() * self.config.pixelScale / 3600.0

            # but this does not (for the delta)
            boxXform = xform.applyForward(afwGeom.Point2D(bbox.getMaxX(), bbox.getMaxY()))
            ccdOffsets['RA_SIZE'][i] = 2. * np.abs(boxXform.getX() -
                                                   pointXform.getX()) / 3600.0
            ccdOffsets['DEC_SIZE'][i] = 2. * np.abs(boxXform.getY() -
                                                    pointXform.getY()) / 3600.0

            # old version below
            # point = detector.getCenter(afwCameraGeom.FOCAL_PLANE)
            # bbox = detector.getBBox()

        noFitsDict = {'lutIndex': lutIndexVals,
                      'lutStd': lutStd,
                      'expInfo': fgcmExpInfo,
                      'ccdOffsets': ccdOffsets}

        # set up the fitter object
        fgcmFitCycle = fgcm.FgcmFitCycle(configDict, useFits=False,
                                         noFitsDict=noFitsDict)

        # create the parameter object
        if (fgcmFitCycle.initialCycle):
            # cycle = 0, initial cycle
            fgcmPars = fgcm.FgcmParameters.newParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                             fgcmLut,
                                                             fgcmExpInfo)
        else:
            # note that we already checked that this is available
            parCat = butler.get('fgcmFitParameters', fgcmcycle=self.config.cycleNumber-1)

            parLutFilterNames = np.array(parCat[0]['lutfilternames'].split(','))
            parFitBands = np.array(parCat[0]['fitbands'].split(','))
            parExtraBands = np.array(parCat[0]['extrabands'].split(','))

            # FIXME: check that these are the same as in the config, to be sure

            inParInfo = np.zeros(1, dtype=[('NCCD', 'i4'),
                                           ('LUTFILTERNAMES', parLutFilterNames.dtype.str,
                                            parLutFilterNames.size),
                                           ('FITBANDS', parFitBands.dtype.str, parFitBands.size),
                                           ('EXTRABANDS', parExtraBands.dtype.str, parExtraBands.size),
                                           ('TAUUNIT', 'f8'),
                                           ('TAUPERSLOPEUNIT', 'f8'),
                                           ('ALPHAUNIT', 'f8'),
                                           ('PWVUNIT', 'f8'),
                                           ('PWVPERSLOPEUNIT', 'f8'),
                                           ('O3UNIT', 'f8'),
                                           ('QESYSUNIT', 'f8'),
                                           ('QESYSSLOPEUNIT', 'f8'),
                                           ('HASEXTERNALPWV', 'i2'),
                                           ('HASEXTERNALTAU', 'i2')])
            inParInfo['NCCD'] = parCat['nccd']
            inParInfo['LUTFILTERNAMES'][:] = parLutFilterNames
            inParInfo['FITBANDS'][:] = parFitBands
            inParInfo['EXTRABANDS'][:] = parExtraBands
            inParInfo['TAUUNIT'] = parCat['tauunit']
            inParInfo['TAUPERSLOPEUNIT'] = parCat['tauperslopeunit']
            inParInfo['ALPHAUNIT'] = parCat['alphaunit']
            inParInfo['PWVUNIT'] = parCat['pwvunit']
            inParInfo['PWVPERSLOPEUNIT'] = parCat['pwvperslopeunit']
            inParInfo['O3UNIT'] = parCat['o3unit']
            inParInfo['QESYSUNIT'] = parCat['qesysunit']
            inParInfo['QESYSSLOPEUNIT'] = parCat['qesysslopeunit']
            inParInfo['HASEXTERNALPWV'] = parCat['hasexternalpwv']
            inParInfo['HASEXTERNALTAU'] = parCat['hasexternaltau']

            inParams = np.zeros(1, dtype=[('PARALPHA', 'f8', parCat['paralpha'].size),
                                          ('PARO3', 'f8', parCat['paro3'].size),
                                          ('PARTAUINTERCEPT', 'f8',
                                           parCat['partauintercept'].size),
                                          ('PARTAUPERSLOPE', 'f8',
                                           parCat['partauperslope'].size),
                                          ('PARPWVINTERCEPT', 'f8',
                                           parCat['parpwvintercept'].size),
                                          ('PARPWVPERSLOPE', 'f8',
                                           parCat['parpwvperslope'].size),
                                          ('PARQESYSINTERCEPT', 'f8',
                                           parCat['parqesysintercept'].size),
                                          ('PARQESYSSLOPE', 'f8',
                                           parCat['parqesysslope'].size),
                                          ('COMPAPERCORRPIVOT', 'f8',
                                           parCat['compapercorrpivot'].size),
                                          ('COMPAPERCORRSLOPE', 'f8',
                                           parCat['compapercorrslope'].size),
                                          ('COMPAPERCORRSLOPEERR', 'f8',
                                           parCat['compapercorrslopeerr'].size),
                                          ('COMPAPERCORRRANGE', 'f8',
                                           parCat['compapercorrrange'].size),
                                          ('COMPEXPGRAY', 'f8',
                                           parCat['compexpgray'].size),
                                          ('COMPVARGRAY', 'f8',
                                           parCat['compvargray'].size),
                                          ('COMPNGOODSTARPEREXP', 'i4',
                                           parCat['compngoodstarperexp'].size),
                                          ('COMPSIGFGCM', 'f8',
                                           parCat['compsigfgcm'].size)])

            inParams['PARALPHA'][:] = parCat['paralpha'][0, :]
            inParams['PARO3'][:] = parCat['paro3'][0, :]
            inParams['PARTAUINTERCEPT'][:] = parCat['partauintercept'][0, :]
            inParams['PARTAUPERSLOPE'][:] = parCat['partauperslope'][0, :]
            inParams['PARPWVINTERCEPT'][:] = parCat['parpwvintercept'][0, :]
            inParams['PARQESYSINTERCEPT'][:] = parCat['parqesysintercept'][0, :]
            inParams['PARQESYSSLOPE'][:] = parCat['parqesysslope'][0, :]
            inParams['COMPAPERCORRPIVOT'][:] = parCat['compapercorrpivot'][0, :]
            inParams['COMPAPERCORRSLOPE'][:] = parCat['compapercorrslope'][0, :]
            inParams['COMPAPERCORRSLOPEERR'][:] = parCat['compapercorrslopeerr'][0, :]
            inParams['COMPAPERCORRRANGE'][:] = parCat['compapercorrrange'][0, :]
            inParams['COMPEXPGRAY'][:] = parCat['compexpgray'][0, :]
            inParams['COMPVARGRAY'][:] = parCat['compvargray'][0, :]
            inParams['COMPNGOODSTARPEREXP'][:] = parCat['compngoodstarperexp'][0, :]
            inParams['COMPSIGFGCM'][:] = parCat['compsigfgcm'][0, :]

            inSuperStar = np.zeros(parCat['superstarsize'][0, :], dtype='f8')
            inSuperStar[:, :, :] = parCat['superstar'][0, :].reshape(inSuperStar.shape)

            fgcmPars = fgcm.FgcmParameters.loadParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                              fgcmExpInfo,
                                                              inParInfo,
                                                              inParams,
                                                              inSuperStar)

        lastCycle = configDict['cycleNumber'] - 1

        # set up the stars...
        fgcmStars = fgcm.FgcmStars(fgcmFitCycle.fgcmConfig)

        starObs = butler.get('fgcmStarObservations')
        starIds = butler.get('fgcmStarIds')
        starIndices = butler.get('fgcmStarIndices')

        # grab the flagged stars if available
        if butler.datasetExists('fgcmFlaggedStars', fgcmcycle=lastCycle):
            flaggedStars = butler.get('fgcmFlaggedStars', fgcmcycle=lastCycle)
            flagId = flaggedStars['objid'][:]
            flagFlag = flaggedStars['objflag'][:]
        else:
            flagId = None
            flagFlag = None

        # match star observations to visits
        visitIndex = np.searchsorted(fgcmExpInfo['VISIT'], starObs['visit'][starIndices['obsindex']])

        # note that we only need the star observations from specific indices

        fgcmStars.loadStars(fgcmPars,
                            starObs['visit'][starIndices['obsindex']],
                            starObs['ccd'][starIndices['obsindex']],
                            np.rad2deg(starObs['ra'][starIndices['obsindex']]),
                            np.rad2deg(starObs['dec'][starIndices['obsindex']]),
                            starObs['mag'][starIndices['obsindex']],
                            starObs['magerr'][starIndices['obsindex']],
                            fgcmExpInfo['FILTERNAME'][visitIndex],
                            starIds['fgcm_id'][:],
                            starIds['ra'][:],
                            starIds['dec'][:],
                            starIds['obsarrindex'][:],
                            starIds['nobs'][:],
                            obsX=starObs['x'][starIndices['obsindex']],
                            obsY=starObs['y'][starIndices['obsindex']],
                            flagID=flagId,
                            flagFlag=flagFlag,
                            computeNobs=True)

        # clear star memory
        starObs = None
        starIds = None
        starIndices = None
        flagId = None
        flagFlag = None
        flaggedStars = None

        # and set the bits in the cycle object
        fgcmFitCycle.setLUT(fgcmLut)
        fgcmFitCycle.setStars(fgcmStars)
        fgcmFitCycle.setPars(fgcmPars)

        # finish the setup
        fgcmFitCycle.finishSetup()

        # and run
        fgcmFitCycle.run()

        ##################
        # Persistance
        ##################

        # parameters
        parInfo, pars = fgcmFitCycle.fgcmPars.parsToArrays()

        parSchema = afwTable.Schema()

        comma = ','
        lutFilterNameString = comma.join(parInfo['LUTFILTERNAMES'][0])
        fitBandString = comma.join(parInfo['FITBANDS'][0])
        extraBandString = comma.join(parInfo['EXTRABANDS'][0])

        # parameter info section
        parSchema.addField('nccd', type=np.int32, doc='Number of CCDs')
        parSchema.addField('lutfilternames', type=str, doc='LUT Filter names in parameter file',
                           size=len(lutFilterNameString))
        parSchema.addField('fitbands', type=str, doc='Bands that were fit',
                           size=len(fitBandString))
        parSchema.addField('extrabands', type=str, doc='Bands that were not fit',
                           size=len(extraBandString))
        parSchema.addField('tauunit', type=np.float64, doc='Step units for AOD')
        parSchema.addField('tauperslopeunit', type=np.float64,
                           doc='Step units for AOD percent slope')
        parSchema.addField('alphaunit', type=np.float64, doc='Step units for alpha')
        parSchema.addField('pwvunit', type=np.float64, doc='Step units for pwv')
        parSchema.addField('pwvperslopeunit', type=np.float64,
                           doc='Step units for PWV percent slope')
        parSchema.addField('o3unit', type=np.float64, doc='Step units for O3')
        parSchema.addField('qesysunit', type=np.float64, doc='Step units for mirror gray')
        parSchema.addField('qesysslopeunit', type=np.float64, doc='Step units for mirror gray slope')
        parSchema.addField('hasexternalpwv', type=np.int32, doc='Parameters fit using external pwv')
        parSchema.addField('hasexternaltau', type=np.int32, doc='Parameters fit using external tau')

        # parameter section
        parSchema.addField('paralpha', type='ArrayD', doc='Alpha parameter vector',
                           size=pars['PARALPHA'].size)
        parSchema.addField('paro3', type='ArrayD', doc='O3 parameter vector',
                           size=pars['PARO3'].size)
        parSchema.addField('partauintercept', type='ArrayD', doc='Tau intercept parameter vector',
                           size=pars['PARTAUINTERCEPT'].size)
        parSchema.addField('partauperslope', type='ArrayD', doc='Tau percent slope parameter vector',
                           size=pars['PARTAUPERSLOPE'].size)
        parSchema.addField('parpwvintercept', type='ArrayD', doc='PWV intercept parameter vector',
                           size=pars['PARPWVINTERCEPT'].size)
        parSchema.addField('parpwvperslope', type='ArrayD', doc='PWV percent slope parameter vector',
                           size=pars['PARPWVPERSLOPE'].size)
        parSchema.addField('parqesysintercept', type='ArrayD', doc='Mirror gray intercept parameter vector',
                           size=pars['PARQESYSINTERCEPT'].size)
        parSchema.addField('parqesysslope', type='ArrayD', doc='Mirror gray slope parameter vector',
                           size=pars['PARQESYSSLOPE'].size)
        parSchema.addField('compapercorrpivot', type='ArrayD', doc='Aperture correction pivot',
                           size=pars['COMPAPERCORRPIVOT'].size)
        parSchema.addField('compapercorrslope', type='ArrayD', doc='Aperture correction slope',
                           size=pars['COMPAPERCORRSLOPE'].size)
        parSchema.addField('compapercorrslopeerr', type='ArrayD', doc='Aperture correction slope error',
                           size=pars['COMPAPERCORRSLOPEERR'].size)
        parSchema.addField('compapercorrrange', type='ArrayD', doc='Aperture correction range',
                           size=pars['COMPAPERCORRRANGE'].size)
        parSchema.addField('compexpgray', type='ArrayD', doc='Computed exposure gray',
                           size=pars['COMPEXPGRAY'].size)
        parSchema.addField('compvargray', type='ArrayD', doc='Computed exposure variance',
                           size=pars['COMPVARGRAY'].size)
        parSchema.addField('compngoodstarperexp', type='ArrayI',
                           doc='Computed number of good stars per exposure',
                           size=pars['COMPNGOODSTARPEREXP'].size)
        parSchema.addField('compsigfgcm', type='ArrayD', doc='Computed sigma_fgcm',
                           size=pars['COMPSIGFGCM'].size)

        # superstarflat section
        parSchema.addField('superstarsize', type='ArrayI', doc='Superstar matrix size',
                           size=3)
        parSchema.addField('superstar', type='ArrayD', doc='Superstar matrix (flattened)',
                           size=fgcmFitCycle.fgcmPars.parSuperStarFlat.size)

        parCat = afwTable.BaseCatalog(parSchema)
        parCat.table.preallocate(1)

        rec = parCat.addNew()

        # info section
        rec['nccd'] = parInfo['NCCD']
        rec['lutfilternames'] = lutFilterNameString
        rec['fitbands'] = fitBandString
        rec['extrabands'] = extraBandString
        rec['tauunit'] = parInfo['TAUUNIT']
        rec['tauperslopeunit'] = parInfo['TAUPERSLOPEUNIT']
        rec['alphaunit'] = parInfo['ALPHAUNIT']
        rec['pwvunit'] = parInfo['PWVUNIT']
        rec['pwvperslopeunit'] = parInfo['PWVPERSLOPEUNIT']
        rec['o3unit'] = parInfo['O3UNIT']
        rec['qesysunit'] = parInfo['QESYSUNIT']
        rec['qesysslopeunit'] = parInfo['QESYSSLOPEUNIT']
        # note these are not currently supported here.
        rec['hasexternalpwv'] = 0
        rec['hasexternaltau'] = 0

        # parameter section
        rec['paralpha'][:] = pars['PARALPHA'][0, :]
        rec['paro3'][:] = pars['PARO3'][0, :]
        rec['partauintercept'][:] = pars['PARTAUINTERCEPT'][0, :]
        rec['partauperslope'][:] = pars['PARTAUPERSLOPE'][0, :]
        rec['parpwvintercept'][:] = pars['PARPWVINTERCEPT'][0, :]
        rec['parpwvperslope'][:] = pars['PARPWVPERSLOPE'][0, :]
        rec['parqesysintercept'][:] = pars['PARQESYSINTERCEPT'][0, :]
        rec['parqesysslope'][:] = pars['PARQESYSSLOPE'][0, :]
        rec['compapercorrpivot'][:] = pars['COMPAPERCORRPIVOT'][0, :]
        rec['compapercorrslope'][:] = pars['COMPAPERCORRSLOPE'][0, :]
        rec['compapercorrslopeerr'][:] = pars['COMPAPERCORRSLOPEERR'][0, :]
        rec['compapercorrrange'][:] = pars['COMPAPERCORRRANGE'][0, :]
        rec['compexpgray'][:] = pars['COMPEXPGRAY'][0, :]
        rec['compvargray'][:] = pars['COMPVARGRAY'][0, :]
        rec['compngoodstarperexp'][:] = pars['COMPNGOODSTARPEREXP'][0, :]
        rec['compsigfgcm'][:] = pars['COMPSIGFGCM'][0, :]

        # superstar section
        rec['superstarsize'][:] = fgcmFitCycle.fgcmPars.parSuperStarFlat.shape
        rec['superstar'][:] = fgcmFitCycle.fgcmPars.parSuperStarFlat.flatten()

        butler.put(parCat, 'fgcmFitParameters', fgcmcycle=self.config.cycleNumber)

        # Save the flagged stars
        flagStarSchema = afwTable.Schema()

        flagStarSchema.addField('objid', type=np.int32, doc='FGCM object id')
        flagStarSchema.addField('objflag', type=np.int32, doc='FGCM object flag')

        flagStarCat = afwTable.BaseCatalog(flagStarSchema)
        flagStarStruct = fgcmFitCycle.fgcmStars.getFlagStarIndices()
        flagStarCat.table.preallocate(flagStarStruct.size)
        for i in xrange(flagStarStruct.size):
            rec = flagStarCat.addNew()

        flagStarCat = flagStarCat.copy(deep=True)

        flagStarCat['objid'][:] = flagStarStruct['OBJID']
        flagStarCat['objflag'][:] = flagStarStruct['OBJFLAG']

        butler.put(flagStarCat, 'fgcmFlaggedStars', fgcmcycle=self.config.cycleNumber)

        # Save zeropoints
        zptSchema = afwTable.Schema()

        zptSchema.addField('visit', type=np.int32, doc='Visit number')
        zptSchema.addField('ccd', type=np.int32, doc='CCD number')
        zptSchema.addField('fgcmflag', type=np.int32, doc='FGCM flag value')
        zptSchema.addField('fgcmzpt', type=np.float32, doc='FGCM zeropoint')
        zptSchema.addField('fgcmzpterr', type=np.float32,
                           doc='Error on zeropoint, estimated from repeatability + number of obs')
        zptSchema.addField('fgcmi0', type=np.float32, doc='Integral of the passband')
        zptSchema.addField('fgcmi10', type=np.float32, doc='Normalized chromatic integral')
        zptSchema.addField('fgcmr0', type=np.float32,
                           doc='Retrieved i0 integral, estimated from stars (only for flag 1)')
        zptSchema.addField('fgcmr10', type=np.float32,
                           doc='Retrieved i10 integral, estimated from stars (only for flag 1)')
        zptSchema.addField('fgcmgry', type=np.float32,
                           doc='Estimated gray extinction relative to atmospheric solution; '
                           'only for flag <= 4')
        zptSchema.addField('fgcmzptvar', type=np.float32, doc='Variance of zeropoint over ccd')
        zptSchema.addField('fgcmtilings', type=np.float32,
                           doc='Number of photometric tilings used for solution for ccd')
        zptSchema.addField('fgcmfpgry', type=np.float32,
                           doc='Average gray extinction over the full focal plane '
                           '(same for all ccds in a visit)')
        zptSchema.addField('fgcmfpvar', type=np.float32,
                           doc='Variance of gray extinction over the full focal plane '
                           '(same for all ccds in a visit)')
        zptSchema.addField('fgcmdust', type=np.float32,
                           doc='Gray dust extinction from the primary/corrector'
                           'at the time of the exposure')
        zptSchema.addField('fgcmflat', type=np.float32, doc='Superstarflat illumination correction')
        zptSchema.addField('fgcmapercorr', type=np.float32, doc='Aperture correction estimated by fgcm')
        zptSchema.addField('exptime', type=np.float32, doc='Exposure time')
        zptSchema.addField('filtername', type=str, size=2, doc='Filter name')

        zptCat = afwTable.BaseCatalog(zptSchema)
        zptCat.table.preallocate(fgcmFitCycle.fgcmZpts.zpStruct.size)
        for filterName in fgcmFitCycle.fgcmZpts.zpStruct['FILTERNAME']:
            rec = zptCat.addNew()
            rec['filtername'] = filterName

        zptCat = zptCat.copy(deep=True)

        zptCat['visit'][:] = fgcmFitCycle.fgcmZpts.zpStruct['VISIT']
        zptCat['ccd'][:] = fgcmFitCycle.fgcmZpts.zpStruct['CCD']
        zptCat['fgcmflag'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FLAG']
        zptCat['fgcmzpt'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_ZPT']
        zptCat['fgcmzpterr'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_ZPTERR']
        zptCat['fgcmi0'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_I0']
        zptCat['fgcmi10'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_I10']
        zptCat['fgcmr0'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_R0']
        zptCat['fgcmr10'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_R10']
        zptCat['fgcmgry'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_GRY']
        zptCat['fgcmzptvar'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_ZPTVAR']
        zptCat['fgcmtilings'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_TILINGS']
        zptCat['fgcmfpgry'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FPGRY']
        zptCat['fgcmfpvar'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FPVAR']
        zptCat['fgcmdust'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_DUST']
        zptCat['fgcmflat'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FLAT']
        zptCat['fgcmapercorr'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_APERCORR']
        zptCat['exptime'][:] = fgcmFitCycle.fgcmZpts.zpStruct['EXPTIME']

        butler.put(zptCat, 'fgcmZeropoints', fgcmcycle=self.config.cycleNumber)

        # Save atmosphere values

        atmSchema = afwTable.Schema()

        atmSchema.addField('visit', type=np.int32, doc='Visit number')
        atmSchema.addField('pmb', type=np.float64, doc='Barometric pressure (mb)')
        atmSchema.addField('pwv', type=np.float64, doc='Water vapor (mm)')
        atmSchema.addField('tau', type=np.float64, doc='Aerosol optical depth')
        atmSchema.addField('alpha', type=np.float64, doc='Aerosol slope')
        atmSchema.addField('o3', type=np.float64, doc='Ozone (dobson)')
        atmSchema.addField('seczenith', type=np.float64, doc='Secant(zenith) (~ airmass)')

        atmCat = afwTable.BaseCatalog(atmSchema)
        atmCat.table.preallocate(fgcmFitCycle.fgcmZpts.atmStruct.size)
        for i in xrange(fgcmFitCycle.fgcmZpts.atmStruct.size):
            rec = atmCat.addNew()

        atmCat = atmCat.copy(deep=True)

        atmCat['visit'][:] = fgcmFitCycle.fgcmZpts.atmStruct['VISIT']
        atmCat['pmb'][:] = fgcmFitCycle.fgcmZpts.atmStruct['PMB']
        atmCat['pwv'][:] = fgcmFitCycle.fgcmZpts.atmStruct['PWV']
        atmCat['tau'][:] = fgcmFitCycle.fgcmZpts.atmStruct['TAU']
        atmCat['alpha'][:] = fgcmFitCycle.fgcmZpts.atmStruct['ALPHA']
        atmCat['o3'][:] = fgcmFitCycle.fgcmZpts.atmStruct['O3']
        atmCat['seczenith'][:] = fgcmFitCycle.fgcmZpts.atmStruct['SECZENITH']

        butler.put(atmCat, 'fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        # Output the config for the next cycle
        # We need to make a copy since the input one has been frozen

        outConfig = FgcmFitCycleConfig()
        outConfig.update(**self.config.toDict())

        outConfig.cycleNumber += 1
        outConfig.precomputeSuperStarInitialCycle = False
        outConfig.freezeStdAtmosphere = False
        configFileName = '%s_cycle%02d_config.py' % (outConfig.outfileBase,
                                                     outConfig.cycleNumber)
        outConfig.save(configFileName)

        self.log.info("Saved config for next cycle to %s" % (configFileName))
        self.log.info("Be sure to look at:")
        self.log.info("   config.freezeStdAtmosphere")
        self.log.info("   config.expGrayPhotometricCut")

        # tear down and clear memory
