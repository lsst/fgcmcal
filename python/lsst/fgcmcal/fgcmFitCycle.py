# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from past.builtins import xrange

import matplotlib
matplotlib.use("Agg")  # noqa

import sys
import traceback

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.geom

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
    requiredFlag = pexConfig.ListField(
        doc="Flag for bands to require to be a calibration star",
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
    superStarSubCcd = pexConfig.Field(
        doc="Compute superstar flat on sub-ccd scale",
        dtype=bool,
        default=True,
    )
    superStarSubCcdChebyshevOrder = pexConfig.Field(
        doc="Order of chebyshev polynomials for sub-ccd superstar fit",
        dtype=int,
        default=1,
    )
    superStarSigmaClip = pexConfig.Field(
        doc="Number of sigma to clip outliers when selecting for superstar flats",
        dtype=float,
        default=5.0,
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
    colorSplitIndices = pexConfig.ListField(
        doc="Band indices to use to split stars by color",
        dtype=int,
        default=None,
    )
    modelMagErrors = pexConfig.Field(
        doc="Should FGCM model the magnitude errors from sky/fwhm? (False means trust inputs)",
        dtype=bool,
        default=True,
    )
    outputStandards = pexConfig.Field(
        doc="Output standard stars? (Usually only for final iteration)",
        dtype=bool,
        default=False,
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

    # This overrides the pipe_base config saving which would fail because it
    # requires the %(fgcmcycle)d dataId key.
    # def precall(self, parsedCmd):
    #    return True

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
            # Note that there's not much results to return (empty list)
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
    def runDataRef(self, butler):
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

    def writeConfig(self, butler, clobber=False, doBackup=True):
        """Write the configuration used for processing the data, or check that an existing
        one is equal to the new one if present.  This is an override of the regular
        version from pipe_base that knows about fgcmcycle.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
            Data butler used to write the config. The config is written to dataset type
            `CmdLineTask._getConfigName`.
        clobber : `bool`, optional
            A boolean flag that controls what happens if a config already has been saved:
            - `True`: overwrite or rename the existing config, depending on ``doBackup``.
            - `False`: raise `TaskError` if this config does not match the existing config.
        doBackup : bool, optional
            Set to `True` to backup the config files if clobbering.
        """
        configName = self._getConfigName()
        if configName is None:
            return
        if clobber:
            butler.put(self.config, configName, doBackup=doBackup, fgcmcycle=self.config.cycleNumber)
        elif butler.datasetExists(configName, write=True, fgcmcycle=self.config.cycleNumber):
            # this may be subject to a race condition; see #2789
            try:
                oldConfig = butler.get(configName, immediate=True, fgcmcycle=self.config.cycleNumber)
            except Exception as exc:
                raise type(exc)("Unable to read stored config file %s (%s); consider using --clobber-config" %
                                (configName, exc))

            def logConfigMismatch(msg):
                self.log.fatal("Comparing configuration: %s", msg)

            if not self.config.compare(oldConfig, shortcut=False, output=logConfigMismatch):
                raise pipeBase.TaskError(
                    ("Config does not match existing task config %r on disk; tasks configurations " +
                     "must be consistent within the same output repo (override with --clobber-config)") %
                    (configName,))
        else:
            butler.put(self.config, configName, fgcmcycle=self.config.cycleNumber)

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

        fitFlag = np.array(self.config.fitFlag, dtype=np.bool)
        requiredFlag = np.array(self.config.requiredFlag, dtype=np.bool)

        fitBands = [b for i, b in enumerate(self.config.bands) if fitFlag[i]]
        notFitBands = [b for i, b in enumerate(self.config.bands) if not fitFlag[i]]
        requiredBands = [b for i, b in enumerate(self.config.bands) if requiredFlag[i]]

        camera = butler.get('camera')

        # process the starColorCuts
        starColorCutList = []
        for ccut in self.config.starColorCuts:
            parts = ccut.split(',')
            starColorCutList.append([parts[0], parts[1], float(parts[2]), float(parts[3])])

        if self.config.maxIter == 0:
            resetParameters = False
        else:
            resetParameters = True

        # Mirror area in cm**2
        mirrorArea = np.pi*(camera.telescopeDiameter*100./2.)**2.

        # create a configuration dictionary for fgcmFitCycle
        configDict = {'outfileBase': self.config.outfileBase,
                      'logger': self.log,
                      'exposureFile': None,
                      'obsFile': None,
                      'indexFile': None,
                      'lutFile': None,
                      'mirrorArea': mirrorArea,
                      'cameraGain': self.config.cameraGain,
                      'ccdStartIndex': camera[0].getId(),
                      'expField': 'VISIT',
                      'ccdField': 'CCD',
                      'seeingField': 'DELTA_APER',
                      'fwhmField': 'PSFSIGMA',
                      'skyBrightnessField': 'SKYBACKGROUND',
                      'deepFlag': 'DEEPFLAG',  # unused
                      'bands': list(self.config.bands),
                      'fitBands': list(fitBands),
                      'notFitBands': list(notFitBands),
                      'requiredBands': list(requiredBands),
                      'filterToBand': dict(self.config.filterToBand),
                      'logLevel': 'INFO',  # FIXME
                      'nCore': self.config.nCore,
                      'nStarPerRun': self.config.nStarPerRun,
                      'nExpPerRun': self.config.nExpPerRun,
                      'reserveFraction': self.config.reserveFraction,
                      'freezeStdAtmosphere': self.config.freezeStdAtmosphere,
                      'precomputeSuperStarInitialCycle': self.config.precomputeSuperStarInitialCycle,
                      'superStarSubCCD': self.config.superStarSubCcd,
                      'superStarSubCCDChebyshevOrder': self.config.superStarSubCcdChebyshevOrder,
                      'superStarSigmaClip': self.config.superStarSigmaClip,
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
                      'sedFudgeFactors': np.array(self.config.sedFudgeFactors),
                      'colorSplitIndices': np.array(self.config.colorSplitIndices),
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
                      'modelMagErrors': self.config.modelMagErrors,
                      'printOnly': False,
                      'outputStars': False,
                      'clobber': True,
                      'useSedLUT': False,
                      'resetParameters': resetParameters}

        fgcmLut, lutIndexVals, lutStd = self._loadFgcmLut(butler,
                                                          filterToBand=self.config.filterToBand)

        # next we need the exposure/visit information

        fgcmExpInfo = self._loadVisitCatalog(butler)

        ccdOffsets = self._loadCcdOffsets(butler)

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
            parNotFitBands = np.array(parCat[0]['notfitbands'].split(','))

            # FIXME: check that these are the same as in the config, to be sure

            inParInfo = np.zeros(1, dtype=[('NCCD', 'i4'),
                                           ('LUTFILTERNAMES', parLutFilterNames.dtype.str,
                                            parLutFilterNames.size),
                                           ('FITBANDS', parFitBands.dtype.str, parFitBands.size),
                                           ('NOTFITBANDS', parNotFitBands.dtype.str, parNotFitBands.size),
                                           ('LNTAUUNIT', 'f8'),
                                           ('LNTAUSLOPEUNIT', 'f8'),
                                           ('ALPHAUNIT', 'f8'),
                                           ('PWVUNIT', 'f8'),
                                           ('PWVPERSLOPEUNIT', 'f8'),
                                           ('PWVGLOBALUNIT', 'f8'),
                                           ('O3UNIT', 'f8'),
                                           ('QESYSUNIT', 'f8'),
                                           ('QESYSSLOPEUNIT', 'f8'),

                                           ('HASEXTERNALPWV', 'i2'),
                                           ('HASEXTERNALTAU', 'i2')])
            inParInfo['NCCD'] = parCat['nccd']
            inParInfo['LUTFILTERNAMES'][:] = parLutFilterNames
            inParInfo['FITBANDS'][:] = parFitBands
            inParInfo['NOTFITBANDS'][:] = parNotFitBands
            inParInfo['LNTAUUNIT'] = parCat['lntauunit']
            inParInfo['LNTAUSLOPEUNIT'] = parCat['lntauslopeunit']
            inParInfo['ALPHAUNIT'] = parCat['alphaunit']
            inParInfo['PWVUNIT'] = parCat['pwvunit']
            inParInfo['PWVPERSLOPEUNIT'] = parCat['pwvperslopeunit']
            inParInfo['PWVGLOBALUNIT'] = parCat['pwvglobalunit']
            inParInfo['O3UNIT'] = parCat['o3unit']
            inParInfo['QESYSUNIT'] = parCat['qesysunit']
            inParInfo['QESYSSLOPEUNIT'] = parCat['qesysslopeunit']
            inParInfo['HASEXTERNALPWV'] = parCat['hasexternalpwv']
            inParInfo['HASEXTERNALTAU'] = parCat['hasexternaltau']

            inParams = np.zeros(1, dtype=[('PARALPHA', 'f8', parCat['paralpha'].size),
                                          ('PARO3', 'f8', parCat['paro3'].size),
                                          ('PARLNTAUINTERCEPT', 'f8',
                                           parCat['parlntauintercept'].size),
                                          ('PARLNTAUSLOPE', 'f8',
                                           parCat['parlntauslope'].size),
                                          ('PARPWVINTERCEPT', 'f8',
                                           parCat['parpwvintercept'].size),
                                          ('PARPWVPERSLOPE', 'f8',
                                           parCat['parpwvperslope'].size),
                                          ('PARQESYSINTERCEPT', 'f8',
                                           parCat['parqesysintercept'].size),
                                          ('PARQESYSSLOPE', 'f8',
                                           parCat['parqesysslope'].size),
                                          ('PARRETRIEVEDPWVSCALE', 'f8'),
                                          ('PARRETRIEVEDPWVOFFSET', 'f8'),
                                          ('PARRETRIEVEDPWVNIGHTLYOFFSET', 'f8',
                                           parCat['parretrievedpwvnightlyoffset'].size),
                                          ('COMPAPERCORRPIVOT', 'f8',
                                           parCat['compapercorrpivot'].size),
                                          ('COMPAPERCORRSLOPE', 'f8',
                                           parCat['compapercorrslope'].size),
                                          ('COMPAPERCORRSLOPEERR', 'f8',
                                           parCat['compapercorrslopeerr'].size),
                                          ('COMPAPERCORRRANGE', 'f8',
                                           parCat['compapercorrrange'].size),
                                          ('COMPMODELERREXPTIMEPIVOT', 'f8',
                                           parCat['compmodelerrexptimepivot'].size),
                                          ('COMPMODELERRFWHMPIVOT', 'f8',
                                           parCat['compmodelerrfwhmpivot'].size),
                                          ('COMPMODELERRSKYPIVOT', 'f8',
                                           parCat['compmodelerrskypivot'].size),
                                          ('COMPMODELERRPARS', 'f8',
                                           parCat['compmodelerrpars'].size),
                                          ('COMPEXPGRAY', 'f8',
                                           parCat['compexpgray'].size),
                                          ('COMPVARGRAY', 'f8',
                                           parCat['compvargray'].size),
                                          ('COMPNGOODSTARPEREXP', 'i4',
                                           parCat['compngoodstarperexp'].size),
                                          ('COMPSIGFGCM', 'f8',
                                           parCat['compsigfgcm'].size),
                                          ('COMPRETRIEVEDPWV', 'f8',
                                           parCat['compretrievedpwv'].size),
                                          ('COMPRETRIEVEDPWVRAW', 'f8',
                                           parCat['compretrievedpwvraw'].size),
                                          ('COMPRETRIEVEDPWVFLAG', 'i2',
                                           parCat['compretrievedpwvflag'].size),
                                          ('COMPRETRIEVEDTAUNIGHT', 'f8',
                                           parCat['compretrievedtaunight'].size)])

            inParams['PARALPHA'][:] = parCat['paralpha'][0, :]
            inParams['PARO3'][:] = parCat['paro3'][0, :]
            inParams['PARLNTAUINTERCEPT'][:] = parCat['parlntauintercept'][0, :]
            inParams['PARLNTAUSLOPE'][:] = parCat['parlntauslope'][0, :]
            inParams['PARPWVINTERCEPT'][:] = parCat['parpwvintercept'][0, :]
            inParams['PARQESYSINTERCEPT'][:] = parCat['parqesysintercept'][0, :]
            inParams['PARQESYSSLOPE'][:] = parCat['parqesysslope'][0, :]
            inParams['PARRETRIEVEDPWVSCALE'] = parCat['parretrievedpwvscale']
            inParams['PARRETRIEVEDPWVOFFSET'] = parCat['parretrievedpwvoffset']
            inParams['PARRETRIEVEDPWVNIGHTLYOFFSET'][:] = parCat['parretrievedpwvnightlyoffset'][0, :]
            inParams['COMPAPERCORRPIVOT'][:] = parCat['compapercorrpivot'][0, :]
            inParams['COMPAPERCORRSLOPE'][:] = parCat['compapercorrslope'][0, :]
            inParams['COMPAPERCORRSLOPEERR'][:] = parCat['compapercorrslopeerr'][0, :]
            inParams['COMPAPERCORRRANGE'][:] = parCat['compapercorrrange'][0, :]
            inParams['COMPMODELERREXPTIMEPIVOT'][:] = parCat['compmodelerrexptimepivot'][0, :]
            inParams['COMPMODELERRFWHMPIVOT'][:] = parCat['compmodelerrfwhmpivot'][0, :]
            inParams['COMPMODELERRSKYPIVOT'][:] = parCat['compmodelerrskypivot'][0, :]
            inParams['COMPMODELERRPARS'][:] = parCat['compmodelerrpars'][0, :]
            inParams['COMPEXPGRAY'][:] = parCat['compexpgray'][0, :]
            inParams['COMPVARGRAY'][:] = parCat['compvargray'][0, :]
            inParams['COMPNGOODSTARPEREXP'][:] = parCat['compngoodstarperexp'][0, :]
            inParams['COMPSIGFGCM'][:] = parCat['compsigfgcm'][0, :]
            inParams['COMPRETRIEVEDPWV'][:] = parCat['compretrievedpwv'][0, :]
            inParams['COMPRETRIEVEDPWVRAW'][:] = parCat['compretrievedpwvraw'][0, :]
            inParams['COMPRETRIEVEDPWVFLAG'][:] = parCat['compretrievedpwvflag'][0, :]
            inParams['COMPRETRIEVEDTAUNIGHT'][:] = parCat['compretrievedtaunight'][0, :]

            inSuperStar = np.zeros(parCat['superstarsize'][0, :], dtype='f8')
            inSuperStar[:, :, :, :] = parCat['superstar'][0, :].reshape(inSuperStar.shape)

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
        lutFilterNameString = comma.join([n.decode('utf-8')
                                          for n in parInfo['LUTFILTERNAMES'][0]])
        fitBandString = comma.join([n.decode('utf-8')
                                    for n in parInfo['FITBANDS'][0]])
        notFitBandString = comma.join([n.decode('utf-8')
                                       for n in parInfo['NOTFITBANDS'][0]])

        # parameter info section
        parSchema.addField('nccd', type=np.int32, doc='Number of CCDs')
        parSchema.addField('lutfilternames', type=str, doc='LUT Filter names in parameter file',
                           size=len(lutFilterNameString))
        parSchema.addField('fitbands', type=str, doc='Bands that were fit',
                           size=len(fitBandString))
        parSchema.addField('notfitbands', type=str, doc='Bands that were not fit',
                           size=len(notFitBandString))
        parSchema.addField('lntauunit', type=np.float64, doc='Step units for ln(AOD)')
        parSchema.addField('lntauslopeunit', type=np.float64,
                           doc='Step units for ln(AOD) slope')
        parSchema.addField('alphaunit', type=np.float64, doc='Step units for alpha')
        parSchema.addField('pwvunit', type=np.float64, doc='Step units for pwv')
        parSchema.addField('pwvperslopeunit', type=np.float64,
                           doc='Step units for PWV percent slope')
        parSchema.addField('pwvglobalunit', type=np.float64,
                           doc='Step units for global PWV parameters')
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
        parSchema.addField('parlntauintercept', type='ArrayD',
                           doc='ln(Tau) intercept parameter vector',
                           size=pars['PARLNTAUINTERCEPT'].size)
        parSchema.addField('parlntauslope', type='ArrayD',
                           doc='ln(Tau) slope parameter vector',
                           size=pars['PARLNTAUSLOPE'].size)
        parSchema.addField('parpwvintercept', type='ArrayD', doc='PWV intercept parameter vector',
                           size=pars['PARPWVINTERCEPT'].size)
        parSchema.addField('parpwvperslope', type='ArrayD', doc='PWV percent slope parameter vector',
                           size=pars['PARPWVPERSLOPE'].size)
        parSchema.addField('parqesysintercept', type='ArrayD', doc='Mirror gray intercept parameter vector',
                           size=pars['PARQESYSINTERCEPT'].size)
        parSchema.addField('parqesysslope', type='ArrayD', doc='Mirror gray slope parameter vector',
                           size=pars['PARQESYSSLOPE'].size)
        parSchema.addField('parretrievedpwvscale', type=np.float64,
                           doc='Global scale for retrieved PWV')
        parSchema.addField('parretrievedpwvoffset', type=np.float64,
                           doc='Global offset for retrieved PWV')
        parSchema.addField('parretrievedpwvnightlyoffset', type='ArrayD',
                           doc='Nightly offset for retrieved PWV',
                           size=pars['PARRETRIEVEDPWVNIGHTLYOFFSET'].size)
        parSchema.addField('compapercorrpivot', type='ArrayD', doc='Aperture correction pivot',
                           size=pars['COMPAPERCORRPIVOT'].size)
        parSchema.addField('compapercorrslope', type='ArrayD', doc='Aperture correction slope',
                           size=pars['COMPAPERCORRSLOPE'].size)
        parSchema.addField('compapercorrslopeerr', type='ArrayD', doc='Aperture correction slope error',
                           size=pars['COMPAPERCORRSLOPEERR'].size)
        parSchema.addField('compapercorrrange', type='ArrayD', doc='Aperture correction range',
                           size=pars['COMPAPERCORRRANGE'].size)
        parSchema.addField('compmodelerrexptimepivot', type='ArrayD', doc='Model error exptime pivot',
                           size=pars['COMPMODELERREXPTIMEPIVOT'].size)
        parSchema.addField('compmodelerrfwhmpivot', type='ArrayD', doc='Model error fwhm pivot',
                           size=pars['COMPMODELERRFWHMPIVOT'].size)
        parSchema.addField('compmodelerrskypivot', type='ArrayD', doc='Model error sky pivot',
                           size=pars['COMPMODELERRSKYPIVOT'].size)
        parSchema.addField('compmodelerrpars', type='ArrayD', doc='Model error parameters',
                           size=pars['COMPMODELERRPARS'].size)
        parSchema.addField('compexpgray', type='ArrayD', doc='Computed exposure gray',
                           size=pars['COMPEXPGRAY'].size)
        parSchema.addField('compvargray', type='ArrayD', doc='Computed exposure variance',
                           size=pars['COMPVARGRAY'].size)
        parSchema.addField('compngoodstarperexp', type='ArrayI',
                           doc='Computed number of good stars per exposure',
                           size=pars['COMPNGOODSTARPEREXP'].size)
        parSchema.addField('compsigfgcm', type='ArrayD', doc='Computed sigma_fgcm',
                           size=pars['COMPSIGFGCM'].size)
        parSchema.addField('compretrievedpwv', type='ArrayD', doc='Retrieved PWV (smoothed)',
                           size=pars['COMPRETRIEVEDPWV'].size)
        parSchema.addField('compretrievedpwvraw', type='ArrayD', doc='Retrieved PWV (raw)',
                           size=pars['COMPRETRIEVEDPWVRAW'].size)
        parSchema.addField('compretrievedpwvflag', type='ArrayI', doc='Retrieved PWV Flag',
                           size=pars['COMPRETRIEVEDPWVFLAG'].size)
        parSchema.addField('compretrievedtaunight', type='ArrayD', doc='Retrieved tau (per night)',
                           size=pars['COMPRETRIEVEDTAUNIGHT'].size)

        # superstarflat section
        parSchema.addField('superstarsize', type='ArrayI', doc='Superstar matrix size',
                           size=4)
        parSchema.addField('superstar', type='ArrayD', doc='Superstar matrix (flattened)',
                           size=fgcmFitCycle.fgcmPars.parSuperStarFlat.size)

        parCat = afwTable.BaseCatalog(parSchema)
        parCat.reserve(1)

        rec = parCat.addNew()

        # info section
        rec['nccd'] = parInfo['NCCD']
        rec['lutfilternames'] = lutFilterNameString
        rec['fitbands'] = fitBandString
        rec['notfitbands'] = notFitBandString
        rec['lntauunit'] = parInfo['LNTAUUNIT']
        rec['lntauslopeunit'] = parInfo['LNTAUSLOPEUNIT']
        rec['alphaunit'] = parInfo['ALPHAUNIT']
        rec['pwvunit'] = parInfo['PWVUNIT']
        rec['pwvperslopeunit'] = parInfo['PWVPERSLOPEUNIT']
        rec['pwvglobalunit'] = parInfo['PWVGLOBALUNIT']
        rec['o3unit'] = parInfo['O3UNIT']
        rec['qesysunit'] = parInfo['QESYSUNIT']
        rec['qesysslopeunit'] = parInfo['QESYSSLOPEUNIT']
        # note these are not currently supported here.
        rec['hasexternalpwv'] = 0
        rec['hasexternaltau'] = 0

        # parameter section

        scalarNames = ['parretrievedpwvscale', 'parretrievedpwvoffset']

        arrNames = ['paralpha', 'paro3', 'parlntauintercept', 'parlntauslope',
                    'parpwvintercept', 'parpwvperslope', 'parqesysintercept',
                    'parqesysslope', 'parretrievedpwvnightlyoffset', 'compapercorrpivot',
                    'compapercorrslope', 'compapercorrslopeerr', 'compapercorrrange',
                    'compmodelerrexptimepivot', 'compmodelerrfwhmpivot',
                    'compmodelerrskypivot', 'compmodelerrpars',
                    'compexpgray', 'compvargray', 'compngoodstarperexp', 'compsigfgcm',
                    'compretrievedpwv', 'compretrievedpwvraw', 'compretrievedpwvflag',
                    'compretrievedtaunight']

        for scalarName in scalarNames:
            rec[scalarName] = pars[scalarName.upper()]

        for arrName in arrNames:
            rec[arrName][:] = np.atleast_1d(pars[0][arrName.upper()])[:]

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
        flagStarCat.reserve(flagStarStruct.size)
        for i in xrange(flagStarStruct.size):
            rec = flagStarCat.addNew()

        if not flagStarCat.isContiguous():
            flagStarCat = flagStarCat.copy(deep=True)

        flagStarCat['objid'][:] = flagStarStruct['OBJID']
        flagStarCat['objflag'][:] = flagStarStruct['OBJFLAG']

        butler.put(flagStarCat, 'fgcmFlaggedStars', fgcmcycle=self.config.cycleNumber)

        # Save zeropoints
        zptSchema = afwTable.Schema()

        zptSchema.addField('visit', type=np.int32, doc='Visit number')
        zptSchema.addField('ccd', type=np.int32, doc='CCD number')
        zptSchema.addField('fgcmflag', type=np.int32, doc='FGCM flag value')
        zptSchema.addField('fgcmzpt', type=np.float32, doc='FGCM zeropoint (center of CCD)')
        zptSchema.addField('fgcmzpterr', type=np.float32,
                           doc='Error on zeropoint, estimated from repeatability + number of obs')
        if self.config.superStarSubCcd:
            zptSchema.addField('fgcmfzptcheb', type='ArrayD',
                               size=fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB'].shape[1],
                               doc='Chebyshev parameters (flattened) for zeropoint')
            zptSchema.addField('fgcmfzptchebxymax', type='ArrayD', size=2,
                               doc='maximum x/maximum y to scale to apply chebyshev parameters')
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
        zptCat.reserve(fgcmFitCycle.fgcmZpts.zpStruct.size)
        for filterName in fgcmFitCycle.fgcmZpts.zpStruct['FILTERNAME']:
            rec = zptCat.addNew()
            rec['filtername'] = filterName.decode('utf-8')

        if not zptCat.isContiguous():
            zptCat = zptCat.copy(deep=True)

        zptCat['visit'][:] = fgcmFitCycle.fgcmZpts.zpStruct['VISIT']
        zptCat['ccd'][:] = fgcmFitCycle.fgcmZpts.zpStruct['CCD']
        zptCat['fgcmflag'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FLAG']
        zptCat['fgcmzpt'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_ZPT']
        zptCat['fgcmzpterr'][:] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_ZPTERR']
        if self.config.superStarSubCcd:
            zptCat['fgcmfzptcheb'][:, :] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB']
            zptCat['fgcmfzptchebxymax'][:, :] = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB_XYMAX']
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
        atmCat.reserve(fgcmFitCycle.fgcmZpts.atmStruct.size)
        for i in xrange(fgcmFitCycle.fgcmZpts.atmStruct.size):
            rec = atmCat.addNew()

        if not atmCat.isContiguous():
            atmCat = atmCat.copy(deep=True)

        atmCat['visit'][:] = fgcmFitCycle.fgcmZpts.atmStruct['VISIT']
        atmCat['pmb'][:] = fgcmFitCycle.fgcmZpts.atmStruct['PMB']
        atmCat['pwv'][:] = fgcmFitCycle.fgcmZpts.atmStruct['PWV']
        atmCat['tau'][:] = fgcmFitCycle.fgcmZpts.atmStruct['TAU']
        atmCat['alpha'][:] = fgcmFitCycle.fgcmZpts.atmStruct['ALPHA']
        atmCat['o3'][:] = fgcmFitCycle.fgcmZpts.atmStruct['O3']
        atmCat['seczenith'][:] = fgcmFitCycle.fgcmZpts.atmStruct['SECZENITH']

        butler.put(atmCat, 'fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        if self.config.outputStandards:
            stdSchema = afwTable.SimpleTable.makeMinimalSchema()
            stdSchema.addField('ngood', type='ArrayI', doc='Number of good observations',
                               size=len(self.config.bands))
            stdSchema.addField('mag_std_noabs', type='ArrayF',
                               doc='Standard magnitude (no absolute calibration)',
                               size=len(self.config.bands))
            stdSchema.addField('magerr_std', type='ArrayF',
                               doc='Standard magnitude error',
                               size=len(self.config.bands))

            outCat = fgcmFitCycle.fgcmStars.retrieveStdStarCatalog(fgcmFitCycle.fgcmPars)
            stdCat = afwTable.SimpleCatalog(stdSchema)

            stdCat.reserve(outCat.size)
            for i in range(outCat.size):
                rec = stdCat.addNew()

            # Sometimes the reserve doesn't actually make a contiguous catalog (sigh)
            if not stdCat.isContiguous():
                stdCat = stdCat.copy(deep=True)

            stdCat['id'][:] = outCat['FGCM_ID']
            stdCat['coord_ra'][:] = outCat['RA'] * lsst.geom.degrees
            stdCat['coord_dec'][:] = outCat['DEC'] * lsst.geom.degrees
            stdCat['ngood'][:, :] = outCat['NGOOD'][:, :]
            stdCat['mag_std_noabs'][:, :] = outCat['MAG_STD'][:, :]
            stdCat['magerr_std'][:, :] = outCat['MAGERR_STD'][:, :]

            butler.put(stdCat, 'fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

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

        if self.config.maxIter == 0 and self.config.outputStandards:
            # We are done, there is no more warning
            self.log.info("Everything is in place to run fgcmOutputProducts.py")
        else:
            self.log.info("Saved config for next cycle to %s" % (configFileName))
            self.log.info("Be sure to look at:")
            self.log.info("   config.expGrayPhotometricCut")
            self.log.info("   config.expGrayHighCut")
            self.log.info("If you are satisfied with the fit, please set:")
            self.log.info("   config.maxIter = 0")
            self.log.info("   config.outputStandards = True")

    def _loadFgcmLut(self, butler, filterToBand=None):
        """
        """

        # set up the look-up-table
        lutCat = butler.get('fgcmLookUpTable')

        # first we need the lutIndexVals
        # dtype is set for py2/py3/fits/fgcm compatibility
        lutFilterNames = np.array(lutCat[0]['filternames'].split(','), dtype='a')
        lutStdFilterNames = np.array(lutCat[0]['stdfilternames'].split(','), dtype='a')

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

        fgcmLut = fgcm.FgcmLUT(lutIndexVals, lutFlat, lutDerivFlat, lutStd,
                               filterToBand=self.config.filterToBand)

        # and clear out the memory of the large temporary objects
        lutFlat = None
        lutDerivFlat = None

        return fgcmLut, lutIndexVals, lutStd

    def _loadVisitCatalog(self, butler):
        """
        """

        # next we need the exposure/visit information
        visitCat = butler.get('fgcmVisitCatalog')

        fgcmExpInfo = np.zeros(len(visitCat), dtype=[('VISIT', 'i8'),
                                                     ('MJD', 'f8'),
                                                     ('EXPTIME', 'f8'),
                                                     ('PSFSIGMA', 'f8'),
                                                     ('DELTA_APER', 'f8'),
                                                     ('SKYBACKGROUND', 'f8'),
                                                     ('DEEPFLAG', 'i2'),
                                                     ('TELHA', 'f8'),
                                                     ('TELRA', 'f8'),
                                                     ('TELDEC', 'f8'),
                                                     ('PMB', 'f8'),
                                                     ('FILTERNAME', 'a2')])
        fgcmExpInfo['VISIT'][:] = visitCat['visit']
        fgcmExpInfo['MJD'][:] = visitCat['mjd']
        fgcmExpInfo['EXPTIME'][:] = visitCat['exptime']
        fgcmExpInfo['DEEPFLAG'][:] = visitCat['deepflag']
        fgcmExpInfo['TELHA'][:] = visitCat['telha']
        fgcmExpInfo['TELRA'][:] = visitCat['telra']
        fgcmExpInfo['TELDEC'][:] = visitCat['teldec']
        fgcmExpInfo['PMB'][:] = visitCat['pmb']
        fgcmExpInfo['PSFSIGMA'][:] = visitCat['psfsigma']
        fgcmExpInfo['DELTA_APER'][:] = visitCat['deltaaper']
        fgcmExpInfo['SKYBACKGROUND'][:] = visitCat['skybackground']
        # Note that we have to go through asAstropy() to get a string
        #  array out of an afwTable
        fgcmExpInfo['FILTERNAME'][:] = visitCat.asAstropy()['filtername']

        return fgcmExpInfo

    def _loadCcdOffsets(self, butler):
        """
        """
        camera = butler.get('camera')

        # and we need to know the ccd offsets from the camera geometry
        ccdOffsets = np.zeros(len(camera), dtype=[('CCDNUM', 'i4'),
                                                  ('DELTA_RA', 'f8'),
                                                  ('DELTA_DEC', 'f8'),
                                                  ('RA_SIZE', 'f8'),
                                                  ('DEC_SIZE', 'f8'),
                                                  ('X_SIZE', 'i4'),
                                                  ('Y_SIZE', 'i4')])

        extent = afwGeom.Extent2D(self.config.pixelScale, self.config.pixelScale)

        for i, detector in enumerate(camera):
            # new version, using proper rotations
            #  but I worry this only works with HSC, as there's a unit inconsistency

            camPoint = detector.getCenter(afwCameraGeom.PIXELS)
            bbox = detector.getBBox()
            orient = detector.getOrientation()

            ccdOffsets['CCDNUM'][i] = detector.getId()

            xform = orient.makePixelFpTransform(extent)
            pointXform = xform.applyForward(camPoint)
            # this requires a pixelScale
            # NOTE that this now works properly with HSC, but I need to work on
            # generalizing this properly
            ccdOffsets['DELTA_RA'][i] = -pointXform.getY() * self.config.pixelScale / 3600.0
            ccdOffsets['DELTA_DEC'][i] = -pointXform.getX() * self.config.pixelScale / 3600.0

            # but this does not (for the delta)
            boxXform = xform.applyForward(afwGeom.Point2D(bbox.getMaxX(), bbox.getMaxY()))
            ccdOffsets['RA_SIZE'][i] = 2. * np.abs(boxXform.getY() -
                                                   pointXform.getY()) / 3600.0
            ccdOffsets['DEC_SIZE'][i] = 2. * np.abs(boxXform.getX() -
                                                    pointXform.getX()) / 3600.0

            ccdOffsets['X_SIZE'][i] = bbox.getMaxX()
            ccdOffsets['Y_SIZE'][i] = bbox.getMaxY()

        return ccdOffsets
