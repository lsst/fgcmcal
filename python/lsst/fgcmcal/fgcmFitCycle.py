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
"""Perform a single fit cycle of FGCM.

This task runs a single "fit cycle" of fgcm.  Prior to running this task
one must run both fgcmMakeLut (to construct the atmosphere and instrumental
look-up-table) and fgcmBuildStars (to extract visits and star observations
for the global fit).

The fgcmFitCycle is meant to be run multiple times, and is tracked by the
'cycleNumber'.  After each run of the fit cycle, diagnostic plots should
be inspected to set parameters for outlier rejection on the following
cycle.  Please see the fgcmcal Cookbook for details.
"""

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

__all__ = ['FgcmFitCycleConfig', 'FgcmFitCycleTask', 'FgcmFitCycleRunner']


class FgcmFitCycleConfig(pexConfig.Config):
    """Config for FgcmFitCycle"""

    bands = pexConfig.ListField(
        doc="Bands to run calibration (in wavelength order)",
        dtype=str,
        default=("NO_DATA",),
    )
    fitFlag = pexConfig.ListField(
        doc=("Flag for which bands are directly constrained in the FGCM fit. "
             "Bands set to 0 will have the atmosphere constrained from observations "
             "in other bands on the same night."),
        dtype=int,
        default=(0,),
    )
    requiredFlag = pexConfig.ListField(
        doc=("Flag for which bands are required for a star to be considered a calibration "
             "star in the FGCM fit.  Typically this should be the same as fitFlag."),
        dtype=int,
        default=(0,),
    )
    filterToBand = pexConfig.DictField(
        doc=("Dictionary to map filterName (e.g. physical filter) to band (e.g. abstract filter). "
             "With this mapping different filters (e.g. HSC r and r2) can be calibrated to the same "
             "'r' band."),
        keytype=str,
        itemtype=str,
        default={},
    )
    doReferenceCalibration = pexConfig.Field(
        doc="Use reference catalog as additional constraint on calibration",
        dtype=bool,
        default=False,
    )
    refStarSnMin = pexConfig.Field(
        doc="Reference star signal-to-noise minimum to use in calibration.  Set to <=0 for no cut.",
        dtype=float,
        default=20.0,
    )
    refStarOutlierNSig = pexConfig.Field(
        doc=("Number of sigma compared to average mag for reference star to be considered an outlier. "
             "Computed per-band, and if it is an outlier in any band it is rejected from fits."),
        dtype=float,
        default=4.0,
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
        doc=("Order of the 2D chebyshev polynomials for sub-ccd superstar fit. "
             "Global default is first-order polynomials, and should be overridden "
             "on a camera-by-camera basis depending on the ISR."),
        dtype=int,
        default=1,
    )
    superStarSubCcdTriangular = pexConfig.Field(
        doc=("Should the sub-ccd superstar chebyshev matrix be triangular to "
             "suppress high-order cross terms?"),
        dtype=bool,
        default=False,
    )
    superStarSigmaClip = pexConfig.Field(
        doc="Number of sigma to clip outliers when selecting for superstar flats",
        dtype=float,
        default=5.0,
    )
    ccdGraySubCcd = pexConfig.Field(
        doc="Compute CCD gray terms on sub-ccd scale",
        dtype=bool,
        default=False,
    )
    ccdGraySubCcdChebyshevOrder = pexConfig.Field(
        doc="Order of the 2D chebyshev polynomials for sub-ccd gray fit.",
        dtype=int,
        default=1,
    )
    ccdGraySubCcdTriangular = pexConfig.Field(
        doc=("Should the sub-ccd gray chebyshev matrix be triangular to "
             "suppress high-order cross terms?"),
        dtype=bool,
        default=True,
    )
    cycleNumber = pexConfig.Field(
        doc=("FGCM fit cycle number.  This is automatically incremented after each run "
             "and stage of outlier rejection.  See cookbook for details."),
        dtype=int,
        default=None,
    )
    isFinalCycle = pexConfig.Field(
        doc=("Is this the final cycle of the fitting?  Will automatically compute final "
             "selection of stars and photometric exposures, and will output zeropoints "
             "and standard stars for use in fgcmOutputProducts"),
        dtype=bool,
        default=False,
    )
    maxIterBeforeFinalCycle = pexConfig.Field(
        doc=("Maximum fit iterations, prior to final cycle.  The number of iterations "
             "will always be 0 in the final cycle for cleanup and final selection."),
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
    # TODO: When DM-16511 is done, it will be possible to get the
    # telescope latitude directly from the camera.
    latitude = pexConfig.Field(
        doc="Observatory latitude",
        dtype=float,
        default=None,
    )
    # TODO: DM-16490 will make this unneccessary
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
        doc=("Minimum number of good stars per CCD to be used in calibration fit. "
             "CCDs with fewer stars will have their calibration estimated from other "
             "CCDs in the same visit, with zeropoint error increased accordingly."),
        dtype=int,
        default=5,
    )
    minCcdPerExp = pexConfig.Field(
        doc=("Minimum number of good CCDs per exposure/visit to be used in calibration fit. "
             "Visits with fewer good CCDs will have CCD zeropoints estimated where possible."),
        dtype=int,
        default=5,
    )
    maxCcdGrayErr = pexConfig.Field(
        doc="Maximum error on CCD gray offset to be considered photometric",
        dtype=float,
        default=0.05,
    )
    minStarPerExp = pexConfig.Field(
        doc=("Minimum number of good stars per exposure/visit to be used in calibration fit. "
             "Visits with fewer good stars will have CCD zeropoints estimated where possible."),
        dtype=int,
        default=600,
    )
    minExpPerNight = pexConfig.Field(
        doc="Minimum number of good exposures/visits to consider a partly photometric night",
        dtype=int,
        default=10,
    )
    expGrayInitialCut = pexConfig.Field(
        doc=("Maximum exposure/visit gray value for initial selection of possible photometric "
             "observations."),
        dtype=float,
        default=-0.25,
    )
    expGrayPhotometricCut = pexConfig.ListField(
        doc=("Maximum (negative) exposure gray for a visit to be considered photometric. "
             "There will be one value per band."),
        dtype=float,
        default=(0.0,),
    )
    expGrayHighCut = pexConfig.ListField(
        doc=("Maximum (positive) exposure gray for a visit to be considered photometric. "
             "There will be one value per band."),
        dtype=float,
        default=(0.0,),
    )
    expGrayRecoverCut = pexConfig.Field(
        doc=("Maximum (negative) exposure gray to be able to recover bad ccds via interpolation. "
             "Visits with more gray extinction will only get CCD zeropoints if there are "
             "sufficient star observations (minStarPerCcd) on that CCD."),
        dtype=float,
        default=-1.0,
    )
    expVarGrayPhotometricCut = pexConfig.Field(
        doc="Maximum exposure variance to be considered possibly photometric",
        dtype=float,
        default=0.0005,
    )
    expGrayErrRecoverCut = pexConfig.Field(
        doc=("Maximum exposure gray error to be able to recover bad ccds via interpolation. "
             "Visits with more gray variance will only get CCD zeropoints if there are "
             "sufficient star observations (minStarPerCcd) on that CCD."),
        dtype=float,
        default=0.05,
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
    approxThroughput = pexConfig.ListField(
        doc="Approximate overall throughput at start of calibration observations",
        dtype=float,
        default=(1.0, ),
    )
    sigmaCalRange = pexConfig.ListField(
        doc="Allowed range for systematic error floor estimation",
        dtype=float,
        default=(0.001, 0.003),
    )
    sigmaCalFitPercentile = pexConfig.ListField(
        doc="Magnitude percentile range to fit systematic error floor",
        dtype=float,
        default=(0.05, 0.15),
    )
    sigmaCalPlotPercentile = pexConfig.ListField(
        doc="Magnitude percentile range to plot systematic error floor",
        dtype=float,
        default=(0.05, 0.95),
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
    useQuadraticPwv = pexConfig.Field(
        doc="Model PWV with a quadratic term for variation through the night?",
        dtype=bool,
        default=False,
    )
    outputStandardsBeforeFinalCycle = pexConfig.Field(
        doc="Output standard stars prior to final cycle?  Used in debugging.",
        dtype=bool,
        default=False,
    )
    outputZeropointsBeforeFinalCycle = pexConfig.Field(
        doc="Output standard stars prior to final cycle?  Used in debugging.",
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
    This Runner does not perform any dataRef parallelization, but the FGCM
    code called by the Task uses python multiprocessing (see the "ncores"
    config option).
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
        """

        task = self.TaskClass(config=self.config, log=self.log)

        exitStatus = 0
        if self.doRaise:
            task.runDataRef(butler)
        else:
            try:
                task.runDataRef(butler)
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
        butler : `lsst.daf.persistence.Butler`
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def runDataRef(self, butler):
        """
        Run a single fit cycle for FGCM

        Parameters
        ----------
        butler:  `lsst.daf.persistence.Butler`
        """

        self._fgcmFitCycle(butler)

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
        doBackup : `bool`, optional
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
        butler: `lsst.daf.persistence.Butler`
        """

        self._checkDatasetsExist(butler)

        # Set defaults on whether to output standards and zeropoints
        self.maxIter = self.config.maxIterBeforeFinalCycle
        self.outputStandards = self.config.outputStandardsBeforeFinalCycle
        self.outputZeropoints = self.config.outputZeropointsBeforeFinalCycle
        self.resetFitParameters = True

        if self.config.isFinalCycle:
            # This is the final fit cycle, so we do not want to reset fit
            # parameters, we want to run a final "clean-up" with 0 fit iterations,
            # and we always want to output standards and zeropoints
            self.maxIter = 0
            self.outputStandards = True
            self.outputZeropoints = True
            self.resetFitParameters = False

        camera = butler.get('camera')
        configDict = self._makeConfigDict(camera)

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
            inParInfo, inParams, inSuperStar = self._loadParameters(butler)
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
            flagId = flaggedStars['objId'][:]
            flagFlag = flaggedStars['objFlag'][:]
        else:
            flagId = None
            flagFlag = None

        if self.config.doReferenceCalibration:
            refStars = butler.get('fgcmReferenceStars')
            refId = refStars['fgcm_id'][:]
            refMag = refStars['refMag'][:, :]
            refMagErr = refStars['refMagErr'][:, :]
        else:
            refId = None
            refMag = None
            refMagErr = None

        # match star observations to visits
        # Only those star observations that match visits from fgcmExpInfo['VISIT'] will
        # actually be transferred into fgcm using the indexing below.
        visitIndex = np.searchsorted(fgcmExpInfo['VISIT'], starObs['visit'][starIndices['obsIndex']])

        # The fgcmStars.loadStars method will copy all the star information into
        # special shared memory objects that will not blow up the memory usage when
        # used with python multiprocessing.  Once all the numbers are copied,
        # it is necessary to release all references to the objects that previously
        # stored the data to ensure that the garbage collector can clear the memory,
        # and ensure that this memory is not copied when multiprocessing kicks in.

        # We determine the conversion from the native units (typically radians) to
        # degrees for the first star.  This allows us to treat coord_ra/coord_dec as
        # numpy arrays rather than Angles, which would we approximately 600x slower.
        conv = starObs[0]['ra'].asDegrees() / float(starObs[0]['ra'])

        fgcmStars.loadStars(fgcmPars,
                            starObs['visit'][starIndices['obsIndex']],
                            starObs['ccd'][starIndices['obsIndex']],
                            starObs['ra'][starIndices['obsIndex']] * conv,
                            starObs['dec'][starIndices['obsIndex']] * conv,
                            starObs['instMag'][starIndices['obsIndex']],
                            starObs['instMagErr'][starIndices['obsIndex']],
                            fgcmExpInfo['FILTERNAME'][visitIndex],
                            starIds['fgcm_id'][:],
                            starIds['ra'][:],
                            starIds['dec'][:],
                            starIds['obsArrIndex'][:],
                            starIds['nObs'][:],
                            obsX=starObs['x'][starIndices['obsIndex']],
                            obsY=starObs['y'][starIndices['obsIndex']],
                            refID=refId,
                            refMag=refMag,
                            refMagErr=refMagErr,
                            flagID=flagId,
                            flagFlag=flagFlag,
                            computeNobs=True)

        # Release all references to temporary objects holding star data (see above)
        starObs = None
        starIds = None
        starIndices = None
        flagId = None
        flagFlag = None
        flaggedStars = None
        refStars = None

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

        self._persistFgcmDatasets(butler, fgcmFitCycle)

        # Output the config for the next cycle
        # We need to make a copy since the input one has been frozen

        outConfig = FgcmFitCycleConfig()
        outConfig.update(**self.config.toDict())

        outConfig.cycleNumber += 1
        outConfig.precomputeSuperStarInitialCycle = False
        outConfig.freezeStdAtmosphere = False
        outConfig.expGrayPhotometricCut[:] = fgcmFitCycle.updatedPhotometricCut
        outConfig.expGrayHighCut[:] = fgcmFitCycle.updatedHighCut
        configFileName = '%s_cycle%02d_config.py' % (outConfig.outfileBase,
                                                     outConfig.cycleNumber)
        outConfig.save(configFileName)

        if self.config.isFinalCycle == 1:
            # We are done, ready to output products
            self.log.info("Everything is in place to run fgcmOutputProducts.py")
        else:
            self.log.info("Saved config for next cycle to %s" % (configFileName))
            self.log.info("Be sure to look at:")
            self.log.info("   config.expGrayPhotometricCut")
            self.log.info("   config.expGrayHighCut")
            self.log.info("If you are satisfied with the fit, please set:")
            self.log.info("   config.isFinalCycle = True")

    def _checkDatasetsExist(self, butler):
        """
        Check if necessary datasets exist to run fgcmFitCycle

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`

        Raises
        ------
        RuntimeError
           If any of fgcmVisitCatalog, fgcmStarObservations, fgcmStarIds,
           fgcmStarIndices, fgcmLookUpTable datasets do not exist.
           If cycleNumber > 0, then also checks for fgcmFitParameters,
           fgcmFlaggedStars.
        """

        if not butler.datasetExists('fgcmVisitCatalog'):
            raise RuntimeError("Could not find fgcmVisitCatalog in repo!")
        if not butler.datasetExists('fgcmStarObservations'):
            raise RuntimeError("Could not find fgcmStarObservations in repo!")
        if not butler.datasetExists('fgcmStarIds'):
            raise RuntimeError("Could not find fgcmStarIds in repo!")
        if not butler.datasetExists('fgcmStarIndices'):
            raise RuntimeError("Could not find fgcmStarIndices in repo!")
        if not butler.datasetExists('fgcmLookUpTable'):
            raise RuntimeError("Could not find fgcmLookUpTable in repo!")

        # Need additional datasets if we are not the initial cycle
        if (self.config.cycleNumber > 0):
            if not butler.datasetExists('fgcmFitParameters',
                                        fgcmcycle=self.config.cycleNumber-1):
                raise RuntimeError("Could not find fgcmFitParameters for previous cycle (%d) in repo!" %
                                   (self.config.cycleNumber-1))
            if not butler.datasetExists('fgcmFlaggedStars',
                                        fgcmcycle=self.config.cycleNumber-1):
                raise RuntimeError("Could not find fgcmFlaggedStars for previous cycle (%d) in repo!" %
                                   (self.config.cycleNumber-1))

        # And additional dataset if we want reference calibration
        if self.config.doReferenceCalibration:
            if not butler.datasetExists('fgcmReferenceStars'):
                raise RuntimeError("Could not find fgcmReferenceStars in repo, and "
                                   "doReferenceCalibration is True.")

    def _loadFgcmLut(self, butler, filterToBand=None):
        """
        Load the FGCM look-up-table

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        filterToBand: `dict`
           Dictionary mapping filters to bands (see self.config.filterToBand)

        Returns
        -------
        fgcmLut: `lsst.fgcm.FgcmLut`
           Lookup table for FGCM
        lutIndexVals: `np.ndarray`
           Numpy array with LUT index information for FGCM
        lutStd: `np.ndarray`
           Numpy array with LUT standard throughput values for FGCM
        """

        # set up the look-up-table
        lutCat = butler.get('fgcmLookUpTable')

        # first we need the lutIndexVals
        # dtype is set for py2/py3/fits/fgcm compatibility
        lutFilterNames = np.array(lutCat[0]['filterNames'].split(','), dtype='a')
        lutStdFilterNames = np.array(lutCat[0]['stdFilterNames'].split(','), dtype='a')

        # Note that any discrepancies between config values will raise relevant
        # exceptions in the FGCM code.

        lutIndexVals = np.zeros(1, dtype=[('FILTERNAMES', lutFilterNames.dtype.str,
                                           lutFilterNames.size),
                                          ('STDFILTERNAMES', lutStdFilterNames.dtype.str,
                                           lutStdFilterNames.size),
                                          ('PMB', 'f8', lutCat[0]['pmb'].size),
                                          ('PMBFACTOR', 'f8', lutCat[0]['pmbFactor'].size),
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
        lutIndexVals['PMBFACTOR'][:] = lutCat[0]['pmbFactor']
        lutIndexVals['PMBELEVATION'] = lutCat[0]['pmbElevation']
        lutIndexVals['LAMBDANORM'] = lutCat[0]['lambdaNorm']
        lutIndexVals['PWV'][:] = lutCat[0]['pwv']
        lutIndexVals['O3'][:] = lutCat[0]['o3']
        lutIndexVals['TAU'][:] = lutCat[0]['tau']
        lutIndexVals['ALPHA'][:] = lutCat[0]['alpha']
        lutIndexVals['ZENITH'][:] = lutCat[0]['zenith']
        lutIndexVals['NCCD'] = lutCat[0]['nCcd']

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
                                    ('ATMLAMBDA', 'f8', lutCat[0]['atmLambda'].size),
                                    ('ATMSTDTRANS', 'f8', lutCat[0]['atmStdTrans'].size)])
        lutStd['PMBSTD'] = lutCat[0]['pmbStd']
        lutStd['PWVSTD'] = lutCat[0]['pwvStd']
        lutStd['O3STD'] = lutCat[0]['o3Std']
        lutStd['TAUSTD'] = lutCat[0]['tauStd']
        lutStd['ALPHASTD'] = lutCat[0]['alphaStd']
        lutStd['ZENITHSTD'] = lutCat[0]['zenithStd']
        lutStd['LAMBDARANGE'][:] = lutCat[0]['lambdaRange'][:]
        lutStd['LAMBDASTEP'] = lutCat[0]['lambdaStep']
        lutStd['LAMBDASTD'][:] = lutCat[0]['lambdaStd']
        lutStd['LAMBDASTDFILTER'][:] = lutCat[0]['lambdaStdFilter']
        lutStd['I0STD'][:] = lutCat[0]['i0Std']
        lutStd['I1STD'][:] = lutCat[0]['i1Std']
        lutStd['I10STD'][:] = lutCat[0]['i10Std']
        lutStd['LAMBDAB'][:] = lutCat[0]['lambdaB']
        lutStd['ATMLAMBDA'][:] = lutCat[0]['atmLambda'][:]
        lutStd['ATMSTDTRANS'][:] = lutCat[0]['atmStdTrans'][:]

        lutTypes = [row['luttype'] for row in lutCat]

        # And the flattened look-up-table
        lutFlat = np.zeros(lutCat[0]['lut'].size, dtype=[('I0', 'f4'),
                                                         ('I1', 'f4')])

        lutFlat['I0'][:] = lutCat[lutTypes.index('I0')]['lut'][:]
        lutFlat['I1'][:] = lutCat[lutTypes.index('I1')]['lut'][:]

        lutDerivFlat = np.zeros(lutCat[0]['lut'].size, dtype=[('D_LNPWV', 'f4'),
                                                              ('D_O3', 'f4'),
                                                              ('D_LNTAU', 'f4'),
                                                              ('D_ALPHA', 'f4'),
                                                              ('D_SECZENITH', 'f4'),
                                                              ('D_LNPWV_I1', 'f4'),
                                                              ('D_O3_I1', 'f4'),
                                                              ('D_LNTAU_I1', 'f4'),
                                                              ('D_ALPHA_I1', 'f4'),
                                                              ('D_SECZENITH_I1', 'f4')])

        for name in lutDerivFlat.dtype.names:
            lutDerivFlat[name][:] = lutCat[lutTypes.index(name)]['lut'][:]

        # The fgcm.FgcmLUT() class copies all the LUT information into special
        # shared memory objects that will not blow up the memory usage when used
        # with python multiprocessing.  Once all the numbers are copied, the
        # references to the temporary objects (lutCat, lutFlat, lutDerivFlat)
        # will fall out of scope and can be cleaned up by the garbage collector.
        fgcmLut = fgcm.FgcmLUT(lutIndexVals, lutFlat, lutDerivFlat, lutStd,
                               filterToBand=self.config.filterToBand)

        return fgcmLut, lutIndexVals, lutStd

    def _loadVisitCatalog(self, butler):
        """
        Load the FGCM visit catalog

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`

        Returns
        -------
        fgcmExpInfo: `np.ndarray`
           Numpy array for visit information for FGCM
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
        fgcmExpInfo['DEEPFLAG'][:] = visitCat['deepFlag']
        fgcmExpInfo['TELHA'][:] = visitCat['telha']
        fgcmExpInfo['TELRA'][:] = visitCat['telra']
        fgcmExpInfo['TELDEC'][:] = visitCat['teldec']
        fgcmExpInfo['PMB'][:] = visitCat['pmb']
        fgcmExpInfo['PSFSIGMA'][:] = visitCat['psfSigma']
        fgcmExpInfo['DELTA_APER'][:] = visitCat['deltaAper']
        fgcmExpInfo['SKYBACKGROUND'][:] = visitCat['skyBackground']
        # Note that we have to go through asAstropy() to get a string
        #  array out of an afwTable.  This is faster than a row-by-row loop.
        fgcmExpInfo['FILTERNAME'][:] = visitCat.asAstropy()['filtername']

        return fgcmExpInfo

    def _loadCcdOffsets(self, butler):
        """
        Load the CCD offsets in ra/dec and x/y space

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`

        Returns
        -------
        ccdOffsets: `np.ndarray`
           Numpy array with ccd offset information for input to FGCM
        """
        # TODO: DM-16490 will simplify and generalize the math.
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
            # Note that this now works properly with HSC, but I need to work on
            # generalizing this properly.  I expect the updates in DM-16490 will
            # generalize these computations, and appropriate tests can be added
            # on that ticket.
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

    def _makeConfigDict(self, camera):
        """
        Make the FGCM configuration dict

        Parameters
        ----------
        camera: 'cameraGeom.Camera`
            Camera from the butler

        Returns
        -------
        configDict: `dict`
        """

        fitFlag = np.array(self.config.fitFlag, dtype=np.bool)
        requiredFlag = np.array(self.config.requiredFlag, dtype=np.bool)

        fitBands = [b for i, b in enumerate(self.config.bands) if fitFlag[i]]
        notFitBands = [b for i, b in enumerate(self.config.bands) if not fitFlag[i]]
        requiredBands = [b for i, b in enumerate(self.config.bands) if requiredFlag[i]]

        # process the starColorCuts
        starColorCutList = []
        for ccut in self.config.starColorCuts:
            parts = ccut.split(',')
            starColorCutList.append([parts[0], parts[1], float(parts[2]), float(parts[3])])

        # TODO: Having direct access to the mirror area from the camera would be
        #  useful.  See DM-16489.
        # Mirror area in cm**2
        mirrorArea = np.pi*(camera.telescopeDiameter*100./2.)**2.

        # Get approximate average camera gain:
        gains = [amp.getGain() for detector in camera for amp in detector.getAmpInfoCatalog()]
        cameraGain = float(np.median(gains))

        # create a configuration dictionary for fgcmFitCycle
        configDict = {'outfileBase': self.config.outfileBase,
                      'logger': self.log,
                      'exposureFile': None,
                      'obsFile': None,
                      'indexFile': None,
                      'lutFile': None,
                      'mirrorArea': mirrorArea,
                      'cameraGain': cameraGain,
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
                      'superStarSubCCDTriangular': self.config.superStarSubCcdTriangular,
                      'superStarSigmaClip': self.config.superStarSigmaClip,
                      'ccdGraySubCCD': self.config.ccdGraySubCcd,
                      'ccdGraySubCCDChebyshevOrder': self.config.ccdGraySubCcdChebyshevOrder,
                      'ccdGraySubCCDTriangular': self.config.ccdGraySubCcdTriangular,
                      'cycleNumber': self.config.cycleNumber,
                      'maxIter': self.maxIter,
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
                      'refStarSnMin': self.config.refStarSnMin,
                      'refStarOutlierNSig': self.config.refStarOutlierNSig,
                      'illegalValue': -9999.0,  # internally used by fgcm.
                      'starColorCuts': starColorCutList,
                      'aperCorrFitNBins': self.config.aperCorrFitNBins,
                      'sedFudgeFactors': np.array(self.config.sedFudgeFactors),
                      'colorSplitIndices': np.array(self.config.colorSplitIndices),
                      'sigFgcmMaxErr': self.config.sigFgcmMaxErr,
                      'sigFgcmMaxEGray': self.config.sigFgcmMaxEGray,
                      'ccdGrayMaxStarErr': self.config.ccdGrayMaxStarErr,
                      'approxThroughput': list(self.config.approxThroughput),
                      'sigmaCalRange': list(self.config.sigmaCalRange),
                      'sigmaCalFitPercentile': list(self.config.sigmaCalFitPercentile),
                      'sigmaCalPlotPercentile': list(self.config.sigmaCalPlotPercentile),
                      'sigma0Phot': self.config.sigma0Phot,
                      'mapLongitudeRef': self.config.mapLongitudeRef,
                      'mapNSide': self.config.mapNSide,
                      'varNSig': 100.0,  # Turn off 'variable star selection' which doesn't work yet
                      'varMinBand': 2,
                      'useRetrievedPwv': False,
                      'useNightlyRetrievedPwv': False,
                      'pwvRetrievalSmoothBlock': 25,
                      'useQuadraticPwv': self.config.useQuadraticPwv,
                      'useRetrievedTauInit': False,
                      'tauRetrievalMinCCDPerNight': 500,
                      'modelMagErrors': self.config.modelMagErrors,
                      'printOnly': False,
                      'outputStars': False,
                      'clobber': True,
                      'useSedLUT': False,
                      'resetParameters': self.resetFitParameters,
                      'outputZeropoints': self.outputZeropoints}

        return configDict

    def _loadParameters(self, butler):
        """
        Load FGCM parameters from a previous fit cycle

        Parameters
        ----------
        butler:  `lsst.daf.persistence.Butler`

        Returns
        -------
        inParInfo: `np.ndarray`
           Numpy array parameter information formatted for input to fgcm
        inParameters: `np.ndarray`
           Numpy array parameter values formatted for input to fgcm
        inSuperStar: `np.array`
           Superstar flat formatted for input to fgcm
        """

        # note that we already checked that this is available
        parCat = butler.get('fgcmFitParameters', fgcmcycle=self.config.cycleNumber-1)

        parLutFilterNames = np.array(parCat[0]['lutFilterNames'].split(','))
        parFitBands = np.array(parCat[0]['fitBands'].split(','))
        parNotFitBands = np.array(parCat[0]['notFitBands'].split(','))

        inParInfo = np.zeros(1, dtype=[('NCCD', 'i4'),
                                       ('LUTFILTERNAMES', parLutFilterNames.dtype.str,
                                        parLutFilterNames.size),
                                       ('FITBANDS', parFitBands.dtype.str, parFitBands.size),
                                       ('NOTFITBANDS', parNotFitBands.dtype.str, parNotFitBands.size),
                                       ('LNTAUUNIT', 'f8'),
                                       ('LNTAUSLOPEUNIT', 'f8'),
                                       ('ALPHAUNIT', 'f8'),
                                       ('LNPWVUNIT', 'f8'),
                                       ('LNPWVSLOPEUNIT', 'f8'),
                                       ('LNPWVQUADRATICUNIT', 'f8'),
                                       ('LNPWVGLOBALUNIT', 'f8'),
                                       ('O3UNIT', 'f8'),
                                       ('QESYSUNIT', 'f8'),
                                       ('QESYSSLOPEUNIT', 'f8'),
                                       ('FILTEROFFSETUNIT', 'f8'),
                                       ('HASEXTERNALPWV', 'i2'),
                                       ('HASEXTERNALTAU', 'i2')])
        inParInfo['NCCD'] = parCat['nCcd']
        inParInfo['LUTFILTERNAMES'][:] = parLutFilterNames
        inParInfo['FITBANDS'][:] = parFitBands
        inParInfo['NOTFITBANDS'][:] = parNotFitBands
        inParInfo['LNTAUUNIT'] = parCat['lnTauUnit']
        inParInfo['LNTAUSLOPEUNIT'] = parCat['lnTauSlopeUnit']
        inParInfo['ALPHAUNIT'] = parCat['alphaUnit']
        inParInfo['LNPWVUNIT'] = parCat['lnPwvUnit']
        inParInfo['LNPWVSLOPEUNIT'] = parCat['lnPwvSlopeUnit']
        inParInfo['LNPWVQUADRATICUNIT'] = parCat['lnPwvQuadraticUnit']
        inParInfo['LNPWVGLOBALUNIT'] = parCat['lnPwvGlobalUnit']
        inParInfo['O3UNIT'] = parCat['o3Unit']
        inParInfo['QESYSUNIT'] = parCat['qeSysUnit']
        inParInfo['QESYSSLOPEUNIT'] = parCat['qeSysSlopeUnit']
        inParInfo['FILTEROFFSETUNIT'] = parCat['filterOffsetUnit']
        inParInfo['HASEXTERNALPWV'] = parCat['hasExternalPwv']
        inParInfo['HASEXTERNALTAU'] = parCat['hasExternalTau']

        inParams = np.zeros(1, dtype=[('PARALPHA', 'f8', parCat['parAlpha'].size),
                                      ('PARO3', 'f8', parCat['parO3'].size),
                                      ('PARLNTAUINTERCEPT', 'f8',
                                       parCat['parLnTauIntercept'].size),
                                      ('PARLNTAUSLOPE', 'f8',
                                       parCat['parLnTauSlope'].size),
                                      ('PARLNPWVINTERCEPT', 'f8',
                                       parCat['parLnPwvIntercept'].size),
                                      ('PARLNPWVSLOPE', 'f8',
                                       parCat['parLnPwvSlope'].size),
                                      ('PARLNPWVQUADRATIC', 'f8',
                                       parCat['parLnPwvQuadratic'].size),
                                      ('PARQESYSINTERCEPT', 'f8',
                                       parCat['parQeSysIntercept'].size),
                                      ('PARQESYSSLOPE', 'f8',
                                       parCat['parQeSysSlope'].size),
                                      ('PARFILTEROFFSET', 'f8',
                                       parCat['parFilterOffset'].size),
                                      ('PARFILTEROFFSETFITFLAG', 'i2',
                                       parCat['parFilterOffsetFitFlag'].size),
                                      ('PARRETRIEVEDLNPWVSCALE', 'f8'),
                                      ('PARRETRIEVEDLNPWVOFFSET', 'f8'),
                                      ('PARRETRIEVEDLNPWVNIGHTLYOFFSET', 'f8',
                                       parCat['parRetrievedLnPwvNightlyOffset'].size),
                                      ('COMPABSTHROUGHPUT', 'f8',
                                       parCat['compAbsThroughput'].size),
                                      ('COMPREFOFFSET', 'f8',
                                       parCat['compRefOffset'].size),
                                      ('COMPREFSIGMA', 'f8',
                                       parCat['compRefSigma'].size),
                                      ('COMPAPERCORRPIVOT', 'f8',
                                       parCat['compAperCorrPivot'].size),
                                      ('COMPAPERCORRSLOPE', 'f8',
                                       parCat['compAperCorrSlope'].size),
                                      ('COMPAPERCORRSLOPEERR', 'f8',
                                       parCat['compAperCorrSlopeErr'].size),
                                      ('COMPAPERCORRRANGE', 'f8',
                                       parCat['compAperCorrRange'].size),
                                      ('COMPMODELERREXPTIMEPIVOT', 'f8',
                                       parCat['compModelErrExptimePivot'].size),
                                      ('COMPMODELERRFWHMPIVOT', 'f8',
                                       parCat['compModelErrFwhmPivot'].size),
                                      ('COMPMODELERRSKYPIVOT', 'f8',
                                       parCat['compModelErrSkyPivot'].size),
                                      ('COMPMODELERRPARS', 'f8',
                                       parCat['compModelErrPars'].size),
                                      ('COMPEXPGRAY', 'f8',
                                       parCat['compExpGray'].size),
                                      ('COMPVARGRAY', 'f8',
                                       parCat['compVarGray'].size),
                                      ('COMPNGOODSTARPEREXP', 'i4',
                                       parCat['compNGoodStarPerExp'].size),
                                      ('COMPSIGFGCM', 'f8',
                                       parCat['compSigFgcm'].size),
                                      ('COMPSIGMACAL', 'f8',
                                       parCat['compSigmaCal'].size),
                                      ('COMPRETRIEVEDLNPWV', 'f8',
                                       parCat['compRetrievedLnPwv'].size),
                                      ('COMPRETRIEVEDLNPWVRAW', 'f8',
                                       parCat['compRetrievedLnPwvRaw'].size),
                                      ('COMPRETRIEVEDLNPWVFLAG', 'i2',
                                       parCat['compRetrievedLnPwvFlag'].size),
                                      ('COMPRETRIEVEDTAUNIGHT', 'f8',
                                       parCat['compRetrievedTauNight'].size)])

        inParams['PARALPHA'][:] = parCat['parAlpha'][0, :]
        inParams['PARO3'][:] = parCat['parO3'][0, :]
        inParams['PARLNTAUINTERCEPT'][:] = parCat['parLnTauIntercept'][0, :]
        inParams['PARLNTAUSLOPE'][:] = parCat['parLnTauSlope'][0, :]
        inParams['PARLNPWVINTERCEPT'][:] = parCat['parLnPwvIntercept'][0, :]
        inParams['PARLNPWVSLOPE'][:] = parCat['parLnPwvSlope'][0, :]
        inParams['PARLNPWVQUADRATIC'][:] = parCat['parLnPwvQuadratic'][0, :]
        inParams['PARQESYSINTERCEPT'][:] = parCat['parQeSysIntercept'][0, :]
        inParams['PARQESYSSLOPE'][:] = parCat['parQeSysSlope'][0, :]
        inParams['PARFILTEROFFSET'][:] = parCat['parFilterOffset'][0, :]
        inParams['PARFILTEROFFSETFITFLAG'][:] = parCat['parFilterOffsetFitFlag'][0, :]
        inParams['PARRETRIEVEDLNPWVSCALE'] = parCat['parRetrievedLnPwvScale']
        inParams['PARRETRIEVEDLNPWVOFFSET'] = parCat['parRetrievedLnPwvOffset']
        inParams['PARRETRIEVEDLNPWVNIGHTLYOFFSET'][:] = parCat['parRetrievedLnPwvNightlyOffset'][0, :]
        inParams['COMPABSTHROUGHPUT'][:] = parCat['compAbsThroughput'][0, :]
        inParams['COMPREFOFFSET'][:] = parCat['compRefOffset'][0, :]
        inParams['COMPREFSIGMA'][:] = parCat['compRefSigma'][0, :]
        inParams['COMPAPERCORRPIVOT'][:] = parCat['compAperCorrPivot'][0, :]
        inParams['COMPAPERCORRSLOPE'][:] = parCat['compAperCorrSlope'][0, :]
        inParams['COMPAPERCORRSLOPEERR'][:] = parCat['compAperCorrSlopeErr'][0, :]
        inParams['COMPAPERCORRRANGE'][:] = parCat['compAperCorrRange'][0, :]
        inParams['COMPMODELERREXPTIMEPIVOT'][:] = parCat['compModelErrExptimePivot'][0, :]
        inParams['COMPMODELERRFWHMPIVOT'][:] = parCat['compModelErrFwhmPivot'][0, :]
        inParams['COMPMODELERRSKYPIVOT'][:] = parCat['compModelErrSkyPivot'][0, :]
        inParams['COMPMODELERRPARS'][:] = parCat['compModelErrPars'][0, :]
        inParams['COMPEXPGRAY'][:] = parCat['compExpGray'][0, :]
        inParams['COMPVARGRAY'][:] = parCat['compVarGray'][0, :]
        inParams['COMPNGOODSTARPEREXP'][:] = parCat['compNGoodStarPerExp'][0, :]
        inParams['COMPSIGFGCM'][:] = parCat['compSigFgcm'][0, :]
        inParams['COMPSIGMACAL'][:] = parCat['compSigmaCal'][0, :]
        inParams['COMPRETRIEVEDLNPWV'][:] = parCat['compRetrievedLnPwv'][0, :]
        inParams['COMPRETRIEVEDLNPWVRAW'][:] = parCat['compRetrievedLnPwvRaw'][0, :]
        inParams['COMPRETRIEVEDLNPWVFLAG'][:] = parCat['compRetrievedLnPwvFlag'][0, :]
        inParams['COMPRETRIEVEDTAUNIGHT'][:] = parCat['compRetrievedTauNight'][0, :]

        inSuperStar = np.zeros(parCat['superstarSize'][0, :], dtype='f8')
        inSuperStar[:, :, :, :] = parCat['superstar'][0, :].reshape(inSuperStar.shape)

        return (inParInfo, inParams, inSuperStar)

    def _persistFgcmDatasets(self, butler, fgcmFitCycle):
        """
        Persist FGCM datasets through the butler.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        fgcmFitCycle: `lsst.fgcm.FgcmFitCycle`
           Fgcm Fit cycle object
        """

        # Save the parameters
        parInfo, pars = fgcmFitCycle.fgcmPars.parsToArrays()

        parSchema = afwTable.Schema()

        comma = ','
        lutFilterNameString = comma.join([n.decode('utf-8')
                                          for n in parInfo['LUTFILTERNAMES'][0]])
        fitBandString = comma.join([n.decode('utf-8')
                                    for n in parInfo['FITBANDS'][0]])
        notFitBandString = comma.join([n.decode('utf-8')
                                       for n in parInfo['NOTFITBANDS'][0]])

        parSchema = self._makeParSchema(parInfo, pars, fgcmFitCycle.fgcmPars.parSuperStarFlat,
                                        lutFilterNameString, fitBandString, notFitBandString)
        parCat = self._makeParCatalog(parSchema, parInfo, pars,
                                      fgcmFitCycle.fgcmPars.parSuperStarFlat,
                                      lutFilterNameString, fitBandString, notFitBandString)

        butler.put(parCat, 'fgcmFitParameters', fgcmcycle=self.config.cycleNumber)

        # Save the indices of the flagged stars
        # (stars that have been (a) reserved from the fit for testing and
        # (b) bad stars that have failed quality checks.)
        flagStarSchema = self._makeFlagStarSchema()
        flagStarStruct = fgcmFitCycle.fgcmStars.getFlagStarIndices()
        flagStarCat = self._makeFlagStarCat(flagStarSchema, flagStarStruct)

        butler.put(flagStarCat, 'fgcmFlaggedStars', fgcmcycle=self.config.cycleNumber)

        # Save the zeropoint information and atmospheres only if desired
        if self.outputZeropoints:
            if self.config.superStarSubCcd or self.config.ccdGraySubCcd:
                chebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB'].shape[1]
            else:
                chebSize = 0
            zptSchema = self._makeZptSchema(chebSize)
            zptCat = self._makeZptCat(zptSchema, fgcmFitCycle.fgcmZpts.zpStruct)

            butler.put(zptCat, 'fgcmZeropoints', fgcmcycle=self.config.cycleNumber)

            # Save atmosphere values
            # These are generated by the same code that generates zeropoints
            atmSchema = self._makeAtmSchema()
            atmCat = self._makeAtmCat(atmSchema, fgcmFitCycle.fgcmZpts.atmStruct)

            butler.put(atmCat, 'fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        # Save the standard stars (if configured)
        if self.outputStandards:
            stdSchema = self._makeStdSchema()
            stdStruct = fgcmFitCycle.fgcmStars.retrieveStdStarCatalog(fgcmFitCycle.fgcmPars)
            stdCat = self._makeStdCat(stdSchema, stdStruct)

            butler.put(stdCat, 'fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

    def _makeParSchema(self, parInfo, pars, parSuperStarFlat,
                       lutFilterNameString, fitBandString, notFitBandString):
        """
        Make the parameter persistence schema

        Parameters
        ----------
        parInfo: `np.ndarray`
           Parameter information returned by fgcm
        pars: `np.ndarray`
           Parameter values returned by fgcm
        parSuperStarFlat: `np.array`
           Superstar flat values returned by fgcm
        lutFilterNameString: `str`
           Combined string of all the lutFilterNames
        fitBandString: `str`
           Combined string of all the fitBands
        notFitBandString: `str`
           Combined string of all the bands not used in the fit

        Returns
        -------
        parSchema: `afwTable.schema`
        """

        parSchema = afwTable.Schema()

        # parameter info section
        parSchema.addField('nCcd', type=np.int32, doc='Number of CCDs')
        parSchema.addField('lutFilterNames', type=str, doc='LUT Filter names in parameter file',
                           size=len(lutFilterNameString))
        parSchema.addField('fitBands', type=str, doc='Bands that were fit',
                           size=len(fitBandString))
        parSchema.addField('notFitBands', type=str, doc='Bands that were not fit',
                           size=len(notFitBandString))
        parSchema.addField('lnTauUnit', type=np.float64, doc='Step units for ln(AOD)')
        parSchema.addField('lnTauSlopeUnit', type=np.float64,
                           doc='Step units for ln(AOD) slope')
        parSchema.addField('alphaUnit', type=np.float64, doc='Step units for alpha')
        parSchema.addField('lnPwvUnit', type=np.float64, doc='Step units for ln(pwv)')
        parSchema.addField('lnPwvSlopeUnit', type=np.float64,
                           doc='Step units for ln(pwv) slope')
        parSchema.addField('lnPwvQuadraticUnit', type=np.float64,
                           doc='Step units for ln(pwv) quadratic term')
        parSchema.addField('lnPwvGlobalUnit', type=np.float64,
                           doc='Step units for global ln(pwv) parameters')
        parSchema.addField('o3Unit', type=np.float64, doc='Step units for O3')
        parSchema.addField('qeSysUnit', type=np.float64, doc='Step units for mirror gray')
        parSchema.addField('qeSysSlopeUnit', type=np.float64, doc='Step units for mirror gray slope')
        parSchema.addField('filterOffsetUnit', type=np.float64, doc='Step units for filter offset')
        parSchema.addField('hasExternalPwv', type=np.int32, doc='Parameters fit using external pwv')
        parSchema.addField('hasExternalTau', type=np.int32, doc='Parameters fit using external tau')

        # parameter section
        parSchema.addField('parAlpha', type='ArrayD', doc='Alpha parameter vector',
                           size=pars['PARALPHA'].size)
        parSchema.addField('parO3', type='ArrayD', doc='O3 parameter vector',
                           size=pars['PARO3'].size)
        parSchema.addField('parLnTauIntercept', type='ArrayD',
                           doc='ln(Tau) intercept parameter vector',
                           size=pars['PARLNTAUINTERCEPT'].size)
        parSchema.addField('parLnTauSlope', type='ArrayD',
                           doc='ln(Tau) slope parameter vector',
                           size=pars['PARLNTAUSLOPE'].size)
        parSchema.addField('parLnPwvIntercept', type='ArrayD', doc='ln(pwv) intercept parameter vector',
                           size=pars['PARLNPWVINTERCEPT'].size)
        parSchema.addField('parLnPwvSlope', type='ArrayD', doc='ln(pwv) slope parameter vector',
                           size=pars['PARLNPWVSLOPE'].size)
        parSchema.addField('parLnPwvQuadratic', type='ArrayD', doc='ln(pwv) quadratic parameter vector',
                           size=pars['PARLNPWVQUADRATIC'].size)
        parSchema.addField('parQeSysIntercept', type='ArrayD', doc='Mirror gray intercept parameter vector',
                           size=pars['PARQESYSINTERCEPT'].size)
        parSchema.addField('parQeSysSlope', type='ArrayD', doc='Mirror gray slope parameter vector',
                           size=pars['PARQESYSSLOPE'].size)
        parSchema.addField('parFilterOffset', type='ArrayD', doc='Filter offset parameter vector',
                           size=pars['PARFILTEROFFSET'].size)
        parSchema.addField('parFilterOffsetFitFlag', type='ArrayI', doc='Filter offset parameter fit flag',
                           size=pars['PARFILTEROFFSETFITFLAG'].size)
        parSchema.addField('parRetrievedLnPwvScale', type=np.float64,
                           doc='Global scale for retrieved ln(pwv)')
        parSchema.addField('parRetrievedLnPwvOffset', type=np.float64,
                           doc='Global offset for retrieved ln(pwv)')
        parSchema.addField('parRetrievedLnPwvNightlyOffset', type='ArrayD',
                           doc='Nightly offset for retrieved ln(pwv)',
                           size=pars['PARRETRIEVEDLNPWVNIGHTLYOFFSET'].size)
        parSchema.addField('compAbsThroughput', type='ArrayD',
                           doc='Absolute throughput (relative to transmission curves)',
                           size=pars['COMPABSTHROUGHPUT'].size)
        parSchema.addField('compRefOffset', type='ArrayD',
                           doc='Offset between reference stars and calibrated stars',
                           size=pars['COMPREFOFFSET'].size)
        parSchema.addField('compRefSigma', type='ArrayD',
                           doc='Width of reference star/calibrated star distribution',
                           size=pars['COMPREFSIGMA'].size)
        parSchema.addField('compAperCorrPivot', type='ArrayD', doc='Aperture correction pivot',
                           size=pars['COMPAPERCORRPIVOT'].size)
        parSchema.addField('compAperCorrSlope', type='ArrayD', doc='Aperture correction slope',
                           size=pars['COMPAPERCORRSLOPE'].size)
        parSchema.addField('compAperCorrSlopeErr', type='ArrayD', doc='Aperture correction slope error',
                           size=pars['COMPAPERCORRSLOPEERR'].size)
        parSchema.addField('compAperCorrRange', type='ArrayD', doc='Aperture correction range',
                           size=pars['COMPAPERCORRRANGE'].size)
        parSchema.addField('compModelErrExptimePivot', type='ArrayD', doc='Model error exptime pivot',
                           size=pars['COMPMODELERREXPTIMEPIVOT'].size)
        parSchema.addField('compModelErrFwhmPivot', type='ArrayD', doc='Model error fwhm pivot',
                           size=pars['COMPMODELERRFWHMPIVOT'].size)
        parSchema.addField('compModelErrSkyPivot', type='ArrayD', doc='Model error sky pivot',
                           size=pars['COMPMODELERRSKYPIVOT'].size)
        parSchema.addField('compModelErrPars', type='ArrayD', doc='Model error parameters',
                           size=pars['COMPMODELERRPARS'].size)
        parSchema.addField('compExpGray', type='ArrayD', doc='Computed exposure gray',
                           size=pars['COMPEXPGRAY'].size)
        parSchema.addField('compVarGray', type='ArrayD', doc='Computed exposure variance',
                           size=pars['COMPVARGRAY'].size)
        parSchema.addField('compNGoodStarPerExp', type='ArrayI',
                           doc='Computed number of good stars per exposure',
                           size=pars['COMPNGOODSTARPEREXP'].size)
        parSchema.addField('compSigFgcm', type='ArrayD', doc='Computed sigma_fgcm (intrinsic repeatability)',
                           size=pars['COMPSIGFGCM'].size)
        parSchema.addField('compSigmaCal', type='ArrayD', doc='Computed sigma_cal (systematic error floor)',
                           size=pars['COMPSIGMACAL'].size)
        parSchema.addField('compRetrievedLnPwv', type='ArrayD', doc='Retrieved ln(pwv) (smoothed)',
                           size=pars['COMPRETRIEVEDLNPWV'].size)
        parSchema.addField('compRetrievedLnPwvRaw', type='ArrayD', doc='Retrieved ln(pwv) (raw)',
                           size=pars['COMPRETRIEVEDLNPWVRAW'].size)
        parSchema.addField('compRetrievedLnPwvFlag', type='ArrayI', doc='Retrieved ln(pwv) Flag',
                           size=pars['COMPRETRIEVEDLNPWVFLAG'].size)
        parSchema.addField('compRetrievedTauNight', type='ArrayD', doc='Retrieved tau (per night)',
                           size=pars['COMPRETRIEVEDTAUNIGHT'].size)
        # superstarflat section
        parSchema.addField('superstarSize', type='ArrayI', doc='Superstar matrix size',
                           size=4)
        parSchema.addField('superstar', type='ArrayD', doc='Superstar matrix (flattened)',
                           size=parSuperStarFlat.size)

        return parSchema

    def _makeParCatalog(self, parSchema, parInfo, pars, parSuperStarFlat,
                        lutFilterNameString, fitBandString, notFitBandString):
        """
        Make the FGCM parameter catalog for persistence

        Parameters
        ----------
        parSchema: `afwTable.schema`
           Parameter catalog schema
        pars: `np.ndarray`
           FGCM parameters to put into parCat
        parSuperStarFlat: `np.array`
           FGCM superstar flat array to put into parCat
        lutFilterNameString: `str`
           Combined string of all the lutFilterNames
        fitBandString: `str`
           Combined string of all the fitBands
        notFitBandString: `str`
           Combined string of all the bands not used in the fit

        Returns
        -------
        parCat: `afwTable.BasicCatalog`
           Atmosphere and instrumental model parameter catalog for persistence
        """

        parCat = afwTable.BaseCatalog(parSchema)
        parCat.reserve(1)

        # The parameter catalog just has one row, with many columns for all the
        # atmosphere and instrument fit parameters
        rec = parCat.addNew()

        # info section
        rec['nCcd'] = parInfo['NCCD']
        rec['lutFilterNames'] = lutFilterNameString
        rec['fitBands'] = fitBandString
        rec['notFitBands'] = notFitBandString
        rec['lnTauUnit'] = parInfo['LNTAUUNIT']
        rec['lnTauSlopeUnit'] = parInfo['LNTAUSLOPEUNIT']
        rec['alphaUnit'] = parInfo['ALPHAUNIT']
        rec['lnPwvUnit'] = parInfo['LNPWVUNIT']
        rec['lnPwvSlopeUnit'] = parInfo['LNPWVSLOPEUNIT']
        rec['lnPwvQuadraticUnit'] = parInfo['LNPWVQUADRATICUNIT']
        rec['lnPwvGlobalUnit'] = parInfo['LNPWVGLOBALUNIT']
        rec['o3Unit'] = parInfo['O3UNIT']
        rec['qeSysUnit'] = parInfo['QESYSUNIT']
        rec['qeSysSlopeUnit'] = parInfo['QESYSSLOPEUNIT']
        rec['filterOffsetUnit'] = parInfo['FILTEROFFSETUNIT']
        # note these are not currently supported here.
        rec['hasExternalPwv'] = 0
        rec['hasExternalTau'] = 0

        # parameter section

        scalarNames = ['parRetrievedLnPwvScale', 'parRetrievedLnPwvOffset']

        arrNames = ['parAlpha', 'parO3', 'parLnTauIntercept', 'parLnTauSlope',
                    'parLnPwvIntercept', 'parLnPwvSlope', 'parLnPwvQuadratic',
                    'parQeSysIntercept',
                    'parQeSysSlope', 'parRetrievedLnPwvNightlyOffset', 'compAperCorrPivot',
                    'parFilterOffset', 'parFilterOffsetFitFlag',
                    'compAbsThroughput', 'compRefOffset', 'compRefSigma',
                    'compAperCorrSlope', 'compAperCorrSlopeErr', 'compAperCorrRange',
                    'compModelErrExptimePivot', 'compModelErrFwhmPivot',
                    'compModelErrSkyPivot', 'compModelErrPars',
                    'compExpGray', 'compVarGray', 'compNGoodStarPerExp', 'compSigFgcm',
                    'compSigmaCal',
                    'compRetrievedLnPwv', 'compRetrievedLnPwvRaw', 'compRetrievedLnPwvFlag',
                    'compRetrievedTauNight']

        for scalarName in scalarNames:
            rec[scalarName] = pars[scalarName.upper()]

        for arrName in arrNames:
            rec[arrName][:] = np.atleast_1d(pars[0][arrName.upper()])[:]

        # superstar section
        rec['superstarSize'][:] = parSuperStarFlat.shape
        rec['superstar'][:] = parSuperStarFlat.flatten()

        return parCat

    def _makeFlagStarSchema(self):
        """
        Make the flagged-stars schema

        Returns
        -------
        flagStarSchema: `afwTable.schema`
        """

        flagStarSchema = afwTable.Schema()

        flagStarSchema.addField('objId', type=np.int32, doc='FGCM object id')
        flagStarSchema.addField('objFlag', type=np.int32, doc='FGCM object flag')

        return flagStarSchema

    def _makeFlagStarCat(self, flagStarSchema, flagStarStruct):
        """
        Make the flagged star catalog for persistence

        Parameters
        ----------
        flagStarSchema: `afwTable.schema`
           Flagged star schema
        flagStarStruct: `np.ndarray`
           Flagged star structure from fgcm

        Returns
        -------
        flagStarCat: `afwTable.BaseCatalog`
           Flagged star catalog for persistence
        """

        flagStarCat = afwTable.BaseCatalog(flagStarSchema)
        flagStarCat.reserve(flagStarStruct.size)
        for i in range(flagStarStruct.size):
            flagStarCat.addNew()

        flagStarCat['objId'][:] = flagStarStruct['OBJID']
        flagStarCat['objFlag'][:] = flagStarStruct['OBJFLAG']

        return flagStarCat

    def _makeZptSchema(self, chebyshevSize):
        """
        Make the zeropoint schema

        Parameters
        ----------
        chebyshevSize: `int`
           Length of the zeropoint chebyshev array

        Returns
        -------
        zptSchema: `afwTable.schema`
        """

        zptSchema = afwTable.Schema()

        zptSchema.addField('visit', type=np.int32, doc='Visit number')
        zptSchema.addField('ccd', type=np.int32, doc='CCD number')
        zptSchema.addField('fgcmFlag', type=np.int32, doc=('FGCM flag value: '
                                                           '1: Photometric, used in fit; '
                                                           '2: Photometric, not used in fit; '
                                                           '4: Non-photometric, on partly photometric night; '
                                                           '8: Non-photometric, on non-photometric night; '
                                                           '16: No zeropoint could be determined; '
                                                           '32: Too few stars for reliable gray computation'))
        zptSchema.addField('fgcmZpt', type=np.float32, doc='FGCM zeropoint (center of CCD)')
        zptSchema.addField('fgcmZptErr', type=np.float32,
                           doc='Error on zeropoint, estimated from repeatability + number of obs')
        if self.config.superStarSubCcd:
            zptSchema.addField('fgcmfZptCheb', type='ArrayD',
                               size=chebyshevSize,
                               doc='Chebyshev parameters (flattened) for zeropoint')
            zptSchema.addField('fgcmfZptChebXyMax', type='ArrayD', size=2,
                               doc='maximum x/maximum y to scale to apply chebyshev parameters')
        zptSchema.addField('fgcmI0', type=np.float32, doc='Integral of the passband')
        zptSchema.addField('fgcmI10', type=np.float32, doc='Normalized chromatic integral')
        zptSchema.addField('fgcmR0', type=np.float32,
                           doc='Retrieved i0 integral, estimated from stars (only for flag 1)')
        zptSchema.addField('fgcmR10', type=np.float32,
                           doc='Retrieved i10 integral, estimated from stars (only for flag 1)')
        zptSchema.addField('fgcmGry', type=np.float32,
                           doc='Estimated gray extinction relative to atmospheric solution; '
                           'only for flag <= 4')
        zptSchema.addField('fgcmZptVar', type=np.float32, doc='Variance of zeropoint over ccd')
        zptSchema.addField('fgcmTilings', type=np.float32,
                           doc='Number of photometric tilings used for solution for ccd')
        zptSchema.addField('fgcmFpGry', type=np.float32,
                           doc='Average gray extinction over the full focal plane '
                           '(same for all ccds in a visit)')
        zptSchema.addField('fgcmFpVar', type=np.float32,
                           doc='Variance of gray extinction over the full focal plane '
                           '(same for all ccds in a visit)')
        zptSchema.addField('fgcmDust', type=np.float32,
                           doc='Gray dust extinction from the primary/corrector'
                           'at the time of the exposure')
        zptSchema.addField('fgcmFlat', type=np.float32, doc='Superstarflat illumination correction')
        zptSchema.addField('fgcmAperCorr', type=np.float32, doc='Aperture correction estimated by fgcm')
        zptSchema.addField('exptime', type=np.float32, doc='Exposure time')
        zptSchema.addField('filtername', type=str, size=2, doc='Filter name')

        return zptSchema

    def _makeZptCat(self, zptSchema, zpStruct):
        """
        Make the zeropoint catalog for persistence

        Parameters
        ----------
        zptSchema: `afwTable.schema`
           Zeropoint catalog schema
        zpStruct: `np.ndarray`
           Zeropoint structure from fgcm

        Returns
        -------
        zptCat: `afwTable.BaseCatalog`
           Zeropoint catalog for persistence
        """

        zptCat = afwTable.BaseCatalog(zptSchema)
        zptCat.reserve(zpStruct.size)

        for filterName in zpStruct['FILTERNAME']:
            rec = zptCat.addNew()
            rec['filtername'] = filterName.decode('utf-8')

        zptCat['visit'][:] = zpStruct['VISIT']
        zptCat['ccd'][:] = zpStruct['CCD']
        zptCat['fgcmFlag'][:] = zpStruct['FGCM_FLAG']
        zptCat['fgcmZpt'][:] = zpStruct['FGCM_ZPT']
        zptCat['fgcmZptErr'][:] = zpStruct['FGCM_ZPTERR']
        if self.config.superStarSubCcd:
            zptCat['fgcmfZptCheb'][:, :] = zpStruct['FGCM_FZPT_CHEB']
            zptCat['fgcmfZptChebXyMax'][:, :] = zpStruct['FGCM_FZPT_CHEB_XYMAX']
        zptCat['fgcmI0'][:] = zpStruct['FGCM_I0']
        zptCat['fgcmI10'][:] = zpStruct['FGCM_I10']
        zptCat['fgcmR0'][:] = zpStruct['FGCM_R0']
        zptCat['fgcmR10'][:] = zpStruct['FGCM_R10']
        zptCat['fgcmGry'][:] = zpStruct['FGCM_GRY']
        zptCat['fgcmZptVar'][:] = zpStruct['FGCM_ZPTVAR']
        zptCat['fgcmTilings'][:] = zpStruct['FGCM_TILINGS']
        zptCat['fgcmFpGry'][:] = zpStruct['FGCM_FPGRY']
        zptCat['fgcmFpVar'][:] = zpStruct['FGCM_FPVAR']
        zptCat['fgcmDust'][:] = zpStruct['FGCM_DUST']
        zptCat['fgcmFlat'][:] = zpStruct['FGCM_FLAT']
        zptCat['fgcmAperCorr'][:] = zpStruct['FGCM_APERCORR']
        zptCat['exptime'][:] = zpStruct['EXPTIME']

        return zptCat

    def _makeAtmSchema(self):
        """
        Make the atmosphere schema

        Returns
        -------
        atmSchema: `afwTable.schema`
        """

        atmSchema = afwTable.Schema()

        atmSchema.addField('visit', type=np.int32, doc='Visit number')
        atmSchema.addField('pmb', type=np.float64, doc='Barometric pressure (mb)')
        atmSchema.addField('pwv', type=np.float64, doc='Water vapor (mm)')
        atmSchema.addField('tau', type=np.float64, doc='Aerosol optical depth')
        atmSchema.addField('alpha', type=np.float64, doc='Aerosol slope')
        atmSchema.addField('o3', type=np.float64, doc='Ozone (dobson)')
        atmSchema.addField('secZenith', type=np.float64, doc='Secant(zenith) (~ airmass)')

        return atmSchema

    def _makeAtmCat(self, atmSchema, atmStruct):
        """
        Make the atmosphere catalog for persistence

        Parameters
        ----------
        atmSchema: `afwTable.schema`
           Atmosphere catalog schema
        atmStruct: `np.ndarray`
           Atmosphere structure from fgcm

        Returns
        -------
        atmCat: `afwTable.BaseCatalog`
           Atmosphere catalog for persistence
        """

        atmCat = afwTable.BaseCatalog(atmSchema)
        atmCat.reserve(atmStruct.size)
        for i in range(atmStruct.size):
            atmCat.addNew()

        atmCat['visit'][:] = atmStruct['VISIT']
        atmCat['pmb'][:] = atmStruct['PMB']
        atmCat['pwv'][:] = atmStruct['PWV']
        atmCat['tau'][:] = atmStruct['TAU']
        atmCat['alpha'][:] = atmStruct['ALPHA']
        atmCat['o3'][:] = atmStruct['O3']
        atmCat['secZenith'][:] = atmStruct['SECZENITH']

        return atmCat

    def _makeStdSchema(self):
        """
        Make the standard star schema

        Returns
        -------
        stdSchema: `afwTable.schema`
        """

        stdSchema = afwTable.SimpleTable.makeMinimalSchema()
        stdSchema.addField('ngood', type='ArrayI', doc='Number of good observations',
                           size=len(self.config.bands))
        stdSchema.addField('mag_std_noabs', type='ArrayF',
                           doc='Standard magnitude (no absolute calibration)',
                           size=len(self.config.bands))
        stdSchema.addField('magErr_std', type='ArrayF',
                           doc='Standard magnitude error',
                           size=len(self.config.bands))

        return stdSchema

    def _makeStdCat(self, stdSchema, stdStruct):
        """
        Make the standard star catalog for persistence

        Parameters
        ----------
        stdSchema: `afwTable.schema`
           Standard star catalog schema
        stdStruct: `np.ndarray`
           Standard star structure in FGCM format

        Returns
        -------
        stdCat: `afwTable.BaseCatalog`
           Standard star catalog for persistence
        """

        stdCat = afwTable.SimpleCatalog(stdSchema)

        stdCat.reserve(stdStruct.size)
        for i in range(stdStruct.size):
            stdCat.addNew()

        stdCat['id'][:] = stdStruct['FGCM_ID']
        stdCat['coord_ra'][:] = stdStruct['RA'] * lsst.geom.degrees
        stdCat['coord_dec'][:] = stdStruct['DEC'] * lsst.geom.degrees
        stdCat['ngood'][:, :] = stdStruct['NGOOD'][:, :]
        stdCat['mag_std_noabs'][:, :] = stdStruct['MAG_STD'][:, :]
        stdCat['magErr_std'][:, :] = stdStruct['MAGERR_STD'][:, :]

        return stdCat
