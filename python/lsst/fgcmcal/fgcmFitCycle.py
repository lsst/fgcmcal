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
import copy

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from .utilities import makeConfigDict, translateFgcmLut, translateVisitCatalog
from .utilities import extractReferenceMags
from .utilities import computeCcdOffsets, makeZptSchema, makeZptCat
from .utilities import makeAtmSchema, makeAtmCat, makeStdSchema, makeStdCat
from .sedterms import SedboundarytermDict, SedtermDict

import fgcm

__all__ = ['FgcmFitCycleConfig', 'FgcmFitCycleTask', 'FgcmFitCycleRunner']


class FgcmFitCycleConfig(pexConfig.Config):
    """Config for FgcmFitCycle"""

    bands = pexConfig.ListField(
        doc="Bands to run calibration",
        dtype=str,
        default=[],
    )
    fitFlag = pexConfig.ListField(
        doc=("Flag for which bands are directly constrained in the FGCM fit. "
             "Bands set to 0 will have the atmosphere constrained from observations "
             "in other bands on the same night.  Must be same length as config.bands, "
             "and matched band-by-band."),
        dtype=int,
        default=(0,),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use fitBands instead."),
    )
    fitBands = pexConfig.ListField(
        doc=("Bands to use in atmospheric fit. The bands not listed here will have "
             "the atmosphere constrained from the 'fitBands' on the same night. "
             "Must be a subset of `config.bands`"),
        dtype=str,
        default=[],
    )
    requiredFlag = pexConfig.ListField(
        doc=("Flag for which bands are required for a star to be considered a calibration "
             "star in the FGCM fit.  Typically this should be the same as fitFlag.  Must "
             "be same length as config.bands, and matched band-by-band."),
        dtype=int,
        default=(0,),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use requiredBands instead."),
    )
    requiredBands = pexConfig.ListField(
        doc=("Bands that are required for a star to be considered a calibration star. "
             "Must be a subset of `config.bands`"),
        dtype=str,
        default=[],
    )
    filterMap = pexConfig.DictField(
        doc="Mapping from 'filterName' to band.",
        keytype=str,
        itemtype=str,
        default={},
    )
    doReferenceCalibration = pexConfig.Field(
        doc="Use reference catalog as additional constraint on calibration",
        dtype=bool,
        default=True,
    )
    refStarSnMin = pexConfig.Field(
        doc="Reference star signal-to-noise minimum to use in calibration.  Set to <=0 for no cut.",
        dtype=float,
        default=50.0,
    )
    refStarOutlierNSig = pexConfig.Field(
        doc=("Number of sigma compared to average mag for reference star to be considered an outlier. "
             "Computed per-band, and if it is an outlier in any band it is rejected from fits."),
        dtype=float,
        default=4.0,
    )
    applyRefStarColorCuts = pexConfig.Field(
        doc="Apply color cuts to reference stars?",
        dtype=bool,
        default=True,
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
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use superStarSubCcdDict instead."),
    )
    superStarSubCcdDict = pexConfig.DictField(
        doc=("Per-band specification on whether to compute superstar flat on sub-ccd scale. "
             "Must have one entry per band."),
        keytype=str,
        itemtype=bool,
        default={},
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
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use ccdGraySubCcdDict instead."),
    )
    ccdGraySubCcdDict = pexConfig.DictField(
        doc=("Per-band specification on whether to compute achromatic per-ccd residual "
             "('ccd gray') on a sub-ccd scale."),
        keytype=str,
        itemtype=bool,
        default={},
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
             "Must be same length as config.bands, and matched band-by-band."),
        dtype=float,
        default=(0.0,),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use expGrayPhotometricCutDict instead."),
    )
    expGrayPhotometricCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum (negative) achromatic exposure residual "
             "('gray term') for a visit to be considered photometric.  Must have one "
             "entry per band. Broad-band filters should be -0.05."),
        keytype=str,
        itemtype=float,
        default={},
    )
    expGrayHighCut = pexConfig.ListField(
        doc=("Maximum (positive) exposure gray for a visit to be considered photometric. "
             "Must be same length as config.bands, and matched band-by-band."),
        dtype=float,
        default=(0.0,),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use expGrayHighCutDict instead."),
    )
    expGrayHighCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum (positive) achromatic exposure residual "
             "('gray term') for a visit to be considered photometric.  Must have one "
             "entry per band.  Broad-band filters should be 0.2."),
        keytype=str,
        itemtype=float,
        default={},
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
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use expVarGrayPhotometricCutDict instead."),
    )
    expVarGrayPhotometricCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum exposure variance to be considered possibly "
             "photometric.  Must have one entry per band.  Broad-band filters should be "
             "0.0005."),
        keytype=str,
        itemtype=float,
        default={},
    )
    expGrayErrRecoverCut = pexConfig.Field(
        doc=("Maximum exposure gray error to be able to recover bad ccds via interpolation. "
             "Visits with more gray variance will only get CCD zeropoints if there are "
             "sufficient star observations (minStarPerCcd) on that CCD."),
        dtype=float,
        default=0.05,
    )
    aperCorrFitNBins = pexConfig.Field(
        doc=("Number of aperture bins used in aperture correction fit.  When set to 0"
             "no fit will be performed, and the config.aperCorrInputSlopes will be "
             "used if available."),
        dtype=int,
        default=10,
    )
    aperCorrInputSlopes = pexConfig.ListField(
        doc=("Aperture correction input slope parameters.  These are used on the first "
             "fit iteration, and aperture correction parameters will be updated from "
             "the data if config.aperCorrFitNBins > 0.  It is recommended to set this"
             "when there is insufficient data to fit the parameters (e.g. tract mode). "
             "If set, must be same length as config.bands, and matched band-by-band."),
        dtype=float,
        default=[],
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use aperCorrInputSlopeDict instead."),
    )
    aperCorrInputSlopeDict = pexConfig.DictField(
        doc=("Per-band specification of aperture correction input slope parameters.  These "
             "are used on the first fit iteration, and aperture correction parameters will "
             "be updated from the data if config.aperCorrFitNBins > 0.  It is recommended "
             "to set this when there is insufficient data to fit the parameters (e.g. "
             "tract mode)."),
        keytype=str,
        itemtype=float,
        default={},
    )
    sedFudgeFactors = pexConfig.ListField(
        doc=("Fudge factors for computing linear SED from colors.  Must be same length as "
             "config.bands, and matched band-by-band."),
        dtype=float,
        default=(0,),
        optional=True,
        deprecated=("This field has been deprecated and will be removed after v20. "
                    "Please use sedSlopeTermMap and sedSlopeMap."),
    )
    sedboundaryterms = pexConfig.ConfigField(
        doc="Mapping from bands to SED boundary term names used is sedterms.",
        dtype=SedboundarytermDict,
    )
    sedterms = pexConfig.ConfigField(
        doc="Mapping from terms to bands for fgcm linear SED approximations.",
        dtype=SedtermDict,
    )
    sigFgcmMaxErr = pexConfig.Field(
        doc="Maximum mag error for fitting sigma_FGCM",
        dtype=float,
        default=0.01,
    )
    sigFgcmMaxEGray = pexConfig.ListField(
        doc=("Maximum (absolute) gray value for observation in sigma_FGCM. "
             "May be 1 element (same for all bands) or the same length as config.bands."),
        dtype=float,
        default=(0.05,),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use sigFgcmMaxEGrayDict instead."),
    )
    sigFgcmMaxEGrayDict = pexConfig.DictField(
        doc=("Per-band specification for maximum (absolute) achromatic residual (gray value) "
             "for observations in sigma_fgcm (raw repeatability).  Broad-band filters "
             "should be 0.05."),
        keytype=str,
        itemtype=float,
        default={},
    )
    ccdGrayMaxStarErr = pexConfig.Field(
        doc=("Maximum error on a star observation to use in ccd gray (achromatic residual) "
             "computation"),
        dtype=float,
        default=0.10,
    )
    approxThroughput = pexConfig.ListField(
        doc=("Approximate overall throughput at start of calibration observations. "
             "May be 1 element (same for all bands) or the same length as config.bands."),
        dtype=float,
        default=(1.0, ),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use approxThroughputDict instead."),
    )
    approxThroughputDict = pexConfig.DictField(
        doc=("Per-band specification of the approximate overall throughput at the start of "
             "calibration observations.  Must have one entry per band.  Typically should "
             "be 1.0."),
        keytype=str,
        itemtype=float,
        default={},
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
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use colorSplitBands instead."),
    )
    colorSplitBands = pexConfig.ListField(
        doc="Band names to use to split stars by color.  Must have 2 entries.",
        dtype=str,
        length=2,
        default=('g', 'i'),
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
    instrumentParsPerBand = pexConfig.Field(
        doc=("Model instrumental parameters per band? "
             "Otherwise, instrumental parameters (QE changes with time) are "
             "shared among all bands."),
        dtype=bool,
        default=False,
    )
    instrumentSlopeMinDeltaT = pexConfig.Field(
        doc=("Minimum time change (in days) between observations to use in constraining "
             "instrument slope."),
        dtype=float,
        default=20.0,
    )
    fitMirrorChromaticity = pexConfig.Field(
        doc="Fit (intraband) mirror chromatic term?",
        dtype=bool,
        default=False,
    )
    coatingMjds = pexConfig.ListField(
        doc="Mirror coating dates in MJD",
        dtype=float,
        default=(0.0,),
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
    useRepeatabilityForExpGrayCuts = pexConfig.ListField(
        doc=("Use star repeatability (instead of exposures) for computing photometric "
             "cuts? Recommended for tract mode or bands with few exposures. "
             "May be 1 element (same for all bands) or the same length as config.bands."),
        dtype=bool,
        default=(False,),
        optional=True,
        deprecated=("This field is no longer used, and has been deprecated by DM-23699. "
                    "It will be removed after v20.  Use useRepeatabilityForExpGrayCutsDict instead."),
    )
    useRepeatabilityForExpGrayCutsDict = pexConfig.DictField(
        doc=("Per-band specification on whether to use star repeatability (instead of exposures) "
             "for computing photometric cuts. Recommended for tract mode or bands with few visits."),
        keytype=str,
        itemtype=bool,
        default={},
    )
    autoPhotometricCutNSig = pexConfig.Field(
        doc=("Number of sigma for automatic computation of (low) photometric cut. "
             "Cut is based on exposure gray width (per band), unless "
             "useRepeatabilityForExpGrayCuts is set, in which case the star "
             "repeatability is used (also per band)."),
        dtype=float,
        default=3.0,
    )
    autoHighCutNSig = pexConfig.Field(
        doc=("Number of sigma for automatic computation of (high) outlier cut. "
             "Cut is based on exposure gray width (per band), unless "
             "useRepeatabilityForExpGrayCuts is set, in which case the star "
             "repeatability is used (also per band)."),
        dtype=float,
        default=4.0,
    )
    quietMode = pexConfig.Field(
        doc="Be less verbose with logging.",
        dtype=bool,
        default=False,
    )

    def setDefaults(self):
        pass

    def validate(self):
        super().validate()

        for band in self.fitBands:
            if band not in self.bands:
                msg = 'fitBand %s not in bands' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.fitBands, self, msg)
        for band in self.requiredBands:
            if band not in self.bands:
                msg = 'requiredBand %s not in bands' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.requiredBands, self, msg)
        for band in self.colorSplitBands:
            if band not in self.bands:
                msg = 'colorSplitBand %s not in bands' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.colorSplitBands, self, msg)
        for band in self.bands:
            if band not in self.superStarSubCcdDict:
                msg = 'band %s not in superStarSubCcdDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.superStarSubCcdDict,
                                                     self, msg)
            if band not in self.ccdGraySubCcdDict:
                msg = 'band %s not in ccdGraySubCcdDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.ccdGraySubCcdDict,
                                                     self, msg)
            if band not in self.expGrayPhotometricCutDict:
                msg = 'band %s not in expGrayPhotometricCutDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.expGrayPhotometricCutDict,
                                                     self, msg)
            if band not in self.expGrayHighCutDict:
                msg = 'band %s not in expGrayHighCutDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.expGrayHighCutDict,
                                                     self, msg)
            if band not in self.expVarGrayPhotometricCutDict:
                msg = 'band %s not in expVarGrayPhotometricCutDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.expVarGrayPhotometricCutDict,
                                                     self, msg)
            if band not in self.sigFgcmMaxEGrayDict:
                msg = 'band %s not in sigFgcmMaxEGrayDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.sigFgcmMaxEGrayDict,
                                                     self, msg)
            if band not in self.approxThroughputDict:
                msg = 'band %s not in approxThroughputDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.approxThroughputDict,
                                                     self, msg)
            if band not in self.useRepeatabilityForExpGrayCutsDict:
                msg = 'band %s not in useRepeatabilityForExpGrayCutsDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.useRepeatabilityForExpGrayCutsDict,
                                                     self, msg)


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
        configDict = makeConfigDict(self.config, self.log, camera,
                                    self.maxIter, self.resetFitParameters,
                                    self.outputZeropoints)

        lutCat = butler.get('fgcmLookUpTable')
        fgcmLut, lutIndexVals, lutStd = translateFgcmLut(lutCat, dict(self.config.filterMap))
        del lutCat

        # next we need the exposure/visit information

        # fgcmExpInfo = self._loadVisitCatalog(butler)
        visitCat = butler.get('fgcmVisitCatalog')
        fgcmExpInfo = translateVisitCatalog(visitCat)
        del visitCat

        # Use the first orientation.
        # TODO: DM-21215 will generalize to arbitrary camera orientations
        ccdOffsets = computeCcdOffsets(camera, fgcmExpInfo['TELROT'][0])

        noFitsDict = {'lutIndex': lutIndexVals,
                      'lutStd': lutStd,
                      'expInfo': fgcmExpInfo,
                      'ccdOffsets': ccdOffsets}

        # set up the fitter object
        fgcmFitCycle = fgcm.FgcmFitCycle(configDict, useFits=False,
                                         noFitsDict=noFitsDict, noOutput=True)

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
            flaggedStars = None
            flagId = None
            flagFlag = None

        if self.config.doReferenceCalibration:
            refStars = butler.get('fgcmReferenceStars')

            refMag, refMagErr = extractReferenceMags(refStars,
                                                     self.config.bands,
                                                     self.config.filterMap)
            refId = refStars['fgcm_id'][:]
        else:
            refStars = None
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
                            psfCandidate=starObs['psf_candidate'][starIndices['obsIndex']],
                            refID=refId,
                            refMag=refMag,
                            refMagErr=refMagErr,
                            flagID=flagId,
                            flagFlag=flagFlag,
                            computeNobs=True)

        # Release all references to temporary objects holding star data (see above)
        del starObs
        del starIds
        del starIndices
        del flagId
        del flagFlag
        del flaggedStars
        del refStars
        del refId
        del refMag
        del refMagErr

        # and set the bits in the cycle object
        fgcmFitCycle.setLUT(fgcmLut)
        fgcmFitCycle.setStars(fgcmStars, fgcmPars)
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

        updatedPhotometricCutDict = {b: float(fgcmFitCycle.updatedPhotometricCut[i]) for
                                     i, b in enumerate(self.config.bands)}
        updatedHighCutDict = {band: float(fgcmFitCycle.updatedHighCut[i]) for
                              i, band in enumerate(self.config.bands)}

        outConfig = copy.copy(self.config)
        outConfig.update(cycleNumber=(self.config.cycleNumber + 1),
                         precomputeSuperStarInitialCycle=False,
                         freezeStdAtmosphere=False,
                         expGrayPhotometricCutDict=updatedPhotometricCutDict,
                         expGrayHighCutDict=updatedHighCutDict)
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

    def _loadParameters(self, butler):
        """
        Load FGCM parameters from a previous fit cycle

        Parameters
        ----------
        butler:  `lsst.daf.persistence.Butler`

        Returns
        -------
        inParInfo: `numpy.ndarray`
           Numpy array parameter information formatted for input to fgcm
        inParameters: `numpy.ndarray`
           Numpy array parameter values formatted for input to fgcm
        inSuperStar: `numpy.array`
           Superstar flat formatted for input to fgcm
        """

        # note that we already checked that this is available
        parCat = butler.get('fgcmFitParameters', fgcmcycle=self.config.cycleNumber-1)

        parLutFilterNames = np.array(parCat[0]['lutFilterNames'].split(','))
        parFitBands = np.array(parCat[0]['fitBands'].split(','))

        inParInfo = np.zeros(1, dtype=[('NCCD', 'i4'),
                                       ('LUTFILTERNAMES', parLutFilterNames.dtype.str,
                                        (parLutFilterNames.size, )),
                                       ('FITBANDS', parFitBands.dtype.str, (parFitBands.size, )),
                                       ('LNTAUUNIT', 'f8'),
                                       ('LNTAUSLOPEUNIT', 'f8'),
                                       ('ALPHAUNIT', 'f8'),
                                       ('LNPWVUNIT', 'f8'),
                                       ('LNPWVSLOPEUNIT', 'f8'),
                                       ('LNPWVQUADRATICUNIT', 'f8'),
                                       ('LNPWVGLOBALUNIT', 'f8'),
                                       ('O3UNIT', 'f8'),
                                       ('QESYSUNIT', 'f8'),
                                       ('FILTEROFFSETUNIT', 'f8'),
                                       ('HASEXTERNALPWV', 'i2'),
                                       ('HASEXTERNALTAU', 'i2')])
        inParInfo['NCCD'] = parCat['nCcd']
        inParInfo['LUTFILTERNAMES'][:] = parLutFilterNames
        inParInfo['FITBANDS'][:] = parFitBands
        inParInfo['HASEXTERNALPWV'] = parCat['hasExternalPwv']
        inParInfo['HASEXTERNALTAU'] = parCat['hasExternalTau']

        inParams = np.zeros(1, dtype=[('PARALPHA', 'f8', (parCat['parAlpha'].size, )),
                                      ('PARO3', 'f8', (parCat['parO3'].size, )),
                                      ('PARLNTAUINTERCEPT', 'f8',
                                       (parCat['parLnTauIntercept'].size, )),
                                      ('PARLNTAUSLOPE', 'f8',
                                       (parCat['parLnTauSlope'].size, )),
                                      ('PARLNPWVINTERCEPT', 'f8',
                                       (parCat['parLnPwvIntercept'].size, )),
                                      ('PARLNPWVSLOPE', 'f8',
                                       (parCat['parLnPwvSlope'].size, )),
                                      ('PARLNPWVQUADRATIC', 'f8',
                                       (parCat['parLnPwvQuadratic'].size, )),
                                      ('PARQESYSINTERCEPT', 'f8',
                                       (parCat['parQeSysIntercept'].size, )),
                                      ('COMPQESYSSLOPE', 'f8',
                                       (parCat['compQeSysSlope'].size, )),
                                      ('PARFILTEROFFSET', 'f8',
                                       (parCat['parFilterOffset'].size, )),
                                      ('PARFILTEROFFSETFITFLAG', 'i2',
                                       (parCat['parFilterOffsetFitFlag'].size, )),
                                      ('PARRETRIEVEDLNPWVSCALE', 'f8'),
                                      ('PARRETRIEVEDLNPWVOFFSET', 'f8'),
                                      ('PARRETRIEVEDLNPWVNIGHTLYOFFSET', 'f8',
                                       (parCat['parRetrievedLnPwvNightlyOffset'].size, )),
                                      ('COMPABSTHROUGHPUT', 'f8',
                                       (parCat['compAbsThroughput'].size, )),
                                      ('COMPREFOFFSET', 'f8',
                                       (parCat['compRefOffset'].size, )),
                                      ('COMPREFSIGMA', 'f8',
                                       (parCat['compRefSigma'].size, )),
                                      ('COMPMIRRORCHROMATICITY', 'f8',
                                       (parCat['compMirrorChromaticity'].size, )),
                                      ('MIRRORCHROMATICITYPIVOT', 'f8',
                                       (parCat['mirrorChromaticityPivot'].size, )),
                                      ('COMPAPERCORRPIVOT', 'f8',
                                       (parCat['compAperCorrPivot'].size, )),
                                      ('COMPAPERCORRSLOPE', 'f8',
                                       (parCat['compAperCorrSlope'].size, )),
                                      ('COMPAPERCORRSLOPEERR', 'f8',
                                       (parCat['compAperCorrSlopeErr'].size, )),
                                      ('COMPAPERCORRRANGE', 'f8',
                                       (parCat['compAperCorrRange'].size, )),
                                      ('COMPMODELERREXPTIMEPIVOT', 'f8',
                                       (parCat['compModelErrExptimePivot'].size, )),
                                      ('COMPMODELERRFWHMPIVOT', 'f8',
                                       (parCat['compModelErrFwhmPivot'].size, )),
                                      ('COMPMODELERRSKYPIVOT', 'f8',
                                       (parCat['compModelErrSkyPivot'].size, )),
                                      ('COMPMODELERRPARS', 'f8',
                                       (parCat['compModelErrPars'].size, )),
                                      ('COMPEXPGRAY', 'f8',
                                       (parCat['compExpGray'].size, )),
                                      ('COMPVARGRAY', 'f8',
                                       (parCat['compVarGray'].size, )),
                                      ('COMPNGOODSTARPEREXP', 'i4',
                                       (parCat['compNGoodStarPerExp'].size, )),
                                      ('COMPSIGFGCM', 'f8',
                                       (parCat['compSigFgcm'].size, )),
                                      ('COMPSIGMACAL', 'f8',
                                       (parCat['compSigmaCal'].size, )),
                                      ('COMPRETRIEVEDLNPWV', 'f8',
                                       (parCat['compRetrievedLnPwv'].size, )),
                                      ('COMPRETRIEVEDLNPWVRAW', 'f8',
                                       (parCat['compRetrievedLnPwvRaw'].size, )),
                                      ('COMPRETRIEVEDLNPWVFLAG', 'i2',
                                       (parCat['compRetrievedLnPwvFlag'].size, )),
                                      ('COMPRETRIEVEDTAUNIGHT', 'f8',
                                       (parCat['compRetrievedTauNight'].size, ))])

        inParams['PARALPHA'][:] = parCat['parAlpha'][0, :]
        inParams['PARO3'][:] = parCat['parO3'][0, :]
        inParams['PARLNTAUINTERCEPT'][:] = parCat['parLnTauIntercept'][0, :]
        inParams['PARLNTAUSLOPE'][:] = parCat['parLnTauSlope'][0, :]
        inParams['PARLNPWVINTERCEPT'][:] = parCat['parLnPwvIntercept'][0, :]
        inParams['PARLNPWVSLOPE'][:] = parCat['parLnPwvSlope'][0, :]
        inParams['PARLNPWVQUADRATIC'][:] = parCat['parLnPwvQuadratic'][0, :]
        inParams['PARQESYSINTERCEPT'][:] = parCat['parQeSysIntercept'][0, :]
        inParams['COMPQESYSSLOPE'][:] = parCat['compQeSysSlope'][0, :]
        inParams['PARFILTEROFFSET'][:] = parCat['parFilterOffset'][0, :]
        inParams['PARFILTEROFFSETFITFLAG'][:] = parCat['parFilterOffsetFitFlag'][0, :]
        inParams['PARRETRIEVEDLNPWVSCALE'] = parCat['parRetrievedLnPwvScale']
        inParams['PARRETRIEVEDLNPWVOFFSET'] = parCat['parRetrievedLnPwvOffset']
        inParams['PARRETRIEVEDLNPWVNIGHTLYOFFSET'][:] = parCat['parRetrievedLnPwvNightlyOffset'][0, :]
        inParams['COMPABSTHROUGHPUT'][:] = parCat['compAbsThroughput'][0, :]
        inParams['COMPREFOFFSET'][:] = parCat['compRefOffset'][0, :]
        inParams['COMPREFSIGMA'][:] = parCat['compRefSigma'][0, :]
        inParams['COMPMIRRORCHROMATICITY'][:] = parCat['compMirrorChromaticity'][0, :]
        inParams['MIRRORCHROMATICITYPIVOT'][:] = parCat['mirrorChromaticityPivot'][0, :]
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

        parSchema = self._makeParSchema(parInfo, pars, fgcmFitCycle.fgcmPars.parSuperStarFlat,
                                        lutFilterNameString, fitBandString)
        parCat = self._makeParCatalog(parSchema, parInfo, pars,
                                      fgcmFitCycle.fgcmPars.parSuperStarFlat,
                                      lutFilterNameString, fitBandString)

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
            superStarChebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_SSTAR_CHEB'].shape[1]
            zptChebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB'].shape[1]

            zptSchema = makeZptSchema(superStarChebSize, zptChebSize)
            zptCat = makeZptCat(zptSchema, fgcmFitCycle.fgcmZpts.zpStruct)

            butler.put(zptCat, 'fgcmZeropoints', fgcmcycle=self.config.cycleNumber)

            # Save atmosphere values
            # These are generated by the same code that generates zeropoints
            atmSchema = makeAtmSchema()
            atmCat = makeAtmCat(atmSchema, fgcmFitCycle.fgcmZpts.atmStruct)

            butler.put(atmCat, 'fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        # Save the standard stars (if configured)
        if self.outputStandards:
            stdStruct, goodBands = fgcmFitCycle.fgcmStars.retrieveStdStarCatalog(fgcmFitCycle.fgcmPars)
            stdSchema = makeStdSchema(len(goodBands))
            stdCat = makeStdCat(stdSchema, stdStruct, goodBands)

            butler.put(stdCat, 'fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

    def _makeParSchema(self, parInfo, pars, parSuperStarFlat,
                       lutFilterNameString, fitBandString):
        """
        Make the parameter persistence schema

        Parameters
        ----------
        parInfo: `numpy.ndarray`
           Parameter information returned by fgcm
        pars: `numpy.ndarray`
           Parameter values returned by fgcm
        parSuperStarFlat: `numpy.array`
           Superstar flat values returned by fgcm
        lutFilterNameString: `str`
           Combined string of all the lutFilterNames
        fitBandString: `str`
           Combined string of all the fitBands

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
        parSchema.addField('compQeSysSlope', type='ArrayD', doc='Mirror gray slope parameter vector',
                           size=pars[0]['COMPQESYSSLOPE'].size)
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
        parSchema.addField('compMirrorChromaticity', type='ArrayD',
                           doc='Computed mirror chromaticity terms',
                           size=pars['COMPMIRRORCHROMATICITY'].size)
        parSchema.addField('mirrorChromaticityPivot', type='ArrayD',
                           doc='Mirror chromaticity pivot mjd',
                           size=pars['MIRRORCHROMATICITYPIVOT'].size)
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
                        lutFilterNameString, fitBandString):
        """
        Make the FGCM parameter catalog for persistence

        Parameters
        ----------
        parSchema: `lsst.afw.table.Schema`
           Parameter catalog schema
        pars: `numpy.ndarray`
           FGCM parameters to put into parCat
        parSuperStarFlat: `numpy.array`
           FGCM superstar flat array to put into parCat
        lutFilterNameString: `str`
           Combined string of all the lutFilterNames
        fitBandString: `str`
           Combined string of all the fitBands

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
        # note these are not currently supported here.
        rec['hasExternalPwv'] = 0
        rec['hasExternalTau'] = 0

        # parameter section

        scalarNames = ['parRetrievedLnPwvScale', 'parRetrievedLnPwvOffset']

        arrNames = ['parAlpha', 'parO3', 'parLnTauIntercept', 'parLnTauSlope',
                    'parLnPwvIntercept', 'parLnPwvSlope', 'parLnPwvQuadratic',
                    'parQeSysIntercept', 'compQeSysSlope',
                    'parRetrievedLnPwvNightlyOffset', 'compAperCorrPivot',
                    'parFilterOffset', 'parFilterOffsetFitFlag',
                    'compAbsThroughput', 'compRefOffset', 'compRefSigma',
                    'compMirrorChromaticity', 'mirrorChromaticityPivot',
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
        flagStarSchema: `lsst.afw.table.Schema`
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
        flagStarSchema: `lsst.afw.table.Schema`
           Flagged star schema
        flagStarStruct: `numpy.ndarray`
           Flagged star structure from fgcm

        Returns
        -------
        flagStarCat: `lsst.afw.table.BaseCatalog`
           Flagged star catalog for persistence
        """

        flagStarCat = afwTable.BaseCatalog(flagStarSchema)
        flagStarCat.resize(flagStarStruct.size)

        flagStarCat['objId'][:] = flagStarStruct['OBJID']
        flagStarCat['objFlag'][:] = flagStarStruct['OBJFLAG']

        return flagStarCat
