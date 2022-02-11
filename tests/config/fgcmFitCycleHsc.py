# All camera defaults were copied from obs_subaru/config/fgcmFitCycle.py
# on 07/21/21, weekly w_2021_29.

import lsst.fgcmcal as fgcmcal

config.outfileBase = 'TestFgcm'
# Use these bands and fit them.
# g band does not have any observations in tests, this ensures
# that the code performs properly when there are missing bands.
config.bands = ['g', 'r', 'i']
config.fitBands = ['g', 'r', 'i']
from lsst.obs.hsc.hscFilters import HSC_FILTER_DEFINITIONS
config.physicalFilterMap = HSC_FILTER_DEFINITIONS.physical_to_band
# Only require r, i observations for a star to be a calibration star.
config.requiredBands = ['r', 'i']
# Do 5 iterations in multi-cycle run mode.
config.maxIterBeforeFinalCycle = 5
config.nCore = 1
config.cycleNumber = 0
config.utBoundary = 0.0
config.washMjds = (0.0, )
# For tests, define 1 observing epoch that encompasses everything.
config.epochMjds = (0.0, 100000.0)
config.coatingMjds = [56650.0, 58050.0]
config.latitude = 19.8256
config.defaultCameraOrientation = 270.0
config.expGrayPhotometricCutDict = {'g': -0.1, 'r': -0.1, 'i': -0.1}
config.expGrayHighCutDict = {'g': 0.1, 'r': 0.1, 'i': 0.1}
config.expVarGrayPhotometricCutDict = {'g': 0.05**2.,
                                       'r': 0.05**2.,
                                       'i': 0.05**2.}
# For tests, make a broad cut for outliers.
config.autoPhotometricCutNSig = 5.0
config.autoHighCutNSig = 5.0
# Fit aperture corrections with only 2 bins to exercise the code.
config.aperCorrFitNBins = 2
config.aperCorrInputSlopeDict = {'g': -1.0,
                                 'r': -0.9694,
                                 'i': -1.7229}
# Define the band to SED constants approximately so they work
# for data that only has r, i observations.
config.sedboundaryterms = fgcmcal.SedboundarytermDict()
config.sedboundaryterms.data['ri'] = fgcmcal.Sedboundaryterm(primary='r',
                                                             secondary='i')
config.sedterms = fgcmcal.SedtermDict()
config.sedterms.data['g'] = fgcmcal.Sedterm(primaryTerm='ri', secondaryTerm=None,
                                            extrapolated=False, constant=0.0)
config.sedterms.data['r'] = fgcmcal.Sedterm(primaryTerm='ri', secondaryTerm=None,
                                            extrapolated=False, constant=1.0)
config.sedterms.data['i'] = fgcmcal.Sedterm(primaryTerm='ri', secondaryTerm=None,
                                            extrapolated=False, constant=0.75)
# Define good stars with an r-i color cut.
config.starColorCuts = ('r,i,-0.50,2.25',)
config.precomputeSuperStarInitialCycle = False
config.superStarSubCcdDict = {'g': True,
                              'r': True,
                              'i': True}
# Allow calibration to work with just 1 exposure on a night.
config.minExpPerNight = 1
# Allow calibration to work with very few stars per exposure.
config.minStarPerExp = 50
# Allow calibration to work with small number of stars in processing batches.
config.nStarPerRun = 50
config.nExpPerRun = 2
# Define r-i color as the primary way to split by color.
config.colorSplitBands = ['r', 'i']
config.freezeStdAtmosphere = True
# For tests, do low-order per-ccd polynomial.
config.superStarSubCcdChebyshevOrder = 1
config.ccdGraySubCcdDict = {'g': False,
                            'r': False,
                            'i': False}
config.ccdGrayFocalPlaneDict = {'g': True,
                                'r': True,
                                'i': True}
config.ccdGrayFocalPlaneFitMinCcd = 50
config.ccdGrayFocalPlaneChebyshevOrder = 1
# Do not model the magnitude errors (use errors as reported).
config.modelMagErrors = False
# Fix the sigma_cal calibration noise to 0.003 mag.
config.sigmaCalRange = (0.003, 0.003)
# Do not fit instrumental parameters (mirror decay) per band.
config.instrumentParsPerBand = False
# Set the random seed for repeatability in fits.
config.randomSeed = 12345
# Do not use star repeatability metrics for selecting exposures.
# (Instead, use exposure repeatability metrics).
config.useRepeatabilityForExpGrayCutsDict = {'g': False,
                                             'r': False,
                                             'i': False}
config.sigFgcmMaxEGrayDict = {'g': 0.05,
                              'r': 0.05,
                              'i': 0.05}
config.approxThroughputDict = {'g': 1.0,
                              'r': 1.0,
                              'i': 1.0}

config.deltaAperFitPerCcdNx = 2
config.deltaAperFitPerCcdNy = 4
config.deltaAperInnerRadiusArcsec = 2.04
config.deltaAperOuterRadiusArcsec = 2.89
config.doComputeDeltaAperPerVisit = True
config.doComputeDeltaAperMap = True
config.doComputeDeltaAperPerCcd = True
