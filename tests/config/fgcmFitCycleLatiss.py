import lsst.fgcmcal as fgcmcal

config.outfileBase = "TestFgcm"
# Use these bands and fit them.
# i band does not have any observations in tests, this ensures
# that the code performs properly when there are missing bands.
config.bands = ["g", "r", "i"]
config.fitBands = ["g", "r", "i"]
from lsst.obs.lsst.filters import LATISS_FILTER_DEFINITIONS
config.physicalFilterMap = LATISS_FILTER_DEFINITIONS.physical_to_band
# Only require g, r observations for a star to be a calibration star.
config.requiredBands = ["g", "r"]
# Do 5 iterations in multi-cycle run mode.
config.maxIterBeforeFinalCycle = 5
config.nCore = 1
config.cycleNumber = 0
config.utBoundary = 0.0
config.washMjds = (0.0, )
# For tests, define 1 observing epoch that encompasses everything.
config.epochMjds = (0.0, 100000.0)
config.coatingMjds = []
config.latitude = -30.2333
# This is pi*(1.2/2.)**2.
config.mirrorArea = 1.13097
config.defaultCameraOrientation = 0.0
config.brightObsGrayMax = 0.5
config.expGrayInitialCut = -0.5
config.expGrayPhotometricCutDict = {"g": -0.5, "r": -0.5, "i": -0.5}
config.expGrayHighCutDict = {"g": 0.2, "r": 0.2, "i": 0.2}
config.expVarGrayPhotometricCutDict = {"g": 0.1**2.,
                                       "r": 0.1**2.,
                                       "i": 0.1**2.}
# For tests, make a broad cut for outliers.
config.autoPhotometricCutNSig = 5.0
config.autoHighCutNSig = 5.0
# Fit aperture corrections with only 2 bins to exercise the code.
config.aperCorrFitNBins = 2
config.aperCorrInputSlopeDict = {"g": -1.0,
                                 "r": -1.0,
                                 "i": -1.0}
# Define the band to SED constants approximately so they work
# for data that only has r, r observations.
config.sedboundaryterms = fgcmcal.SedboundarytermDict()
config.sedboundaryterms.data["gr"] = fgcmcal.Sedboundaryterm(primary="g",
                                                             secondary="r")
config.sedterms = fgcmcal.SedtermDict()
config.sedterms.data["g"] = fgcmcal.Sedterm(primaryTerm="gr", secondaryTerm=None,
                                            extrapolated=False, constant=0.0)
config.sedterms.data["r"] = fgcmcal.Sedterm(primaryTerm="gr", secondaryTerm=None,
                                            extrapolated=False, constant=1.0)
config.sedterms.data["i"] = fgcmcal.Sedterm(primaryTerm="gr", secondaryTerm=None,
                                            extrapolated=False, constant=0.75)
# Define good stars with an g-r color cut.
config.starColorCuts = ("g, r, -0.50, 2.25",)
config.refStarColorCuts = ("g, r, -0.50, 2.25",)
config.useExposureReferenceOffset = True
config.precomputeSuperStarInitialCycle = False
config.superStarSubCcdDict = {"g": True,
                              "r": True,
                              "i": True}
config.superStarPlotCcdResiduals = True
# Allow calibration to work with just 1 exposure on a night.
config.minExpPerNight = 1
# Allow calibration to work with very few stars per exposure.
config.minStarPerExp = 5
# Allow calibration to work with small number of stars in processing batches.
config.nStarPerRun = 50
config.nExpPerRun = 2
# Define g-r color as the primary way to split by color.
config.colorSplitBands = ["g", "r"]
config.freezeStdAtmosphere = True
# For tests, do low-order per-ccd polynomial.
config.superStarSubCcdChebyshevOrder = 1
config.ccdGraySubCcdDict = {"g": False,
                            "r": False,
                            "i": False}
config.ccdGrayFocalPlaneDict = {"g": False,
                                "r": False,
                                "i": False}
config.ccdGrayFocalPlaneFitMinCcd = 1
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
config.useRepeatabilityForExpGrayCutsDict = {"g": False,
                                             "r": False,
                                             "i": False}
config.sigFgcmMaxEGrayDict = {"g": 0.1,
                              "r": 0.1,
                              "i": 0.1}
config.approxThroughputDict = {"g": 1.0,
                              "r": 1.0,
                              "i": 1.0}

config.deltaAperFitPerCcdNx = 2
config.deltaAperFitPerCcdNy = 2
config.deltaAperInnerRadiusArcsec = 2.04
config.deltaAperOuterRadiusArcsec = 2.89
config.doComputeDeltaAperPerVisit = True
config.doComputeDeltaAperMap = True
config.doComputeDeltaAperPerCcd = True
