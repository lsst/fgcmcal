"""
HSC-specific overrides for FgcmFitCycle
"""

from lsst.fgcmcal import Sedterm, Sedboundaryterm

# Output file base name for diagnostic plots
config.outfileBase = 'fgcmFitCycleHscCookbook'
# Bands to be used in the fit
config.bands = ('g', 'r', 'i', 'z', 'y')
# Flag to specify if these should be part of the fit, or generated from the atmosphere
# parameters.
config.fitFlag = (1, 1, 1, 1, 1)
# Flag for bands that are required to be a calibration star.
config.requiredFlag = (1, 1, 1, 1, 1)
# Dictionary that maps "filters" (instrumental configurations) to "bands"
# (abstract names).  All filters must be listed in the LUT.
config.filterMap = {'g':'g', 'r':'r', 'i':'i', 'z':'z', 'y':'y'}
# Maximum number of fit iterations (15 for testing, 50+ for a full run.)
config.maxIterBeforeFinalCycle = 30
# Number of cores to run with python multiprocessing
config.nCore = 4
# Cycle number (should start at 0)
config.cycleNumber = 0
# Value to add to MJD to ensure that different MJDs fall on different nights
# This value will depend on your longitude/time zone!
config.utBoundary = 0.0
# MJD dates on which the mirror was washed
config.washMjds = (56700.0, 57500.0, 57700.0, 58050.0)
# Dividing point between observing epochs (years, camera events, etc.)
config.epochMjds = (56700., 57420., 57606.)
# Latitude of the observatory
config.latitude = 19.8256
# Amount of gray extinction to be considered "photometric".  This will
# get ratcheded down in further cycles.
config.expGrayPhotometricCut = (-0.05, -0.05, -0.05, -0.05, -0.05)
# Amount of gray "positive extinction" to be considered "photometric".  Used
# to cull out-of-model exposures.
config.expGrayHighCut = (0.2, 0.2, 0.2, 0.2, 0.2)
# Number of bins to do aperture correction.  Not currently supported in LSST stack.
config.aperCorrFitNBins = 0
# Aperture correction input slope parameters.  There should be one slope ber band.
# This is used when there is insufficient data to fit the parameters from the data
# itself (e.g. tract mode or RC2).
config.aperCorrInputSlopes = (-1.0150, -0.9694, -1.7229, -1.4549, -1.1998)
# Mapping from bands to SED boundary term names used is sedterms.
config.sedboundaryterms.data = {'gr': Sedboundaryterm(primary='g', secondary='r'),
                                'ri': Sedboundaryterm(primary='r', secondary='i'),
                                'iz': Sedboundaryterm(primary='i', secondary='z'),
                                'zy': Sedboundaryterm(primary='z', secondary='y')}
# Mapping from terms to bands for fgcm linear SED approximations.
config.sedterms.data = {'g': Sedterm(primaryTerm='gr', secondaryTerm='ri', constant=1.6),
                        'r': Sedterm(primaryTerm='gr', secondaryTerm='ri', constant=0.9),
                        'i': Sedterm(primaryTerm='ri', secondaryTerm='iz', constant=1.0),
                        'z': Sedterm(primaryTerm='iz', secondaryTerm='zy', constant=1.0),
                        'y': Sedterm(primaryTerm='zy', secondaryTerm='iz', constant=0.25,
                                     extrapolated=True, primaryBand='y', secondaryBand='z',
                                     tertiaryBand='i')}
# Color cuts for stars to use for calibration.  Each element is a string with
# band1, band2, range_low, range_high such that range_low < (band1 - band2) < range_high
config.starColorCuts = ('g,r,-0.25,2.25',
                        'r,i,-0.50,2.25',
                        'i,z,-0.50,1.00',
                        'g,i,0.0,3.5')
# Which band indices are used to do color splits? (Default g-i)
config.colorSplitIndices = (0, 2)
# Freeze atmosphere to standard values?  Recommended for first fit cycle.
config.freezeStdAtmosphere = True
# Precompute "superstar" in initial cycle (==00) based on bright star observations?  Recommended for HSC.
config.precomputeSuperStarInitialCycle = True
# Do sub-ccd calibration of the superstar flat (recommended for HSC)
config.superStarSubCcd = True
# Chebyshev order of sub-ccd superstar fits
config.superStarSubCcdChebyshevOrder = 2
# Compute CCD gray terms on sub-ccd scale
config.ccdGraySubCcd = True
# Order of the 2D chebyshev polynomials for sub-ccd gray fit
config.ccdGraySubCcdChebyshevOrder = 1
# Model instrumental variation over time per band
config.instrumentParsPerBand = True
# Use reference catalog as additional constraint on calibration
config.doReferenceCalibration = True
# Reference star signal-to-noise minimum to use in calibration
config.refStarSnMin = 50.0
# Number of sigma compared to average mag for reference star to be considered an outlier
config.refStarOutlierNSig = 4.0

config.useRepeatabilityForExpGrayCuts = (False, False, False, False, False)
config.sigFgcmMaxEGray = (0.05, 0.05, 0.05, 0.05, 0.05)
