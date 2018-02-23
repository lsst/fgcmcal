"""
HSC-specific overrides for FgcmFitCycle
"""

# Output file base name for diagnostic plots
config.outfileBase = 'fgcmFitCycleHscCookbook'
# Bands to be used in the fit
config.bands = ('g', 'r', 'i', 'z', 'y')
# Flag to specify if these should be part of the fit, or generated from the atmosphere
# parameters.
config.fitFlag = (1, 1, 1, 1, 1)
# Dictionary that maps "filters" (instrumental configurations) to "bands"
# (abstract names).  All filters must be listed in the LUT.
config.filterToBand = {'g':'g', 'r':'r', 'i':'i', 'i2':'i', 'z':'z', 'y':'y'}
# Maximum number of fit iterations (15 for testing, 50+ for a full run.)
config.maxIter = 15
# Number of cores to run with python multiprocessing
config.nCore = 4
# Cycle number (should start at 0)
config.cycleNumber = 0
# Value to add to MJD to ensure that different MJDs fall on different nights
# This value will depend on your longitude/time zone!
config.utBoundary = 0.0
# MJD dates on which the mirror was washed
config.washMjds = (0.0, )
# Dividing point between observing epochs (years, camera events, etc.)
config.epochMjds = (0.0, 100000.0)
# Latitude of the observatory
config.latitude = 19.8256
# Gain of the camera (used to dead-reckon absolute scale)
config.cameraGain = 3.0
# Pixel scale (arcseconds)
config.pixelScale = 0.17
# Amount of gray extinction to be considered "photometric".  This will
# get ratcheded down in further cycles
config.expGrayPhotometricCut = (-0.05, -0.05, -0.05, -0.05, -0.05)
# Number of bins to do aperture correction.  Not currently supported in LSST stack.
config.aperCorrFitNBins = 0
# "Fudge factors" for computing SED slope (best values for HSC not determined yet)
config.sedFudgeFactors = (1.0, 1.0, 1.0, 1.0, 1.0)
# Color cuts for stars to use for calibration.  Each element is a string with
# band1, band2, range_low, range_high such that range_low < (band1 - band2) < range_high
config.starColorCuts = ('g,r,-0.25,2.25',
                        'r,i,-0.50,2.25',
                        'i,z,-0.50,1.00',
                        'g,i,0.0,3.5')
# Freeze atmosphere to standard values?  Recommended for first fit cycle.
config.freezeStdAtmosphere = True
# Precompute "superstar" in initial cycle (==00) based on bright star observations?  Recommended for HSC.
config.precomputeSuperStarInitialCycle = True
