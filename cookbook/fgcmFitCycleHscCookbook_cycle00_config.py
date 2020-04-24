"""
HSC-specific overrides for FgcmFitCycle
"""
# Output file base name for diagnostic plots
config.outfileBase = 'fgcmFitCycleHscCookbook'
# Bands to be used in the fit
config.bands = ['g', 'r', 'i', 'z', 'y']
config.fitBands = ['g', 'r', 'i', 'z', 'y']
config.requiredBands = ['g', 'r', 'i', 'z', 'y']
# Maximum number of fit iterations (15 for testing, 50+ for a full run.)
config.maxIterBeforeFinalCycle = 30
# Dividing point between observing epochs (years, camera events, etc.)
config.epochMjds = [56700., 57420., 57606.]
# Number of bins to do aperture correction.  Not currently supported in LSST stack.
config.aperCorrFitNBins = 0
