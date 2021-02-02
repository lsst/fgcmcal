"""
HSC-specific overrides for FgcmFitCycle
"""

from lsst.fgcmcal import Sedterm, Sedboundaryterm

# Output file base name for diagnostic plots
config.outfileBase = 'fgcmFitCycleHscCookbook'
# Bands to run calibration
config.bands = ('g', 'r', 'i', 'z', 'y')
# Bands to use in atmospheric fit.
config.fitBands = ('g', 'r', 'i', 'z', 'y')
# Bands that are required for a star to be considered a calibration star.
config.requiredBands = ('g', 'r', 'i', 'z', 'y')


# Maximum number of fit iterations (15 for testing, 50+ for a full run.)
config.maxIterBeforeFinalCycle = 30
# Number of cores to run with python multiprocessing
config.nCore = 4
# MJD dates on which the mirror was washed
config.washMjds = (56700.0, 57500.0, 57700.0, 58050.0)
# Dividing point between observing epochs (years, camera events, etc.)
config.epochMjds = (56700., 57420., 57606.)
# Number of bins to do aperture correction.  Not currently supported in LSST stack.
config.aperCorrFitNBins = 0

config.sedboundaryterms.data = {'gr': Sedboundaryterm(primary='g', secondary='r'),
                                'ri': Sedboundaryterm(primary='r', secondary='i'),
                                'iz': Sedboundaryterm(primary='i', secondary='z'),
                                'zy': Sedboundaryterm(primary='z', secondary='y')}
config.sedterms.data = {'g': Sedterm(primaryTerm='gr', secondaryTerm='ri', constant=1.6),
                        'r': Sedterm(primaryTerm='gr', secondaryTerm='ri', constant=0.9),
                        'i': Sedterm(primaryTerm='ri', secondaryTerm='iz', constant=1.0),
                        'z': Sedterm(primaryTerm='iz', secondaryTerm='zy', constant=1.0),
                        'y': Sedterm(primaryTerm='zy', secondaryTerm='iz', constant=0.25,
                                     extrapolated=True, primaryBand='y', secondaryBand='z',
                                     tertiaryBand='i')}
