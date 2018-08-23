"""
HSC-specific overrides for FgcmBuildStars
"""

# Minimum number of observations per band for a star to be considered for calibration
config.minPerBand = 2
# Match radius to associate stars from src catalogs (arcseconds)
config.matchRadius = 1.0
# Isolation radius: stars must be at least this far from a neighbor to be considered (arcseconds)
config.isolationRadius = 2.0
# Measure the stellar density with healpix nside=densityCutNside
config.densityCutNside = 128
# If there are more than densityCutMaxPerPixel stars per pixel, sample them
config.densityCutMaxPerPixel = 1500
# Dictionary that maps "filters" (instrumental configurations) to "bands"
# (abstract names).  All filters must be listed in the LUT.
config.filterToBand = {'g':'g', 'r':'r', 'i':'i', 'z':'z', 'y':'y'}
# Which bands are required to be observed to be considered a calibration star
config.requiredBands = ['g','r','i','z']
# The reference band is used for initial star selection
config.referenceBands = ['i']
# The reference CCD is a good CCD used to select visit to speed up the scanning
config.referenceCCD = 13
# If smatch matching is available, use this nside.  Not used with default LSST stack.
config.matchNside = 4096
