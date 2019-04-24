"""
HSC-specific overrides for FgcmBuildStars
"""

import os.path

from lsst.utils import getPackageDir

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
# Which bands are required to be observed to be considered a calibration star
config.requiredBands = ['g', 'r', 'i', 'z']
# The reference band is used for initial star selection
config.primaryBands = ['i']
# The reference CCD is a good CCD used to select visit to speed up the scanning
config.referenceCCD = 13

hscConfigDir = os.path.join(getPackageDir("obs_subaru"), "config", "hsc")

# The filter map which goes from filter name to (abstract) reference band
config.load(os.path.join(hscConfigDir, 'filterMap.py'))

config.referenceBands = ['g', 'r', 'i', 'z', 'y']
config.doReferenceMatches = True
config.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name = 'ps1_pv3_3pi_20170110'
config.fgcmLoadReferenceCatalog.refObjLoader.load(os.path.join(hscConfigDir, 'filterMap.py'))
config.fgcmLoadReferenceCatalog.applyColorTerms = True
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(hscConfigDir, 'colorterms.py'))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 10.0
