"""
HSC-specific overrides for FgcmBuildStars
"""

import os.path

from lsst.utils import getPackageDir

# Check repo for all CCDs for each visit
config.checkAllCcds = True
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
# The reference CCD is a good CCD used to select visit to speed up the scanning
config.referenceCCD = 13
# If smatch matching is available, use this nside.  Not used with default LSST stack.
config.matchNside = 4096
# A star must be observed in one of these bands to be considered as a calibration star
config.primaryBands = ['i']
# Match reference catalog as additional constraint on calibration
config.doReferenceMatches = True
# Subtract the local background before performing calibration?
config.doSubtractLocalBackground = True
# Number of visits read between checkpoints
config.nVisitsPerCheckpoint = 100

# Reference object loader configuration parameters
config.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name = 'ps1_pv3_3pi_20170110'
config.fgcmLoadReferenceCatalog.refFilterMap = {'g': 'g', 'r': 'r', 'i': 'i', 'z': 'z', 'y': 'y'}
config.fgcmLoadReferenceCatalog.applyColorTerms = True
hscConfigDir = os.path.join(getPackageDir('obs_subaru'), 'config')
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(hscConfigDir, 'colorterms.py'))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 10.0
