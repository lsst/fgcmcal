"""
HSC cookbook specific overrides for fgcmOutputProducts
"""

import os.path

from lsst.utils import getPackageDir

# Do not do the post-fit reference catalog calibration
config.doReferenceCalibration = False
# Do the standard star catalog output
config.doRefcatOutput = True
# Do the atmosphere output
config.doAtmosphereOutput = True
# Do output the zeropoints in jointcal_photoCalib format (though uses a lot of space)
config.doZeropointOutput = True

# Last cycle number that was run, preferably with outputStandards == True
config.cycleNumber = 3
# For quicker runs, we cut down the area used for cross-calibration
config.referencePixelizationNPixels = 10

# Reference object info
config.refObjLoader.ref_dataset_name = 'ps1_pv3_3pi_20170110'

# Photometric calibration information
config.photoCal.photoCatName = 'ps1'
hscConfigDir = os.path.join(getPackageDir("obs_subaru"), "config", "hsc")
config.photoCal.colorterms.load(os.path.join(hscConfigDir, 'colorterms.py'))



