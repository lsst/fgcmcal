import os.path


config.connections.cycleNumber = 2

config.doReferenceCalibration = True
config.referencePixelizationMinStars = 20
config.referenceMinMatch = 15

from lsst.obs.lsst.filters import LATISS_FILTER_DEFINITIONS

config.physicalFilterMap = {
    'SDSSg_65mm~empty': 'g',
    'SDSSr_65mm~empty': 'r',
    'SDSSi_65mm~empty': 'i',
}

config.photoCal.applyColorTerms = False
config.photoCal.photoCatName = 'atlas_refcat2_20220201'

config.connections.refCat = 'atlas_refcat2_20220201'
