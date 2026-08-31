config.connections.cycleNumber = 2

config.doReferenceCalibration = True
config.referencePixelizationMinStars = 20
config.referenceMinMatch = 15
config.doTractStars = False

from lsst.obs.lsst.filters import LATISS_FILTER_DEFINITIONS

config.physicalFilterMap = {
    'SDSSg_65mm~empty': 'g',
    'SDSSr_65mm~empty': 'r',
    'SDSSi_65mm~empty': 'i',
}

config.photoCal.applyColorTerms = True
config.photoCal.photoCatName = 'atlas_refcat2_20220201'

config.photoCal.colorterms.load('colortermsLatiss.py')
config.connections.refCat = 'atlas_refcat2_20220201'
