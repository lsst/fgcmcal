# All camera defaults were copied from obs_subaru/config/fgcmOutputProducts.py
# on 07/21/21, weekly w_2021_29.

config.connections.cycleNumber = 2

config.doReferenceCalibration = True

from lsst.obs.hsc.hscFilters import HSC_FILTER_DEFINITIONS
config.physicalFilterMap = HSC_FILTER_DEFINITIONS.physical_to_band

config.photoCal.applyColorTerms = True
config.photoCal.photoCatName = 'ps1_pv3_3pi_20170110'

config.photoCal.colorterms.load('colortermsHsc.py')
config.connections.refCat = 'ps1_pv3_3pi_20170110'
