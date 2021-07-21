# All camera defaults were copied from obs_subaru/config/fgcmBuildStars.py
# on 07/21/21, weekly w_2021_29.

import os
from lsst.obs.hsc.hscFilters import HSC_FILTER_DEFINITIONS


config.requiredBands = ['r', 'i']
config.primaryBands = ['i']
config.minPerBand = 2
config.checkAllCcds = False
config.coarseNside = 64
config.visitDataRefName = 'visit'
config.ccdDataRefName = 'ccd'
config.doReferenceMatches = True
config.nVisitsPerCheckpoint = 5
# The test data do not have persisted backgrounds, so don't use them
config.doModelErrorsWithBackground = False
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0
config.densityCutMaxPerPixel = 2000

configDir = os.path.join(os.path.dirname(__file__))
config.physicalFilterMap = HSC_FILTER_DEFINITIONS.physical_to_band
config.doSubtractLocalBackground = True
config.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name = 'ps1_pv3_3pi_20170110'
config.fgcmLoadReferenceCatalog.load(os.path.join(configDir, 'filterMap.py'))
config.fgcmLoadReferenceCatalog.applyColorTerms = True
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(configDir, 'colorterms.py'))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
