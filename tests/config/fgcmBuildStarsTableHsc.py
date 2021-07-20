# All camera defaults were copied from obs_subaru/config/fgcmBuildStarsTable.py
# on 07/21/21, weekly w_2021_29.

import os
from lsst.obs.hsc.hscFilters import HSC_FILTER_DEFINITIONS

# This override is to be consistent with the old tests
config.referenceCCD = 13
# The filterMap and bands are for the small subset of bands used in the tests
config.physicalFilterMap = {'HSC-G': 'g', 'HSC-R': 'r', 'HSC-I': 'i'}
config.requiredBands = ['r', 'i']
config.primaryBands = ['i']
# The coarseNside is set appropriate to the area of the test data
config.coarseNside = 64
# We have only a few visits, so checkpointing is more frequent
config.nVisitsPerCheckpoint = 5
# The tests are done with only the brightest reference stars
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0
# The test data do not have persisted backgrounds, so don't use them
config.doModelErrorsWithBackground = False

config.minPerBand = 2
config.connections.refCat = "ps1_pv3_3pi_20170110"
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
