import os
from lsst.obs.hsc.hscFilters import HSC_FILTER_DEFINITIONS

config.instFluxField = 'apFlux_12_0_instFlux'
config.sourceSelector["science"].signalToNoise.fluxField = 'apFlux_12_0_instFlux'
config.sourceSelector["science"].signalToNoise.errField = 'apFlux_12_0_instFluxErr'
config.sourceSelector["science"].signalToNoise.minimum = 11.0

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
config.connections.ref_cat = "ps1_pv3_3pi_20170110"
# This is set low enough to ensure that this code path is triggered.
config.densityCutMaxPerPixel = 200

configDir = os.path.join(os.path.dirname(__file__))
config.physicalFilterMap = HSC_FILTER_DEFINITIONS.physical_to_band
config.doSubtractLocalBackground = True
config.sourceSelector["science"].flags.bad.append("localBackground_flag")
config.fgcmLoadReferenceCatalog.load(os.path.join(configDir, 'filterMapHsc.py'))
config.fgcmLoadReferenceCatalog.applyColorTerms = True
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(configDir, 'colortermsHsc.py'))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
