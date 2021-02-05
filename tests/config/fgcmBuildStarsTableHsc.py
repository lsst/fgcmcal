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

config.connections.refCat = "ps1_pv3_3pi_20170110"
