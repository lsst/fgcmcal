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
