import os.path

from lsst.utils import getPackageDir

config.filterMap = {'g': 'g', 'r': 'r', 'i': 'i'}
config.requiredBands = ['r', 'i']
config.primaryBands = ['i']
config.minPerBand = 2
config.checkAllCcds = False
config.coarseNside = 64
config.visitDataRefName = 'visit'
config.ccdDataRefName = 'ccd'
config.doReferenceMatches = True
# The testdata do not have local background information
config.doSubtractLocalBackground = True
config.nVisitsPerCheckpoint = 5
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0
