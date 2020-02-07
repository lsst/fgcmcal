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
config.doSubtractLocalBackground = False
# The testdata catalogs have the old name
config.psfCandidateName = 'calib_psfCandidate'
config.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name = 'sdss-dr9-fink-v5b'
config.fgcmLoadReferenceCatalog.refFilterMap = {'g': 'g', 'r': 'r', 'i': 'i'}
config.fgcmLoadReferenceCatalog.applyColorTerms = True
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(getPackageDir('fgcmcal'),
                                                                           'tests',
                                                                           'config',
                                                                           'sdssColortermsHsc.py'))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0
