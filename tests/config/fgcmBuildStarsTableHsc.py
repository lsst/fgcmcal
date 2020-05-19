import os.path

import os
from lsst.utils import getPackageDir

config.filterMap = {'g': 'g', 'r': 'r', 'i': 'i'}
config.requiredBands = ['r', 'i']
config.primaryBands = ['i']
config.minPerBand = 2
config.coarseNside = 64
config.visitDataRefName = 'visit'
config.ccdDataRefName = 'ccd'
config.doReferenceMatches = True
# The testdata do not have local background information
config.doSubtractLocalBackground = True
config.nVisitsPerCheckpoint = 5
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0
config.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name = 'ps1_pv3_3pi_20170110'
config.fgcmLoadReferenceCatalog.refFilterMap = {'g': 'g', 'r': 'r', 'r2': 'r',
                                                'i': 'i', 'i2': 'i', 'z': 'z', 'y': 'y',
                                                'N387': 'g', 'N816': 'i', 'N921': 'z',
                                                'N1010': 'y'}
config.fgcmLoadReferenceCatalog.applyColorTerms = True
hscConfigDir = os.path.join(getPackageDir('obs_subaru'), 'config', 'hsc')
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(hscConfigDir, 'colorterms.py'))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 10.0
