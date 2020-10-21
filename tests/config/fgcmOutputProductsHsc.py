import os.path

from lsst.utils import getPackageDir

config.cycleNumber = 2

config.doReferenceCalibration = True

config.photoCal.applyColorTerms = True

config.photoCal.photoCatName = 'ps1_pv3_3pi_20170110'
config.photoCal.colorterms.load(os.path.join(getPackageDir('obs_subaru'),
                                             'config',
                                             'colorterms.py'))
config.refObjLoader.ref_dataset_name = 'ps1_pv3_3pi_20170110'
config.connections.refCat = 'ps1_pv3_3pi_20170110'
