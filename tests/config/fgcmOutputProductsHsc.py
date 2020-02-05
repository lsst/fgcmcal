import os.path

from lsst.utils import getPackageDir

config.cycleNumber = 2

config.doReferenceCalibration = True

config.photoCal.applyColorTerms = True
config.photoCal.photoCatName = 'sdss-dr9-fink-v5b'
config.photoCal.colorterms.load(os.path.join(getPackageDir('fgcmcal'),
                                             'tests',
                                             'config',
                                             'sdssColortermsHsc.py'))
config.refObjLoader.ref_dataset_name = "sdss-dr9-fink-v5b"
