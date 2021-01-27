import os.path

from lsst.utils import getPackageDir

config.fgcmBuildStars.load(os.path.join(getPackageDir('fgcmcal'),
                                        'tests',
                                        'config',
                                        'fgcmBuildStarsTableHsc.py'))
config.fgcmFitCycle.load(os.path.join(getPackageDir('fgcmcal'),
                                      'tests',
                                      'config',
                                      'fgcmFitCycleHsc.py'))
config.maxFitCycles = 3
config.fgcmOutputProducts.doRefcatOutput = True

config.connections.refCat = 'ps1_pv3_3pi_20170110'
