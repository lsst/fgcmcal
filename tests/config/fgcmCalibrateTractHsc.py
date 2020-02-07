import os.path

from lsst.utils import getPackageDir

config.fgcmBuildStars.load(os.path.join(getPackageDir('fgcmcal'),
                                        'tests',
                                        'config',
                                        'fgcmBuildStarsHsc.py'))
config.fgcmFitCycle.load(os.path.join(getPackageDir('fgcmcal'),
                                      'tests',
                                      'config',
                                      'fgcmFitCycleHsc.py'))
config.maxFitCycles = 3
config.fgcmOutputProducts.doRefcatOutput = True
