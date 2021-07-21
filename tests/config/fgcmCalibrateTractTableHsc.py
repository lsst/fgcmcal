# All camera defaults were copied from obs_subaru/config/fgcmCalibrateTractTableHsc.py
# on 07/21/21, weekly w_2021_29.

import os.path


configDir = os.path.join(os.path.dirname(__file__))

config.fgcmBuildStars.load(os.path.join(configDir, 'fgcmBuildStarsTableHsc.py'))
config.fgcmFitCycle.load(os.path.join(configDir, 'fgcmFitCycleHsc.py'))

config.maxFitCycles = 3
config.fgcmOutputProducts.doRefcatOutput = True
config.fgcmFitCycle.aperCorrFitNBins = 0
config.fgcmFitCycle.useRepeatabilityForExpGrayCutsDict = {'g': True,
                                                          'r': True,
                                                          'i': True}

config.connections.refCat = 'ps1_pv3_3pi_20170110'
