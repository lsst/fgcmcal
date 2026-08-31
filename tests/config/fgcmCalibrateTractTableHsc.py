# All camera defaults were copied from obs_subaru/config/fgcmCalibrateTractTableHsc.py
# on 07/21/21, weekly w_2021_29.

config.fgcmBuildStars.load('fgcmBuildStarsTableHsc.py')
config.fgcmFitCycle.load('fgcmFitCycleHsc.py')

config.maxFitCycles = 3
config.fgcmFitCycle.aperCorrFitNBins = 0
config.fgcmFitCycle.useRepeatabilityForExpGrayCutsDict = {'g': True,
                                                          'r': True,
                                                          'i': True}
config.fgcmFitCycle.maxIterBeforeFinalCycle = 20

config.connections.refCat = 'ps1_pv3_3pi_20170110'
