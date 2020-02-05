import lsst.fgcmcal as fgcmcal

config.outfileBase = 'TestFgcm'
config.filterMap = {'r': 'r', 'i': 'i'}
config.bands = ['r', 'i']
config.fitFlag = (1, 1)
config.requiredFlag = (1, 1)
config.maxIterBeforeFinalCycle = 5
config.nCore = 1
config.washMjds = (0.0, )
config.epochMjds = (0.0, 100000.0)
config.expGrayPhotometricCut = (-0.1, -0.1)
config.expGrayHighCut = (0.1, 0.1)
config.expVarGrayPhotometricCut = 0.05**2.
config.autoPhotometricCutNSig = 5.0
config.autoHighCutNSig = 5.0
config.aperCorrFitNBins = 0
config.aperCorrInputSlopes = (-0.9694, -1.7229)
config.sedboundaryterms = fgcmcal.SedboundarytermDict()
config.sedboundaryterms.data['ri'] = fgcmcal.Sedboundaryterm(primary='r',
                                                             secondary='i')
config.sedterms = fgcmcal.SedtermDict()
config.sedterms.data['r'] = fgcmcal.Sedterm(primaryTerm='ri', secondaryTerm=None,
                                            extrapolated=False, constant=1.0)
config.sedterms.data['i'] = fgcmcal.Sedterm(primaryTerm='ri', secondaryTerm=None,
                                            extrapolated=False, constant=0.75)
config.starColorCuts = ('r,i,-0.50,2.25',)
config.precomputeSuperStarInitialCycle = False
config.minExpPerNight = 1
config.minStarPerExp = 50
config.nStarPerRun = 50
config.nExpPerRun = 2
config.colorSplitIndices = (0, 1)
config.superStarSubCcdChebyshevOrder = 1
config.ccdGraySubCcd = False
config.modelMagErrors = False
config.sigmaCalRange = (0.003, 0.003)
config.instrumentParsPerBand = False
config.useRepeatabilityForExpGrayCuts = (False, False)
config.sigFgcmMaxEGray = (0.05, 0.05)
