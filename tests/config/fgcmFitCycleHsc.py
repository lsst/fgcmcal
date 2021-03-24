import lsst.fgcmcal as fgcmcal

config.outfileBase = 'TestFgcm'
# The unused z-band is here to test the case that extra filters
# are in the map that are not used in the calibrations.
config.physicalFilterMap = {'HSC-G': 'g',
                            'HSC-R': 'r',
                            'HSC-I': 'i',
                            'HSC-Z': 'z'}
config.bands = ['g', 'r', 'i']
config.fitBands = ['g', 'r', 'i']
config.requiredBands = ['r', 'i']
config.maxIterBeforeFinalCycle = 5
config.nCore = 1
config.washMjds = (0.0, )
config.epochMjds = (0.0, 100000.0)
config.expGrayPhotometricCutDict = {'g': -0.1, 'r': -0.1, 'i': -0.1}
config.expGrayHighCutDict = {'g': 0.1, 'r': 0.1, 'i': 0.1}
config.expVarGrayPhotometricCutDict = {'g': 0.05**2.,
                                       'r': 0.05**2.,
                                       'i': 0.05**2.}
config.autoPhotometricCutNSig = 5.0
config.autoHighCutNSig = 5.0
config.aperCorrFitNBins = 2
config.aperCorrInputSlopeDict = {'g': -1.0,
                                 'r': -0.9694,
                                 'i': -1.7229}
config.sedboundaryterms = fgcmcal.SedboundarytermDict()
config.sedboundaryterms.data['ri'] = fgcmcal.Sedboundaryterm(primary='r',
                                                             secondary='i')
config.sedterms = fgcmcal.SedtermDict()
config.sedterms.data['g'] = fgcmcal.Sedterm(primaryTerm='ri', secondaryTerm=None,
                                            extrapolated=False, constant=0.0)
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
config.colorSplitBands = ['r', 'i']
config.superStarSubCcdChebyshevOrder = 1
config.ccdGraySubCcdDict = {'g': False,
                            'r': False,
                            'i': False}
config.ccdGrayFocalPlaneDict = {'g': True,
                                'r': True,
                                'i': True}
config.ccdGrayFocalPlaneChebyshevOrder = 1
config.modelMagErrors = False
config.sigmaCalRange = (0.003, 0.003)
config.instrumentParsPerBand = False
config.randomSeed = 12345
