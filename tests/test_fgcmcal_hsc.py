# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import inspect
import unittest
import os
import tempfile
import numpy as np

import lsst.utils
import lsst.pex.exceptions

import fgcmcalTestBase

import lsst.fgcmcal as fgcmcal

ROOT = os.path.abspath(os.path.dirname(__file__))

class FgcmcalTestHSC(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        inputDir = os.path.join(self.dataDir, 'hsc')

        testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")

        self.setUp_base(inputDir=inputDir, testDir=testDir)

        lsst.log.setLevel("HscMapper", lsst.log.FATAL)

    def test_fgcmcalTasks(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        # First test making the LUT
        self.config = fgcmcal.FgcmMakeLutConfig()
        self.config.filterNames = ['r', 'i']
        self.config.stdFilterNames = ['r', 'i']
        self.config.atmosphereTableName = 'fgcm_atm_subaru1_test'
        self.otherArgs = []

        nBand = 2
        i0Std = [0.07877351, 0.06464688]
        i10Std = [-0.00061516, -0.00063434]
        i0Recon = [0.06897538, 0.05616964]
        i10Recon = [-6.97094875, 3.69335364]

        self._runFgcmMakeLut(nBand, i0Std, i0Recon, i10Std, i10Recon)

        # Now the star building
        self.config = fgcmcal.FgcmBuildStarsConfig()
        self.config.filterToBand = {'r':'r', 'i':'i'}
        self.config.requiredBands = ['r', 'i']
        self.config.referenceBand = 'i'
        self.config.checkAllCcds = True
        self.otherArgs = []

        nVisit = 11
        nStar = 472
        nObs = 5431

        self._runFgcmBuildStars(nVisit, nStar, nObs)

        # And the fit cycle
        self.config = fgcmcal.FgcmFitCycleConfig()
        self.config.outfileBase = 'TestFgcm'
        self.config.bands = ['r', 'i']
        self.config.fitFlag = (1, 1)
        self.config.filterToBand = {'r':'r', 'i':'i'}
        self.config.maxIter = 1
        self.config.nCore = 1
        self.config.cycleNumber = 0
        self.config.utBoundary = 0.0
        self.config.washMjds = (0.0, )
        self.config.epochMjds = (0.0, 100000.0)
        self.config.latitude = 19.8256
        self.config.cameraGain = 3.0
        self.config.pixelScale = 0.17
        self.config.expGrayPhotometricCut = (-0.05, -0.05)
        self.config.expGrayHighCut = (0.2, 0.2)
        self.config.aperCorrFitNBins = 0
        self.config.sedFudgeFactors = (1.0, 1.0)
        self.config.starColorCuts = ('r,i,-0.50,2.25',)
        self.config.freezeStdAtmosphere = True
        self.config.precomputeSuperStarInitialCycle = False
        self.config.minExpPerNight = 1
        self.config.minStarPerExp = 50
        self.config.nStarPerRun = 50
        self.config.nExpPerRun = 2
        self.config.colorSplitIndices = (0, 1)
        self.otherArgs = []

        nZp = 1232
        nGoodZp = 27

        self._runFgcmFitCycle(nZp, nGoodZp)

        # And output the products
        self.config = fgcmcal.FgcmOutputProductsConfig()
        self.config.cycleNumber = 0

        self._runFgcmOutputProducts()



class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
