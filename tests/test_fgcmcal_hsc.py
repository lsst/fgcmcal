# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os
import tempfile
import numpy as np

import lsst.utils
import lsst.pex.exceptions

from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask

import fgcmcalTestBase

import lsst.fgcmcal as fgcmcal

ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmcalTestHSC(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.dataDir, 'hsc_and_index')
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

        visitDataRefName = 'visit'
        ccdDataRefName = 'ccd'

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
        self.config.filterToBand = {'r': 'r', 'i': 'i'}
        self.config.requiredBands = ['r', 'i']
        self.config.referenceBands = ['i']
        self.config.checkAllCcds = True
        self.config.visitDataRefName = visitDataRefName
        self.config.ccdDataRefName = ccdDataRefName
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
        self.config.requiredFlag = (1, 1)
        self.config.filterToBand = {'r': 'r', 'i': 'i'}
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
        self.config.outputStandards = True
        self.config.superStarSubCcd = True
        self.config.superStarSubCcdChebyshevOrder = 1
        self.config.modelMagErrors = False
        self.otherArgs = []

        nZp = 1232
        nGoodZp = 27
        nStdStars = 472

        self._runFgcmFitCycle(nZp, nGoodZp, nStdStars)

        # And output the products
        self.config = fgcmcal.FgcmOutputProductsConfig()
        self.config.cycleNumber = 0
        self.config.photoCal.photoCatName = 'sdss-dr9-fink-v5b'
        self.config.photoCal.colorterms.data = {}
        self.config.photoCal.colorterms.data['sdss*'] = lsst.pipe.tasks.colorterms.ColortermDict()
        self.config.photoCal.colorterms.data['sdss*'].data = {}
        self.config.photoCal.colorterms.data['sdss*'].data['g'] = lsst.pipe.tasks.colorterms.Colorterm()
        self.config.photoCal.colorterms.data['sdss*'].data['g'].primary = 'g'
        self.config.photoCal.colorterms.data['sdss*'].data['g'].secondary = 'r'
        self.config.photoCal.colorterms.data['sdss*'].data['g'].c0 = -0.00816446
        self.config.photoCal.colorterms.data['sdss*'].data['g'].c1 = -0.08366937
        self.config.photoCal.colorterms.data['sdss*'].data['g'].c2 = -0.00726883
        self.config.photoCal.colorterms.data['sdss*'].data['r'] = lsst.pipe.tasks.colorterms.Colorterm()
        self.config.photoCal.colorterms.data['sdss*'].data['r'].primary = 'r'
        self.config.photoCal.colorterms.data['sdss*'].data['r'].secondary = 'i'
        self.config.photoCal.colorterms.data['sdss*'].data['r'].c0 = 0.0013181
        self.config.photoCal.colorterms.data['sdss*'].data['r'].c1 = 0.01284177
        self.config.photoCal.colorterms.data['sdss*'].data['r'].c2 = -0.03068248
        self.config.photoCal.colorterms.data['sdss*'].data['i'] = lsst.pipe.tasks.colorterms.Colorterm()
        self.config.photoCal.colorterms.data['sdss*'].data['i'].primary = 'i'
        self.config.photoCal.colorterms.data['sdss*'].data['i'].secondary = 'z'
        self.config.photoCal.colorterms.data['sdss*'].data['i'].c0 = 0.00130204
        self.config.photoCal.colorterms.data['sdss*'].data['i'].c1 = -0.16922042
        self.config.photoCal.colorterms.data['sdss*'].data['i'].c2 = -0.01374245
        self.config.refObjLoader.retarget(target=LoadAstrometryNetObjectsTask)

        filterMapping = {'r': 'HSC-R', 'i': 'HSC-I'}
        zpOffsets = np.array([8.685752, 8.971653])

        self._runFgcmOutputProducts(visitDataRefName, ccdDataRefName,
                                    filterMapping, zpOffsets,
                                    904014, 12, 'i', 1)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
