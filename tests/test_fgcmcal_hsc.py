# See COPYRIGHT file at the top of the source tree.
#
# This file is part of fgcmcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Test the fgcmcal code with testdata_jointcal/hsc.

Run test suite on fgcmcal using HSC data from testdata_jointcal.
"""

import matplotlib
matplotlib.use("Agg")  # noqa E402

import unittest
import os
import tempfile
import numpy as np

import lsst.utils
import lsst.pipe.tasks
from lsst.pipe.tasks.colorterms import ColortermLibrary

import fgcmcalTestBase

import lsst.fgcmcal as fgcmcal

ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmcalTestHSC(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.dataDir, 'hsc_and_index')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        inputDir = os.path.join(self.dataDir, 'hsc')

        self.testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")

        self.setUp_base(inputDir=inputDir, testDir=self.testDir)

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
        self.config.atmosphereTableName = 'fgcm_atm_subaru2_test'
        self.otherArgs = []

        nBand = 2
        i0Std = np.array([0.07877351, 0.06464688])
        i10Std = np.array([-0.00061516, -0.00063434])
        i0Recon = np.array([0.0689530429, 0.05600673])
        i10Recon = np.array([-7.01847144, 3.62675740])

        self._testFgcmMakeLut(nBand, i0Std, i0Recon, i10Std, i10Recon)

        # Build the stars, adding in the reference stars
        self.config = fgcmcal.FgcmBuildStarsConfig()
        self.fillDefaultBuildStarsConfig(self.config, visitDataRefName, ccdDataRefName)
        self.otherArgs = []

        nVisit = 11
        nStar = 472
        nObs = 5431

        self._testFgcmBuildStars(nVisit, nStar, nObs)

        # Perform the fit cycle
        self.config = fgcmcal.FgcmFitCycleConfig()
        self.fillDefaultFitCycleConfig(self.config)
        self.otherArgs = []

        nZp = 1232
        nGoodZp = 26
        nOkZp = 26
        nBadZp = 1206
        nStdStars = 390
        nPlots = 34

        self._testFgcmFitCycle(nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)

        # Test the second fit cycle -- need to copy to unfreeze config
        newConfig = fgcmcal.FgcmFitCycleConfig()
        newConfig.update(**self.config.toDict())
        newConfig.cycleNumber = 1
        self.config = newConfig

        self._testFgcmFitCycle(nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)

        # Test the "final" fit cycle
        newConfig = fgcmcal.FgcmFitCycleConfig()
        newConfig.update(**self.config.toDict())
        newConfig.cycleNumber = 2
        newConfig.ccdGraySubCcd = True
        newConfig.isFinalCycle = True
        self.config = newConfig

        self._testFgcmFitCycle(nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots)

        # Output the products

        self.config = fgcmcal.FgcmOutputProductsConfig()
        self.config.cycleNumber = 2
        # Turning on "reference calibration" is redundant in practice if
        # reference stars have been used in the fit, but this needs to
        # be exercised in testing.
        self.config.doReferenceCalibration = True
        self.config.photoCal.applyColorTerms = True
        self.config.photoCal.photoCatName = 'sdss-dr9-fink-v5b'
        self.config.photoCal.colorterms = self.sdssColorterms()
        self.config.refObjLoader.ref_dataset_name = "sdss-dr9-fink-v5b"

        filterMapping = {'r': 'HSC-R', 'i': 'HSC-I'}
        zpOffsets = np.array([-0.0013903317740, -0.0020539460238])

        self._testFgcmOutputProducts(visitDataRefName, ccdDataRefName,
                                     filterMapping, zpOffsets,
                                     904014, 12, 'i', 1)

    def test_fgcmcalTract(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        visitDataRefName = 'visit'
        ccdDataRefName = 'ccd'

        # First need to make the LUT
        self.config = fgcmcal.FgcmMakeLutConfig()
        self.config.filterNames = ['r', 'i']
        self.config.stdFilterNames = ['r', 'i']
        self.config.atmosphereTableName = 'fgcm_atm_subaru2_test'
        self.otherArgs = []

        nBand = 2
        i0Std = np.array([0.07877351, 0.06464688])
        i10Std = np.array([-0.00061516, -0.00063434])
        i0Recon = np.array([0.0689530429, 0.05600673])
        i10Recon = np.array([-7.01847144, 3.62675740])

        self._testFgcmMakeLut(nBand, i0Std, i0Recon, i10Std, i10Recon)

        self.config = fgcmcal.FgcmCalibrateTractConfig()
        self.fillDefaultBuildStarsConfig(self.config.fgcmBuildStars, visitDataRefName, ccdDataRefName)
        self.config.fgcmBuildStars.checkAllCcds = False
        self.fillDefaultFitCycleConfig(self.config.fgcmFitCycle)
        self.config.maxFitCycles = 2

        self.config.fgcmOutputProducts.doRefcatOutput = True

        rawRepeatability = np.array([0.007070288705, 0.0074971053995])
        filterNCalibMap = {'HSC-R': 13,
                           'HSC-I': 13}

        visits = [903334, 903336, 903338, 903342, 903344, 903346,
                  903986, 903988, 903990, 904010, 904014]
        tract = 0

        self._testFgcmCalibrateTract(visits, tract,
                                     rawRepeatability, filterNCalibMap)

    def sdssColorterms(self):
        """
        Return the SDSS to HSC color terms
        """

        colorterms = ColortermLibrary()
        colorterms.data = {}
        colorterms.data['sdss*'] = lsst.pipe.tasks.colorterms.ColortermDict()
        colorterms.data['sdss*'].data = {}
        colorterms.data['sdss*'].data['g'] = lsst.pipe.tasks.colorterms.Colorterm()
        colorterms.data['sdss*'].data['g'].primary = 'g'
        colorterms.data['sdss*'].data['g'].secondary = 'r'
        colorterms.data['sdss*'].data['g'].c0 = -0.00816446
        colorterms.data['sdss*'].data['g'].c1 = -0.08366937
        colorterms.data['sdss*'].data['g'].c2 = -0.00726883
        colorterms.data['sdss*'].data['r'] = lsst.pipe.tasks.colorterms.Colorterm()
        colorterms.data['sdss*'].data['r'].primary = 'r'
        colorterms.data['sdss*'].data['r'].secondary = 'i'
        colorterms.data['sdss*'].data['r'].c0 = 0.0013181
        colorterms.data['sdss*'].data['r'].c1 = 0.01284177
        colorterms.data['sdss*'].data['r'].c2 = -0.03068248
        colorterms.data['sdss*'].data['i'] = lsst.pipe.tasks.colorterms.Colorterm()
        colorterms.data['sdss*'].data['i'].primary = 'i'
        colorterms.data['sdss*'].data['i'].secondary = 'z'
        colorterms.data['sdss*'].data['i'].c0 = 0.00130204
        colorterms.data['sdss*'].data['i'].c1 = -0.16922042
        colorterms.data['sdss*'].data['i'].c2 = -0.01374245

        return colorterms

    def fillDefaultBuildStarsConfig(self, config, visitDataRefName, ccdDataRefName):
        """
        Fill the config parameters for a build stars configuration

        Parameters
        ----------
        config: `lsst.fgcmcal.FgcmBuildStarsConfig`
        visitDataRefName: `str`
           Name of the dataRef key for the visit
        ccdDataRefName: `str`
           Name of the dataRef key for the ccd
        """

        config.filterMap = {'r': 'r', 'i': 'i'}
        config.requiredBands = ['r', 'i']
        config.primaryBands = ['i']
        config.checkAllCcds = True
        config.coarseNside = 64
        config.visitDataRefName = visitDataRefName
        config.ccdDataRefName = ccdDataRefName
        config.doReferenceMatches = True
        # The testdata catalogs have the old name
        config.psfCandidateName = 'calib_psfCandidate'
        config.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name = 'sdss-dr9-fink-v5b'
        config.fgcmLoadReferenceCatalog.applyColorTerms = True
        config.fgcmLoadReferenceCatalog.colorterms = self.sdssColorterms()
        config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
        config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = 'i_flux'
        config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = 'i_fluxErr'
        config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0

    def fillDefaultFitCycleConfig(self, config):
        """
        Fill the config parameters for a fit cycle configuration.

        Parameters
        ----------
        config: `lsst.fgcmcal.FgcmFitCycleConfig`
        """

        config.outfileBase = 'TestFgcm'
        config.bands = ['r', 'i']
        config.fitFlag = (1, 1)
        config.requiredFlag = (1, 1)
        config.filterMap = {'r': 'r', 'i': 'i'}
        config.doReferenceCalibration = True
        config.maxIterBeforeFinalCycle = 5
        config.nCore = 1
        config.cycleNumber = 0
        config.utBoundary = 0.0
        config.washMjds = (0.0, )
        config.epochMjds = (0.0, 100000.0)
        config.latitude = 19.8256
        config.expGrayPhotometricCut = (-0.05, -0.05)
        config.expGrayHighCut = (0.2, 0.2)
        config.aperCorrFitNBins = 0
        config.aperCorrInputSlopes = (-0.9694, -1.7229)
        config.sedFudgeFactors = (1.0, 1.0)
        config.starColorCuts = ('r,i,-0.50,2.25',)
        config.freezeStdAtmosphere = True
        config.precomputeSuperStarInitialCycle = False
        config.minExpPerNight = 1
        config.minStarPerExp = 50
        config.nStarPerRun = 50
        config.nExpPerRun = 2
        config.colorSplitIndices = (0, 1)
        config.superStarSubCcd = True
        config.superStarSubCcdChebyshevOrder = 1
        config.ccdGraySubCcd = False
        config.ccdGraySubCcdChebyshevOrder = 1
        config.ccdGraySubCcdTriangular = True
        config.modelMagErrors = False
        config.outputZeropointsBeforeFinalCycle = False
        config.outputStandardsBeforeFinalCycle = False
        config.sigmaCalRange = (0.003, 0.003)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
