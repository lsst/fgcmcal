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

Run test suite on fgcmcal using Gen2 HSC data from testdata_jointcal.
"""

import matplotlib
matplotlib.use("Agg")  # noqa E402

import unittest
import os
import copy
import tempfile
import numpy as np

import lsst.utils
import lsst.pipe.tasks

import fgcmcalTestBaseGen2

import lsst.fgcmcal as fgcmcal

ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmcalTestHSCGen2(fgcmcalTestBaseGen2.FgcmcalTestBaseGen2, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")
        try:
            lsst.utils.getPackageDir('obs_subaru')
        except LookupError:
            raise unittest.SkipTest("obs_subaru not setup")

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
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmMakeLutHsc.py')
        self.configfiles = [testConfigFile]

        self.otherArgs = []

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut(nBand, i0Std, i0Recon, i10Std, i10Recon)

        # Build the stars, adding in the reference stars
        self.config = fgcmcal.FgcmBuildStarsTableConfig()
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmBuildStarsTableHsc.py')
        self.configfiles = [testConfigFile]
        self.otherArgs = []

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        nStar = 305
        nObs = 3789

        self._testFgcmBuildStarsTable(visits, nStar, nObs)

        self.config = fgcmcal.FgcmBuildStarsConfig()
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmBuildStarsHsc.py')
        self.configfiles = [testConfigFile]
        self.otherArgs = []

        self._testFgcmBuildStarsAndCompare(visits)

        # Perform the fit cycle
        self.config = fgcmcal.FgcmFitCycleConfig()
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmFitCycleHsc.py')
        self.config.load(testConfigFile)
        self.configfiles = [testConfigFile]
        self.otherArgs = []

        nZp = 1120
        nGoodZp = 27
        nOkZp = 27
        nBadZp = 1093
        nStdStars = 237
        nPlots = 35

        self._testFgcmFitCycle(nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)

        # Test the second fit cycle -- need to copy to unfreeze config
        newConfig = copy.copy(self.config)
        newConfig.update(cycleNumber=1)
        newConfig.connections.update(cycleNumber='1',
                                     previousCycleNumber='0')
        self.config = newConfig

        newConfigFile = os.path.join(self.testDir,
                                     f'fgcmFitCycle_cycle{newConfig.cycleNumber}.py')
        newConfig.save(newConfigFile)
        self.configfiles.append(newConfigFile)

        self._testFgcmFitCycle(nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)

        # Test the "final" fit cycle
        newConfig = copy.copy(self.config)
        newConfig.update(cycleNumber=2,
                         ccdGraySubCcdDict={'g': True, 'r': True, 'i': True},
                         isFinalCycle=True)
        newConfig.connections.update(cycleNumber='2',
                                     previousCycleNumber='1')
        self.config = newConfig

        newConfigFile = os.path.join(self.testDir,
                                     f'fgcmFitCycle_cycle{newConfig.cycleNumber}.py')
        newConfig.save(newConfigFile)
        self.configfiles.append(newConfigFile)

        self._testFgcmFitCycle(nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots)

        # Output the products

        self.config = fgcmcal.FgcmOutputProductsConfig()
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmOutputProductsHsc.py')
        self.configfiles = [testConfigFile]
        self.otherArgs = []

        filterMapping = {'r': 'HSC-R', 'i': 'HSC-I'}
        zpOffsets = np.array([0.0010470541892573237, 0.005398149602115154])

        self._testFgcmOutputProducts(visitDataRefName, ccdDataRefName,
                                     filterMapping, zpOffsets,
                                     36236, 87, 'i', 1)

    def test_fgcmcalTract(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        # First need to make the LUT
        self.config = fgcmcal.FgcmMakeLutConfig()
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmMakeLutHsc.py')
        self.configfiles = [testConfigFile]
        self.otherArgs = []

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut(nBand, i0Std, i0Recon, i10Std, i10Recon)

        self.config = fgcmcal.FgcmCalibrateTractTableConfig()
        testConfigFile = os.path.join(ROOT, 'config', 'fgcmCalibrateTractTableHsc.py')
        self.configfiles = [testConfigFile]
        self.otherArgs = []

        rawRepeatability = np.array([0.0, 0.008282480993703009, 0.006739350255884648])
        filterNCalibMap = {'HSC-R': 14,
                           'HSC-I': 15}

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]
        tract = 9697

        self._testFgcmCalibrateTract(visits, tract,
                                     rawRepeatability, filterNCalibMap)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
