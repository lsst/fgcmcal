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

Run test suite on fgcmcal using Gen3 HSC data from testdata_jointcal.
"""
import matplotlib
matplotlib.use("Agg")  # noqa E402

import unittest
import os
import tempfile
import numpy as np

import lsst.utils
import lsst.pipe.tasks

import fgcmcalTestBase

ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmcalTestHSC(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
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
        testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")
        self.setUp_base(testDir)

        self._importRepository('lsst.obs.subaru.HyperSuprimeCam',
                               os.path.join(self.dataDir, 'hsc'),
                               os.path.join(self.dataDir, 'hsc', 'exports.yaml'))

    def test_fgcmcalPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut('HSC', nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        nStar = 305
        nObs = 3789

        self._testFgcmBuildStarsTable('HSC', "physical_filter IN ('HSC-G', 'HSC-R', 'HSC-I')",
                                      visits, nStar, nObs)

        nZp = 1120
        nGoodZp = 27
        nOkZp = 27
        nBadZp = 1093
        nStdStars = 237
        nPlots = 35

        self._testFgcmFitCycle('HSC', 0, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)
        self._testFgcmFitCycle('HSC', 1, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)

        # We need to create an extra config file to turn on "sub-ccd gray" for testing.
        extraConfigFile = os.path.join(self.testDir, "cycle03_patch_config.py")
        with open(extraConfigFile, "w") as f:
            f.write("config.isFinalCycle = True\n")
            f.write("config.ccdGraySubCcdDict = {'g': True, 'r': True, 'i': True}\n")

        self._testFgcmFitCycle('HSC', 2, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots,
                               extraConfig=extraConfigFile)

        zpOffsets = np.array([0.0010470541892573237, 0.005398149602115154])

        self._testFgcmOutputProducts('HSC', zpOffsets, 36236, 87, 'i', 1)

    def test_fgcmcalTractPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut('HSC', nBand, i0Std, i0Recon, i10Std, i10Recon)

        rawRepeatability = np.array([0.0, 0.008282480993703009, 0.006739350255884648])
        filterNCalibMap = {'HSC-R': 14,
                           'HSC-I': 15}

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]
        tract = 9697

        self._testFgcmCalibrateTract('HSC', visits, tract, 'hsc_rings_v1',
                                     rawRepeatability, filterNCalibMap)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
