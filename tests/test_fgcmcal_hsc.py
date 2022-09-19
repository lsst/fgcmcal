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
import unittest
import os
import tempfile
import numpy as np

import matplotlib
matplotlib.use("Agg")

import lsst.utils  # noqa: E402
import lsst.pipe.tasks  # noqa: E402
import lsst.daf.butler  # noqa: E402

import fgcmcalTestBase  # noqa: E402


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

        lsst.daf.butler.cli.cliLog.CliLog.initLog(longlog=False)

        cls.testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")

        cls._importRepository('lsst.obs.subaru.HyperSuprimeCam',
                              os.path.join(cls.dataDir, 'hsc/repo'),
                              os.path.join(cls.dataDir, 'hsc', 'exports.yaml'))

    def test_fgcmcalPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcalpipe'

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        nStar = 305
        nObs = 3789

        self._testFgcmBuildStarsTable(instName, testName,
                                      "physical_filter IN ('HSC-G', 'HSC-R', 'HSC-I')",
                                      visits, nStar, nObs)

        nZp = 1120
        nGoodZp = 27
        nOkZp = 27
        nBadZp = 1093
        nStdStars = 235
        nPlots = 47

        self._testFgcmFitCycle(instName, testName,
                               0, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)
        self._testFgcmFitCycle(instName, testName,
                               1, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True)

        # We need to create an extra config file to turn on "sub-ccd gray" for testing.
        # We also want to exercise the code path setting useExposureReferenceOffset = False.
        extraConfigFile = os.path.join(self.testDir, "cycle03_patch_config.py")
        with open(extraConfigFile, "w") as f:
            f.write("config.isFinalCycle = True\n")
            f.write("config.ccdGraySubCcdDict = {'g': True, 'r': True, 'i': True}\n")
            f.write("config.useExposureReferenceOffset = False")

        self._testFgcmFitCycle(instName, testName,
                               2, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots,
                               extraConfig=extraConfigFile)

        zpOffsets = np.array([-0.0008051003096625209,
                              0.0072303167544305325])

        self._testFgcmOutputProducts(instName, testName,
                                     zpOffsets, 36236, 87, 'i', 1)

    def test_fgcmcalMultipleFitPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcalmultiple'

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        # These are slightly different from above due to the configuration change
        # mid-way in the separate fits.
        zpOffsets = np.array([-0.0006988655077293515,
                              0.004102597013115883])

        self._testFgcmMultiFit(instName, testName,
                               "physical_filter IN ('HSC-G', 'HSC-R', 'HSC-I')",
                               visits, zpOffsets)

    def test_fgcmcalTractPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcaltract'

        nBand = 3
        i0Std = np.array([0.08294534, 0.07877351, 0.06464688])
        i10Std = np.array([-0.000091981, -0.00061516, -0.00063434])
        i0Recon = np.array([0.07322632, 0.0689530429, 0.05600673])
        i10Recon = np.array([-5.89816122, -7.01847144, 3.62675740])

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        rawRepeatability = np.array([0.0,
                                     0.0025195920941720683,
                                     0.004095912225403857])
        filterNCalibMap = {'HSC-R': 12,
                           'HSC-I': 15}

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]
        tract = 9697

        self._testFgcmCalibrateTract(instName, testName,
                                     visits, tract, 'hsc_rings_v1',
                                     rawRepeatability, filterNCalibMap)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
