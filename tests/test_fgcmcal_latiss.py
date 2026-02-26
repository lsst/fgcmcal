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
"""Test the fgcmcal code with testdata_jointcal/latiss.

Run test suite on fgcmcal using LATISS data from testdata_jointcal.
"""
import unittest
import os
import tempfile
import numpy as np

# Need to import pyproj to prevent file handle leakage since importing
# pyproj automatically opens proj.db and never closes it. We can not wait
# for some dependent code to import it whilst the test is running since then
# the leak checker will think it is a leak.
import pyproj  # noqa: F401

# Ensure that matplotlib doesn't try to open a display during testing.
import matplotlib
matplotlib.use("Agg")

import lsst.utils  # noqa: E402
import lsst.pipe.tasks  # noqa: E402
import lsst.daf.butler.cli.cliLog  # noqa: E402

import fgcmcalTestBase  # noqa: E402


ROOT = os.path.abspath(os.path.dirname(__file__))

I0STD = [0.0, 0.0, 0.0]
I10STD = [0.0, 0.0, 0.0]
I0RECON = [0.14882544833566255, 0.1082095952782785, 0.10403867395420009]
I10RECON = [-4.1568350444828335, -1.0020576281832512, -1.2088515357441083]


class FgcmcalTestLatiss(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")
        try:
            lsst.utils.getPackageDir('obs_lsst')
        except LookupError:
            raise unittest.SkipTest("obs_lsst not setup")

        lsst.daf.butler.cli.cliLog.CliLog.initLog(longlog=False)

        cls.testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")

        cls._importRepository('lsst.obs.lsst.Latiss',
                              os.path.join(cls.dataDir, 'latiss/testdata'),
                              os.path.join(cls.dataDir, 'latiss', 'exports.yaml'))

    def test_fgcmcalPipeline(self):
        """Test running the full pipeline, using isolated star association code.
        """
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'LATISS'
        testName = 'testfgcmcalpipe'

        nBand = 3
        i0Std = np.array(I0STD)
        i10Std = np.array(I10STD)
        i0Recon = np.array(I0RECON)
        i10Recon = np.array(I10RECON)

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [
            2023051100320,
            2023051100357,
            2023051100390,
            2023051100406,
            2023051100448,
            2023051100454,
            2023051100278,
            2023051100473,
            2023051100263,
            2023051100509,
            2023051100304,
            2023051100431,
            2023051100547,
            2023051100379,
            2023051100495,
            2023051100489,
            2023051100401,
            2023051100280,
            2023051100303,
            2023051100508,
        ]

        nStar = 54
        nObs = 301

        self._testFgcmBuildFromIsolatedStars(
            instName,
            testName,
            "band IN ('g', 'r', 'i')",
            visits,
            nStar,
            nObs,
            refcatCollection="refcats/DM-33444",
        )

        nZp = 20
        nGoodZp = 13
        nOkZp = 13
        nBadZp = 7
        nStdStars = 48
        nPlots = 52

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

        zpOffsets = np.array([0.025642290711402893,
                              -0.001271035522222519])

        self._testFgcmOutputProducts(
            instName,
            testName,
            zpOffsets,
            2023051100320,
            0,
            'r',
            1,
            testSrc=False,
        )

        self._testFgcmOutputIlluminationCorrection(instName, testName, 0)

    def test_fgcmcalMultipleFitPipeline(self):
        np.random.seed(seed=1000)

        instName = 'LATISS'
        testName = 'testfgcmcalmultiple'

        nBand = 3
        i0Std = np.array(I0STD)
        i10Std = np.array(I10STD)
        i0Recon = np.array(I0RECON)
        i10Recon = np.array(I10RECON)

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [
            2023051100320,
            2023051100357,
            2023051100390,
            2023051100406,
            2023051100448,
            2023051100454,
            2023051100278,
            2023051100473,
            2023051100263,
            2023051100509,
            2023051100304,
            2023051100431,
            2023051100547,
            2023051100379,
            2023051100495,
            2023051100489,
            2023051100401,
            2023051100280,
            2023051100303,
            2023051100508,
        ]

        # These are slightly different from above due to the configuration change
        # mid-way in the separate fits.
        zpOffsets = np.array([0.00999188981950283,
                              -0.009526489302515984])

        self._testFgcmMultiFit(
            instName,
            testName,
            "band IN ('g', 'r', 'i')",
            visits,
            zpOffsets,
            58,
            50,
            refcatCollection="refcats/DM-33444",
        )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
