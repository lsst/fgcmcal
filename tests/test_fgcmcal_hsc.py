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

# Ensure that matplotlib doesn't try to open a display during testing.
import matplotlib
matplotlib.use("Agg")

import lsst.utils  # noqa: E402
import lsst.pipe.tasks  # noqa: E402
import lsst.daf.butler.cli.cliLog  # noqa: E402

import fgcmcalTestBase  # noqa: E402


ROOT = os.path.abspath(os.path.dirname(__file__))

I0STD = [0.08294534, 0.07877351, 0.06464688]
I10STD = [-0.000091981, -0.00061516, -0.00063434]
I0RECON = [0.07322179342588758, 0.0689530429, 0.05600673]
I10RECON = [-4.490243571049125, -7.01786443508, 3.62738180611]


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

    def test_fgcmcalPipelineBuildFromTable(self):
        """Test running the full pipeline, using older association code.

        This test uses the FgcmBuildStarsFromTableTask instead of the new
        FgcmBuildFromIsolatedStarsTask.
        """
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcalpipe'

        nBand = 3
        i0Std = np.array(I0STD)
        i10Std = np.array(I10STD)
        i0Recon = np.array(I0RECON)
        i10Recon = np.array(I10RECON)

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
        nStdStars = 237
        nPlots = 59

        # We need an extra config file to turn off parquet format.
        extraConfigFile = os.path.join(self.testDir, "turn_off_parquet.py")
        with open(extraConfigFile, "w") as f:
            f.write("config.useParquetCatalogFormat = False\n")

        self._testFgcmFitCycle(instName, testName,
                               0, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True,
                               extraConfig=extraConfigFile)
        self._testFgcmFitCycle(instName, testName,
                               1, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=True,
                               extraConfig=extraConfigFile)

        # We need to create an extra config file to turn on "sub-ccd gray" for testing.
        # We also want to exercise the code path setting useExposureReferenceOffset = False.
        extraConfigFile = os.path.join(self.testDir, "cycle03_patch_config.py")
        with open(extraConfigFile, "w") as f:
            f.write("config.useParquetCatalogFormat = False\n")
            f.write("config.isFinalCycle = True\n")
            f.write("config.ccdGraySubCcdDict = {'g': True, 'r': True, 'i': True}\n")
            f.write("config.useExposureReferenceOffset = False")

        self._testFgcmFitCycle(instName, testName,
                               2, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots,
                               extraConfig=extraConfigFile)

    def test_fgcmcalPipeline(self):
        """Test running the full pipeline, using new isolated star association code.

        This test uses the FgcmBuildFromIsolatedStarsTask instead of the old
        FgcmBuildStarsFromTableTask.
        """
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcalpipe'

        nBand = 3
        i0Std = np.array(I0STD)
        i10Std = np.array(I10STD)
        i0Recon = np.array(I0RECON)
        i10Recon = np.array(I10RECON)

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        nStar = 295
        nObs = 1808

        self._testFgcmBuildFromIsolatedStars(
            instName,
            testName,
            "physical_filter IN ('HSC-G', 'HSC-R', 'HSC-I')",
            visits,
            nStar,
            nObs,
        )

        nZp = 1120
        nGoodZp = 27
        nOkZp = 27
        nBadZp = 1093
        nStdStars = 227
        nPlots = 59

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

        zpOffsets = np.array([-0.001339632086455822,
                              0.005496968515217304])

        self._testFgcmOutputProducts(instName, testName,
                                     zpOffsets, 36236, 87, 'i', 1, 'hsc_rings_v1')

        # Test a single detector illumination correction.
        self._testFgcmOutputIlluminationCorrection(instName, testName, 51)

    def test_fgcmcalMultipleFitPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcalmultiple'

        nBand = 3
        i0Std = np.array(I0STD)
        i10Std = np.array(I10STD)
        i0Recon = np.array(I0RECON)
        i10Recon = np.array(I10RECON)

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        # These are slightly different from above due to the configuration change
        # mid-way in the separate fits.
        zpOffsets = np.array([-0.0015136072179302573,
                              0.0019038696773350239])

        self._testFgcmMultiFit(instName, testName,
                               "physical_filter IN ('HSC-G', 'HSC-R', 'HSC-I') and skymap='hsc_rings_v1'",
                               visits, zpOffsets, 118, 57)

    def test_fgcmcalTractPipeline(self):
        # Set numpy seed for stability
        np.random.seed(seed=1000)

        instName = 'HSC'
        testName = 'testfgcmcaltract'

        nBand = 3
        i0Std = np.array(I0STD)
        i10Std = np.array(I10STD)
        i0Recon = np.array(I0RECON)
        i10Recon = np.array(I10RECON)

        self._testFgcmMakeLut(instName, testName,
                              nBand, i0Std, i0Recon, i10Std, i10Recon)

        rawRepeatability = np.array([0.0,
                                     0.013999312239597281,
                                     0.005144248063188913])
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
