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
"""Test the fgcmcal computeApertureRadius code with testdata_jointcal.
"""

import unittest
import os
import tempfile

import lsst.daf.butler
import lsst.utils

from lsst.fgcmcal.utilities import computeApertureRadiusFromDataRef, computeApertureRadiusFromName

import fgcmcalTestBase


ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmApertureTestHsc(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
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

    def test_fgcmAperture(self):
        """
        Test computeApertureRadius for HSC.
        """
        butler = lsst.daf.butler.Butler(self.repo, instrument='HSC', collections=['HSC/testdata'])

        dataHandle = butler.getDeferred('src', visit=34648, detector=51)

        self.assertRaises(RuntimeError, computeApertureRadiusFromDataRef, dataHandle, 'base_PsfFlux_instFlux')
        self.assertRaises(RuntimeError, computeApertureRadiusFromDataRef, dataHandle, 'not_a_field')
        self.assertEqual(computeApertureRadiusFromDataRef(dataHandle, 'slot_CalibFlux_instFlux'), 12.0)
        self.assertEqual(computeApertureRadiusFromDataRef(dataHandle,
                                                          'base_CircularApertureFlux_12_0_instFlux'),
                         12.0)
        self.assertEqual(computeApertureRadiusFromDataRef(dataHandle,
                                                          'base_CircularApertureFlux_4_5_instFlux'),
                         4.5)

        self.assertEqual(computeApertureRadiusFromName('ApFlux_12_0_instFlux'), 12.0)
        self.assertEqual(computeApertureRadiusFromName('ApFlux_4_5_instFlux'), 4.5)
        self.assertEqual(computeApertureRadiusFromName('apFlux_12_0_instFlux'), 12.0)
        self.assertEqual(computeApertureRadiusFromName('apFlux_4_5_instFlux'), 4.5)
        self.assertRaises(RuntimeError, computeApertureRadiusFromName, 'not_a_field')


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
