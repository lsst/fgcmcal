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

import lsst.utils
import lsst.daf.persistence as dafPersist

from lsst.fgcmcal.utilities import computeApertureRadiusFromSchema, computeApertureRadiusFromName


class FgcmApertureTest(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def test_fgcmApertureHsc(self):
        """
        Test computeApertureRadius for HSC.
        """
        lsst.log.setLevel("HscMapper", lsst.log.FATAL)

        butler = dafPersist.Butler(os.path.join(self.dataDir, 'hsc'))

        schema = butler.get('src_schema').schema

        self.assertRaises(RuntimeError, computeApertureRadiusFromSchema, schema, 'base_PsfFlux_instFlux')
        self.assertRaises(LookupError, computeApertureRadiusFromSchema, schema, 'not_a_field')
        self.assertEqual(computeApertureRadiusFromSchema(schema, 'slot_CalibFlux_instFlux'), 12.0)
        self.assertEqual(computeApertureRadiusFromSchema(schema,
                                                         'base_CircularApertureFlux_12_0_instFlux'), 12.0)
        self.assertEqual(computeApertureRadiusFromSchema(schema,
                                                         'base_CircularApertureFlux_4_5_instFlux'), 4.5)
        self.assertEqual(computeApertureRadiusFromName('ApFlux_12_0_instFlux'), 12.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
