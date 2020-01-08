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
"""Test the fgcmcal computeCcdOffsets code with testdata_jointcal.
"""

import unittest
import os

import lsst.utils
import lsst.daf.persistence as dafPersist

import lsst.fgcmcal as fgcmcal


class FgcmCcdOffsetsTest(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def test_fgcmCcdOffsetsHsc(self):
        """
        Test computation of ccd offsets for HSC.
        """

        lsst.log.setLevel("HscMapper", lsst.log.FATAL)

        butler = dafPersist.Butler(os.path.join(self.dataDir, 'hsc'))

        visit = 903986
        ccd = 16

        visitInfo = butler.get('calexp_visitInfo', visit=visit, ccd=ccd)
        rotAngle = visitInfo.getBoresightRotAngle().asDegrees()

        camera = butler.get('camera')

        ccdOffsets = fgcmcal.utilities.computeCcdOffsets(camera, rotAngle)

        # Spot check relative orientations of some of the ccds
        # This is based on
        # https://subarutelescope.org/Observing/Instruments/HSC/CCDPosition_20170212.png
        # The goal here is to check that North is Up, South is Down,
        # East is Left, and West is Right.
        self.assertLess(ccdOffsets['DELTA_RA'][15], ccdOffsets['DELTA_RA'][10])
        self.assertLess(ccdOffsets['DELTA_RA'][95], ccdOffsets['DELTA_RA'][90])
        self.assertGreater(ccdOffsets['DELTA_RA'][46], 0.0)

        self.assertLess(ccdOffsets['DELTA_DEC'][15], ccdOffsets['DELTA_DEC'][95])
        self.assertLess(ccdOffsets['DELTA_DEC'][10], ccdOffsets['DELTA_DEC'][90])
        self.assertGreater(ccdOffsets['DELTA_DEC'][97], 0.0)

        # Check the sizes
        # The projected size of the ccds varies with radius over the large
        # HSC field-of-view. Empirically, the x size is between 0.07 and 0.10 deg
        # and the y size is between 0.17 and 0.20 deg for the non-rotated CCDs.
        # This test checks that the orientations of the CCDs are as expected for
        # rotated/non-rotated CCDs (relative size of RA/DEC), and that the size
        # is roughly correct.  Because these values are only used for visualization
        # in the fgcm code, this does not need to be perfect.

        rotatedCcds = [100, 101, 102, 103]

        for i, ccd in enumerate(ccdOffsets['CCDNUM']):
            if ccd in rotatedCcds:
                self.assertLess(ccdOffsets['RA_SIZE'][i], ccdOffsets['DEC_SIZE'][i])
                self.assertGreater(ccdOffsets['DEC_SIZE'][i], 0.17)
                self.assertLess(ccdOffsets['DEC_SIZE'][i], 0.20)
                self.assertGreater(ccdOffsets['RA_SIZE'][i], 0.07)
                self.assertLess(ccdOffsets['RA_SIZE'][i], 0.10)
            else:
                self.assertGreater(ccdOffsets['RA_SIZE'][i], ccdOffsets['DEC_SIZE'][i])
                self.assertGreater(ccdOffsets['RA_SIZE'][i], 0.17)
                self.assertLess(ccdOffsets['RA_SIZE'][i], 0.20)
                self.assertGreater(ccdOffsets['DEC_SIZE'][i], 0.07)
                self.assertLess(ccdOffsets['DEC_SIZE'][i], 0.10)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
