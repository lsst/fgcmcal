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
"""Test the fgcmcal FocalPlaneProjector code with testdata_jointcal/hsc.
"""
import unittest
import os
import tempfile
import numpy as np
import warnings

import lsst.utils
import lsst.daf.butler as dafButler
import lsst.fgcmcal as fgcmcal

import fgcmcalTestBase


ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmFocalPlaneProjectorTestHsc(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest('testdata_jointcal not setup.')
        try:
            lsst.utils.getPackageDir('obs_subaru')
        except LookupError:
            raise unittest.SkipTest("obs_subaru not setup")

        cls.testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")

        cls._importRepository('lsst.obs.subaru.HyperSuprimeCam',
                              os.path.join(cls.dataDir, 'hsc/repo'),
                              os.path.join(cls.dataDir, 'hsc', 'exports.yaml'))

    def test_focalPlaneProjector(self):
        """
        Test focal plane projector code for HSC.
        """
        butler = dafButler.Butler(os.path.join(self.testDir, 'testrepo'),
                                  instrument='HSC',
                                  collections=['HSC/calib/unbounded', 'HSC/testdata'])
        camera = butler.get('camera', instrument='HSC')

        visit = 36236
        detector = 87

        visitInfo = butler.get('calexp.visitInfo', visit=visit, detector=detector)
        wcs = butler.get('calexp.wcs', visit=visit, detector=detector)
        rotAngle = visitInfo.getBoresightRotAngle().asDegrees()

        defaultRotation = 270.0
        self.assertAlmostEqual(rotAngle, defaultRotation)

        focalPlaneProjector = fgcmcal.FocalPlaneProjector(camera, defaultRotation)
        self.assertAlmostEqual(focalPlaneProjector.defaultOrientation, defaultRotation)

        focalPlaneProjector2 = fgcmcal.FocalPlaneProjector(camera, -90.0)
        self.assertAlmostEqual(focalPlaneProjector2.defaultOrientation, defaultRotation)

        # Load a delta mapper at the default angle
        deltaMapper = focalPlaneProjector(int(defaultRotation))

        # Spot check relative orientations of some of the ccds
        # This is based on
        # https://subarutelescope.org/Observing/Instruments/HSC/CCDPosition_20170212.png
        # The goal here is to check that North is Up, South is Down,
        # East is Left, and West is Right.
        self.assertLess(deltaMapper['delta_ra_cent'][15], deltaMapper['delta_ra_cent'][10])
        self.assertLess(deltaMapper['delta_ra_cent'][95], deltaMapper['delta_ra_cent'][90])
        self.assertGreater(deltaMapper['delta_ra_cent'][46], 0.0)

        # Check that loading at NaN gives the same default (with a warning)
        self.assertWarns(UserWarning, focalPlaneProjector, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            deltaMapperDefault = focalPlaneProjector(np.nan)

        np.testing.assert_array_almost_equal(deltaMapperDefault['delta_ra_cent'],
                                             deltaMapper['delta_ra_cent'])
        np.testing.assert_array_almost_equal(deltaMapperDefault['delta_dec_cent'],
                                             deltaMapper['delta_dec_cent'])

        # Compare the ra/dec to x/y mapping for the one detector from the wcs.
        wcsRa, wcsDec = wcs.pixelToSkyArray(deltaMapper['x'][detector, :],
                                            deltaMapper['y'][detector, :],
                                            degrees=True)
        boresightRa = visitInfo.getBoresightRaDec().getRa().asDegrees()
        boresightDec = visitInfo.getBoresightRaDec().getDec().asDegrees()

        wcsDeltaRa = (wcsRa - boresightRa)*np.cos(np.deg2rad(boresightDec))
        wcsDeltaDec = wcsDec - boresightDec

        meanDeltaRa = np.mean(wcsDeltaRa - deltaMapper['delta_ra'][detector, :])
        meanDeltaDec = np.mean(wcsDeltaDec - deltaMapper['delta_dec'][detector, :])

        # Check that these offsets are reasonable (they will not be zero
        # because the boresight is not updated after the WCS is fit...)
        self.assertLess(np.abs(meanDeltaRa), 0.002)
        self.assertLess(np.abs(meanDeltaDec), 0.002)

        self.assertLess(np.abs(wcsDeltaRa - deltaMapper['delta_ra'][detector, :] - meanDeltaRa).max(),
                        0.0005)
        self.assertLess(np.abs(wcsDeltaDec - deltaMapper['delta_dec'][detector, :] - meanDeltaDec).max(),
                        0.0005)

        # Check the sizes
        # The projected size of the ccds varies with radius over the large
        # HSC field-of-view. Empirically, the x size is between 0.07 and 0.10 deg
        # and the y size is between 0.17 and 0.20 deg for the non-rotated CCDs.
        # This test checks that the orientations of the CCDs are as expected for
        # rotated/non-rotated CCDs (relative size of RA/DEC), and that the size
        # is roughly correct.

        rotatedDetectors = [100, 101, 102, 103]

        for i, detectorId in enumerate(deltaMapper['id']):
            ra_size = np.max(deltaMapper['delta_ra'][i, :]) - np.min(deltaMapper['delta_ra'][i, :])
            dec_size = np.max(deltaMapper['delta_dec'][i, :]) - np.min(deltaMapper['delta_dec'][i, :])
            if detectorId in rotatedDetectors:
                self.assertLess(ra_size, dec_size)
                self.assertGreater(dec_size, 0.17)
                self.assertLess(dec_size, 0.20)
                self.assertGreater(ra_size, 0.07)
                self.assertLess(ra_size, 0.10)
            else:
                self.assertGreater(ra_size, dec_size)
                self.assertGreater(ra_size, 0.17)
                self.assertLess(ra_size, 0.20)
                self.assertGreater(dec_size, 0.07)
                self.assertLess(dec_size, 0.10)

        # And spin it around, testing the reverse of the test above.
        deltaMapperRot = focalPlaneProjector(defaultRotation + 180.0)
        self.assertGreater(deltaMapperRot['delta_ra_cent'][15], deltaMapperRot['delta_ra_cent'][10])
        self.assertGreater(deltaMapperRot['delta_ra_cent'][95], deltaMapperRot['delta_ra_cent'][90])
        self.assertLess(deltaMapperRot['delta_ra_cent'][46], 0.0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
