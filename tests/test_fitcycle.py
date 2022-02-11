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
"""Test the fitcycle configurations
"""

import unittest
import copy

import lsst.fgcmcal as fgcmcal
import lsst.utils


class FgcmcalTestFitCycleConfig(lsst.utils.tests.TestCase):
    """
    Test FgcmFitCycleConfig validation.
    """
    def test_fgcmfitcycle_config_validation(self):
        """
        Test the FgcmFitCycleConfig validation.
        """
        config = fgcmcal.FgcmFitCycleConfig()

        config.cycleNumber = 0
        config.utBoundary = 0.0
        config.latitude = 0.0
        config.outfileBase = 'None'
        config.bands = ['r', 'i']
        config.fitBands = ['r', 'i']
        config.requiredBands = ['r']
        config.colorSplitBands = ['r', 'i']
        config.superStarSubCcdDict = {'r': True, 'i': True}
        config.ccdGraySubCcdDict = {'r': True, 'i': True}
        config.expGrayPhotometricCutDict = {'r': -0.05, 'i': -0.05}
        config.expGrayHighCutDict = {'r': 0.10, 'i': 0.10}
        config.expVarGrayPhotometricCutDict = {'r': 0.0005, 'i': 0.0005}
        config.sigFgcmMaxEGrayDict = {'r': 0.05, 'i': 0.05}
        config.approxThroughputDict = {'r': 1.0, 'i': 1.0}
        config.useRepeatabilityForExpGrayCutsDict = {'r': False, 'i': False}
        config.defaultCameraOrientation = 0.0

        # Ensure that it validates
        config.validate()

        self._test_misconfig(config, 'fitBands', ['r', 'i', 'z'])
        self._test_misconfig(config, 'requiredBands', ['r', 'i', 'z'])
        self._test_misconfig(config, 'colorSplitBands', ['r', 'z'])
        self._test_misconfig(config, 'superStarSubCcdDict', {'r': True})
        self._test_misconfig(config, 'ccdGraySubCcdDict', {'r': True})
        self._test_misconfig(config, 'expGrayPhotometricCutDict', {'r': -0.05})
        self._test_misconfig(config, 'expGrayHighCutDict', {'r': 0.10})
        self._test_misconfig(config, 'expVarGrayPhotometricCutDict', {'r': 0.0005})
        self._test_misconfig(config, 'sigFgcmMaxEGrayDict', {'r': 0.05})
        self._test_misconfig(config, 'approxThroughputDict', {'r': 1.0})
        self._test_misconfig(config, 'useRepeatabilityForExpGrayCutsDict', {'r': False})

        config.doComputeDeltaAperMap = True
        self._test_misconfig(config, 'deltaAperInnerRadiusArcsec', 0.0)
        config.deltaAperInnerRadiusArcsec = 1.0
        self._test_misconfig(config, 'deltaAperOuterRadiusArcsec', 0.0)
        self._test_misconfig(config, 'deltaAperOuterRadiusArcsec', 1.0)

    def _test_misconfig(self, config, field, value):
        """
        Test misconfigured field.
        """
        config2 = copy.copy(config)
        config2.update(**{field: value})
        with self.assertRaises(ValueError):
            config2.validate()


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
