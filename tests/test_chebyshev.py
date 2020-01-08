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
"""Test the fgcm and fgcmcal (stack) chebyshev polynomial code.

Run small tests to ensure that the two chebyshev implementations are
equivalent.  The fgcm chebyshev polynomial code
(`fgcm.fgcmUtilities.Cheb2dField`) was developed independently of the stack
`lsst.afw.math.ChebyshevBoundedField`, so as to avoid dependencies.
"""

import matplotlib
matplotlib.use("Agg")  # noqa E402

import unittest
import numpy as np

import lsst.afw.math as afwMath
import lsst.geom

import fgcm


class FgcmChebyshevTest(lsst.utils.tests.TestCase):
    """
    Test the fgcm and fgcmcal (stack) chebyshev polynomial code.
    """
    def setUp(self):
        """
        Set up test class.
        """
        self.xSize = 2000
        self.ySize = 4000
        self.nStar = 1000

        self.order = 2
        self.pars = np.zeros((self.order + 1, self.order + 1))
        self.pars[0, 0] = 1.0
        self.pars[0, 1] = 3e-3
        self.pars[0, 2] = 1e-5
        self.pars[1, 0] = 5e-3
        self.pars[1, 1] = 7e-6
        self.pars[1, 2] = 0.0
        self.pars[2, 0] = 3e-5
        self.pars[2, 1] = 0.0
        self.pars[2, 2] = 0.0

        self.order2 = 1
        self.pars2 = np.zeros((self.order2 + 1, self.order2 + 1))
        self.pars2[0, 0] = 0.98
        self.pars2[0, 1] = 5e-3
        self.pars2[1, 0] = 2e-2

    def test_chebyshev_evaluate(self, seed=1000):
        """
        Test the evaluation of chebyshev polynomials.

        Parameters
        ----------
        seed: `int`, optional
           Numpy random seed
        """
        # Set numpy seed for stability
        np.random.seed(seed=seed)

        xPos = self.xSize * np.random.rand(self.nStar)
        yPos = self.ySize * np.random.rand(self.nStar)

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Point2I(self.xSize - 1, self.ySize - 1))

        # Compute the chebyshev values using the fgcm code
        fgcmField = fgcm.fgcmUtilities.Cheb2dField(self.xSize, self.ySize,
                                                   self.pars)
        fgcmValues = fgcmField.evaluate(xPos, yPos)

        # Compute the chebyshev values using the afw code
        field = afwMath.ChebyshevBoundedField(bbox, self.pars)
        fieldValues = field.evaluate(xPos, yPos)

        self.assertFloatsAlmostEqual(fieldValues, fgcmValues, rtol=5e-15)

    def test_chebyshev_fit(self, seed=1000):
        """
        Test the fitting of chebyshev polynomials.

        Parameters
        ----------
        seed: `int`, optional
           Numpy random seed
        """
        # Set numpy seed for stability
        np.random.seed(seed=seed)

        # Generate some points to fit
        xPos = self.xSize * np.random.rand(self.nStar)
        yPos = self.ySize * np.random.rand(self.nStar)
        fgcmField = fgcm.fgcmUtilities.Cheb2dField(self.xSize, self.ySize,
                                                   self.pars)
        fgcmValues = fgcmField.evaluate(xPos, yPos)

        # Fit the points using the fgcm code
        fgcmField = fgcm.fgcmUtilities.Cheb2dField.fit(self.xSize, self.ySize,
                                                       self.order, xPos, yPos,
                                                       fgcmValues)

        # Fit the points using the afw code
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Point2I(self.xSize - 1, self.ySize - 1))

        ctrl = afwMath.ChebyshevBoundedFieldControl()
        ctrl.orderX = self.order
        ctrl.orderY = self.order
        ctrl.triangular = True
        field = afwMath.ChebyshevBoundedField.fit(bbox, xPos, yPos,
                                                  fgcmValues, ctrl)

        # Compare the fit parameters
        # The tolerance here must be looser than the application, I believe
        # because of rounding errors in the fit implementations.  But the
        # good news is that a tolerance of 1e-9 in parameters in these
        # tests yields a recovered tolerance of < 5e-15.
        self.assertFloatsAlmostEqual(fgcmField.pars, field.getCoefficients(),
                                     rtol=1e-9)

        # And compare the input and output
        fgcmValues2 = fgcmField.evaluate(xPos, yPos)
        fieldValues2 = field.evaluate(xPos, yPos)

        self.assertFloatsAlmostEqual(fgcmValues, fgcmValues2, rtol=5e-15)
        self.assertFloatsAlmostEqual(fgcmValues2, fieldValues2, rtol=5e-15)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
