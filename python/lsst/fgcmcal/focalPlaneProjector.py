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
"""A class to project the focal plane in arbitrary rotations for fgcm.

This file contains a class used by fgcm ...
"""
from functools import lru_cache
import warnings
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.geom as geom
from lsst.obs.base import createInitialSkyWcs

__all__ = ['FocalPlaneProjector']


class FocalPlaneProjector(object):
    """
    Class to project the focal plane onto the sky.

    Parameters
    ----------
    camera : `lsst.afw.cameraGeom.Camera`
        Camera from the butler.
    defaultOrientation : `int`
        Default camera orientation in degrees.  This angle is the position
        angle of the focal plane +Y with respect to north.
    useScienceDetectors : `bool`, optional
        Use only science detectors in projector?
    """
    def __init__(self, camera, defaultOrientation, useScienceDetectors=False):
        self.camera = camera

        # Put the reference boresight at the equator to avoid cos(dec) problems.
        self.boresight = geom.SpherePoint(180.0*geom.degrees, 0.0*geom.degrees)
        self.flipX = False
        self.defaultOrientation = int(defaultOrientation) % 360
        self.useScienceDetectors = useScienceDetectors

    def _makeWcsDict(self, orientation):
        """
        Make a dictionary of WCSs at the reference boresight position.

        Parameters
        ----------
        orientation : `int`
            Orientation in degrees.  This angle is the position
            angle of the focal plane +Y with respect to north.

        Returns
        -------
        wcsDict : `dict`
            Dictionary of WCS, with the detector id as the key.
        """
        _orientation = orientation*geom.degrees

        visitInfo = afwImage.VisitInfo(boresightRaDec=self.boresight,
                                       boresightRotAngle=_orientation,
                                       rotType=afwImage.RotType.SKY)

        wcsDict = {}

        for detector in self.camera:
            if self.useScienceDetectors:
                if not detector.getType() == afwCameraGeom.DetectorType.SCIENCE:
                    continue

            detectorId = detector.getId()
            wcsDict[detectorId] = createInitialSkyWcs(visitInfo, detector, self.flipX)

        return wcsDict

    def __call__(self, orientation, nstep=100, use_cache=True):
        """
        Make a focal plane projection mapping for use with fgcm.

        Parameters
        ----------
        orientation : `float` or `int`
            Camera orientation in degrees.  This angle is the position
            angle of the focal plane +Y with respect to north.
        nstep : `int`
            Number of steps in x/y per detector for the mapping.
        use_cache : `bool`, optional
            Use integerized cached lookup.

        Returns
        -------
        projectionMapping : `np.ndarray`
            A projection mapping object with x, y, x_size, y_size,
            delta_ra_cent, delta_dec_cent, delta_ra, delta_dec for
            each detector id.
        """
        if not np.isfinite(orientation):
            warnings.warn('Encountered non-finite orientation; using default.')
            _orientation = self.defaultOrientation
        else:
            _orientation = orientation % 360

        if use_cache:
            _orientation = int(_orientation)

            return self._compute_cached_projection(int(_orientation), nstep=nstep)
        else:
            return self._compute_projection(_orientation, nstep=nstep)

    @lru_cache(maxsize=360)
    def _compute_cached_projection(self, orientation, nstep=50):
        """
        Compute the focal plane projection, with caching.

        Parameters
        ----------
        orientation : `int`
            Camera orientation in degrees. This angle is the position
            angle of the focal plane +Y with respect to north.
        nstep : `int`
            Number of steps in x/y per detector for the mapping.

        Returns
        -------
        projectionMapping : `np.ndarray`
            A projection mapping object with x, y, x_size, y_size,
            delta_ra_cent, delta_dec_cent, delta_ra, delta_dec for
            each detector id.
        """
        return self._compute_projection(orientation, nstep=nstep)

    def _compute_projection(self, orientation, nstep=50):
        """
        Compute the focal plane projection.

        Parameters
        ----------
        orientation : `float` or `int`
            Camera orientation in degrees. This angle is the position
            angle of the focal plane +Y with respect to north.
        nstep : `int`
            Number of steps in x/y per detector for the mapping.

        Returns
        -------
        projectionMapping : `np.ndarray`
            A projection mapping object with x, y, x_size, y_size,
            delta_ra_cent, delta_dec_cent, delta_ra, delta_dec for
            each detector id.
        """
        wcsDict = self._makeWcsDict(orientation)

        # Need something for the max detector ...
        deltaMapper = np.zeros(
            len(wcsDict),
            dtype=[
                ('id', 'i4'),
                ('x', 'f8', nstep**2),
                ('y', 'f8', nstep**2),
                ('x_size', 'i4'),
                ('y_size', 'i4'),
                ('delta_ra_cent', 'f8'),
                ('delta_dec_cent', 'f8'),
                ('delta_ra', 'f8', nstep**2),
                ('delta_dec', 'f8', nstep**2)
            ],
        )

        for detector in self.camera:
            if self.useScienceDetectors:
                if not detector.getType() == afwCameraGeom.DetectorType.SCIENCE:
                    continue

            detectorId = detector.getId()

            deltaMapper['id'][detectorId] = detectorId

            xSize = detector.getBBox().getMaxX()
            ySize = detector.getBBox().getMaxY()

            xValues = np.linspace(0.0, xSize, nstep)
            yValues = np.linspace(0.0, ySize, nstep)

            deltaMapper['x'][detectorId, :] = np.repeat(xValues, yValues.size)
            deltaMapper['y'][detectorId, :] = np.tile(yValues, xValues.size)
            deltaMapper['x_size'][detectorId] = xSize
            deltaMapper['y_size'][detectorId] = ySize

            radec = wcsDict[detector.getId()].pixelToSkyArray(deltaMapper['x'][detectorId, :],
                                                              deltaMapper['y'][detectorId, :],
                                                              degrees=True)

            deltaMapper['delta_ra'][detectorId, :] = radec[0] - self.boresight.getRa().asDegrees()
            deltaMapper['delta_dec'][detectorId, :] = radec[1] - self.boresight.getDec().asDegrees()

            detCenter = wcsDict[detector.getId()].pixelToSky(detector.getCenter(afwCameraGeom.PIXELS))
            deltaMapper['delta_ra_cent'][detectorId] = (detCenter.getRa()
                                                        - self.boresight.getRa()).asDegrees()
            deltaMapper['delta_dec_cent'][detectorId] = (detCenter.getDec()
                                                         - self.boresight.getDec()).asDegrees()

        return deltaMapper
