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
"""Configuration for SED terms in fgcmcal.
Analogous to colorterms.
"""

from lsst.pex.config import Config, Field, ConfigDictField

__all__ = ["Sedterm", "SedtermDict", "Sedboundaryterm", "SedboundarytermDict"]


class Sedboundaryterm(Config):
    """SED boundary term for a pair of bands.

    The SED slope (in flux units) at the boundary between two bands is given by:

        S = -0.921 * (primary - secondary) / (lambda_primary - lambda_secondary)

    To construct a Sedboundaryterm, use keyword arguments:
    Sedboundaryterm(primary=primaryBandName, secondary=secondaryBandName)

    This is a subclass of Config.  This follows the form of
    `lsst.pipe.tasks.Colorterm`.
    """
    primary = Field(dtype=str, doc="name of primary band")
    secondary = Field(dtype=str, doc="name of secondary band")


class SedboundarytermDict(Config):
    """A mapping of Sedboundaryterm name to Sedterm.

    To construct a SedboundarytermDict use keyword arguments:
    SedboundarytermDict(data=dataDict)
    where dataDict is a Python dict of name: Sedterm
    For example::

        SedboundarytermDict(data={
            'gr': Sedboundaryterm(primary="g", secondary="r"),
            'ri': Sedboundaryterm(primary="r", secondary="i"),
        })

    This is a subclass of Config.  This follows the form of
    `lsst.pipe.tasks.ColortermDict`.
    """
    data = ConfigDictField(
        doc="Mapping of Sedboundary term name to Sedboundaryterm",
        keytype=str,
        itemtype=Sedboundaryterm,
        default={},
    )


class Sedterm(Config):
    """SED term for a single band.

    The SED slope (in flux units) in the middle of a band is computed either
    as an "interpolated" or "extrapolated" computation.  See Burke et al. 2018
    Appendix A (https://ui.adsabs.harvard.edu/abs/2018AJ....155...41B).

    For interpolation, with a secondary term::

       F'_nu ~ constant * (primaryTerm + secondaryTerm) / 2.0

    For interpolation, without a secondary term::

       F'_nu ~ constant * primaryTerm

    For extrapolation::

       F'_nu ~ primaryTerm + constant * (((lambda_primaryBand - lambda_secondaryBand) /
                                         (lambda_primaryBand - lambda_tertiaryBand)) *
                                         (primaryTerm - secondaryTerm))

    where primaryTerm and secondaryTerm are names from a `SedboundarytermDict`, and
    primaryBand, secondaryBand, and tertiaryBand are band names.

    To construct a Sedterm, use keyword arguments::

        Sedterm(primaryTerm=primaryTermName, secondaryTerm=secondaryTermName,
                extrapolated=False, constant=1.0)

    or::

        Sedterm(primaryTerm=primaryTermName, secondaryTerm=secondaryTermName,
                extrapolated=True, constant=1.0, primaryBand=primaryBandName,
                secondaryBand=secondaryBandName, tertiaryBand=tertiaryBandName)

    This is a subclass of Config.  This follows the form of
    `lsst.pipe.tasks.Colorterm`.
    """
    primaryTerm = Field(dtype=str, doc="Name of primary Sedboundaryterm")
    secondaryTerm = Field(dtype=str, default=None, optional=True,
                          doc="Name of secondary Sedboundaryterm")
    extrapolated = Field(dtype=bool, default=False, doc="Extrapolate to compute SED slope")
    constant = Field(dtype=float, default=1.0, doc="Adjustment constant for SED slope")
    primaryBand = Field(dtype=str, default=None, optional=True,
                        doc="Primary band name for extrapolation")
    secondaryBand = Field(dtype=str, default=None, optional=True,
                          doc="Secondary band name for extrapolation")
    tertiaryBand = Field(dtype=str, default=None, optional=True,
                         doc="Tertiary band name for extrapolation")

    def validate(self):
        Config.validate(self)
        if self.extrapolated:
            if self.primaryBand is None or \
                    self.secondaryBand is None or \
                    self.tertiaryBand is None:
                raise RuntimeError("extrapolated requires primaryBand, secondaryBand, and "
                                   "tertiaryBand are provided.")


class SedtermDict(Config):
    """A mapping of bands to Sedterms.

    To construct a SedtermDict use keyword arguments::

        SedtermDict(data=dataDict)

    where dataDict is a Python dict of band to Sedterm
    For example::

        SedtermDict(data={
            'g': Sedterm(primaryTerm='gr', secondaryTerm='ri', extrapolated=True, constant=0.25,
                         primaryBand='g', secondaryBand='r', tertiaryBand='i'),
            'r': Sedterm(primaryTerm='gr', secondaryTerm='ri', extrapolated=False)
        })

    This is a subclass of Config.  This follows the form of
    `lsst.pipe.tasks.ColortermDict`.
    """
    data = ConfigDictField(
        doc="Mapping of band name to Sedterm",
        keytype=str,
        itemtype=Sedterm,
        default={},
    )
