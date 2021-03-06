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
"""Class for running fgcmcal on a single tract using src tables.
"""

import lsst.pipe.base as pipeBase
from lsst.jointcal.dataIds import PerTractCcdDataIdContainer

from .fgcmCalibrateTractBase import (FgcmCalibrateTractConfigBase, FgcmCalibrateTractRunner,
                                     FgcmCalibrateTractBaseTask)

__all__ = ['FgcmCalibrateTractConfig', 'FgcmCalibrateTractTask']


class FgcmCalibrateTractConfig(FgcmCalibrateTractConfigBase):
    """Config for FgcmCalibrateTract task"""
    def setDefaults(self):
        super().setDefaults()

        # For tract mode, turn off checkAllCcds optimization for
        # FgcmBuildStarsTask.
        self.fgcmBuildStars.checkAllCcds = False
        # For tract mode, we set a very high effective density cut.
        self.fgcmBuildStars.densityCutMaxPerPixel = 10000


class FgcmCalibrateTractTask(FgcmCalibrateTractBaseTask):
    """
    Calibrate a single tract using fgcmcal
    """
    ConfigClass = FgcmCalibrateTractConfig
    RunnerClass = FgcmCalibrateTractRunner
    _DefaultName = "fgcmCalibrateTract"

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "src", help="Data ID, e.g. --id visit=6789",
                               ContainerClass=PerTractCcdDataIdContainer)

        return parser
