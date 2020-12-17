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
"""Class for running fgcmcal on a single tract using sourceTable_visit tables.
"""

import lsst.pipe.base as pipeBase

from .dataIds import TractCheckDataIdContainer
from .fgcmBuildStarsTable import FgcmBuildStarsTableTask
from .fgcmCalibrateTractBase import (FgcmCalibrateTractConfigBase, FgcmCalibrateTractRunner,
                                     FgcmCalibrateTractBaseTask)

__all__ = ['FgcmCalibrateTractTableConfig', 'FgcmCalibrateTractTableTask']


class FgcmCalibrateTractTableConfig(FgcmCalibrateTractConfigBase):
    """Config for FgcmCalibrateTractTable task"""
    def setDefaults(self):
        super().setDefaults()

        # For the Table version of CalibrateTract, use the associated
        # Table version of the BuildStars task.
        self.fgcmBuildStars.retarget(FgcmBuildStarsTableTask)
        # For tract mode, we set a very high effective density cut.
        self.fgcmBuildStars.densityCutMaxPerPixel = 10000


class FgcmCalibrateTractTableTask(FgcmCalibrateTractBaseTask):
    """
    Calibrate a single tract using fgcmcal, using sourceTable_visit
    input catalogs.
    """
    ConfigClass = FgcmCalibrateTractTableConfig
    RunnerClass = FgcmCalibrateTractRunner
    _DefaultName = "fgcmCalibrateTractTable"

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "sourceTable_visit",
                               help="Data ID, e.g. --id visit=6789 tract=9617",
                               ContainerClass=TractCheckDataIdContainer)

        return parser
