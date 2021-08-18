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
"""Class for checking sourceTable_visit dataIds
"""

import lsst.pipe.base as pipeBase
import lsst.log as lsstLog


class TractCheckDataIdContainer(pipeBase.DataIdContainer):
    """A version of lsst.pipe.base.DataIdContainer that ensures that tract
    is additionally set with sourceTable_visit catalogs.
    """
    # TODO: Move to pipe_tasks with DM-25778

    def castDataIds(self, butler):
        """Validate data IDs and cast them to the correct type
        (modify self.idList in place).

        This code casts the values in the data IDs dicts in `self.idList`
        to the type required by the butler. Data IDs are read from the
        command line as `str`, but the butler requires some values to be
        other types. For example "visit" values should be `int` (which
        is defined by the templates).

        This is taken from lsst.pipe.base.ArgumentParser.castDataIds(),
        adding in a check that the tract is an int, which must be done
        explicitly because it is not in the template.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
            Data butler.
        """
        datasetType = 'sourceTable_visit'
        try:
            idKeyTypeDict = butler.getKeys(datasetType=datasetType, level=self.level)
        except KeyError as e:
            msg = f"Cannot get keys for datasetType {datasetType} at level {self.level}"
            raise KeyError(msg) from e

        idKeyTypeDict = idKeyTypeDict.copy()
        idKeyTypeDict['tract'] = int

        log = None

        for dataDict in self.idList:
            for key, strVal in dataDict.items():
                try:
                    keyType = idKeyTypeDict[key]
                except KeyError:
                    # OK, assume that it's a valid key and guess that it's a string
                    keyType = str

                    if log is None:
                        log = lsstLog.Log.getDefaultLogger()
                    log.warning("Unexpected ID %s; guessing type is \"%s\"",
                                key, 'str' if keyType == str else keyType)
                    idKeyTypeDict[key] = keyType

                if keyType != str:
                    try:
                        castVal = keyType(strVal)
                    except Exception:
                        raise TypeError(f"Cannot cast value {strVal!r} to {keyType} for ID key {key}")
                    dataDict[key] = castVal

    def makeDataRefList(self, namespace):
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")

        tracts = set()
        for dataId in self.idList:
            if "tract" not in dataId:
                raise RuntimeError("Must set tract for dataId")
            tracts.add(dataId['tract'])

            self.refList.append(namespace.butler.dataRef(datasetType=self.datasetType,
                                                         dataId=dataId))

        if len(tracts) > 1:
            raise RuntimeError("Can only run a single tract")
