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
"""Convert gen2 to gen3 photocalib files.
"""
# import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.utils as utils


__all__ = ['PhotoCalibConsolidateGen2Gen3Config',
           'PhotoCalibConsolidateGen2Gen3Connections',
           'PhotoCalibConsolidateGen2Gen3Task',
           'SkyWcsConsolidateGen2Gen3Config',
           'SkyWcsConsolidateGen2Gen3Connections',
           'SkyWcsConsolidateGen2Gen3Task']


class PhotoCalibConsolidateGen2Gen3Connections(pipeBase.PipelineTaskConnections,
                                               dimensions=("instrument", "visit")):
    photoCalibList = connectionTypes.Input(
        doc="Per-detector photocalibs from fgcm",
        name="fgcm_photoCalib",
        storageClass="PhotoCalib",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
        multiple=True,
    )
    photoCalibGlobalCatalog = connectionTypes.Output(
        doc="Global per-visit photocalibs.",
        name="fgcmPhotoCalibCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
    )


class PhotoCalibConsolidateGen2Gen3Config(pipeBase.PipelineTaskConfig,
                                          pipelineConnections=PhotoCalibConsolidateGen2Gen3Connections):
    def validate(self):
        super().validate()


class PhotoCalibConsolidateGen2Gen3Task(pipeBase.PipelineTask):
    """Consolidate gen2 photocalibs into gen3 photocalibs."""
    ConfigClass = PhotoCalibConsolidateGen2Gen3Config
    _DefaultName = "photoCalibConsolidateGen2Gen3"

    def __init__(self, butler=None, **kwargs):
        super().__init__(**kwargs)

    @utils.inheritDoc(pipeBase.PipelineTask)
    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        visit = butlerQC.quantum.dataId['visit']

        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField('visit', type='L', doc='visit number')

        metadata = dafBase.PropertyList()
        metadata.add("COMMENT", "Catalog id is detector id, sorted")
        metadata.add("COMMENT", "Only detectors with data have entries")

        photoCalibCat = afwTable.ExposureCatalog(schema)
        photoCalibCat.setMetadata(metadata)
        photoCalibCat.reserve(len(inputRefs.photoCalibList))

        photoCalibList = butlerQC.get(inputRefs.photoCalibList)
        for dataRef in photoCalibList:
            detector = dataRef.dataId['detector']
            photoCalib = dataRef.get()
            rec = photoCalibCat.addNew()
            rec['id'] = detector
            rec['visit'] = visit
            rec.setPhotoCalib(photoCalib)

        photoCalibCat.sort()

        butlerQC.put(photoCalibCat, outputRefs.photoCalibGlobalCatalog)


class SkyWcsConsolidateGen2Gen3Connections(pipeBase.PipelineTaskConnections,
                                           dimensions=("instrument", "visit",
                                                       "skymap", "tract")):
    skyWcsList = connectionTypes.Input(
        doc="Per-tract, per-detector wcs calibrations.",
        name="jointcal_wcs",
        storageClass="Wcs",
        dimensions=("instrument", "visit", "detector", "tract", "skymap"),
        deferLoad=True,
        multiple=True,
    )
    skyWcsTractCatalog = connectionTypes.Output(
        doc="Per-tract, per-visit wcs calibrations.",
        name="jointcalSkyWcsCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit", "tract"),
    )


class SkyWcsConsolidateGen2Gen3Config(pipeBase.PipelineTaskConfig,
                                      pipelineConnections=SkyWcsConsolidateGen2Gen3Connections):
    def validate(self):
        super().validate()


class SkyWcsConsolidateGen2Gen3Task(pipeBase.PipelineTask):
    """Consolidate gen2 skywcss into gen3 skywcss."""
    ConfigClass = SkyWcsConsolidateGen2Gen3Config
    _DefaultName = "skyWcsConsolidateGen2Gen3"

    def __init__(self, butler=None, **kwargs):
        super().__init__(**kwargs)

    @utils.inheritDoc(pipeBase.PipelineTask)
    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        visit = butlerQC.quantum.dataId['visit']

        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField('visit', type='L', doc='visit number')

        metadata = dafBase.PropertyList()
        metadata.add("COMMENT", "Catalog id is detector id, sorted")
        metadata.add("COMMENT", "Only detectors with data have entries")

        skyWcsCat = afwTable.ExposureCatalog(schema)
        skyWcsCat.setMetadata(metadata)
        skyWcsCat.reserve(len(inputRefs.skyWcsList))

        skyWcsList = butlerQC.get(inputRefs.skyWcsList)
        for dataRef in skyWcsList:
            detector = dataRef.dataId['detector']
            skyWcs = dataRef.get()
            rec = skyWcsCat.addNew()
            rec['id'] = detector
            rec['visit'] = visit
            rec.setWcs(skyWcs)

        skyWcsCat.sort()

        butlerQC.put(skyWcsCat, outputRefs.skyWcsTractCatalog)
