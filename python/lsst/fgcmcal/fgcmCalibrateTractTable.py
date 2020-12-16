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
import numpy as np

import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
from lsst.meas.algorithms import ReferenceObjectLoader
import lsst.afw.table as afwTable

from .dataIds import TractCheckDataIdContainer
from .fgcmBuildStarsTable import FgcmBuildStarsTableTask
from .fgcmCalibrateTractBase import (FgcmCalibrateTractConfigBase, FgcmCalibrateTractRunner,
                                     FgcmCalibrateTractBaseTask)
from .utilities import lookupStaticCalibrations

__all__ = ['FgcmCalibrateTractTableConfig', 'FgcmCalibrateTractTableTask']


class FgcmCalibrateTractTableConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument",
                                                     "tract",)):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibrations,
        isCalibration=True,
    )

    fgcmLookUpTable = connectionTypes.PrerequisiteInput(
        doc=("Atmosphere + instrument look-up-table for FGCM throughput and "
             "chromatic corrections."),
        name="fgcmLookUpTable",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    sourceSchema = connectionTypes.PrerequisiteInput(
        doc="Schema for source catalogs",
        name="src_schema",
        storageClass="SourceCatalog",
        deferLoad=True,
    )

    refCat = connectionTypes.PrerequisiteInput(
        doc="Reference catalog to use for photometric calibration",
        name="cal_ref_cat",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    source_catalogs = connectionTypes.Input(
        doc="Source table in parquet format, per visit",
        name="sourceTable_visit",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
        deferLoad=True,
        multiple=True,
    )

    calexp = connectionTypes.Input(
        doc="Calibrated exposures used for psf and metadata",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
        multiple=True,
    )

    background = connectionTypes.Input(
        doc="Calexp background model",
        name="calexpBackground",
        storageClass="Background",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
        multiple=True,
    )

    fgcmPhotoCalib = connectionTypes.Output(
        doc="Per-tract, per-visit photoCalib exposure catalogs produced from fgcm calibration",
        name="fgcmPhotoCalibTractCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "tract", "visit",),
        multiple=True,
    )

    fgcmTransmissionAtmosphere = connectionTypes.Output(
        doc="Per-visit atmosphere transmission files produced from fgcm calibration",
        name="transmission_atmosphere_fgcm_tract",
        storageClass="TransmissionCurve",
        dimensions=("instrument", "tract", "visit",),
        multiple=True,
    )

    fgcmRepeatability = connectionTypes.Output(
        doc="Per-band raw repeatability numbers in the fgcm tract calibration",
        name="fgcmRawRepeatability",
        storageClass="Catalog",
        dimensions=("instrument", "tract",),
        multiple=False,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        # The ref_dataset_name will be deprecated with Gen2
        loaderName = config.fgcmBuildStars.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name
        if config.connections.refCat != loaderName:
            raise ValueError("connections.refCat must be the same as "
                             "config.fgcmBuildStars.fgcmLoadReferenceCatalog.refObjLoader.ref_dataset_name")
        if config.fgcmOutputProducts.doReferenceCalibration:
            loaderName = config.fgcmOutputProducts.refObjLoader.ref_dataset_name
            if config.connections.refCat != loaderName:
                raise ValueError("connections.refCat must be the same as "
                                 "config.fgcmOutputProducts.refObjLoader.ref_dataset_name")

        if not config.fgcmBuildStars.doModelErrorsWithBackground:
            self.inputs.remove("background")

        if config.fgcmOutputProducts.doRefcatOutput:
            raise ValueError("FgcmCalibrateTractTableTask (Gen3) does not support doRefcatOutput")
        if not config.fgcmOutputProducts.doAtmosphereOutput:
            self.prerequisiteInputs.remove("fgcmAtmosphereParameters")
        if not config.fgcmOutputProducts.doZeropointOutput:
            self.prerequisiteInputs.remove("fgcmZeropoints")


class FgcmCalibrateTractTableConfig(FgcmCalibrateTractConfigBase, pipeBase.PipelineTaskConfig,
                                    pipelineConnections=FgcmCalibrateTractTableConnections):
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
    Calibrate a single tract using fgcmcal, using sourceTable_visit (parquet)
    input catalogs.
    """
    ConfigClass = FgcmCalibrateTractTableConfig
    RunnerClass = FgcmCalibrateTractRunner
    _DefaultName = "fgcmCalibrateTractTable"

    canMultiprocess = False

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataRefDict = butlerQC.get(inputRefs)

        self.log.info("Running with %d sourceTable_visit dataRefs", (len(dataRefDict['source_catalogs'])))

        # Run the build stars tasks
        tract = butlerQC.quantum.dataId['tract']

        calexpRefs = dataRefDict['calexp']
        calexpRefDict = {(calexpRef.dataId.byName()['visit'],
                          calexpRef.dataId.byName()['detector']):
                         calexpRef for calexpRef in calexpRefs}
        dataRefDict['calexps'] = calexpRefDict

        # And the outputs
        if self.config.fgcmOutputProducts.doZeropointOutput:
            photoCalibRefDict = {(photoCalibRef.dataId.byName()['tract'],
                                  photoCalibRef.dataId.byName()['visit']):
                                 photoCalibRef for photoCalibRef in outputRefs.fgcmPhotoCalib}
            dataRefDict['fgcmPhotoCalibs'] = photoCalibRefDict

        if self.config.fgcmOutputProducts.doAtmosphereOutput:
            atmRefDict = {(atmRef.dataId.byName()['tract'],
                           atmRef.dataId.byName()['visit']): atmRef for
                          atmRef in outputRefs.fgcmTransmissionAtmosphere}
            dataRefDict['fgcmTransmissionAtmospheres'] = atmRefDict

        if self.config.fgcmBuildStars.doReferenceMatches:
            refConfig = self.config.fgcmBuildStars.fgcmLoadReferenceCatalog.refObjLoader
            loader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                    for ref in inputRefs.refCat],
                                           refCats=butlerQC.get(inputRefs.refCat),
                                           config=refConfig,
                                           log=self.log)
            buildStarsRefObjLoader = loader
        else:
            buildStarsRefObjLoader = None

        if self.config.fgcmOutputProducts.doReferenceCalibration:
            refConfig = self.config.fgcmOutputProducts.refObjLoader
            loader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                    for ref in inputRefs.refCat],
                                           refCats=butlerQC.get(inputRefs.refCat),
                                           config=refConfig,
                                           log=self.log)
            self.fgcmOutputProducts.refObjLoader = loader

        struct = self.run(dataRefDict, tract,
                          buildStarsRefObjLoader=buildStarsRefObjLoader, butler=butlerQC)

        # Turn raw repeatability into simple catalog for persistence
        schema = afwTable.Schema()
        schema.addField('rawRepeatability', type=np.float64,
                        doc="Per-band raw repeatability in FGCM calibration.")
        repeatabilityCat = afwTable.BaseCatalog(schema)
        repeatabilityCat.resize(len(struct.repeatability))
        repeatabilityCat['rawRepeatability'][:] = struct.repeatability

        butlerQC.put(repeatabilityCat, outputRefs.fgcmRepeatability)

        return

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "sourceTable_visit",
                               help="Data ID, e.g. --id visit=6789 tract=9617",
                               ContainerClass=TractCheckDataIdContainer)

        return parser
