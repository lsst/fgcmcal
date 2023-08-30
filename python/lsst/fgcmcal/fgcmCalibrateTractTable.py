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
from lsst.meas.algorithms import ReferenceObjectLoader, LoadReferenceObjectsConfig
import lsst.afw.table as afwTable

from .fgcmBuildStarsTable import FgcmBuildStarsTableTask
from .fgcmCalibrateTractBase import (FgcmCalibrateTractConfigBase,
                                     FgcmCalibrateTractBaseTask)

__all__ = ['FgcmCalibrateTractTableConfig', 'FgcmCalibrateTractTableTask']


class FgcmCalibrateTractTableConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument",
                                                     "tract",)):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
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

    sourceSchema = connectionTypes.InitInput(
        doc="Schema for source catalogs",
        name="src_schema",
        storageClass="SourceCatalog",
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

    visitSummary = connectionTypes.Input(
        doc="Per-visit summary statistics table",
        name="visitSummary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
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

        if not config.fgcmBuildStars.doModelErrorsWithBackground:
            self.inputs.remove("background")

        if not config.fgcmOutputProducts.doAtmosphereOutput:
            self.prerequisiteInputs.remove("fgcmAtmosphereParameters")
        if not config.fgcmOutputProducts.doZeropointOutput:
            self.prerequisiteInputs.remove("fgcmZeropoints")

    def getSpatialBoundsConnections(self):
        return ("visitSummary",)


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
    _DefaultName = "fgcmCalibrateTractTable"

    canMultiprocess = False

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(initInputs=initInputs, **kwargs)
        if initInputs is not None:
            self.sourceSchema = initInputs["sourceSchema"].schema

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        handleDict = butlerQC.get(inputRefs)

        self.log.info("Running with %d sourceTable_visit handles", (len(handleDict['source_catalogs'])))

        # Run the build stars tasks
        tract = butlerQC.quantum.dataId['tract']

        handleDict['sourceSchema'] = self.sourceSchema

        sourceTableHandles = handleDict['source_catalogs']
        sourceTableHandleDict = {sourceTableHandle.dataId['visit']: sourceTableHandle for
                                 sourceTableHandle in sourceTableHandles}

        visitSummaryHandles = handleDict['visitSummary']
        visitSummaryHandleDict = {visitSummaryHandle.dataId['visit']: visitSummaryHandle for
                                  visitSummaryHandle in visitSummaryHandles}

        handleDict['sourceTableHandleDict'] = sourceTableHandleDict
        handleDict['visitSummaryHandleDict'] = visitSummaryHandleDict

        # And the outputs
        if self.config.fgcmOutputProducts.doZeropointOutput:
            photoCalibRefDict = {photoCalibRef.dataId.byName()['visit']:
                                 photoCalibRef for photoCalibRef in outputRefs.fgcmPhotoCalib}
            handleDict['fgcmPhotoCalibs'] = photoCalibRefDict

        if self.config.fgcmOutputProducts.doAtmosphereOutput:
            atmRefDict = {atmRef.dataId.byName()['visit']: atmRef for
                          atmRef in outputRefs.fgcmTransmissionAtmosphere}
            handleDict['fgcmTransmissionAtmospheres'] = atmRefDict

        if self.config.fgcmBuildStars.doReferenceMatches:
            refConfig = LoadReferenceObjectsConfig()
            refConfig.filterMap = self.config.fgcmBuildStars.fgcmLoadReferenceCatalog.filterMap
            loader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                    for ref in inputRefs.refCat],
                                           refCats=butlerQC.get(inputRefs.refCat),
                                           name=self.config.connections.refCat,
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
                                           name=self.config.connections.refCat,
                                           config=refConfig,
                                           log=self.log)
            self.fgcmOutputProducts.refObjLoader = loader

        struct = self.run(handleDict, tract,
                          buildStarsRefObjLoader=buildStarsRefObjLoader)

        if struct.photoCalibCatalogs is not None:
            self.log.info("Outputting photoCalib catalogs.")
            for visit, expCatalog in struct.photoCalibCatalogs:
                butlerQC.put(expCatalog, photoCalibRefDict[visit])
            self.log.info("Done outputting photoCalib catalogs.")

        if struct.atmospheres is not None:
            self.log.info("Outputting atmosphere transmission files.")
            for visit, atm in struct.atmospheres:
                butlerQC.put(atm, atmRefDict[visit])
            self.log.info("Done outputting atmosphere files.")

        # Turn raw repeatability into simple catalog for persistence
        schema = afwTable.Schema()
        schema.addField('rawRepeatability', type=np.float64,
                        doc="Per-band raw repeatability in FGCM calibration.")
        repeatabilityCat = afwTable.BaseCatalog(schema)
        repeatabilityCat.resize(len(struct.repeatability))
        repeatabilityCat['rawRepeatability'][:] = struct.repeatability

        butlerQC.put(repeatabilityCat, outputRefs.fgcmRepeatability)

        return
