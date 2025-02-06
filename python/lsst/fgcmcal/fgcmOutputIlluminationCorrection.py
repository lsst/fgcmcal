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
"""Make an illumination correction from the fgcmcal outputs.
"""
from datetime import datetime, UTC
import numpy as np
from astropy.time import Time, TimeDelta
import dataclasses

from lsst.afw.image import ExposureF, FilterLabel
from lsst.daf.butler import Timespan
from lsst.pipe.base import PipelineTaskConnections, PipelineTaskConfig, connectionTypes, PipelineTask, Struct
import lsst.pex.config

from .fgcmOutputProducts import FgcmOutputProductsTask
from .utilities import computePixelAreaFieldDetector, computeReferencePixelScale

__all__ = ["FgcmOutputIlluminationCorrectionConfig", "FgcmOutputIlluminationCorrectionTask"]


class FgcmOutputIlluminationCorrectionConnections(
    PipelineTaskConnections,
    dimensions=("instrument", "detector"),
    defaultTemplates={"cycleNumber": "0"},
):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )
    fgcm_visit_catalog = connectionTypes.Input(
        doc="Catalog of visit information for fgcm.",
        name="fgcmVisitCatalog",
        storageClass="Catalog",
        dimensions=("instrument",),
    )
    fgcm_fit_parameters_catalog = connectionTypes.Input(
        doc="Catalog of fgcm fit parameters.",
        name="fgcm_Cycle{cycleNumber}_FitParameters",
        storageClass="Catalog",
        dimensions=("instrument",),
    )
    flat = connectionTypes.PrerequisiteInput(
        doc="Flat fields associated with illumination correction epoch.",
        name="flat",
        storageClass="ExposureF",
        dimensions=["instrument", "detector", "physical_filter"],
        isCalibration=True,
        multiple=True,
    )
    illumination_corrections = connectionTypes.Output(
        doc="Illumination corrections from fgcm fit.",
        name="illuminationCorrection",
        storageClass="ExposureF",
        dimensions=["instrument", "detector", "physical_filter"],
        isCalibration=True,
        multiple=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if str(int(config.connections.cycleNumber)) != config.connections.cycleNumber:
            raise ValueError("cycleNumber must be of integer format")

        if not config.use_flat_metadata:
            del self.flat
        else:
            def lookup_flat(dataset_type, registry, data_id, collections):

                time = Time(config.epoch_time, format=config.epoch_format)

                return [
                    registry.findDataset(
                        dataset_type,
                        data_id,
                        physical_filter=physical_filter,
                        timespan=Timespan(
                            time - TimeDelta(1, format="sec"),
                            time + TimeDelta(1, format="sec"),
                        ),
                        collections=collections,
                    )
                    for physical_filter in config.physical_filters
                ]

            self.flat = dataclasses.replace(self.flat, lookupFunction=lookup_flat)


class FgcmOutputIlluminationCorrectionConfig(
    PipelineTaskConfig,
    pipelineConnections=FgcmOutputIlluminationCorrectionConnections,
):
    """Configuration for FgcmOutputIlluminationCorrectionTask."""

    use_flat_metadata = lsst.pex.config.Field(
        doc="Use flat-field metadata for illumination correction metadata?",
        dtype=bool,
        default=True,
    )
    epoch_time = lsst.pex.config.Field(
        doc="Time string (UTC) that corresponds to a date in the desired epoch.",
        dtype=str,
        default=None,
        optional=False,
    )
    epoch_format = lsst.pex.config.Field(
        doc="Format for time string (e.g. iso, mjd, etc.), used by astropy.time.Time()",
        dtype=str,
        default="iso",
    )
    physical_filters = lsst.pex.config.ListField(
        doc="List of physical filters to produce illumination corrections.",
        dtype=str,
        default=[],
    )
    include_wcs_jacobian = lsst.pex.config.Field(
        doc="Include WCS jacobian in illumination correction?",
        dtype=bool,
        default=True,
    )
    approximate_wcs_jacobian = lsst.pex.config.Field(
        doc="Use a Chebyshev approximation of the WCS jacobian in illumination correction?",
        dtype=bool,
        default=True,
    )

    def validate(self):
        try:
            _ = Time(self.epoch_time, format=self.epoch_format)
        except Exception as e:
            raise ValueError(
                "Could not parse epoch_time/epoch_format ", e)


class FgcmOutputIlluminationCorrectionTask(PipelineTask):
    """Output illumination corrections from fgcm."""
    ConfigClass = FgcmOutputIlluminationCorrectionConfig
    _DefaultName = "fgcmOutputIlluminationCorrection"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):

        inputs = butlerQC.get(inputRefs)

        detector_id = butlerQC.quantum.dataId["detector"]

        filter_label_dict = {ref.dataId["physical_filter"]:
                             FilterLabel(physical=ref.dataId["physical_filter"], band=ref.dataId["band"])
                             for ref in outputRefs.illumination_corrections}

        flat_dict = {}
        if self.config.use_flat_metadata:
            for i, flat in enumerate(inputs["flat"]):
                ref = inputRefs.flat[i]
                flat_dict[ref.dataId["physical_filter"]] = (ref.id, flat)

        retval = self.run(
            camera=inputs["camera"],
            detector_id=detector_id,
            fgcm_fit_parameters_catalog=inputs["fgcm_fit_parameters_catalog"],
            filter_label_dict=filter_label_dict,
            flat_dict=flat_dict,
        )

        # And put the outputs.
        illum_corr_ref_dict = {ref.dataId["physical_filter"]:
                               ref for ref in outputRefs.illumination_corrections}
        for physical_filter in retval.illum_corr_dict:
            if physical_filter in illum_corr_ref_dict:
                self.log.debug(
                    "Serializing illumination correction for detector %d, physical_filter %s",
                    detector_id,
                    physical_filter,
                )
                butlerQC.put(retval.illum_corr_dict[physical_filter], illum_corr_ref_dict[physical_filter])

    def run(self, *, camera, detector_id, fgcm_fit_parameters_catalog, filter_label_dict, flat_dict={}):
        """Run the illumination correction output task.

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.camera`
            The camera with camera geometry
        detector_id : `int`
            The id of the detector.
        fgcm_fit_parameters_catalog : `lsst.afw.SimpleCatalog`
            Catalog of fgcm fit parameters.
        filter_label_dict : `dict` [`str`: `lsst.afw.image.FilterLabel`]
            Dictionary of filter labels, keyed by physical_filter.
        flat_dict : `dict` [`str`: (`uuid.UUID`, `lsst.afw.image.ExposureF`]
            Dictionary of UUIDs and flats, keyed by physical_filter.

        Returns
        -------
        struct : `lsst.pipe.base.Struct`
            Output structure with keys:

                ``illum_corr_dict``: dictionary keyed by physical_filter,
                                     with illumination correction ExposureF.
        """
        epoch_time = Time(self.config.epoch_time, format=self.config.epoch_format)
        epoch_mjd = epoch_time.mjd

        detector_index = detector_id - camera[0].getId()

        # This is the illumination correction array from fgcm.
        fgcm_illum_corr = np.zeros(fgcm_fit_parameters_catalog["superstarSize"][0, :], dtype="f8")
        fgcm_illum_corr[:, :, :, :] = fgcm_fit_parameters_catalog["superstar"][0, :].reshape(
            fgcm_illum_corr.shape,
        )

        # These are the filter names associated with the illumination
        # corrections.
        fgcm_filter_names = np.asarray(fgcm_fit_parameters_catalog[0]["lutFilterNames"].split(","))

        epoch_mjd_start = fgcm_fit_parameters_catalog[0]["epochMjdStart"]
        epoch_mjd_end = fgcm_fit_parameters_catalog[0]["epochMjdEnd"]

        epoch_index, = np.where((epoch_mjd > epoch_mjd_start) & (epoch_mjd < epoch_mjd_end))

        if len(epoch_index) == 0:
            raise RuntimeError(f"Could not find epoch at {epoch_time} in fgcm epochs.")

        detector = camera[detector_id]
        xymax = np.array([detector.getBBox().getMaxX(), detector.getBBox().getMaxY()])
        area_scaling = 1. / computeReferencePixelScale(camera)**2.

        illum_corr_dict = {}

        count = 0
        for physical_filter in self.config.physical_filters:
            if physical_filter not in filter_label_dict:
                # There are no data associated, so we do not make an
                # illumination correction.
                continue

            if physical_filter not in fgcm_filter_names:
                self.log.warning(
                    "FgcmOutputIlluminationCorrectionTask configured to generate correction for "
                    f"physical filter {physical_filter} but that filter was not calibrated.",
                )
                continue

            filter_index, = np.where(fgcm_filter_names == physical_filter)

            illum_corr_pars = fgcm_illum_corr[epoch_index, filter_index, detector_index, :].ravel()

            illum_corr = ExposureF(detector.getBBox())
            illum_corr.image.array[:, :] = 1.0

            # Get the flat uuid and units if available.
            if physical_filter in flat_dict:
                flat_uuid = str(flat_dict[physical_filter][0])
                flat_md = flat_dict[physical_filter][1].metadata
                units = flat_md["LSST ISR UNITS"]
            else:
                # This is used for testdata_jointcal testing only.
                flat_uuid = "Unknown"
                # Assume adu unless otherwise specified.
                units = "adu"

            illum_corr.info.setFilter(filter_label_dict[physical_filter])
            illum_corr.metadata["LSST CALIB UUID FLAT"] = flat_uuid
            illum_corr.metadata["LSST ISR UNITS"] = units

            # Creation date. Calibration team standard is for local time to be
            # available. Also form UTC (not TAI) version for easier comparisons
            # across multiple processing sites.
            now = datetime.now(tz=UTC)
            illum_corr.metadata.set(
                "CALIB_CREATION_DATETIME",
                now.strftime("%Y-%m-%dT%T"),
                comment="UTC of processing",
            )
            local_time = now.astimezone()
            calib_date = local_time.strftime("%Y-%m-%d")
            calib_time = local_time.strftime("%X %Z")
            illum_corr.metadata.set(
                "CALIB_CREATION_DATE",
                calib_date,
                comment="Local time day of creation",
            )
            illum_corr.metadata.set(
                "CALIB_CREATION_TIME",
                calib_time,
                comment="Local time in day of creation",
            )

            # Make sure this is a legal illumination correction; fgcm
            # uses 100.0 as a sentinel for unfit detectors.
            if illum_corr_pars[0] < 100.0:
                illum_corr_field = FgcmOutputProductsTask._getChebyshevBoundedField(
                    illum_corr_pars,
                    xymax,
                )

                # Check if this is the correct operation!
                illum_corr_field.multiplyImage(illum_corr.image)
                # fgcm includes an additional clipping for strongly vignetted regions.
                illum_corr.image.array[:, :] = np.clip(illum_corr.image.array[:, :], 0.1, None)

            else:
                self.log.warning(
                    f"Invalid illumination correction found for detector {physical_filter} {detector_id}; "
                    "setting to all 1.0s.",
                )

            if self.config.include_wcs_jacobian:
                # This code is here (per-filter) because in the future
                # we plan to allow for a different field per band.
                pixel_area_field = computePixelAreaFieldDetector(
                    camera[detector_id],
                    areaScaling=area_scaling,
                    approximate=self.config.approximate_wcs_jacobian,
                )

                # Check if this is the correct operation!
                pixel_area_field.multiplyImage(illum_corr.image)

            count += 1
            illum_corr_dict[physical_filter] = illum_corr

        self.log.info("Successfully created %d illumination corrections for detector %d", count, detector_id)

        return Struct(illum_corr_dict=illum_corr_dict)
