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
"""Build star observations for input to FGCM using sourceTable_visit.

This task finds all the visits and sourceTable_visits in a repository (or a
subset based on command line parameters) and extracts all the potential
calibration stars for input into fgcm.  This task additionally uses fgcm to
match star observations into unique stars, and performs as much cleaning of the
input catalog as possible.
"""
import warnings
import numpy as np
import esutil
import hpgeom as hpg
from smatch.matcher import Matcher
from astropy.table import Table, vstack

from fgcm.fgcmUtilities import objFlagDict, histogram_rev_sorted

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
from lsst.meas.algorithms import ReferenceObjectLoader, LoadReferenceObjectsConfig
from lsst.pipe.tasks.reserveIsolatedStars import ReserveIsolatedStarsTask

from .fgcmBuildStarsBase import FgcmBuildStarsConfigBase, FgcmBuildStarsBaseTask
from .utilities import computeApproxPixelAreaFields, computeApertureRadiusFromName

__all__ = ["FgcmBuildFromIsolatedStarsConfig", "FgcmBuildFromIsolatedStarsTask"]


class FgcmBuildFromIsolatedStarsConnections(pipeBase.PipelineTaskConnections,
                                            dimensions=("instrument",),
                                            defaultTemplates={}):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )
    fgcm_lookup_table = connectionTypes.PrerequisiteInput(
        doc=("Atmosphere + instrument look-up-table for FGCM throughput and "
             "chromatic corrections."),
        name="fgcmLookUpTable",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )
    ref_cat = connectionTypes.PrerequisiteInput(
        doc="Reference catalog to use for photometric calibration.",
        name="cal_ref_cat",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )
    isolated_star_cats = pipeBase.connectionTypes.Input(
        doc=("Catalog of isolated stars with average positions, number of associated "
             "sources, and indexes to the isolated_star_sources catalogs."),
        name="isolated_star_cat",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "tract", "skymap"),
        deferLoad=True,
        multiple=True,
    )
    isolated_star_sources = pipeBase.connectionTypes.Input(
        doc=("Catalog of isolated star sources with sourceIds, and indexes to the "
             "isolated_star_cats catalogs."),
        name="isolated_star_sources",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "tract", "skymap"),
        deferLoad=True,
        multiple=True,
    )
    visit_summaries = connectionTypes.Input(
        doc=("Per-visit consolidated exposure metadata.  These catalogs use "
             "detector id for the id and must be sorted for fast lookups of a "
             "detector."),
        name="visitSummary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
        deferLoad=True,
        multiple=True,
    )
    fgcm_visit_catalog = connectionTypes.Output(
        doc="Catalog of visit information for fgcm.",
        name="fgcmVisitCatalog",
        storageClass="Catalog",
        dimensions=("instrument",),
    )
    fgcm_star_observations = connectionTypes.Output(
        doc="Catalog of star observations for fgcm.",
        name="fgcm_star_observations",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
    )
    fgcm_star_ids = connectionTypes.Output(
        doc="Catalog of fgcm calibration star IDs.",
        name="fgcm_star_ids",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
    )
    fgcm_reference_stars = connectionTypes.Output(
        doc="Catalog of fgcm-matched reference stars.",
        name="fgcm_reference_stars",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doReferenceMatches:
            self.prerequisiteInputs.remove("ref_cat")
            self.prerequisiteInputs.remove("fgcm_lookup_table")
            self.outputs.remove("fgcm_reference_stars")

    def getSpatialBoundsConnections(self):
        return ("isolated_star_cats", "visit_summaries")


class FgcmBuildFromIsolatedStarsConfig(FgcmBuildStarsConfigBase, pipeBase.PipelineTaskConfig,
                                       pipelineConnections=FgcmBuildFromIsolatedStarsConnections):
    """Config for FgcmBuildFromIsolatedStarsTask."""
    referenceCCD = pexConfig.Field(
        doc="Reference detector for checking PSF and background.",
        dtype=int,
        default=40,
    )
    reserve_selection = pexConfig.ConfigurableField(
        target=ReserveIsolatedStarsTask,
        doc="Task to select reserved stars.",
    )

    def setDefaults(self):
        super().setDefaults()

        self.reserve_selection.reserve_name = "fgcmcal"
        self.reserve_selection.reserve_fraction = 0.1

        # The names here correspond to the isolated_star_sources.
        self.instFluxField = 'normCompGaussianFlux_instFlux'
        self.localBackgroundFluxField = 'localBackground_instFlux'
        self.apertureInnerInstFluxField = 'apFlux_12_0_instFlux'
        self.apertureOuterInstFluxField = 'apFlux_17_0_instFlux'

        source_selector = self.sourceSelector["science"]
        source_selector.setDefaults()

        source_selector.doFlags = False
        source_selector.doSignalToNoise = True
        source_selector.doUnresolved = False
        source_selector.doIsolated = False
        source_selector.doRequireFiniteRaDec = False

        source_selector.flags.bad = []

        source_selector.signalToNoise.minimum = 11.0
        source_selector.signalToNoise.maximum = 1000.0
        source_selector.signalToNoise.fluxField = self.instFluxField
        source_selector.signalToNoise.errField = self.instFluxField + 'Err'


class FgcmBuildFromIsolatedStarsTask(FgcmBuildStarsBaseTask):
    """Build star catalog for FGCM global calibration, using the isolated star catalogs.
    """
    ConfigClass = FgcmBuildFromIsolatedStarsConfig
    _DefaultName = "fgcmBuildFromIsolatedStars"

    canMultiprocess = False

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask('reserve_selection')

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        input_ref_dict = butlerQC.get(inputRefs)

        isolated_star_cat_handles = input_ref_dict["isolated_star_cats"]
        isolated_star_source_handles = input_ref_dict["isolated_star_sources"]

        isolated_star_cat_handle_dict = {
            handle.dataId["tract"]: handle for handle in isolated_star_cat_handles
        }
        isolated_star_source_handle_dict = {
            handle.dataId["tract"]: handle for handle in isolated_star_source_handles
        }

        if len(isolated_star_cat_handle_dict) != len(isolated_star_source_handle_dict):
            raise RuntimeError("isolated_star_cats and isolate_star_sources must have same length.")

        for tract in isolated_star_cat_handle_dict:
            if tract not in isolated_star_source_handle_dict:
                raise RuntimeError(f"tract {tract} in isolated_star_cats but not isolated_star_sources")

        if self.config.doReferenceMatches:
            lookup_table_handle = input_ref_dict["fgcm_lookup_table"]

            # Prepare the reference catalog loader
            ref_config = LoadReferenceObjectsConfig()
            ref_config.filterMap = self.config.fgcmLoadReferenceCatalog.filterMap
            ref_obj_loader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                            for ref in inputRefs.ref_cat],
                                                   refCats=butlerQC.get(inputRefs.ref_cat),
                                                   name=self.config.connections.ref_cat,
                                                   log=self.log,
                                                   config=ref_config)
            self.makeSubtask('fgcmLoadReferenceCatalog',
                             refObjLoader=ref_obj_loader,
                             refCatName=self.config.connections.ref_cat)
        else:
            lookup_table_handle = None

        # The visit summary handles for use with fgcmMakeVisitCatalog must be keyed with
        # visit, and values are a list with the first value being the visit_summary_handle,
        # the second is the source catalog (which is not used, but a place-holder is
        # left for compatibility reasons.)
        visit_summary_handle_dict = {handle.dataId['visit']: [handle, None] for
                                     handle in input_ref_dict['visit_summaries']}

        camera = input_ref_dict["camera"]

        struct = self.run(
            camera=camera,
            visit_summary_handle_dict=visit_summary_handle_dict,
            isolated_star_cat_handle_dict=isolated_star_cat_handle_dict,
            isolated_star_source_handle_dict=isolated_star_source_handle_dict,
            lookup_table_handle=lookup_table_handle,
        )

        butlerQC.put(struct.fgcm_visit_catalog, outputRefs.fgcm_visit_catalog)
        butlerQC.put(struct.fgcm_star_observations, outputRefs.fgcm_star_observations)
        butlerQC.put(struct.fgcm_star_ids, outputRefs.fgcm_star_ids)
        if self.config.doReferenceMatches:
            butlerQC.put(struct.fgcm_reference_stars, outputRefs.fgcm_reference_stars)

    def run(self, *, camera, visit_summary_handle_dict, isolated_star_cat_handle_dict,
            isolated_star_source_handle_dict, lookup_table_handle=None):
        """Run the fgcmBuildFromIsolatedStarsTask.

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`
            Camera object.
        visit_summary_handle_dict : `dict` [`int`, [`lsst.daf.butler.DeferredDatasetHandle`]]
            Visit summary dataset handles, with the visit as key.
        isolated_star_cat_handle_dict : `dict` [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Isolated star catalog dataset handles, with the tract as key.
        isolated_star_source_handle_dict : `dict` [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Isolated star source dataset handles, with the tract as key.
        lookup_table_handle : `lsst.daf.butler.DeferredDatasetHandle`, optional
            Data reference to fgcm look-up table (used if matching reference stars).

        Returns
        -------
        struct : `lsst.pipe.base.struct`
            Catalogs for persistence, with attributes:

            ``fgcm_visit_catalog``
                Catalog of per-visit data (`lsst.afw.table.ExposureCatalog`).
            ``fgcm_star_observations``
                Catalog of individual star observations (`astropy.table.Table`).
            ``fgcm_star_ids``
                Catalog of per-star ids and positions (`astropy.table.Table`).
            ``fgcm_reference_stars``
                Catalog of reference stars matched to fgcm stars (`astropy.table.Table`).
        """
        # Compute aperture radius if necessary.  This is useful to do now before
        # any heave lifting has happened (fail early).
        calib_flux_aperture_radius = None
        if self.config.doSubtractLocalBackground:
            try:
                calib_flux_aperture_radius = computeApertureRadiusFromName(self.config.instFluxField)
            except RuntimeError as e:
                raise RuntimeError("Could not determine aperture radius from %s. "
                                   "Cannot use doSubtractLocalBackground." %
                                   (self.config.instFluxField)) from e

        # Check that we have the lookup_table_handle if we are doing reference matches.
        if self.config.doReferenceMatches:
            if lookup_table_handle is None:
                raise RuntimeError("Must supply lookup_table_handle if config.doReferenceMatches is True.")

        visit_cat = self.fgcmMakeVisitCatalog(camera, visit_summary_handle_dict)

        # Select and concatenate the isolated stars and sources.
        fgcm_obj, star_obs = self._make_all_star_obs_from_isolated_stars(
            isolated_star_cat_handle_dict,
            isolated_star_source_handle_dict,
            visit_cat,
            camera,
            calib_flux_aperture_radius=calib_flux_aperture_radius,
        )

        self._compute_delta_aper_summary_statistics(
            visit_cat,
            star_obs,
        )

        # Do density down-sampling.
        self._density_downsample(fgcm_obj, star_obs)

        # Mark reserve stars
        self._mark_reserve_stars(fgcm_obj)

        # Do reference association.
        if self.config.doReferenceMatches:
            fgcm_ref = self._associate_reference_stars(lookup_table_handle, visit_cat, fgcm_obj)
        else:
            fgcm_ref = None

        return pipeBase.Struct(
            fgcm_visit_catalog=visit_cat,
            fgcm_star_observations=star_obs,
            fgcm_star_ids=fgcm_obj,
            fgcm_reference_stars=fgcm_ref,
        )

    def _make_all_star_obs_from_isolated_stars(
            self,
            isolated_star_cat_handle_dict,
            isolated_star_source_handle_dict,
            visit_cat,
            camera,
            calib_flux_aperture_radius=None,
    ):
        """Make all star observations from isolated star catalogs.

        Parameters
        ----------
        isolated_star_cat_handle_dict : `dict` [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Isolated star catalog dataset handles, with the tract as key.
        isolated_star_source_handle_dict : `dict` [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Isolated star source dataset handles, with the tract as key.
        visit_cat : `lsst.afw.table.ExposureCatalog`
            Catalog of per-visit data.
        camera : `lsst.afw.cameraGeom.Camera`
            The camera object.
        calib_flux_aperture_radius : `float`, optional
            Radius for the calibration flux aperture.

        Returns
        -------
        fgcm_obj : `astropy.table.Table`
            Catalog of ids and positions for each star.
        star_obs : `astropy.table.Table`
            Catalog of individual star observations.
        """
        source_columns = [
            "sourceId",
            "visit",
            "detector",
            "ra",
            "dec",
            "x",
            "y",
            "physical_filter",
            "band",
            "obj_index",
            self.config.instFluxField,
            self.config.instFluxField + "Err",
            self.config.apertureInnerInstFluxField,
            self.config.apertureInnerInstFluxField + "Err",
            self.config.apertureOuterInstFluxField,
            self.config.apertureOuterInstFluxField + "Err",
        ]

        if self.config.doSubtractLocalBackground:
            source_columns.append(self.config.localBackgroundFluxField)
            local_background_flag_name = self.config.localBackgroundFluxField[0: -len('instFlux')] + 'flag'
            source_columns.append(local_background_flag_name)

        if self.sourceSelector.config.doFlags:
            source_columns.extend(self.sourceSelector.config.flags.bad)

        if self.config.doSubtractLocalBackground:
            local_background_area = np.pi*calib_flux_aperture_radius**2.

        # Compute the approximate pixel area over the full focal plane
        # from the WCS jacobian using the camera model.
        approx_pixel_area_fields = computeApproxPixelAreaFields(camera)

        # Construct mapping from detector number to index.
        detector_mapping = {}
        for detector_index, detector in enumerate(camera):
            detector_mapping[detector.getId()] = detector_index

        star_obs_dtype = [
            ("ra", "f8"),
            ("dec", "f8"),
            ("x", "f8"),
            ("y", "f8"),
            ("visit", "i8"),
            ("detector", "i4"),
            ("inst_mag", "f4"),
            ("inst_mag_err", "f4"),
            ("jacobian", "f4"),
            ("delta_mag_bkg", "f4"),
            ("delta_mag_aper", "f4"),
            ("delta_mag_err_aper", "f4"),
        ]

        fgcm_obj_dtype = [
            ("fgcm_id", "i8"),
            ("isolated_star_id", "i8"),
            ("ra", "f8"),
            ("dec", "f8"),
            ("obs_arr_index", "i4"),
            ("n_obs", "i4"),
            ("obj_flag", "i4"),
        ]

        fgcm_objs = []
        star_obs_cats = []
        merge_source_counter = 0

        k = 2.5/np.log(10.)

        visit_cat_table = visit_cat.asAstropy()

        for tract in sorted(isolated_star_cat_handle_dict):
            stars = isolated_star_cat_handle_dict[tract].get()
            sources = isolated_star_source_handle_dict[tract].get(parameters={"columns": source_columns})

            # Down-select sources.
            good_sources = self.sourceSelector.selectSources(sources).selected
            if self.config.doSubtractLocalBackground:
                good_sources &= (~sources[local_background_flag_name])
                local_background = local_background_area*sources[self.config.localBackgroundFluxField]
                good_sources &= ((sources[self.config.instFluxField] - local_background) > 0)

            if good_sources.sum() == 0:
                self.log.info("No good sources found in tract %d", tract)
                continue

            # Need to count the observations of each star after cuts, per band.
            # If we have requiredBands specified, we must make sure that each star
            # has the minumum number of observations in each of thos bands.
            # Otherwise, we must make sure that each star has at least the minimum
            # number of observations in _any_ band.
            if len(self.config.requiredBands) > 0:
                loop_bands = self.config.requiredBands
            else:
                loop_bands = np.unique(sources["band"])

            n_req = np.zeros((len(loop_bands), len(stars)), dtype=np.int32)
            for i, band in enumerate(loop_bands):
                (band_use,) = (sources[good_sources]["band"] == band).nonzero()
                np.add.at(
                    n_req,
                    (i, sources[good_sources]["obj_index"][band_use]),
                    1,
                )

            if len(self.config.requiredBands) > 0:
                # The min gives us the band with the fewest observations, which must be
                # above the limit.
                (good_stars,) = (n_req.min(axis=0) >= self.config.minPerBand).nonzero()
            else:
                # The max gives us the band with the most observations, which must be
                # above the limit.
                (good_stars,) = (n_req.max(axis=0) >= self.config.minPerBand).nonzero()

            if len(good_stars) == 0:
                self.log.info("No good stars found in tract %d", tract)
                continue

            # With the following matching:
            #   sources[good_sources][b] <-> stars[good_stars[a]]
            obj_index = sources["obj_index"][good_sources]
            a, b = esutil.numpy_util.match(good_stars, obj_index)

            # Update indexes and cut to selected stars/sources
            _, index_new = np.unique(a, return_index=True)
            stars["source_cat_index"][good_stars] = index_new
            sources = sources[good_sources][b]
            sources["obj_index"][:] = a
            stars = stars[good_stars]

            nsource = np.zeros(len(stars), dtype=np.int32)
            np.add.at(
                nsource,
                sources["obj_index"],
                1,
            )
            stars["nsource"][:] = nsource

            # After these cuts, the catalogs have the following properties:
            # - ``stars`` only contains isolated stars that have the minimum number of good
            #   sources in the required bands.
            # - ``sources`` has been cut to the good sources.
            # - The slice [stars["source_cat_index"]: stars["source_cat_index"]
            #                                         + stars["nsource"]]
            #   applied to ``sources`` will give all the sources associated with the star.
            # - For each source, sources["obj_index"] points to the index of the associated
            #   isolated star.

            # We now reformat the sources and compute the ``observations`` that fgcm expects.
            star_obs = Table(data=np.zeros(len(sources), dtype=star_obs_dtype))
            star_obs["ra"] = sources["ra"]
            star_obs["dec"] = sources["dec"]
            star_obs["x"] = sources["x"]
            star_obs["y"] = sources["y"]
            star_obs["visit"] = sources["visit"]
            star_obs["detector"] = sources["detector"]

            visit_match, obs_match = esutil.numpy_util.match(visit_cat_table["visit"], sources["visit"])

            exp_time = np.zeros(len(star_obs))
            exp_time[obs_match] = visit_cat_table["exptime"][visit_match]

            with warnings.catch_warnings():
                # Ignore warnings, we will filter infinities and nans below.
                warnings.simplefilter("ignore")

                inst_mag_inner = -2.5*np.log10(sources[self.config.apertureInnerInstFluxField])
                inst_mag_err_inner = k*(sources[self.config.apertureInnerInstFluxField + "Err"]
                                        / sources[self.config.apertureInnerInstFluxField])
                inst_mag_outer = -2.5*np.log10(sources[self.config.apertureOuterInstFluxField])
                inst_mag_err_outer = k*(sources[self.config.apertureOuterInstFluxField + "Err"]
                                        / sources[self.config.apertureOuterInstFluxField])
                star_obs["delta_mag_aper"] = inst_mag_inner - inst_mag_outer
                star_obs["delta_mag_err_aper"] = np.sqrt(inst_mag_err_inner**2. + inst_mag_err_outer**2.)
                # Set bad values to sentinel values for fgcm.
                bad = ~np.isfinite(star_obs["delta_mag_aper"])
                star_obs["delta_mag_aper"][bad] = 99.0
                star_obs["delta_mag_err_aper"][bad] = 99.0

            if self.config.doSubtractLocalBackground:
                # At the moment we only adjust the flux and not the flux
                # error by the background because the error on
                # base_LocalBackground_instFlux is the rms error in the
                # background annulus, not the error on the mean in the
                # background estimate (which is much smaller, by sqrt(n)
                # pixels used to estimate the background, which we do not
                # have access to in this task).  In the default settings,
                # the annulus is sufficiently large such that these
                # additional errors are negligibly small (much less
                # than a mmag in quadrature).

                # This is the difference between the mag with local background correction
                # and the mag without local background correction.
                local_background = local_background_area*sources[self.config.localBackgroundFluxField]
                star_obs["delta_mag_bkg"] = (-2.5*np.log10(sources[self.config.instFluxField]
                                                           - local_background) -
                                             -2.5*np.log10(sources[self.config.instFluxField]))

            # Need to loop over detectors here.
            for detector in camera:
                detector_id = detector.getId()
                # used index for all observations with a given detector
                (use,) = (star_obs["detector"][obs_match] == detector_id).nonzero()
                # Prior to running the calibration, we want to remove the effect
                # of the jacobian of the WCS because that is a known quantity.
                # Ideally, this would be done for each individual WCS, but this
                # is extremely slow and makes small differences that are much
                # smaller than the variation in the throughput due to other issues.
                # Therefore, we use the approximate jacobian estimated from the
                # camera model.
                jac = approx_pixel_area_fields[detector_id].evaluate(
                    star_obs["x"][obs_match][use],
                    star_obs["y"][obs_match][use],
                )
                star_obs["jacobian"][obs_match[use]] = jac
                scaled_inst_flux = (sources[self.config.instFluxField][obs_match[use]]
                                    * visit_cat_table["scaling"][visit_match[use],
                                                                 detector_mapping[detector_id]])
                star_obs["inst_mag"][obs_match[use]] = (-2.5 * np.log10(scaled_inst_flux
                                                                        / exp_time[use]))

            # Compute instMagErr from inst_flux_err/inst_flux; scaling will cancel out.
            star_obs["inst_mag_err"] = k*(sources[self.config.instFluxField + "Err"]
                                          / sources[self.config.instFluxField])

            # Apply the jacobian if configured to do so.
            if self.config.doApplyWcsJacobian:
                star_obs["inst_mag"] -= 2.5*np.log10(star_obs["jacobian"])

            # We now reformat the stars and compute the ''objects'' that fgcm expects.
            fgcm_obj = Table(data=np.zeros(len(stars), dtype=fgcm_obj_dtype))
            fgcm_obj["isolated_star_id"] = stars["isolated_star_id"]
            fgcm_obj["ra"] = stars["ra"]
            fgcm_obj["dec"] = stars["dec"]
            fgcm_obj["obs_arr_index"] = stars["source_cat_index"]
            fgcm_obj["n_obs"] = stars["nsource"]

            # Offset indexes to account for tract merging
            fgcm_obj["obs_arr_index"] += merge_source_counter

            fgcm_objs.append(fgcm_obj)
            star_obs_cats.append(star_obs)

            merge_source_counter += len(star_obs)

        fgcm_obj = vstack(fgcm_objs)

        # Set the fgcm_id to a unique 64-bit integer for easier sorting.
        fgcm_obj["fgcm_id"][:] = np.arange(len(fgcm_obj)) + 1

        return fgcm_obj, vstack(star_obs_cats)

    def _density_downsample(self, fgcm_obj, star_obs):
        """Downsample stars according to density.

        Catalogs are modified in-place.

        Parameters
        ----------
        fgcm_obj : `astropy.table.Table`
            Catalog of per-star ids and positions.
        star_obs : `astropy.table.Table`
            Catalog of individual star observations.
        """
        if self.config.randomSeed is not None:
            np.random.seed(seed=self.config.randomSeed)

        ipnest = hpg.angle_to_pixel(self.config.densityCutNside, fgcm_obj["ra"], fgcm_obj["dec"])
        # Use the esutil.stat.histogram function to pull out the histogram of
        # grouped pixels, and the rev_indices which describes which inputs
        # are grouped together.  The fgcm histogram_rev_sorted shim
        # ensures that the indices are properly sorted.
        hist, rev_indices = histogram_rev_sorted(ipnest)

        obj_use = np.ones(len(fgcm_obj), dtype=bool)

        (high,) = (hist > self.config.densityCutMaxPerPixel).nonzero()
        (ok,) = (hist > 0).nonzero()
        self.log.info("There are %d/%d pixels with high stellar density.", high.size, ok.size)
        for i in range(high.size):
            # The pix_indices are the indices of every star in the pixel.
            pix_indices = rev_indices[rev_indices[high[i]]: rev_indices[high[i] + 1]]
            # Cut down to the maximum number of stars in the pixel.
            cut = np.random.choice(
                pix_indices,
                size=pix_indices.size - self.config.densityCutMaxPerPixel,
                replace=False,
            )
            obj_use[cut] = False

        fgcm_obj = fgcm_obj[obj_use]

        obs_index = np.zeros(np.sum(fgcm_obj["n_obs"]), dtype=np.int32)
        ctr = 0
        for i in range(len(fgcm_obj)):
            n_obs = fgcm_obj["n_obs"][i]
            obs_index[ctr: ctr + n_obs] = (
                np.arange(fgcm_obj["obs_arr_index"][i], fgcm_obj["obs_arr_index"][i] + n_obs)
            )
            fgcm_obj["obs_arr_index"][i] = ctr
            ctr += n_obs

        star_obs = star_obs[obs_index]

    def _mark_reserve_stars(self, fgcm_obj):
        """Run the star reservation task to flag reserved stars.

        Parameters
        ----------
        fgcm_obj : `astropy.table.Table`
            Catalog of per-star ids and positions.
        """
        reserved = self.reserve_selection.run(
            len(fgcm_obj),
            extra=','.join(self.config.requiredBands),
        )
        fgcm_obj["obj_flag"][reserved] |= objFlagDict["RESERVED"]

    def _associate_reference_stars(self, lookup_table_handle, visit_cat, fgcm_obj):
        """Associate fgcm isolated stars with reference stars.

        Parameters
        ----------
        lookup_table_handle : `lsst.daf.butler.DeferredDatasetHandle`, optional
            Data reference to fgcm look-up table (used if matching reference stars).
        visit_cat : `lsst.afw.table.ExposureCatalog`
            Catalog of per-visit data.
        fgcm_obj : `astropy.table.Table`
            Catalog of per-star ids and positions

        Returns
        -------
        ref_cat : `astropy.table.Table`
            Catalog of reference stars matched to fgcm stars.
        """
        # Figure out the correct bands/filters that we need based on the data.
        lut_cat = lookup_table_handle.get()

        std_filter_dict = {filter_name: std_filter for (filter_name, std_filter) in
                           zip(lut_cat[0]["physicalFilters"].split(","),
                               lut_cat[0]["stdPhysicalFilters"].split(","))}
        std_lambda_dict = {std_filter: std_lambda for (std_filter, std_lambda) in
                           zip(lut_cat[0]["stdPhysicalFilters"].split(","),
                               lut_cat[0]["lambdaStdFilter"])}
        del lut_cat

        reference_filter_names = self._getReferenceFilterNames(
            visit_cat,
            std_filter_dict,
            std_lambda_dict,
        )
        self.log.info("Using the following reference filters: %s", (", ".join(reference_filter_names)))

        # Split into healpix pixels for more efficient matching.
        ipnest = hpg.angle_to_pixel(self.config.coarseNside, fgcm_obj["ra"], fgcm_obj["dec"])
        hpix, revpix = histogram_rev_sorted(ipnest)

        pixel_cats = []

        # Compute the dtype from the filter names.
        dtype = [("fgcm_id", "i4"),
                 ("refMag", "f4", (len(reference_filter_names), )),
                 ("refMagErr", "f4", (len(reference_filter_names), ))]

        (gdpix,) = (hpix > 0).nonzero()
        for ii, gpix in enumerate(gdpix):
            p1a = revpix[revpix[gpix]: revpix[gpix + 1]]

            # We do a simple wrapping of RA if we need to.
            ra_wrap = fgcm_obj["ra"][p1a]
            if (ra_wrap.min() < 10.0) and (ra_wrap.max() > 350.0):
                ra_wrap[ra_wrap > 180.0] -= 360.0
                mean_ra = np.mean(ra_wrap) % 360.0
            else:
                mean_ra = np.mean(ra_wrap)
            mean_dec = np.mean(fgcm_obj["dec"][p1a])

            dist = esutil.coords.sphdist(mean_ra, mean_dec, fgcm_obj["ra"][p1a], fgcm_obj["dec"][p1a])
            rad = dist.max()

            if rad < hpg.nside_to_resolution(self.config.coarseNside)/2.:
                # Small radius, read just the circle.
                ref_cat = self.fgcmLoadReferenceCatalog.getFgcmReferenceStarsSkyCircle(
                    mean_ra,
                    mean_dec,
                    rad,
                    reference_filter_names,
                )
            else:
                # Otherwise the less efficient but full coverage.
                ref_cat = self.fgcmLoadReferenceCatalog.getFgcmReferenceStarsHealpix(
                    self.config.coarseNside,
                    hpg.nest_to_ring(self.config.coarseNside, ipnest[p1a[0]]),
                    reference_filter_names,
                )
            if ref_cat.size == 0:
                # No stars in this pixel; that's okay.
                continue

            with Matcher(fgcm_obj["ra"][p1a], fgcm_obj["dec"][p1a]) as matcher:
                idx, i1, i2, d = matcher.query_radius(
                    ref_cat["ra"],
                    ref_cat["dec"],
                    self.config.matchRadius/3600.,
                    return_indices=True,
                )

            if len(i1) == 0:
                # No matched stars in this pixel; that's okay.
                continue

            pixel_cat = Table(data=np.zeros(i1.size, dtype=dtype))
            pixel_cat["fgcm_id"] = fgcm_obj["fgcm_id"][p1a[i1]]
            pixel_cat["refMag"][:, :] = ref_cat["refMag"][i2, :]
            pixel_cat["refMagErr"][:, :] = ref_cat["refMagErr"][i2, :]

            pixel_cats.append(pixel_cat)

            self.log.info(
                "Found %d reference matches in pixel %d (%d of %d).",
                len(pixel_cat),
                ipnest[p1a[0]],
                ii,
                gdpix.size - 1,
            )

        ref_cat_full = vstack(pixel_cats)
        ref_cat_full.meta.update({'FILTERNAMES': reference_filter_names})

        return ref_cat_full

    def _compute_delta_aper_summary_statistics(self, visit_cat, star_obs):
        """Compute delta aperture summary statistics (per visit).

        Parameters
        ----------
        visit_cat : `lsst.afw.table.ExposureCatalog`
            Catalog of per-visit data.
        star_obs : `astropy.table.Table`
            Catalog of individual star observations.
        """
        (ok,) = ((star_obs["delta_mag_aper"] < 99.0)
                 & (star_obs["delta_mag_err_aper"] < 99.0)).nonzero()

        visit_index = np.zeros(len(star_obs[ok]), dtype=np.int32)
        a, b = esutil.numpy_util.match(visit_cat["visit"], star_obs["visit"][ok])
        visit_index[b] = a

        h, rev = histogram_rev_sorted(visit_index)

        visit_cat["deltaAper"] = -9999.0
        h_use, = np.where(h >= 3)
        for index in h_use:
            i1a = rev[rev[index]: rev[index + 1]]
            visit_cat["deltaAper"][index] = np.median(star_obs["delta_mag_aper"][ok[i1a]])
