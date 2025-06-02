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
"""Perform a single fit cycle of FGCM.

This task runs a single "fit cycle" of fgcm.  Prior to running this task
one must run both fgcmMakeLut (to construct the atmosphere and instrumental
look-up-table) and fgcmBuildStars (to extract visits and star observations
for the global fit).

The fgcmFitCycle is meant to be run multiple times, and is tracked by the
'cycleNumber'.  After each run of the fit cycle, diagnostic plots should
be inspected to set parameters for outlier rejection on the following
cycle.  Please see the fgcmcal Cookbook for details.
"""

import copy

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
import lsst.afw.table as afwTable

from .utilities import makeConfigDict, translateFgcmLut, translateVisitCatalog
from .utilities import extractReferenceMags
from .utilities import makeZptSchema, makeZptCat
from .utilities import makeAtmSchema, makeAtmCat, makeStdSchema, makeStdCat
from .sedterms import SedboundarytermDict, SedtermDict
from .focalPlaneProjector import FocalPlaneProjector

import fgcm

__all__ = ['FgcmFitCycleConfig', 'FgcmFitCycleTask']

MULTIPLE_CYCLES_MAX = 10


class FgcmFitCycleConnections(pipeBase.PipelineTaskConnections,
                              dimensions=("instrument",),
                              defaultTemplates={"previousCycleNumber": "-1",
                                                "cycleNumber": "0"}):
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

    fgcmVisitCatalog = connectionTypes.Input(
        doc="Catalog of visit information for fgcm",
        name="fgcmVisitCatalog",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStarObservationsParquet = connectionTypes.Input(
        doc=("Catalog of star observations for fgcm, in parquet format. "
             "Used if useParquetCatalogFormat is True."),
        name="fgcm_star_observations",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStarIdsParquet = connectionTypes.Input(
        doc=("Catalog of fgcm calibration star IDs, in parquet format. "
             "Used if useParquetCatalogFormat is True."),
        name="fgcm_star_ids",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmReferenceStarsParquet = connectionTypes.Input(
        doc=("Catalog of fgcm-matched reference stars, in parquet format. "
             "Used if useParquetCatalogFormat is True."),
        name="fgcm_reference_stars",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStarObservations = connectionTypes.Input(
        doc=("Catalog of star observations for fgcm; old format. "
             "Used if useParquetCatalogFormat is False."),
        name="fgcmStarObservations",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStarIds = connectionTypes.Input(
        doc=("Catalog of fgcm calibration star IDs. "
             "Used if useParquetCatalogFormat is False."),
        name="fgcmStarIds",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmStarIndices = connectionTypes.Input(
        doc=("Catalog of fgcm calibration star indices; old format."
             "Used if useParquetCatalogFormat is False."),
        name="fgcmStarIndices",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    fgcmReferenceStars = connectionTypes.Input(
        doc=("Catalog of fgcm-matched reference stars; old format."
             "Used if useParquetCatalogFormat is False."),
        name="fgcmReferenceStars",
        storageClass="Catalog",
        dimensions=("instrument",),
        deferLoad=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        # The connection parameters ``cycleNumber`` and ``previousCycleNumber``
        # always comes in as strings, and doing a round-trip via
        # integer-to-string ensures that they actually are integers.
        if str(int(config.connections.cycleNumber)) != config.connections.cycleNumber:
            raise ValueError("cycleNumber must be of integer format")
        if str(int(config.connections.previousCycleNumber)) != config.connections.previousCycleNumber:
            raise ValueError("previousCycleNumber must be of integer format")
        # We additionally confirm that there is a delta of 1 between them.
        # Note that this math cannot be done in the parameters because they are
        # strings.
        if int(config.connections.previousCycleNumber) != (int(config.connections.cycleNumber) - 1):
            raise ValueError("previousCycleNumber must be 1 less than cycleNumber")

        instDims = ("instrument",)
        bandDims = ("instrument", "band")
        filterDims = ("instrument", "physical_filter")

        inputAndOutputConnections = [
            ("FitParameters", "Catalog", "Catalog of fgcm fit parameters.", instDims),
            ("FlaggedStars", "Catalog", "Catalog of flagged stars for fgcm calibration.", instDims),
        ]
        multicycleOutputConnections = [
            ("OutputConfig", "Config", "Configuration for next fgcm fit cycle.", instDims),
        ]
        optionalZpOutputConnections = [
            ("Zeropoints", "Catalog", "Catalog of fgcm zeropoint data.", instDims),
            ("AtmosphereParameters", "Catalog", "Catalog of atmospheric fit parameters.", instDims),
        ]
        optionalStarOutputConnections = [
            ("StandardStars", "SimpleCatalog", "Catalog of standard star magnitudes.", instDims),
        ]

        epochs = [f"epoch{i}" for i in range(len(config.epochMjds))]

        # All names will be preceeded with ``fgcm_CycleN_``.
        plotConnections = [
            ("Zeropoints_Plot", "Plot", "Plot of fgcm zeropoints.", instDims),
            ("ExpgrayDeep_Plot",
             "Plot",
             "Plot of gray term per exposure per time for deep fields.",
             instDims),
            ("NightlyAlpha_Plot", "Plot", "Plot of nightly AOD alpha term.", instDims),
            ("NightlyTau_Plot", "Plot", "Plot of nightly aerosol optical depth (tau).", instDims),
            ("NightlyPwv_Plot", "Plot", "Plot of nightly water vapor.", instDims),
            ("NightlyO3_Plot", "Plot", "Plot of nightly ozone.", instDims),
            ("FilterOffsets_Plot", "Plot", "Plot of in-band filter offsets.", instDims),
            ("AbsThroughputs_Plot", "Plot", "Plot of absolute throughput fractions.", instDims),
            ("QESysWashesInitial_Plot", "Plot", "Plot of initial system QE with mirror washes.", instDims),
            ("QESysWashesFinal_Plot", "Plot", "Plot of final system QE with mirror washes.", instDims),
            ("RpwvVsRpwvInput_Plot",
             "Plot",
             "Plot of change in per-visit ``retrieved`` PWV from previous fit cycle.",
             instDims),
            ("RpwvVsRpwvSmooth_Plot",
             "Plot",
             "Plot of per-visit ``retrieved`` PWV vs. smoothed PWV.",
             instDims),
            ("ModelPwvVsRpwv_Plot",
             "Plot",
             "Plot of model PWV vs. per-visit ``retrieved`` PWV.",
             instDims),
            ("ChisqFit_Plot",
             "Plot",
             "Plot of chisq as a function of iteration.",
             instDims),
        ]

        plotConnections.extend(
            [
                ("Apercorr_Plot", "Plot", "Plot of fgcm aperture corrections.", bandDims),
                ("EpsilonGlobal_Plot",
                 "Plot",
                 "Plot of global background over/undersubtraction.",
                 bandDims),
                ("EpsilonMap_Plot",
                 "Plot",
                 "Map of spatially varying background over/undersubtraction.",
                 bandDims),
                ("ExpgrayInitial_Plot",
                 "Plot",
                 "Histogram of initial gray term per exposure.",
                 bandDims),
                ("CompareRedblueExpgray_Plot",
                 "Plot",
                 "Plot of red/blue split gray term per exposure",
                 bandDims),
                ("Expgray_Plot", "Plot", "Histogram of gray term per exposure.", bandDims),
                ("ExpgrayAirmass_Plot",
                 "Plot",
                 "Plot of exposure gray term as a function of airmass.",
                 bandDims),
                ("ExpgrayCompareMjdRedblue_Plot",
                 "Plot",
                 "Plot of red/blue split gray term per exposure as a function of time.",
                 bandDims),
                ("ExpgrayUT_Plot",
                 "Plot",
                 "Plot of grey term per exposure as a function of time of night.",
                 bandDims),
                ("ExpgrayCompareBands_Plot",
                 "Plot",
                 "Plot of gray term per exposure between bands nearby in time.",
                 bandDims),
                ("ExpgrayReference_Plot",
                 "Plot",
                 "Histogram of gray term per exposure compared to reference mags.",
                 bandDims),
                ("QESysRefstarsStdInitial_Plot",
                 "Plot",
                 "Plot of reference mag - calibrated (standard) mag vs. time (before fit).",
                 bandDims),
                ("QESysRefstarsStdFinal_Plot",
                 "Plot",
                 "Plot of reference mag - calibrated (standard) mag vs. time (after fit).",
                 bandDims),
                ("QESysRefstarsObsInitial_Plot",
                 "Plot",
                 "Plot of reference mag - observed (instrumental) mag vs. time (before fit).",
                 bandDims),
                ("QESysRefstarsObsFinal_Plot",
                 "Plot",
                 "Plot of reference mag - observed (instrumental) mag vs. time (after fit).",
                 bandDims),
                ("ModelMagerrInitial_Plot",
                 "Plot",
                 "Plots for magnitude error model, initial estimate.",
                 bandDims),
                ("ModelMagerrPostfit_Plot",
                 "Plot",
                 "Plots for magnitude error model, after fitting.",
                 bandDims),
                ("SigmaFgcmAllStars_Plot",
                 "Plot",
                 "Histograms for intrinsic scatter for all bright stars.",
                 bandDims),
                ("SigmaFgcmReservedStars_Plot",
                 "Plot",
                 "Histograms for intrinsic scatter for reserved bright stars.",
                 bandDims),
                ("SigmaFgcmReservedStarsCrunched_Plot",
                 "Plot",
                 "Histograms for intrinsic scatter for reserved bright stars (after gray correction).",
                 bandDims),
                ("SigmaFgcmPullsAllStars_Plot",
                 "Plot",
                 "Histograms for pulls for all bright stars.",
                 bandDims),
                ("SigmaFgcmPullsReservedStars_Plot",
                 "Plot",
                 "Histograms for pulls for reserved bright stars.",
                 bandDims),
                ("SigmaFgcmPullsReservedStarsCrunched_Plot",
                 "Plot",
                 "Histograms for pulls for reserved bright stars (after gray correction).",
                 bandDims),
                ("SigmaCal_Plot",
                 "Plot",
                 "Plot showing scatter as a function of systematic error floor.",
                 bandDims),
                ("SigmaRef_Plot",
                 "Plot",
                 "Histograms of scatter with respect to reference stars.",
                 bandDims),
                ("RefResidVsColorAll_Plot",
                 "Plot",
                 "Plot of reference star residuals vs. color (all stars).",
                 bandDims),
                ("RefResidVsColorCut_Plot",
                 "Plot",
                 "Plot of reference star residuals vs. color (reference star color cuts).",
                 bandDims),
            ]
        )

        plotConnections.extend(
            [
                ("I1R1_Plot", "Plot", "Plot of fgcm R1 vs. I1.", filterDims),
                ("I1_Plot", "Plot", "Focal plane map of fgcm I1.", filterDims),
                ("R1_Plot", "Plot", "Focal plane map of fgcm R1.", filterDims),
                ("R1mI1Matchscale_Plot", "Plot", "Focal plane map of fgcm R1 - I1.", filterDims),
                ("R1mI1_Plot", "Plot", "Focal plane map of fgcm R1 - I1 (rescaled).", filterDims),
                ("R1mI1_vs_mjd_Plot", "Plot", "R1 - I1 residuals vs. mjd.", filterDims),
                ("CompareRedblueMirrorchrom_Plot",
                 "Plot",
                 "Comparison of mirror chromaticity residuals for red/blue stars.",
                 filterDims),
                ("CcdChromaticity_Plot",
                 "Plot",
                 "Focal plane map of fgcm ccd chromaticity.",
                 filterDims),
                ("EpsilonDetector_Plot",
                 "Plot",
                 "Focal plane map of background over/undersubtraction.",
                 filterDims),
                ("EpsilonDetectorMatchscale_Plot",
                 "Plot",
                 "Focal plane map of background over/undersubtraction.",
                 filterDims),
            ]
        )
        for epoch in epochs:
            plotConnections.extend(
                [
                    (
                        f"Superstar_{epoch}_Plot",
                        "Plot",
                        "Plot of illumination Correction.",
                        filterDims,
                    )
                ]
            )

        if config.doMultipleCycles:
            # Multiple cycle run.

            # All but the final cycle get appended here.
            for cycle in range(config.multipleCyclesFinalCycleNumber):
                outputConnections = copy.copy(inputAndOutputConnections)
                outputConnections.extend(multicycleOutputConnections)
                if config.outputZeropointsBeforeFinalCycle:
                    outputConnections.extend(optionalZpOutputConnections)
                if config.outputStandardsBeforeFinalCycle:
                    outputConnections.extend(optionalStarOutputConnections)

                if config.doPlots:
                    # We will make the plots if doPlots is True and either
                    # it is the next-to-last cycle or the
                    # doPlotsBeforeFinalCycles configuration is set which
                    # means we want plots for all cycles.
                    # The final cycle doPlots is configured below.
                    if cycle == (config.multipleCyclesFinalCycleNumber - 1) \
                       or config.doPlotsBeforeFinalCycles:
                        outputConnections.extend(plotConnections)

                for (name, storageClass, doc, dims) in outputConnections:
                    connectionName = f"fgcm_Cycle{cycle}_{name}"
                    storageName = connectionName
                    outConnection = connectionTypes.Output(
                        name=storageName,
                        storageClass=storageClass,
                        doc=doc,
                        dimensions=dims,
                        multiple=(len(dims) > 1),
                    )
                    setattr(self, connectionName, outConnection)

            # The final cycle has everything.
            outputConnections = copy.copy(inputAndOutputConnections)
            outputConnections.extend(multicycleOutputConnections)
            outputConnections.extend(optionalZpOutputConnections)
            outputConnections.extend(optionalStarOutputConnections)
            if config.doPlots:
                outputConnections.extend(plotConnections)
            for (name, storageClass, doc, dims) in outputConnections:
                connectionName = f"fgcm_Cycle{config.multipleCyclesFinalCycleNumber}_{name}"
                storageName = connectionName
                outConnection = connectionTypes.Output(
                    name=storageName,
                    storageClass=storageClass,
                    doc=doc,
                    dimensions=dims,
                    multiple=(len(dims) > 1),
                )
                setattr(self, connectionName, outConnection)
        else:
            # Single cycle run.
            if config.cycleNumber > 0:
                inputConnections = copy.copy(inputAndOutputConnections)
            else:
                inputConnections = []
            outputConnections = copy.copy(inputAndOutputConnections)
            # The following configurations are also useful for runs
            # where fit cycles are run one-at-a-time since it is
            # not typical to look at these outputs for intermediate
            # steps.
            if config.isFinalCycle or config.outputZeropointsBeforeFinalCycle:
                outputConnections.extend(optionalZpOutputConnections)
            if config.isFinalCycle or config.outputStandardsBeforeFinalCycle:
                outputConnections.extend(optionalStarOutputConnections)

            if config.doPlots:
                outputConnections.extend(plotConnections)

            for (name, storageClass, doc, dims) in inputConnections:
                connectionName = f"fgcm{name}Input"
                storageName = f"fgcm_Cycle{config.cycleNumber - 1}_{name}"
                inConnection = connectionTypes.PrerequisiteInput(
                    name=storageName,
                    storageClass=storageClass,
                    doc=doc,
                    dimensions=dims,
                )
                setattr(self, connectionName, inConnection)

            for (name, storageClass, doc, dims) in outputConnections:
                connectionName = f"fgcm{name}"
                storageName = f"fgcm_Cycle{config.cycleNumber}_{name}"
                # Plots have unique names as well.
                if storageClass == "Plot":
                    connectionName = storageName
                outConnection = connectionTypes.Output(
                    name=storageName,
                    storageClass=storageClass,
                    doc=doc,
                    dimensions=dims,
                    multiple=(len(dims) > 1),
                )
                setattr(self, connectionName, outConnection)

        if not config.doReferenceCalibration:
            self.inputs.remove("fgcmReferenceStars")
            self.inputs.remove("fgcmReferenceStarsParquet")

        if config.useParquetCatalogFormat:
            self.inputs.remove("fgcmStarObservations")
            self.inputs.remove("fgcmStarIds")
            self.inputs.remove("fgcmStarIndices")
            if config.doReferenceCalibration:
                self.inputs.remove("fgcmReferenceStars")
        else:
            self.inputs.remove("fgcmStarObservationsParquet")
            self.inputs.remove("fgcmStarIdsParquet")
            if config.doReferenceCalibration:
                self.inputs.remove("fgcmReferenceStarsParquet")


class FgcmFitCycleConfig(pipeBase.PipelineTaskConfig,
                         pipelineConnections=FgcmFitCycleConnections):
    """Config for FgcmFitCycle"""

    doMultipleCycles = pexConfig.Field(
        doc="Run multiple fit cycles in one task",
        dtype=bool,
        default=False,
    )
    useParquetCatalogFormat = pexConfig.Field(
        doc="Use parquet catalog format?",
        dtype=bool,
        default=True,
    )
    multipleCyclesFinalCycleNumber = pexConfig.RangeField(
        doc=("Final cycle number in multiple cycle mode.  The initial cycle "
             "is 0, with limited parameters fit.  The next cycle is 1 with "
             "full parameter fit.  The final cycle is a clean-up with no "
             "parameters fit.  There will be a total of "
             "(multipleCycleFinalCycleNumber + 1) cycles run, and the final "
             "cycle number cannot be less than 2."),
        dtype=int,
        default=5,
        min=2,
        max=MULTIPLE_CYCLES_MAX,
        inclusiveMax=True,
    )
    bands = pexConfig.ListField(
        doc="Bands to run calibration",
        dtype=str,
        default=[],
    )
    fitBands = pexConfig.ListField(
        doc=("Bands to use in atmospheric fit. The bands not listed here will have "
             "the atmosphere constrained from the 'fitBands' on the same night. "
             "Must be a subset of `config.bands`"),
        dtype=str,
        default=[],
    )
    requiredBands = pexConfig.ListField(
        doc=("Bands that are required for a star to be considered a calibration star. "
             "Must be a subset of `config.bands`"),
        dtype=str,
        default=[],
    )
    physicalFilterMap = pexConfig.DictField(
        doc="Mapping from 'physicalFilter' to band.",
        keytype=str,
        itemtype=str,
        default={},
    )
    doReferenceCalibration = pexConfig.Field(
        doc="Use reference catalog as additional constraint on calibration",
        dtype=bool,
        default=True,
    )
    refStarSnMin = pexConfig.Field(
        doc="Reference star signal-to-noise minimum to use in calibration.  Set to <=0 for no cut.",
        dtype=float,
        default=50.0,
    )
    refStarOutlierNSig = pexConfig.Field(
        doc=("Number of sigma compared to average mag for reference star to be considered an outlier. "
             "Computed per-band, and if it is an outlier in any band it is rejected from fits."),
        dtype=float,
        default=4.0,
    )
    applyRefStarColorCuts = pexConfig.Field(
        doc=("Apply color cuts defined in ``starColorCuts`` to reference stars? "
             "These cuts are in addition to any cuts defined in ``refStarColorCuts``"),
        dtype=bool,
        default=True,
    )
    refStarMaxFracUse = pexConfig.Field(
        doc=("Maximum fraction of reference stars to use in the fit. Remainder will "
             "be used only for validation."),
        dtype=float,
        default=0.5,
    )
    useExposureReferenceOffset = pexConfig.Field(
        doc=("Use per-exposure (visit) offsets between calibrated stars and reference stars "
             "for final zeropoints? This may help uniformity for disjoint surveys."),
        dtype=bool,
        default=False,
    )
    nCore = pexConfig.Field(
        doc="Number of cores to use",
        dtype=int,
        default=4,
        deprecated="Number of cores is deprecated as a config, and will be removed after v27. "
                   "Please use ``pipetask run --cores-per-quantum`` instead.",
    )
    nStarPerRun = pexConfig.Field(
        doc="Number of stars to run in each chunk",
        dtype=int,
        default=200000,
    )
    nExpPerRun = pexConfig.Field(
        doc="Number of exposures to run in each chunk",
        dtype=int,
        default=1000,
    )
    reserveFraction = pexConfig.Field(
        doc="Fraction of stars to reserve for testing",
        dtype=float,
        default=0.1,
    )
    freezeStdAtmosphere = pexConfig.Field(
        doc="Freeze atmosphere parameters to standard (for testing)",
        dtype=bool,
        default=False,
    )
    precomputeSuperStarInitialCycle = pexConfig.Field(
        doc="Precompute superstar flat for initial cycle",
        dtype=bool,
        default=False,
    )
    superStarSubCcdDict = pexConfig.DictField(
        doc=("Per-band specification on whether to compute superstar flat on sub-ccd scale. "
             "Must have one entry per band."),
        keytype=str,
        itemtype=bool,
        default={},
    )
    superStarSubCcdChebyshevOrder = pexConfig.Field(
        doc=("Order of the 2D chebyshev polynomials for sub-ccd superstar fit. "
             "Global default is first-order polynomials, and should be overridden "
             "on a camera-by-camera basis depending on the ISR."),
        dtype=int,
        default=1,
    )
    superStarSubCcdTriangular = pexConfig.Field(
        doc=("Should the sub-ccd superstar chebyshev matrix be triangular to "
             "suppress high-order cross terms?"),
        dtype=bool,
        default=False,
    )
    superStarSigmaClip = pexConfig.Field(
        doc="Number of sigma to clip outliers when selecting for superstar flats",
        dtype=float,
        default=5.0,
    )
    superStarPlotCcdResiduals = pexConfig.Field(
        doc="If plotting is enabled, should per-detector residuals be plotted? "
            "This may produce a lot of output, and should be used only for "
            "debugging purposes.",
        dtype=bool,
        default=False,
    )
    superStarForceZeroMean = pexConfig.Field(
        doc="When computing the super-star flat, force the focal-plane mean to "
            "zero (per band)? This should only be used when computing stand-alone "
            "illumination corrections.",
        dtype=bool,
        default=False,
    )
    focalPlaneSigmaClip = pexConfig.Field(
        doc="Number of sigma to clip outliers per focal-plane.",
        dtype=float,
        default=4.0,
    )
    ccdGraySubCcdDict = pexConfig.DictField(
        doc=("Per-band specification on whether to compute achromatic per-ccd residual "
             "('ccd gray') on a sub-ccd scale."),
        keytype=str,
        itemtype=bool,
        default={},
    )
    ccdGraySubCcdChebyshevOrder = pexConfig.Field(
        doc="Order of the 2D chebyshev polynomials for sub-ccd gray fit.",
        dtype=int,
        default=1,
    )
    ccdGraySubCcdTriangular = pexConfig.Field(
        doc=("Should the sub-ccd gray chebyshev matrix be triangular to "
             "suppress high-order cross terms?"),
        dtype=bool,
        default=True,
    )
    ccdGrayFocalPlaneDict = pexConfig.DictField(
        doc=("Per-band specification on whether to compute focal-plane residual "
             "('ccd gray') corrections."),
        keytype=str,
        itemtype=bool,
        default={},
    )
    ccdGrayFocalPlaneFitMinCcd = pexConfig.Field(
        doc=("Minimum number of 'good' CCDs required to perform focal-plane "
             "gray corrections.  If there are fewer good CCDs then the gray "
             "correction is computed per-ccd."),
        dtype=int,
        default=1,
    )
    ccdGrayFocalPlaneChebyshevOrder = pexConfig.Field(
        doc="Order of the 2D chebyshev polynomials for focal plane fit.",
        dtype=int,
        default=3,
    )
    ccdGrayFocalPlaneMaxStars = pexConfig.Field(
        doc="Maximum number of stars to use for focal plane fit. Required to keep "
            "matrix memory usage from running away. If there are more stars than "
            "this will be down-sampled.",
        dtype=int,
        default=50_000,
    )
    cycleNumber = pexConfig.Field(
        doc=("FGCM fit cycle number.  This is automatically incremented after each run "
             "and stage of outlier rejection.  See cookbook for details."),
        dtype=int,
        default=None,
    )
    isFinalCycle = pexConfig.Field(
        doc=("Is this the final cycle of the fitting?  Will automatically compute final "
             "selection of stars and photometric exposures, and will output zeropoints "
             "and standard stars for use in fgcmOutputProducts"),
        dtype=bool,
        default=False,
    )
    maxIterBeforeFinalCycle = pexConfig.Field(
        doc=("Maximum fit iterations, prior to final cycle.  The number of iterations "
             "will always be 0 in the final cycle for cleanup and final selection."),
        dtype=int,
        default=50,
    )
    deltaMagBkgOffsetPercentile = pexConfig.Field(
        doc=("Percentile brightest stars on a visit/ccd to use to compute net "
             "offset from local background subtraction."),
        dtype=float,
        default=0.25,
    )
    deltaMagBkgPerCcd = pexConfig.Field(
        doc=("Compute net offset from local background subtraction per-ccd? "
             "Otherwise, use computation per visit."),
        dtype=bool,
        default=False,
    )
    utBoundary = pexConfig.Field(
        doc="Boundary (in UTC) from day-to-day",
        dtype=float,
        default=None,
    )
    washMjds = pexConfig.ListField(
        doc="Mirror wash MJDs",
        dtype=float,
        default=(0.0,),
    )
    epochMjds = pexConfig.ListField(
        doc="Epoch boundaries in MJD",
        dtype=float,
        default=(0.0,),
    )
    minObsPerBand = pexConfig.Field(
        doc="Minimum good observations per band",
        dtype=int,
        default=2,
    )
    # TODO: When DM-16511 is done, it will be possible to get the
    # telescope latitude directly from the camera.
    latitude = pexConfig.Field(
        doc="Observatory latitude",
        dtype=float,
        default=None,
    )
    mirrorArea = pexConfig.Field(
        doc="Mirror area (square meters) of telescope.  If not set, will "
            "try to estimate from camera.telescopeDiameter.",
        dtype=float,
        default=None,
        optional=True,
    )
    cameraGain = pexConfig.Field(
        doc="Average camera gain. If not set, will use the median of the "
            "camera model/detector/amplifier gains.",
        dtype=float,
        default=None,
        optional=True,
    )
    defaultCameraOrientation = pexConfig.Field(
        doc="Default camera orientation for QA plots.",
        dtype=float,
        default=None,
    )
    brightObsGrayMax = pexConfig.Field(
        doc="Maximum gray extinction to be considered bright observation",
        dtype=float,
        default=0.15,
    )
    minStarPerCcd = pexConfig.Field(
        doc=("Minimum number of good stars per CCD to be used in calibration fit. "
             "CCDs with fewer stars will have their calibration estimated from other "
             "CCDs in the same visit, with zeropoint error increased accordingly."),
        dtype=int,
        default=5,
    )
    minCcdPerExp = pexConfig.Field(
        doc=("Minimum number of good CCDs per exposure/visit to be used in calibration fit. "
             "Visits with fewer good CCDs will have CCD zeropoints estimated where possible."),
        dtype=int,
        default=5,
    )
    maxCcdGrayErr = pexConfig.Field(
        doc="Maximum error on CCD gray offset to be considered photometric",
        dtype=float,
        default=0.05,
    )
    minStarPerExp = pexConfig.Field(
        doc=("Minimum number of good stars per exposure/visit to be used in calibration fit. "
             "Visits with fewer good stars will have CCD zeropoints estimated where possible."),
        dtype=int,
        default=600,
    )
    minExpPerNight = pexConfig.Field(
        doc="Minimum number of good exposures/visits to consider a partly photometric night",
        dtype=int,
        default=10,
    )
    expGrayInitialCut = pexConfig.Field(
        doc=("Maximum exposure/visit gray value for initial selection of possible photometric "
             "observations."),
        dtype=float,
        default=-0.25,
    )
    expFwhmCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum exposure FWHM (arcseconds) that will "
             "be considered for the model fit. Exposures with median FWHM larger "
             "than this threshold will get zeropoints based on matching to good "
             "stars."),
        keytype=str,
        itemtype=float,
        default={},
    )
    expGrayPhotometricCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum (negative) achromatic exposure residual "
             "('gray term') for a visit to be considered photometric.  Must have one "
             "entry per band. Broad-band filters should be -0.05."),
        keytype=str,
        itemtype=float,
        default={},
    )
    expGrayHighCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum (positive) achromatic exposure residual "
             "('gray term') for a visit to be considered photometric.  Must have one "
             "entry per band.  Broad-band filters should be 0.2."),
        keytype=str,
        itemtype=float,
        default={},
    )
    expGrayRecoverCut = pexConfig.Field(
        doc=("Maximum (negative) exposure gray to be able to recover bad ccds via interpolation. "
             "Visits with more gray extinction will only get CCD zeropoints if there are "
             "sufficient star observations (minStarPerCcd) on that CCD."),
        dtype=float,
        default=-1.0,
    )
    expVarGrayPhotometricCutDict = pexConfig.DictField(
        doc=("Per-band specification on maximum exposure variance to be considered possibly "
             "photometric.  Must have one entry per band.  Broad-band filters should be "
             "0.0005."),
        keytype=str,
        itemtype=float,
        default={},
    )
    expGrayErrRecoverCut = pexConfig.Field(
        doc=("Maximum exposure gray error to be able to recover bad ccds via interpolation. "
             "Visits with more gray variance will only get CCD zeropoints if there are "
             "sufficient star observations (minStarPerCcd) on that CCD."),
        dtype=float,
        default=0.05,
    )
    aperCorrPerCcd = pexConfig.Field(
        doc="Use aperture corrections per-ccd (detector) instead of per-visit?",
        dtype=bool,
        default=False,
    )
    aperCorrFitNBins = pexConfig.Field(
        doc=("Number of aperture bins used in aperture correction fit.  When set to 0"
             "no fit will be performed, and the config.aperCorrInputSlopes will be "
             "used if available."),
        dtype=int,
        default=10,
    )
    aperCorrInputSlopeDict = pexConfig.DictField(
        doc=("Per-band specification of aperture correction input slope parameters.  These "
             "are used on the first fit iteration, and aperture correction parameters will "
             "be updated from the data if config.aperCorrFitNBins > 0.  It is recommended "
             "to set this when there is insufficient data to fit the parameters (e.g. "
             "tract mode)."),
        keytype=str,
        itemtype=float,
        default={},
    )
    sedboundaryterms = pexConfig.ConfigField(
        doc="Mapping from bands to SED boundary term names used is sedterms.",
        dtype=SedboundarytermDict,
    )
    sedterms = pexConfig.ConfigField(
        doc="Mapping from terms to bands for fgcm linear SED approximations.",
        dtype=SedtermDict,
    )
    sigFgcmMaxErr = pexConfig.Field(
        doc="Maximum mag error for fitting sigma_FGCM",
        dtype=float,
        default=0.01,
    )
    sigFgcmMaxEGrayDict = pexConfig.DictField(
        doc=("Per-band specification for maximum (absolute) achromatic residual (gray value) "
             "for observations in sigma_fgcm (raw repeatability).  Broad-band filters "
             "should be 0.05."),
        keytype=str,
        itemtype=float,
        default={},
    )
    ccdGrayMaxStarErr = pexConfig.Field(
        doc=("Maximum error on a star observation to use in ccd gray (achromatic residual) "
             "computation"),
        dtype=float,
        default=0.10,
    )
    approxThroughputDict = pexConfig.DictField(
        doc=("Per-band specification of the approximate overall throughput at the start of "
             "calibration observations.  Must have one entry per band.  Typically should "
             "be 1.0."),
        keytype=str,
        itemtype=float,
        default={},
    )
    sigmaCalRange = pexConfig.ListField(
        doc="Allowed range for systematic error floor estimation",
        dtype=float,
        default=(0.001, 0.003),
    )
    sigmaCalFitPercentile = pexConfig.ListField(
        doc="Magnitude percentile range to fit systematic error floor",
        dtype=float,
        default=(0.05, 0.15),
    )
    sigmaCalPlotPercentile = pexConfig.ListField(
        doc="Magnitude percentile range to plot systematic error floor",
        dtype=float,
        default=(0.05, 0.95),
    )
    sigma0Phot = pexConfig.Field(
        doc="Systematic error floor for all zeropoints",
        dtype=float,
        default=0.003,
    )
    mapLongitudeRef = pexConfig.Field(
        doc="Reference longitude for plotting maps",
        dtype=float,
        default=0.0,
    )
    mapNSide = pexConfig.Field(
        doc="Healpix nside for plotting maps",
        dtype=int,
        default=256,
    )
    outfileBase = pexConfig.Field(
        doc="Filename start for plot output files",
        dtype=str,
        default=None,
    )
    starColorCuts = pexConfig.ListField(
        doc=("Encoded star-color cuts (using calibration star colors). "
             "This is a list with each entry a string of the format "
             "``band1,band2,low,high`` such that only stars of color "
             "low < band1 - band2 < high will be used for calibration."),
        dtype=str,
        default=("NO_DATA",),
    )
    refStarColorCuts = pexConfig.ListField(
        doc=("Encoded star color cuts specifically to apply to reference stars. "
             "This is a list with each entry a string of the format "
             "``band1,band2,low,high`` such that only stars of color "
             "low < band1 - band2 < high will be used as reference stars."),
        dtype=str,
        default=("NO_DATA",),
    )
    colorSplitBands = pexConfig.ListField(
        doc="Band names to use to split stars by color.  Must have 2 entries.",
        dtype=str,
        length=2,
        default=('g', 'i'),
    )
    modelMagErrors = pexConfig.Field(
        doc="Should FGCM model the magnitude errors from sky/fwhm? (False means trust inputs)",
        dtype=bool,
        default=True,
    )
    useQuadraticPwv = pexConfig.Field(
        doc="Model PWV with a quadratic term for variation through the night?",
        dtype=bool,
        default=False,
    )
    instrumentParsPerBand = pexConfig.Field(
        doc=("Model instrumental parameters per band? "
             "Otherwise, instrumental parameters (QE changes with time) are "
             "shared among all bands."),
        dtype=bool,
        default=False,
    )
    instrumentSlopeMinDeltaT = pexConfig.Field(
        doc=("Minimum time change (in days) between observations to use in constraining "
             "instrument slope."),
        dtype=float,
        default=20.0,
    )
    fitMirrorChromaticity = pexConfig.Field(
        doc="Fit (intraband) mirror chromatic term?",
        dtype=bool,
        default=False,
    )
    fitCcdChromaticityDict = pexConfig.DictField(
        doc="Specification on whether to compute first-order quantum efficiency (QE) "
            "adjustments. Key is band, and value will be True or False. Any band "
            "not explicitly specified will default to False.",
        keytype=str,
        itemtype=bool,
        default={},
    )
    coatingMjds = pexConfig.ListField(
        doc="Mirror coating dates in MJD",
        dtype=float,
        default=(0.0,),
    )
    outputStandardsBeforeFinalCycle = pexConfig.Field(
        doc="Output standard stars prior to final cycle?  Used in debugging.",
        dtype=bool,
        default=False,
    )
    outputZeropointsBeforeFinalCycle = pexConfig.Field(
        doc="Output standard stars prior to final cycle?  Used in debugging.",
        dtype=bool,
        default=False,
    )
    useRepeatabilityForExpGrayCutsDict = pexConfig.DictField(
        doc=("Per-band specification on whether to use star repeatability (instead of exposures) "
             "for computing photometric cuts. Recommended for tract mode or bands with few visits."),
        keytype=str,
        itemtype=bool,
        default={},
    )
    autoPhotometricCutNSig = pexConfig.Field(
        doc=("Number of sigma for automatic computation of (low) photometric cut. "
             "Cut is based on exposure gray width (per band), unless "
             "useRepeatabilityForExpGrayCuts is set, in which case the star "
             "repeatability is used (also per band)."),
        dtype=float,
        default=3.0,
    )
    autoHighCutNSig = pexConfig.Field(
        doc=("Number of sigma for automatic computation of (high) outlier cut. "
             "Cut is based on exposure gray width (per band), unless "
             "useRepeatabilityForExpGrayCuts is set, in which case the star "
             "repeatability is used (also per band)."),
        dtype=float,
        default=4.0,
    )
    quietMode = pexConfig.Field(
        doc="Be less verbose with logging.",
        dtype=bool,
        default=False,
    )
    doPlots = pexConfig.Field(
        doc="Make fgcm QA plots.",
        dtype=bool,
        default=True,
    )
    doPlotsBeforeFinalCycles = pexConfig.Field(
        doc="Make fgcm QA plots before the final two fit cycles? This only applies in"
            "multi-cycle mode, and if doPlots is True. These are typically the most"
            "important QA plots to inspect.",
        dtype=bool,
        default=False,
    )
    randomSeed = pexConfig.Field(
        doc="Random seed for fgcm for consistency in tests.",
        dtype=int,
        default=None,
        optional=True,
    )
    deltaAperFitMinNgoodObs = pexConfig.Field(
        doc="Minimum number of good observations to use mean delta-aper values in fits.",
        dtype=int,
        default=2,
    )
    deltaAperFitPerCcdNx = pexConfig.Field(
        doc=("Number of x bins per ccd when computing delta-aper background offsets. "
             "Only used when ``doComputeDeltaAperPerCcd`` is True."),
        dtype=int,
        default=10,
    )
    deltaAperFitPerCcdNy = pexConfig.Field(
        doc=("Number of y bins per ccd when computing delta-aper background offsets. "
             "Only used when ``doComputeDeltaAperPerCcd`` is True."),
        dtype=int,
        default=10,
    )
    deltaAperFitSpatialNside = pexConfig.Field(
        doc="Healpix nside to compute spatial delta-aper background offset maps.",
        dtype=int,
        default=64,
    )
    deltaAperInnerRadiusArcsec = pexConfig.Field(
        doc=("Inner radius used to compute deltaMagAper (arcseconds). "
             "Must be positive and less than ``deltaAperOuterRadiusArcsec`` if "
             "any of ``doComputeDeltaAperPerVisit``, ``doComputeDeltaAperPerStar``, "
             "``doComputeDeltaAperMap``, ``doComputeDeltaAperPerCcd`` are set."),
        dtype=float,
        default=0.0,
    )
    deltaAperOuterRadiusArcsec = pexConfig.Field(
        doc=("Outer radius used to compute deltaMagAper (arcseconds). "
             "Must be positive and greater than ``deltaAperInnerRadiusArcsec`` if "
             "any of ``doComputeDeltaAperPerVisit``, ``doComputeDeltaAperPerStar``, "
             "``doComputeDeltaAperMap``, ``doComputeDeltaAperPerCcd`` are set."),
        dtype=float,
        default=0.0,
    )
    doComputeDeltaAperPerVisit = pexConfig.Field(
        doc=("Do the computation of delta-aper background offsets per visit? "
             "Note: this option can be very slow when there are many visits."),
        dtype=bool,
        default=False,
    )
    doComputeDeltaAperPerStar = pexConfig.Field(
        doc="Do the computation of delta-aper mean values per star?",
        dtype=bool,
        default=True,
    )
    doComputeDeltaAperMap = pexConfig.Field(
        doc=("Do the computation of delta-aper spatial maps? "
             "This is only used if ``doComputeDeltaAperPerStar`` is True,"),
        dtype=bool,
        default=False,
    )
    doComputeDeltaAperPerCcd = pexConfig.Field(
        doc="Do the computation of per-ccd delta-aper background offsets?",
        dtype=bool,
        default=False,
    )

    def validate(self):
        super().validate()

        if self.connections.previousCycleNumber != str(self.cycleNumber - 1):
            msg = "cycleNumber in template must be connections.previousCycleNumber + 1"
            raise RuntimeError(msg)
        if self.connections.cycleNumber != str(self.cycleNumber):
            msg = "cycleNumber in template must be equal to connections.cycleNumber"
            raise RuntimeError(msg)

        for band in self.fitBands:
            if band not in self.bands:
                msg = 'fitBand %s not in bands' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.fitBands, self, msg)
        for band in self.requiredBands:
            if band not in self.bands:
                msg = 'requiredBand %s not in bands' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.requiredBands, self, msg)
        for band in self.colorSplitBands:
            if band not in self.bands:
                msg = 'colorSplitBand %s not in bands' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.colorSplitBands, self, msg)
        for band in self.bands:
            if band not in self.superStarSubCcdDict:
                msg = 'band %s not in superStarSubCcdDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.superStarSubCcdDict,
                                                     self, msg)
            if band not in self.ccdGraySubCcdDict:
                msg = 'band %s not in ccdGraySubCcdDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.ccdGraySubCcdDict,
                                                     self, msg)
            if band not in self.expGrayPhotometricCutDict:
                msg = 'band %s not in expGrayPhotometricCutDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.expGrayPhotometricCutDict,
                                                     self, msg)
            if band not in self.expGrayHighCutDict:
                msg = 'band %s not in expGrayHighCutDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.expGrayHighCutDict,
                                                     self, msg)
            if band not in self.expVarGrayPhotometricCutDict:
                msg = 'band %s not in expVarGrayPhotometricCutDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.expVarGrayPhotometricCutDict,
                                                     self, msg)
            if band not in self.sigFgcmMaxEGrayDict:
                msg = 'band %s not in sigFgcmMaxEGrayDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.sigFgcmMaxEGrayDict,
                                                     self, msg)
            if band not in self.approxThroughputDict:
                msg = 'band %s not in approxThroughputDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.approxThroughputDict,
                                                     self, msg)
            if band not in self.useRepeatabilityForExpGrayCutsDict:
                msg = 'band %s not in useRepeatabilityForExpGrayCutsDict' % (band)
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.useRepeatabilityForExpGrayCutsDict,
                                                     self, msg)

        if self.doComputeDeltaAperPerVisit or self.doComputeDeltaAperMap \
           or self.doComputeDeltaAperPerCcd:
            if self.deltaAperInnerRadiusArcsec <= 0.0:
                msg = 'deltaAperInnerRadiusArcsec must be positive if deltaAper computations are turned on.'
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.deltaAperInnerRadiusArcsec,
                                                     self, msg)
            if self.deltaAperOuterRadiusArcsec <= 0.0:
                msg = 'deltaAperOuterRadiusArcsec must be positive if deltaAper computations are turned on.'
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.deltaAperOuterRadiusArcsec,
                                                     self, msg)
            if self.deltaAperOuterRadiusArcsec <= self.deltaAperInnerRadiusArcsec:
                msg = ('deltaAperOuterRadiusArcsec must be greater than deltaAperInnerRadiusArcsec if '
                       'deltaAper computations are turned on.')
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.deltaAperOuterRadiusArcsec,
                                                     self, msg)


class FgcmFitCycleTask(pipeBase.PipelineTask):
    """
    Run Single fit cycle for FGCM global calibration
    """

    ConfigClass = FgcmFitCycleConfig
    _DefaultName = "fgcmFitCycle"

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        camera = butlerQC.get(inputRefs.camera)

        nCore = butlerQC.resources.num_cores

        handleDict = {}

        handleDict['fgcmLookUpTable'] = butlerQC.get(inputRefs.fgcmLookUpTable)
        handleDict['fgcmVisitCatalog'] = butlerQC.get(inputRefs.fgcmVisitCatalog)

        if self.config.useParquetCatalogFormat:
            handleDict['fgcmStarObservations'] = butlerQC.get(inputRefs.fgcmStarObservationsParquet)
            handleDict['fgcmStarIds'] = butlerQC.get(inputRefs.fgcmStarIdsParquet)
            if self.config.doReferenceCalibration:
                handleDict['fgcmReferenceStars'] = butlerQC.get(inputRefs.fgcmReferenceStarsParquet)
        else:
            handleDict['fgcmStarObservations'] = butlerQC.get(inputRefs.fgcmStarObservations)
            handleDict['fgcmStarIds'] = butlerQC.get(inputRefs.fgcmStarIds)
            handleDict['fgcmStarIndices'] = butlerQC.get(inputRefs.fgcmStarIndices)
            if self.config.doReferenceCalibration:
                handleDict['fgcmReferenceStars'] = butlerQC.get(inputRefs.fgcmReferenceStars)
        if self.config.cycleNumber > 0:
            handleDict['fgcmFlaggedStars'] = butlerQC.get(inputRefs.fgcmFlaggedStarsInput)
            handleDict['fgcmFitParameters'] = butlerQC.get(inputRefs.fgcmFitParametersInput)

        fgcmDatasetDict = None
        if self.config.doMultipleCycles:
            # Run multiple cycles at once.
            config = copy.copy(self.config)
            config.update(cycleNumber=0)
            for cycle in range(self.config.multipleCyclesFinalCycleNumber + 1):
                if cycle == self.config.multipleCyclesFinalCycleNumber:
                    config.update(isFinalCycle=True)

                if cycle > 0:
                    handleDict['fgcmFlaggedStars'] = fgcmDatasetDict['fgcmFlaggedStars']
                    handleDict['fgcmFitParameters'] = fgcmDatasetDict['fgcmFitParameters']

                # Set up plot outputs.
                # Note that nothing will go in the dict if doPlots is False.
                plotHandleDict = {}
                for outputRefName in outputRefs.keys():
                    if outputRefName.endswith("Plot") and f"Cycle{cycle}" in outputRefName:
                        ref = getattr(outputRefs, outputRefName)
                        if isinstance(ref, (tuple, list)):
                            if "physical_filter" in ref[0].dimensions:
                                for filterRef in ref:
                                    handleDictKey = f"{outputRefName}_{filterRef.dataId['physical_filter']}"
                                    plotHandleDict[handleDictKey] = filterRef
                            if "band" in ref[0].dimensions:
                                for bandRef in ref:
                                    handleDictKey = f"{outputRefName}_{bandRef.dataId['band']}"
                                    plotHandleDict[handleDictKey] = bandRef
                        else:
                            plotHandleDict[outputRefName] = ref

                fgcmDatasetDict, config = self._fgcmFitCycle(
                    camera,
                    handleDict,
                    butlerQC=butlerQC,
                    plotHandleDict=plotHandleDict,
                    config=config,
                    nCore=nCore,
                )
                butlerQC.put(fgcmDatasetDict['fgcmFitParameters'],
                             getattr(outputRefs, f'fgcm_Cycle{cycle}_FitParameters'))
                butlerQC.put(fgcmDatasetDict['fgcmFlaggedStars'],
                             getattr(outputRefs, f'fgcm_Cycle{cycle}_FlaggedStars'))
                butlerQC.put(config,
                             getattr(outputRefs, f'fgcm_Cycle{cycle}_OutputConfig'))
                if self.outputZeropoints:
                    butlerQC.put(fgcmDatasetDict['fgcmZeropoints'],
                                 getattr(outputRefs, f'fgcm_Cycle{cycle}_Zeropoints'))
                    butlerQC.put(fgcmDatasetDict['fgcmAtmosphereParameters'],
                                 getattr(outputRefs, f'fgcm_Cycle{cycle}_AtmosphereParameters'))
                if self.outputStandards:
                    butlerQC.put(fgcmDatasetDict['fgcmStandardStars'],
                                 getattr(outputRefs, f'fgcm_Cycle{cycle}_StandardStars'))
        else:
            # Run a single cycle

            # Set up plot outputs.
            # Note that nothing will go in the dict if doPlots is False.
            plotHandleDict = {}
            for outputRefName in outputRefs.keys():
                if outputRefName.endswith("Plot") and f"Cycle{self.config.cycleNumber}" in outputRefName:
                    ref = getattr(outputRefs, outputRefName)
                    if isinstance(ref, (tuple, list)):
                        if "physical_filter" in ref[0].dimensions:
                            for filterRef in ref:
                                handleDictKey = f"{outputRefName}_{filterRef.dataId['physical_filter']}"
                                plotHandleDict[handleDictKey] = filterRef
                        if "band" in ref[0].dimensions:
                            for bandRef in ref:
                                handleDictKey = f"{outputRefName}_{bandRef.dataId['band']}"
                                plotHandleDict[handleDictKey] = bandRef
                    else:
                        plotHandleDict[outputRefName] = ref

            fgcmDatasetDict, _ = self._fgcmFitCycle(
                camera,
                handleDict,
                nCore=nCore,
                butlerQC=butlerQC,
                plotHandleDict=plotHandleDict,
                multiCycle=False,
            )

            butlerQC.put(fgcmDatasetDict['fgcmFitParameters'], outputRefs.fgcmFitParameters)
            butlerQC.put(fgcmDatasetDict['fgcmFlaggedStars'], outputRefs.fgcmFlaggedStars)
            if self.outputZeropoints:
                butlerQC.put(fgcmDatasetDict['fgcmZeropoints'], outputRefs.fgcmZeropoints)
                butlerQC.put(fgcmDatasetDict['fgcmAtmosphereParameters'], outputRefs.fgcmAtmosphereParameters)
            if self.outputStandards:
                butlerQC.put(fgcmDatasetDict['fgcmStandardStars'], outputRefs.fgcmStandardStars)

    def _fgcmFitCycle(
        self,
        camera,
        handleDict,
        butlerQC=None,
        plotHandleDict=None,
        config=None,
        nCore=1,
        multiCycle=True,
    ):
        """
        Run the fit cycle

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`
        handleDict : `dict`
            All handles are `lsst.daf.butler.DeferredDatasetHandle`
            handle dictionary with keys:

            ``"fgcmLookUpTable"``
                handle for the FGCM look-up table.
            ``"fgcmVisitCatalog"``
                handle for visit summary catalog.
            ``"fgcmStarObservations"``
                handle for star observation catalog.
            ``"fgcmStarIds"``
                handle for star id catalog.
            ``"fgcmStarIndices"``
                handle for star index catalog.
            ``"fgcmReferenceStars"``
                handle for matched reference star catalog.
            ``"fgcmFlaggedStars"``
                handle for flagged star catalog.
            ``"fgcmFitParameters"``
                handle for fit parameter catalog.
        butlerQC : `lsst.pipe.base.QuantumContext`, optional
            Quantum context used for serializing plots.
        plotHandleDict : `dict` [`str`, `lsst.daf.butler.DatasetRef`], optional
            Dictionary of plot dataset refs, keyed by plot name.
        config : `lsst.pex.config.Config`, optional
            Configuration to use to override self.config.
        nCore : `int`, optional
            Number of cores to use during fitting.
        multiCycle : `bool`, optional
            Is this part of a multicycle run?

        Returns
        -------
        fgcmDatasetDict : `dict`
            Dictionary of datasets to persist.
        """
        if config is not None:
            _config = config
        else:
            _config = self.config

        # Set defaults on whether to output standards and zeropoints
        self.maxIter = _config.maxIterBeforeFinalCycle
        self.outputStandards = _config.outputStandardsBeforeFinalCycle
        self.outputZeropoints = _config.outputZeropointsBeforeFinalCycle
        self.resetFitParameters = True

        if _config.isFinalCycle:
            # This is the final fit cycle, so we do not want to reset fit
            # parameters, we want to run a final "clean-up" with 0 fit iterations,
            # and we always want to output standards and zeropoints
            self.maxIter = 0
            self.outputStandards = True
            self.outputZeropoints = True
            self.resetFitParameters = False

        lutCat = handleDict['fgcmLookUpTable'].get()
        fgcmLut, lutIndexVals, lutStd = translateFgcmLut(lutCat,
                                                         dict(_config.physicalFilterMap))
        del lutCat

        # Check if we want to do plots.
        doPlots = _config.doPlots
        if doPlots and multiCycle:
            if _config.cycleNumber < (_config.multipleCyclesFinalCycleNumber - 1) \
               and not _config.doPlotsBeforeFinalCycles:
                doPlots = False

        configDict = makeConfigDict(_config, self.log, camera,
                                    self.maxIter, self.resetFitParameters,
                                    self.outputZeropoints,
                                    lutIndexVals[0]['FILTERNAMES'],
                                    nCore=nCore,
                                    doPlots=doPlots)

        # next we need the exposure/visit information
        visitCat = handleDict['fgcmVisitCatalog'].get()
        fgcmExpInfo = translateVisitCatalog(visitCat)
        del visitCat

        if len(camera) == fgcmLut.nCCD:
            useScienceDetectors = False
        else:
            # If the LUT has a different number of detectors than
            # the camera, then we only want to use science detectors
            # in the focal plane projector.
            useScienceDetectors = True

        focalPlaneProjector = FocalPlaneProjector(
            camera,
            self.config.defaultCameraOrientation,
            useScienceDetectors=useScienceDetectors,
        )

        noFitsDict = {'lutIndex': lutIndexVals,
                      'lutStd': lutStd,
                      'expInfo': fgcmExpInfo,
                      'focalPlaneProjector': focalPlaneProjector}

        # set up the fitter object
        fgcmFitCycle = fgcm.FgcmFitCycle(
            configDict,
            useFits=False,
            noFitsDict=noFitsDict,
            noOutput=True,
            butlerQC=butlerQC,
            plotHandleDict=plotHandleDict,
        )

        # create the parameter object
        if (fgcmFitCycle.initialCycle):
            # cycle = 0, initial cycle
            fgcmPars = fgcm.FgcmParameters.newParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                             fgcmLut,
                                                             fgcmExpInfo,
                                                             butlerQC=butlerQC,
                                                             plotHandleDict=plotHandleDict)
        else:
            if isinstance(handleDict['fgcmFitParameters'], afwTable.BaseCatalog):
                parCat = handleDict['fgcmFitParameters']
            else:
                parCat = handleDict['fgcmFitParameters'].get()
            inParInfo, inParams, inSuperStar = self._loadParameters(parCat)
            del parCat
            fgcmPars = fgcm.FgcmParameters.loadParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                              fgcmExpInfo,
                                                              inParInfo,
                                                              inParams,
                                                              inSuperStar,
                                                              butlerQC=butlerQC,
                                                              plotHandleDict=plotHandleDict)

        # set up the stars...
        fgcmStars = fgcm.FgcmStars(fgcmFitCycle.fgcmConfig, butlerQC=butlerQC, plotHandleDict=plotHandleDict)

        starObs = handleDict['fgcmStarObservations'].get()
        starIds = handleDict['fgcmStarIds'].get()
        if not self.config.useParquetCatalogFormat:
            starIndices = handleDict['fgcmStarIndices'].get()
        else:
            starIndices = None

        # grab the flagged stars if available
        if 'fgcmFlaggedStars' in handleDict:
            if isinstance(handleDict['fgcmFlaggedStars'], afwTable.BaseCatalog):
                flaggedStars = handleDict['fgcmFlaggedStars']
            else:
                flaggedStars = handleDict['fgcmFlaggedStars'].get()
            flagId = flaggedStars['objId'][:]
            flagFlag = flaggedStars['objFlag'][:]

            del flaggedStars
        elif self.config.useParquetCatalogFormat:
            # If we are using the parquet catalog format, then that means that
            # reserved stars have already been flagged.  We extract the flags here
            # to input to fgcm, which will then be persisted (with additional
            # quality flags) as the fgcmFlaggedStars datatype in subsequent
            # fit cycles.
            (flagged,) = (starIds['obj_flag'] > 0).nonzero()
            flagId = starIds['fgcm_id'][flagged]
            flagFlag = starIds['obj_flag'][flagged]
        else:
            flagId = None
            flagFlag = None

        if _config.doReferenceCalibration:
            refStars = handleDict['fgcmReferenceStars'].get()

            refMag, refMagErr = extractReferenceMags(refStars,
                                                     _config.bands,
                                                     _config.physicalFilterMap)

            refId = refStars['fgcm_id'][:]
        else:
            refStars = None
            refId = None
            refMag = None
            refMagErr = None

        # match star observations to visits
        # Only those star observations that match visits from fgcmExpInfo['VISIT'] will
        # actually be transferred into fgcm using the indexing below.
        if self.config.useParquetCatalogFormat:
            visitIndex = np.searchsorted(fgcmExpInfo['VISIT'], starObs['visit'])
        else:
            visitIndex = np.searchsorted(fgcmExpInfo['VISIT'], starObs['visit'][starIndices['obsIndex']])

        # The fgcmStars.loadStars method will copy all the star information into
        # special shared memory objects that will not blow up the memory usage when
        # used with python multiprocessing.  Once all the numbers are copied,
        # it is necessary to release all references to the objects that previously
        # stored the data to ensure that the garbage collector can clear the memory,
        # and ensure that this memory is not copied when multiprocessing kicks in.

        if self.config.useParquetCatalogFormat:
            # Note that the ra/dec coordinates for the parquet format are in
            # degrees, which is what fgcm expects.
            fgcmStars.loadStars(fgcmPars,
                                starObs['visit'][:],
                                starObs['detector'][:],
                                starObs['ra'][:],
                                starObs['dec'][:],
                                starObs['inst_mag'][:],
                                starObs['inst_mag_err'][:],
                                fgcmExpInfo['FILTERNAME'][visitIndex],
                                starIds['fgcm_id'][:],
                                starIds['ra'][:],
                                starIds['dec'][:],
                                starIds['obs_arr_index'][:],
                                starIds['n_obs'][:],
                                obsX=starObs['x'][:],
                                obsY=starObs['y'][:],
                                obsDeltaMagBkg=starObs['delta_mag_bkg'][:],
                                obsDeltaAper=starObs['delta_mag_aper'][:],
                                refID=refId,
                                refMag=refMag,
                                refMagErr=refMagErr,
                                flagID=flagId,
                                flagFlag=flagFlag,
                                computeNobs=True,
                                objIDAlternate=starIds['isolated_star_id'])
        else:
            # We determine the conversion from the native units (typically radians) to
            # degrees for the first star.  This allows us to treat coord_ra/coord_dec as
            # numpy arrays rather than Angles, which would we approximately 600x slower.
            conv = starObs[0]['ra'].asDegrees() / float(starObs[0]['ra'])

            fgcmStars.loadStars(fgcmPars,
                                starObs['visit'][starIndices['obsIndex']],
                                starObs['ccd'][starIndices['obsIndex']],
                                starObs['ra'][starIndices['obsIndex']] * conv,
                                starObs['dec'][starIndices['obsIndex']] * conv,
                                starObs['instMag'][starIndices['obsIndex']],
                                starObs['instMagErr'][starIndices['obsIndex']],
                                fgcmExpInfo['FILTERNAME'][visitIndex],
                                starIds['fgcm_id'][:],
                                starIds['ra'][:],
                                starIds['dec'][:],
                                starIds['obsArrIndex'][:],
                                starIds['nObs'][:],
                                obsX=starObs['x'][starIndices['obsIndex']],
                                obsY=starObs['y'][starIndices['obsIndex']],
                                obsDeltaMagBkg=starObs['deltaMagBkg'][starIndices['obsIndex']],
                                obsDeltaAper=starObs['deltaMagAper'][starIndices['obsIndex']],
                                psfCandidate=starObs['psf_candidate'][starIndices['obsIndex']],
                                refID=refId,
                                refMag=refMag,
                                refMagErr=refMagErr,
                                flagID=flagId,
                                flagFlag=flagFlag,
                                computeNobs=True)

        # Release all references to temporary objects holding star data (see above)
        del starObs
        del starIds
        del starIndices
        del flagId
        del flagFlag
        del refStars
        del refId
        del refMag
        del refMagErr

        # and set the bits in the cycle object
        fgcmFitCycle.setLUT(fgcmLut)
        fgcmFitCycle.setStars(fgcmStars, fgcmPars)
        fgcmFitCycle.setPars(fgcmPars)

        # finish the setup
        fgcmFitCycle.finishSetup()

        # and run
        fgcmFitCycle.run()

        ##################
        # Persistance
        ##################

        fgcmDatasetDict = self._makeFgcmOutputDatasets(fgcmFitCycle)

        # Output the config for the next cycle
        # We need to make a copy since the input one has been frozen

        updatedPhotometricCutDict = {b: float(fgcmFitCycle.updatedPhotometricCut[i]) for
                                     i, b in enumerate(_config.bands)}
        updatedHighCutDict = {band: float(fgcmFitCycle.updatedHighCut[i]) for
                              i, band in enumerate(_config.bands)}

        outConfig = copy.copy(_config)
        outConfig.update(cycleNumber=(_config.cycleNumber + 1),
                         precomputeSuperStarInitialCycle=False,
                         freezeStdAtmosphere=False,
                         expGrayPhotometricCutDict=updatedPhotometricCutDict,
                         expGrayHighCutDict=updatedHighCutDict)

        outConfig.connections.update(previousCycleNumber=str(_config.cycleNumber),
                                     cycleNumber=str(_config.cycleNumber + 1))

        if not multiCycle:
            configFileName = '%s_cycle%02d_config.py' % (outConfig.outfileBase,
                                                         outConfig.cycleNumber)
            outConfig.save(configFileName)

            if _config.isFinalCycle == 1:
                # We are done, ready to output products
                self.log.info("Everything is in place to run fgcmOutputProducts.py")
            else:
                self.log.info("Saved config for next cycle to %s" % (configFileName))
                self.log.info("Be sure to look at:")
                self.log.info("   config.expGrayPhotometricCut")
                self.log.info("   config.expGrayHighCut")
                self.log.info("If you are satisfied with the fit, please set:")
                self.log.info("   config.isFinalCycle = True")

        fgcmFitCycle.freeSharedMemory()

        return fgcmDatasetDict, outConfig

    def _loadParameters(self, parCat):
        """
        Load FGCM parameters from a previous fit cycle

        Parameters
        ----------
        parCat : `lsst.afw.table.BaseCatalog`
            Parameter catalog in afw table form.

        Returns
        -------
        inParInfo: `numpy.ndarray`
           Numpy array parameter information formatted for input to fgcm
        inParameters: `numpy.ndarray`
           Numpy array parameter values formatted for input to fgcm
        inSuperStar: `numpy.array`
           Superstar flat formatted for input to fgcm
        """
        parLutFilterNames = np.array(parCat[0]['lutFilterNames'].split(','))
        parFitBands = np.array(parCat[0]['fitBands'].split(','))

        inParInfo = np.zeros(1, dtype=[('NCCD', 'i4'),
                                       ('LUTFILTERNAMES', parLutFilterNames.dtype.str,
                                        (parLutFilterNames.size, )),
                                       ('FITBANDS', parFitBands.dtype.str, (parFitBands.size, )),
                                       ('LNTAUUNIT', 'f8'),
                                       ('LNTAUSLOPEUNIT', 'f8'),
                                       ('ALPHAUNIT', 'f8'),
                                       ('LNPWVUNIT', 'f8'),
                                       ('LNPWVSLOPEUNIT', 'f8'),
                                       ('LNPWVQUADRATICUNIT', 'f8'),
                                       ('LNPWVGLOBALUNIT', 'f8'),
                                       ('O3UNIT', 'f8'),
                                       ('QESYSUNIT', 'f8'),
                                       ('FILTEROFFSETUNIT', 'f8'),
                                       ('HASEXTERNALPWV', 'i2'),
                                       ('HASEXTERNALTAU', 'i2')])
        inParInfo['NCCD'] = parCat['nCcd']
        inParInfo['LUTFILTERNAMES'][:] = parLutFilterNames
        inParInfo['FITBANDS'][:] = parFitBands
        inParInfo['HASEXTERNALPWV'] = parCat['hasExternalPwv']
        inParInfo['HASEXTERNALTAU'] = parCat['hasExternalTau']

        inParams = np.zeros(1, dtype=[('PARALPHA', 'f8', (parCat['parAlpha'].size, )),
                                      ('PARO3', 'f8', (parCat['parO3'].size, )),
                                      ('PARLNTAUINTERCEPT', 'f8',
                                       (parCat['parLnTauIntercept'].size, )),
                                      ('PARLNTAUSLOPE', 'f8',
                                       (parCat['parLnTauSlope'].size, )),
                                      ('PARLNPWVINTERCEPT', 'f8',
                                       (parCat['parLnPwvIntercept'].size, )),
                                      ('PARLNPWVSLOPE', 'f8',
                                       (parCat['parLnPwvSlope'].size, )),
                                      ('PARLNPWVQUADRATIC', 'f8',
                                       (parCat['parLnPwvQuadratic'].size, )),
                                      ('PARQESYSINTERCEPT', 'f8',
                                       (parCat['parQeSysIntercept'].size, )),
                                      ('COMPQESYSSLOPE', 'f8',
                                       (parCat['compQeSysSlope'].size, )),
                                      ('PARFILTEROFFSET', 'f8',
                                       (parCat['parFilterOffset'].size, )),
                                      ('PARFILTEROFFSETFITFLAG', 'i2',
                                       (parCat['parFilterOffsetFitFlag'].size, )),
                                      ('PARRETRIEVEDLNPWVSCALE', 'f8'),
                                      ('PARRETRIEVEDLNPWVOFFSET', 'f8'),
                                      ('PARRETRIEVEDLNPWVNIGHTLYOFFSET', 'f8',
                                       (parCat['parRetrievedLnPwvNightlyOffset'].size, )),
                                      ('COMPABSTHROUGHPUT', 'f8',
                                       (parCat['compAbsThroughput'].size, )),
                                      ('COMPREFOFFSET', 'f8',
                                       (parCat['compRefOffset'].size, )),
                                      ('COMPREFSIGMA', 'f8',
                                       (parCat['compRefSigma'].size, )),
                                      ('COMPMIRRORCHROMATICITY', 'f8',
                                       (parCat['compMirrorChromaticity'].size, )),
                                      ('MIRRORCHROMATICITYPIVOT', 'f8',
                                       (parCat['mirrorChromaticityPivot'].size, )),
                                      ('COMPCCDCHROMATICITY', 'f8',
                                       (parCat['compCcdChromaticity'].size, )),
                                      ('COMPMEDIANSEDSLOPE', 'f8',
                                       (parCat['compMedianSedSlope'].size, )),
                                      ('COMPAPERCORRPIVOT', 'f8',
                                       (parCat['compAperCorrPivot'].size, )),
                                      ('COMPAPERCORRSLOPE', 'f8',
                                       (parCat['compAperCorrSlope'].size, )),
                                      ('COMPAPERCORRSLOPEERR', 'f8',
                                       (parCat['compAperCorrSlopeErr'].size, )),
                                      ('COMPAPERCORRRANGE', 'f8',
                                       (parCat['compAperCorrRange'].size, )),
                                      ('COMPMODELERREXPTIMEPIVOT', 'f8',
                                       (parCat['compModelErrExptimePivot'].size, )),
                                      ('COMPMODELERRFWHMPIVOT', 'f8',
                                       (parCat['compModelErrFwhmPivot'].size, )),
                                      ('COMPMODELERRSKYPIVOT', 'f8',
                                       (parCat['compModelErrSkyPivot'].size, )),
                                      ('COMPMODELERRPARS', 'f8',
                                       (parCat['compModelErrPars'].size, )),
                                      ('COMPEXPGRAY', 'f8',
                                       (parCat['compExpGray'].size, )),
                                      ('COMPVARGRAY', 'f8',
                                       (parCat['compVarGray'].size, )),
                                      ('COMPEXPDELTAMAGBKG', 'f8',
                                       (parCat['compExpDeltaMagBkg'].size, )),
                                      ('COMPNGOODSTARPEREXP', 'i4',
                                       (parCat['compNGoodStarPerExp'].size, )),
                                      ('COMPEXPREFOFFSET', 'f8',
                                       (parCat['compExpRefOffset'].size, )),
                                      ('COMPSIGFGCM', 'f8',
                                       (parCat['compSigFgcm'].size, )),
                                      ('COMPSIGMACAL', 'f8',
                                       (parCat['compSigmaCal'].size, )),
                                      ('COMPRETRIEVEDLNPWV', 'f8',
                                       (parCat['compRetrievedLnPwv'].size, )),
                                      ('COMPRETRIEVEDLNPWVRAW', 'f8',
                                       (parCat['compRetrievedLnPwvRaw'].size, )),
                                      ('COMPRETRIEVEDLNPWVFLAG', 'i2',
                                       (parCat['compRetrievedLnPwvFlag'].size, )),
                                      ('COMPRETRIEVEDTAUNIGHT', 'f8',
                                       (parCat['compRetrievedTauNight'].size, )),
                                      ('COMPEPSILON', 'f8',
                                       (parCat['compEpsilon'].size, )),
                                      ('COMPMEDDELTAAPER', 'f8',
                                       (parCat['compMedDeltaAper'].size, )),
                                      ('COMPGLOBALEPSILON', 'f4',
                                       (parCat['compGlobalEpsilon'].size, )),
                                      ('COMPEPSILONMAP', 'f4',
                                       (parCat['compEpsilonMap'].size, )),
                                      ('COMPEPSILONNSTARMAP', 'i4',
                                       (parCat['compEpsilonNStarMap'].size, )),
                                      ('COMPEPSILONCCDMAP', 'f4',
                                       (parCat['compEpsilonCcdMap'].size, )),
                                      ('COMPEPSILONCCDNSTARMAP', 'i4',
                                       (parCat['compEpsilonCcdNStarMap'].size, ))])

        inParams['PARALPHA'][:] = parCat['parAlpha'][0, :]
        inParams['PARO3'][:] = parCat['parO3'][0, :]
        inParams['PARLNTAUINTERCEPT'][:] = parCat['parLnTauIntercept'][0, :]
        inParams['PARLNTAUSLOPE'][:] = parCat['parLnTauSlope'][0, :]
        inParams['PARLNPWVINTERCEPT'][:] = parCat['parLnPwvIntercept'][0, :]
        inParams['PARLNPWVSLOPE'][:] = parCat['parLnPwvSlope'][0, :]
        inParams['PARLNPWVQUADRATIC'][:] = parCat['parLnPwvQuadratic'][0, :]
        inParams['PARQESYSINTERCEPT'][:] = parCat['parQeSysIntercept'][0, :]
        inParams['COMPQESYSSLOPE'][:] = parCat['compQeSysSlope'][0, :]
        inParams['PARFILTEROFFSET'][:] = parCat['parFilterOffset'][0, :]
        inParams['PARFILTEROFFSETFITFLAG'][:] = parCat['parFilterOffsetFitFlag'][0, :]
        inParams['PARRETRIEVEDLNPWVSCALE'] = parCat['parRetrievedLnPwvScale']
        inParams['PARRETRIEVEDLNPWVOFFSET'] = parCat['parRetrievedLnPwvOffset']
        inParams['PARRETRIEVEDLNPWVNIGHTLYOFFSET'][:] = parCat['parRetrievedLnPwvNightlyOffset'][0, :]
        inParams['COMPABSTHROUGHPUT'][:] = parCat['compAbsThroughput'][0, :]
        inParams['COMPREFOFFSET'][:] = parCat['compRefOffset'][0, :]
        inParams['COMPREFSIGMA'][:] = parCat['compRefSigma'][0, :]
        inParams['COMPMIRRORCHROMATICITY'][:] = parCat['compMirrorChromaticity'][0, :]
        inParams['MIRRORCHROMATICITYPIVOT'][:] = parCat['mirrorChromaticityPivot'][0, :]
        inParams['COMPCCDCHROMATICITY'][:] = parCat['compCcdChromaticity'][0, :]
        inParams['COMPMEDIANSEDSLOPE'][:] = parCat['compMedianSedSlope'][0, :]
        inParams['COMPAPERCORRPIVOT'][:] = parCat['compAperCorrPivot'][0, :]
        inParams['COMPAPERCORRSLOPE'][:] = parCat['compAperCorrSlope'][0, :]
        inParams['COMPAPERCORRSLOPEERR'][:] = parCat['compAperCorrSlopeErr'][0, :]
        inParams['COMPAPERCORRRANGE'][:] = parCat['compAperCorrRange'][0, :]
        inParams['COMPMODELERREXPTIMEPIVOT'][:] = parCat['compModelErrExptimePivot'][0, :]
        inParams['COMPMODELERRFWHMPIVOT'][:] = parCat['compModelErrFwhmPivot'][0, :]
        inParams['COMPMODELERRSKYPIVOT'][:] = parCat['compModelErrSkyPivot'][0, :]
        inParams['COMPMODELERRPARS'][:] = parCat['compModelErrPars'][0, :]
        inParams['COMPEXPGRAY'][:] = parCat['compExpGray'][0, :]
        inParams['COMPVARGRAY'][:] = parCat['compVarGray'][0, :]
        inParams['COMPEXPDELTAMAGBKG'][:] = parCat['compExpDeltaMagBkg'][0, :]
        inParams['COMPNGOODSTARPEREXP'][:] = parCat['compNGoodStarPerExp'][0, :]
        inParams['COMPEXPREFOFFSET'][:] = parCat['compExpRefOffset'][0, :]
        inParams['COMPSIGFGCM'][:] = parCat['compSigFgcm'][0, :]
        inParams['COMPSIGMACAL'][:] = parCat['compSigmaCal'][0, :]
        inParams['COMPRETRIEVEDLNPWV'][:] = parCat['compRetrievedLnPwv'][0, :]
        inParams['COMPRETRIEVEDLNPWVRAW'][:] = parCat['compRetrievedLnPwvRaw'][0, :]
        inParams['COMPRETRIEVEDLNPWVFLAG'][:] = parCat['compRetrievedLnPwvFlag'][0, :]
        inParams['COMPRETRIEVEDTAUNIGHT'][:] = parCat['compRetrievedTauNight'][0, :]
        inParams['COMPEPSILON'][:] = parCat['compEpsilon'][0, :]
        inParams['COMPMEDDELTAAPER'][:] = parCat['compMedDeltaAper'][0, :]
        inParams['COMPGLOBALEPSILON'][:] = parCat['compGlobalEpsilon'][0, :]
        inParams['COMPEPSILONMAP'][:] = parCat['compEpsilonMap'][0, :]
        inParams['COMPEPSILONNSTARMAP'][:] = parCat['compEpsilonNStarMap'][0, :]
        inParams['COMPEPSILONCCDMAP'][:] = parCat['compEpsilonCcdMap'][0, :]
        inParams['COMPEPSILONCCDNSTARMAP'][:] = parCat['compEpsilonCcdNStarMap'][0, :]

        inSuperStar = np.zeros(parCat['superstarSize'][0, :], dtype='f8')
        inSuperStar[:, :, :, :] = parCat['superstar'][0, :].reshape(inSuperStar.shape)

        return (inParInfo, inParams, inSuperStar)

    def _makeFgcmOutputDatasets(self, fgcmFitCycle):
        """
        Persist FGCM datasets through the butler.

        Parameters
        ----------
        fgcmFitCycle: `lsst.fgcm.FgcmFitCycle`
           Fgcm Fit cycle object
        """
        fgcmDatasetDict = {}

        # Save the parameters
        parInfo, pars = fgcmFitCycle.fgcmPars.parsToArrays()

        parSchema = afwTable.Schema()

        comma = ','
        lutFilterNameString = comma.join([n.decode('utf-8')
                                          for n in parInfo['LUTFILTERNAMES'][0]])
        fitBandString = comma.join([n.decode('utf-8')
                                    for n in parInfo['FITBANDS'][0]])

        parSchema = self._makeParSchema(parInfo, pars, fgcmFitCycle.fgcmPars.parSuperStarFlat,
                                        lutFilterNameString, fitBandString)
        parCat = self._makeParCatalog(parSchema, parInfo, pars,
                                      fgcmFitCycle.fgcmPars.parSuperStarFlat,
                                      lutFilterNameString, fitBandString)

        fgcmDatasetDict['fgcmFitParameters'] = parCat

        # Save the indices of the flagged stars
        # (stars that have been (a) reserved from the fit for testing and
        # (b) bad stars that have failed quality checks.)
        flagStarSchema = self._makeFlagStarSchema()
        flagStarStruct = fgcmFitCycle.fgcmStars.getFlagStarIndices()
        flagStarCat = self._makeFlagStarCat(flagStarSchema, flagStarStruct)

        fgcmDatasetDict['fgcmFlaggedStars'] = flagStarCat

        # Save the zeropoint information and atmospheres only if desired
        if self.outputZeropoints:
            superStarChebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_SSTAR_CHEB'].shape[1]
            zptChebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB'].shape[1]

            zptSchema = makeZptSchema(superStarChebSize, zptChebSize)
            zptCat = makeZptCat(zptSchema, fgcmFitCycle.fgcmZpts.zpStruct)

            fgcmDatasetDict['fgcmZeropoints'] = zptCat

            # Save atmosphere values
            # These are generated by the same code that generates zeropoints
            atmSchema = makeAtmSchema()
            atmCat = makeAtmCat(atmSchema, fgcmFitCycle.fgcmZpts.atmStruct)

            fgcmDatasetDict['fgcmAtmosphereParameters'] = atmCat

        # Save the standard stars (if configured)
        if self.outputStandards:
            stdStruct, goodBands = fgcmFitCycle.fgcmStars.retrieveStdStarCatalog(fgcmFitCycle.fgcmPars)
            stdSchema = makeStdSchema(len(goodBands))
            stdCat = makeStdCat(stdSchema, stdStruct, goodBands)

            fgcmDatasetDict['fgcmStandardStars'] = stdCat

        return fgcmDatasetDict

    def _makeParSchema(self, parInfo, pars, parSuperStarFlat,
                       lutFilterNameString, fitBandString):
        """
        Make the parameter persistence schema

        Parameters
        ----------
        parInfo: `numpy.ndarray`
           Parameter information returned by fgcm
        pars: `numpy.ndarray`
           Parameter values returned by fgcm
        parSuperStarFlat: `numpy.array`
           Superstar flat values returned by fgcm
        lutFilterNameString: `str`
           Combined string of all the lutFilterNames
        fitBandString: `str`
           Combined string of all the fitBands

        Returns
        -------
        parSchema: `afwTable.schema`
        """

        parSchema = afwTable.Schema()

        # parameter info section
        parSchema.addField('nCcd', type=np.int32, doc='Number of CCDs')
        parSchema.addField('lutFilterNames', type=str, doc='LUT Filter names in parameter file',
                           size=len(lutFilterNameString))
        parSchema.addField('fitBands', type=str, doc='Bands that were fit',
                           size=len(fitBandString))
        parSchema.addField('lnTauUnit', type=np.float64, doc='Step units for ln(AOD)')
        parSchema.addField('lnTauSlopeUnit', type=np.float64,
                           doc='Step units for ln(AOD) slope')
        parSchema.addField('alphaUnit', type=np.float64, doc='Step units for alpha')
        parSchema.addField('lnPwvUnit', type=np.float64, doc='Step units for ln(pwv)')
        parSchema.addField('lnPwvSlopeUnit', type=np.float64,
                           doc='Step units for ln(pwv) slope')
        parSchema.addField('lnPwvQuadraticUnit', type=np.float64,
                           doc='Step units for ln(pwv) quadratic term')
        parSchema.addField('lnPwvGlobalUnit', type=np.float64,
                           doc='Step units for global ln(pwv) parameters')
        parSchema.addField('o3Unit', type=np.float64, doc='Step units for O3')
        parSchema.addField('qeSysUnit', type=np.float64, doc='Step units for mirror gray')
        parSchema.addField('filterOffsetUnit', type=np.float64, doc='Step units for filter offset')
        parSchema.addField('hasExternalPwv', type=np.int32, doc='Parameters fit using external pwv')
        parSchema.addField('hasExternalTau', type=np.int32, doc='Parameters fit using external tau')

        # parameter section
        parSchema.addField('parAlpha', type='ArrayD', doc='Alpha parameter vector',
                           size=pars['PARALPHA'].size)
        parSchema.addField('parO3', type='ArrayD', doc='O3 parameter vector',
                           size=pars['PARO3'].size)
        parSchema.addField('parLnTauIntercept', type='ArrayD',
                           doc='ln(Tau) intercept parameter vector',
                           size=pars['PARLNTAUINTERCEPT'].size)
        parSchema.addField('parLnTauSlope', type='ArrayD',
                           doc='ln(Tau) slope parameter vector',
                           size=pars['PARLNTAUSLOPE'].size)
        parSchema.addField('parLnPwvIntercept', type='ArrayD', doc='ln(pwv) intercept parameter vector',
                           size=pars['PARLNPWVINTERCEPT'].size)
        parSchema.addField('parLnPwvSlope', type='ArrayD', doc='ln(pwv) slope parameter vector',
                           size=pars['PARLNPWVSLOPE'].size)
        parSchema.addField('parLnPwvQuadratic', type='ArrayD', doc='ln(pwv) quadratic parameter vector',
                           size=pars['PARLNPWVQUADRATIC'].size)
        parSchema.addField('parQeSysIntercept', type='ArrayD', doc='Mirror gray intercept parameter vector',
                           size=pars['PARQESYSINTERCEPT'].size)
        parSchema.addField('compQeSysSlope', type='ArrayD', doc='Mirror gray slope parameter vector',
                           size=pars[0]['COMPQESYSSLOPE'].size)
        parSchema.addField('parFilterOffset', type='ArrayD', doc='Filter offset parameter vector',
                           size=pars['PARFILTEROFFSET'].size)
        parSchema.addField('parFilterOffsetFitFlag', type='ArrayI', doc='Filter offset parameter fit flag',
                           size=pars['PARFILTEROFFSETFITFLAG'].size)
        parSchema.addField('parRetrievedLnPwvScale', type=np.float64,
                           doc='Global scale for retrieved ln(pwv)')
        parSchema.addField('parRetrievedLnPwvOffset', type=np.float64,
                           doc='Global offset for retrieved ln(pwv)')
        parSchema.addField('parRetrievedLnPwvNightlyOffset', type='ArrayD',
                           doc='Nightly offset for retrieved ln(pwv)',
                           size=pars['PARRETRIEVEDLNPWVNIGHTLYOFFSET'].size)
        parSchema.addField('compAbsThroughput', type='ArrayD',
                           doc='Absolute throughput (relative to transmission curves)',
                           size=pars['COMPABSTHROUGHPUT'].size)
        parSchema.addField('compRefOffset', type='ArrayD',
                           doc='Offset between reference stars and calibrated stars',
                           size=pars['COMPREFOFFSET'].size)
        parSchema.addField('compRefSigma', type='ArrayD',
                           doc='Width of reference star/calibrated star distribution',
                           size=pars['COMPREFSIGMA'].size)
        parSchema.addField('compMirrorChromaticity', type='ArrayD',
                           doc='Computed mirror chromaticity terms',
                           size=pars['COMPMIRRORCHROMATICITY'].size)
        parSchema.addField('mirrorChromaticityPivot', type='ArrayD',
                           doc='Mirror chromaticity pivot mjd',
                           size=pars['MIRRORCHROMATICITYPIVOT'].size)
        parSchema.addField('compCcdChromaticity', type='ArrayD',
                           doc='Computed CCD chromaticity terms',
                           size=pars['COMPCCDCHROMATICITY'].size)
        parSchema.addField('compMedianSedSlope', type='ArrayD',
                           doc='Computed median SED slope (per band)',
                           size=pars['COMPMEDIANSEDSLOPE'].size)
        parSchema.addField('compAperCorrPivot', type='ArrayD', doc='Aperture correction pivot',
                           size=pars['COMPAPERCORRPIVOT'].size)
        parSchema.addField('compAperCorrSlope', type='ArrayD', doc='Aperture correction slope',
                           size=pars['COMPAPERCORRSLOPE'].size)
        parSchema.addField('compAperCorrSlopeErr', type='ArrayD', doc='Aperture correction slope error',
                           size=pars['COMPAPERCORRSLOPEERR'].size)
        parSchema.addField('compAperCorrRange', type='ArrayD', doc='Aperture correction range',
                           size=pars['COMPAPERCORRRANGE'].size)
        parSchema.addField('compModelErrExptimePivot', type='ArrayD', doc='Model error exptime pivot',
                           size=pars['COMPMODELERREXPTIMEPIVOT'].size)
        parSchema.addField('compModelErrFwhmPivot', type='ArrayD', doc='Model error fwhm pivot',
                           size=pars['COMPMODELERRFWHMPIVOT'].size)
        parSchema.addField('compModelErrSkyPivot', type='ArrayD', doc='Model error sky pivot',
                           size=pars['COMPMODELERRSKYPIVOT'].size)
        parSchema.addField('compModelErrPars', type='ArrayD', doc='Model error parameters',
                           size=pars['COMPMODELERRPARS'].size)
        parSchema.addField('compExpGray', type='ArrayD', doc='Computed exposure gray',
                           size=pars['COMPEXPGRAY'].size)
        parSchema.addField('compVarGray', type='ArrayD', doc='Computed exposure variance',
                           size=pars['COMPVARGRAY'].size)
        parSchema.addField('compExpDeltaMagBkg', type='ArrayD',
                           doc='Computed exposure offset due to background',
                           size=pars['COMPEXPDELTAMAGBKG'].size)
        parSchema.addField('compNGoodStarPerExp', type='ArrayI',
                           doc='Computed number of good stars per exposure',
                           size=pars['COMPNGOODSTARPEREXP'].size)
        parSchema.addField('compExpRefOffset', type='ArrayD',
                           doc='Computed per-visit median offset between standard stars and ref stars.',
                           size=pars['COMPEXPREFOFFSET'].size)
        parSchema.addField('compSigFgcm', type='ArrayD', doc='Computed sigma_fgcm (intrinsic repeatability)',
                           size=pars['COMPSIGFGCM'].size)
        parSchema.addField('compSigmaCal', type='ArrayD', doc='Computed sigma_cal (systematic error floor)',
                           size=pars['COMPSIGMACAL'].size)
        parSchema.addField('compRetrievedLnPwv', type='ArrayD', doc='Retrieved ln(pwv) (smoothed)',
                           size=pars['COMPRETRIEVEDLNPWV'].size)
        parSchema.addField('compRetrievedLnPwvRaw', type='ArrayD', doc='Retrieved ln(pwv) (raw)',
                           size=pars['COMPRETRIEVEDLNPWVRAW'].size)
        parSchema.addField('compRetrievedLnPwvFlag', type='ArrayI', doc='Retrieved ln(pwv) Flag',
                           size=pars['COMPRETRIEVEDLNPWVFLAG'].size)
        parSchema.addField('compRetrievedTauNight', type='ArrayD', doc='Retrieved tau (per night)',
                           size=pars['COMPRETRIEVEDTAUNIGHT'].size)
        parSchema.addField('compEpsilon', type='ArrayD',
                           doc='Computed epsilon background offset per visit (nJy/arcsec2)',
                           size=pars['COMPEPSILON'].size)
        parSchema.addField('compMedDeltaAper', type='ArrayD',
                           doc='Median delta mag aper per visit',
                           size=pars['COMPMEDDELTAAPER'].size)
        parSchema.addField('compGlobalEpsilon', type='ArrayD',
                           doc='Computed epsilon bkg offset (global) (nJy/arcsec2)',
                           size=pars['COMPGLOBALEPSILON'].size)
        parSchema.addField('compEpsilonMap', type='ArrayD',
                           doc='Computed epsilon maps (nJy/arcsec2)',
                           size=pars['COMPEPSILONMAP'].size)
        parSchema.addField('compEpsilonNStarMap', type='ArrayI',
                           doc='Number of stars per pixel in computed epsilon maps',
                           size=pars['COMPEPSILONNSTARMAP'].size)
        parSchema.addField('compEpsilonCcdMap', type='ArrayD',
                           doc='Computed epsilon ccd maps (nJy/arcsec2)',
                           size=pars['COMPEPSILONCCDMAP'].size)
        parSchema.addField('compEpsilonCcdNStarMap', type='ArrayI',
                           doc='Number of stars per ccd bin in epsilon ccd maps',
                           size=pars['COMPEPSILONCCDNSTARMAP'].size)
        parSchema.addField('epochMjdStart', type='ArrayD',
                           doc='Epoch MJD start times',
                           size=pars['EPOCHMJDSTART'].size)
        parSchema.addField('epochMjdEnd', type='ArrayD',
                           doc='EpochMJD end times',
                           size=pars['EPOCHMJDEND'].size)
        # superstarflat section
        parSchema.addField('superstarSize', type='ArrayI', doc='Superstar matrix size',
                           size=4)
        parSchema.addField('superstar', type='ArrayD', doc='Superstar matrix (flattened)',
                           size=parSuperStarFlat.size)

        return parSchema

    def _makeParCatalog(self, parSchema, parInfo, pars, parSuperStarFlat,
                        lutFilterNameString, fitBandString):
        """
        Make the FGCM parameter catalog for persistence

        Parameters
        ----------
        parSchema: `lsst.afw.table.Schema`
           Parameter catalog schema
        pars: `numpy.ndarray`
           FGCM parameters to put into parCat
        parSuperStarFlat: `numpy.array`
           FGCM superstar flat array to put into parCat
        lutFilterNameString: `str`
           Combined string of all the lutFilterNames
        fitBandString: `str`
           Combined string of all the fitBands

        Returns
        -------
        parCat: `afwTable.BasicCatalog`
           Atmosphere and instrumental model parameter catalog for persistence
        """

        parCat = afwTable.BaseCatalog(parSchema)
        parCat.reserve(1)

        # The parameter catalog just has one row, with many columns for all the
        # atmosphere and instrument fit parameters
        rec = parCat.addNew()

        # info section
        rec['nCcd'] = parInfo['NCCD'][0]
        rec['lutFilterNames'] = lutFilterNameString
        rec['fitBands'] = fitBandString
        # note these are not currently supported here.
        rec['hasExternalPwv'] = 0
        rec['hasExternalTau'] = 0

        # parameter section

        scalarNames = ['parRetrievedLnPwvScale', 'parRetrievedLnPwvOffset']

        arrNames = ['parAlpha', 'parO3', 'parLnTauIntercept', 'parLnTauSlope',
                    'parLnPwvIntercept', 'parLnPwvSlope', 'parLnPwvQuadratic',
                    'parQeSysIntercept', 'compQeSysSlope',
                    'parRetrievedLnPwvNightlyOffset', 'compAperCorrPivot',
                    'parFilterOffset', 'parFilterOffsetFitFlag',
                    'compAbsThroughput', 'compRefOffset', 'compRefSigma',
                    'compMirrorChromaticity', 'mirrorChromaticityPivot', 'compCcdChromaticity',
                    'compAperCorrSlope', 'compAperCorrSlopeErr', 'compAperCorrRange',
                    'compModelErrExptimePivot', 'compModelErrFwhmPivot',
                    'compModelErrSkyPivot', 'compModelErrPars',
                    'compExpGray', 'compVarGray', 'compNGoodStarPerExp', 'compSigFgcm',
                    'compSigmaCal', 'compExpDeltaMagBkg', 'compMedianSedSlope',
                    'compRetrievedLnPwv', 'compRetrievedLnPwvRaw', 'compRetrievedLnPwvFlag',
                    'compRetrievedTauNight', 'compEpsilon', 'compMedDeltaAper',
                    'compGlobalEpsilon', 'compEpsilonMap', 'compEpsilonNStarMap',
                    'compEpsilonCcdMap', 'compEpsilonCcdNStarMap', 'compExpRefOffset',
                    'epochMjdStart', 'epochMjdEnd']

        for scalarName in scalarNames:
            rec[scalarName] = pars[scalarName.upper()][0]

        for arrName in arrNames:
            rec[arrName][:] = np.atleast_1d(pars[0][arrName.upper()])[:]

        # superstar section
        rec['superstarSize'][:] = parSuperStarFlat.shape
        rec['superstar'][:] = parSuperStarFlat.ravel()

        return parCat

    def _makeFlagStarSchema(self):
        """
        Make the flagged-stars schema

        Returns
        -------
        flagStarSchema: `lsst.afw.table.Schema`
        """

        flagStarSchema = afwTable.Schema()

        flagStarSchema.addField('objId', type=np.int32, doc='FGCM object id')
        flagStarSchema.addField('objFlag', type=np.int32, doc='FGCM object flag')

        return flagStarSchema

    def _makeFlagStarCat(self, flagStarSchema, flagStarStruct):
        """
        Make the flagged star catalog for persistence

        Parameters
        ----------
        flagStarSchema: `lsst.afw.table.Schema`
           Flagged star schema
        flagStarStruct: `numpy.ndarray`
           Flagged star structure from fgcm

        Returns
        -------
        flagStarCat: `lsst.afw.table.BaseCatalog`
           Flagged star catalog for persistence
        """

        flagStarCat = afwTable.BaseCatalog(flagStarSchema)
        flagStarCat.resize(flagStarStruct.size)

        flagStarCat['objId'][:] = flagStarStruct['OBJID']
        flagStarCat['objFlag'][:] = flagStarStruct['OBJFLAG']

        return flagStarCat
