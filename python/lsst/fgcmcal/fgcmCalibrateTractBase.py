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
"""Base class for running fgcmcal on a single tract using src tables
or sourceTable_visit tables.
"""

import abc

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .fgcmBuildStarsTable import FgcmBuildStarsTableTask
from .fgcmFitCycle import FgcmFitCycleConfig
from .fgcmOutputProducts import FgcmOutputProductsTask
from .utilities import makeConfigDict, translateFgcmLut, translateVisitCatalog
from .utilities import computeApertureRadiusFromName, extractReferenceMags
from .utilities import makeZptSchema, makeZptCat
from .utilities import makeAtmSchema, makeAtmCat
from .utilities import makeStdSchema, makeStdCat
from .focalPlaneProjector import FocalPlaneProjector

import fgcm

__all__ = ['FgcmCalibrateTractConfigBase', 'FgcmCalibrateTractBaseTask']


class FgcmCalibrateTractConfigBase(pexConfig.Config):
    """Config for FgcmCalibrateTract"""

    fgcmBuildStars = pexConfig.ConfigurableField(
        target=FgcmBuildStarsTableTask,
        doc="Task to load and match stars for fgcm",
    )
    fgcmFitCycle = pexConfig.ConfigField(
        dtype=FgcmFitCycleConfig,
        doc="Config to run a single fgcm fit cycle",
    )
    fgcmOutputProducts = pexConfig.ConfigurableField(
        target=FgcmOutputProductsTask,
        doc="Task to output fgcm products",
    )
    convergenceTolerance = pexConfig.Field(
        doc="Tolerance on repeatability convergence (per band)",
        dtype=float,
        default=0.005,
    )
    maxFitCycles = pexConfig.Field(
        doc="Maximum number of fit cycles",
        dtype=int,
        default=5,
    )
    doDebuggingPlots = pexConfig.Field(
        doc="Make plots for debugging purposes?",
        dtype=bool,
        default=False,
    )

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

        self.fgcmFitCycle.quietMode = True
        self.fgcmFitCycle.doPlots = False
        self.fgcmOutputProducts.doReferenceCalibration = False
        self.fgcmOutputProducts.doRefcatOutput = False
        self.fgcmOutputProducts.cycleNumber = 0
        self.fgcmOutputProducts.photoCal.applyColorTerms = False

    def validate(self):
        super().validate()

        for band in self.fgcmFitCycle.bands:
            if not self.fgcmFitCycle.useRepeatabilityForExpGrayCutsDict[band]:
                msg = 'Must set useRepeatabilityForExpGrayCutsDict[band]=True for all bands'
                raise pexConfig.FieldValidationError(FgcmFitCycleConfig.useRepeatabilityForExpGrayCutsDict,
                                                     self, msg)


class FgcmCalibrateTractBaseTask(pipeBase.PipelineTask, abc.ABC):
    """Base class to calibrate a single tract using fgcmcal
    """
    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("fgcmBuildStars", initInputs=initInputs)
        self.makeSubtask("fgcmOutputProducts")

    def run(self, handleDict, tract,
            buildStarsRefObjLoader=None, returnCatalogs=True):
        """Run the calibrations for a single tract with fgcm.

        Parameters
        ----------
        handleDict : `dict`
            All handles are `lsst.daf.butler.DeferredDatasetHandle`
            handle dictionary with the following keys.  Note that all
            keys need not be set based on config parameters.

            ``"camera"``
                Camera object (`lsst.afw.cameraGeom.Camera`)
            ``"source_catalogs"``
                `list` of handles for input source catalogs.
            ``"sourceSchema"``
                Schema for the source catalogs.
            ``"fgcmLookUpTable"``
                handle for the FGCM look-up table.
            ``"calexps"``
                `list` of handles for the input calexps
            ``"fgcmPhotoCalibs"``
                `dict` of output photoCalib handles.  Key is
                (tract, visit, detector).
                Present if doZeropointOutput is True.
            ``"fgcmTransmissionAtmospheres"``
                `dict` of output atmosphere transmission handles.
                Key is (tract, visit).
                Present if doAtmosphereOutput is True.
        tract : `int`
            Tract number
        buildStarsRefObjLoader : `lsst.meas.algorithms.ReferenceObjectLoader`, optional
            Reference object loader object for fgcmBuildStars.
        returnCatalogs : `bool`, optional
            Return photoCalibs as per-visit exposure catalogs.

        Returns
        -------
        outstruct : `lsst.pipe.base.Struct`
            Output structure with keys:

            offsets : `np.ndarray`
                Final reference offsets, per band.
            repeatability : `np.ndarray`
                Raw fgcm repeatability for bright stars, per band.
            atmospheres : `generator` [(`int`, `lsst.afw.image.TransmissionCurve`)]
                Generator that returns (visit, transmissionCurve) tuples.
            photoCalibs : `generator` [(`int`, `int`, `str`, `lsst.afw.image.PhotoCalib`)]
                Generator that returns (visit, ccd, filtername, photoCalib) tuples.
                (returned if returnCatalogs is False).
            photoCalibCatalogs : `generator` [(`int`, `lsst.afw.table.ExposureCatalog`)]
                Generator that returns (visit, exposureCatalog) tuples.
                (returned if returnCatalogs is True).
        """
        self.log.info("Running on tract %d", (tract))

        # Compute the aperture radius if necessary.  This is useful to do now before
        # any heavy lifting has happened (fail early).
        calibFluxApertureRadius = None
        if self.config.fgcmBuildStars.doSubtractLocalBackground:
            try:
                field = self.config.fgcmBuildStars.instFluxField
                calibFluxApertureRadius = computeApertureRadiusFromName(field)
            except RuntimeError:
                raise RuntimeError("Could not determine aperture radius from %s. "
                                   "Cannot use doSubtractLocalBackground." %
                                   (field))

        # Run the build stars tasks

        # Note that we will need visitCat at the end of the procedure for the outputs
        groupedHandles = self.fgcmBuildStars._groupHandles(handleDict['sourceTableHandleDict'],
                                                           handleDict['visitSummaryHandleDict'])
        visitCat = self.fgcmBuildStars.fgcmMakeVisitCatalog(handleDict['camera'], groupedHandles)
        rad = calibFluxApertureRadius
        fgcmStarObservationCat = self.fgcmBuildStars.fgcmMakeAllStarObservations(groupedHandles,
                                                                                 visitCat,
                                                                                 handleDict['sourceSchema'],
                                                                                 handleDict['camera'],
                                                                                 calibFluxApertureRadius=rad)

        if self.fgcmBuildStars.config.doReferenceMatches:
            lutHandle = handleDict['fgcmLookUpTable']
            self.fgcmBuildStars.makeSubtask("fgcmLoadReferenceCatalog",
                                            refObjLoader=buildStarsRefObjLoader)
        else:
            lutHandle = None

        fgcmStarIdCat, fgcmStarIndicesCat, fgcmRefCat = \
            self.fgcmBuildStars.fgcmMatchStars(visitCat,
                                               fgcmStarObservationCat,
                                               lutHandle=lutHandle)

        # Load the LUT
        lutCat = handleDict['fgcmLookUpTable'].get()
        fgcmLut, lutIndexVals, lutStd = translateFgcmLut(lutCat,
                                                         dict(self.config.fgcmFitCycle.physicalFilterMap))
        del lutCat

        # Translate the visit catalog into fgcm format
        fgcmExpInfo = translateVisitCatalog(visitCat)

        configDict = makeConfigDict(self.config.fgcmFitCycle, self.log, handleDict['camera'],
                                    self.config.fgcmFitCycle.maxIterBeforeFinalCycle,
                                    True, False, lutIndexVals[0]['FILTERNAMES'],
                                    tract=tract)

        focalPlaneProjector = FocalPlaneProjector(handleDict['camera'],
                                                  self.config.fgcmFitCycle.defaultCameraOrientation)

        # Set up the fit cycle task

        noFitsDict = {'lutIndex': lutIndexVals,
                      'lutStd': lutStd,
                      'expInfo': fgcmExpInfo,
                      'focalPlaneProjector': focalPlaneProjector}

        fgcmFitCycle = fgcm.FgcmFitCycle(configDict, useFits=False,
                                         noFitsDict=noFitsDict, noOutput=True)

        # We determine the conversion from the native units (typically radians) to
        # degrees for the first star.  This allows us to treat coord_ra/coord_dec as
        # numpy arrays rather than Angles, which would we approximately 600x slower.
        conv = fgcmStarObservationCat[0]['ra'].asDegrees() / float(fgcmStarObservationCat[0]['ra'])

        # To load the stars, we need an initial parameter object
        fgcmPars = fgcm.FgcmParameters.newParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                         fgcmLut,
                                                         fgcmExpInfo)

        # Match star observations to visits
        # Only those star observations that match visits from fgcmExpInfo['VISIT'] will
        # actually be transferred into fgcm using the indexing below.

        obsIndex = fgcmStarIndicesCat['obsIndex']
        visitIndex = np.searchsorted(fgcmExpInfo['VISIT'],
                                     fgcmStarObservationCat['visit'][obsIndex])

        refMag, refMagErr = extractReferenceMags(fgcmRefCat,
                                                 self.config.fgcmFitCycle.bands,
                                                 self.config.fgcmFitCycle.physicalFilterMap)
        refId = fgcmRefCat['fgcm_id'][:]

        fgcmStars = fgcm.FgcmStars(fgcmFitCycle.fgcmConfig)
        fgcmStars.loadStars(fgcmPars,
                            fgcmStarObservationCat['visit'][obsIndex],
                            fgcmStarObservationCat['ccd'][obsIndex],
                            fgcmStarObservationCat['ra'][obsIndex] * conv,
                            fgcmStarObservationCat['dec'][obsIndex] * conv,
                            fgcmStarObservationCat['instMag'][obsIndex],
                            fgcmStarObservationCat['instMagErr'][obsIndex],
                            fgcmExpInfo['FILTERNAME'][visitIndex],
                            fgcmStarIdCat['fgcm_id'][:],
                            fgcmStarIdCat['ra'][:],
                            fgcmStarIdCat['dec'][:],
                            fgcmStarIdCat['obsArrIndex'][:],
                            fgcmStarIdCat['nObs'][:],
                            obsX=fgcmStarObservationCat['x'][obsIndex],
                            obsY=fgcmStarObservationCat['y'][obsIndex],
                            obsDeltaMagBkg=fgcmStarObservationCat['deltaMagBkg'][obsIndex],
                            obsDeltaAper=fgcmStarObservationCat['deltaMagAper'][obsIndex],
                            psfCandidate=fgcmStarObservationCat['psf_candidate'][obsIndex],
                            refID=refId,
                            refMag=refMag,
                            refMagErr=refMagErr,
                            flagID=None,
                            flagFlag=None,
                            computeNobs=True)

        # Clear out some memory
        del fgcmStarIdCat
        del fgcmStarIndicesCat
        del fgcmRefCat

        fgcmFitCycle.setLUT(fgcmLut)
        fgcmFitCycle.setStars(fgcmStars, fgcmPars)

        converged = False
        cycleNumber = 0

        previousReservedRawRepeatability = np.zeros(fgcmPars.nBands) + 1000.0
        previousParInfo = None
        previousParams = None
        previousSuperStar = None

        while (not converged and cycleNumber < self.config.maxFitCycles):

            fgcmFitCycle.fgcmConfig.updateCycleNumber(cycleNumber)

            if cycleNumber > 0:
                # Use parameters from previous cycle
                fgcmPars = fgcm.FgcmParameters.loadParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                                  fgcmExpInfo,
                                                                  previousParInfo,
                                                                  previousParams,
                                                                  previousSuperStar)
                # We need to reset the star magnitudes and errors for the next
                # cycle
                fgcmFitCycle.fgcmStars.reloadStarMagnitudes(fgcmStarObservationCat['instMag'][obsIndex],
                                                            fgcmStarObservationCat['instMagErr'][obsIndex])
                fgcmFitCycle.initialCycle = False

            fgcmFitCycle.setPars(fgcmPars)
            fgcmFitCycle.finishSetup()

            fgcmFitCycle.run()

            # Grab the parameters for the next cycle
            previousParInfo, previousParams = fgcmFitCycle.fgcmPars.parsToArrays()
            previousSuperStar = fgcmFitCycle.fgcmPars.parSuperStarFlat.copy()

            self.log.info("Raw repeatability after cycle number %d is:" % (cycleNumber))
            for i, band in enumerate(fgcmFitCycle.fgcmPars.bands):
                if not fgcmFitCycle.fgcmPars.hasExposuresInBand[i]:
                    continue
                rep = fgcmFitCycle.fgcmPars.compReservedRawRepeatability[i] * 1000.0
                self.log.info("  Band %s, repeatability: %.2f mmag" % (band, rep))

            # Check for convergence
            if np.alltrue((previousReservedRawRepeatability
                           - fgcmFitCycle.fgcmPars.compReservedRawRepeatability)
                          < self.config.convergenceTolerance):
                self.log.info("Raw repeatability has converged after cycle number %d." % (cycleNumber))
                converged = True
            else:
                fgcmFitCycle.fgcmConfig.expGrayPhotometricCut[:] = fgcmFitCycle.updatedPhotometricCut
                fgcmFitCycle.fgcmConfig.expGrayHighCut[:] = fgcmFitCycle.updatedHighCut
                fgcmFitCycle.fgcmConfig.precomputeSuperStarInitialCycle = False
                fgcmFitCycle.fgcmConfig.freezeStdAtmosphere = False
                previousReservedRawRepeatability[:] = fgcmFitCycle.fgcmPars.compReservedRawRepeatability
                self.log.info("Setting exposure gray photometricity cuts to:")
                for i, band in enumerate(fgcmFitCycle.fgcmPars.bands):
                    if not fgcmFitCycle.fgcmPars.hasExposuresInBand[i]:
                        continue
                    cut = fgcmFitCycle.updatedPhotometricCut[i] * 1000.0
                    self.log.info("  Band %s, photometricity cut: %.2f mmag" % (band, cut))

            cycleNumber += 1

        # Log warning if not converged
        if not converged:
            self.log.warning("Maximum number of fit cycles exceeded (%d) without convergence.", cycleNumber)

        # Do final clean-up iteration
        fgcmFitCycle.fgcmConfig.freezeStdAtmosphere = False
        fgcmFitCycle.fgcmConfig.resetParameters = False
        fgcmFitCycle.fgcmConfig.maxIter = 0
        fgcmFitCycle.fgcmConfig.outputZeropoints = True
        fgcmFitCycle.fgcmConfig.outputStandards = True
        fgcmFitCycle.fgcmConfig.doPlots = self.config.doDebuggingPlots
        fgcmFitCycle.fgcmConfig.updateCycleNumber(cycleNumber)
        fgcmFitCycle.initialCycle = False

        fgcmPars = fgcm.FgcmParameters.loadParsWithArrays(fgcmFitCycle.fgcmConfig,
                                                          fgcmExpInfo,
                                                          previousParInfo,
                                                          previousParams,
                                                          previousSuperStar)
        fgcmFitCycle.fgcmStars.reloadStarMagnitudes(fgcmStarObservationCat['instMag'][obsIndex],
                                                    fgcmStarObservationCat['instMagErr'][obsIndex])
        fgcmFitCycle.setPars(fgcmPars)
        fgcmFitCycle.finishSetup()

        self.log.info("Running final clean-up fit cycle...")
        fgcmFitCycle.run()

        self.log.info("Raw repeatability after clean-up cycle is:")
        for i, band in enumerate(fgcmFitCycle.fgcmPars.bands):
            if not fgcmFitCycle.fgcmPars.hasExposuresInBand[i]:
                continue
            rep = fgcmFitCycle.fgcmPars.compReservedRawRepeatability[i] * 1000.0
            self.log.info("  Band %s, repeatability: %.2f mmag" % (band, rep))

        # Do the outputs.  Need to keep track of tract.

        superStarChebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_SSTAR_CHEB'].shape[1]
        zptChebSize = fgcmFitCycle.fgcmZpts.zpStruct['FGCM_FZPT_CHEB'].shape[1]

        zptSchema = makeZptSchema(superStarChebSize, zptChebSize)
        zptCat = makeZptCat(zptSchema, fgcmFitCycle.fgcmZpts.zpStruct)

        atmSchema = makeAtmSchema()
        atmCat = makeAtmCat(atmSchema, fgcmFitCycle.fgcmZpts.atmStruct)

        stdStruct, goodBands = fgcmFitCycle.fgcmStars.retrieveStdStarCatalog(fgcmFitCycle.fgcmPars)
        stdSchema = makeStdSchema(len(goodBands))
        stdCat = makeStdCat(stdSchema, stdStruct, goodBands)

        outStruct = self.fgcmOutputProducts.generateTractOutputProducts(handleDict,
                                                                        tract,
                                                                        visitCat,
                                                                        zptCat, atmCat, stdCat,
                                                                        self.config.fgcmBuildStars)

        outStruct.repeatability = fgcmFitCycle.fgcmPars.compReservedRawRepeatability

        fgcmFitCycle.freeSharedMemory()

        return outStruct
