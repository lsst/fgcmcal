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
"""Utility functions for fgcmcal.

This file contains utility functions that are used by more than one task,
and do not need to be part of a task.
"""

import numpy as np
import os
import re

from lsst.daf.base import PropertyList
from lsst.daf.butler import Timespan
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom as geom
from lsst.obs.base import createInitialSkyWcs

import fgcm


FGCM_EXP_FIELD = 'VISIT'
FGCM_CCD_FIELD = 'DETECTOR'
FGCM_ILLEGAL_VALUE = -9999.0


def makeConfigDict(config, log, camera, maxIter,
                   resetFitParameters, outputZeropoints,
                   lutFilterNames, tract=None, nCore=1, doPlots=False):
    """
    Make the FGCM fit cycle configuration dict

    Parameters
    ----------
    config : `lsst.fgcmcal.FgcmFitCycleConfig`
        Configuration object
    log : `logging.Logger`
        Log object.
    camera : `lsst.afw.cameraGeom.Camera`
        Camera from the butler
    maxIter : `int`
        Maximum number of iterations
    resetFitParameters: `bool`
        Reset fit parameters before fitting?
    outputZeropoints : `bool`
        Compute zeropoints for output?
    lutFilterNames : array-like, `str`
        Array of physical filter names in the LUT.
    tract : `int`, optional
        Tract number for extending the output file name for debugging.
        Default is None.
    nCore : `int`, optional
        Number of cores to use.
    doPlots : `bool`, optional
        Make FGCM QA plots?

    Returns
    -------
    configDict : `dict`
        Configuration dictionary for fgcm
    """
    # Extract the bands that are _not_ being fit for fgcm configuration
    notFitBands = [b for b in config.bands if b not in config.fitBands]

    # process the starColorCuts
    starColorCutList = []
    for ccut in config.starColorCuts:
        if ccut == 'NO_DATA':
            # No color cuts to apply.
            break
        parts = ccut.split(',')
        starColorCutList.append([parts[0], parts[1], float(parts[2]), float(parts[3])])

    # process the refStarColorCuts
    refStarColorCutList = []
    for ccut in config.refStarColorCuts:
        if ccut == 'NO_DATA':
            # No color cuts to apply.
            break
        parts = ccut.split(',')
        refStarColorCutList.append([parts[0], parts[1], float(parts[2]), float(parts[3])])

    # TODO: Having direct access to the mirror area from the camera would be
    #  useful.  See DM-16489.
    # Mirror area in cm**2
    if config.mirrorArea is None:
        mirrorArea = np.pi*(camera.telescopeDiameter*100./2.)**2.
    else:
        # Convert to square cm.
        mirrorArea = config.mirrorArea * 100.**2.

    # Get approximate average camera gain:
    gains = [amp.getGain() for detector in camera for amp in detector.getAmplifiers()]
    cameraGain = float(np.median(gains))

    # Cut down the filter map to those that are in the LUT
    filterToBand = {filterName: config.physicalFilterMap[filterName] for
                    filterName in lutFilterNames}

    if tract is None:
        outfileBase = config.outfileBase
    else:
        outfileBase = '%s-%06d' % (config.outfileBase, tract)

    # create a configuration dictionary for fgcmFitCycle
    configDict = {'outfileBase': outfileBase,
                  'logger': log,
                  'exposureFile': None,
                  'obsFile': None,
                  'indexFile': None,
                  'lutFile': None,
                  'mirrorArea': mirrorArea,
                  'cameraGain': cameraGain,
                  'ccdStartIndex': camera[0].getId(),
                  'expField': FGCM_EXP_FIELD,
                  'ccdField': FGCM_CCD_FIELD,
                  'seeingField': 'DELTA_APER',
                  'fwhmField': 'PSFSIGMA',
                  'skyBrightnessField': 'SKYBACKGROUND',
                  'deepFlag': 'DEEPFLAG',  # unused
                  'bands': list(config.bands),
                  'fitBands': list(config.fitBands),
                  'notFitBands': notFitBands,
                  'requiredBands': list(config.requiredBands),
                  'filterToBand': filterToBand,
                  'logLevel': 'INFO',
                  'nCore': nCore,
                  'nStarPerRun': config.nStarPerRun,
                  'nExpPerRun': config.nExpPerRun,
                  'reserveFraction': config.reserveFraction,
                  'freezeStdAtmosphere': config.freezeStdAtmosphere,
                  'precomputeSuperStarInitialCycle': config.precomputeSuperStarInitialCycle,
                  'superStarSubCCDDict': dict(config.superStarSubCcdDict),
                  'superStarSubCCDChebyshevOrder': config.superStarSubCcdChebyshevOrder,
                  'superStarSubCCDTriangular': config.superStarSubCcdTriangular,
                  'superStarSigmaClip': config.superStarSigmaClip,
                  'superStarPlotCCDResiduals': config.superStarPlotCcdResiduals,
                  'focalPlaneSigmaClip': config.focalPlaneSigmaClip,
                  'ccdGraySubCCDDict': dict(config.ccdGraySubCcdDict),
                  'ccdGraySubCCDChebyshevOrder': config.ccdGraySubCcdChebyshevOrder,
                  'ccdGraySubCCDTriangular': config.ccdGraySubCcdTriangular,
                  'ccdGrayFocalPlaneDict': dict(config.ccdGrayFocalPlaneDict),
                  'ccdGrayFocalPlaneChebyshevOrder': config.ccdGrayFocalPlaneChebyshevOrder,
                  'ccdGrayFocalPlaneFitMinCcd': config.ccdGrayFocalPlaneFitMinCcd,
                  'cycleNumber': config.cycleNumber,
                  'maxIter': maxIter,
                  'deltaMagBkgOffsetPercentile': config.deltaMagBkgOffsetPercentile,
                  'deltaMagBkgPerCcd': config.deltaMagBkgPerCcd,
                  'UTBoundary': config.utBoundary,
                  'washMJDs': config.washMjds,
                  'epochMJDs': config.epochMjds,
                  'coatingMJDs': config.coatingMjds,
                  'minObsPerBand': config.minObsPerBand,
                  'latitude': config.latitude,
                  'defaultCameraOrientation': config.defaultCameraOrientation,
                  'brightObsGrayMax': config.brightObsGrayMax,
                  'minStarPerCCD': config.minStarPerCcd,
                  'minCCDPerExp': config.minCcdPerExp,
                  'maxCCDGrayErr': config.maxCcdGrayErr,
                  'minStarPerExp': config.minStarPerExp,
                  'minExpPerNight': config.minExpPerNight,
                  'expGrayInitialCut': config.expGrayInitialCut,
                  'expGrayPhotometricCutDict': dict(config.expGrayPhotometricCutDict),
                  'expGrayHighCutDict': dict(config.expGrayHighCutDict),
                  'expGrayRecoverCut': config.expGrayRecoverCut,
                  'expVarGrayPhotometricCutDict': dict(config.expVarGrayPhotometricCutDict),
                  'expGrayErrRecoverCut': config.expGrayErrRecoverCut,
                  'refStarSnMin': config.refStarSnMin,
                  'refStarOutlierNSig': config.refStarOutlierNSig,
                  'applyRefStarColorCuts': config.applyRefStarColorCuts,
                  'refStarMaxFracUse': config.refStarMaxFracUse,
                  'useExposureReferenceOffset': config.useExposureReferenceOffset,
                  'illegalValue': FGCM_ILLEGAL_VALUE,  # internally used by fgcm.
                  'starColorCuts': starColorCutList,
                  'refStarColorCuts': refStarColorCutList,
                  'aperCorrFitNBins': config.aperCorrFitNBins,
                  'aperCorrInputSlopeDict': dict(config.aperCorrInputSlopeDict),
                  'sedBoundaryTermDict': config.sedboundaryterms.toDict()['data'],
                  'sedTermDict': config.sedterms.toDict()['data'],
                  'colorSplitBands': list(config.colorSplitBands),
                  'sigFgcmMaxErr': config.sigFgcmMaxErr,
                  'sigFgcmMaxEGrayDict': dict(config.sigFgcmMaxEGrayDict),
                  'ccdGrayMaxStarErr': config.ccdGrayMaxStarErr,
                  'approxThroughputDict': dict(config.approxThroughputDict),
                  'sigmaCalRange': list(config.sigmaCalRange),
                  'sigmaCalFitPercentile': list(config.sigmaCalFitPercentile),
                  'sigmaCalPlotPercentile': list(config.sigmaCalPlotPercentile),
                  'sigma0Phot': config.sigma0Phot,
                  'mapLongitudeRef': config.mapLongitudeRef,
                  'mapNSide': config.mapNSide,
                  'varNSig': 100.0,  # Turn off 'variable star selection' which doesn't work yet
                  'varMinBand': 2,
                  'useRetrievedPwv': False,
                  'useNightlyRetrievedPwv': False,
                  'pwvRetrievalSmoothBlock': 25,
                  'useQuadraticPwv': config.useQuadraticPwv,
                  'useRetrievedTauInit': False,
                  'tauRetrievalMinCCDPerNight': 500,
                  'modelMagErrors': config.modelMagErrors,
                  'instrumentParsPerBand': config.instrumentParsPerBand,
                  'instrumentSlopeMinDeltaT': config.instrumentSlopeMinDeltaT,
                  'fitMirrorChromaticity': config.fitMirrorChromaticity,
                  'fitCCDChromaticityDict': dict(config.fitCcdChromaticityDict),
                  'useRepeatabilityForExpGrayCutsDict': dict(config.useRepeatabilityForExpGrayCutsDict),
                  'autoPhotometricCutNSig': config.autoPhotometricCutNSig,
                  'autoHighCutNSig': config.autoHighCutNSig,
                  'deltaAperInnerRadiusArcsec': config.deltaAperInnerRadiusArcsec,
                  'deltaAperOuterRadiusArcsec': config.deltaAperOuterRadiusArcsec,
                  'deltaAperFitMinNgoodObs': config.deltaAperFitMinNgoodObs,
                  'deltaAperFitPerCcdNx': config.deltaAperFitPerCcdNx,
                  'deltaAperFitPerCcdNy': config.deltaAperFitPerCcdNy,
                  'deltaAperFitSpatialNside': config.deltaAperFitSpatialNside,
                  'doComputeDeltaAperExposures': config.doComputeDeltaAperPerVisit,
                  'doComputeDeltaAperStars': config.doComputeDeltaAperPerStar,
                  'doComputeDeltaAperMap': config.doComputeDeltaAperMap,
                  'doComputeDeltaAperPerCcd': config.doComputeDeltaAperPerCcd,
                  'printOnly': False,
                  'quietMode': config.quietMode,
                  'randomSeed': config.randomSeed,
                  'outputStars': False,
                  'outputPath': os.path.abspath('.'),
                  'clobber': True,
                  'useSedLUT': False,
                  'resetParameters': resetFitParameters,
                  'doPlots': doPlots,
                  'outputFgcmcalZpts': True,  # when outputting zpts, use fgcmcal format
                  'outputZeropoints': outputZeropoints}

    return configDict


def translateFgcmLut(lutCat, physicalFilterMap):
    """
    Translate the FGCM look-up-table into an fgcm-compatible object

    Parameters
    ----------
    lutCat: `lsst.afw.table.BaseCatalog`
       Catalog describing the FGCM look-up table
    physicalFilterMap: `dict`
       Physical filter to band mapping

    Returns
    -------
    fgcmLut: `lsst.fgcm.FgcmLut`
       Lookup table for FGCM
    lutIndexVals: `numpy.ndarray`
       Numpy array with LUT index information for FGCM
    lutStd: `numpy.ndarray`
       Numpy array with LUT standard throughput values for FGCM

    Notes
    -----
    After running this code, it is wise to `del lutCat` to clear the memory.
    """

    # first we need the lutIndexVals
    lutFilterNames = np.array(lutCat[0]['physicalFilters'].split(','), dtype='U')
    lutStdFilterNames = np.array(lutCat[0]['stdPhysicalFilters'].split(','), dtype='U')

    # Note that any discrepancies between config values will raise relevant
    # exceptions in the FGCM code.

    lutIndexVals = np.zeros(1, dtype=[('FILTERNAMES', lutFilterNames.dtype.str,
                                       lutFilterNames.size),
                                      ('STDFILTERNAMES', lutStdFilterNames.dtype.str,
                                       lutStdFilterNames.size),
                                      ('PMB', 'f8', lutCat[0]['pmb'].size),
                                      ('PMBFACTOR', 'f8', lutCat[0]['pmbFactor'].size),
                                      ('PMBELEVATION', 'f8'),
                                      ('LAMBDANORM', 'f8'),
                                      ('PWV', 'f8', lutCat[0]['pwv'].size),
                                      ('O3', 'f8', lutCat[0]['o3'].size),
                                      ('TAU', 'f8', lutCat[0]['tau'].size),
                                      ('ALPHA', 'f8', lutCat[0]['alpha'].size),
                                      ('ZENITH', 'f8', lutCat[0]['zenith'].size),
                                      ('NCCD', 'i4')])

    lutIndexVals['FILTERNAMES'][:] = lutFilterNames
    lutIndexVals['STDFILTERNAMES'][:] = lutStdFilterNames
    lutIndexVals['PMB'][:] = lutCat[0]['pmb']
    lutIndexVals['PMBFACTOR'][:] = lutCat[0]['pmbFactor']
    lutIndexVals['PMBELEVATION'] = lutCat[0]['pmbElevation']
    lutIndexVals['LAMBDANORM'] = lutCat[0]['lambdaNorm']
    lutIndexVals['PWV'][:] = lutCat[0]['pwv']
    lutIndexVals['O3'][:] = lutCat[0]['o3']
    lutIndexVals['TAU'][:] = lutCat[0]['tau']
    lutIndexVals['ALPHA'][:] = lutCat[0]['alpha']
    lutIndexVals['ZENITH'][:] = lutCat[0]['zenith']
    lutIndexVals['NCCD'] = lutCat[0]['nCcd']

    # now we need the Standard Values
    lutStd = np.zeros(1, dtype=[('PMBSTD', 'f8'),
                                ('PWVSTD', 'f8'),
                                ('O3STD', 'f8'),
                                ('TAUSTD', 'f8'),
                                ('ALPHASTD', 'f8'),
                                ('ZENITHSTD', 'f8'),
                                ('LAMBDARANGE', 'f8', 2),
                                ('LAMBDASTEP', 'f8'),
                                ('LAMBDASTD', 'f8', lutFilterNames.size),
                                ('LAMBDASTDFILTER', 'f8', lutStdFilterNames.size),
                                ('I0STD', 'f8', lutFilterNames.size),
                                ('I1STD', 'f8', lutFilterNames.size),
                                ('I10STD', 'f8', lutFilterNames.size),
                                ('I2STD', 'f8', lutFilterNames.size),
                                ('LAMBDAB', 'f8', lutFilterNames.size),
                                ('ATMLAMBDA', 'f8', lutCat[0]['atmLambda'].size),
                                ('ATMSTDTRANS', 'f8', lutCat[0]['atmStdTrans'].size)])
    lutStd['PMBSTD'] = lutCat[0]['pmbStd']
    lutStd['PWVSTD'] = lutCat[0]['pwvStd']
    lutStd['O3STD'] = lutCat[0]['o3Std']
    lutStd['TAUSTD'] = lutCat[0]['tauStd']
    lutStd['ALPHASTD'] = lutCat[0]['alphaStd']
    lutStd['ZENITHSTD'] = lutCat[0]['zenithStd']
    lutStd['LAMBDARANGE'][:] = lutCat[0]['lambdaRange'][:]
    lutStd['LAMBDASTEP'] = lutCat[0]['lambdaStep']
    lutStd['LAMBDASTD'][:] = lutCat[0]['lambdaStd']
    lutStd['LAMBDASTDFILTER'][:] = lutCat[0]['lambdaStdFilter']
    lutStd['I0STD'][:] = lutCat[0]['i0Std']
    lutStd['I1STD'][:] = lutCat[0]['i1Std']
    lutStd['I10STD'][:] = lutCat[0]['i10Std']
    lutStd['I2STD'][:] = lutCat[0]['i2Std']
    lutStd['LAMBDAB'][:] = lutCat[0]['lambdaB']
    lutStd['ATMLAMBDA'][:] = lutCat[0]['atmLambda'][:]
    lutStd['ATMSTDTRANS'][:] = lutCat[0]['atmStdTrans'][:]

    lutTypes = [row['luttype'] for row in lutCat]

    # And the flattened look-up-table
    lutFlat = np.zeros(lutCat[0]['lut'].size, dtype=[('I0', 'f4'),
                                                     ('I1', 'f4')])

    lutFlat['I0'][:] = lutCat[lutTypes.index('I0')]['lut'][:]
    lutFlat['I1'][:] = lutCat[lutTypes.index('I1')]['lut'][:]

    lutDerivFlat = np.zeros(lutCat[0]['lut'].size, dtype=[('D_LNPWV', 'f4'),
                                                          ('D_O3', 'f4'),
                                                          ('D_LNTAU', 'f4'),
                                                          ('D_ALPHA', 'f4'),
                                                          ('D_SECZENITH', 'f4'),
                                                          ('D_LNPWV_I1', 'f4'),
                                                          ('D_O3_I1', 'f4'),
                                                          ('D_LNTAU_I1', 'f4'),
                                                          ('D_ALPHA_I1', 'f4'),
                                                          ('D_SECZENITH_I1', 'f4')])

    for name in lutDerivFlat.dtype.names:
        lutDerivFlat[name][:] = lutCat[lutTypes.index(name)]['lut'][:]

    # The fgcm.FgcmLUT() class copies all the LUT information into special
    # shared memory objects that will not blow up the memory usage when used
    # with python multiprocessing.  Once all the numbers are copied, the
    # references to the temporary objects (lutCat, lutFlat, lutDerivFlat)
    # will fall out of scope and can be cleaned up by the garbage collector.
    fgcmLut = fgcm.FgcmLUT(lutIndexVals, lutFlat, lutDerivFlat, lutStd,
                           filterToBand=physicalFilterMap)

    return fgcmLut, lutIndexVals, lutStd


def translateVisitCatalog(visitCat):
    """
    Translate the FGCM visit catalog to an fgcm-compatible object

    Parameters
    ----------
    visitCat: `lsst.afw.table.BaseCatalog`
       FGCM visitCat from `lsst.fgcmcal.FgcmBuildStarsTask`

    Returns
    -------
    fgcmExpInfo: `numpy.ndarray`
       Numpy array for visit information for FGCM

    Notes
    -----
    After running this code, it is wise to `del visitCat` to clear the memory.
    """

    fgcmExpInfo = np.zeros(len(visitCat), dtype=[('VISIT', 'i8'),
                                                 ('MJD', 'f8'),
                                                 ('EXPTIME', 'f8'),
                                                 ('PSFSIGMA', 'f8'),
                                                 ('DELTA_APER', 'f8'),
                                                 ('SKYBACKGROUND', 'f8'),
                                                 ('DEEPFLAG', 'i2'),
                                                 ('TELHA', 'f8'),
                                                 ('TELRA', 'f8'),
                                                 ('TELDEC', 'f8'),
                                                 ('TELROT', 'f8'),
                                                 ('PMB', 'f8'),
                                                 ('FILTERNAME', 'a50')])
    fgcmExpInfo['VISIT'][:] = visitCat['visit']
    fgcmExpInfo['MJD'][:] = visitCat['mjd']
    fgcmExpInfo['EXPTIME'][:] = visitCat['exptime']
    fgcmExpInfo['DEEPFLAG'][:] = visitCat['deepFlag']
    fgcmExpInfo['TELHA'][:] = visitCat['telha']
    fgcmExpInfo['TELRA'][:] = visitCat['telra']
    fgcmExpInfo['TELDEC'][:] = visitCat['teldec']
    fgcmExpInfo['TELROT'][:] = visitCat['telrot']
    fgcmExpInfo['PMB'][:] = visitCat['pmb']
    fgcmExpInfo['PSFSIGMA'][:] = visitCat['psfSigma']
    fgcmExpInfo['DELTA_APER'][:] = visitCat['deltaAper']
    fgcmExpInfo['SKYBACKGROUND'][:] = visitCat['skyBackground']
    # Note that we have to go through asAstropy() to get a string
    #  array out of an afwTable.  This is faster than a row-by-row loop.
    fgcmExpInfo['FILTERNAME'][:] = visitCat.asAstropy()['physicalFilter']

    return fgcmExpInfo


def computeReferencePixelScale(camera):
    """
    Compute the median pixel scale in the camera

    Returns
    -------
    pixelScale: `float`
       Average pixel scale (arcsecond) over the camera
    """

    boresight = geom.SpherePoint(180.0*geom.degrees, 0.0*geom.degrees)
    orientation = 0.0*geom.degrees
    flipX = False

    # Create a temporary visitInfo for input to createInitialSkyWcs
    visitInfo = afwImage.VisitInfo(boresightRaDec=boresight,
                                   boresightRotAngle=orientation,
                                   rotType=afwImage.RotType.SKY)

    pixelScales = np.zeros(len(camera))
    for i, detector in enumerate(camera):
        wcs = createInitialSkyWcs(visitInfo, detector, flipX)
        pixelScales[i] = wcs.getPixelScale(detector.getBBox().getCenter()).asArcseconds()

    ok, = np.where(pixelScales > 0.0)
    return np.median(pixelScales[ok])


def computeApproxPixelAreaFields(camera):
    """
    Compute the approximate pixel area bounded fields from the camera
    geometry.

    Parameters
    ----------
    camera: `lsst.afw.cameraGeom.Camera`

    Returns
    -------
    approxPixelAreaFields: `dict`
       Dictionary of approximate area fields, keyed with detector ID
    """

    areaScaling = 1. / computeReferencePixelScale(camera)**2.

    # Generate fake WCSs centered at 180/0 to avoid the RA=0/360 problem,
    # since we are looking for relative scales
    boresight = geom.SpherePoint(180.0*geom.degrees, 0.0*geom.degrees)

    flipX = False
    # Create a temporary visitInfo for input to createInitialSkyWcs
    # The orientation does not matter for the area computation
    visitInfo = afwImage.VisitInfo(boresightRaDec=boresight,
                                   boresightRotAngle=0.0*geom.degrees,
                                   rotType=afwImage.RotType.SKY)

    approxPixelAreaFields = {}

    for i, detector in enumerate(camera):
        key = detector.getId()

        wcs = createInitialSkyWcs(visitInfo, detector, flipX)
        bbox = detector.getBBox()

        areaField = afwMath.PixelAreaBoundedField(bbox, wcs,
                                                  unit=geom.arcseconds, scaling=areaScaling)
        approxAreaField = afwMath.ChebyshevBoundedField.approximate(areaField)

        approxPixelAreaFields[key] = approxAreaField

    return approxPixelAreaFields


def makeZptSchema(superStarChebyshevSize, zptChebyshevSize):
    """
    Make the zeropoint schema

    Parameters
    ----------
    superStarChebyshevSize: `int`
       Length of the superstar chebyshev array
    zptChebyshevSize: `int`
       Length of the zeropoint chebyshev array

    Returns
    -------
    zptSchema: `lsst.afw.table.schema`
    """

    zptSchema = afwTable.Schema()

    zptSchema.addField('visit', type=np.int64, doc='Visit number')
    zptSchema.addField('detector', type=np.int32, doc='Detector ID number')
    zptSchema.addField('fgcmFlag', type=np.int32, doc=('FGCM flag value: '
                                                       '1: Photometric, used in fit; '
                                                       '2: Photometric, not used in fit; '
                                                       '4: Non-photometric, on partly photometric night; '
                                                       '8: Non-photometric, on non-photometric night; '
                                                       '16: No zeropoint could be determined; '
                                                       '32: Too few stars for reliable gray computation'))
    zptSchema.addField('fgcmZpt', type=np.float64, doc='FGCM zeropoint (center of CCD)')
    zptSchema.addField('fgcmZptErr', type=np.float64,
                       doc='Error on zeropoint, estimated from repeatability + number of obs')
    zptSchema.addField('fgcmfZptChebXyMax', type='ArrayD', size=2,
                       doc='maximum x/maximum y to scale to apply chebyshev parameters')
    zptSchema.addField('fgcmfZptCheb', type='ArrayD',
                       size=zptChebyshevSize,
                       doc='Chebyshev parameters (flattened) for zeropoint')
    zptSchema.addField('fgcmfZptSstarCheb', type='ArrayD',
                       size=superStarChebyshevSize,
                       doc='Chebyshev parameters (flattened) for superStarFlat')
    zptSchema.addField('fgcmI0', type=np.float64, doc='Integral of the passband')
    zptSchema.addField('fgcmI10', type=np.float64, doc='Normalized chromatic integral')
    zptSchema.addField('fgcmR0', type=np.float64,
                       doc='Retrieved i0 integral, estimated from stars (only for flag 1)')
    zptSchema.addField('fgcmR10', type=np.float64,
                       doc='Retrieved i10 integral, estimated from stars (only for flag 1)')
    zptSchema.addField('fgcmGry', type=np.float64,
                       doc='Estimated gray extinction relative to atmospheric solution; '
                       'only for fgcmFlag <= 4 (see fgcmFlag) ')
    zptSchema.addField('fgcmDeltaChrom', type=np.float64,
                       doc='Mean chromatic correction for stars in this ccd; '
                       'only for fgcmFlag <= 4 (see fgcmFlag)')
    zptSchema.addField('fgcmZptVar', type=np.float64, doc='Variance of zeropoint over ccd')
    zptSchema.addField('fgcmTilings', type=np.float64,
                       doc='Number of photometric tilings used for solution for ccd')
    zptSchema.addField('fgcmFpGry', type=np.float64,
                       doc='Average gray extinction over the full focal plane '
                       '(same for all ccds in a visit)')
    zptSchema.addField('fgcmFpGryBlue', type=np.float64,
                       doc='Average gray extinction over the full focal plane '
                       'for 25% bluest stars')
    zptSchema.addField('fgcmFpGryBlueErr', type=np.float64,
                       doc='Error on Average gray extinction over the full focal plane '
                       'for 25% bluest stars')
    zptSchema.addField('fgcmFpGryRed', type=np.float64,
                       doc='Average gray extinction over the full focal plane '
                       'for 25% reddest stars')
    zptSchema.addField('fgcmFpGryRedErr', type=np.float64,
                       doc='Error on Average gray extinction over the full focal plane '
                       'for 25% reddest stars')
    zptSchema.addField('fgcmFpVar', type=np.float64,
                       doc='Variance of gray extinction over the full focal plane '
                       '(same for all ccds in a visit)')
    zptSchema.addField('fgcmDust', type=np.float64,
                       doc='Gray dust extinction from the primary/corrector'
                       'at the time of the exposure')
    zptSchema.addField('fgcmFlat', type=np.float64, doc='Superstarflat illumination correction')
    zptSchema.addField('fgcmAperCorr', type=np.float64, doc='Aperture correction estimated by fgcm')
    zptSchema.addField('fgcmDeltaMagBkg', type=np.float64,
                       doc=('Local background correction from brightest percentile '
                            '(value set by deltaMagBkgOffsetPercentile) calibration '
                            'stars.'))
    zptSchema.addField('exptime', type=np.float32, doc='Exposure time')
    zptSchema.addField('filtername', type=str, size=30, doc='Filter name')

    return zptSchema


def makeZptCat(zptSchema, zpStruct):
    """
    Make the zeropoint catalog for persistence

    Parameters
    ----------
    zptSchema: `lsst.afw.table.Schema`
       Zeropoint catalog schema
    zpStruct: `numpy.ndarray`
       Zeropoint structure from fgcm

    Returns
    -------
    zptCat: `afwTable.BaseCatalog`
       Zeropoint catalog for persistence
    """

    zptCat = afwTable.BaseCatalog(zptSchema)
    zptCat.reserve(zpStruct.size)

    for filterName in zpStruct['FILTERNAME']:
        rec = zptCat.addNew()
        rec['filtername'] = filterName.decode('utf-8')

    zptCat['visit'][:] = zpStruct[FGCM_EXP_FIELD]
    zptCat['detector'][:] = zpStruct[FGCM_CCD_FIELD]
    zptCat['fgcmFlag'][:] = zpStruct['FGCM_FLAG']
    zptCat['fgcmZpt'][:] = zpStruct['FGCM_ZPT']
    zptCat['fgcmZptErr'][:] = zpStruct['FGCM_ZPTERR']
    zptCat['fgcmfZptChebXyMax'][:, :] = zpStruct['FGCM_FZPT_XYMAX']
    zptCat['fgcmfZptCheb'][:, :] = zpStruct['FGCM_FZPT_CHEB']
    zptCat['fgcmfZptSstarCheb'][:, :] = zpStruct['FGCM_FZPT_SSTAR_CHEB']
    zptCat['fgcmI0'][:] = zpStruct['FGCM_I0']
    zptCat['fgcmI10'][:] = zpStruct['FGCM_I10']
    zptCat['fgcmR0'][:] = zpStruct['FGCM_R0']
    zptCat['fgcmR10'][:] = zpStruct['FGCM_R10']
    zptCat['fgcmGry'][:] = zpStruct['FGCM_GRY']
    zptCat['fgcmDeltaChrom'][:] = zpStruct['FGCM_DELTACHROM']
    zptCat['fgcmZptVar'][:] = zpStruct['FGCM_ZPTVAR']
    zptCat['fgcmTilings'][:] = zpStruct['FGCM_TILINGS']
    zptCat['fgcmFpGry'][:] = zpStruct['FGCM_FPGRY']
    zptCat['fgcmFpGryBlue'][:] = zpStruct['FGCM_FPGRY_CSPLIT'][:, 0]
    zptCat['fgcmFpGryBlueErr'][:] = zpStruct['FGCM_FPGRY_CSPLITERR'][:, 0]
    zptCat['fgcmFpGryRed'][:] = zpStruct['FGCM_FPGRY_CSPLIT'][:, 2]
    zptCat['fgcmFpGryRedErr'][:] = zpStruct['FGCM_FPGRY_CSPLITERR'][:, 2]
    zptCat['fgcmFpVar'][:] = zpStruct['FGCM_FPVAR']
    zptCat['fgcmDust'][:] = zpStruct['FGCM_DUST']
    zptCat['fgcmFlat'][:] = zpStruct['FGCM_FLAT']
    zptCat['fgcmAperCorr'][:] = zpStruct['FGCM_APERCORR']
    zptCat['fgcmDeltaMagBkg'][:] = zpStruct['FGCM_DELTAMAGBKG']
    zptCat['exptime'][:] = zpStruct['EXPTIME']

    return zptCat


def makeAtmSchema():
    """
    Make the atmosphere schema

    Returns
    -------
    atmSchema: `lsst.afw.table.Schema`
    """

    atmSchema = afwTable.Schema()

    atmSchema.addField('visit', type=np.int64, doc='Visit number')
    atmSchema.addField('pmb', type=np.float64, doc='Barometric pressure (mb)')
    atmSchema.addField('pwv', type=np.float64, doc='Water vapor (mm)')
    atmSchema.addField('tau', type=np.float64, doc='Aerosol optical depth')
    atmSchema.addField('alpha', type=np.float64, doc='Aerosol slope')
    atmSchema.addField('o3', type=np.float64, doc='Ozone (dobson)')
    atmSchema.addField('secZenith', type=np.float64, doc='Secant(zenith) (~ airmass)')
    atmSchema.addField('cTrans', type=np.float64, doc='Transmission correction factor')
    atmSchema.addField('lamStd', type=np.float64, doc='Wavelength for transmission correction')

    return atmSchema


def makeAtmCat(atmSchema, atmStruct):
    """
    Make the atmosphere catalog for persistence

    Parameters
    ----------
    atmSchema: `lsst.afw.table.Schema`
       Atmosphere catalog schema
    atmStruct: `numpy.ndarray`
       Atmosphere structure from fgcm

    Returns
    -------
    atmCat: `lsst.afw.table.BaseCatalog`
       Atmosphere catalog for persistence
    """

    atmCat = afwTable.BaseCatalog(atmSchema)
    atmCat.resize(atmStruct.size)

    atmCat['visit'][:] = atmStruct['VISIT']
    atmCat['pmb'][:] = atmStruct['PMB']
    atmCat['pwv'][:] = atmStruct['PWV']
    atmCat['tau'][:] = atmStruct['TAU']
    atmCat['alpha'][:] = atmStruct['ALPHA']
    atmCat['o3'][:] = atmStruct['O3']
    atmCat['secZenith'][:] = atmStruct['SECZENITH']
    atmCat['cTrans'][:] = atmStruct['CTRANS']
    atmCat['lamStd'][:] = atmStruct['LAMSTD']

    return atmCat


def makeStdSchema(nBands):
    """
    Make the standard star schema

    Parameters
    ----------
    nBands: `int`
       Number of bands in standard star catalog

    Returns
    -------
    stdSchema: `lsst.afw.table.Schema`
    """

    stdSchema = afwTable.SimpleTable.makeMinimalSchema()
    stdSchema.addField('ngood', type='ArrayI', doc='Number of good observations',
                       size=nBands)
    stdSchema.addField('ntotal', type='ArrayI', doc='Number of total observations',
                       size=nBands)
    stdSchema.addField('mag_std_noabs', type='ArrayF',
                       doc='Standard magnitude (no absolute calibration)',
                       size=nBands)
    stdSchema.addField('magErr_std', type='ArrayF',
                       doc='Standard magnitude error',
                       size=nBands)
    stdSchema.addField('npsfcand', type='ArrayI',
                       doc='Number of observations flagged as psf candidates',
                       size=nBands)
    stdSchema.addField('delta_aper', type='ArrayF',
                       doc='Delta mag (small - large aperture)',
                       size=nBands)

    return stdSchema


def makeStdCat(stdSchema, stdStruct, goodBands):
    """
    Make the standard star catalog for persistence

    Parameters
    ----------
    stdSchema: `lsst.afw.table.Schema`
       Standard star catalog schema
    stdStruct: `numpy.ndarray`
       Standard star structure in FGCM format
    goodBands: `list`
       List of good band names used in stdStruct

    Returns
    -------
    stdCat: `lsst.afw.table.BaseCatalog`
       Standard star catalog for persistence
    """

    stdCat = afwTable.SimpleCatalog(stdSchema)
    stdCat.resize(stdStruct.size)

    stdCat['id'][:] = stdStruct['FGCM_ID']
    stdCat['coord_ra'][:] = stdStruct['RA'] * geom.degrees
    stdCat['coord_dec'][:] = stdStruct['DEC'] * geom.degrees
    stdCat['ngood'][:, :] = stdStruct['NGOOD'][:, :]
    stdCat['ntotal'][:, :] = stdStruct['NTOTAL'][:, :]
    stdCat['mag_std_noabs'][:, :] = stdStruct['MAG_STD'][:, :]
    stdCat['magErr_std'][:, :] = stdStruct['MAGERR_STD'][:, :]
    if 'NPSFCAND' in stdStruct.dtype.names:
        stdCat['npsfcand'][:, :] = stdStruct['NPSFCAND'][:, :]
    stdCat['delta_aper'][:, :] = stdStruct['DELTA_APER'][:, :]

    md = PropertyList()
    md.set("BANDS", list(goodBands))
    stdCat.setMetadata(md)

    return stdCat


def computeApertureRadiusFromName(fluxField):
    """
    Compute the radius associated with a CircularApertureFlux or ApFlux field.

    Parameters
    ----------
    fluxField : `str`
       CircularApertureFlux or ApFlux

    Returns
    -------
    apertureRadius : `float`
        Radius of the aperture field, in pixels.

    Raises
    ------
     RuntimeError: Raised if flux field is not a CircularApertureFlux,
       ApFlux, or apFlux.
    """
    # TODO: Move this method to more general stack method in DM-25775
    m = re.search(r'(CircularApertureFlux|ApFlux|apFlux)_(\d+)_(\d+)_', fluxField)

    if m is None:
        raise RuntimeError(f"Flux field {fluxField} does not correspond to a CircularApertureFlux or ApFlux")

    apertureRadius = float(m.groups()[1]) + float(m.groups()[2])/10.

    return apertureRadius


def extractReferenceMags(refStars, bands, filterMap):
    """
    Extract reference magnitudes from refStars for given bands and
    associated filterMap.

    Parameters
    ----------
    refStars : `astropy.table.Table` or `lsst.afw.table.BaseCatalog`
        FGCM reference star catalog.
    bands : `list`
        List of bands for calibration.
    filterMap: `dict`
        FGCM mapping of filter to band.

    Returns
    -------
    refMag : `np.ndarray`
        nstar x nband array of reference magnitudes.
    refMagErr : `np.ndarray`
        nstar x nband array of reference magnitude errors.
    """
    hasAstropyMeta = False
    try:
        meta = refStars.meta
        hasAstropyMeta = True
    except AttributeError:
        meta = refStars.getMetadata()

    if 'FILTERNAMES' in meta:
        if hasAstropyMeta:
            filternames = meta['FILTERNAMES']
        else:
            filternames = meta.getArray('FILTERNAMES')

        # The reference catalog that fgcm wants has one entry per band
        # in the config file
        refMag = np.zeros((len(refStars), len(bands)),
                          dtype=refStars['refMag'].dtype) + 99.0
        refMagErr = np.zeros_like(refMag) + 99.0
        for i, filtername in enumerate(filternames):
            # We are allowed to run the fit configured so that we do not
            # use every column in the reference catalog.
            try:
                band = filterMap[filtername]
            except KeyError:
                continue
            try:
                ind = bands.index(band)
            except ValueError:
                continue

            refMag[:, ind] = refStars['refMag'][:, i]
            refMagErr[:, ind] = refStars['refMagErr'][:, i]
    else:
        raise RuntimeError("FGCM reference stars missing FILTERNAMES metadata.")

    return refMag, refMagErr


def lookupStaticCalibrations(datasetType, registry, quantumDataId, collections):
    # For static calibrations, we search with a timespan that has unbounded
    # begin and end; we'll get an error if there's more than one match (because
    # then it's not static).
    timespan = Timespan(begin=None, end=None)
    result = []
    # First iterate over all of the data IDs for this dataset type that are
    # consistent with the quantum data ID.
    for dataId in registry.queryDataIds(datasetType.dimensions, dataId=quantumDataId):
        # Find the dataset with this data ID using the unbounded timespan.
        if ref := registry.findDataset(datasetType, dataId, collections=collections, timespan=timespan):
            result.append(ref)
    return result
