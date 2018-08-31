# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import os
import shutil
import numpy as np

import lsst.daf.persistence as dafPersistence
import lsst.afw.geom as afwGeom
import lsst.log
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask, LoadIndexedReferenceObjectsConfig
import lsst.afw.image as afwImage

import lsst.fgcmcal as fgcmcal


class FgcmcalTestBase(object):
    """
    Base class for fgcmcal tests, to genericize some test running and setup.

    Derive from this first, then from TestCase.
    """

    def setUp_base(self, inputDir=None, testDir=None, logLevel=None, otherArgs=[]):
        """
        Call from your child class's setUp() to get variables built.

        Parameters
        ----------
        None

        """

        self.config = None
        self.inputDir = inputDir
        self.testDir = testDir
        self.logLevel = logLevel
        self.otherArgs = otherArgs

        lsst.log.setLevel("daf.persistence.butler", lsst.log.FATAL)
        lsst.log.setLevel("CameraMapper", lsst.log.FATAL)

    def _runFgcmMakeLut(self, nBand, i0Std, i0Recon, i10Std, i10Recon):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])
        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmMakeLutTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        butler = dafPersistence.butler.Butler(self.testDir)
        tempTask = fgcmcal.FgcmFitCycleTask()
        fgcmLut, lutIndexVals, lutStd = tempTask._loadFgcmLut(butler)

        # Check that we got the requested number of bands...
        self.assertEqual(nBand, len(lutIndexVals[0]['FILTERNAMES']))

        self.assertFloatsAlmostEqual(i0Std, lutStd[0]['I0STD'], msg='I0Std', rtol=1e-5)
        self.assertFloatsAlmostEqual(i10Std, lutStd[0]['I10STD'], msg='I10Std', rtol=1e-5)

        indices = fgcmLut.getIndices(np.arange(nBand, dtype=np.int32),
                                     np.zeros(nBand) + lutStd[0]['PWVSTD'],
                                     np.zeros(nBand) + lutStd[0]['O3STD'],
                                     np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                                     np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                                     np.zeros(nBand) + 1. / np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                                     np.zeros(nBand, dtype=np.int32),
                                     np.zeros(nBand) + lutStd[0]['PMBSTD'])
        i0 = fgcmLut.computeI0(np.zeros(nBand) + lutStd[0]['PWVSTD'],
                               np.zeros(nBand) + lutStd[0]['O3STD'],
                               np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                               np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                               np.zeros(nBand) + 1. / np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                               np.zeros(nBand) + lutStd[0]['PMBSTD'],
                               indices)

        self.assertFloatsAlmostEqual(i0Recon, i0, msg='i0Recon', rtol=1e-5)

        i1 = fgcmLut.computeI1(np.zeros(nBand) + lutStd[0]['PWVSTD'],
                               np.zeros(nBand) + lutStd[0]['O3STD'],
                               np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                               np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                               np.zeros(nBand) + 1. / np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                               np.zeros(nBand) + lutStd[0]['PMBSTD'],
                               indices)

        self.assertFloatsAlmostEqual(i10Recon, i1 / i0, msg='i10Recon', rtol=1e-5)

    def _runFgcmBuildStars(self, nVisit, nStar, nObs):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])
        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmBuildStarsTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        butler = dafPersistence.butler.Butler(self.testDir)

        visitCat = butler.get('fgcmVisitCatalog')
        self.assertEqual(nVisit, len(visitCat))

        starIds = butler.get('fgcmStarIds')
        self.assertEqual(nStar, len(starIds))

        starObs = butler.get('fgcmStarObservations')
        self.assertEqual(nObs, len(starObs))

    def _runFgcmFitCycle(self, nZp, nGoodZp, nStdStars):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])

        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        # Move into the test directory so the plots will get cleaned in tearDown
        cwd = os.getcwd()
        os.chdir(self.testDir)

        result = fgcmcal.FgcmFitCycleTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        # Move back to the previous directory
        os.chdir(cwd)

        butler = dafPersistence.butler.Butler(self.testDir)

        zps = butler.get('fgcmZeropoints', fgcmcycle=0)

        self.assertEqual(nZp, len(zps))

        gd, = np.where(zps['fgcmflag'] == 1)
        self.assertEqual(nGoodZp, len(gd))

        stds = butler.get('fgcmStandardStars', fgcmcycle=0)

        self.assertEqual(nStdStars, len(stds))

    def _runFgcmOutputProducts(self, visitDataRefName, ccdDataRefName, filterMapping,
                               zpOffsets, testVisit, testCcd, testFilter, testBandIndex):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])

        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmOutputProductsTask.parseAndRun(args=args, config=self.config,
                                                            doReturnResults=True)
        self._checkResult(result)

        # Extract the offsets from the results
        offsets = result.resultList[0].results.offsets

        self.assertFloatsAlmostEqual(offsets[0], zpOffsets[0], rtol=1e-6)
        self.assertFloatsAlmostEqual(offsets[1], zpOffsets[1], rtol=1e-6)

        butler = dafPersistence.butler.Butler(self.testDir)

        # Test the reference catalog stars

        # Read in the raw stars...
        rawStars = butler.get('fgcmStandardStars', fgcmcycle=0)

        # Read in the new reference catalog...
        config = LoadIndexedReferenceObjectsConfig()
        config.ref_dataset_name = 'fgcm_stars'
        task = LoadIndexedReferenceObjectsTask(butler, config=config)
        # Read in a giant radius to get them all
        refStruct = task.loadSkyCircle(rawStars[0].getCoord(), 5.0 * lsst.geom.degrees,
                                       filterName='r')

        # Make sure all the stars are there
        self.assertEqual(len(rawStars), len(refStruct.refCat))

        # And make sure the numbers are consistent
        test, = np.where(rawStars['id'][0] == refStruct.refCat['id'])

        mag = rawStars['mag_std_noabs'][0, 0] + offsets[0]
        flux = afwImage.fluxFromABMag(mag)
        fluxErr = afwImage.fluxErrFromABMagErr(rawStars['magerr_std'][0, 0], mag)
        self.assertFloatsAlmostEqual(flux, refStruct.refCat['r_flux'][test[0]], rtol=1e-6)
        self.assertFloatsAlmostEqual(fluxErr, refStruct.refCat['r_fluxErr'][test[0]], rtol=1e-6)

        # Test the joincal_photoCalib output

        zptCat = butler.get('fgcmZeropoints', fgcmcycle=0)
        selected = (zptCat['fgcmflag'] < 16)

        # Read in all the calibrations, these should all be there
        for rec in zptCat[selected]:
            testCal = butler.get('jointcal_photoCalib',
                                 dataId={visitDataRefName: int(rec['visit']),
                                         ccdDataRefName: int(rec['ccd']),
                                         'filter': filterMapping[rec['filtername']],
                                         'tract': 0})

        # Our round-trip tests will be on this final one which is still loaded
        testCal = butler.get('jointcal_photoCalib',
                             dataId={visitDataRefName: int(testVisit),
                                     ccdDataRefName: int(testCcd),
                                     'filter': filterMapping[testFilter],
                                     'tract': 0})

        src = butler.get('src', dataId={visitDataRefName: int(testVisit),
                                        ccdDataRefName: int(testCcd)})

        # Only test sources with positive flux
        gdSrc = (src['slot_CalibFlux_flux'] > 0.0)

        # We need to apply the calibration offset to the fgcmzpt (which is internal
        # and doesn't know about that yet)
        testZpInd, = np.where((zptCat['visit'] == testVisit) &
                              (zptCat['ccd'] == testCcd))
        fgcmZpt = zptCat['fgcmzpt'][testZpInd] + offsets[testBandIndex]

        # This is the magnitude through the mean calibration
        photoCalMeanCalMags = np.zeros(gdSrc.sum())
        # This is the magnitude through the full focal-plane variable mags
        photoCalMags = np.zeros_like(photoCalMeanCalMags)
        # This is the magnitude with the FGCM (central-ccd) zeropoint
        zptMeanCalMags = np.zeros_like(photoCalMeanCalMags)

        for i, rec in enumerate(src[gdSrc]):
            photoCalMeanCalMags[i] = testCal.instFluxToMagnitude(rec['slot_CalibFlux_flux'])
            photoCalMags[i] = testCal.instFluxToMagnitude(rec['slot_CalibFlux_flux'],
                                                          rec.getCentroid())
            zptMeanCalMags[i] = fgcmZpt - 2.5*np.log10(rec['slot_CalibFlux_flux'])

        # These should be very close but some tiny differences because the fgcm value
        # is defined at the center of the bbox, and the photoCal is the mean over the box
        self.assertFloatsAlmostEqual(photoCalMeanCalMags,
                                     zptMeanCalMags, rtol=1e-6)
        # These should be roughly equal, but not precisely because of the focal-plane
        # variation.  However, this is a useful sanity check for something going totally
        # wrong.
        self.assertFloatsAlmostEqual(photoCalMeanCalMags,
                                     photoCalMags, rtol=1e-2)

        # Test the transmission output

        visitCatalog = butler.get('fgcmVisitCatalog')
        lutCat = butler.get('fgcmLookUpTable')

        testTrans = butler.get('transmission_atmosphere_fgcm',
                               dataId={visitDataRefName: visitCatalog[0]['visit']})
        testResp = testTrans.sampleAt(position=afwGeom.Point2D(0, 0),
                                      wavelengths=lutCat[0]['atmlambda'])

        # The fit to be roughly consistent with the standard, although the
        # airmass is taken into account even with the "frozen" atmosphere.
        # This is also a rough comparison, because the interpolation does
        # not work well with such a coarse look-up table used for the test.
        self.assertFloatsAlmostEqual(testResp, lutCat[0]['atmstdtrans'], atol=0.06)

        # The second should be close to the first, but there is the airmass
        # difference so they aren't identical
        testTrans2 = butler.get('transmission_atmosphere_fgcm',
                                dataId={visitDataRefName: visitCatalog[1]['visit']})
        testResp2 = testTrans2.sampleAt(position=afwGeom.Point2D(0, 0),
                                        wavelengths=lutCat[0]['atmlambda'])
        self.assertFloatsAlmostEqual(testResp, testResp2, atol=1e-4)

    def _checkResult(self, result):
        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.assertEqual(result.resultList[0].exitStatus, 0)

    def tearDown(self):
        if getattr(self, 'config', None) is not None:
            del self.config

        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir, True)
