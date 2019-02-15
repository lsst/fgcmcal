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
"""General fgcmcal testing class.

This class is used as the basis for individual obs package tests using
data from testdata_jointcal.
"""

import os
import shutil
import numpy as np
import glob

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
        inputDir: `str`, optional
           Input directory
        testDir: `str`, optional
           Test directory
        logLevel: `str`, optional
           Override loglevel for command-line tasks
        otherArgs: `list`, default=[]
           List of additional arguments to send to command-line tasks
        """

        self.config = None
        self.inputDir = inputDir
        self.testDir = testDir
        self.logLevel = logLevel
        self.otherArgs = otherArgs

        lsst.log.setLevel("daf.persistence.butler", lsst.log.FATAL)
        lsst.log.setLevel("CameraMapper", lsst.log.FATAL)

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])

    def _testFgcmMakeLut(self, nBand, i0Std, i0Recon, i10Std, i10Recon):
        """
        Test running of FgcmMakeLutTask

        Parameters
        ----------
        nBand: `int`
           Number of bands tested
        i0Std: `np.array', size nBand
           Values of i0Std to compare to
        i10Std: `np.array`, size nBand
           Values of i10Std to compare to
        i0Recon: `np.array`, size nBand
           Values of reconstructed i0 to compare to
        i10Recon: `np.array`, size nBand
           Values of reconsntructed i10 to compare to

        Raises
        ------
        Exceptions on test failures
        """

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
                                     np.zeros(nBand) + np.log(lutStd[0]['PWVSTD']),
                                     np.zeros(nBand) + lutStd[0]['O3STD'],
                                     np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                                     np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                                     np.zeros(nBand) + 1. / np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                                     np.zeros(nBand, dtype=np.int32),
                                     np.zeros(nBand) + lutStd[0]['PMBSTD'])
        i0 = fgcmLut.computeI0(np.zeros(nBand) + np.log(lutStd[0]['PWVSTD']),
                               np.zeros(nBand) + lutStd[0]['O3STD'],
                               np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                               np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                               np.zeros(nBand) + 1. / np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                               np.zeros(nBand) + lutStd[0]['PMBSTD'],
                               indices)

        self.assertFloatsAlmostEqual(i0Recon, i0, msg='i0Recon', rtol=1e-5)

        i1 = fgcmLut.computeI1(np.zeros(nBand) + np.log(lutStd[0]['PWVSTD']),
                               np.zeros(nBand) + lutStd[0]['O3STD'],
                               np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                               np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                               np.zeros(nBand) + 1. / np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                               np.zeros(nBand) + lutStd[0]['PMBSTD'],
                               indices)

        self.assertFloatsAlmostEqual(i10Recon, i1 / i0, msg='i10Recon', rtol=1e-5)

    def _testFgcmBuildStars(self, nVisit, nStar, nObs):
        """
        Test running of FgcmBuildStarsTask

        Parameters
        ----------
        nVisit: `int`
           Number of visits expected
        nStar: `int`
           Number of stars expected
        nObs: `int`
           Number of observations of stars expected

        Raises
        ------
        Exceptions on test failures
        """

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

    def _testFgcmFitCycle(self, nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots, skipChecks=False):
        """
        Test running of FgcmFitCycleTask

        Parameters
        ----------
        nZp: `int`
           Number of zeropoints created by the task
        nGoodZp: `int`
           Number of good (photometric) zeropoints created
        nOkZp: `int`
           Number of constrained zeropoints (photometric or not)
        nBadZp: `int`
           Number of unconstrained (bad) zeropoints
        nStdStars: `int`
           Number of standard stars produced
        nPlots: `int`
           Number of plots produced
        skipChecks: `bool`, optional
           Skip number checks, when running less-than-final cycle.
           Default is False.
        """

        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        # Move into the test directory so the plots will get cleaned in tearDown
        # In the future, with Gen3, we will probably have a better way of managing
        # non-data output such as plots.
        cwd = os.getcwd()
        os.chdir(self.testDir)

        result = fgcmcal.FgcmFitCycleTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        # Move back to the previous directory
        os.chdir(cwd)

        if skipChecks:
            return

        # Check that the expected number of plots are there.
        plots = glob.glob(os.path.join(self.testDir, self.config.outfileBase +
                                       '_cycle%02d_plots/' % (self.config.cycleNumber) +
                                       '*.png'))
        self.assertEqual(nPlots, len(plots))

        butler = dafPersistence.butler.Butler(self.testDir)

        zps = butler.get('fgcmZeropoints', fgcmcycle=self.config.cycleNumber)

        # Check the numbers of zeropoints in all, good, okay, and bad
        self.assertEqual(nZp, len(zps))

        gd, = np.where(zps['fgcmFlag'] == 1)
        self.assertEqual(nGoodZp, len(gd))

        ok, = np.where(zps['fgcmFlag'] < 16)
        self.assertEqual(nOkZp, len(ok))

        bd, = np.where(zps['fgcmFlag'] >= 16)
        self.assertEqual(nBadZp, len(bd))

        # Check that there are no illegal values with the ok zeropoints
        test, = np.where(zps['fgcmZpt'][gd] < -9000.0)
        self.assertEqual(0, len(test))

        stds = butler.get('fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

        self.assertEqual(nStdStars, len(stds))

    def _testFgcmOutputProducts(self, visitDataRefName, ccdDataRefName, filterMapping,
                                zpOffsets, testVisit, testCcd, testFilter, testBandIndex):
        """
        Test running of FgcmOutputProductsTask

        Parameters
        ----------
        visitDataRefName: `str`
           Name of column in dataRef to get the visit
        ccdDataRefName: `str`
           Name of column in dataRef to get the ccd
        filterMapping: `dict`
           Mapping of filterName to dataRef filter names
        zpOffsets: `np.array`
           Zeropoint offsets expected
        testVisit: `int`
           Visit id to check for round-trip computations
        testCcd: `int`
           Ccd id to check for round-trip computations
        testFilter: `str`
           Filtername for testVisit/testCcd
        testBandIndex: `int`
           Band index for testVisit/testCcd
        """

        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmOutputProductsTask.parseAndRun(args=args, config=self.config,
                                                            doReturnResults=True)
        self._checkResult(result)

        # Extract the offsets from the results
        offsets = result.resultList[0].results.offsets

        self.assertFloatsAlmostEqual(offsets[0], zpOffsets[0], atol=1e-7)
        self.assertFloatsAlmostEqual(offsets[1], zpOffsets[1], atol=1e-7)

        butler = dafPersistence.butler.Butler(self.testDir)

        # Test the reference catalog stars

        # Read in the raw stars...
        rawStars = butler.get('fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

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
        fluxErr = afwImage.fluxErrFromABMagErr(rawStars['magErr_std'][0, 0], mag)
        self.assertFloatsAlmostEqual(flux, refStruct.refCat['r_flux'][test[0]], rtol=1e-6)
        self.assertFloatsAlmostEqual(fluxErr, refStruct.refCat['r_fluxErr'][test[0]], rtol=1e-6)

        # Test the joincal_photoCalib output

        zptCat = butler.get('fgcmZeropoints', fgcmcycle=self.config.cycleNumber)
        selected = (zptCat['fgcmFlag'] < 16)

        # Read in all the calibrations, these should all be there
        # This test is simply to ensure that all the photoCalib files exist
        for rec in zptCat[selected]:
            testCal = butler.get('fgcm_photoCalib',
                                 dataId={visitDataRefName: int(rec['visit']),
                                         ccdDataRefName: int(rec['ccd']),
                                         'filter': filterMapping[rec['filtername']]})
            self.assertIsNotNone(testCal)

        # We do round-trip value checking on just the final one (chosen arbitrarily)
        testCal = butler.get('fgcm_photoCalib',
                             dataId={visitDataRefName: int(testVisit),
                                     ccdDataRefName: int(testCcd),
                                     'filter': filterMapping[testFilter]})
        self.assertIsNotNone(testCal)

        src = butler.get('src', dataId={visitDataRefName: int(testVisit),
                                        ccdDataRefName: int(testCcd)})

        # Only test sources with positive flux
        gdSrc = (src['slot_CalibFlux_flux'] > 0.0)

        # We need to apply the calibration offset to the fgcmzpt (which is internal
        # and doesn't know about that yet)
        testZpInd, = np.where((zptCat['visit'] == testVisit) &
                              (zptCat['ccd'] == testCcd))
        fgcmZpt = zptCat['fgcmZpt'][testZpInd] + offsets[testBandIndex]

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
                                      wavelengths=lutCat[0]['atmLambda'])

        # The test fit is performed with the atmosphere parameters frozen
        # (freezeStdAtmosphere = True).  Thus the only difference between
        # these output atmospheres and the standard is the different
        # airmass.  Furthermore, this is a very rough comparison because
        # the look-up table is computed with very coarse sampling for faster
        # testing.
        # Therefore, this rough comparison can only be seen as a sanity check
        # and is not high precision.
        self.assertFloatsAlmostEqual(testResp, lutCat[0]['atmStdTrans'], atol=0.06)

        # The second should be close to the first, but there is the airmass
        # difference so they aren't identical
        testTrans2 = butler.get('transmission_atmosphere_fgcm',
                                dataId={visitDataRefName: visitCatalog[1]['visit']})
        testResp2 = testTrans2.sampleAt(position=afwGeom.Point2D(0, 0),
                                        wavelengths=lutCat[0]['atmLambda'])
        self.assertFloatsAlmostEqual(testResp, testResp2, atol=1e-4)

    def _checkResult(self, result):
        """
        Check the result output from the task

        Parameters
        ----------
        result: `pipeBase.struct`
           Result structure output from a task

        Raises
        ------
        Exceptions on test failures
        """

        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.assertEqual(result.resultList[0].exitStatus, 0)

    def tearDown(self):
        """
        Tear down and clear directories
        """

        if getattr(self, 'config', None) is not None:
            del self.config

        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir, True)
