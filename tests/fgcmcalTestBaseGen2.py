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
data from testdata_jointcal for Gen2 repos.
"""

import os
import shutil
import numpy as np
import numpy.testing as testing
import glob
import esutil

import lsst.daf.persistence as dafPersist
import lsst.geom as geom
import lsst.log
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask, LoadIndexedReferenceObjectsConfig
from astropy import units

import lsst.fgcmcal as fgcmcal


class FgcmcalTestBaseGen2(object):
    """
    Base class for gen2 fgcmcal tests, to genericize some test running and setup.

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

        self.inputDir = inputDir
        self.testDir = testDir
        self.logLevel = logLevel
        self.otherArgs = otherArgs

        self.config = None
        self.configfiles = []

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
        """

        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        if len(self.configfiles) > 0:
            args.extend(['--configfile', *self.configfiles])
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmMakeLutTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        butler = dafPersist.butler.Butler(self.testDir)
        tempTask = fgcmcal.FgcmFitCycleTask()
        lutCat = butler.get('fgcmLookUpTable')
        fgcmLut, lutIndexVals, lutStd = fgcmcal.utilities.translateFgcmLut(lutCat,
                                                                           dict(tempTask.config.filterMap))

        # Check that we got the requested number of bands...
        self.assertEqual(nBand, len(lutIndexVals[0]['FILTERNAMES']))

        self.assertFloatsAlmostEqual(i0Std, lutStd[0]['I0STD'], msg='I0Std', rtol=1e-5)
        self.assertFloatsAlmostEqual(i10Std, lutStd[0]['I10STD'], msg='I10Std', rtol=1e-5)

        indices = fgcmLut.getIndices(np.arange(nBand, dtype=np.int32),
                                     np.zeros(nBand) + np.log(lutStd[0]['PWVSTD']),
                                     np.zeros(nBand) + lutStd[0]['O3STD'],
                                     np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                                     np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                                     np.zeros(nBand) + 1./np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                                     np.zeros(nBand, dtype=np.int32),
                                     np.zeros(nBand) + lutStd[0]['PMBSTD'])
        i0 = fgcmLut.computeI0(np.zeros(nBand) + np.log(lutStd[0]['PWVSTD']),
                               np.zeros(nBand) + lutStd[0]['O3STD'],
                               np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                               np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                               np.zeros(nBand) + 1./np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                               np.zeros(nBand) + lutStd[0]['PMBSTD'],
                               indices)

        self.assertFloatsAlmostEqual(i0Recon, i0, msg='i0Recon', rtol=1e-5)

        i1 = fgcmLut.computeI1(np.zeros(nBand) + np.log(lutStd[0]['PWVSTD']),
                               np.zeros(nBand) + lutStd[0]['O3STD'],
                               np.zeros(nBand) + np.log(lutStd[0]['TAUSTD']),
                               np.zeros(nBand) + lutStd[0]['ALPHASTD'],
                               np.zeros(nBand) + 1./np.cos(np.radians(lutStd[0]['ZENITHSTD'])),
                               np.zeros(nBand) + lutStd[0]['PMBSTD'],
                               indices)

        self.assertFloatsAlmostEqual(i10Recon, i1/i0, msg='i10Recon', rtol=1e-5)

    def _testFgcmBuildStarsTable(self, visits, nStar, nObs):
        """
        Test running of FgcmBuildStarsTableTask

        Parameters
        ----------
        visits: `list`
           List of visits to calibrate
        nStar: `int`
           Number of stars expected
        nObs: `int`
           Number of observations of stars expected
        """

        args = [self.inputDir, '--output', self.testDir,
                '--id', 'visit='+'^'.join([str(visit) for visit in visits]),
                '--doraise']
        if len(self.configfiles) > 0:
            args.extend(['--configfile', *self.configfiles])
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmBuildStarsTableTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        butler = dafPersist.butler.Butler(self.testDir)

        visitCat = butler.get('fgcmVisitCatalog')
        self.assertEqual(len(visits), len(visitCat))

        starIds = butler.get('fgcmStarIds')
        self.assertEqual(nStar, len(starIds))

        starObs = butler.get('fgcmStarObservations')
        self.assertEqual(nObs, len(starObs))

    def _testFgcmBuildStarsAndCompare(self, visits):
        """
        Test running of FgcmBuildStarsTask and compare to Table run

        Parameters
        ----------
        visits: `list`
           List of visits to calibrate
        """
        args = [self.testDir, '--output', os.path.join(self.testDir, 'rerun', 'src'),
                '--id', 'visit='+'^'.join([str(visit) for visit in visits]),
                '--doraise']
        if len(self.configfiles) > 0:
            args.extend(['--configfile', *self.configfiles])
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmBuildStarsTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        butlerSrc = dafPersist.Butler(os.path.join(self.testDir, 'rerun', 'src'))
        butlerTable = dafPersist.Butler(os.path.join(self.testDir))

        # We compare the two catalogs to ensure they contain the same data.  They will
        # not be identical in ordering because the input data was ingested in a different
        # order (hence the stars are rearranged).
        self._compareBuildStars(butlerSrc, butlerTable)

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
        if len(self.configfiles) > 0:
            args.extend(['--configfile', *self.configfiles])
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
        self.assertEqual(len(plots), nPlots)

        butler = dafPersist.butler.Butler(self.testDir)

        zps = butler.get('fgcmZeropoints', fgcmcycle=self.config.cycleNumber)

        # Check the numbers of zeropoints in all, good, okay, and bad
        self.assertEqual(len(zps), nZp)

        gd, = np.where(zps['fgcmFlag'] == 1)
        self.assertEqual(len(gd), nGoodZp)

        ok, = np.where(zps['fgcmFlag'] < 16)
        self.assertEqual(len(ok), nOkZp)

        bd, = np.where(zps['fgcmFlag'] >= 16)
        self.assertEqual(len(bd), nBadZp)

        # Check that there are no illegal values with the ok zeropoints
        test, = np.where(zps['fgcmZpt'][gd] < -9000.0)
        self.assertEqual(len(test), 0)

        stds = butler.get('fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

        self.assertEqual(len(stds), nStdStars)

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
        if len(self.configfiles) > 0:
            args.extend(['--configfile', *self.configfiles])
        args.extend(self.otherArgs)

        result = fgcmcal.FgcmOutputProductsTask.parseAndRun(args=args, config=self.config,
                                                            doReturnResults=True)
        self._checkResult(result)

        # Extract the offsets from the results
        offsets = result.resultList[0].results.offsets

        self.assertFloatsAlmostEqual(offsets, zpOffsets, atol=1e-6)

        butler = dafPersist.butler.Butler(self.testDir)

        # Test the reference catalog stars

        # Read in the raw stars...
        rawStars = butler.get('fgcmStandardStars', fgcmcycle=self.config.cycleNumber)

        # Read in the new reference catalog...
        config = LoadIndexedReferenceObjectsConfig()
        config.ref_dataset_name = 'fgcm_stars'
        task = LoadIndexedReferenceObjectsTask(butler, config=config)

        # Read in a giant radius to get them all
        refStruct = task.loadSkyCircle(rawStars[0].getCoord(), 5.0*geom.degrees,
                                       filterName='r')

        # Make sure all the stars are there
        self.assertEqual(len(rawStars), len(refStruct.refCat))

        # And make sure the numbers are consistent
        test, = np.where(rawStars['id'][0] == refStruct.refCat['id'])

        # Perform math on numpy arrays to maintain datatypes
        mags = rawStars['mag_std_noabs'][:, 0].astype(np.float64) + offsets[0]
        fluxes = (mags*units.ABmag).to_value(units.nJy)
        fluxErrs = (np.log(10.)/2.5)*fluxes*rawStars['magErr_std'][:, 0].astype(np.float64)
        # Only check the first one
        self.assertFloatsAlmostEqual(fluxes[0], refStruct.refCat['r_flux'][test[0]])
        self.assertFloatsAlmostEqual(fluxErrs[0], refStruct.refCat['r_fluxErr'][test[0]])

        # Test the psf candidate counting, ratio should be between 0.0 and 1.0
        candRatio = (refStruct.refCat['r_nPsfCandidate'].astype(np.float64) /
                     refStruct.refCat['r_nTotal'].astype(np.float64))
        self.assertFloatsAlmostEqual(candRatio.min(), 0.0)
        self.assertFloatsAlmostEqual(candRatio.max(), 1.0)

        # Test the fgcm_photoCalib output

        zptCat = butler.get('fgcmZeropoints', fgcmcycle=self.config.cycleNumber)
        selected = (zptCat['fgcmFlag'] < 16)

        # Read in all the calibrations, these should all be there
        # This test is simply to ensure that all the photoCalib files exist
        for rec in zptCat[selected]:
            testCal = butler.get('fgcm_photoCalib',
                                 dataId={visitDataRefName: int(rec['visit']),
                                         ccdDataRefName: int(rec['detector']),
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
        gdSrc = (src['slot_CalibFlux_instFlux'] > 0.0)

        # We need to apply the calibration offset to the fgcmzpt (which is internal
        # and doesn't know about that yet)
        testZpInd, = np.where((zptCat['visit'] == testVisit) &
                              (zptCat['detector'] == testCcd))
        fgcmZpt = (zptCat['fgcmZpt'][testZpInd] + offsets[testBandIndex] +
                   zptCat['fgcmDeltaChrom'][testZpInd])
        fgcmZptGrayErr = np.sqrt(zptCat['fgcmZptVar'][testZpInd])

        if self.config.doComposeWcsJacobian:
            # The raw zeropoint needs to be modified to know about the wcs jacobian
            camera = butler.get('camera')
            approxPixelAreaFields = fgcmcal.utilities.computeApproxPixelAreaFields(camera)
            center = approxPixelAreaFields[testCcd].getBBox().getCenter()
            pixAreaCorr = approxPixelAreaFields[testCcd].evaluate(center)
            fgcmZpt += -2.5*np.log10(pixAreaCorr)

        # This is the magnitude through the mean calibration
        photoCalMeanCalMags = np.zeros(gdSrc.sum())
        # This is the magnitude through the full focal-plane variable mags
        photoCalMags = np.zeros_like(photoCalMeanCalMags)
        # This is the magnitude with the FGCM (central-ccd) zeropoint
        zptMeanCalMags = np.zeros_like(photoCalMeanCalMags)

        for i, rec in enumerate(src[gdSrc]):
            photoCalMeanCalMags[i] = testCal.instFluxToMagnitude(rec['slot_CalibFlux_instFlux'])
            photoCalMags[i] = testCal.instFluxToMagnitude(rec['slot_CalibFlux_instFlux'],
                                                          rec.getCentroid())
            zptMeanCalMags[i] = fgcmZpt - 2.5*np.log10(rec['slot_CalibFlux_instFlux'])

        # These should be very close but some tiny differences because the fgcm value
        # is defined at the center of the bbox, and the photoCal is the mean over the box
        self.assertFloatsAlmostEqual(photoCalMeanCalMags,
                                     zptMeanCalMags, rtol=1e-6)
        # These should be roughly equal, but not precisely because of the focal-plane
        # variation.  However, this is a useful sanity check for something going totally
        # wrong.
        self.assertFloatsAlmostEqual(photoCalMeanCalMags,
                                     photoCalMags, rtol=1e-2)

        # The next test compares the "FGCM standard magnitudes" (which are output
        # from the fgcm code itself) to the "calibrated magnitudes" that are
        # obtained from running photoCalib.calibrateCatalog() on the original
        # src catalogs.  This summary comparison ensures that using photoCalibs
        # yields the same results as what FGCM is computing internally.
        # Note that we additionally need to take into account the post-processing
        # offsets used in the tests.

        # For decent statistics, we are matching all the sources from one visit
        # (multiple ccds)

        subset = butler.subset('src', dataId={visitDataRefName: int(testVisit)})

        matchMag, matchDelta = self._getMatchedVisitCat(rawStars, subset, testBandIndex, offsets)

        st = np.argsort(matchMag)
        # Compare the brightest 25% of stars.  No matter the setting of
        # deltaMagBkgOffsetPercentile, we want to ensure that these stars
        # match on average.
        brightest, = np.where(matchMag < matchMag[st[int(0.25*st.size)]])
        self.assertFloatsAlmostEqual(np.median(matchDelta[brightest]), 0.0, atol=0.002)

        # And the photoCal error is just the zeropoint gray error
        self.assertFloatsAlmostEqual(testCal.getCalibrationErr(),
                                     (np.log(10.0)/2.5)*testCal.getCalibrationMean()*fgcmZptGrayErr)

        # Test the transmission output

        visitCatalog = butler.get('fgcmVisitCatalog')
        lutCat = butler.get('fgcmLookUpTable')

        testTrans = butler.get('transmission_atmosphere_fgcm',
                               dataId={visitDataRefName: visitCatalog[0]['visit']})
        testResp = testTrans.sampleAt(position=geom.Point2D(0, 0),
                                      wavelengths=lutCat[0]['atmLambda'])

        # The test fit is performed with the atmosphere parameters frozen
        # (freezeStdAtmosphere = True).  Thus the only difference between
        # these output atmospheres and the standard is the different
        # airmass.  Furthermore, this is a very rough comparison because
        # the look-up table is computed with very coarse sampling for faster
        # testing.

        # To account for overall throughput changes, we scale by the median ratio,
        # we only care about the shape
        ratio = np.median(testResp/lutCat[0]['atmStdTrans'])
        self.assertFloatsAlmostEqual(testResp/ratio, lutCat[0]['atmStdTrans'], atol=0.04)

        # The second should be close to the first, but there is the airmass
        # difference so they aren't identical.
        testTrans2 = butler.get('transmission_atmosphere_fgcm',
                                dataId={visitDataRefName: visitCatalog[1]['visit']})
        testResp2 = testTrans2.sampleAt(position=geom.Point2D(0, 0),
                                        wavelengths=lutCat[0]['atmLambda'])

        # As above, we scale by the ratio to compare the shape of the curve.
        ratio = np.median(testResp/testResp2)
        self.assertFloatsAlmostEqual(testResp/ratio, testResp2, atol=0.04)

    def _testFgcmCalibrateTract(self, visits, tract,
                                rawRepeatability, filterNCalibMap):
        """
        Test running of FgcmCalibrateTractTask

        Parameters
        ----------
        visits: `list`
           List of visits to calibrate
        tract: `int`
           Tract number
        rawRepeatability: `np.array`
           Expected raw repeatability after convergence.
           Length should be number of bands.
        filterNCalibMap: `dict`
           Mapping from filter name to number of photoCalibs created.
        """

        args = [self.inputDir, '--output', self.testDir,
                '--id', 'visit='+'^'.join([str(visit) for visit in visits]),
                'tract=%d' % (tract),
                '--doraise']
        if len(self.configfiles) > 0:
            args.extend(['--configfile', *self.configfiles])
        args.extend(self.otherArgs)

        # Move into the test directory so the plots will get cleaned in tearDown
        # In the future, with Gen3, we will probably have a better way of managing
        # non-data output such as plots.
        cwd = os.getcwd()
        os.chdir(self.testDir)

        result = fgcmcal.FgcmCalibrateTractTableTask.parseAndRun(args=args, config=self.config,
                                                                 doReturnResults=True)
        self._checkResult(result)

        # Move back to the previous directory
        os.chdir(cwd)

        # Check that the converged repeatability is what we expect
        repeatability = result.resultList[0].results.repeatability
        self.assertFloatsAlmostEqual(repeatability, rawRepeatability, atol=4e-6)

        butler = dafPersist.butler.Butler(self.testDir)

        # Check that the number of photoCalib objects in each filter are what we expect
        for filterName in filterNCalibMap.keys():
            subset = butler.subset('fgcm_tract_photoCalib', tract=tract, filter=filterName)
            tot = 0
            for dataRef in subset:
                if butler.datasetExists('fgcm_tract_photoCalib', dataId=dataRef.dataId):
                    tot += 1
            self.assertEqual(tot, filterNCalibMap[filterName])

        # Check that every visit got a transmission
        visits = butler.queryMetadata('fgcm_tract_photoCalib', ('visit'), tract=tract)
        for visit in visits:
            self.assertTrue(butler.datasetExists('transmission_atmosphere_fgcm_tract',
                                                 tract=tract, visit=visit))

        # Check that we got the reference catalog output.
        # This will raise an exception if the catalog is not there.
        config = LoadIndexedReferenceObjectsConfig()
        config.ref_dataset_name = 'fgcm_stars_%d' % (tract)
        task = LoadIndexedReferenceObjectsTask(butler, config=config)

        coord = geom.SpherePoint(337.656174*geom.degrees, 0.823595*geom.degrees)

        refStruct = task.loadSkyCircle(coord, 5.0*geom.degrees, filterName='r')

        # Test the psf candidate counting, ratio should be between 0.0 and 1.0
        candRatio = (refStruct.refCat['r_nPsfCandidate'].astype(np.float64) /
                     refStruct.refCat['r_nTotal'].astype(np.float64))
        self.assertFloatsAlmostEqual(candRatio.min(), 0.0)
        self.assertFloatsAlmostEqual(candRatio.max(), 1.0)

        # Test that temporary files aren't stored
        self.assertFalse(butler.datasetExists('fgcmVisitCatalog'))
        self.assertFalse(butler.datasetExists('fgcmStarObservations'))
        self.assertFalse(butler.datasetExists('fgcmStarIndices'))
        self.assertFalse(butler.datasetExists('fgcmReferenceStars'))

    def _compareBuildStars(self, butler1, butler2):
        """
        Compare the full set of BuildStars outputs with files from two
        repos.

        Parameters
        ----------
        butler1, butler2 : `lsst.daf.persistence.Butler`
        """
        # Check the visit catalogs are identical
        visitCat1 = butler1.get('fgcmVisitCatalog').asAstropy()
        visitCat2 = butler2.get('fgcmVisitCatalog').asAstropy()

        for col in visitCat1.columns:
            if isinstance(visitCat1[col][0], str):
                testing.assert_array_equal(visitCat1[col], visitCat2[col])
            else:
                testing.assert_array_almost_equal(visitCat1[col], visitCat2[col])

        # Check that the observation catalogs have the same length
        # Detailed comparisons of the contents are below.
        starObs1 = butler1.get('fgcmStarObservations')
        starObs2 = butler2.get('fgcmStarObservations')
        self.assertEqual(len(starObs1), len(starObs2))

        # Check that the number of stars is the same and all match.
        starIds1 = butler1.get('fgcmStarIds')
        starIds2 = butler2.get('fgcmStarIds')
        self.assertEqual(len(starIds1), len(starIds2))
        matcher = esutil.htm.Matcher(11, starIds1['ra'], starIds1['dec'])
        matches = matcher.match(starIds2['ra'], starIds2['dec'], 1./3600., maxmatch=1)
        self.assertEqual(len(matches[0]), len(starIds1))

        # Check that the number of observations of each star is the same.
        testing.assert_array_equal(starIds1['nObs'][matches[1]],
                                   starIds2['nObs'][matches[0]])

        # And to test the contents, we need to unravel the observations and make
        # sure that they are matched individually, because the two catalogs
        # are constructed in a different order.
        starIndices1 = butler1.get('fgcmStarIndices')
        starIndices2 = butler2.get('fgcmStarIndices')

        test1 = np.zeros(len(starIndices1), dtype=[('ra', 'f8'),
                                                   ('dec', 'f8'),
                                                   ('x', 'f8'),
                                                   ('y', 'f8'),
                                                   ('psf_candidate', 'b1'),
                                                   ('visit', 'i4'),
                                                   ('ccd', 'i4'),
                                                   ('instMag', 'f4'),
                                                   ('instMagErr', 'f4'),
                                                   ('jacobian', 'f4')])
        test2 = np.zeros_like(test1)

        # Fill the test1 numpy recarray with sorted and unpacked data from starObs1.
        # Note that each star has a different number of observations, leading to
        # a "ragged" array that is packed in here.
        counter = 0
        obsIndex = starIndices1['obsIndex']
        for i in range(len(starIds1)):
            ind = starIds1['obsArrIndex'][matches[1][i]]
            nObs = starIds1['nObs'][matches[1][i]]
            for name in test1.dtype.names:
                test1[name][counter: counter + nObs] = starObs1[name][obsIndex][ind: ind + nObs]
            counter += nObs

        # Fill the test2 numpy recarray with sorted and unpacked data from starObs2.
        # Note that we have to match these observations per star by matching "visit"
        # (implicitly assuming each star is observed only once per visit) to ensure
        # that the observations in test2 are in the same order as test1.
        counter = 0
        obsIndex = starIndices2['obsIndex']
        for i in range(len(starIds2)):
            ind = starIds2['obsArrIndex'][matches[0][i]]
            nObs = starIds2['nObs'][matches[0][i]]
            a, b = esutil.numpy_util.match(test1['visit'][counter: counter + nObs],
                                           starObs2['visit'][obsIndex][ind: ind + nObs])
            for name in test2.dtype.names:
                test2[name][counter: counter + nObs][a] = starObs2[name][obsIndex][ind: ind + nObs][b]
            counter += nObs

        for name in test1.dtype.names:
            testing.assert_array_almost_equal(test1[name], test2[name])

    def _getMatchedVisitCat(self, rawStars, dataRefs, bandIndex, offsets):
        """
        Get a list of matched magnitudes and deltas from calibrated src catalogs.

        Parameters
        ----------
        rawStars : `lsst.afw.table.SourceCatalog`
           Fgcm standard stars
        dataRefs : `list` or `lsst.daf.persist.ButlerSubset`
           Data references for source catalogs to match
        bandIndex : `int`
           Index of the band for the source catalogs
        offsets : `np.ndarray`
           Testing calibration offsets to apply to rawStars

        Returns
        -------
        matchMag : `np.ndarray`
           Array of matched magnitudes
        matchDelta : `np.ndarray`
           Array of matched deltas between src and standard stars.
        """
        matcher = esutil.htm.Matcher(11, np.rad2deg(rawStars['coord_ra']),
                                     np.rad2deg(rawStars['coord_dec']))

        matchDelta = None
        for dataRef in dataRefs:
            src = dataRef.get()
            photoCal = dataRef.get('fgcm_photoCalib')
            src = photoCal.calibrateCatalog(src)

            gdSrc, = np.where(np.nan_to_num(src['slot_CalibFlux_flux']) > 0.0)

            matches = matcher.match(np.rad2deg(src['coord_ra'][gdSrc]),
                                    np.rad2deg(src['coord_dec'][gdSrc]),
                                    1./3600., maxmatch=1)

            srcMag = src['slot_CalibFlux_mag'][gdSrc][matches[0]]
            # Apply offset here to the catalog mag
            catMag = rawStars['mag_std_noabs'][matches[1]][:, bandIndex] + offsets[bandIndex]
            delta = srcMag - catMag
            if matchDelta is None:
                matchDelta = delta
                matchMag = catMag
            else:
                matchDelta = np.append(matchDelta, delta)
                matchMag = np.append(matchMag, catMag)

        return matchMag, matchDelta

    def _checkResult(self, result):
        """
        Check the result output from the task

        Parameters
        ----------
        result: `pipeBase.struct`
           Result structure output from a task
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
