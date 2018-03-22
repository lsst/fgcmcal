# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import os
import inspect
import shutil
import numpy as np

import lsst.daf.persistence as dafPersistence
import lsst.log

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

    def tearDown(self):
        if getattr(self, 'config', None) is not None:
            del self.config

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

    def _runFgcmFitCycle(self, nZp, nGoodZp):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])

        args = [self.inputDir, '--output', self.testDir,
                '--doraise']
        args.extend(self.otherArgs)

        args = [self.inputDir, '--output', self.testDir, '--doraise']
        result = fgcmcal.FgcmFitCycleTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

        butler = dafPersistence.butler.Butler(self.testDir)

        zps = butler.get('fgcmZeropoints')

        self.assertEqual(nZp, len(zps))

        gd, = np.where(zps['fgcmflag'] == 1)
        self.assertEqual(nGoodZp, len(gd))

    def _checkResult(self, result):
        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.assertEqual(result.resultList[0].exitStatus, 0)

    def tearDown(self):
        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir, True)


