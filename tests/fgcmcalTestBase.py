# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import os
import inspect

from lsst.fgcmcal import fgcmcal

class FgcmcalTestBase(object):
    """
    Base class for fgcmcal tests, to genericize some test running and setup.

    Derive from this first, then from TestCase.
    """

    def setUp_base(self):
        """
        Call from your child class's setUp() to get variables built.

        Parameters
        ----------
        None

        """

        self.config = None
        self.inputDir = None
        self.logLevel = None

    def tearDown(self):
        if getattr(self, 'config', None) is not None:
            del self.config

    def _runFgcmMakeLut(self, outputDir):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])
        args = [self.inputDir, '--output', outputDir,
                '--doraise']
        args.extend(self.otherArgs)

        result = fgcmcal.fgcmMakeLutTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

    def _runFgcmBuildStars(self, outputDir):
        """
        """

        if self.logLevel is not None:
            self.otherArgs.extend(['--loglevel', 'fgcmcal=%s'%self.logLevel])
        args = [self.inputDir, '--output', outputDir,
                '--doraise']
        args.extend(self.otherArgs)

        result = fgcmcal.fgcmBuildStarsTask.parseAndRun(args=args, config=self.config)
        self._checkResult(result)

    def _checkResult(self, result):
        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.assertEqual(result.resultList[0].exitStatus, 0)




