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
data from testdata_jointcal for Gen3 repos.
"""

import os
import shutil
import numpy as np
import esutil

import lsst.daf.butler as dafButler
import lsst.pipe.base as pipeBase
import lsst.geom as geom
from lsst.pipe.base import Pipeline, ExecutionResources
from lsst.ctrl.mpexec import SimplePipelineExecutor

import lsst.fgcmcal as fgcmcal

ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmcalTestBase(object):
    """Base class for gen3 fgcmcal tests, to genericize some test running and setup.

    Derive from this first, then from TestCase.
    """
    @classmethod
    def _importRepository(cls, instrument, exportPath, exportFile):
        """Import a test repository into self.testDir

        Parameters
        ----------
        instrument : `str`
            Full string name for the instrument.
        exportPath : `str`
            Path to location of repository to export.
        exportFile : `str`
            Filename of export data.
        """
        cls.repo = os.path.join(cls.testDir, 'testrepo')

        # Make the repo and retrieve a writeable Butler
        _ = dafButler.Butler.makeRepo(cls.repo)
        butler = dafButler.Butler(cls.repo, writeable=True)
        # Register the instrument
        instrInstance = pipeBase.Instrument.from_string(instrument)
        instrInstance.register(butler.registry)
        # Import the exportFile
        butler.import_(directory=exportPath, filename=exportFile,
                       transfer='symlink',
                       skip_dimensions={'instrument', 'detector', 'physical_filter'})

    def _runPipeline(self, repo, pipelineFile, queryString='',
                     inputCollections=None, outputCollection=None,
                     configFiles={}, configOptions={},
                     registerDatasetTypes=False):
        """Run a pipeline via pipetask.

        Parameters
        ----------
        repo : `str`
            Gen3 repository yaml file.
        pipelineFile : `str`
            Pipeline definition file.
        queryString : `str`, optional
            Where query that defines the data to use.
        inputCollections : `list` [`str`], optional
            Input collections list.
        outputCollection : `str`, optional
            Output collection name.
        configFiles : `dict` [`list`], optional
            Dictionary of config files.  The key of the ``configFiles``
            dict is the relevant task label.  The value of ``configFiles``
            is a list of config files to apply (in order) to that task.
        configOptions : `dict` [`dict`], optional
            Dictionary of individual config options.  The key of the
            ``configOptions`` dict is the relevant task label.  The value
            of ``configOptions`` is another dict that contains config
            key/value overrides to apply.
        configOptions : `list` [`str`], optional
            List of individual config options to use.  Each string will
            be of the form ``taskName:configField=value``.
        registerDatasetTypes : `bool`, optional
            Register new dataset types?

        Returns
        -------
        exit_code : `int`
            Exit code for pipetask run.

        Raises
        ------
        RuntimeError : Raised if the "pipetask" call fails.
        """
        butler = SimplePipelineExecutor.prep_butler(repo,
                                                    inputs=inputCollections,
                                                    output=outputCollection)
        pipeline = Pipeline.fromFile(pipelineFile)
        for taskName, fileList in configFiles.items():
            for fileName in fileList:
                pipeline.addConfigFile(taskName, fileName)
        for taskName, configDict in configOptions.items():
            for option, value in configDict.items():
                pipeline.addConfigOverride(taskName, option, value)

        resources = ExecutionResources(num_cores=1)

        executor = SimplePipelineExecutor.from_pipeline(pipeline,
                                                        where=queryString,
                                                        butler=butler,
                                                        resources=resources)
        quanta = executor.run(register_dataset_types=registerDatasetTypes)

        return len(quanta)

    def _testFgcmMakeLut(self, instName, testName, nBand, i0Std, i0Recon, i10Std, i10Recon):
        """Test running of FgcmMakeLutTask

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        nBand : `int`
            Number of bands tested.
        i0Std : `np.ndarray'
            Values of i0Std to compare to.
        i10Std : `np.ndarray`
            Values of i10Std to compare to.
        i0Recon : `np.ndarray`
            Values of reconstructed i0 to compare to.
        i10Recon : `np.ndarray`
            Values of reconsntructed i10 to compare to.
        """
        instCamel = instName.title()

        configFiles = {'fgcmMakeLut': [os.path.join(ROOT,
                                                    'config',
                                                    f'fgcmMakeLut{instCamel}.py')]}
        outputCollection = f'{instName}/{testName}/lut'

        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       f'fgcmMakeLut{instCamel}.yaml'),
                          configFiles=configFiles,
                          inputCollections=[f'{instName}/calib', f'{instName}/testdata'],
                          outputCollection=outputCollection,
                          registerDatasetTypes=True)

        # Check output values
        butler = dafButler.Butler(self.repo)
        lutCat = butler.get('fgcmLookUpTable',
                            collections=[outputCollection],
                            instrument=instName)
        fgcmLut, lutIndexVals, lutStd = fgcmcal.utilities.translateFgcmLut(lutCat, {})

        self.assertEqual(nBand, len(lutIndexVals[0]['FILTERNAMES']))

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

        # Check that the standard atmosphere was output and non-zero.
        atmStd = butler.get('fgcm_standard_atmosphere',
                            collections=[outputCollection],
                            instrument=instName)
        bounds = atmStd.getWavelengthBounds()
        lambdas = np.linspace(bounds[0], bounds[1], 1000)
        tputs = atmStd.sampleAt(position=geom.Point2D(0.0, 0.0), wavelengths=lambdas)
        self.assertGreater(np.min(tputs), 0.0)

        # Check that the standard passbands were output and non-zero.
        for physical_filter in fgcmLut.filterNames:
            passband = butler.get('fgcm_standard_passband',
                                  collections=[outputCollection],
                                  instrument=instName,
                                  physical_filter=physical_filter)
            tputs = passband.sampleAt(position=geom.Point2D(0.0, 0.0), wavelengths=lambdas)
            self.assertEqual(np.min(tputs), 0.0)
            self.assertGreater(np.max(tputs), 0.0)

    def _testFgcmBuildStarsTable(self, instName, testName, queryString, visits, nStar, nObs,
                                 refcatCollection="refcats/gen2"):
        """Test running of FgcmBuildStarsTableTask

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        queryString : `str`
            Query to send to the pipetask.
        visits : `list`
            List of visits to calibrate.
        nStar : `int`
            Number of stars expected.
        nObs : `int`
            Number of observations of stars expected.
        refcatCollection : `str`, optional
            Name of reference catalog collection.
        """
        instCamel = instName.title()

        configFiles = {'fgcmBuildStarsTable': [os.path.join(ROOT,
                                                            'config',
                                                            f'fgcmBuildStarsTable{instCamel}.py')]}
        outputCollection = f'{instName}/{testName}/buildstars'

        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       'fgcmBuildStarsTable%s.yaml' % (instCamel)),
                          configFiles=configFiles,
                          inputCollections=[f'{instName}/{testName}/lut',
                                            refcatCollection],
                          outputCollection=outputCollection,
                          queryString=queryString,
                          registerDatasetTypes=True)

        butler = dafButler.Butler(self.repo)

        visitCat = butler.get('fgcmVisitCatalog', collections=[outputCollection],
                              instrument=instName)
        self.assertEqual(len(visits), len(visitCat))

        starIds = butler.get('fgcmStarIds', collections=[outputCollection],
                             instrument=instName)
        self.assertEqual(len(starIds), nStar)

        starObs = butler.get('fgcmStarObservations', collections=[outputCollection],
                             instrument=instName)
        self.assertEqual(len(starObs), nObs)

    def _testFgcmBuildFromIsolatedStars(self, instName, testName, queryString, visits, nStar, nObs,
                                        refcatCollection="refcats/gen2"):
        """Test running of FgcmBuildFromIsolatedStarsTask.

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        queryString : `str`
            Query to send to the pipetask.
        visits : `list`
            List of visits to calibrate.
        nStar : `int`
            Number of stars expected.
        nObs : `int`
            Number of observations of stars expected.
        refcatCollection : `str`, optional
            Name of reference catalog collection.
        """
        instCamel = instName.title()

        configFiles = {'fgcmBuildFromIsolatedStars': [
            os.path.join(ROOT,
                         'config',
                         f'fgcmBuildFromIsolatedStars{instCamel}.py')
        ]}
        outputCollection = f'{instName}/{testName}/buildstars'

        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       'fgcmBuildFromIsolatedStars%s.yaml' % (instCamel)),
                          configFiles=configFiles,
                          inputCollections=[f'{instName}/{testName}/lut',
                                            refcatCollection],
                          outputCollection=outputCollection,
                          queryString=queryString,
                          registerDatasetTypes=True)

        butler = dafButler.Butler(self.repo)

        visitCat = butler.get('fgcmVisitCatalog', collections=[outputCollection],
                              instrument=instName)
        self.assertEqual(len(visits), len(visitCat))

        starIds = butler.get('fgcm_star_ids', collections=[outputCollection],
                             instrument=instName)
        self.assertEqual(len(starIds), nStar)

        starObs = butler.get('fgcm_star_observations', collections=[outputCollection],
                             instrument=instName)
        self.assertEqual(len(starObs), nObs)

    def _testFgcmFitCycle(self, instName, testName, cycleNumber,
                          nZp, nGoodZp, nOkZp, nBadZp, nStdStars, nPlots,
                          skipChecks=False, extraConfig=None):
        """Test running of FgcmFitCycleTask

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        cycleNumber : `int`
            Fit cycle number.
        nZp : `int`
            Number of zeropoints created by the task.
        nGoodZp : `int`
            Number of good (photometric) zeropoints created.
        nOkZp : `int`
            Number of constrained zeropoints (photometric or not).
        nBadZp : `int`
            Number of unconstrained (bad) zeropoints.
        nStdStars : `int`
            Number of standard stars produced.
        nPlots : `int`
            Number of plots produced.
        skipChecks : `bool`, optional
            Skip number checks, when running less-than-final cycle.
        extraConfig : `str`, optional
            Name of an extra config file to apply.
        """
        instCamel = instName.title()

        configFiles = {'fgcmFitCycle': [os.path.join(ROOT,
                                                     'config',
                                                     f'fgcmFitCycle{instCamel}.py')]}
        if extraConfig is not None:
            configFiles['fgcmFitCycle'].append(extraConfig)

        outputCollection = f'{instName}/{testName}/fit'

        if cycleNumber == 0:
            inputCollections = [f'{instName}/{testName}/buildstars']
        else:
            # In these tests we are running the fit cycle task multiple
            # times into the same output collection.  This code allows
            # us to find the correct chained input collections to use
            # so that we can both read from previous runs in the output
            # collection and write to a new run in the output collection.
            # Note that this behavior is handled automatically by the
            # pipetask command-line interface, but not by the python
            # API.
            butler = dafButler.Butler(self.repo)
            inputCollections = list(butler.registry.getCollectionChain(outputCollection))

        cwd = os.getcwd()
        runDir = os.path.join(self.testDir, testName)
        os.makedirs(runDir, exist_ok=True)
        os.chdir(runDir)

        configOptions = {'fgcmFitCycle':
                         {'cycleNumber': f'{cycleNumber}',
                          'connections.previousCycleNumber': f'{cycleNumber - 1}',
                          'connections.cycleNumber': f'{cycleNumber}'}}
        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       f'fgcmFitCycle{instCamel}.yaml'),
                          configFiles=configFiles,
                          inputCollections=inputCollections,
                          outputCollection=outputCollection,
                          configOptions=configOptions,
                          registerDatasetTypes=True)

        os.chdir(cwd)

        if skipChecks:
            return

        butler = dafButler.Butler(self.repo)

        # Check that the expected number of plots are there.
        plotDatasets = list(butler.registry.queryDatasets(
            f"fgcm_Cycle{cycleNumber}_*Plot",
            collections=[outputCollection],
        ))
        self.assertEqual(len(plotDatasets), nPlots)

        zps = butler.get(f'fgcm_Cycle{cycleNumber}_Zeropoints',
                         collections=[outputCollection],
                         instrument=instName)
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

        stds = butler.get(f'fgcm_Cycle{cycleNumber}_StandardStars',
                          collections=[outputCollection],
                          instrument=instName)

        self.assertEqual(len(stds), nStdStars)

    def _testFgcmOutputProducts(self, instName, testName,
                                zpOffsets, testVisit, testCcd, testFilter, testBandIndex,
                                testSrc=True):
        """Test running of FgcmOutputProductsTask.

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        zpOffsets : `np.ndarray`
            Zeropoint offsets expected.
        testVisit : `int`
            Visit id to check for round-trip computations.
        testCcd : `int`
            Ccd id to check for round-trip computations.
        testFilter : `str`
            Filtername for testVisit/testCcd.
        testBandIndex : `int`
            Band index for testVisit/testCcd.
        testSrc : `bool`, optional
            Test the source catalogs?  (Only if available in dataset.)
        """
        instCamel = instName.title()

        configFiles = {'fgcmOutputProducts': [os.path.join(ROOT,
                                                           'config',
                                                           f'fgcmOutputProducts{instCamel}.py')]}
        inputCollection = f'{instName}/{testName}/fit'
        outputCollection = f'{instName}/{testName}/fit/output'

        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       'fgcmOutputProducts%s.yaml' % (instCamel)),
                          configFiles=configFiles,
                          inputCollections=[inputCollection],
                          outputCollection=outputCollection,
                          registerDatasetTypes=True)

        butler = dafButler.Butler(self.repo)
        offsetCat = butler.get('fgcmReferenceCalibrationOffsets',
                               collections=[outputCollection], instrument=instName)
        offsets = offsetCat['offset'][:]
        self.assertFloatsAlmostEqual(offsets, zpOffsets, atol=1e-6)

        config = butler.get('fgcmOutputProducts_config',
                            collections=[outputCollection], instrument=instName)

        rawStars = butler.get(f'fgcm_Cycle{config.connections.cycleNumber}_StandardStars',
                              collections=[inputCollection], instrument=instName)

        # Test the fgcm_photoCalib output
        zptCat = butler.get(f'fgcm_Cycle{config.connections.cycleNumber}_Zeropoints',
                            collections=[inputCollection], instrument=instName)

        good = (zptCat['fgcmFlag'] < 16)
        bad = (zptCat['fgcmFlag'] >= 16)

        # Read in all the calibrations, these should all be there
        # This test is simply to ensure that all the photoCalib files exist
        visits = np.unique(zptCat['visit'][good])
        photoCalibDict = {}
        for visit in visits:
            expCat = butler.get('fgcmPhotoCalibCatalog',
                                visit=visit,
                                collections=[outputCollection], instrument=instName)
            for row in expCat:
                if row['visit'] == visit:
                    photoCalibDict[(visit, row['id'])] = row.getPhotoCalib()

        # Check that all of the good photocalibs are there.
        for rec in zptCat[good]:
            self.assertTrue((rec['visit'], rec['detector']) in photoCalibDict)

        # Check that none of the bad photocalibs are there.
        for rec in zptCat[bad]:
            self.assertFalse((rec['visit'], rec['detector']) in photoCalibDict)

        # We do round-trip value checking on just the final one (chosen arbitrarily)
        testCal = photoCalibDict[(testVisit, testCcd)]

        if testSrc:
            src = butler.get('src', visit=int(testVisit), detector=int(testCcd),
                             collections=[outputCollection], instrument=instName)

            # Only test sources with positive flux
            gdSrc = (src['slot_CalibFlux_instFlux'] > 0.0)

            # We need to apply the calibration offset to the fgcmzpt (which is internal
            # and doesn't know about that yet)
            testZpInd, = np.where((zptCat['visit'] == testVisit)
                                  & (zptCat['detector'] == testCcd))
            fgcmZpt = (zptCat['fgcmZpt'][testZpInd] + offsets[testBandIndex]
                       + zptCat['fgcmDeltaChrom'][testZpInd])
            fgcmZptGrayErr = np.sqrt(zptCat['fgcmZptVar'][testZpInd])

            if config.doComposeWcsJacobian:
                # The raw zeropoint needs to be modified to know about the wcs jacobian
                refs = butler.registry.queryDatasets('camera', dimensions=['instrument'],
                                                     collections=...)
                camera = butler.get(list(refs)[0])
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
                zptMeanCalMags[i] = fgcmZpt[0] - 2.5*np.log10(rec['slot_CalibFlux_instFlux'])

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
            whereClause = f"instrument='{instName:s}' and visit={testVisit:d}"
            srcRefs = butler.registry.queryDatasets('src', dimensions=['visit'],
                                                    collections='%s/testdata' % (instName),
                                                    where=whereClause,
                                                    findFirst=True)
            photoCals = []
            for srcRef in srcRefs:
                photoCals.append(photoCalibDict[(testVisit, srcRef.dataId['detector'])])

            matchMag, matchDelta = self._getMatchedVisitCat(butler, srcRefs, photoCals,
                                                            rawStars, testBandIndex, offsets)

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
        visitCatalog = butler.get('fgcmVisitCatalog', collections=[inputCollection],
                                  instrument=instName)
        lutCat = butler.get('fgcmLookUpTable', collections=[inputCollection],
                            instrument=instName)

        testTrans = butler.get('transmission_atmosphere_fgcm',
                               visit=visitCatalog[0]['visit'],
                               collections=[outputCollection], instrument=instName)
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
        self.assertFloatsAlmostEqual(testResp/ratio, lutCat[0]['atmStdTrans'], atol=0.2)

        # The second should be close to the first, but there is the airmass
        # difference so they aren't identical.
        testTrans2 = butler.get('transmission_atmosphere_fgcm',
                                visit=visitCatalog[1]['visit'],
                                collections=[outputCollection], instrument=instName)
        testResp2 = testTrans2.sampleAt(position=geom.Point2D(0, 0),
                                        wavelengths=lutCat[0]['atmLambda'])

        # As above, we scale by the ratio to compare the shape of the curve.
        ratio = np.median(testResp/testResp2)
        self.assertFloatsAlmostEqual(testResp/ratio, testResp2, atol=0.04)

    def _testFgcmMultiFit(self, instName, testName, queryString, visits, zpOffsets,
                          refcatCollection="refcats/gen2"):
        """Test running the full pipeline with multiple fit cycles.

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        queryString : `str`
            Query to send to the pipetask.
        visits : `list`
            List of visits to calibrate.
        zpOffsets : `np.ndarray`
            Zeropoint offsets expected.
        refcatCollection : `str`, optional
            Name of reference catalog collection.
        """
        instCamel = instName.title()

        configFiles = {'fgcmBuildFromIsolatedStars': [
            os.path.join(ROOT,
                         'config',
                         f'fgcmBuildFromIsolatedStars{instCamel}.py'
                         )],
                       'fgcmFitCycle': [os.path.join(ROOT,
                                                     'config',
                                                     f'fgcmFitCycle{instCamel}.py')],
                       'fgcmOutputProducts': [os.path.join(ROOT,
                                                           'config',
                                                           f'fgcmOutputProducts{instCamel}.py')]}
        outputCollection = f'{instName}/{testName}/unified'

        cwd = os.getcwd()
        runDir = os.path.join(self.testDir, testName)
        os.makedirs(runDir)
        os.chdir(runDir)

        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       f'fgcmFullPipeline{instCamel}.yaml'),
                          configFiles=configFiles,
                          inputCollections=[f'{instName}/{testName}/lut',
                                            refcatCollection],
                          outputCollection=outputCollection,
                          queryString=queryString,
                          registerDatasetTypes=True)

        os.chdir(cwd)

        butler = dafButler.Butler(self.repo)

        offsetCat = butler.get('fgcmReferenceCalibrationOffsets',
                               collections=[outputCollection], instrument=instName)
        offsets = offsetCat['offset'][:]
        self.assertFloatsAlmostEqual(offsets, zpOffsets, atol=1e-6)

    def _getMatchedVisitCat(self, butler, srcHandles, photoCals,
                            rawStars, bandIndex, offsets):
        """
        Get a list of matched magnitudes and deltas from calibrated src catalogs.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
        srcHandles : `list`
           Handles of source catalogs.
        photoCals : `list`
           photoCalib objects, matched to srcHandles.
        rawStars : `lsst.afw.table.SourceCatalog`
           Fgcm standard stars.
        bandIndex : `int`
           Index of the band for the source catalogs.
        offsets : `np.ndarray`
           Testing calibration offsets to apply to rawStars.

        Returns
        -------
        matchMag : `np.ndarray`
           Array of matched magnitudes.
        matchDelta : `np.ndarray`
           Array of matched deltas between src and standard stars.
        """
        matcher = esutil.htm.Matcher(11, np.rad2deg(rawStars['coord_ra']),
                                     np.rad2deg(rawStars['coord_dec']))

        matchDelta = None
        for srcHandle, photoCal in zip(srcHandles, photoCals):
            src = butler.get(srcHandle)
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

    def _testFgcmCalibrateTract(self, instName, testName, visits, tract, skymapName,
                                rawRepeatability, filterNCalibMap):
        """Test running of FgcmCalibrateTractTask

        Parameters
        ----------
        instName : `str`
            Short name of the instrument.
        testName : `str`
            Base name of the test collection.
        visits : `list`
            List of visits to calibrate.
        tract : `int`
            Tract number.
        skymapName : `str`
            Name of the sky map.
        rawRepeatability : `np.array`
            Expected raw repeatability after convergence.
            Length should be number of bands.
        filterNCalibMap : `dict`
            Mapping from filter name to number of photoCalibs created.
        """
        instCamel = instName.title()

        configFiles = {'fgcmCalibrateTractTable':
                       [os.path.join(ROOT,
                                     'config',
                                     f'fgcmCalibrateTractTable{instCamel}.py')]}

        outputCollection = f'{instName}/{testName}/tract'

        inputCollections = [f'{instName}/{testName}/lut',
                            'refcats/gen2']

        queryString = f"tract={tract:d} and skymap='{skymapName:s}'"

        self._runPipeline(self.repo,
                          os.path.join(ROOT,
                                       'pipelines',
                                       f'fgcmCalibrateTractTable{instCamel:s}.yaml'),
                          queryString=queryString,
                          configFiles=configFiles,
                          inputCollections=inputCollections,
                          outputCollection=outputCollection,
                          registerDatasetTypes=True)

        butler = dafButler.Butler(self.repo)

        whereClause = f"instrument='{instName:s}' and tract={tract:d} and skymap='{skymapName:s}'"

        repRefs = butler.registry.queryDatasets('fgcmRawRepeatability',
                                                dimensions=['tract'],
                                                collections=outputCollection,
                                                where=whereClause)

        repeatabilityCat = butler.get(list(repRefs)[0])
        repeatability = repeatabilityCat['rawRepeatability'][:]
        self.assertFloatsAlmostEqual(repeatability, rawRepeatability, atol=4e-6)

        # Check that the number of photoCalib objects in each filter are what we expect
        for filterName in filterNCalibMap.keys():
            whereClause = (f"instrument='{instName:s}' and tract={tract:d} and "
                           f"physical_filter='{filterName:s}' and skymap='{skymapName:s}'")

            refs = butler.registry.queryDatasets('fgcmPhotoCalibTractCatalog',
                                                 dimensions=['tract', 'physical_filter'],
                                                 collections=outputCollection,
                                                 where=whereClause)

            count = 0
            for ref in set(refs):
                expCat = butler.get(ref)
                test, = np.where((expCat['visit'] > 0) & (expCat['id'] >= 0))
                count += test.size

            self.assertEqual(count, filterNCalibMap[filterName])

        # Check that every visit got a transmission
        for visit in visits:
            whereClause = (f"instrument='{instName:s}' and tract={tract:d} and "
                           f"visit={visit:d} and skymap='{skymapName:s}'")
            refs = butler.registry.queryDatasets('transmission_atmosphere_fgcm_tract',
                                                 dimensions=['tract', 'visit'],
                                                 collections=outputCollection,
                                                 where=whereClause)
            self.assertEqual(len(set(refs)), 1)

    @classmethod
    def tearDownClass(cls):
        """Tear down and clear directories.
        """
        if os.path.exists(cls.testDir):
            shutil.rmtree(cls.testDir, True)
