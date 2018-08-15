# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import sys
import traceback

import numpy as np
import healpy as hp
import esutil

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import TransmissionCurve
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.pipe.tasks.photoCal import PhotoCalConfig, PhotoCalTask
from .fgcmFitCycle import FgcmFitCycleConfig, FgcmFitCycleTask
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable


import fgcm

__all__ = ['FgcmOutputProductsConfig', 'FgcmOutputProductsTask']

class FgcmOutputProductsConfig(pexConfig.Config):
    """Config for FgcmOutputProductsTask"""

    cycleNumber = pexConfig.Field(
        doc="Final fit cycle from FGCM fit",
        dtype=int,
        default=None,
    )

    # The following fields refer to calibrating from a reference
    # catalog, but in the future this might need to be expanded
    doReferenceCalibration = pexConfig.Field(
        doc="Apply 'absolute' calibration from reference catalog",
        dtype=bool,
        default=True,
    )
    refObjLoader = pexConfig.ConfigurableField(
        target=LoadIndexedReferenceObjectsTask,
        doc="reference object loader for 'absolute' photometric calibration",
    )
    photoCal = pexConfig.ConfigurableField(
        target=PhotoCalTask,
        doc="task to perform 'absolute' calibration",
    )
    referencePixelizationNside = pexConfig.Field(
        doc="Healpix nside to pixelize catalog to compare to reference",
        dtype=int,
        default=64,
    )
    referencePixelizationMinStars = pexConfig.Field(
        doc="Minimum number of stars per pixel to select for comparison",
        dtype=int,
        default=200,
    )
    referenceMinMatch = pexConfig.Field(
        doc="Minimum number of stars matched to reference catalog to be used in statistics",
        dtype=int,
        default=50,
    )
    referencePixelizationNPixels = pexConfig.Field(
        doc="Number of pixels to sample to do comparison",
        dtype=int,
        default=100,
    )

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

        self.photoCal.applyColorTerms = True
        self.photoCal.fluxField = 'flux'
        self.photoCal.magErrFloor = 0.003
        self.photoCal.match.referenceSelection.doSignalToNoise = True
        self.photoCal.match.referenceSelection.signalToNoise.minimum = 10.0
        self.photoCal.match.sourceSelection.doSignalToNoise = True
        self.photoCal.match.sourceSelection.signalToNoise.minimum = 10.0
        self.photoCal.match.sourceSelection.doFlags = True
        self.photoCal.match.sourceSelection.flags.good = ['flag_goodStar']
        self.photoCal.match.sourceSelection.flags.bad = []
        self.photoCal.match.sourceSelection.doUnresolved = False


class FgcmOutputProductsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmOutputProductsTask

    fgcmOutputProductsTask.run() takes one argument, the butler, and
    does not run on any data in the repository.
    This runner does not use any parallelization.
    """

    @staticmethod
    def getTargetList(parsedCmd):
        """
        Return a list with one element, the butler.
        """
        return [parsedCmd.butler]

    def precall(self, parsedCmd):
        return True

    def __call__(self, butler):
        """
        Parameters
        ----------
        butler: lsst.daf.persistence.Butler

        Returns
        -------
        None if self.doReturnResults is False
        An empty list if self.doReturnResults is True
        """
        task = self.TaskClass(butler=butler, config=self.config, log=self.log)

        exitStatus = 0
        if self.doRaise:
            results = task.runDataRef(butler)
        else:
            try:
                results = task.runDataRef(butler)
            except Exception as e:
                exitStatus = 1
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)

        if self.doReturnResults:
            # Not much in terms of results to return (empty list)
            return [pipeBase.Struct(exitStatus=exitStatus,
                                    results=results)]
        else:
            return [pipeBase.Struct(exitStatus=exitStatus)]

    # turn off any multiprocessing

    def run(self, parsedCmd):
        """
        Run the task, with no multiprocessing

        Parameters
        ----------
        parsedCmd: ArgumentParser parsed command line
        """

        resultList = []

        if self.precall(parsedCmd):
            # profileName = parsedCmd.profile if hasattr(parsedCmd, "profile") else None
            # log = parsedCmd.log
            targetList = self.getTargetList(parsedCmd)
            # make sure that we only get 1
            resultList = self(targetList[0])

        return resultList


class FgcmOutputProductsTask(pipeBase.CmdLineTask):
    """
    Output products from FGCM global calibration.

    Currently this is atmosphere transmission tables, in the future will
    include zeropoint files.
    """

    ConfigClass = FgcmOutputProductsConfig
    RunnerClass = FgcmOutputProductsRunner
    _DefaultName = "fgcmOutputProducts"

    def __init__(self, butler=None, **kwargs):
        """
        Instantiate an fgcmOutputProductsTask.

        Parameters
        ----------
        butler : lsst.daf.persistence.Butler
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

        if self.config.doReferenceCalibration:
            # We need the ref obj loader to get the flux field
            self.makeSubtask("refObjLoader", butler=butler)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)

        return parser

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def runDataRef(self, butler):
        """
        Make FGCM output products for use in the stack

        Parameters
        ----------
        butler:  lsst.daf.persistence.Butler

        Returns
        -------
        Empty list
        """

        # Check to make sure that the fgcmBuildStars config exists, to retrieve
        # the visit and ccd dataset tags
        if not butler.datasetExists('fgcmBuildStars_config'):
            raise RuntimeError("Cannot find fgcmBuildStars_config, which is prereq for fgcmOutputProducts")

        fgcmBuildStarsConfig = butler.get('fgcmBuildStars_config')
        self.visitDataRefName = fgcmBuildStarsConfig.visitDataRefName
        self.ccdDataRefName = fgcmBuildStarsConfig.ccdDataRefName

        # Make sure that the fit config exists, to retrieve bands and other info
        if not butler.datasetExists('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber):
            raise RuntimeError("Cannot find fgcmFitCycle_config from cycle %d " % (self.config.cycleNumber) +
                               "which is required for fgcmOutputProducts.")

        fitCycleConfig = butler.get('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber)
        self.bands = fitCycleConfig.bands

        # And make sure that the atmosphere was output properly
        if not butler.datasetExists('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber):
            raise RuntimeError("The specified fgcm cycle (%d) has not been run on this repo." %
                               (self.config.cycleNumber))

        # And make sure this is the last cycle
        if butler.datasetExists('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber + 1):
            raise RuntimeError("The task fgcmOutputProducts should only be run on the final fit cycle products")

        # Check if we need to run one final cycle to generate standards
        if not butler.datasetExists('fgcmStandardStars', fgcmcycle=self.config.cycleNumber):
            self.log.info("Standard stars will be generated using final fgcmFitCycle settings")
            self._generateFgcmStandardStars(butler)
            # Increment the cycle number
            self.useCycle = self.config.cycleNumber + 1
        else:
            self.useCycle = self.config.cycleNumber

        # And make sure that the atmosphere was output properly
        if not butler.datasetExists('fgcmAtmosphereParameters', fgcmcycle=self.useCycle):
            raise RuntimeError("The latest fgcm cycle (%d) has not been run on this repo." %
                               (self.useCycle))

        if self.config.doReferenceCalibration:
            offsets = self._computeReferenceOffsets(butler)
        else:
            offsets = np.zeros(len(self.bands))

        # Output the standard stars in stack format
        self._outputStandardStars(butler, offsets)

        # Output the atmospheres
        self._outputAtmospheres(butler)

        # Output the gray zeropoints

        return []

    def _generateFgcmStandardStars(self, butler):
        """
        Call fgcmFitCycle with settings to output the standard stars.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
        """

        fitCycleConfig = butler.get('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber)
        fitCycleConfig.cycleNumber += 1
        fitCycleConfig.maxIter = 0
        fitCycleConfig.outputStandards = True

        fitCycleTask = FgcmFitCycleTask(butler=butler, config=fitCycleConfig)
        fitCycleTask.runDataRef(butler)

    def _computeReferenceOffsets(self, butler):
        """
        Compute offsets relative to a reference catalog.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler

        Returns
        -------
        offsets: np.array of floats
           Per band zeropoint offsets
        """

        # Load the stars
        stars = butler.get('fgcmStandardStars', fgcmcycle=self.useCycle)

        # Only use stars that are observed in all the bands
        minObs = stars['ngood'].min(axis=1)

        # Depending on how things work, this might need to be configurable
        # However, I think that best results will occur when we use the same
        # pixels to do the calibration for all the bands, so I think this
        # is the right choice.

        goodStars = (minObs >= 1)
        stars = stars[goodStars]

        self.log.info("Found %d stars with at least 1 good observation in each band" %
                      (len(stars)))

        # We have to make a table for each pixel with flux/fluxErr
        sourceMapper = afwTable.SchemaMapper(stars.schema)
        sourceMapper.addMapping(stars.schema.find('id'))
        sourceMapper.addMapping(stars.schema.find('coord_ra'))
        sourceMapper.addMapping(stars.schema.find('coord_dec'))
        sourceMapper.editOutputSchema().addField('flux', type=np.float64, doc="flux")
        sourceMapper.editOutputSchema().addField('fluxErr', type=np.float64, doc="flux error")
        sourceMapper.editOutputSchema().addField('flag_GoodStar', type=bool, doc="Good flag")

        # The exposure is used to record the filter name
        exposure = afwImage.ExposureF()

        # Split up the stars (units are radians)
        theta = np.pi / 2. - stars['coord_dec']
        phi = stars['coord_ra']

        ipring = hp.ang2pix(self.config.referencePixelizationNside, theta, phi)
        h, rev = esutil.stat.histogram(ipring, rev=True)

        gdpix, = np.where(h >= self.config.referencePixelizationMinStars)

        self.log.info("Found %d pixels (nside=%d) with at least %d good stars" %
                      (gdpix.size, self.config.referencePixelizationMinStars))

        if gdpix.size < self.config.referencePixelizationNPixels:
            self.log.warn("Found fewer good pixels (%d) than preferred in configuration (%d)" %
                          (gdpix.size, self.config.referencePixelizationNPixels))
        else:
            # Sample out the pixels we want to use
            gdpix = np.random.choice(gdpix, size=self.config.referencePixelizationNPixels, replace=False)

        results = np.zeros(gdpix.size, dtype=[('hpix', 'i4'),
                                              ('nstar', 'i4', len(bands)),
                                              ('nmatch', 'i4', len(bands)),
                                              ('zp', 'f4', len(bands)),
                                              ('zpErr', 'f4', len(bands))])
        results['hpix'] = ipring[rev[rev[gdpix]]]

        # We need a boolean index to deal with catalogs...
        selected = np.zeros(len(stars), dtype=np.bool)

        refFluxFields = [None] * len(bands)

        for p, pix in enumerate(gdpix):
            i1a = rev[rev[pix]: rev[pix + 1]]

            # Need to convert to a boolean array
            selected[:] = False
            selected[i1a] = True

            for b, band in enumerate(bands):

                sourceCat = afwTable.SourceCatalog(sourceMapper.getOutputSchema())
                sourceCat.reserve(len(use))
                sourceCat.extend(sources[selected], mapper=sourceMapper)
                sourceCat['flux'] = afwImage.fluxFromABMag(stars['mag_std_noabs'][selected, b])
                sourceCat['fluxErr'] = afwImage.fluxErrFromABMagErr(stars['magerr_std'][selected, b], stars['mag_std_noabs'][selected, b])

                # Make sure we only use stars that have valid measurements
                # (This is perhaps redundant with requirements above that the
                # stars be observed in all bands, but it can't hurt)
                goodStar = (stars['mag_std_noabs'][selected, b] < 99.0)
                sourceCat['flag_goodStar'][goodStar] = True

                exposure.setFilter(afwImage.Filter(band))

                if refFluxFields[b] is None:
                    # Need to find the flux field in the reference catalog
                    # to work around limitations of DirectMatch in PhotoCal
                    ctr = lsst.geom.SpherePoint(sourceCat['coord_ra'][0],
                                                sourceCat['coord_dec'][0])
                    rad = 0.05 * lsst.geom.degrees
                    refDataTest = self.refObjLoader.loadSkyCircle(ctr, rad, band)
                    refFluxFields[b] = refDataTest.fluxField

                calConfig = self.config.photoCal.config
                calConfig.match.referenceSelection.signalToNoise.fluxField = refFluxFields[b]
                calConfig.match.referenceSelection.signalToNoise.errField = refFluxField + 'Err'

                calTask = self.config.photoCal(self.config.refObjLoader,
                                               config=calConfig,
                                               schema=sourceCat.getSchema())
                try:
                    struct = calTask.run(exposure, sourceCat)
                except:
                    self.log.warn("Calibration task failed on pixel %d for %s with %d stars." %
                                  (ipring[i1a[0]], band, len(sourceCat)))
                    continue

                results['nstar'][p, b] = len(use)
                results['nmatch'][p, b] = len(struct.arrays.refMag)
                results['zp'][p, b] = struct.zp
                results['zpErr'][p, b] = struct.sigma

        # And compute the summary statistics
        offsets = np.zeros(len(self.bands))

        for b, band in enumerate(self.bands):
            # make configurable
            ok, = np.where(results['nmatch'][:, b] >= self.config.referenceMinMatch)
            offsets[b] = np.median(results['zp'][ok, b])
            madSigma = 1.4826 * np.median(np.abs(results['zp'][ok, b] - offsets[b]))
            self.log.info("Reference catalog offset for %s band: %.6f +/- %.6f" %
                          (offsets[b], madSigma))

        return offsets

    def _outputStandardStars(self, butler, offsets):
        """
        Output standard stars in indexed reference catalog format.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
        offsets: np.array of floats
           Per band zeropoint offsets
        """

        pass


    def _outputAtmospheres(self, butler):
        """
        Output the atmospheres.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
           (used for mapper information)
        """

        # First, we need to grab the look-up table and key info
        lutCat = butler.get('fgcmLookUpTable')

        atmosphereTableName = lutCat[0]['tablename']
        elevation = lutCat[0]['elevation']
        atmLambda = lutCat[0]['atmlambda']
        lutCat = None

        # Get the fgcmBuildStarsConfig for visit name
        fgcmBuildStarsConfig = butler.get("fgcmBuildStars_config")
        visitDataRefName = fgcmBuildStarsConfig.visitDataRefName

        # Make the atmosphere table if possible
        try:
            atmTable = fgcm.FgcmAtmosphereTable.initWithTableName(atmosphereTableName)
            atmTable.loadTable()
        except IOError:
            atmTable = None

        if atmTable is None:
            # Try to use MODTRAN instead
            try:
                modGen = fgcm.ModtranGenerator(elevation)
                lambdaRange = np.array([atmLambda[0], atmLambda[-1]]) / 10.
                lambdaStep = (atmLambda[1] - atmLambda[0]) / 10.
            except (ValueError, IOError):
                raise RuntimeError("FGCM look-up-table generated without modtran, "
                                   "but modtran not configured to run.")

        # Next, we need to grab the atmosphere parameters
        atmCat = butler.get('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        zenith = np.degrees(np.arccos(1. / atmCat['seczenith']))

        for i, visit in enumerate(atmCat['visit']):
            if atmTable is not None:
                # Interpolate the atmosphere table
                atmVals = atmTable.interpolateAtmosphere(pmb=atmCat[i]['pmb'],
                                                         pwv=atmCat[i]['pwv'],
                                                         o3=atmCat[i]['o3'],
                                                         tau=atmCat[i]['tau'],
                                                         alpha=atmCat[i]['alpha'],
                                                         zenith=zenith[i])
            else:
                # Run modtran
                modAtm = modGen(pmb=atmCat[i]['pmb'],
                                pwv=atmCat[i]['pwv'],
                                o3=atmCat[i]['o3'],
                                tau=atmCat[i]['tau'],
                                alpha=atmCat[i]['alpha'],
                                zenith=zenith[i],
                                lambdaRange=lambdaRange,
                                lambdaStep=lambdaStep)
                atmVals = modAtm['COMBINED']

            # Now need to create something to persist...
            curve = TransmissionCurve.makeSpatiallyConstant(throughput=atmVals,
                                                            wavelengths=atmLambda,
                                                            throughputAtMin=atmVals[0],
                                                            throughputAtMax=atmVals[-1])

            butler.put(curve, "transmission_atmosphere_fgcm",
                       dataId={visitDataRefName: visit})
