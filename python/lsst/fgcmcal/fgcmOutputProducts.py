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
"""Make the final fgcmcal output products.

This task takes the final output from fgcmFitCycle and produces the following
outputs for use in the DM stack: the FGCM standard stars in a reference
catalog format; the model atmospheres in "transmission_atmosphere_fgcm"
format; and the zeropoints in "fgcm_photoCalib" format.  Optionally, the
task can transfer the 'absolute' calibration from a reference catalog
to put the fgcm standard stars in units of Jansky.  This is accomplished
by matching stars in a sample of healpix pixels, and applying the median
offset per band.
"""

import sys
import traceback
import copy

import numpy as np
import healpy as hp
import esutil
from astropy import units

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import TransmissionCurve
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.pipe.tasks.photoCal import PhotoCalTask
from .fgcmFitCycle import FgcmFitCycleTask
import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.meas.algorithms import IndexerRegistry
from lsst.meas.algorithms import DatasetConfig

import fgcm

__all__ = ['FgcmOutputProductsConfig', 'FgcmOutputProductsTask', 'FgcmOutputProductsRunner']


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
        doc="Transfer 'absolute' calibration from reference catalog",
        dtype=bool,
        default=True,
    )
    doRefcatOutput = pexConfig.Field(
        doc="Output standard stars in reference catalog format",
        dtype=bool,
        default=True,
    )
    doAtmosphereOutput = pexConfig.Field(
        doc="Output atmospheres in transmission_atmosphere_fgcm format",
        dtype=bool,
        default=True,
    )
    doZeropointOutput = pexConfig.Field(
        doc="Output zeropoints in fgcm_photoCalib format",
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
        doc="Healpix nside to pixelize catalog to compare to reference catalog",
        dtype=int,
        default=64,
    )
    referencePixelizationMinStars = pexConfig.Field(
        doc=("Minimum number of stars per healpix pixel to select for comparison"
             "to the specified reference catalog"),
        dtype=int,
        default=200,
    )
    referenceMinMatch = pexConfig.Field(
        doc="Minimum number of stars matched to reference catalog to be used in statistics",
        dtype=int,
        default=50,
    )
    referencePixelizationNPixels = pexConfig.Field(
        doc=("Number of healpix pixels to sample to do comparison. "
             "Doing too many will take a long time and not yield any more "
             "precise results because the final number is the median offset "
             "(per band) from the set of pixels."),
        dtype=int,
        default=100,
    )
    datasetConfig = pexConfig.ConfigField(
        dtype=DatasetConfig,
        doc="Configuration for writing/reading ingested catalog",
    )

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

        # In order to transfer the "absolute" calibration from a reference
        # catalog to the relatively calibrated FGCM standard stars (one number
        # per band), we use the PhotoCalTask to match stars in a sample of healpix
        # pixels.  These basic settings ensure that only well-measured, good stars
        # from the source and reference catalogs are used for the matching.
        self.photoCal.applyColorTerms = True
        self.photoCal.fluxField = 'instFlux'
        self.photoCal.magErrFloor = 0.003
        self.photoCal.match.referenceSelection.doSignalToNoise = True
        self.photoCal.match.referenceSelection.signalToNoise.minimum = 10.0
        self.photoCal.match.sourceSelection.doSignalToNoise = True
        self.photoCal.match.sourceSelection.signalToNoise.minimum = 10.0
        self.photoCal.match.sourceSelection.signalToNoise.fluxField = 'instFlux'
        self.photoCal.match.sourceSelection.signalToNoise.errField = 'instFluxErr'
        self.photoCal.match.sourceSelection.doFlags = True
        self.photoCal.match.sourceSelection.flags.good = []
        self.photoCal.match.sourceSelection.flags.bad = ['flag_badStar']
        self.photoCal.match.sourceSelection.doUnresolved = False
        self.datasetConfig.ref_dataset_name = 'fgcm_stars'


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
        butler: `lsst.daf.persistence.Butler`

        Returns
        -------
        exitStatus: `list` with `pipeBase.Struct`
           exitStatus (0: success; 1: failure)
           if self.doReturnResults also
           results (`np.array` with absolute zeropoint offsets)
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
            # The results here are the zeropoint offsets for each band
            return [pipeBase.Struct(exitStatus=exitStatus,
                                    results=results)]
        else:
            return [pipeBase.Struct(exitStatus=exitStatus)]

    def run(self, parsedCmd):
        """
        Run the task, with no multiprocessing

        Parameters
        ----------
        parsedCmd: ArgumentParser parsed command line
        """

        resultList = []

        if self.precall(parsedCmd):
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
        butler : `lsst.daf.persistence.Butler`
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

        if self.config.doReferenceCalibration:
            # We need the ref obj loader to get the flux field
            self.makeSubtask("refObjLoader", butler=butler)

        if self.config.doRefcatOutput:
            self.indexer = IndexerRegistry[self.config.datasetConfig.indexer.name](
                self.config.datasetConfig.indexer.active)

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def runDataRef(self, butler):
        """
        Make FGCM output products for use in the stack

        Parameters
        ----------
        butler:  `lsst.daf.persistence.Butler`

        Returns
        -------
        offsets: `pipeBase.Struct`
           A structure with array of zeropoint offsets
        """

        # Check to make sure that the fgcmBuildStars config exists, to retrieve
        # the visit and ccd dataset tags
        if not butler.datasetExists('fgcmBuildStars_config'):
            raise RuntimeError("Cannot find fgcmBuildStars_config, which is prereq for fgcmOutputProducts")

        fgcmBuildStarsConfig = butler.get('fgcmBuildStars_config')
        self.visitDataRefName = fgcmBuildStarsConfig.visitDataRefName
        self.ccdDataRefName = fgcmBuildStarsConfig.ccdDataRefName
        self.filterMap = fgcmBuildStarsConfig.filterMap

        # Make sure that the fit config exists, to retrieve bands and other info
        if not butler.datasetExists('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber):
            raise RuntimeError("Cannot find fgcmFitCycle_config from cycle %d " % (self.config.cycleNumber) +
                               "which is required for fgcmOutputProducts.")

        fitCycleConfig = butler.get('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber)
        self.bands = fitCycleConfig.bands
        self.superStarSubCcd = fitCycleConfig.superStarSubCcd
        self.chebyshevOrder = fitCycleConfig.superStarSubCcdChebyshevOrder

        # And make sure that the atmosphere was output properly
        if not butler.datasetExists('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber):
            raise RuntimeError("The specified fgcm cycle (%d) has not been run on this repo." %
                               (self.config.cycleNumber))

        # And make sure this is the last cycle
        if butler.datasetExists('fgcmFitCycle_config', fgcmcycle=self.config.cycleNumber + 1):
            raise RuntimeError("The task fgcmOutputProducts should only be run"
                               "on the final fit cycle products")

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
        if self.config.doRefcatOutput:
            self._outputStandardStars(butler, offsets)

        # Output the gray zeropoints
        if self.config.doZeropointOutput:
            self._outputZeropoints(butler, offsets)

        # Output the atmospheres
        if self.config.doAtmosphereOutput:
            self._outputAtmospheres(butler)

        # We return the zp offsets
        return pipeBase.Struct(offsets=offsets)

    def _generateFgcmStandardStars(self, butler):
        """
        Call fgcmFitCycle with settings to output the standard stars.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
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

        This method splits the star catalog into healpix pixels
        and computes the calibration transfer for a sample of
        these pixels to approximate the 'absolute' calibration
        values (on for each band) to apply to transfer the
        absolute scale.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`

        Returns
        -------
        offsets: `np.array` of floats
           Per band zeropoint offsets
        """

        # Load the stars
        stars = butler.get('fgcmStandardStars', fgcmcycle=self.useCycle)

        # Only use stars that are observed in all the bands
        # This will ensure that we use the same healpix pixels for the absolute
        # calibration of each band.
        minObs = stars['ngood'].min(axis=1)

        goodStars = (minObs >= 1)
        stars = stars[goodStars]

        self.log.info("Found %d stars with at least 1 good observation in each band" %
                      (len(stars)))

        # We have to make a table for each pixel with flux/fluxErr
        # This is a temporary table generated for input to the photoCal task.
        # These fluxes are not instFlux (they are top-of-the-atmosphere approximate and
        # have had chromatic corrections applied to get to the standard system
        # specified by the atmosphere/instrumental parameters), nor are they
        # in Jansky (since they don't have a proper absolute calibration: the overall
        # zeropoint is estimated from the telescope size, etc.)
        sourceMapper = afwTable.SchemaMapper(stars.schema)
        sourceMapper.addMinimalSchema(afwTable.SimpleTable.makeMinimalSchema())
        sourceMapper.editOutputSchema().addField('instFlux', type=np.float64,
                                                 doc="instrumental flux (counts)")
        sourceMapper.editOutputSchema().addField('instFluxErr', type=np.float64,
                                                 doc="instrumental flux error (counts)")
        badStarKey = sourceMapper.editOutputSchema().addField('flag_badStar',
                                                              type='Flag',
                                                              doc="bad flag")

        # Split up the stars
        # Note that there is an assumption here that the ra/dec coords stored
        # on-disk are in radians, and therefore that starObs['coord_ra'] /
        # starObs['coord_dec'] return radians when used as an array of numpy float64s.
        theta = np.pi / 2. - stars['coord_dec']
        phi = stars['coord_ra']

        ipring = hp.ang2pix(self.config.referencePixelizationNside, theta, phi)
        h, rev = esutil.stat.histogram(ipring, rev=True)

        gdpix, = np.where(h >= self.config.referencePixelizationMinStars)

        self.log.info("Found %d pixels (nside=%d) with at least %d good stars" %
                      (gdpix.size,
                       self.config.referencePixelizationNside,
                       self.config.referencePixelizationMinStars))

        if gdpix.size < self.config.referencePixelizationNPixels:
            self.log.warn("Found fewer good pixels (%d) than preferred in configuration (%d)" %
                          (gdpix.size, self.config.referencePixelizationNPixels))
        else:
            # Sample out the pixels we want to use
            gdpix = np.random.choice(gdpix, size=self.config.referencePixelizationNPixels, replace=False)

        results = np.zeros(gdpix.size, dtype=[('hpix', 'i4'),
                                              ('nstar', 'i4', len(self.bands)),
                                              ('nmatch', 'i4', len(self.bands)),
                                              ('zp', 'f4', len(self.bands)),
                                              ('zpErr', 'f4', len(self.bands))])
        results['hpix'] = ipring[rev[rev[gdpix]]]

        # We need a boolean index to deal with catalogs...
        selected = np.zeros(len(stars), dtype=np.bool)

        refFluxFields = [None] * len(self.bands)

        for p, pix in enumerate(gdpix):
            i1a = rev[rev[pix]: rev[pix + 1]]

            # the stars afwTable can only be indexed with boolean arrays,
            # and not numpy index arrays (see DM-16497).  This little trick
            # converts the index array into a boolean array
            selected[:] = False
            selected[i1a] = True

            for b, band in enumerate(self.bands):

                struct = self._computeOffsetOneBand(sourceMapper, badStarKey, b, band, stars,
                                                    selected, refFluxFields)
                results['nstar'][p, b] = len(i1a)
                results['nmatch'][p, b] = len(struct.arrays.refMag)
                results['zp'][p, b] = struct.zp
                results['zpErr'][p, b] = struct.sigma

        # And compute the summary statistics
        offsets = np.zeros(len(self.bands))

        for b, band in enumerate(self.bands):
            # make configurable
            ok, = np.where(results['nmatch'][:, b] >= self.config.referenceMinMatch)
            offsets[b] = np.median(results['zp'][ok, b])
            # use median absolute deviation to estimate Normal sigma
            # see https://en.wikipedia.org/wiki/Median_absolute_deviation
            madSigma = 1.4826 * np.median(np.abs(results['zp'][ok, b] - offsets[b]))
            self.log.info("Reference catalog offset for %s band: %.12f +/- %.12f" %
                          (band, offsets[b], madSigma))

        return offsets

    def _computeOffsetOneBand(self, sourceMapper, badStarKey,
                              b, band, stars, selected, refFluxFields):
        """
        Compute the zeropoint offset between the fgcm stars and the reference
        stars for one pixel in one band

        Parameters
        ----------
        sourceMapper: `afwTable.SchemaMapper`
           Mapper to go from stars to calibratable catalog
        badStarKey: `afwTable.Key`
           Key for the field with bad stars
        b: `int`
           Index of the band in the star catalog
        band: `str`
           Name of band for reference catalog
        stars: `afwTable.SimpleCatalog`
           FGCM standard stars
        selected: `np.array(dtype=np.bool)`
           Boolean array of which stars are in the pixel
        refFluxFields: `list`
           List of names of flux fields for reference catalog
        """

        sourceCat = afwTable.SimpleCatalog(sourceMapper.getOutputSchema())
        sourceCat.reserve(selected.sum())
        sourceCat.extend(stars[selected], mapper=sourceMapper)
        sourceCat['instFlux'] = 10.**(stars['mag_std_noabs'][selected, b] / (-2.5))
        sourceCat['instFluxErr'] = (np.log(10.) / 2.5) * (stars['magErr_std'][selected, b] *
                                                          sourceCat['instFlux'])
        # Make sure we only use stars that have valid measurements
        # (This is perhaps redundant with requirements above that the
        # stars be observed in all bands, but it can't hurt)
        badStar = (stars['mag_std_noabs'][selected, b] > 90.0)
        for rec in sourceCat[badStar]:
            rec.set(badStarKey, True)

        exposure = afwImage.ExposureF()
        exposure.setFilter(afwImage.Filter(band))

        if refFluxFields[b] is None:
            # Need to find the flux field in the reference catalog
            # to work around limitations of DirectMatch in PhotoCal
            ctr = stars[0].getCoord()
            rad = 0.05 * lsst.geom.degrees
            refDataTest = self.refObjLoader.loadSkyCircle(ctr, rad, band)
            refFluxFields[b] = refDataTest.fluxField

        # Make a copy of the config so that we can modify it
        calConfig = copy.copy(self.config.photoCal.value)
        calConfig.match.referenceSelection.signalToNoise.fluxField = refFluxFields[b]
        calConfig.match.referenceSelection.signalToNoise.errField = refFluxFields[b] + 'Err'
        calTask = self.config.photoCal.target(refObjLoader=self.refObjLoader,
                                              config=calConfig,
                                              schema=sourceCat.getSchema())

        struct = calTask.run(exposure, sourceCat)

        return struct

    def _outputStandardStars(self, butler, offsets):
        """
        Output standard stars in indexed reference catalog format.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        offsets: `np.array` of floats
           Per band zeropoint offsets
        """

        self.log.info("Outputting standard stars to %s" % (self.config.datasetConfig.ref_dataset_name))

        # Load the stars (this is the full set of stars, no cuts)
        stars = butler.get('fgcmStandardStars', fgcmcycle=self.useCycle)

        # We determine the conversion from the native units (typically radians) to
        # degrees for the first star.  This allows us to treat coord_ra/coord_dec as
        # numpy arrays rather than Angles, which would we approximately 600x slower.
        # TODO: Fix this after DM-16524 (HtmIndexer.indexPoints should take coords
        # (as Angles) for input
        conv = stars[0]['coord_ra'].asDegrees() / float(stars[0]['coord_ra'])
        indices = np.array(self.indexer.indexPoints(stars['coord_ra'] * conv,
                                                    stars['coord_dec'] * conv))

        formattedCat = self._formatCatalog(stars, offsets)

        # Write the master schema
        dataId = self.indexer.makeDataId('master_schema',
                                         self.config.datasetConfig.ref_dataset_name)
        butler.put(afwTable.SimpleCatalog(formattedCat.schema), 'ref_cat', dataId=dataId)

        # Break up the pixels using a histogram
        h, rev = esutil.stat.histogram(indices, rev=True)
        gd, = np.where(h > 0)
        selected = np.zeros(len(formattedCat), dtype=np.bool)
        for i in gd:
            i1a = rev[rev[i]: rev[i + 1]]

            # the formattedCat afwTable can only be indexed with boolean arrays,
            # and not numpy index arrays (see DM-16497).  This little trick
            # converts the index array into a boolean array
            selected[:] = False
            selected[i1a] = True

            # Write the individual pixel
            dataId = self.indexer.makeDataId(indices[i1a[0]],
                                             self.config.datasetConfig.ref_dataset_name)
            butler.put(formattedCat[selected], 'ref_cat', dataId=dataId)

        # And save the dataset configuration
        dataId = self.indexer.makeDataId(None, self.config.datasetConfig.ref_dataset_name)
        butler.put(self.config.datasetConfig, 'ref_cat_config', dataId=dataId)

        self.log.info("Done outputting standard stars.")

    def _formatCatalog(self, fgcmStarCat, offsets):
        """
        Turn an FGCM-formatted star catalog, applying zeropoint offsets.

        Parameters
        ----------
        fgcmStarCat: `afwTable.SimpleCatalog`
           SimpleCatalog as output by fgcmcal
        offsets: `list` with len(self.bands) entries
           Zeropoint offsets to apply

        Returns
        -------
        formattedCat: `afwTable.SimpleCatalog`
           SimpleCatalog suitable for using as a reference catalog
        """

        sourceMapper = afwTable.SchemaMapper(fgcmStarCat.schema)
        minSchema = LoadIndexedReferenceObjectsTask.makeMinimalSchema(self.bands,
                                                                      addCentroid=False,
                                                                      addIsResolved=True,
                                                                      coordErrDim=0)
        sourceMapper.addMinimalSchema(minSchema)
        for band in self.bands:
            sourceMapper.editOutputSchema().addField('%s_nGood' % (band), type=np.int32)

        formattedCat = afwTable.SimpleCatalog(sourceMapper.getOutputSchema())
        formattedCat.reserve(len(fgcmStarCat))
        formattedCat.extend(fgcmStarCat, mapper=sourceMapper)

        # Note that we don't have to set `resolved` because the default is False

        for b, band in enumerate(self.bands):
            mag = fgcmStarCat['mag_std_noabs'][:, b] + offsets[b]
            # We want fluxes in Jy from calibrated AB magnitudes
            # (after applying offset)
            # TODO: Full implementation of RFC-549 will have all reference
            # catalogs in nJy instead of Jy.
            flux = afwImage.fluxFromABMag(mag)
            fluxErr = afwImage.fluxErrFromABMagErr(fgcmStarCat['magErr_std'][:, b], mag)
            formattedCat['%s_flux' % (band)][:] = flux
            formattedCat['%s_fluxErr' % (band)][:] = fluxErr
            formattedCat['%s_nGood' % (band)][:] = fgcmStarCat['ngood'][:, b]

        return formattedCat

    def _outputZeropoints(self, butler, offsets):
        """
        Output the zeropoints in fgcm_photoCalib format.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        offsets: `np.array`
           Float array of absolute calibration offsets, one for each filter
        """

        self.log.info("Outputting fgcm_photoCalib objects")

        zptCat = butler.get('fgcmZeropoints', fgcmcycle=self.useCycle)
        visitCat = butler.get('fgcmVisitCatalog')

        # Only output those that we have a calibration
        # See fgcmFitCycle._makeZptSchema for flag definitions
        selected = (zptCat['fgcmFlag'] < 16)

        # Get the mapping from filtername to dataId filter name
        filterMapping = {}
        nFound = 0
        for rec in zptCat[selected]:
            if rec['filtername'] in filterMapping:
                continue
            dataId = {self.visitDataRefName: int(rec['visit']),
                      self.ccdDataRefName: int(rec['ccd'])}
            dataRef = butler.dataRef('raw', dataId=dataId)
            filterMapping[rec['filtername']] = dataRef.dataId['filter']
            nFound += 1
            if nFound == len(self.filterMap):
                break

        # Get a mapping from filtername to the offsets
        offsetMapping = {}
        for f in self.filterMap:
            offsetMapping[f] = offsets[self.bands.index(self.filterMap[f])]

        # Get a mapping from "ccd" to the ccd index used for the scaling
        camera = butler.get('camera')
        ccdMapping = {}
        for ccdIndex, detector in enumerate(camera):
            ccdMapping[detector.getId()] = ccdIndex

        # And a mapping to get the flat-field scaling values
        scalingMapping = {}
        for rec in visitCat:
            scalingMapping[rec['visit']] = rec['scaling']

        for rec in zptCat[selected]:

            if self.superStarSubCcd:
                # Spatially varying zeropoint

                scaling = scalingMapping[rec['visit']][ccdMapping[rec['ccd']]]
                photoCalib = self._getChebyshevPhotoCalib(rec['fgcmfZptCheb'],
                                                          rec['fgcmZptErr'],
                                                          rec['fgcmfZptChebXyMax'],
                                                          offsetMapping[rec['filtername']],
                                                          scaling)
            else:
                # Spatially constant zeropoint

                scaling = scalingMapping[rec['visit']][ccdMapping[rec['ccd']]]
                photoCalib = self._getConstantPhotoCalib(rec['fgcmZpt'], rec['fgcmZptErr'],
                                                         offsetMapping[rec['filtername']],
                                                         scaling)

            butler.put(photoCalib, 'fgcm_photoCalib',
                       dataId={self.visitDataRefName: int(rec['visit']),
                               self.ccdDataRefName: int(rec['ccd']),
                               'filter': filterMapping[rec['filtername']]})

        self.log.info("Done outputting fgcm_photoCalib objects")

    def _getChebyshevPhotoCalib(self, coefficients, err, xyMax, offset, scaling):
        """
        Get the PhotoCalib object from a chebyshev polynomial zeropoint.

        Parameters
        ----------
        coefficients: `np.array`
           Flattened array of chebyshev coefficients
        err: `float`
           Error on zeropoint
        xyMax: `list` of length 2
           Maximum x and y of the chebyshev bounding box
        offset: `float`
           Absolute calibration offset
        scaling: `float`
           Flat scaling value from fgcmBuildStars

        Returns
        -------
        photoCalib: `afwImage.PhotoCalib`
        """

        orderPlus1 = self.chebyshevOrder + 1
        pars = np.zeros((orderPlus1, orderPlus1))

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0.0, 0.0),
                               lsst.geom.Point2I(*xyMax))
        # Take the zeropoint, apply the absolute relative calibration offset,
        # and whatever flat-field scaling was applied
        pars[:, :] = (coefficients.reshape(orderPlus1, orderPlus1) *
                      (offset*units.ABmag).to_value(units.nJy) * scaling)

        field = afwMath.ChebyshevBoundedField(bbox, pars)
        calibMean = field.mean()

        calibErr = (np.log(10.) / 2.5) * calibMean * err

        photoCalib = afwImage.PhotoCalib(field, calibErr)

        return photoCalib

    def _getConstantPhotoCalib(self, zeropoint, err, offset, scaling):
        """
        Get the PhotoCalib object from a constant zeropoint.

        Parameters
        ----------
        zeropoint: `float`
           Zeropoint value (mag)
        err: `float`
           Error on zeropoint
        offset: `float`
           Absolute calibration offset
        scaling: `float`
           Flat scaling value from fgcmBuildStars

        Returns
        -------
        photoCalib: `afwImage.PhotoCalib`
        """

        # Take the zeropoint, apply the absolute relative calibration offset,
        # and whatever flat-field scaling was applied

        calibMean = ((zeropoint + offset)*units.ABmag).to_value(units.nJy) * scaling
        calibErr = (np.log(10.) / 2.5) * calibMean * err
        photoCalib = afwImage.PhotoCalib(calibMean, calibErr)

        return photoCalib

    def _outputAtmospheres(self, butler):
        """
        Output the atmospheres.

        Parameters
        ----------
        butler: `lsst.daf.persistence.Butler`
        """

        self.log.info("Outputting atmosphere transmissions")

        # First, we need to grab the look-up table and key info
        lutCat = butler.get('fgcmLookUpTable')

        atmosphereTableName = lutCat[0]['tablename']
        elevation = lutCat[0]['elevation']
        atmLambda = lutCat[0]['atmLambda']
        lutCat = None

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
            except (ValueError, IOError) as e:
                raise RuntimeError("FGCM look-up-table generated with modtran, "
                                   "but modtran not configured to run.") from e

        # Next, we need to grab the atmosphere parameters
        atmCat = butler.get('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        zenith = np.degrees(np.arccos(1. / atmCat['secZenith']))

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
                       dataId={self.visitDataRefName: visit})

        self.log.info("Done outputting atmosphere transmissions")
