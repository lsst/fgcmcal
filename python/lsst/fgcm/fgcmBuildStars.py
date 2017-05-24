# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import sys
import traceback

import numpy as np

import lsst.utils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pex.exceptions as pexExceptions
import lsst.afw.table as afwTable
from lsst.daf.base.dateTime import DateTime
import lsst.afw.geom as afwGeom
import lsst.daf.persistence.butlerExceptions as butlerExceptions

import time

import lsst.fgcm as lsstFgcm

import fgcm


__all__ = ['FgcmBuildStarsConfig','FgcmBuildStarsTask']

class FgcmBuildStarsConfig(pexConfig.Config):
    """Config for FgcmBuildStarsTask"""

    minPerBand = pexConfig.Field(
        doc="Minimum observations per band",
        dtype=int,
        default=2,
        )
    matchRadius = pexConfig.Field(
        doc="Match radius (arcseconds)",
        dtype=float,
        default=1.0,
        )
    isolationRadius = pexConfig.Field(
        doc="Isolation radius (arcseconds)",
        dtype=float,
        default=2.0,
        )
    densityCutNside = pexConfig.Field(
        doc="Density cut healpix nside",
        dtype=int,
        default=128,
        )
    densityCutMaxPerPixel = pexConfig.Field(
        doc="Density cut number of stars per pixel",
        dtype=int,
        default=1000,
        )
    zeropointDefault = pexConfig.Field(
        doc="Zeropoint default (arbitrary?)",
        dtype=float,
        default=25.0,
        )
    bands = pexConfig.ListField(
        doc="Bands to run calibration",
        dtype=str,
        default=("NO_DATA",),
        )
    requiredFlag = pexConfig.ListField(
        doc="Flag for required bands",
        dtype=int,
        default=(0,),
        )
    referenceCCD = pexConfig.Field(
        doc="Reference CCD for scanning visits",
        dtype=int,
        default=13,
        )

    def setDefaults(self):
        pass

class FgcmBuildStarsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmBuildStarsTask

    """

    #TaskClass = FgcmBuildStarsTask

    # only need a single butler instance to run on
    @staticmethod
    def getTargetList(parsedCmd):
        return [parsedCmd.butler]

    def precall(self, parsedCmd):
        return True

    def __call__(self, butler):
        print("In taskrunner __call__")
        task = self.TaskClass(config=self.config, log=self.log)
        if self.doRaise:
            results = task.run(butler)
        else:
            try:
                results = task.run(butler)
            except Exception as e:
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)
        print("Done with taskrunner")
        if self.doReturnResults:
            return results

    # turn off any multiprocessing

    def run(self, parsedCmd):
        """ runs the task, but doesn't do multiprocessing"""

        print("in taskrunner run")
        resultList = []

        if self.precall(parsedCmd):
            profileName = parsedCmd.profile if hasattr(parsedCmd, "profile") else None
            log = parsedCmd.log
            #targetList = self.getTargetList(parsedCmd)
            # I think targetList should just have 1 entry?
            #if len(targetList) > 0:
            #    with profile(profileName, log):
            #        resultList = list(map(self, targetList))
            #else:
            #    log.warn("not running the task because there is no data to process")
            targetList = self.getTargetList(parsedCmd)
            # make sure that we only get 1
            resultList = self(targetList[0])

        return resultList

    # and override getTargetList ... want just one?
    #@staticmethod

class FgcmBuildStarsTask(pipeBase.CmdLineTask):
    """
    Build stars for the FGCM global calibration
    """

    ConfigClass = FgcmBuildStarsConfig
    RunnerClass = FgcmBuildStarsRunner
    _DefaultName = "fgcmBuildStars"

    def __init__(self, butler=None, **kwargs):
        """
        Instantiate an FgcmBuildStarsTask.

        Parameters
        ----------
        butler : lsst.daf.persistence.Butler
          Something about the butler
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)

        return parser

    # no saving of the config for now
    #def _getConfigName(self):
    #    return None

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def run(self, butler):
        """
        Cross-match and make star list for FGCM

        Parameters
        ----------
        butler:  a butler.  try to run all from the rerun?  is that crazy?
        dataRefs: list of lsst.daf.persistence.ButlerDataRef
            List of data references to the exposures to be fit

        Returns
        -------
        pipe.base.Struct
            struct containing:
            * dataRefs: the provided data references consolidated
            (others)
        """

        #if len(dataRefs) == 0:
        #    raise ValueError("Need a list of data references!")

        print("Run")
        print("min obs: %d" % (self.config.minPerBand))
        print("bands:")
        print(self.config.bands)


        # make the visit catalog if necessary
        #  question: what's the propper clobber interface?
        if (butler.datasetExists('fgcmVisitCatalog')):
            visitCat = butler.get('fgcmVisitCatalog')
        else:
            # we need to build visitCat
            visitCat = self._fgcmMakeVisitCatalog(butler)

        # and compile all the stars
        #  this will put this dataset out.
        if (not butler.datasetExists('fgcmStarObservations')):
            self._fgcmMakeAllStarObservations(butler, visitCat)

        if (not butler.datasetExists('fgcmStarIds') or
            not butler.datasetExists('fgcmStarIndices')):
            self._fgcmMatchStars(butler, visitCat)

        # next: need to get a list of source catalogs, etc.
        #  just a few would be fine.  Then I could see the formatting of things.
        # how to get into interactive as well?

        # a-ha!
        # first, need to compile all the visits
        # second, need to compile all the observations

        return None

    def _fgcmMakeVisitCatalog(self,butler):
        """
        """

        # check to see if this already exists...
        startTime = time.time()

        allVisits = butler.queryMetadata('src',
                                         format=['visit','filter'],
                                         dataId={'CCD':self.config.referenceCCD})

        srcVisits = []
        for dataset in allVisits:
            if (butler.datasetExists('src', dataId={'visit':dataset[0],
                                                    'ccd':self.config.referenceCCD})):
                srcVisits.append(dataset[0])

        print("Found all visits in %.2f s" % (time.time()-startTime))

        schema = afwTable.Schema()
        schema.addField('visit', type=np.int32, doc="Visit number")
        schema.addField('band', type=str,size=2, doc="Filter band")
        schema.addField('telra', type=np.float64, doc="Pointing RA (deg)")
        schema.addField('teldec', type=np.float64, doc="Pointing Dec (deg)")
        schema.addField('telha', type=np.float64, doc="Pointing Hour Angle (deg)")
        schema.addField('mjd', type=np.float64, doc="MJD of visit")
        schema.addField('exptime', type=np.float32, doc="Exposure time")
        schema.addField('pmb', type=np.float32, doc="Pressure (millibar)")
        schema.addField('fwhm', type=np.float32, doc="Seeing FWHM?")
        schema.addField('deepflag', type=np.int32, doc="Deep observation")

        visitCat = afwTable.BaseCatalog(schema)
        visitCat.table.preallocate(len(srcVisits))

        startTime = time.time()
        # reading in a small bbox is marginally faster in the scan
        bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.PointI(1, 1))

        # now loop over visits and get the information
        for srcVisit in srcVisits:
            calexp = butler.get('calexp_sub', dataId={'visit':srcVisit,
                                                      'ccd':self.config.referenceCCD},
                                bbox=bbox)

            visitInfo = calexp.getInfo().getVisitInfo()

            rec=visitCat.addNew()
            rec['visit'] = srcVisit
            rec['band'] = calexp.getInfo().getFilter().getName()
            radec = visitInfo.getBoresightRaDec()
            rec['telra'] = radec.getRa().asDegrees()
            rec['teldec'] = radec.getDec().asDegrees()
            rec['telha'] = visitInfo.getBoresightHourAngle().asDegrees()
            rec['mjd'] = visitInfo.getDate().get(system=DateTime.MJD)
            rec['exptime'] = visitInfo.getExposureTime()
            # convert from Pa to millibar
            rec['pmb'] = visitInfo.getWeather().getAirPressure() / 100
            rec['fwhm'] = 0.0
            rec['deepflag'] = 0

        print("Found all VisitInfo in %.2f s" % (time.time() - startTime))

        # and now persist it
        butler.put(visitCat, 'fgcmVisitCatalog')

        return visitCat

    def _fgcmMakeAllStarObservations(self, butler, visitCat):
        """
        """

        startTime=time.time()

        # create our source schema
        sourceSchema = butler.get('src_schema', immediate=True).schema

        # create a mapper to the preferred output
        sourceMapper = afwTable.SchemaMapper(sourceSchema)

        # map to ra/dec
        sourceMapper.addMapping(sourceSchema.find('coord_ra').key, 'ra')
        sourceMapper.addMapping(sourceSchema.find('coord_dec').key, 'dec')

        # and add the fields we want
        sourceMapper.editOutputSchema().addField(
            "visit", type=np.int32, doc="Visit number")
        sourceMapper.editOutputSchema().addField(
            "ccd", type=np.int32, doc="CCD number")
        sourceMapper.editOutputSchema().addField(
            "mag", type=np.float32, doc="Raw magnitude")
        sourceMapper.editOutputSchema().addField(
            "magerr", type=np.float32, doc="Raw magnitude error")

        # create the stub of the full catalog
        fullCatalog = afwTable.BaseCatalog(sourceMapper.getOutputSchema())

        # we need to know the ccds...
        camera = butler.get('camera')

        # loop over visits
        for visit in visitCat:
            print("Reading sources from visit %d" % (visit['visit']))
            # loop over CCDs
            for detector in camera:
                ccdId = detector.getId()

                # get the dataref
                ref = butler.dataRef('raw', dataId={'visit':visit['visit'],
                                                    'ccd':ccdId})
                try:
                    sources = ref.get('src',
                                      flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                except butlerExceptions.NoResults:
                    # this ccd does not exist.  That's fine.
                    continue

                # based on ApFlux.  Maybe make this configurable
                magErr = (2.5/np.log(10.)) * (sources['slot_ApFlux_fluxSigma'] /
                                              sources['slot_ApFlux_flux'])
                magErr = np.nan_to_num(magErr)

                # general flag, child/parent/etc cuts
                # will want to make magErr range configurable.
                gdFlag = np.logical_and.reduce([~sources['base_PixelFlags_flag_saturatedCenter'],
                                     ~sources['base_PixelFlags_flag_interpolatedCenter'],
                                     ~sources['base_PixelFlags_flag_edge'],
                                     ~sources['base_PixelFlags_flag_crCenter'],
                                     ~sources['base_PixelFlags_flag_bad'],
                                     ~sources['base_PixelFlags_flag_interpolated'],
                                     ~sources['slot_Centroid_flag'],
                                     ~sources['slot_Centroid_flag_edge'],
                                     ~sources['slot_ApFlux_flag'],
                                     ~sources['base_ClassificationExtendedness_flag'],
                                     sources['deblend_nChild'] == 0,
                                     sources['parent'] == 0,
                                     sources['base_ClassificationExtendedness_value'] < 0.5,
                                     np.isfinite(magErr),
                                     magErr > 0.001,
                                     magErr < 0.1])

                tempCat = afwTable.BaseCatalog(fullCatalog.schema)
                tempCat.table.preallocate(gdFlag.sum())
                tempCat.extend(sources[gdFlag], mapper=sourceMapper)
                tempCat['visit'][:] = visit['visit']
                tempCat['ccd'][:] = ccdId
                tempCat['mag'][:] = 25.0 - 2.5*np.log10(sources['slot_ApFlux_flux'][gdFlag])
                tempCat['magerr'][:] = magErr[gdFlag]

                fullCatalog.extend(tempCat)

        print("Found all good star observations in %.2f s" %
              (time.time() - startTime))

        butler.put(fullCatalog, 'fgcmStarObservations')

        print("Done with all stars in %.2f s" %
              (time.time() - startTime))

    def _fgcmMatchStars(self, butler, visitCat):
        """
        """

        obsCat = butler.get('fgcmStarObservations')

        # get bands into a numpy array...
        visitBands = np.zeros(len(visitCat), dtype='a2')
        for i in xrange(len(visitCat)):
            visitBands[i] = visitCat[i]['band']

        # match to put bands with observations
        visitIndex = np.searchsorted(visitCat['visit'],
                                     obsCat['visit'])

        obsBands = visitBands[visitIndex]

        # make the fgcm starConfig dict
        ## FIXME: make the fgcm configuration dict

        # initialize the FgcmMakeStars object
        fgcmMakeStars = fgcm.FgcmMakeStars(starConfig)

        # make the reference stars
        fgcmMakeStars.makeReferenceStars(obsCat['ra'], obsCat['dec'],
                                         bandArray = obsBands,
                                         bandSelected = False)

        # and match all the stars
        fgcmMakeStars.makeMatchedStars(obsCat['ra'], obsCat['dec'], obsBands)

        # now persist

        # afwTable for objects
        objSchema = afwTable.Schema()
        objSchema.addField('fgcm_id', type=np.int32, doc='FGCM Unique ID')
        objSchema.addField('ra', type=np.float64, doc='Mean object RA')
        objSchema.addField('dec', type=np.float64, doc='Mean object Dec')
        objSchema.addField('obsarrindex', type=np.int32,
                           doc='Index in obsIndexTable for first observation')
        objSchema.addField('nobs', type=np.int32, doc='Total number of observations')

        fgcmStarIdCat = afwTable.BaseCatalog(objSchema)
        fgcmStarIdCat.table.preallocate(fgcmMakeStars.objIndexCat.size)

        fgcmStarIdCat['fgcm_id'][:] = fgcmMakeStars.objIndexCat['FGCM_ID']
        fgcmStarIdCat['ra'][:] = fgcmMakeStars.objIndexCat['RA']
        fgcmStarIdCat['dec'][:] = fgcmMakeStars.objIndexCat['DEC']
        fgcmStarIdCat['obsarrindex'][:] = fgcmMakeStars.objIndexCat['OBSARRINDEX']
        fgcmStarIdCat['nobs'][:] = fgcmMakeStars.objIndexCat['NOBS']

        butler.put(fgcmStarIdCat, 'fgcmStarIds')

        # afwTable for observation indices
        obsSchema = afwTable.Schema()
        obsSchema.addField('obsindex', type=np.int32, doc='Index in observation table')

        fgcmStarIndicesCat = afwTable.BaseCatalog(obsSchema)
        fgcmStarIndicesCat.table.preallocate(fgcmMakeStars.obsIndexCat.size)

        fgcmStarIndicesCat['obsindex'][:] = fgcmMakeStars.obsIndexCat['OBSINDEX']

        butler.put(fgcmStarIndicesCat, 'fgcmStarIndices')

        # and we're done with the stars
