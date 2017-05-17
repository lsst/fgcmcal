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

import time

import lsst.fgcm


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
    def _getConfigName(self):
        return None

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



        visitCat = self._fgcmMakeVisitCatalog(butler)

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
            if (butler.datasetExists('src',dataId={'visit':dataset[0],
                                                   'ccd':self.config.referenceCCD})):
                srcVisits.append(dataset[0])

        print("Found all visits in %.2f" % (time.time()-startTime))

        schema = afwTable.Schema()
        schema.addField('visit',type=np.int32,doc="Visit number")
        schema.addField('band',type=str,size=2,doc="Filter band")
        schema.addField('telra',type=np.float64,doc="Pointing RA (deg)")
        schema.addField('teldec',type=np.float64,doc="Pointing Dec (deg)")
        schema.addField('telha',type=np.float64,doc="Pointing Hour Angle (deg)")
        schema.addField('mjd',type=np.float64,doc="MJD of visit")
        schema.addField('exptime',type=np.float32,doc="Exposure time")
        schema.addField('pmb',type=np.float32,doc="Pressure (millibar)")
        schema.addField('fwhm',type=np.float32,doc="Seeing FWHM?")
        schema.addField('deepflag',type=np.int32,doc="Deep observation")

        visitCat = afwTable.BaseCatalog(schema)
        visitCat.table.preallocate(len(srcVisits))

        startTime = time.time()

        # now loop over visits and get the information
        for srcVisit in srcVisits:
            calexp = butler.get('calexp',dataId={'visit':srcVisit,'ccd':refCCD})

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

        print("Found all VisitInfo in %.2f" % (time.time()-startTime))

        # and now persist it
        butler.put(visitCat, 'fgcmVisitCatalog')

        return visitCat
