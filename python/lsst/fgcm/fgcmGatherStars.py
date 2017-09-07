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
import lsst.daf.persistence


import time

__all__ = ['FgcmGatherStarsConfig', 'FgcmGatherStarsTask']

class FgcmGatherStarsConfig(pexConfig.Config):
    """Config for FgcmGatherStarsTask"""

    def setDefaults(self):
        pass

class FgcmGatherStarsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmGatherStarsTask

    """
    @staticmethod
    def getTargetList(parsedCmd):
        """
        Return a list of tuples per visit, each containing dataRefs

        """
        ## check that this works.
        refListDict = {}
        for ref in parsedCmd.id.refList:
            refListDict.setdefault(ref.dataId['visit'], []).append(ref)

        result = [refListDict[visit] for visit in sorted(refListDict.keys())]

        for r in result:
            sys.log.info("Number of thingies is %d" % (len(r)))
        return result

    def __call__(self, args):
        """
        """
        self.log.info("Called taskrunner")
        raise ValueError("kick out here")
        dataRefList, kwargs = args
        butler = kwargs.pop('butler')
        task = self.TaskClass(config=self.config, butler=butler)
        result = task.run(butler, dataRefList)


class FgcmGatherStarsTask(pipeBase.CmdLineTask):
    """
    Gather visits of stars for the FGCM global calibration star building
    """

    ConfigClass = FgcmGatherStarsConfig
    RunnerClass = FgcmGatherStarsRunner
    _DefaultName = "fgcmGatherStars"

    def __init__(self, butler=None, **kwargs):
        """
        Instantiate an FgcmGatherStarsTask.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
          Something about the butler
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="Data ID, e.g. --id visit=6789 (optional)")

        return parser

    # no saving of the config for now
    def _getConfigName(self):
        return None

    # no saving of metadata for now
    def _getMetadataName(self):
        return None


    @pipeBase.timeMethod
    def run(self, butler, dataRefs):
        """
        """

        startTime=time.time()

        visit=dataRefs[0].dataId['visit']

        if (butler.datasetExists('fgcmVisitObservations',visit=visit)):
            # We already have this, and are done
            return

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

        started=False

        for dataRef in dataRefs:
            print("Reading sources from visit %d/ccd %d" %
                  (visit, dataRef.dataId['ccd']))

            sources = ref.get('src',
                              flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

            if not started:
                # get the keys for quicker look-up

                fluxKey = sources.getApFluxKey()
                fluxErrKey = sources.getApFluxErrKey()
                satCenterKey = sources.schema.find('flag_pixel_saturated_center').key
                intCenterKey = sources.schema.find('flag_pixel_interpolated_center').key
                pixEdgeKey = sources.schema.find('flag_pixel_edge').key
                pixCrCenterKey = sources.schema.find('flag_pixel_cr_center').key
                pixBadKey = sources.schema.find('flag_pixel_bad').key
                pixInterpAnyKey = sources.schema.find('flag_pixel_interpolated_any').key
                centroidFlagKey = sources.schema.find('slot_Centroid_flag').key
                apFluxFlagKey = sources.schema.find('slot_ApFlux_flag').key
                deblendNchildKey = sources.schema.find('deblend_nchild').key
                parentKey = sources.schema.find('parent').key
                extKey = sources.schema.find('classification_extendedness').key

                started=True

            magErr = (2.5/np.log(10.)) * (sources[fluxKey] /
                                          sources[fluxErrKey])
            magErr = np.nan_to_num(magErr)

            gdFlag = np.logical_and.reduce([~sources[satCenterKey],
                                             ~sources[intCenterKey],
                                             ~sources[pixEdgeKey],
                                             ~sources[pixCrCenterKey],
                                             ~sources[pixBadKey],
                                             ~sources[pixInterpAnyKey],
                                             ~sources[centroidFlagKey],
                                             ~sources[apFluxFlagKey],
                                             sources[deblendNchildKey] == 0,
                                             sources[parentKey] == 0,
                                             sources[extKey] < 0.5,
                                             np.isfinite(magErr),
                                             magErr > 0.001,
                                             magErr < 0.1])

            tempCat = afwTable.BaseCatalog(fullCatalog.schema)
            tempCat.table.preallocate(gdFlag.sum())
            tempCat.extend(sources[gdFlag], mapper=sourceMapper)
            tempCat['visit'][:] = visit['visit']
            tempCat['ccd'][:] = ccdId
            tempCat['mag'][:] = 25.0 - 2.5*np.log10(sources[fluxKey][gdFlag])
            tempCat['magerr'][:] = magErr[gdFlag]

            fullCatalog.extend(tempCat)

        print("Found %d good star observations for visit %d in %.2f s" %
              (len(fullCatalog), time.time() - startTime))

        butler.put(fullCatalog, 'fgcmVisitObservations',visit=visit)

