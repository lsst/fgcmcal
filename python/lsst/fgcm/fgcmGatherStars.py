# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import sys
import traceback

import numpy as np
import healpy as hp

import lsst.utils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pex.exceptions as pexExceptions
import lsst.afw.table as afwTable
from lsst.daf.base.dateTime import DateTime
import lsst.afw.geom as afwGeom
import lsst.daf.persistence.butlerExceptions as butlerExceptions
import lsst.daf.persistence
import lsst.afw.image as afwImage


import time

__all__ = ['FgcmGatherStarsConfig', 'FgcmGatherStarsTask']

class FgcmGatherStarsConfig(pexConfig.Config):
    """Config for FgcmGatherStarsTask"""

    nSide = pexConfig.Field(
        doc="healpix nside to group observations",
        dtype=int,
        default=8,
        )

    def setDefaults(self):
        pass

class FgcmGatherStarsRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmGatherStarsTask

    """
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """
        Return a list of tuples per visit, each containing dataRefs

        """
        kwargs['butler'] = parsedCmd.butler

        ## check that this works.
        refListDict = {}
        for ref in parsedCmd.id.refList:
            refListDict.setdefault(ref.dataId['visit'], []).append(ref)

        result = [(refListDict[visit], kwargs) for visit in sorted(refListDict.keys())]

        #for ref in parsedCmd.id.refList:
        #    md = dataRef.get("calexp_md", immediate=True)
        #    wcs = afwImage.makeWcs(md)
        #    center = wcs.pixelToSky(md.get("NAXIS1")/2., md.get("NAXIS2")/2.)
        #    theta = np.pi/2. - center.getDec().asRadians()
        #    phi = center.getRa().asRadians()
        #    ipring = hp.ang2pix(self.config.nside, theta, phi)

        #    refListDict.setdefault(ipring, []).append(ref)

        #result = [(refListDict[ipring], kwargs) for ipring in sorted(refListDict.keys())]

        return result

    def __call__(self, args):
        """
        """
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

        self.log.info("Working on visit %d with %d ccds" % (visit, len(dataRefs)))
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
        outputStarted=False
        #starsSelector = FgcmGatherStarsSelector()

        for dataRef in dataRefs:
            self.log.info("Reading sources from visit %d/ccd %d" %
                          (visit, dataRef.dataId['ccd']))

            sources = dataRef.get('src',
                                  flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

            #cutSources = starsSelector.selectSources(sources)

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

                #outputSchema = sourceMapper.getOutputSchema()
                #visitKey = outputSchema.find('visit')
                #ccdKey = outputSchema.find('ccd')
                #magKey = outputSchema.find('mag')
                #magErrKey = outputSchema.find('magerr')

                started=True

            #magErr = (2.5/np.log(10.)) * (sources.getApFluxErr() /
            #                              sources.getApFlux())
            #magErr = np.nan_to_num(magErr)

            #gdFlag = np.logical_and.reduce([~sources.get('flag_pixel_saturated_center'),
            #                                 ~sources.get('flag_pixel_interpolated_center'),
            #                                 ~sources.get('flag_pixel_edge'),
            #                                 ~sources.get('flag_pixel_cr_center'),
            #                                 ~sources.get('flag_pixel_bad'),
            #                                 ~sources.get('flag_pixel_interpolated_any'),
            #                                 ~sources.get('slot_Centroid_flag'),
            #                                 ~sources.get('slot_ApFlux_flag'),
            #                                 sources.get('deblend_nchild') == 0,
            #                                 sources.get('parent') == 0,
            #                                 sources.get('classification_extendedness') < 0.5,
            #                                 np.isfinite(magErr),
            #                                 magErr > 0.001,
             #                                magErr < 0.1])

            magErr = (2.5/np.log(10.)) * (sources.get(fluxKey) /
                                          sources.get(fluxErrKey))
            magErr = np.nan_to_num(magErr)

            gdFlag = np.logical_and.reduce([~sources.get(satCenterKey),
                                             ~sources.get(intCenterKey),
                                             ~sources.get(pixEdgeKey),
                                             ~sources.get(pixCrCenterKey),
                                             ~sources.get(pixBadKey),
                                             ~sources.get(pixInterpAnyKey),
                                             ~sources.get(centroidFlagKey),
                                             ~sources.get(apFluxFlagKey),
                                             sources.get(deblendNchildKey) == 0,
                                             sources.get(parentKey) == 0,
                                             sources.get(extKey) < 0.5,
                                             np.isfinite(magErr),
                                             magErr > 0.001,
                                             magErr < 0.1])

            tempCat = afwTable.BaseCatalog(fullCatalog.schema)
            tempCat.table.preallocate(gdFlag.sum())
            tempCat.extend(sources[gdFlag], mapper=sourceMapper)
            #tempCat.get('visit')[:] = visit
            #tempCat.get('ccd')[:] = dataRef.dataId['ccd']
            #tempCat.get('mag')[:] = 25.0 - 2.5*np.log10(sources.getApFlux()[gdFlag])
            #tempCat.get('magerr')[:] = magErr[gdFlag]

            if (not outputStarted):
                visitKey = tempCat.schema.find('visit')
                ccdKey = tempCat.schema.find('ccd')
                magKey = tempCat.schema.find('mag')
                magErrKey = tempCat.schema.find('magerr')

                outputStarted = True

            tempCat.get(visitKey)[:] = visit
            tempCat.get(ccdKey)[:] = dataRef.dataId['ccd']
            tempCat.get(magKey)[:] = 25.0 - 2.5*np.log10(sources.getApFlux()[gdFlag])
            tempCat.get(magErrKey)[:] = magErr[gdFlag]


            fullCatalog.extend(tempCat)

        butler.put(fullCatalog, 'fgcmVisitObservations',visit=visit)

        self.log.info("Found %d good star observations for visit %d in %.2f s" %
                      (len(fullCatalog), visit, time.time() - startTime))


class FgcmGatherStarsSelector(object):
    """
    """

    def __init__(self):
        self.badFlags = ['flag_pixel_saturated_center',
                         'flag_pixel_interpolated_center',
                         'flag_pixel_edge',
                         'flag_pixel_cr_center',
                         'flag_pixel_bad',
                         'flag_pixel_interpolated_any',
                         'slot_Centroid_flag',
                         'slot_ApFlux_flag']
        self.magConst = (2.5/np.log(10.))
        # will want to make these configurable
        self.minMagErr = 0.001
        self.maxMagErr = 0.1

    def _isBad(self, source):
        return any(source.get(flag) for flag in self.badFlags)

    def selectSources(self, sourceCat):
        result = afwTable.SourceCatalog(sourceCat.table)

        for source in sourceCat:
            magErr = self.magConst * (source.getApFluxErr() /
                                      source.getApFlux())

            if (not self._isBad(source) and
                source['deblend_nchild'] == 0 and
                source['parent'] == 0 and
                source['classification_extendedness'] < 0.5 and
                np.isfinite(magErr) and
                magErr > 0.001 and
                magErr < 0.1):
                result.append(source)
        # return a copy for contiguous access
        return result.copy(deep=True)
