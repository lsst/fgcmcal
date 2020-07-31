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
"""Build star observations for input to FGCM.

This task finds all the visits and calexps in a repository (or a subset
based on command line parameters) and extract all the potential calibration
stars for input into fgcm.  This task additionally uses fgcm to match
star observations into unique stars, and performs as much cleaning of
the input catalog as possible.
"""

import time

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from .fgcmBuildStarsBase import FgcmBuildStarsConfigBase, FgcmBuildStarsRunner, FgcmBuildStarsBaseTask
from .utilities import computeApproxPixelAreaFields

__all__ = ['FgcmBuildStarsConfig', 'FgcmBuildStarsTask']


class FgcmBuildStarsConfig(FgcmBuildStarsConfigBase):
    """Config for FgcmBuildStarsTask"""

    referenceCCD = pexConfig.Field(
        doc="Reference CCD for scanning visits",
        dtype=int,
        default=13,
    )
    checkAllCcds = pexConfig.Field(
        doc=("Check repo for all CCDs for each visit specified.  To be used when the "
             "full set of ids (visit/ccd) are not specified on the command line.  For "
             "Gen2, specifying one ccd and setting checkAllCcds=True is significantly "
             "faster than the alternatives."),
        dtype=bool,
        default=True,
    )

    def setDefaults(self):
        super().setDefaults()

        sourceSelector = self.sourceSelector["science"]

        # The names here correspond to raw src catalogs, which differ
        # from the post-transformed sourceTable_visit catalogs.
        # Therefore, field and flag names cannot be easily
        # derived from the base config class.
        fluxFlagName = self.instFluxField[0: -len('instFlux')] + 'flag'
        sourceSelector.flags.bad = ['base_PixelFlags_flag_edge',
                                    'base_PixelFlags_flag_interpolatedCenter',
                                    'base_PixelFlags_flag_saturatedCenter',
                                    'base_PixelFlags_flag_crCenter',
                                    'base_PixelFlags_flag_bad',
                                    'base_PixelFlags_flag_interpolated',
                                    'base_PixelFlags_flag_saturated',
                                    'slot_Centroid_flag',
                                    fluxFlagName]

        if self.doSubtractLocalBackground:
            localBackgroundFlagName = self.localBackgroundFluxField[0: -len('instFlux')] + 'flag'
            sourceSelector.flags.bad.append(localBackgroundFlagName)

        sourceSelector.signalToNoise.fluxField = self.instFluxField
        sourceSelector.signalToNoise.errField = self.instFluxField + 'Err'


class FgcmBuildStarsTask(FgcmBuildStarsBaseTask):
    """
    Build stars for the FGCM global calibration, using src catalogs.
    """
    ConfigClass = FgcmBuildStarsConfig
    RunnerClass = FgcmBuildStarsRunner
    _DefaultName = "fgcmBuildStars"

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "src", help="Data ID, e.g. --id visit=6789")

        return parser

    def findAndGroupDataRefs(self, butler, dataRefs):
        self.log.info("Grouping dataRefs by %s" % (self.config.visitDataRefName))

        camera = butler.get('camera')

        ccdIds = []
        for detector in camera:
            ccdIds.append(detector.getId())

        # TODO: related to DM-13730, this dance of looking for source visits
        # will be unnecessary with Gen3 Butler.  This should be part of
        # DM-13730.

        nVisits = 0

        groupedDataRefs = {}
        for dataRef in dataRefs:
            visit = dataRef.dataId[self.config.visitDataRefName]
            # If we don't have the dataset, just continue
            if not dataRef.datasetExists(datasetType='src'):
                continue
            # If we need to check all ccds, do it here
            if self.config.checkAllCcds:
                if visit in groupedDataRefs:
                    # We already have found this visit
                    continue
                dataId = dataRef.dataId.copy()
                # For each ccd we must check that a valid source catalog exists.
                for ccdId in ccdIds:
                    dataId[self.config.ccdDataRefName] = ccdId
                    if butler.datasetExists('src', dataId=dataId):
                        goodDataRef = butler.dataRef('src', dataId=dataId)
                        if visit in groupedDataRefs:
                            if (goodDataRef.dataId[self.config.ccdDataRefName] not in
                               [d.dataId[self.config.ccdDataRefName] for d in groupedDataRefs[visit]]):
                                groupedDataRefs[visit].append(goodDataRef)
                        else:
                            # This is a new visit
                            nVisits += 1
                            groupedDataRefs[visit] = [goodDataRef]
            else:
                # We have already confirmed that the dataset exists, so no need
                # to check here.
                if visit in groupedDataRefs:
                    if (dataRef.dataId[self.config.ccdDataRefName] not in
                       [d.dataId[self.config.ccdDataRefName] for d in groupedDataRefs[visit]]):
                        groupedDataRefs[visit].append(dataRef)
                else:
                    # This is a new visit
                    nVisits += 1
                    groupedDataRefs[visit] = [dataRef]

            if (nVisits % 100) == 0 and nVisits > 0:
                self.log.info("Found %d unique %ss..." % (nVisits,
                                                          self.config.visitDataRefName))

        self.log.info("Found %d unique %ss total." % (nVisits,
                                                      self.config.visitDataRefName))

        # Put them in ccd order, with the reference ccd first (if available)
        def ccdSorter(dataRef):
            ccdId = dataRef.dataId[self.config.ccdDataRefName]
            if ccdId == self.config.referenceCCD:
                return -100
            else:
                return ccdId

        # If we did not check all ccds, put them in ccd order
        if not self.config.checkAllCcds:
            for visit in groupedDataRefs:
                groupedDataRefs[visit] = sorted(groupedDataRefs[visit], key=ccdSorter)

        return groupedDataRefs

    def fgcmMakeAllStarObservations(self, groupedDataRefs, visitCat,
                                    calibFluxApertureRadius=None,
                                    visitCatDataRef=None,
                                    starObsDataRef=None,
                                    inStarObsCat=None):
        startTime = time.time()

        # If both dataRefs are None, then we assume the caller does not
        # want to store checkpoint files.  If both are set, we will
        # do checkpoint files.  And if only one is set, this is potentially
        # unintentional and we will warn.
        if (visitCatDataRef is not None and starObsDataRef is None or
           visitCatDataRef is None and starObsDataRef is not None):
            self.log.warn("Only one of visitCatDataRef and starObsDataRef are set, so "
                          "no checkpoint files will be persisted.")

        if self.config.doSubtractLocalBackground and calibFluxApertureRadius is None:
            raise RuntimeError("Must set calibFluxApertureRadius if doSubtractLocalBackground is True.")

        # create our source schema.  Use the first valid dataRef
        dataRef = groupedDataRefs[list(groupedDataRefs.keys())[0]][0]
        sourceSchema = dataRef.get('src_schema', immediate=True).schema

        # Construct a mapping from ccd number to index
        camera = dataRef.get('camera')
        ccdMapping = {}
        for ccdIndex, detector in enumerate(camera):
            ccdMapping[detector.getId()] = ccdIndex

        approxPixelAreaFields = computeApproxPixelAreaFields(camera)

        sourceMapper = self._makeSourceMapper(sourceSchema)

        # We also have a temporary catalog that will accumulate aperture measurements
        aperMapper = self._makeAperMapper(sourceSchema)

        outputSchema = sourceMapper.getOutputSchema()

        if inStarObsCat is not None:
            fullCatalog = inStarObsCat
            comp1 = fullCatalog.schema.compare(outputSchema, outputSchema.EQUAL_KEYS)
            comp2 = fullCatalog.schema.compare(outputSchema, outputSchema.EQUAL_NAMES)
            if not comp1 or not comp2:
                raise RuntimeError("Existing fgcmStarObservations file found with mismatched schema.")
        else:
            fullCatalog = afwTable.BaseCatalog(outputSchema)

        # FGCM will provide relative calibration for the flux in config.instFluxField

        instFluxKey = sourceSchema[self.config.instFluxField].asKey()
        instFluxErrKey = sourceSchema[self.config.instFluxField + 'Err'].asKey()
        visitKey = outputSchema['visit'].asKey()
        ccdKey = outputSchema['ccd'].asKey()
        instMagKey = outputSchema['instMag'].asKey()
        instMagErrKey = outputSchema['instMagErr'].asKey()
        deltaMagBkgKey = outputSchema['deltaMagBkg'].asKey()

        # Prepare local background if desired
        if self.config.doSubtractLocalBackground:
            localBackgroundFluxKey = sourceSchema[self.config.localBackgroundFluxField].asKey()
            localBackgroundArea = np.pi*calibFluxApertureRadius**2.

        aperOutputSchema = aperMapper.getOutputSchema()

        instFluxAperInKey = sourceSchema[self.config.apertureInnerInstFluxField].asKey()
        instFluxErrAperInKey = sourceSchema[self.config.apertureInnerInstFluxField + 'Err'].asKey()
        instFluxAperOutKey = sourceSchema[self.config.apertureOuterInstFluxField].asKey()
        instFluxErrAperOutKey = sourceSchema[self.config.apertureOuterInstFluxField + 'Err'].asKey()
        instMagInKey = aperOutputSchema['instMag_aper_inner'].asKey()
        instMagErrInKey = aperOutputSchema['instMagErr_aper_inner'].asKey()
        instMagOutKey = aperOutputSchema['instMag_aper_outer'].asKey()
        instMagErrOutKey = aperOutputSchema['instMagErr_aper_outer'].asKey()

        k = 2.5/np.log(10.)

        # loop over visits
        for ctr, visit in enumerate(visitCat):
            if visit['sources_read']:
                continue

            expTime = visit['exptime']

            nStarInVisit = 0

            # Reset the aperture catalog (per visit)
            aperVisitCatalog = afwTable.BaseCatalog(aperOutputSchema)

            for dataRef in groupedDataRefs[visit['visit']]:

                ccdId = dataRef.dataId[self.config.ccdDataRefName]

                sources = dataRef.get(datasetType='src', flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
                goodSrc = self.sourceSelector.selectSources(sources)

                tempCat = afwTable.BaseCatalog(fullCatalog.schema)
                tempCat.reserve(goodSrc.selected.sum())
                tempCat.extend(sources[goodSrc.selected], mapper=sourceMapper)
                tempCat[visitKey][:] = visit['visit']
                tempCat[ccdKey][:] = ccdId

                # Compute "instrumental magnitude" by scaling flux with exposure time.
                scaledInstFlux = (sources[instFluxKey][goodSrc.selected] *
                                  visit['scaling'][ccdMapping[ccdId]])
                tempCat[instMagKey][:] = (-2.5*np.log10(scaledInstFlux) + 2.5*np.log10(expTime))

                # Compute the change in magnitude from the background offset
                if self.config.doSubtractLocalBackground:
                    # At the moment we only adjust the flux and not the flux
                    # error by the background because the error on
                    # base_LocalBackground_instFlux is the rms error in the
                    # background annulus, not the error on the mean in the
                    # background estimate (which is much smaller, by sqrt(n)
                    # pixels used to estimate the background, which we do not
                    # have access to in this task).  In the default settings,
                    # the annulus is sufficiently large such that these
                    # additional errors are are negligibly small (much less
                    # than a mmag in quadrature).

                    localBackground = localBackgroundArea*sources[localBackgroundFluxKey]

                    # This is the difference between the mag with background correction
                    # and the mag without background correction.
                    tempCat[deltaMagBkgKey][:] = (-2.5*np.log10(sources[instFluxKey][goodSrc.selected] -
                                                                localBackground[goodSrc.selected]) -
                                                  -2.5*np.log10(sources[instFluxKey][goodSrc.selected]))
                else:
                    tempCat[deltaMagBkgKey][:] = 0.0

                # Compute instMagErr from instFluxErr/instFlux, any scaling
                # will cancel out.

                tempCat[instMagErrKey][:] = k*(sources[instFluxErrKey][goodSrc.selected] /
                                               sources[instFluxKey][goodSrc.selected])

                # Compute the jacobian from an approximate PixelAreaBoundedField
                tempCat['jacobian'] = approxPixelAreaFields[ccdId].evaluate(tempCat['x'],
                                                                            tempCat['y'])

                # Apply the jacobian if configured
                if self.config.doApplyWcsJacobian:
                    tempCat[instMagKey][:] -= 2.5*np.log10(tempCat['jacobian'][:])

                fullCatalog.extend(tempCat)

                # And the aperture information
                # This does not need the jacobian because it is all locally relative
                tempAperCat = afwTable.BaseCatalog(aperVisitCatalog.schema)
                tempAperCat.reserve(goodSrc.selected.sum())
                tempAperCat.extend(sources[goodSrc.selected], mapper=aperMapper)

                with np.warnings.catch_warnings():
                    # Ignore warnings, we will filter infinities and
                    # nans below.
                    np.warnings.simplefilter("ignore")

                    tempAperCat[instMagInKey][:] = -2.5*np.log10(
                        sources[instFluxAperInKey][goodSrc.selected])
                    tempAperCat[instMagErrInKey][:] = k*(
                        sources[instFluxErrAperInKey][goodSrc.selected] /
                        sources[instFluxAperInKey][goodSrc.selected])
                    tempAperCat[instMagOutKey][:] = -2.5*np.log10(
                        sources[instFluxAperOutKey][goodSrc.selected])
                    tempAperCat[instMagErrOutKey][:] = k*(
                        sources[instFluxErrAperOutKey][goodSrc.selected] /
                        sources[instFluxAperOutKey][goodSrc.selected])

                aperVisitCatalog.extend(tempAperCat)

                nStarInVisit += len(tempCat)

            # Compute the median delta-aper
            if not aperVisitCatalog.isContiguous():
                aperVisitCatalog = aperVisitCatalog.copy(deep=True)

            instMagIn = aperVisitCatalog[instMagInKey]
            instMagErrIn = aperVisitCatalog[instMagErrInKey]
            instMagOut = aperVisitCatalog[instMagOutKey]
            instMagErrOut = aperVisitCatalog[instMagErrOutKey]

            ok = (np.isfinite(instMagIn) & np.isfinite(instMagErrIn) &
                  np.isfinite(instMagOut) & np.isfinite(instMagErrOut))

            visit['deltaAper'] = np.median(instMagIn[ok] - instMagOut[ok])
            visit['sources_read'] = True

            self.log.info("  Found %d good stars in visit %d (deltaAper = %.3f)" %
                          (nStarInVisit, visit['visit'], visit['deltaAper']))

            if ((ctr % self.config.nVisitsPerCheckpoint) == 0 and
               starObsDataRef is not None and visitCatDataRef is not None):
                # We need to persist both the stars and the visit catalog which gets
                # additional metadata from each visit.
                starObsDataRef.put(fullCatalog)
                visitCatDataRef.put(visitCat)

        self.log.info("Found all good star observations in %.2f s" %
                      (time.time() - startTime))

        return fullCatalog

    def _makeAperMapper(self, sourceSchema):
        """
        Make a schema mapper for fgcm aperture measurements

        Parameters
        ----------
        sourceSchema: `afwTable.Schema`
           Default source schema from the butler

        Returns
        -------
        aperMapper: `afwTable.schemaMapper`
           Mapper to the FGCM aperture schema
        """

        aperMapper = afwTable.SchemaMapper(sourceSchema)
        aperMapper.addMapping(sourceSchema['coord_ra'].asKey(), 'ra')
        aperMapper.addMapping(sourceSchema['coord_dec'].asKey(), 'dec')
        aperMapper.editOutputSchema().addField('instMag_aper_inner', type=np.float64,
                                               doc="Magnitude at inner aperture")
        aperMapper.editOutputSchema().addField('instMagErr_aper_inner', type=np.float64,
                                               doc="Magnitude error at inner aperture")
        aperMapper.editOutputSchema().addField('instMag_aper_outer', type=np.float64,
                                               doc="Magnitude at outer aperture")
        aperMapper.editOutputSchema().addField('instMagErr_aper_outer', type=np.float64,
                                               doc="Magnitude error at outer aperture")

        return aperMapper
