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
"""Build star observations for input to FGCM using sourceTable_visit.

This task finds all the visits and sourceTable_visits in a repository (or a
subset based on command line parameters) and extracts all the potential
calibration stars for input into fgcm.  This task additionally uses fgcm to
match star observations into unique stars, and performs as much cleaning of the
input catalog as possible.
"""

import time

import numpy as np
import collections

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from .fgcmBuildStarsBase import FgcmBuildStarsConfigBase, FgcmBuildStarsRunner, FgcmBuildStarsBaseTask
from .utilities import computeApproxPixelAreaFields

__all__ = ['FgcmBuildStarsTableConfig', 'FgcmBuildStarsTableTask']


class FgcmBuildStarsTableConfig(FgcmBuildStarsConfigBase):
    """Config for FgcmBuildStarsTableTask"""

    referenceCCD = pexConfig.Field(
        doc="Reference CCD for checking PSF and background",
        dtype=int,
        default=40,
    )

    def setDefaults(self):
        super().setDefaults()

        # The names here correspond to the post-transformed
        # sourceTable_visit catalogs, which differ from the raw src
        # catalogs.  Therefore, all field and flag names cannot
        # be derived from the base config class.
        self.instFluxField = 'ApFlux_12_0_instFlux'
        self.localBackgroundFluxField = 'LocalBackground_instFlux'
        self.apertureInnerInstFluxField = 'ApFlux_12_0_instFlux'
        self.apertureOuterInstFluxField = 'ApFlux_17_0_instFlux'
        self.psfCandidateName = 'Calib_psf_candidate'

        sourceSelector = self.sourceSelector["science"]

        fluxFlagName = self.instFluxField[0: -len('instFlux')] + 'flag'

        sourceSelector.flags.bad = ['PixelFlags_edge',
                                    'PixelFlags_interpolatedCenter',
                                    'PixelFlags_saturatedCenter',
                                    'PixelFlags_crCenter',
                                    'PixelFlags_bad',
                                    'PixelFlags_interpolated',
                                    'PixelFlags_saturated',
                                    'Centroid_flag',
                                    fluxFlagName]

        if self.doSubtractLocalBackground:
            localBackgroundFlagName = self.localBackgroundFluxField[0: -len('instFlux')] + 'flag'
            sourceSelector.flags.bad.append(localBackgroundFlagName)

        sourceSelector.signalToNoise.fluxField = self.instFluxField
        sourceSelector.signalToNoise.errField = self.instFluxField + 'Err'

        sourceSelector.isolated.parentName = 'parentSourceId'
        sourceSelector.isolated.nChildName = 'Deblend_nChild'

        sourceSelector.unresolved.name = 'extendedness'


class FgcmBuildStarsTableTask(FgcmBuildStarsBaseTask):
    """
    Build stars for the FGCM global calibration, using sourceTable_visit catalogs.
    """
    ConfigClass = FgcmBuildStarsTableConfig
    RunnerClass = FgcmBuildStarsRunner
    _DefaultName = "fgcmBuildStarsTable"

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "sourceTable_visit", help="Data ID, e.g. --id visit=6789")

        return parser

    def findAndGroupDataRefs(self, butler, dataRefs):
        self.log.info("Grouping dataRefs by %s" % (self.config.visitDataRefName))

        camera = butler.get('camera')

        ccdIds = []
        for detector in camera:
            ccdIds.append(detector.getId())
        # Insert our preferred referenceCCD first:
        # It is fine that this is listed twice, because we only need
        # the first calexp that is found.
        ccdIds.insert(0, self.config.referenceCCD)

        # The visitTable building code expects a dictionary of groupedDataRefs
        # keyed by visit, the first element as the "primary" calexp dataRef.
        # We then append the sourceTable_visit dataRef at the end for the
        # code which does the data reading (fgcmMakeAllStarObservations).

        groupedDataRefs = collections.defaultdict(list)
        for dataRef in dataRefs:
            visit = dataRef.dataId[self.config.visitDataRefName]

            # Find an existing calexp (we need for psf and metadata)
            # and make the relevant dataRef
            for ccdId in ccdIds:
                try:
                    calexpRef = butler.dataRef('calexp', dataId={self.config.visitDataRefName: visit,
                                                                 self.config.ccdDataRefName: ccdId})
                except RuntimeError:
                    # Not found
                    continue
                # It was found.  Add and quit out, since we only
                # need one calexp per visit.
                groupedDataRefs[visit].append(calexpRef)
                break

            # And append this dataRef
            groupedDataRefs[visit].append(dataRef)

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

        # To get the correct output schema, we use similar code as fgcmBuildStarsTask
        # We are not actually using this mapper, except to grab the outputSchema
        dataRef = groupedDataRefs[list(groupedDataRefs.keys())[0]][0]
        sourceSchema = dataRef.get('src_schema', immediate=True).schema
        sourceMapper = self._makeSourceMapper(sourceSchema)
        outputSchema = sourceMapper.getOutputSchema()

        # Construct mapping from ccd number to index
        camera = dataRef.get('camera')
        ccdMapping = {}
        for ccdIndex, detector in enumerate(camera):
            ccdMapping[detector.getId()] = ccdIndex

        approxPixelAreaFields = computeApproxPixelAreaFields(camera)

        if inStarObsCat is not None:
            fullCatalog = inStarObsCat
            comp1 = fullCatalog.schema.compare(outputSchema, outputSchema.EQUAL_KEYS)
            comp2 = fullCatalog.schema.compare(outputSchema, outputSchema.EQUAL_NAMES)
            if not comp1 or not comp2:
                raise RuntimeError("Existing fgcmStarObservations file found with mismatched schema.")
        else:
            fullCatalog = afwTable.BaseCatalog(outputSchema)

        visitKey = outputSchema['visit'].asKey()
        ccdKey = outputSchema['ccd'].asKey()
        instMagKey = outputSchema['instMag'].asKey()
        instMagErrKey = outputSchema['instMagErr'].asKey()

        # Prepare local background if desired
        if self.config.doSubtractLocalBackground:
            localBackgroundArea = np.pi*calibFluxApertureRadius**2.
        else:
            localBackground = 0.0

        # Determine which columns we need from the sourceTable_visit catalogs
        columns = self._get_sourceTable_visit_columns()

        k = 2.5/np.log(10.)

        for counter, visit in enumerate(visitCat):
            # Check if these sources have already been read and stored in the checkpoint file
            if visit['sources_read']:
                continue

            expTime = visit['exptime']

            dataRef = groupedDataRefs[visit['visit']][-1]
            srcTable = dataRef.get()

            df = srcTable.toDataFrame(columns)

            if self.config.doSubtractLocalBackground:
                localBackground = localBackgroundArea*df[self.config.localBackgroundFluxField].values
                df[self.config.instFluxField] -= localBackground

            goodSrc = self.sourceSelector.selectSources(df)
            use, = np.where(goodSrc.selected)

            tempCat = afwTable.BaseCatalog(fullCatalog.schema)
            tempCat.resize(use.size)

            tempCat['ra'][:] = np.deg2rad(df['ra'].values[use])
            tempCat['dec'][:] = np.deg2rad(df['decl'].values[use])
            tempCat['x'][:] = df['x'].values[use]
            tempCat['y'][:] = df['y'].values[use]
            tempCat[visitKey][:] = df[self.config.visitDataRefName].values[use]
            tempCat[ccdKey][:] = df[self.config.ccdDataRefName].values[use]
            tempCat['psf_candidate'] = df['Calib_psf_candidate'].values[use]

            # Need to loop over ccds here
            for detector in camera:
                ccdId = detector.getId()
                # used index for all observations with a given ccd
                use2 = (tempCat[ccdKey] == ccdId)
                tempCat['jacobian'][use2] = approxPixelAreaFields[ccdId].evaluate(tempCat['x'][use2],
                                                                                  tempCat['y'][use2])
                scaledInstFlux = (df[self.config.instFluxField].values[use[use2]] *
                                  visit['scaling'][ccdMapping[ccdId]])
                tempCat[instMagKey][use2] = (-2.5*np.log10(scaledInstFlux) + 2.5*np.log10(expTime))

            # Compute instMagErr from instFluxErr/instFlux, any scaling
            # will cancel out.
            tempCat[instMagErrKey][:] = k*(df[self.config.instFluxField + 'Err'].values[use] /
                                           df[self.config.instFluxField].values[use])

            # Apply the jacobian if configured
            if self.config.doApplyWcsJacobian:
                tempCat[instMagKey][:] -= 2.5*np.log10(tempCat['jacobian'][:])

            fullCatalog.extend(tempCat)

            # Now do the aperture information
            instMagIn = -2.5*np.log10(df[self.config.apertureInnerInstFluxField].values[use])
            instMagErrIn = k*(df[self.config.apertureInnerInstFluxField + 'Err'].values[use] /
                              df[self.config.apertureInnerInstFluxField].values[use])
            instMagOut = -2.5*np.log10(df[self.config.apertureOuterInstFluxField].values[use])
            instMagErrOut = k*(df[self.config.apertureOuterInstFluxField + 'Err'].values[use] /
                               df[self.config.apertureOuterInstFluxField].values[use])
            ok = (np.isfinite(instMagIn) & np.isfinite(instMagErrIn) &
                  np.isfinite(instMagOut) & np.isfinite(instMagErrOut))

            visit['deltaAper'] = np.median(instMagIn[ok] - instMagOut[ok])
            visit['sources_read'] = True

            self.log.info("  Found %d good stars in visit %d (deltaAper = %0.3f)",
                          use.size, visit['visit'], visit['deltaAper'])

            if ((counter % self.config.nVisitsPerCheckpoint) == 0 and
               starObsDataRef is not None and visitCatDataRef is not None):
                # We need to persist both the stars and the visit catalog which gets
                # additional metadata from each visit.
                starObsDataRef.put(fullCatalog)
                visitCatDataRef.put(visitCat)

        self.log.info("Found all good star observations in %.2f s" %
                      (time.time() - startTime))

        return fullCatalog

    def _get_sourceTable_visit_columns(self):
        """
        Get the sourceTable_visit columns from the config.

        Returns
        -------
        columns : `list`
           List of columns to read from sourceTable_visit
        """
        columns = [self.config.visitDataRefName, self.config.ccdDataRefName,
                   'ra', 'decl', 'x', 'y', self.config.psfCandidateName,
                   self.config.instFluxField, self.config.instFluxField + 'Err',
                   self.config.apertureInnerInstFluxField, self.config.apertureInnerInstFluxField + 'Err',
                   self.config.apertureOuterInstFluxField, self.config.apertureOuterInstFluxField + 'Err']
        if self.sourceSelector.config.doFlags:
            columns.extend(self.sourceSelector.config.flags.bad)
        if self.sourceSelector.config.doUnresolved:
            columns.append(self.sourceSelector.config.unresolved.name)
        if self.sourceSelector.config.doIsolated:
            columns.append(self.sourceSelector.config.isolated.parentName)
            columns.append(self.sourceSelector.config.isolated.nChildName)
        if self.config.doSubtractLocalBackground:
            columns.append(self.config.localBackgroundFluxField)

        return columns
