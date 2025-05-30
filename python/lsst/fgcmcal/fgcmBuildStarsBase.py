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
"""Base class for BuildStars using src tables or sourceTable_visit tables.
"""

import abc

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
from lsst.daf.base import PropertyList
from lsst.daf.base.dateTime import DateTime
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .fgcmLoadReferenceCatalog import FgcmLoadReferenceCatalogTask
from .utilities import computeReferencePixelScale, countDetectors

import fgcm

REFSTARS_FORMAT_VERSION = 1

__all__ = ['FgcmBuildStarsConfigBase', 'FgcmBuildStarsBaseTask']


class FgcmBuildStarsConfigBase(pexConfig.Config):
    """Base config for FgcmBuildStars tasks"""

    instFluxField = pexConfig.Field(
        doc=("Faull name of the source instFlux field to use, including 'instFlux'. "
             "The associated flag will be implicitly included in badFlags"),
        dtype=str,
        default='slot_CalibFlux_instFlux',
    )
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
    randomSeed = pexConfig.Field(
        doc="Random seed for high density down-sampling.",
        dtype=int,
        default=None,
        optional=True,
    )
    matchNside = pexConfig.Field(
        doc="Healpix Nside for matching",
        dtype=int,
        default=4096,
    )
    coarseNside = pexConfig.Field(
        doc="Healpix coarse Nside for partitioning matches",
        dtype=int,
        default=8,
    )
    physicalFilterMap = pexConfig.DictField(
        doc="Mapping from 'physicalFilter' to band.",
        keytype=str,
        itemtype=str,
        default={},
    )
    requiredBands = pexConfig.ListField(
        doc="Bands required for each star",
        dtype=str,
        default=(),
    )
    primaryBands = pexConfig.ListField(
        doc=("Bands for 'primary' star matches. "
             "A star must be observed in one of these bands to be considered "
             "as a calibration star."),
        dtype=str,
        default=None
    )
    doApplyWcsJacobian = pexConfig.Field(
        doc="Apply the jacobian of the WCS to the star observations prior to fit?",
        dtype=bool,
        default=True
    )
    doModelErrorsWithBackground = pexConfig.Field(
        doc="Model flux errors with background term?",
        dtype=bool,
        default=True
    )
    psfCandidateName = pexConfig.Field(
        doc="Name of field with psf candidate flag for propagation",
        dtype=str,
        default="calib_psf_candidate"
    )
    doSubtractLocalBackground = pexConfig.Field(
        doc=("Subtract the local background before performing calibration? "
             "This is only supported for circular aperture calibration fluxes."),
        dtype=bool,
        default=False
    )
    localBackgroundFluxField = pexConfig.Field(
        doc="Full name of the local background instFlux field to use.",
        dtype=str,
        default='base_LocalBackground_instFlux'
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources",
        default="science"
    )
    apertureInnerInstFluxField = pexConfig.Field(
        doc=("Full name of instFlux field that contains inner aperture "
             "flux for aperture correction proxy"),
        dtype=str,
        default='base_CircularApertureFlux_12_0_instFlux'
    )
    apertureOuterInstFluxField = pexConfig.Field(
        doc=("Full name of instFlux field that contains outer aperture "
             "flux for aperture correction proxy"),
        dtype=str,
        default='base_CircularApertureFlux_17_0_instFlux'
    )
    doReferenceMatches = pexConfig.Field(
        doc="Match reference catalog as additional constraint on calibration",
        dtype=bool,
        default=True,
    )
    fgcmLoadReferenceCatalog = pexConfig.ConfigurableField(
        target=FgcmLoadReferenceCatalogTask,
        doc="FGCM reference object loader",
    )
    nVisitsPerCheckpoint = pexConfig.Field(
        doc="Number of visits read between checkpoints",
        dtype=int,
        default=500,
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["science"]
        sourceSelector.setDefaults()

        sourceSelector.doFlags = True
        sourceSelector.doUnresolved = True
        sourceSelector.doSignalToNoise = True
        sourceSelector.doIsolated = True
        sourceSelector.doRequireFiniteRaDec = True

        sourceSelector.signalToNoise.minimum = 10.0
        sourceSelector.signalToNoise.maximum = 1000.0

        # FGCM operates on unresolved sources, and this setting is
        # appropriate for the current base_ClassificationExtendedness
        sourceSelector.unresolved.maximum = 0.5


class FgcmBuildStarsBaseTask(pipeBase.PipelineTask, abc.ABC):
    """
    Base task to build stars for FGCM global calibration
    """
    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)

        self.makeSubtask("sourceSelector")
        # Only log warning and fatal errors from the sourceSelector
        self.sourceSelector.log.setLevel(self.sourceSelector.log.WARN)

    def fgcmMakeAllStarObservations(self, groupedHandles, visitCat,
                                    sourceSchema,
                                    camera,
                                    calibFluxApertureRadius=None):
        """
        Compile all good star observations from visits in visitCat.

        Parameters
        ----------
        groupedHandles : `dict` [`list` [`lsst.daf.butler.DeferredDatasetHandle`]]
            Dataset handles, grouped by visit.
        visitCat : `afw.table.BaseCatalog`
            Catalog with visit data for FGCM
        sourceSchema : `lsst.afw.table.Schema`
            Schema for the input src catalogs.
        camera : `lsst.afw.cameraGeom.Camera`
        calibFluxApertureRadius : `float`, optional
            Aperture radius for calibration flux.
        inStarObsCat : `afw.table.BaseCatalog`
            Input observation catalog.  If this is incomplete, observations
            will be appended from when it was cut off.

        Returns
        -------
        fgcmStarObservations : `afw.table.BaseCatalog`
            Full catalog of good observations.

        Raises
        ------
        RuntimeError: Raised if doSubtractLocalBackground is True and
           calibFluxApertureRadius is not set.
        """
        raise NotImplementedError("fgcmMakeAllStarObservations not implemented.")

    def fgcmMakeVisitCatalog(self, camera, groupedHandles, useScienceDetectors=False):
        """
        Make a visit catalog with all the keys from each visit

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`
            Camera from the butler
        groupedHandles : `dict` [`list` [`lsst.daf.butler.DeferredDatasetHandle`]]
            Dataset handles, grouped by visit.
        useScienceDetectors : `bool`, optional
            Limit to science detectors?

        Returns
        -------
        visitCat: `afw.table.BaseCatalog`
        """

        self.log.info("Assembling visitCatalog from %d visits", len(groupedHandles))

        nCcd = countDetectors(camera, useScienceDetectors)

        schema = self._makeFgcmVisitSchema(nCcd)

        visitCat = afwTable.BaseCatalog(schema)
        visitCat.reserve(len(groupedHandles))
        visitCat.resize(len(groupedHandles))

        visitCat['visit'] = list(groupedHandles.keys())
        visitCat['used'] = 0
        visitCat['sources_read'] = False

        defaultPixelScale = computeReferencePixelScale(camera, useScienceDetectors=useScienceDetectors)

        # No matter what, fill the catalog. This will check if it was
        # already read.
        self._fillVisitCatalog(visitCat, groupedHandles, defaultPixelScale)

        return visitCat

    def _fillVisitCatalog(self, visitCat, groupedHandles, defaultPixelScale):
        """
        Fill the visit catalog with visit metadata

        Parameters
        ----------
        visitCat : `afw.table.BaseCatalog`
            Visit catalog.  See _makeFgcmVisitSchema() for schema definition.
        groupedHandles : `dict` [`list` [`lsst.daf.butler.DeferredDatasetHandle`]]
            Dataset handles, grouped by visit.
        defaultPixelScale : `float`
            Default pixel scale to use if not in visit summary (arcsecond/pixel).
        """

        # Guarantee that these are sorted.
        for i, visit in enumerate(sorted(groupedHandles)):
            if (i % self.config.nVisitsPerCheckpoint) == 0:
                self.log.info("Retrieving metadata for visit %d (%d/%d)", visit, i, len(groupedHandles))

            handle = groupedHandles[visit][0]
            summary = handle.get()

            summaryRow = summary.find(self.config.referenceCCD)
            if summaryRow is None:
                # Take the first available ccd if reference isn't available
                summaryRow = summary[0]

            visitInfo = summaryRow.getVisitInfo()
            physicalFilter = summaryRow['physical_filter']
            # Compute the median psf sigma and fwhm if possible.
            if 'pixelScale' in summary.schema:
                # This is not available in the older test summaries
                pixelScales = summary['pixelScale']
            else:
                pixelScales = np.full(len(summary['psfSigma']), defaultPixelScale)
            psfSigmas = summary['psfSigma']
            goodSigma, = np.where((np.nan_to_num(psfSigmas) > 0) & (np.nan_to_num(pixelScales) > 0))
            if goodSigma.size > 2:
                psfSigma = np.median(psfSigmas[goodSigma])
                psfFwhm = np.median(psfSigmas[goodSigma] * pixelScales[goodSigma]) * np.sqrt(8.*np.log(2.))
            elif goodSigma.size > 0:
                psfSigma = psfSigmas[goodSigma[0]]
                psfFwhm = psfSigmas[goodSigma[0]] * pixelScales[goodSigma[0]] * np.sqrt(8.)*np.log(2.)
            else:
                self.log.warning("Could not find any good summary psfSigma for visit %d", visit)
                psfSigma = 0.0
                psfFwhm = 0.0
            # Compute median background if possible
            goodBackground, = np.where(np.nan_to_num(summary['skyBg']) > 0.0)
            if goodBackground.size > 2:
                skyBackground = np.median(summary['skyBg'][goodBackground])
            elif goodBackground.size > 0:
                skyBackground = summary['skyBg'][goodBackground[0]]
            else:
                self.log.warning('Could not find any good summary skyBg for visit %d', visit)
                skyBackground = -1.0

            rec = visitCat[i]
            rec['visit'] = visit
            rec['physicalFilter'] = physicalFilter
            # TODO DM-26991: Use the wcs to refine the focal-plane center.
            radec = visitInfo.getBoresightRaDec()
            rec['telra'] = radec.getRa().asDegrees()
            rec['teldec'] = radec.getDec().asDegrees()
            rec['telha'] = visitInfo.getBoresightHourAngle().asDegrees()
            rec['telrot'] = visitInfo.getBoresightRotAngle().asDegrees()
            rec['mjd'] = visitInfo.getDate().get(system=DateTime.MJD)
            rec['exptime'] = visitInfo.getExposureTime()
            # convert from Pa to millibar
            # Note that I don't know if this unit will need to be per-camera config
            rec['pmb'] = visitInfo.getWeather().getAirPressure() / 100
            # Flag to signify if this is a "deep" field.  Not currently used
            rec['deepFlag'] = 0
            # Relative flat scaling (1.0 means no relative scaling)
            rec['scaling'][:] = 1.0
            # Median delta aperture, to be measured from stars
            rec['deltaAper'] = 0.0
            rec['psfSigma'] = psfSigma
            rec['psfFwhm'] = psfFwhm
            rec['skyBackground'] = skyBackground
            rec['used'] = 1

    def _makeSourceMapper(self, sourceSchema):
        """
        Make a schema mapper for fgcm sources

        Parameters
        ----------
        sourceSchema: `afwTable.Schema`
           Default source schema from the butler

        Returns
        -------
        sourceMapper: `afwTable.schemaMapper`
           Mapper to the FGCM source schema
        """

        # create a mapper to the preferred output
        sourceMapper = afwTable.SchemaMapper(sourceSchema)

        # map to ra/dec
        sourceMapper.addMapping(sourceSchema['coord_ra'].asKey(), 'ra')
        sourceMapper.addMapping(sourceSchema['coord_dec'].asKey(), 'dec')
        sourceMapper.addMapping(sourceSchema['slot_Centroid_x'].asKey(), 'x')
        sourceMapper.addMapping(sourceSchema['slot_Centroid_y'].asKey(), 'y')
        # Add the mapping if the field exists in the input catalog.
        # If the field does not exist, simply add it (set to False).
        # This field is not required for calibration, but is useful
        # to collate if available.
        try:
            sourceMapper.addMapping(sourceSchema[self.config.psfCandidateName].asKey(),
                                    'psf_candidate')
        except LookupError:
            sourceMapper.editOutputSchema().addField(
                "psf_candidate", type='Flag',
                doc=("Flag set if the source was a candidate for PSF determination, "
                     "as determined by the star selector."))

        # and add the fields we want
        sourceMapper.editOutputSchema().addField(
            "visit", type=np.int64, doc="Visit number")
        sourceMapper.editOutputSchema().addField(
            "ccd", type=np.int32, doc="CCD number")
        sourceMapper.editOutputSchema().addField(
            "instMag", type=np.float32, doc="Instrumental magnitude")
        sourceMapper.editOutputSchema().addField(
            "instMagErr", type=np.float32, doc="Instrumental magnitude error")
        sourceMapper.editOutputSchema().addField(
            "jacobian", type=np.float32, doc="Relative pixel scale from wcs jacobian")
        sourceMapper.editOutputSchema().addField(
            "deltaMagBkg", type=np.float32, doc="Change in magnitude due to local background offset")
        sourceMapper.editOutputSchema().addField(
            "deltaMagAper", type=np.float32, doc="Change in magnitude from larger to smaller aperture")

        return sourceMapper

    def fgcmMatchStars(self, visitCat, obsCat, lutHandle=None):
        """
        Use FGCM code to match observations into unique stars.

        Parameters
        ----------
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for fgcm
        obsCat: `afw.table.BaseCatalog`
           Full catalog of star observations for fgcm
        lutHandle: `lsst.daf.butler.DeferredDatasetHandle`, optional
           Data reference to fgcm look-up table (used if matching reference stars).

        Returns
        -------
        fgcmStarIdCat: `afw.table.BaseCatalog`
           Catalog of unique star identifiers and index keys
        fgcmStarIndicesCat: `afwTable.BaseCatalog`
           Catalog of unique star indices
        fgcmRefCat: `afw.table.BaseCatalog`
           Catalog of matched reference stars.
           Will be None if `config.doReferenceMatches` is False.
        """
        # get filter names into a numpy array...
        # This is the type that is expected by the fgcm code
        visitFilterNames = np.zeros(len(visitCat), dtype='S30')
        for i in range(len(visitCat)):
            visitFilterNames[i] = visitCat[i]['physicalFilter']

        # match to put filterNames with observations
        visitIndex = np.searchsorted(visitCat['visit'],
                                     obsCat['visit'])

        obsFilterNames = visitFilterNames[visitIndex]

        if self.config.doReferenceMatches:
            # Get the reference filter names, using the LUT
            lutCat = lutHandle.get()

            stdFilterDict = {filterName: stdFilter for (filterName, stdFilter) in
                             zip(lutCat[0]['physicalFilters'].split(','),
                                 lutCat[0]['stdPhysicalFilters'].split(','))}
            stdLambdaDict = {stdFilter: stdLambda for (stdFilter, stdLambda) in
                             zip(lutCat[0]['stdPhysicalFilters'].split(','),
                                 lutCat[0]['lambdaStdFilter'])}

            del lutCat

            referenceFilterNames = self._getReferenceFilterNames(visitCat,
                                                                 stdFilterDict,
                                                                 stdLambdaDict)
            self.log.info("Using the following reference filters: %s" %
                          (', '.join(referenceFilterNames)))

        else:
            # This should be an empty list
            referenceFilterNames = []

        # make the fgcm starConfig dict
        starConfig = {'logger': self.log,
                      'useHtm': True,
                      'filterToBand': self.config.physicalFilterMap,
                      'requiredBands': self.config.requiredBands,
                      'minPerBand': self.config.minPerBand,
                      'matchRadius': self.config.matchRadius,
                      'isolationRadius': self.config.isolationRadius,
                      'matchNSide': self.config.matchNside,
                      'coarseNSide': self.config.coarseNside,
                      'densNSide': self.config.densityCutNside,
                      'densMaxPerPixel': self.config.densityCutMaxPerPixel,
                      'randomSeed': self.config.randomSeed,
                      'primaryBands': self.config.primaryBands,
                      'referenceFilterNames': referenceFilterNames}

        # initialize the FgcmMakeStars object
        fgcmMakeStars = fgcm.FgcmMakeStars(starConfig)

        # make the primary stars
        # note that the ra/dec native Angle format is radians
        # We determine the conversion from the native units (typically
        # radians) to degrees for the first observation.  This allows us
        # to treate ra/dec as numpy arrays rather than Angles, which would
        # be approximately 600x slower.
        conv = obsCat[0]['ra'].asDegrees() / float(obsCat[0]['ra'])
        fgcmMakeStars.makePrimaryStars(obsCat['ra'] * conv,
                                       obsCat['dec'] * conv,
                                       filterNameArray=obsFilterNames,
                                       bandSelected=False)

        # and match all the stars
        fgcmMakeStars.makeMatchedStars(obsCat['ra'] * conv,
                                       obsCat['dec'] * conv,
                                       obsFilterNames)

        if self.config.doReferenceMatches:
            fgcmMakeStars.makeReferenceMatches(self.fgcmLoadReferenceCatalog)

        # now persist

        objSchema = self._makeFgcmObjSchema()

        # make catalog and records
        fgcmStarIdCat = afwTable.BaseCatalog(objSchema)
        fgcmStarIdCat.reserve(fgcmMakeStars.objIndexCat.size)
        for i in range(fgcmMakeStars.objIndexCat.size):
            fgcmStarIdCat.addNew()

        # fill the catalog
        fgcmStarIdCat['fgcm_id'][:] = fgcmMakeStars.objIndexCat['fgcm_id']
        fgcmStarIdCat['ra'][:] = fgcmMakeStars.objIndexCat['ra']
        fgcmStarIdCat['dec'][:] = fgcmMakeStars.objIndexCat['dec']
        fgcmStarIdCat['obsArrIndex'][:] = fgcmMakeStars.objIndexCat['obsarrindex']
        fgcmStarIdCat['nObs'][:] = fgcmMakeStars.objIndexCat['nobs']

        obsSchema = self._makeFgcmObsSchema()

        fgcmStarIndicesCat = afwTable.BaseCatalog(obsSchema)
        fgcmStarIndicesCat.reserve(fgcmMakeStars.obsIndexCat.size)
        for i in range(fgcmMakeStars.obsIndexCat.size):
            fgcmStarIndicesCat.addNew()

        fgcmStarIndicesCat['obsIndex'][:] = fgcmMakeStars.obsIndexCat['obsindex']

        if self.config.doReferenceMatches:
            refSchema = self._makeFgcmRefSchema(len(referenceFilterNames))

            fgcmRefCat = afwTable.BaseCatalog(refSchema)
            fgcmRefCat.reserve(fgcmMakeStars.referenceCat.size)

            for i in range(fgcmMakeStars.referenceCat.size):
                fgcmRefCat.addNew()

            fgcmRefCat['fgcm_id'][:] = fgcmMakeStars.referenceCat['fgcm_id']
            fgcmRefCat['refMag'][:, :] = fgcmMakeStars.referenceCat['refMag']
            fgcmRefCat['refMagErr'][:, :] = fgcmMakeStars.referenceCat['refMagErr']

            md = PropertyList()
            md.set("REFSTARS_FORMAT_VERSION", REFSTARS_FORMAT_VERSION)
            md.set("FILTERNAMES", referenceFilterNames)
            fgcmRefCat.setMetadata(md)

        else:
            fgcmRefCat = None

        return fgcmStarIdCat, fgcmStarIndicesCat, fgcmRefCat

    def _makeFgcmVisitSchema(self, nCcd):
        """
        Make a schema for an fgcmVisitCatalog

        Parameters
        ----------
        nCcd: `int`
           Number of CCDs in the camera

        Returns
        -------
        schema: `afwTable.Schema`
        """

        schema = afwTable.Schema()
        schema.addField('visit', type=np.int64, doc="Visit number")
        schema.addField('physicalFilter', type=str, size=30, doc="Physical filter")
        schema.addField('telra', type=np.float64, doc="Pointing RA (deg)")
        schema.addField('teldec', type=np.float64, doc="Pointing Dec (deg)")
        schema.addField('telha', type=np.float64, doc="Pointing Hour Angle (deg)")
        schema.addField('telrot', type=np.float64, doc="Camera rotation (deg)")
        schema.addField('mjd', type=np.float64, doc="MJD of visit")
        schema.addField('exptime', type=np.float32, doc="Exposure time")
        schema.addField('pmb', type=np.float32, doc="Pressure (millibar)")
        schema.addField('psfSigma', type=np.float32, doc="PSF sigma (median); pixels")
        schema.addField('psfFwhm', type=np.float32, doc="PSF FWHM (median); arcseconds")
        schema.addField('deltaAper', type=np.float32, doc="Delta-aperture")
        schema.addField('skyBackground', type=np.float32, doc="Sky background (ADU) (reference CCD)")
        # the following field is not used yet
        schema.addField('deepFlag', type=np.int32, doc="Deep observation")
        schema.addField('scaling', type='ArrayD', doc="Scaling applied due to flat adjustment",
                        size=nCcd)
        schema.addField('used', type=np.int32, doc="This visit has been ingested.")
        schema.addField('sources_read', type='Flag', doc="This visit had sources read.")

        return schema

    def _makeFgcmObjSchema(self):
        """
        Make a schema for the objIndexCat from fgcmMakeStars

        Returns
        -------
        schema: `afwTable.Schema`
        """

        objSchema = afwTable.Schema()
        objSchema.addField('fgcm_id', type=np.int32, doc='FGCM Unique ID')
        # Will investigate making these angles...
        objSchema.addField('ra', type=np.float64, doc='Mean object RA (deg)')
        objSchema.addField('dec', type=np.float64, doc='Mean object Dec (deg)')
        objSchema.addField('obsArrIndex', type=np.int32,
                           doc='Index in obsIndexTable for first observation')
        objSchema.addField('nObs', type=np.int32, doc='Total number of observations')

        return objSchema

    def _makeFgcmObsSchema(self):
        """
        Make a schema for the obsIndexCat from fgcmMakeStars

        Returns
        -------
        schema: `afwTable.Schema`
        """

        obsSchema = afwTable.Schema()
        obsSchema.addField('obsIndex', type=np.int32, doc='Index in observation table')

        return obsSchema

    def _makeFgcmRefSchema(self, nReferenceBands):
        """
        Make a schema for the referenceCat from fgcmMakeStars

        Parameters
        ----------
        nReferenceBands: `int`
           Number of reference bands

        Returns
        -------
        schema: `afwTable.Schema`
        """

        refSchema = afwTable.Schema()
        refSchema.addField('fgcm_id', type=np.int32, doc='FGCM Unique ID')
        refSchema.addField('refMag', type='ArrayF', doc='Reference magnitude array (AB)',
                           size=nReferenceBands)
        refSchema.addField('refMagErr', type='ArrayF', doc='Reference magnitude error array',
                           size=nReferenceBands)

        return refSchema

    def _getReferenceFilterNames(self, visitCat, stdFilterDict, stdLambdaDict):
        """
        Get the reference filter names, in wavelength order, from the visitCat and
        information from the look-up-table.

        Parameters
        ----------
        visitCat: `afw.table.BaseCatalog`
           Catalog with visit data for FGCM
        stdFilterDict: `dict`
           Mapping of filterName to stdFilterName from LUT
        stdLambdaDict: `dict`
           Mapping of stdFilterName to stdLambda from LUT

        Returns
        -------
        referenceFilterNames: `list`
           Wavelength-ordered list of reference filter names
        """

        # Find the unique list of filter names in visitCat
        filterNames = np.unique(visitCat.asAstropy()['physicalFilter'])

        # Find the unique list of "standard" filters
        stdFilterNames = {stdFilterDict[filterName] for filterName in filterNames}

        # And sort these by wavelength
        referenceFilterNames = sorted(stdFilterNames, key=stdLambdaDict.get)

        return referenceFilterNames
