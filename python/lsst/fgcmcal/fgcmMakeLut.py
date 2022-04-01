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
"""Make a look-up-table (LUT) for FGCM calibration.

This task computes a look-up-table for the range in expected atmosphere
variation and variation in instrumental throughput (as tracked by the
transmission_filter products).  By pre-computing linearized integrals,
the FGCM fit is orders of magnitude faster for stars with a broad range
of colors and observing bands, yielding precision at the 1-2 mmag level.

Computing a LUT requires running MODTRAN or with a pre-generated
atmosphere table packaged with fgcm.
"""

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes
import lsst.afw.table as afwTable
import lsst.afw.cameraGeom as afwCameraGeom
from .utilities import lookupStaticCalibrations

import fgcm

__all__ = ['FgcmMakeLutParametersConfig', 'FgcmMakeLutConfig', 'FgcmMakeLutTask']


class FgcmMakeLutConnections(pipeBase.PipelineTaskConnections,
                             dimensions=('instrument',),
                             defaultTemplates={}):
    camera = connectionTypes.PrerequisiteInput(
        doc="Camera instrument",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibrations,
        isCalibration=True,
    )

    transmission_optics = connectionTypes.PrerequisiteInput(
        doc="Optics transmission curve information",
        name="transmission_optics",
        storageClass="TransmissionCurve",
        dimensions=("instrument",),
        lookupFunction=lookupStaticCalibrations,
        isCalibration=True,
        deferLoad=True,
    )

    transmission_sensor = connectionTypes.PrerequisiteInput(
        doc="Sensor transmission curve information",
        name="transmission_sensor",
        storageClass="TransmissionCurve",
        dimensions=("instrument", "detector",),
        lookupFunction=lookupStaticCalibrations,
        isCalibration=True,
        deferLoad=True,
        multiple=True,
    )

    transmission_filter = connectionTypes.PrerequisiteInput(
        doc="Filter transmission curve information",
        name="transmission_filter",
        storageClass="TransmissionCurve",
        dimensions=("band", "instrument", "physical_filter",),
        lookupFunction=lookupStaticCalibrations,
        isCalibration=True,
        deferLoad=True,
        multiple=True,
    )

    fgcmLookUpTable = connectionTypes.Output(
        doc=("Atmosphere + instrument look-up-table for FGCM throughput and "
             "chromatic corrections."),
        name="fgcmLookUpTable",
        storageClass="Catalog",
        dimensions=("instrument",),
    )


class FgcmMakeLutParametersConfig(pexConfig.Config):
    """Config for parameters if atmosphereTableName not available"""
    # TODO: When DM-16511 is done, it will be possible to get the
    # telescope elevation directly from the camera.
    elevation = pexConfig.Field(
        doc="Telescope elevation (m)",
        dtype=float,
        default=None,
    )
    pmbRange = pexConfig.ListField(
        doc=("Barometric Pressure range (millibar) "
             "Recommended range depends on the site."),
        dtype=float,
        default=None,
    )
    pmbSteps = pexConfig.Field(
        doc="Barometric Pressure number of steps",
        dtype=int,
        default=5,
    )
    pwvRange = pexConfig.ListField(
        doc=("Precipitable Water Vapor range (mm) "
             "Recommended range depends on the site."),
        dtype=float,
        default=None,
    )
    pwvSteps = pexConfig.Field(
        doc="Precipitable Water Vapor number of steps",
        dtype=int,
        default=15,
    )
    o3Range = pexConfig.ListField(
        doc="Ozone range (dob)",
        dtype=float,
        default=[220.0, 310.0],
    )
    o3Steps = pexConfig.Field(
        doc="Ozone number of steps",
        dtype=int,
        default=3,
    )
    tauRange = pexConfig.ListField(
        doc="Aerosol Optical Depth range (unitless)",
        dtype=float,
        default=[0.002, 0.35],
    )
    tauSteps = pexConfig.Field(
        doc="Aerosol Optical Depth number of steps",
        dtype=int,
        default=11,
    )
    alphaRange = pexConfig.ListField(
        doc="Aerosol alpha range (unitless)",
        dtype=float,
        default=[0.0, 2.0],
    )
    alphaSteps = pexConfig.Field(
        doc="Aerosol alpha number of steps",
        dtype=int,
        default=9,
    )
    zenithRange = pexConfig.ListField(
        doc="Zenith angle range (degree)",
        dtype=float,
        default=[0.0, 70.0],
    )
    zenithSteps = pexConfig.Field(
        doc="Zenith angle number of steps",
        dtype=int,
        default=21,
    )
    # Note that the standard atmosphere parameters depend on the observatory
    # and elevation, and so these should be set on a per-camera basis.
    pmbStd = pexConfig.Field(
        doc=("Standard Atmosphere pressure (millibar); "
             "Recommended default depends on the site."),
        dtype=float,
        default=None,
    )
    pwvStd = pexConfig.Field(
        doc=("Standard Atmosphere PWV (mm); "
             "Recommended default depends on the site."),
        dtype=float,
        default=None,
    )
    o3Std = pexConfig.Field(
        doc="Standard Atmosphere O3 (dob)",
        dtype=float,
        default=263.0,
    )
    tauStd = pexConfig.Field(
        doc="Standard Atmosphere aerosol optical depth",
        dtype=float,
        default=0.03,
    )
    alphaStd = pexConfig.Field(
        doc="Standard Atmosphere aerosol alpha",
        dtype=float,
        default=1.0,
    )
    airmassStd = pexConfig.Field(
        doc=("Standard Atmosphere airmass; "
             "Recommended default depends on the survey strategy."),
        dtype=float,
        default=None,
    )
    lambdaNorm = pexConfig.Field(
        doc="Aerosol Optical Depth normalization wavelength (Angstrom)",
        dtype=float,
        default=7750.0,
    )
    lambdaStep = pexConfig.Field(
        doc="Wavelength step for generating atmospheres (nm)",
        dtype=float,
        default=0.5,
    )
    lambdaRange = pexConfig.ListField(
        doc="Wavelength range for LUT (Angstrom)",
        dtype=float,
        default=[3000.0, 11000.0],
    )


class FgcmMakeLutConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=FgcmMakeLutConnections):
    """Config for FgcmMakeLutTask"""
    physicalFilters = pexConfig.ListField(
        doc="List of physicalFilter labels to generate look-up table.",
        dtype=str,
        default=[],
    )
    stdPhysicalFilterOverrideMap = pexConfig.DictField(
        doc=("Override mapping from physical filter labels to 'standard' physical "
             "filter labels. The 'standard' physical filter defines the transmission "
             "curve that the FGCM standard bandpass will be based on. "
             "Any filter not listed here will be mapped to "
             "itself (e.g. g->g or HSC-G->HSC-G).  Use this override for cross-"
             "filter calibration such as HSC-R->HSC-R2 and HSC-I->HSC-I2."),
        keytype=str,
        itemtype=str,
        default={},
    )
    atmosphereTableName = pexConfig.Field(
        doc="FGCM name or filename of precomputed atmospheres",
        dtype=str,
        default=None,
        optional=True,
    )
    parameters = pexConfig.ConfigField(
        doc="Atmosphere parameters (required if no atmosphereTableName)",
        dtype=FgcmMakeLutParametersConfig,
        default=None,
        check=None)

    def validate(self):
        """
        Validate the config parameters.

        This method behaves differently from the parent validate in the case
        that atmosphereTableName is set.  In this case, the config values
        for standard values, step sizes, and ranges are loaded
        directly from the specified atmosphereTableName.
        """
        # check that filterNames and stdFilterNames are okay
        self._fields['physicalFilters'].validate(self)
        self._fields['stdPhysicalFilterOverrideMap'].validate(self)

        if self.atmosphereTableName is None:
            # Validate the parameters
            self._fields['parameters'].validate(self)


class FgcmMakeLutTask(pipeBase.PipelineTask):
    """
    Make Look-Up Table for FGCM.

    This task computes a look-up-table for the range in expected atmosphere
    variation and variation in instrumental throughput (as tracked by the
    transmission_filter products).  By pre-computing linearized integrals,
    the FGCM fit is orders of magnitude faster for stars with a broad range
    of colors and observing bands, yielding precision at the 1-2 mmag level.

    Computing a LUT requires running MODTRAN or with a pre-generated
    atmosphere table packaged with fgcm.
    """

    ConfigClass = FgcmMakeLutConfig
    _DefaultName = "fgcmMakeLut"

    def __init__(self, butler=None, initInputs=None, **kwargs):
        super().__init__(**kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        camera = butlerQC.get(inputRefs.camera)

        opticsHandle = butlerQC.get(inputRefs.transmission_optics)

        sensorHandles = butlerQC.get(inputRefs.transmission_sensor)
        sensorHandleDict = {sensorHandle.dataId.byName()['detector']: sensorHandle for
                            sensorHandle in sensorHandles}

        filterHandles = butlerQC.get(inputRefs.transmission_filter)
        filterHandleDict = {filterHandle.dataId['physical_filter']: filterHandle for
                            filterHandle in filterHandles}

        lutCat = self._fgcmMakeLut(camera,
                                   opticsHandle,
                                   sensorHandleDict,
                                   filterHandleDict)
        butlerQC.put(lutCat, outputRefs.fgcmLookUpTable)

    def _fgcmMakeLut(self, camera, opticsHandle, sensorHandleDict,
                     filterHandleDict):
        """
        Make a FGCM Look-up Table

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`
            Camera from the butler.
        opticsHandle : `lsst.daf.butler.DeferredDatasetHandle`
            Reference to optics transmission curve.
        sensorHandleDict : `dict` of [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Dictionary of references to sensor transmission curves.  Key will
            be detector id.
        filterHandleDict : `dict` of [`str`, `lsst.daf.butler.DeferredDatasetHandle`]
            Dictionary of references to filter transmission curves.  Key will
            be physical filter label.

        Returns
        -------
        fgcmLookUpTable : `BaseCatalog`
            The FGCM look-up table.
        """
        # number of ccds from the length of the camera iterator
        nCcd = len(camera)
        self.log.info("Found %d ccds for look-up table" % (nCcd))

        # Load in optics, etc.
        self._loadThroughputs(camera,
                              opticsHandle,
                              sensorHandleDict,
                              filterHandleDict)

        lutConfig = self._createLutConfig(nCcd)

        # make the lut object
        self.log.info("Making the LUT maker object")
        self.fgcmLutMaker = fgcm.FgcmLUTMaker(lutConfig)

        # generate the throughput dictionary.

        # these will be in Angstroms
        # note that lambdaStep is currently in nm, because of historical
        # reasons in the code.  Convert to Angstroms here.
        throughputLambda = np.arange(self.fgcmLutMaker.lambdaRange[0],
                                     self.fgcmLutMaker.lambdaRange[1]+self.fgcmLutMaker.lambdaStep*10,
                                     self.fgcmLutMaker.lambdaStep*10.)

        self.log.info("Built throughput lambda, %.1f-%.1f, step %.2f" %
                      (throughputLambda[0], throughputLambda[-1],
                       throughputLambda[1] - throughputLambda[0]))

        throughputDict = {}
        for i, physicalFilter in enumerate(self.config.physicalFilters):
            tDict = {}
            tDict['LAMBDA'] = throughputLambda
            for ccdIndex, detector in enumerate(camera):
                tDict[ccdIndex] = self._getThroughputDetector(detector, physicalFilter, throughputLambda)
            throughputDict[physicalFilter] = tDict

        # set the throughputs
        self.fgcmLutMaker.setThroughputs(throughputDict)

        # make the LUT
        self.log.info("Making LUT")
        self.fgcmLutMaker.makeLUT()

        # and save the LUT

        # build the index values
        comma = ','
        physicalFilterString = comma.join(self.config.physicalFilters)
        stdPhysicalFilterString = comma.join(self._getStdPhysicalFilterList())

        atmosphereTableName = 'NoTableWasUsed'
        if self.config.atmosphereTableName is not None:
            atmosphereTableName = self.config.atmosphereTableName

        lutSchema = self._makeLutSchema(physicalFilterString, stdPhysicalFilterString,
                                        atmosphereTableName)

        lutCat = self._makeLutCat(lutSchema, physicalFilterString,
                                  stdPhysicalFilterString, atmosphereTableName)
        return lutCat

    def _getStdPhysicalFilterList(self):
        """Get the standard physical filter lists from config.physicalFilters
        and config.stdPhysicalFilterOverrideMap

        Returns
        -------
        stdPhysicalFilters : `list`
        """
        override = self.config.stdPhysicalFilterOverrideMap
        return [override.get(physicalFilter, physicalFilter) for
                physicalFilter in self.config.physicalFilters]

    def _createLutConfig(self, nCcd):
        """
        Create the fgcmLut config dictionary

        Parameters
        ----------
        nCcd: `int`
           Number of CCDs in the camera
        """

        # create the common stub of the lutConfig
        lutConfig = {}
        lutConfig['logger'] = self.log
        lutConfig['filterNames'] = self.config.physicalFilters
        lutConfig['stdFilterNames'] = self._getStdPhysicalFilterList()
        lutConfig['nCCD'] = nCcd

        # atmosphereTable already validated if available
        if self.config.atmosphereTableName is not None:
            lutConfig['atmosphereTableName'] = self.config.atmosphereTableName
        else:
            # use the regular paramters (also validated if needed)
            lutConfig['elevation'] = self.config.parameters.elevation
            lutConfig['pmbRange'] = self.config.parameters.pmbRange
            lutConfig['pmbSteps'] = self.config.parameters.pmbSteps
            lutConfig['pwvRange'] = self.config.parameters.pwvRange
            lutConfig['pwvSteps'] = self.config.parameters.pwvSteps
            lutConfig['o3Range'] = self.config.parameters.o3Range
            lutConfig['o3Steps'] = self.config.parameters.o3Steps
            lutConfig['tauRange'] = self.config.parameters.tauRange
            lutConfig['tauSteps'] = self.config.parameters.tauSteps
            lutConfig['alphaRange'] = self.config.parameters.alphaRange
            lutConfig['alphaSteps'] = self.config.parameters.alphaSteps
            lutConfig['zenithRange'] = self.config.parameters.zenithRange
            lutConfig['zenithSteps'] = self.config.parameters.zenithSteps
            lutConfig['pmbStd'] = self.config.parameters.pmbStd
            lutConfig['pwvStd'] = self.config.parameters.pwvStd
            lutConfig['o3Std'] = self.config.parameters.o3Std
            lutConfig['tauStd'] = self.config.parameters.tauStd
            lutConfig['alphaStd'] = self.config.parameters.alphaStd
            lutConfig['airmassStd'] = self.config.parameters.airmassStd
            lutConfig['lambdaRange'] = self.config.parameters.lambdaRange
            lutConfig['lambdaStep'] = self.config.parameters.lambdaStep
            lutConfig['lambdaNorm'] = self.config.parameters.lambdaNorm

        return lutConfig

    def _loadThroughputs(self, camera, opticsHandle, sensorHandleDict, filterHandleDict):
        """Internal method to load throughput data for filters

        Parameters
        ----------
        camera: `lsst.afw.cameraGeom.Camera`
            Camera from the butler
        opticsHandle : `lsst.daf.butler.DeferredDatasetHandle`
            Reference to optics transmission curve.
        sensorHandleDict : `dict` of [`int`, `lsst.daf.butler.DeferredDatasetHandle`]
            Dictionary of references to sensor transmission curves.  Key will
            be detector id.
        filterHandleDict : `dict` of [`str`, `lsst.daf.butler.DeferredDatasetHandle`]
            Dictionary of references to filter transmission curves.  Key will
            be physical filter label.

        Raises
        ------
        ValueError : Raised if configured filter name does not match any of the
            available filter transmission curves.
        """
        self._opticsTransmission = opticsHandle.get()

        self._sensorsTransmission = {}
        for detector in camera:
            self._sensorsTransmission[detector.getId()] = sensorHandleDict[detector.getId()].get()

        self._filtersTransmission = {}
        for physicalFilter in self.config.physicalFilters:
            self._filtersTransmission[physicalFilter] = filterHandleDict[physicalFilter].get()

    def _getThroughputDetector(self, detector, physicalFilter, throughputLambda):
        """Internal method to get throughput for a detector.

        Returns the throughput at the center of the detector for a given filter.

        Parameters
        ----------
        detector: `lsst.afw.cameraGeom._detector.Detector`
           Detector on camera
        physicalFilter: `str`
           Physical filter label
        throughputLambda: `np.array(dtype=np.float64)`
           Wavelength steps (Angstrom)

        Returns
        -------
        throughput: `np.array(dtype=np.float64)`
           Throughput (max 1.0) at throughputLambda
        """

        c = detector.getCenter(afwCameraGeom.FOCAL_PLANE)
        c.scale(1.0/detector.getPixelSize()[0])  # Assumes x and y pixel sizes in arcsec are the same

        throughput = self._opticsTransmission.sampleAt(position=c,
                                                       wavelengths=throughputLambda)

        throughput *= self._sensorsTransmission[detector.getId()].sampleAt(position=c,
                                                                           wavelengths=throughputLambda)

        throughput *= self._filtersTransmission[physicalFilter].sampleAt(position=c,
                                                                         wavelengths=throughputLambda)

        # Clip the throughput from 0 to 1
        throughput = np.clip(throughput, 0.0, 1.0)

        return throughput

    def _makeLutSchema(self, physicalFilterString, stdPhysicalFilterString,
                       atmosphereTableName):
        """
        Make the LUT schema

        Parameters
        ----------
        physicalFilterString: `str`
           Combined string of all the physicalFilters
        stdPhysicalFilterString: `str`
           Combined string of all the standard physicalFilters
        atmosphereTableName: `str`
           Name of the atmosphere table used to generate LUT

        Returns
        -------
        lutSchema: `afwTable.schema`
        """

        lutSchema = afwTable.Schema()

        lutSchema.addField('tablename', type=str, doc='Atmosphere table name',
                           size=len(atmosphereTableName))
        lutSchema.addField('elevation', type=float, doc="Telescope elevation used for LUT")
        lutSchema.addField('physicalFilters', type=str, doc='physicalFilters in LUT',
                           size=len(physicalFilterString))
        lutSchema.addField('stdPhysicalFilters', type=str, doc='Standard physicalFilters in LUT',
                           size=len(stdPhysicalFilterString))
        lutSchema.addField('pmb', type='ArrayD', doc='Barometric Pressure',
                           size=self.fgcmLutMaker.pmb.size)
        lutSchema.addField('pmbFactor', type='ArrayD', doc='PMB scaling factor',
                           size=self.fgcmLutMaker.pmb.size)
        lutSchema.addField('pmbElevation', type=np.float64, doc='PMB Scaling at elevation')
        lutSchema.addField('pwv', type='ArrayD', doc='Preciptable Water Vapor',
                           size=self.fgcmLutMaker.pwv.size)
        lutSchema.addField('o3', type='ArrayD', doc='Ozone',
                           size=self.fgcmLutMaker.o3.size)
        lutSchema.addField('tau', type='ArrayD', doc='Aerosol optical depth',
                           size=self.fgcmLutMaker.tau.size)
        lutSchema.addField('lambdaNorm', type=np.float64, doc='AOD wavelength')
        lutSchema.addField('alpha', type='ArrayD', doc='Aerosol alpha',
                           size=self.fgcmLutMaker.alpha.size)
        lutSchema.addField('zenith', type='ArrayD', doc='Zenith angle',
                           size=self.fgcmLutMaker.zenith.size)
        lutSchema.addField('nCcd', type=np.int32, doc='Number of CCDs')

        # and the standard values
        lutSchema.addField('pmbStd', type=np.float64, doc='PMB Standard')
        lutSchema.addField('pwvStd', type=np.float64, doc='PWV Standard')
        lutSchema.addField('o3Std', type=np.float64, doc='O3 Standard')
        lutSchema.addField('tauStd', type=np.float64, doc='Tau Standard')
        lutSchema.addField('alphaStd', type=np.float64, doc='Alpha Standard')
        lutSchema.addField('zenithStd', type=np.float64, doc='Zenith angle Standard')
        lutSchema.addField('lambdaRange', type='ArrayD', doc='Wavelength range',
                           size=2)
        lutSchema.addField('lambdaStep', type=np.float64, doc='Wavelength step')
        lutSchema.addField('lambdaStd', type='ArrayD', doc='Standard Wavelength',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('lambdaStdFilter', type='ArrayD', doc='Standard Wavelength (raw)',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('i0Std', type='ArrayD', doc='I0 Standard',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('i1Std', type='ArrayD', doc='I1 Standard',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('i10Std', type='ArrayD', doc='I10 Standard',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('i2Std', type='ArrayD', doc='I2 Standard',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('lambdaB', type='ArrayD', doc='Wavelength for passband (no atm)',
                           size=len(self.fgcmLutMaker.filterNames))
        lutSchema.addField('atmLambda', type='ArrayD', doc='Atmosphere wavelengths (Angstrom)',
                           size=self.fgcmLutMaker.atmLambda.size)
        lutSchema.addField('atmStdTrans', type='ArrayD', doc='Standard Atmosphere Throughput',
                           size=self.fgcmLutMaker.atmStdTrans.size)

        # and the look-up-tables
        lutSchema.addField('luttype', type=str, size=20, doc='Look-up table type')
        lutSchema.addField('lut', type='ArrayF', doc='Look-up table for luttype',
                           size=self.fgcmLutMaker.lut['I0'].size)

        return lutSchema

    def _makeLutCat(self, lutSchema, physicalFilterString, stdPhysicalFilterString,
                    atmosphereTableName):
        """
        Make the LUT schema

        Parameters
        ----------
        lutSchema: `afwTable.schema`
           Lut catalog schema
        physicalFilterString: `str`
           Combined string of all the physicalFilters
        stdPhysicalFilterString: `str`
           Combined string of all the standard physicalFilters
        atmosphereTableName: `str`
           Name of the atmosphere table used to generate LUT

        Returns
        -------
        lutCat: `afwTable.BaseCatalog`
           Lut catalog for persistence
        """

        # The somewhat strange format is to make sure that
        # the rows of the afwTable do not get too large
        # (see DM-11419)

        lutCat = afwTable.BaseCatalog(lutSchema)
        lutCat.table.preallocate(14)

        # first fill the first index
        rec = lutCat.addNew()

        rec['tablename'] = atmosphereTableName
        rec['elevation'] = self.fgcmLutMaker.atmosphereTable.elevation
        rec['physicalFilters'] = physicalFilterString
        rec['stdPhysicalFilters'] = stdPhysicalFilterString
        rec['pmb'][:] = self.fgcmLutMaker.pmb
        rec['pmbFactor'][:] = self.fgcmLutMaker.pmbFactor
        rec['pmbElevation'] = self.fgcmLutMaker.pmbElevation
        rec['pwv'][:] = self.fgcmLutMaker.pwv
        rec['o3'][:] = self.fgcmLutMaker.o3
        rec['tau'][:] = self.fgcmLutMaker.tau
        rec['lambdaNorm'] = self.fgcmLutMaker.lambdaNorm
        rec['alpha'][:] = self.fgcmLutMaker.alpha
        rec['zenith'][:] = self.fgcmLutMaker.zenith
        rec['nCcd'] = self.fgcmLutMaker.nCCD

        rec['pmbStd'] = self.fgcmLutMaker.pmbStd
        rec['pwvStd'] = self.fgcmLutMaker.pwvStd
        rec['o3Std'] = self.fgcmLutMaker.o3Std
        rec['tauStd'] = self.fgcmLutMaker.tauStd
        rec['alphaStd'] = self.fgcmLutMaker.alphaStd
        rec['zenithStd'] = self.fgcmLutMaker.zenithStd
        rec['lambdaRange'][:] = self.fgcmLutMaker.lambdaRange
        rec['lambdaStep'] = self.fgcmLutMaker.lambdaStep
        rec['lambdaStd'][:] = self.fgcmLutMaker.lambdaStd
        rec['lambdaStdFilter'][:] = self.fgcmLutMaker.lambdaStdFilter
        rec['i0Std'][:] = self.fgcmLutMaker.I0Std
        rec['i1Std'][:] = self.fgcmLutMaker.I1Std
        rec['i10Std'][:] = self.fgcmLutMaker.I10Std
        rec['i2Std'][:] = self.fgcmLutMaker.I2Std
        rec['lambdaB'][:] = self.fgcmLutMaker.lambdaB
        rec['atmLambda'][:] = self.fgcmLutMaker.atmLambda
        rec['atmStdTrans'][:] = self.fgcmLutMaker.atmStdTrans

        rec['luttype'] = 'I0'
        rec['lut'][:] = self.fgcmLutMaker.lut['I0'].flatten()

        # and add the rest
        rec = lutCat.addNew()
        rec['luttype'] = 'I1'
        rec['lut'][:] = self.fgcmLutMaker.lut['I1'].flatten()

        derivTypes = ['D_PMB', 'D_LNPWV', 'D_O3', 'D_LNTAU', 'D_ALPHA', 'D_SECZENITH',
                      'D_PMB_I1', 'D_LNPWV_I1', 'D_O3_I1', 'D_LNTAU_I1', 'D_ALPHA_I1',
                      'D_SECZENITH_I1']
        for derivType in derivTypes:
            rec = lutCat.addNew()
            rec['luttype'] = derivType
            rec['lut'][:] = self.fgcmLutMaker.lutDeriv[derivType].flatten()

        return lutCat
