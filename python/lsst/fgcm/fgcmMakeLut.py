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

from .detectorThroughput import DetectorThroughput

import time

import fgcm

__all__ = ['FgcmMakeLutConfig', 'FgcmMakeLutTask']

class FgcmMakeLutConfig(pexConfig.Config):
    """Config for FgcmMakeLutTask"""

    filterNames = pexConfig.ListField(
        doc="Filter names to build LUT",
        dtype=str,
        default=("NO_DATA",),
        )
    elevation = pexConfig.Field(
        doc="Telescope elevation (m)",
        dtype=float,
        default=None,
        )
    pmbRange = pexConfig.ListField(
        doc="Barometric Pressure range (millibar)",
        dtype=float,
        default=(770.0,790.0,),
        )
    pmbSteps = pexConfig.Field(
        doc="Barometric Pressure number of steps",
        dtype=int,
        default=5,
        )
    pwvRange = pexConfig.ListField(
        doc="Precipitable Water Vapor range (mm)",
        dtype=float,
        default=(0.1,12.0,),
        )
    pwvSteps = pexConfig.Field(
        doc="Precipitable Water Vapor number of steps",
        dtype=int,
        default=21,
        )
    o3Range = pexConfig.ListField(
        doc="Ozone range (dob)",
        dtype=float,
        default=(220.0,310.0,),
        )
    o3Steps = pexConfig.Field(
        doc="Ozone number of steps",
        dtype=int,
        default=3,
        )
    tauRange = pexConfig.ListField(
        doc="Aerosol Optical Depth range (unitless)",
        dtype=float,
        default=(0.002,0.25,),
        )
    tauSteps = pexConfig.Field(
        doc="Aerosol Optical Depth number of steps",
        dtype=int,
        default=11,
        )
    alphaRange = pexConfig.ListField(
        doc="Aerosol alpha range (unitless)",
        dtype=float,
        default=(0.0,2.0),
        )
    alphaSteps = pexConfig.Field(
        doc="Aerosol alpha number of steps",
        dtype=int,
        default=9,
        )
    zenithRange = pexConfig.ListField(
        doc="Zenith angle range (degree)",
        dtype=float,
        default=(0.0,70.0,),
        )
    zenithSteps = pexConfig.Field(
        doc="Zenith angle number of steps",
        dtype=int,
        default=21,
        )
    pmbStd = pexConfig.Field(
        doc="Standard Atmosphere pressure (millibar)",
        dtype=float,
        default=None,
        )
    pwvStd = pexConfig.Field(
        doc="Standard Atmosphere PWV (mm)",
        dtype=float,
        default=None,
        )
    o3Std = pexConfig.Field(
        doc="Standard Atmosphere O3 (dob)",
        dtype=float,
        default=None,
        )
    tauStd = pexConfig.Field(
        doc="Standard Atmosphere aerosol optical depth",
        dtype=float,
        default=None,
        )
    alphaStd = pexConfig.Field(
        doc="Standard Atmosphere aerosol alpha",
        dtype=float,
        default=None,
        )
    airmassStd = pexConfig.Field(
        doc="Standard Atmosphere airmass",
        dtype=float,
        default=None,
        )
    lambdaNorm = pexConfig.Field(
        doc="Aerosol Optical Depth normalization wavelength (A)",
        dtype=float,
        default=None,
        )
    lambdaStep = pexConfig.Field(
        doc="Wavelength step for generating atmospheres (nm)",
        dtype=float,
        default=0.5,
        )
    lambdaRange = pexConfig.ListField(
        doc="Wavelength range for LUT (A)",
        dtype=float,
        default=None,
        )

    def setDefaults(self):
        pass

class FgcmMakeLutRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for fgcmMakeLutTask

    fgcmMakeLutTask.run() takes one argument, the butler, and
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
        if self.doReturnResults:
            return results

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
            profileName = parsedCmd.profile if hasattr(parsedCmd, "profile") else None
            log = parsedCmd.log
            targetList = self.getTargetList(parsedCmd)
            # make sure that we only get 1
            resultList = self(targetList[0])

        return resultList

class FgcmMakeLutTask(pipeBase.CmdLineTask):
    """
    Make Look-Up Table for FGCM global calibration
    """

    ConfigClass = FgcmMakeLutConfig
    RunnerClass = FgcmMakeLutRunner
    _DefaultName = "fgcmMakeLut"

    def __init__(self, butler=None, **kwargs):
        """
        Instantiate an fgcmMakeLutTask.

        Parameters
        ----------
        butler : lsst.daf.persistence.Butler
        """

        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""

        parser = pipeBase.ArgumentParser(name=cls._DefaultName)

        return parser

    # no saving of metadata for now
    def _getMetadataName(self):
        return None

    @pipeBase.timeMethod
    def run(self, butler):
        """
        Make a Look-Up Table for FGCM

        Parameters
        ----------
        butler:  lsst.daf.persistence.Butler

        Returns
        -------
        Empty list
        """

        if (not butler.datasetExists('fgcmLookUpTable')):
            self._fgcmMakeLut(butler)

        return []

    def _fgcmMakeLut(self, butler):
        """
        Make a FGCM Look-up Table

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
           (used for mapper information)
        """

        # need the camera for the detectors
        camera = butler.get('camera')

        # number of ccds from the length of the camera iterator
        nCcd = len(camera)
        self.log.info("Found %d ccds for look-up table" % (nCcd))

        # make the config dictionary
        lutConfig = {}
        lutConfig['logger'] = self.log
        lutConfig['elevation'] = self.config.elevation
        lutConfig['filterNames'] = self.config.filterNames
        lutConfig['nCCD'] = nCcd
        lutConfig['pmbRange'] = self.config.pmbRange
        lutConfig['pmbSteps'] = self.config.pmbSteps
        lutConfig['pwvRange'] = self.config.pwvRange
        lutConfig['pwvSteps'] = self.config.pwvSteps
        lutConfig['o3Range'] = self.config.o3Range
        lutConfig['o3Steps'] = self.config.o3Steps
        lutConfig['tauRange'] = self.config.tauRange
        lutConfig['tauSteps'] = self.config.tauSteps
        lutConfig['alphaRange'] = self.config.alphaRange
        lutConfig['alphaSteps'] = self.config.alphaSteps
        lutConfig['zenithRange'] = self.config.zenithRange
        lutConfig['zenithSteps'] = self.config.zenithSteps
        lutConfig['pmbStd'] = self.config.pmbStd
        lutConfig['pwvStd'] = self.config.pwvStd
        lutConfig['o3Std'] = self.config.o3Std
        lutConfig['tauStd'] = self.config.tauStd
        lutConfig['alphaStd'] = self.config.alphaStd
        lutConfig['airmassStd'] = self.config.airmassStd
        lutConfig['lambdaRange'] = self.config.lambdaRange
        lutConfig['lambdaStep'] = self.config.lambdaStep
        lutConfig['lambdaNorm'] = self.config.lambdaNorm

        # make the lut object
        self.log.info("Making the LUT maker object")
        self.fgcmLutMaker = fgcm.FgcmLUTMaker(lutConfig)

        # generate the throughput dictionary.  Fun!
        # do this internally here at first.  Later, break it into its own thing

        # these will be in Angstroms
        # note that lambdaStep is currently in nm, because dumb.  convert to A
        throughputLambda = np.arange(self.config.lambdaRange[0],
                                     self.config.lambdaRange[1]+self.config.lambdaStep*10,
                                     self.config.lambdaStep*10.)

        self.log.info("Built throughput lambda, %.1f-%.1f, step %.2f" %
              (throughputLambda[0], throughputLambda[-1], throughputLambda[1]-throughputLambda[0]))

        tput = DetectorThroughput()

        throughputDict = {}
        for i,filterName in enumerate(self.config.filterNames):
            tDict = {}
            tDict['LAMBDA'] = throughputLambda
            for ccdIndex,detector in enumerate(camera):
                # make sure we convert the calling units from A to nm
                tDict[ccdIndex] = tput.getThroughputDetector(detector, filterNames,
                                                             throughputLambda/10.)
            throughputDict[filterName] = tDict

        # set the throughputs
        self.fgcmLutMaker.setThroughputs(throughputDict)

        # make the LUT
        self.log.info("Making LUT")
        self.fgcmLutMaker.makeLUT()

        # and save the LUT
        lutSchema = afwTable.Schema()

        # new version, which gets around the afwTable row length limitation
        #  each LUT will be saved in a different row
        #  there is overhead of the arrays that we only need one copy, but this
        #  is going to be insignificant overall

        # build the index values
        comma=','
        filterNameString = comma.join(self.config.filterNames)

        lutSchema.addField('filternames', type=str, doc='filterNames in LUT',
                           size=len(filterNameString))
        lutSchema.addField('pmb', type='ArrayD', doc='Barometric Pressure',
                           size=self.fgcmLutMaker.pmb.size)
        lutSchema.addField('pmbfactor', type='ArrayD', doc='PMB scaling factor',
                           size=self.fgcmLutMaker.pmb.size)
        lutSchema.addField('pmbelevation', type=np.float64, doc='PMB Scaling at elevation')
        lutSchema.addField('pwv', type='ArrayD', doc='Preciptable Water Vapor',
                           size=self.fgcmLutMaker.pwv.size)
        lutSchema.addField('o3', type='ArrayD', doc='Ozone',
                           size=self.fgcmLutMaker.o3.size)
        lutSchema.addField('tau', type='ArrayD', doc='Aerosol optical depth',
                           size=self.fgcmLutMaker.tau.size)
        lutSchema.addField('lambdanorm', type=np.float64, doc='AOD wavelength')
        lutSchema.addField('alpha', type='ArrayD', doc='Aerosol alpha',
                           size=self.fgcmLutMaker.alpha.size)
        lutSchema.addField('zenith', type='ArrayD', doc='Zenith angle',
                           size=self.fgcmLutMaker.zenith.size)
        lutSchema.addField('nccd', type=np.int32, doc='Number of CCDs')

        # and the standard values
        lutSchema.addField('pmbstd', type=np.float64, doc='PMB Standard')
        lutSchema.addField('pwvstd', type=np.float64, doc='PWV Standard')
        lutSchema.addField('o3std', type=np.float64, doc='O3 Standard')
        lutSchema.addField('taustd', type=np.float64, doc='Tau Standard')
        lutSchema.addField('alphastd', type=np.float64, doc='Alpha Standard')
        lutSchema.addField('zenithstd', type=np.float64, doc='Zenith angle Standard')
        lutSchema.addField('lambdarange', type='ArrayD', doc='Wavelength range',
                           size=2)
        lutSchema.addField('lambdastep', type=np.float64, doc='Wavelength step')
        lutSchema.addField('lambdastd', type='ArrayD', doc='Standard Wavelength',
                           size=self.fgcmLutMaker.filterNames.size)
        lutSchema.addField('i0std', type='ArrayD', doc='I0 Standard',
                           size=self.fgcmLutMaker.filterNames.size)
        lutSchema.addField('i1std', type='ArrayD', doc='I1 Standard',
                           size=self.fgcmLutMaker.filterNames.size)
        lutSchema.addField('i10std', type='ArrayD', doc='I10 Standard',
                           size=self.fgcmLutMaker.filterNames.size)
        lutSchema.addField('lambdab', type='ArrayD', doc='Wavelength for passband (no atm)',
                           size=self.fgcmLutMaker.filterNames.size)
        lutSchema.addField('atmlambda', type='ArrayD', doc='Atmosphere wavelengths',
                           size=self.fgcmLutMaker.atmLambda.size)
        lutSchema.addField('atmstdtrans', type='ArrayD', doc='Standard Atmosphere Throughput',
                           size=self.fgcmLutMaker.atmStd.size)

        # and the look-up-tables
        lutSchema.addField('luttype', type=str, size=20, doc='Look-up table type')
        lutSchema.addField('lut', type='ArrayF', doc='Look-up table for luttype',
                           size=self.fgcmLutMaker.lut['I0'].size)

        lutCat = afwTable.BaseCatalog(lutSchema)
        lutCat.table.preallocate(14)

        # first fill the first index
        rec=lutCat.addNew()

        rec['filterNames'] = filterNameString
        rec['pmb'][:] = self.fgcmLutMaker.pmb
        rec['pmbfactor'][:] = self.fgcmLutMaker.pmbFactor
        rec['pmbelevation'] = self.fgcmLutMaker.pmbElevation
        rec['pwv'][:] = self.fgcmLutMaker.pwv
        rec['o3'][:] = self.fgcmLutMaker.o3
        rec['tau'][:] = self.fgcmLutMaker.tau
        rec['lambdanorm'] = self.fgcmLutMaker.lambdaNorm
        rec['alpha'][:] = self.fgcmLutMaker.alpha
        rec['zenith'][:] = self.fgcmLutMaker.zenith
        rec['nccd'] = self.fgcmLutMaker.nCCD

        rec['pmbstd'] = self.fgcmLutMaker.pmbStd
        rec['pwvstd'] = self.fgcmLutMaker.pwvStd
        rec['o3std'] = self.fgcmLutMaker.o3Std
        rec['taustd'] = self.fgcmLutMaker.tauStd
        rec['alphastd'] = self.fgcmLutMaker.alphaStd
        rec['zenithstd'] = self.fgcmLutMaker.zenithStd
        rec['lambdarange'][:] = self.fgcmLutMaker.lambdaRange
        rec['lambdastep'] = self.fgcmLutMaker.lambdaStep
        rec['lambdastd'][:] = self.fgcmLutMaker.lambdaStd
        rec['i0std'][:] = self.fgcmLutMaker.I0Std
        rec['i1std'][:] = self.fgcmLutMaker.I1Std
        rec['i10std'][:] = self.fgcmLutMaker.I10Std
        rec['lambdab'][:] = self.fgcmLutMaker.lambdaB
        rec['atmlambda'][:] = self.fgcmLutMaker.atmLambda
        rec['atmstdtrans'][:] = self.fgcmLutMaker.atmStdTrans

        rec['luttype'] = 'I0'
        rec['lut'][:] = self.fgcmLutMaker.lut['I0'].flatten()

        # and add the rest
        rec = lutCat.addNew()
        rec['luttype'] = 'I1'
        rec['lut'][:] = self.fgcmLutMaker.lut['I1'].flatten()

        derivTypes = ['D_PMB', 'D_PWV', 'D_O3', 'D_LNTAU', 'D_ALPHA', 'D_SECZENITH',
                      'D_PMB_I1', 'D_PWV_I1', 'D_O3_I1', 'D_LNTAU_I1', 'D_ALPHA_I1', 'D_SECZENITH_I1']
        for derivType in derivTypes:
            rec = lutCat.addNew()
            rec['luttype'] = derivType
            rec['lut'][:] = self.fgcmLutMaker.lutDeriv[derivType].flatten()

        butler.put(lutCat, 'fgcmLookUpTable')

        # and we're done
