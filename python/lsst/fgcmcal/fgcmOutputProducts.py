# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import sys
import traceback

import numpy as np

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import TransmissionCurve
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask

import fgcm

__all__ = ['FgcmOutputProductsConfig', 'FgcmOutputProductsTask']


class FgcmOutputProductsConfig(pexConfig.Config):
    """Config for FgcmOutputProductsTask"""

    cycleNumber = pexConfig.Field(
        doc="Fit cycle as basis for products to output",
        dtype=int,
        default=None,
    )
    doReferenceCalibration = pexConfig.Field(
        doc="Apply 'absolute' calibration from reference catalog",
        dtype=bool,
        default=False,
    )
    photoRefObjLoader = pexConfig.ConfigurableField(
        default=LoadIndexedReferenceObjectsTask,
        doc="reference object loader for 'absolute' photometric calibration",
    )
    photoCal = pexConfig.ConfigurableField(
        default=PhotoCalTask,
        doc="perform 'absolute' calibration",
    )

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
        self.photoCal.applyColorTerms = True
        self.photoCal.fluxField = 'slot_calibFlux_flux'  # NOT THIS
        # also need to set _flags for this ...


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
        butler: lsst.daf.persistence.Butler

        Returns
        -------
        None if self.doReturnResults is False
        An empty list if self.doReturnResults is True
        """
        task = self.TaskClass(config=self.config, log=self.log)

        exitStatus = 0
        if self.doRaise:
            results = task.run(butler)
        else:
            try:
                results = task.run(butler)
            except Exception as e:
                exitStatus = 1
                task.log.fatal("Failed: %s" % e)
                if not isinstance(e, pipeBase.TaskError):
                    traceback.print_exc(file=sys.stderr)

        task.writeMetadata(butler)

        if self.doReturnResults:
            # Not much in terms of results to return (empty list)
            return [pipeBase.Struct(exitStatus=exitStatus,
                                    results=results)]
        else:
            return [pipeBase.Struct(exitStatus=exitStatus)]

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
            # profileName = parsedCmd.profile if hasattr(parsedCmd, "profile") else None
            # log = parsedCmd.log
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

        # Check to make sure that the fgcmBuildStars config exists, to retrieve
        # the visit and ccd dataset tags
        if not butler.datasetExists('fgcmBuildStars_config'):
            raise RuntimeError("Cannot find fgcmBuildStars_config, which is prereq for fgcmOutputProducts")
        if not butler.datasetExists('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber):
            raise RuntimeError("The specified fgcm cycle (%d) has not been run on this repo." %
                               (self.config.cycleNumber))

        self._outputAtmospheres(butler)

        return []

    def _outputAtmospheres(self, butler):
        """
        Output the atmospheres.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
           (used for mapper information)
        """

        # First, we need to grab the look-up table and key info
        lutCat = butler.get('fgcmLookUpTable')

        atmosphereTableName = lutCat[0]['tablename']
        elevation = lutCat[0]['elevation']
        atmLambda = lutCat[0]['atmlambda']
        lutCat = None

        # Get the fgcmBuildStarsConfig for visit name
        fgcmBuildStarsConfig = butler.get("fgcmBuildStars_config")
        visitDataRefName = fgcmBuildStarsConfig.visitDataRefName

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
            except (ValueError, IOError):
                raise RuntimeError("FGCM look-up-table generated without modtran, "
                                   "but modtran not configured to run.")

        # Next, we need to grab the atmosphere parameters
        atmCat = butler.get('fgcmAtmosphereParameters', fgcmcycle=self.config.cycleNumber)

        zenith = np.degrees(np.arccos(1. / atmCat['seczenith']))

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
                       dataId={visitDataRefName: visit})
