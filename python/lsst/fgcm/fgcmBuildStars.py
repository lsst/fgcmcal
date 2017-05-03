# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import sys
import traceback

import lsst.utils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pex.exceptions as pexExceptions
import lsst.afw.table

import lsst.fgcm


__all__ = ['FgcmBuildStarsConfig','FgcmBuildStarsTask']

class FgcmBuildStarsConfig(pexConfig.Config):
    """Config for FgcmBuildStarsTask"""

    minPerBand = pexConfig.Field(
        doc="Minimum observations per band",
        dtype=int,
        default=2,
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
    def run(self, dataRefs):
        """
        Cross-match and make star list for FGCM

        Parameters
        ----------
        dataRefs: list of lsst.daf.persistence.ButlerDataRef
            List of data references to the exposures to be fit

        Returns
        -------
        pipe.base.Struct
            struct containing:
            * dataRefs: the provided data references consolidated
            (others)
        """

        if len(dataRefs) == 0:
            raise ValueError("Need a list of data references!")

        print("Run")
        print("min obs: %d" % (self.config.minPerBand))

