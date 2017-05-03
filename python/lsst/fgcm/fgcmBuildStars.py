# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import lsst.utils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pex.exceptions as pexExceptions
import lsst.afw.table

from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

import lsst.fgcm


__all__ = ['FgcmBuildStarsConfig','FgcmBuildStarsTask']

class FgcmBuildStarsConfig(pexConfig.Config):
    """Config for FgcmBuildStarsTask"""

    minPerBand = pexConfig.Field(
        doc="Minimum observations per band"
        dtype=int,
        default=2,
        )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching",
        )

    def setDefaults(self):
        sourceSelector = self.sourceSelector
        sourceSelector.setDefaults()
        sourceSelector.sourceFluxType = 'Calib'


class FgcmBuildStarsTask(pipeBase.CmdLineTask):
    """
    Build stars for the FGCM global calibration
    """

    ConfigClass = FgcmBuildStarsConfig
    #RunnerClass = FgcmBuildStarsRunner
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

