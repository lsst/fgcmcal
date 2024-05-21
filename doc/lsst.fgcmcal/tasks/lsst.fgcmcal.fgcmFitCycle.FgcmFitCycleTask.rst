.. lsst-task-topic:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

################
FgcmFitCycleTask
################

``FgcmFitCycleTask`` is the task to run the ``fgcm`` fit from star observations generated from :doc:`lsst.fgcmcal.fgcmBuildFromIsolatedStars.FgcmBuildFromIsolatedStarsTask` or :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask` using the look-up table generated from :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`.
This code can either be run one "fit cycle" at a time (useful for development/testing) or all together (if ``config.doMultipleCycles = True``.
The final "fit cycle" is a cleanup that does not fit atmosphere or model parameters, but instead generates the final tables required for input to :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

This task will produce a large number of QA plots which will be put into the butler with names following the template ``fgcm_Cycle{N}_{Name}_Plot``.
This template makes it easy to see all the plots from a single fit cycle sorted together.
By default, only QA plots from the final two fit cycles (the next-to-last cycle is a full fit and the last cycle is the post-fit cleanup) will be output, unless ``config.doPlotsBeforeFinalCycles = True``.

This is the third task in a typical ``fgcmcal`` processing chain.  The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the second is :doc:`lsst.fgcmcal.fgcmBuildFromIsolatedStars.FgcmBuildFromIsolatedStarsTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-summary:

Processing summary
==================

``FgcmFitCycleTask`` reads in the star observations and the look-up table, performs an atmosphere and instrument fit, and outputs fit parameters as well as a comprehensive list of QA plots.
If the config option ``config.isFinalCycle = True`` then additional datafiles are output that are used by :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.
This is set automatically when ``config.doMultipleCycles = True``.

In single-cycle mode, ``FgcmFitCycleTask`` will output a new config file in the current working directory with recommended settings for the subsequent fit cycle.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-examples:
