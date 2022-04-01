.. lsst-task-topic:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

################
FgcmFitCycleTask
################

``FgcmFitCycleTask`` runs a single ``fgcm`` fit on the star observations generated from :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask` using the look-up table generated from :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`.  This code is meant to be run multiple times until convergence, and the results output from one "fit cycle" are used as an input to the subsequent fit cycle.  One final cleanup run is performed to generate the tables required for input to :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

This is the third task in a typical ``fgcmcal`` processing chain.  The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the second is :doc:`lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask` or :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-summary:

Processing summary
==================

``FgcmFitCycleTask`` reads in the star observations and the look-up table, performs an atmosphere and instrument fit, and outputs fit parameters as well as a comprehensive list of QA plots.  If the config option ``isFinalCycle`` is set to ``True`` then additional datafiles are output that are used by :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

``FgcmFitCycleTask`` will output a new config file in the current working directory with recommended settings for the subsequent fit cycle.  Furthermore, there are a wide range of diagnostic QA plots that are output by the task.  For details on the contents of these QA plots, please see the  `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ as well as `Burke, Rykoff, et al. 2018 <http://adsabs.harvard.edu/abs/2018AJ....155...41B>`_.

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
