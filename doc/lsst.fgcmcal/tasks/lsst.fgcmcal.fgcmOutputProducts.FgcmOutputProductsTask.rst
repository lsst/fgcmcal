.. lsst-task-topic:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

######################
FgcmOutputProductsTask
######################

``FgcmOutputProductsTask`` uses the output from :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask` to generate a full suite of output products (photometric calibration files and atmosphere transmissions) for use in downstream processing.

This is the fourth and final task in a typical ``fgcmcal`` processing chain. The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the second is :doc:`lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask` or :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`, and the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-summary:

Processing summary
==================

``FgcmOutputProductsTask`` reads in outputs from :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, with the ``cycleNumber`` specified in the config file, translates these tables to formats used by coaddition and other downstream processing.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-examples:
