.. lsst-task-topic:: lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask

###############
FgcmMakeLutTask
###############

``FgcmMakeLutTask`` computes a look-up table tracking atmosphere and instrumental variations for use in the ``fgcmcal`` fits.  The intention is that this is run once for a given observatory/instrument and should only be updated when instrumental parameters change (e.g., updated filter and ccd throughput curves), or in the extremely rare case if an observatory is relocated to a different elevation or latitude.

This is the first task in a typical ``fgcmcal`` processing chain.
The second is :doc:`lsst.fgcmcal.fgcmBuildFromIsolatedStars.FgcmBuildFromIsolatedStarsTask`, the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-summary:

Processing summary
==================

``FgcmMakeLutTask`` uses an input atmosphere table (preferable) or a list of atmosphere parameters, combined with the instrumental throughput as a function of position, to compute a look-up table used in the ``fgcmcal`` fits.

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-butler:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-examples:

Examples
========

See the `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ for worked examples.
