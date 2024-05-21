.. lsst-task-topic:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

#######################
FgcmBuildStarsTableTask
#######################

``FgcmBuildStarsTableTask`` finds all the single-visit sources in a repository (or a subset based on command-line parameters) from ``sourceTable_visit`` parquet tables and extracts all the potential photometric calibration stars for input into fgcm.
This task additionally uses fgcm to match star observations into unique stars, and performs as much cleaning of the input catalog as possible.

This is the second task in a typical ``fgcmcal`` processing chain.
The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.
This task is still usable and tested, but the replacement :doc:`lsst.fgcmcal.fgcmBuildFromIsolatedStars.FgcmBuildFromIsolatedStarsTask` is preferred for modern pipelines.

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTableTask-summary:

Processing summary
==================

``FgcmBuildStarsTableTask`` runs this sequence of operations:

#. Finds unique visits and collates visit metadata, including exposure time, pointing, typical psf size, background level.

#. Reads in all sources, selecting good stars according to the chosen source selector.

#. Matches sources internally across bands to get a unique multi-band list of possible calibration stars.

#. Matches possible calibration stars to a reference catalog.

#. All results are stored in the output repo ``fgcm-process`` directory.

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-examples:
