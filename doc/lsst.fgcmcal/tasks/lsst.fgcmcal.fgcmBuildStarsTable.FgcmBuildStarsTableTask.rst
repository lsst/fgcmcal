.. lsst-task-topic:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

#######################
FgcmBuildStarsTableTask
#######################

``FgcmBuildStarsTableTask`` finds all the single-visit sources in a repository (or a subset based on command-line parameters) from ``sourceTable_visit`` parquet tables and extracts all the potential photometric calibration stars for input into fgcm.
This task additionally uses fgcm to match star observations into unique stars, and performs as much cleaning of the input catalog as possible.
A tutorial on the steps of running ``fgcmcal`` are found in the `cookbook`_.

This ``fgcmcal`` task runs on ``sourceTable_visit`` catalogs from visits constrained by the ``--id`` parameter on the command line.

At the current time, ``fgcmcal`` does not support Gen3.

This is the second task in a typical ``fgcmcal`` processing chain.
The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

``FgcmBuildStarsTableTask`` is available as a :ref:`command-line task <pipe-tasks-command-line-tasks>`, :command:`fgcmBuildStarsTable.py`.

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTableTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``FgcmBuildStarsTableTask`` runs this sequence of operations:

#. Finds unique visits and collates visit metadata, including exposure time, pointing, typical psf size, background level.

#. Reads in all sources, selecting good stars according to the chosen source selector.

#. Matches sources internally across bands to get a unique multi-band list of possible calibration stars.

#. Matches possible calibration stars to a reference catalog.

#. All results are stored in the output repo ``fgcm-process`` directory.


.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-cli:

fgcmBuildStarsTable.py command-line interface
=============================================

.. code-block:: text

   fgcmBuildStarsTable.py REPOPATH [@file [@file2 ...]] [--output OUTPUTREPO | --rerun RERUN] [--id] [other options]

Key arguments:

:option:`REPOPATH`
   The input Butler repository's URI or file path.

Key options:

:option:`--id`:
   The data IDs to process.

.. seealso::

   See :ref:`command-line-task-argument-reference` for details and additional options.

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-butler:

Butler datasets
===============

When run as the ``fgcmBuildStarsTable.py`` command-line task, or directly through the `~lsst.fgcmcal.FgcmBuildStarsTableTask.runDataRef` method, ``FgcmBuildStarsTableTask`` obtains datasets from the input Butler data repository and persists outputs to the output Butler data repository.
Note that configurations for ``FgcmBuildStarsTableTask``, and its subtasks, affect what datasets are persisted and what their content is.

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-butler-inputs:

Input datasets
--------------

``sourceTable_visit``
    Full-depth source catalog, per-visit, in parquet format
``calexp`` (s)
    Calibrated exposures produced by `ProcessCcdTask` (for exposure metadata)
``fgcmLookupTable``
    FGCM look-up table produced by :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-butler-outputs:

Output datasets
---------------

``fgcmVisitCatalog``
    Catalog (`lsst.afw.table`) of visit metadata
``fgcmStarObservations``
    Catalog of star observations
``fgcmStarIds``
    Catalog of unique star ids, positions, and number of observations
``fgcmStarIndices``
    Catalog of indices linking unique star ids to star observations
``fgcmReferenceStars``
    Catalog of reference stars matched to unique star ids.

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask

.. _lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask-examples:

Examples
========

See the `cookbook`_ for worked examples.

.. _cookbook: https://github.com/lsst/fgcmcal/tree/master/cookbook/
