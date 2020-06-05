.. lsst-task-topic:: lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask

##################
FgcmBuildStarsTask
##################

``FgcmBuildStarsTask`` finds all the single-visit sources in a repository (or a subset based on command-line parameters) and extracts all the potential photometric calibration stars for input into fgcm.  This task additionally uses fgcm to match star observations into unique stars, and performs as much cleaning of the input catalog as possible.  A tutorial on the steps of running ``fgcmcal`` are found in the `cookbook`_.

The ``fgcmcal`` code runs on calexp source catalogs from visits constrained by the ``--id`` parameter on the command line.  Best results are obtained when ``fgcmcal`` is run with full visits.

In Gen2, due to limitations of the Gen2 Butler, optimal performance is obtained by specifying a single "reference" ccd on the command line (e.g. ``ccd=13``) and setting the config variable ``checkAllCcds = True`` (which is the default).  The alternative is to specify all the desired CCDs and set ``checkAllCcds = False``, e.g., ``ccd=0..8^10..103``.  However, this is slower than the first option, and the improvement in speed in the first option is greater the more visits are specified.  If instead you want to process all the visits in a rerun selected by filter, field, or some other dataid field, then by using a reference ccd and setting ``checkAllCcds = True`` you can speed things up by a factor of approximately 100 relative to the alternative (naming CCDs specifically).

Be aware that if a visit does not have a ``calexp`` available with the given reference CCD then it will be skipped.  If this is possibly an issue, multiple reference ccds can be specified on the command line, although performance will degrade the more are specified.

At the current time, ``fgcmcal`` does not support Gen3.

This is the second task in a typical ``fgcmcal`` processing chain.  The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

``FgcmBuildStarsTask`` is available as a :ref:`command-line task <pipe-tasks-command-line-tasks>`, :command:`fgcmBuildStars.py`.

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-summary:

Processing summary
==================

.. If the task does not break work down into multiple steps, don't use a list.
.. Instead, summarize the computation itself in a paragraph or two.

``FgcmBuildStarsTask`` runs this sequence of operations:

#. Finds unique visits and collates visit metadata, including exposure time, pointing, typical psf size, background level.

#. Reads in all sources, selecting good stars according to the chosen source selector.

#. Matches sources internally across bands to get a unique multi-band list of possible calibration stars.

#. Matches possible calibration stars to a reference catalog.

#. All results are stored in the output repo ``fgcm-process`` directory.


.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-cli:

fgcmBuildStars.py command-line interface
========================================

.. code-block:: text

   fgcmBuildStars.py REPOPATH [@file [@file2 ...]] [--output OUTPUTREPO | --rerun RERUN] [--id] [other options]

Key arguments:

:option:`REPOPATH`
   The input Butler repository's URI or file path.

Key options:

:option:`--id`:
   The data IDs to process.

.. seealso::

   See :ref:`command-line-task-argument-reference` for details and additional options.

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-butler:

Butler datasets
===============

When run as the ``fgcmBuildStars.py`` command-line task, or directly through the `~lsst.fgcmcal.FgcmBuildStarsTask.runDataRef` method, ``FgcmBuildStarsTask`` obtains datasets from the input Butler data repository and persists outputs to the output Butler data repository.
Note that configurations for ``FgcmBuildStarsTask``, and its subtasks, affect what datasets are persisted and what their content is.

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-butler-inputs:

Input datasets
--------------

``src``
    Full-depth source catalog (`lsst.afw.table`) produced by `ProcessCcdTask`
``calexp``
    Calibrated exposure produced by `ProcessCcdTask` (for exposure metadata)
``fgcmLookupTable``
    FGCM look-up table produced by :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-butler-outputs:

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

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask

.. _lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask-examples:

Examples
========

See the `cookbook`_ for worked examples.

.. _cookbook: https://github.com/lsst/fgcmcal/tree/master/cookbook/
