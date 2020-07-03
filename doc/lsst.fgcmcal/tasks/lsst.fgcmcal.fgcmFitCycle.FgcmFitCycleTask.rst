.. lsst-task-topic:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

################
FgcmFitCycleTask
################

``FgcmFitCycleTask`` runs a single ``fgcm`` fit on the star observations generated from :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask` using the look-up table generated from :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`.  This code is meant to be run multiple times until convergence, and the results output from one "fit cycle" are used as an input to the subsequent fit cycle.  One final cleanup run is performed to generate the tables required for input to :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

This is the third task in a typical ``fgcmcal`` processing chain.  The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the second is :doc:`lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask` or :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

``FgcmFitCycleTask`` is available as a :ref:`command-line task <pipe-tasks-command-line-tasks>`, :command:`fgcmFitCycle.py`.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-summary:

Processing summary
==================

``FgcmFitCycleTask`` reads in the star observations and the look-up table, performs an atmosphere and instrument fit, and outputs fit parameters as well as a comprehensive list of QA plots.  If the config option ``isFinalCycle`` is set to ``True`` then additional datafiles are output that are used by :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.  Please see the `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ for information on deciding when and how to set the ``isFinalCycle`` parameter to ``True``.

Aside from the ``butler`` outputs (listed below), the ``FgcmFitCycleTask`` will output a new config file in the current working directory with recommended settings for the subsequent fit cycle.  Furthermore, there are a wide range of diagnostic QA plots that are output by the task.  For details on the contents of these QA plots, please see the  `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ as well as `Burke, Rykoff, et al. 2018 <http://adsabs.harvard.edu/abs/2018AJ....155...41B>`_.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-cli:

fgcmFitCycle.py command-line interface
======================================

Note that no ``--id`` arguments are used by ``fgcmFitCycle.py``.

.. code-block:: text

   fgcmFitCycle.py REPOPATH [@file [@file2 ...]] [--output OUTPUTREPO | --rerun RERUN] [--configfile configfile] [other options]

Key arguments:

:option:`REPOPATH`
   The input Butler repository's URI or file path.

Key options:

:option:`--configfile`:
   The config file to use.  For all runs except the initial, this should be the config file output from the previous cycle.

.. seealso::

   See :ref:`command-line-task-argument-reference` for details and additional options.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-butler:

Butler datasets
===============

When run as the ``fgcmFitCycle.py`` command-line task, or directly through the `~lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask.runDataRef` method, ``FgcmFitCycleTask`` obtains datasets from the input Butler data repository and persists outputs to the output Butler data repository.
Note that configurations for ``FgcmFitCycleTask``, and its subtasks, affect what datasets are persisted and what their content is.

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-butler-inputs:

Input datasets
--------------

``camera``
    Camera geometry and detector object
``fgcmLookupTable``
    FGCM look-up table produced by :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`
``fgcmVisitCatalog``
    Catalog (`lsst.afw.table`) of visit metadata produced by :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`
``fgcmStarObservations``
    Catalog of star observations produced by :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`
``fgcmStarIds``
    Catalog of unique star ids, positions, and number of observations produced by :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`
``fgcmStarIndices``
    Catalog of indices linking unique star ids to star observations produced by :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`
``fgcmReferenceStars``
    Catalog of reference stars matched to unique star ids produced by :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`
``fgcmFitParameters``
    Catalog of fit parameters from previous fit cycle (if available)

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-butler-outputs:

Output datasets
---------------

``fgcmFitParameters``
    Catalog of fit parameters.  Not output if run as part of :doc:`lsst.fgcmcal.fgcmCalibrateTract.FgcmCalibrateTractTask`
``fgcmFlaggedStars``
    Catalog of flagged star ids, either bad or reserved from fit.   Not output if run as part of :doc:`lsst.fgcmcal.fgcmCalibrateTract.FgcmCalibrateTractTask`
``fgcmZeropoints``
    Catalog of zero-point information.  Only output if ``isFinalCycle`` is ``True``.  Not output if run as part of :doc:`lsst.fgcmcal.fgcmCalibrateTract.FgcmCalibrateTractTask`
``fgcmAtmosphereParameters``
    Catalog of atmosphere parameters.  Only output if ``isFinalCycle`` is ``True``.  Not output if run as part of :doc:`lsst.fgcmcal.fgcmCalibrateTract.FgcmCalibrateTractTask`
``fgcmStandardStars``
    Catalog of standard stars from fit.  Only output if ``isFinalCycle`` is ``True``.  Not output if run as part of :doc:`lsst.fgcmcal.fgcmCalibrateTract.FgcmCalibrateTractTask`

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask

.. _lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask-examples:

Examples
========

See the `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ for worked examples.
