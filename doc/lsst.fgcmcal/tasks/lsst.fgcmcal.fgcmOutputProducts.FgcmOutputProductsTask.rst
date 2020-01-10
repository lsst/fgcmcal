.. lsst-task-topic:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

######################
FgcmOutputProductsTask
######################

``FgcmOutputProductsTask`` uses the output from :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask` to generate a full suite of output products (photometric calibration files, ``fgcm_photoCalib``; atmosphere transmissions, ``transmission_atmosphere_fgcm``; and standard star calibrated reference catalogs ``ref_cat``) for use in downstream processing.

This is the fourth and final task in a typical ``fgcmcal`` processing chain. The first is :doc:`lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`, the second is :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`, and the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`.

``FgcmOutputProductsTask`` is available as a :ref:`command-line task <pipe-tasks-command-line-tasks>`, :command:`fgcmOutputProducts.py`.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-summary:

Processing summary
==================

``FgcmOutputProductsTask`` reads in outputs from :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, with the ``cycleNumber`` specified in the config file, translates these tables to formats used by coaddition and other downstream processing.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-cli:

fgcmOutputProducts.py command-line interface
============================================

Note that no ``--id`` arguments are used by ``fgcmOutputProducts.py``.

.. code-block:: text

   fgcmOutputProducts.py REPOPATH [@file [@file2 ...]] [--output OUTPUTREPO | --rerun RERUN] [--config cycleNumber=finalCycleNumber] [other options]

Key arguments:

:option:`REPOPATH`
   The input Butler repository's URI or file path.

Key options:

:option:`--config cycleNumber=finalCycleNumer`:
   The ``cycleNumber`` used for the ``isFinalCycle`` run of :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`

.. seealso::

   See :ref:`command-line-task-argument-reference` for details and additional options.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-butler:

Butler datasets
===============

When run as the ``fgcmOutputProducts.py`` command-line task, or directly through the `~lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask.runDataRef` method, ``FgcmOutputProductsTask`` obtains datasets from the input Butler data repository and persists outputs to the output Butler data repository.
Note that configurations for ``FgcmOutputProductsTask``, and its subtasks, affect what datasets are persisted and what their content is.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-butler-inputs:

Input datasets
--------------

``fgcmStandardStars``
    Catalog of standard stars from fit.
``fgcmZeropoints``
    Catalog of zero-point information.
``fgcmAtmosphereParameters``
    Catalog of atmosphere parameters.

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-butler-outputs:

Output datasets
---------------

``fgcm_stars``
    Reference catalog of standard stars.
``fgcm_photoCalib``
    One ``fgcm_photoCalib`` photometric calibration file is output for each visit / ccd.
``transmission_atmosphere_fgcm``
    One atmospheric transmission curve is output for each visit.

The reference catalog (``fgcm_stars``) is an ``lsst.afw.table.SimpleCatalog`` that can be used as a reference catalog, including ``band_flux`` and ``band_fluxErr`` measurements per band (in ``nJy``), as well as additional information including ``band_nGood`` (the number of photometric observations per band), ``band_nTotal`` (the total number of observations per band), and ``band_nPsfCandidate`` (the number of observations per band that were flagged as PSF candidates).

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask

.. _lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask-examples:

Examples
========

See the `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ for worked examples.
