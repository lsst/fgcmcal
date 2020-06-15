.. lsst-task-topic:: lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask

###############
FgcmMakeLutTask
###############

``FgcmMakeLutTask`` computes a look-up table tracking atmosphere and instrumental variations for use in the ``fgcmcal`` fits.  The intention is that this is run once for a given observatory/instrument and should only be updated when instrumental parameters change (e.g., updated filter and ccd throughput curves), or in the extremely rare case if an observatory is relocated to a different elevation or latitude.

This is the first task in a typical ``fgcmcal`` processing chain.  The second is :doc:`lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask` or :doc:`lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`, the third is :doc:`lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`, and the fourth is :doc:`lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`.

``FgcmMakeLutTask`` is available as a :ref:`command-line task <pipe-tasks-command-line-tasks>`, :command:`fgcmMakeLut.py`.

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-summary:

Processing summary
==================

``FgcmMakeLutTask`` uses an input atmosphere table (preferable) or a list of atmosphere parameters, combined with the instrumental throughput as a function of position, to compute a look-up table used in the ``fgcmcal`` fits.

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-cli:

fgcmMakeLut.py command-line interface
=====================================

Note that no ``--id`` arguments are used by ``fgcmMakeLut.py``.

.. code-block:: text

   fgcmMakeLut.py REPOPATH [@file [@file2 ...]] [--output OUTPUTREPO | --rerun RERUN] [other options]

Key arguments:

:option:`REPOPATH`
   The input Butler repository's URI or file path.

.. seealso::

   See :ref:`command-line-task-argument-reference` for details and additional options.

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-butler:

Butler datasets
===============

When run as the ``fgcmMakeLut.py`` command-line task, or directly through the `~lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask.runDataRef` method, ``FgcmMakeLutTask`` obtains datasets from the input Butler data repository and persists outputs to the output Butler data repository.
Note that configurations for ``FgcmMakeLutTask``, and its subtasks, affect what datasets are persisted and what their content is.

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-butler-inputs:

Input datasets
--------------

``camera``
   Camera geometry and detector object
``transmission_optics``
   Optics transmission curve for the instrument
``transmission_filter``
   Filter transmission curve (as a function of position) for the instrument

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-butler-outputs:

Output datasets
---------------

``fgcmLookUpTable``
   FGCM atmosphere and instrument look-up table

.. _lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask-subtasks:

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
