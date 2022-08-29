.. py:currentmodule:: lsst.fgcmcal

.. _lsst.fgcmcal:

############
lsst.fgcmcal
############

The ``lsst.fgcmcal`` module runs the Forward Global Calibration Method (FGCM) to perform global photometric survey calibration for the Rubin Observatory LSST.
Please see `Burke, Rykoff, et al. 2018 <http://adsabs.harvard.edu/abs/2018AJ....155...41B>`_ for the paper describing the method.  This ``lsst.fgcmcal`` package wraps the third-party package `fgcm <https://github.com/lsst/fgcm/tree/lsst-dev>`_.

.. _lsst.fgcmcal.pythononly-using:

Using lsst.fgcmcal
==================

Please see the `cookbook <https://github.com/lsst/fgcmcal/tree/master/cookbook/README.md>`_ for a runthrough on how to use ``lsst.fgcmcal``.

There are four tasks to be run in a typical global ``fgcmcal`` processing chain.  They are:

#. Make a look-up table: :doc:`tasks/lsst.fgcmcal.fgcmMakeLut.FgcmMakeLutTask`

#. Build the star lists: :doc:`tasks/lsst.fgcmcal.fgcmBuildStarsTable.FgcmBuildStarsTableTask` (if ``sourceTable_visit`` parquet tables are available) or :doc:`tasks/lsst.fgcmcal.fgcmBuildStars.FgcmBuildStarsTask`

#. Run the fitter: :doc:`tasks/lsst.fgcmcal.fgcmFitCycle.FgcmFitCycleTask`

#. Output the final products: :doc:`tasks/lsst.fgcmcal.fgcmOutputProducts.FgcmOutputProductsTask`

Alternatively, ``fgcmcal`` can be run on a single tract (with multi-band coverage), although the results will not be as robust as a full global calibration.  This can be run with :doc:`tasks/lsst.fgcmcal.fgcmCalibrateTractTable.FgcmCalibrateTractTableTask` or :doc:`tasks/lsst.fgcmcal.fgcmCalibrateTract.FgcmCalibrateTractTask`, which will run all of the tasks above except for the making of the look-up table.

.. _lsst.fgcmcal.pythononly-contributing:

Contributing
============

``lsst.fgcmcal`` is developed at https://github.com/lsst/fgcmcal.  You can find Jira issues for this module under the `fgcmcal <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20fgcmcal>`_ component.

.. _lsst.fgcmcal.pythononly-pipeline-tasks:

Pipeline tasks
--------------

.. lsst-pipelinetasks::
   :root: lsst.fgcmcal
   :toctree: tasks

.. _lsst.example.pythononly-tasks:

Configurations
--------------

.. lsst-configs::
   :root: lsst.fgcmcal
   :toctree: configs

.. _lsst.example.pythononly-pyapi:

Python API reference
====================

.. automodapi:: lsst.fgcmcal
   :no-main-docstr:
   :no-inheritance-diagram:
