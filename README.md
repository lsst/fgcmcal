Forward Global Calibration Method (FGCM)
========================================

This is LSST stack code to interface with the FGCM package.  See
http://adsabs.harvard.edu/abs/2018AJ....155...41B for the paper describing the
method, and https://github.com/erykoff/fgcm for a Python implementation that is
used by this LSST package.

FGCM performs global survey calibration with forward modeling of the atmosphere
and instrument, including chromatic corrections.  Currently, this LSST package
will only run on HSC data within the stack, but this is being actively worked
on.  In general, the FGCM code should be compatible with any survey with the
following requirements:

* Exposure/CCD based observations
* Transmission curves for each filter, preferably for each CCD
* MODTRAN4 is required to generate new atmosphere tables.
* Atmosphere tables for the following telescope locations are included:
    - Blanco telescope at CTIO (DES)
    - Subaru telescope at Mauna Kea (HSC)
    - LSST telescope at Cerro Pachon (LSST)
* Enough memory to hold all the observations in memory at once.
    - A full run of four years of DES survey data can be run on a machine
      with 128 Gb RAM and 32 processors in less than a day.

Installing the FGCM Python Package
----------------------------------

At the moment, the FGCM python package is not distributed with the LSST stack
as a third-party package, but it is my intention that it will be.  Installation
is simple.  First, download the code from https://github.com/erykoff/fgcm .
Note that all of the [required
dependencies](https://github.com/erykoff/fgcm#dependencies) are satisfied with
a simple `setup lsst_apps` or `setup lsst_distrib`.

To install to use with the LSST stack, you simply need to run `setup -j -r
/path/to/thirdparty/fgcm/.` which will read in the associated `eups` table and you should
be good to go.

Installing the FGCM LSST Package (`fgcmcal`)
--------------------------------------------

Also, at the moment, the FGCM LSST package is not distributed with the LSST
stack, but it will be.  Installation is just like above.  First, download the
code from https://github.com/lsst-dm/fgcmcal .  To install, you first need to
run `scons` as:

```
setup lsst_distrib
cd /path/to/lsst-dm/fgcmcal
scons
```

Then, to set it up for use you simply need to run `setup -j -r
/path/to/lsst-dm/fgcmcal/.`  This will read in the associated `eups` table.

FGCM Cookbook
-------------

Please see [the FGCM Cookbook](cookbook/README.md) for doing a test run on HSC
RC data on `lsst-dev`.

