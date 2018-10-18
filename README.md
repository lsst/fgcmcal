fgcmcal: Global Photometric Calibration in LSST with FGCM
=========================================================

This is the LSST stack code to interface with the Forward Global Calibration
Method (FGCM) third-party package to perform global photometric survey
calibration.  Please see [Burke, Rykoff, et
al. 2018](http://adsabs.harvard.edu/abs/2018AJ....155...41B) for the paper
describing the method, and https://github.com/lsst/fgcm/tree/lsst-dev (the LSST
fork of https://github.com/erykoff/fgcm) for the Python implementation that is
used by this `fgcmcal` package.

FGCM performs a global photometric calibration, starting with instrumental
fluxes and producing top-of-the-atmosphere standard fluxes by forward modeling
the atmosphere and instrument, including chromatic corrections.  Overall
absolute calibration is reduced to the computation of one offset per band,
which `fgcmcal` can optionally compute based on a reference catalog.

Currently, the `fgcmcal` package will only run on HSC data within the stack,
but this is being actively worked on.  In general, the FGCM code should be
compatible with any survey with the following requirements:

* Visit/CCD based observations
* Transmission curves for each filter, preferably for each CCD
* MODTRAN4 is required to generate new atmosphere tables.
* Atmosphere tables for the following telescope locations are included:
    - Blanco telescope at CTIO (DES)
    - Subaru telescope at Mauna Kea (HSC)
    - LSST telescope at Cerro Pachon (LSST)
* Enough memory to hold all the observations in memory at once.
    - A full run of four years of DES survey data can be run on a machine
      with 128 Gb RAM and 32 processors in less than a day.


Installing the `fgcm` and `fgcmcal` Packages
--------------------------------------------

As of DM-16128, both `fgcm` (the third-party package) and `fgcmcal` (the stack
interface) are distributed with `lsst_distrib`.  Therefore, with a stack
installation after `w_2018_47` or v17.0 or later, setting up
`fgcmcal` is as simple as:

```
setup lsst_distrib
```

Or, to only get `fgcmcal` and required dependencies:

```
setup fgcmcal
```

FGCM Cookbook
-------------

Please see [the FGCM Cookbook](cookbook/README.md) for doing a test run on HSC
RC data on `lsst-dev` and to learn the control flow of `fgcmcal`.
