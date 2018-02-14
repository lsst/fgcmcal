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
a simple `setup lsst_apps`.

To install to use with the LSST stack, you simply need to run `setup -j -r
/path/to/fgcm/.` which will read in the associated `eups` table and you should
be good to go.

Installing the FGCM LSST Package
--------------------------------

Also, at the moment, the FGCM LSST package is not distributed with the LSST
stack, but it will be.  Installation is just like above.  First, download the
code from https://github.com/lsst-dm/fgcm .  To install, you simply need to run
`setup -j -r /path/to/lsst-dm/fgcm/.` (after doing `setup lsst_apps`) which
will read in the associated `eups` table and you should be good to go.  

Running FGCM
------------

There are three stages to running an FGCM calibration on a dataset:

1. Construct an atmosphere and instrument look-up table (LUT)
2. Ingest and match individual star observations from visit/ccd sets
3. Run the calibration for 3-4 "fit cycles", each of which contains 25-100 fit
iterations.

## Constructing a Look-Up Table (LUT)

The LUT is constructed via the `fgcmMakeLut.py` command line task.
If you are satisfied with the defaults for your telescope, you can use a
precomputed atmosphere table.  If you are not, you either need to install MODTRAN4 (not
free, unfortunately), or contact erykoff@stanford.edu.

The LUT depends on both the atmosphere parameters (which can be precomputed
generically), and the instrumental parameters (which will change as better
metrology becomes available, filters are replaced, etc).  (Though hopefully
filters aren't replaced that often.)

### Constructing a LUT from a precomputed atmosphere table

A [sample config](./examples/fgcmMakeLutHscFromTable.py) for HSC is in the
`examples/` directory.  There are only 3 fields to know about.  `filterNames`
is a list of the "short" filter names (though I don't know the technical
term).  Multiple filters can reference a single band (for example as the r and
i filters have been replaced on HSC).  The `stdFilterNames` field is a list
which matches each filter to a "standard filter".  For example, you can specify
that all observations taken with the i filter should be transformed to the i2
filter.  Finally, `atmosphereTableName` specifies the name of the table
distributed with the third-party FGCM code; see that code for details.

### Constructing a LUT using MODTRAN4

A [sample config](./examples/fgcmMakeLutHscFromModtran.py) for HSC is also in
the `examples/` directory.  This allows you to specify the range of atmosphere
parameters in the table (pressure, water vapor, ozone, aerosol tau and alpha,
and airmass) as well as the "standard atmosphere" parameters which should be
representative of the "typical" conditions.  It is not my expectation that this
will be run very often, as it will be more convenient to decide on parameters
and run

### Running `fgcmMakeLut.py`

This is a very simple command-line task.  E.g.,

`fgcmMakeLut.py /path/to/input/repo/with/all/the/data --output
/path/to/output/repo/for/fgcm/calibrations --configfile fgcmMakeLutHscFromTable.py`

## Ingesting and Matching Star Observations

In order to make the fit cycles tractable without redoing processor-intensive
steps, all data are collated and star observations are matched and indexed before the fit
runs are started.  Depending on how many visits are being processed, and where
the data are located, this step can take quite a while.  (In the future, with a
database backend for the butler, I believe this step will be much faster).

The FGCM code can run on a list of visits specified by the `--id` parameter on
the command line, or it can search the input repository for all visits with
`src` catalogs generated from `calexps`.  Be aware that FGCM uses **entire
visits** even if only one CCD is specified.  Also note that input processing
will be much faster if you specify a single CCD with the visits if `--id` is
used.  E.g., `--id visit=13376^13450 ccd=13`.  For HSC, using a reference CCD
to scan for available `src` catalogs speeds things up by x100, which is
necessary.  A [sample config](./examples/fgcmBuildStarsHsc.py) for HSC is available.

### Running `fgcmBuildStars.py`

This is a very simple command-line task.  E.g.,

`fgcmBuildStars.py /path/to/output/repo/for/fgcm/calibrations --configfile
fgcmBuildStarsHsc.py [--id @visits.txt]`

## Running a Fit Cycle

### The First Cycle (Cycle 0)

See the [sample config](./examples/fgcmFitCycleHsc.py) for a sample config
file.  In this config is a place to describe the observational epochs (divided by
MJD) which should correspond to changes in observing (years, camera warm-ups,
etc), preferably coordinated with when new flat-fields were generated.  There
is also a place to describe when the mirror was washed (in MJD), as the system
throughput usually jumps at these times.

There are two important differences between the first fit cycle and subsequent
cycles.  The first is that a "bright observation" algorithm is employed in the
first cycle to choose approximately photometric observations.  The second is
that it is recommended that you freeze the atmosphere to the "standard
parameters" for the first fit cycle.

At the end of the first (and all subsequent) cycles a whole bunch of diagnostic
plots are made in a subdirectory generated from the `outfileBase` and the
`cycleNumber`.  In addition, a new config file is output for the next fit cycle
that automatically increments the `cycleNumber` and turns off `freezeStdAtmosphere`.

`fgcmFitCycle.py /path/to/output/repo/for/fgcm/calibrations --configfile fgcmFitCycle_cycle00_config.fit`

### Subsequent Cycles

Before running any subsequent cycle, you should look at the diagnostic plots.
The most important plots are `cycleNUMBER_expgray_BAND.png` plots.  Use these
histograms to choose appropriate cuts for each band to select "photometric"
exposures in the next cycle with the `expGrayPhotometricCut` (negative side)
and `expGrayHighCut` (positive side) parameters.  You want to capture the core
and reject outliers, but not too tight that you lose a lot of visits and can't
constrain the fit.  You can also up the number of iterations per fit cycle.  In
my experience, the fit does not improve if you go beyond ~50 iterations.  The
best way to get the fit to improve is to remove non-photometric exposures.

`fgcmFitCycle.py /path/to/output/repo/for/fgcm/calibrations --configfile fgcmFitCycle_cycleXX_config.fit`

## Outputs

Need to write up on outputs.


