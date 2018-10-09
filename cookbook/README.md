FGCM Cookbook
=============

This is a simple cookbook/example for running FGCM on the HSC RC dataset on the
lsst-dev server.  Please see [the FGCM README](../README.md) for installation
instructions.

Requirements
------------

This cookbook assumes access to
[lsst-dev](https://developer.lsst.io/services/lsst-dev.html) and the existence
of the HSC RC dataset from PDR1.  This is regularly reprocessed, starting with
[DM-10404](https://jira.lsstcorp.org/browse/DM-10404).  See Epic "Dataset
Reprocessing Campaigns (FY18b-2)"
[DM-14950](https://jira.lsstcorp.org/browse/DM-14950) for the list of tickets
with the most recent reprocessing.  On `lsst-dev` the recent reprocessing can
be listed at `/datasets/hsc/repo/rerun/RC/`.  However, I do recommend checking
the tickets at the linked Epic to confirm that the most recent processing
actually has been completed.

Setup
-----

The environment should be set up as follows:

```
setup lsst_distrib
setup -j -r /path/to/thirdparty/fgcm/.
setup -j -r /path/to/lsst-dm/fgcmcal/.

export RCRERUN=RC/w_2018_38/DM-15690
export COOKBOOKRERUN=fgcm_cookbook_w_2018_38
```

The `RCRERUN` env variable should be set to the most recent completed rerun
(see above).  The `COOKBOOKRERUN` env variable can be set to something
different if you have already run the cookbook with an older version of the
code and don't want to delete the previous run.

Finally, there is the expectation below that the `USER` environment variable is
set to your username (bash default), and that your private rerun directory is
`/datasets/hsc/repo/rerun/private/${USER}`.


Running FGCM
------------

There are four stages to running an FGCM calibration on a dataset:

1. Construct an atmosphere and instrument look-up table (LUT)
2. Ingest and match individual star observations from visit/ccd sets
3. Run the calibration for 3-4 "fit cycles", each of which contains 25-100 fit
iterations.
4. One more task will output the per-visit atmospheres for use downstream.

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

A [sample config](fgcmMakeLutHscFromTable.py) for HSC in in the `cookbook/`
directory.  There are only 3 fields to know about.  `filterNames` is a list of
the "canonical" filter names.  Multiple "filters" can reference a single "band"
(for example as the `r` and `i` filters have been replaced on HSC).  All filters
that map to a single band will be calibrated together and put on the same
system.  The `stdFilterNames` field is a list which matches each filter to a
"standard filter".  For example, you can specify that all observations taken
with the `i` filter should be transformed to the `i2` filter.  Note that at the
current time the RC dataset on `lsst-dev` does not contain any `i2` data.  Finally,
`atmosphereTableName` specifies the name of the table distributed with the
third-party FGCM code; see that code for details.

### Constructing a LUT using MODTRAN4

A [sample config](fgcmMakeLutHscFromModtran.py) for HSC is also in
the `cookbook/` directory.  This allows you to specify the range of atmosphere
parameters in the table (pressure, water vapor, ozone, aerosol tau and alpha,
and airmass) as well as the "standard atmosphere" parameters which should be
representative of the "typical" conditions.  It is not my expectation that this
will be run very often, as it will be more convenient to decide on parameters
and compute the atmosphere table once.

### Running `fgcmMakeLut.py`

This is a very simple command-line task (assuming that your username is in the
environment variable ${USER}).  E.g.,

```bash
fgcmMakeLut.py /datasets/hsc/repo --rerun \
${RCRERUN}:private/${USER}/${COOKBOOKRERUN}/lut --configfile \
$FGCMCAL_DIR/cookbook/fgcmMakeLutHscFromTable.py
```

## Ingesting and Matching Star Observations

In order to make the fit cycles tractable without redoing processor-intensive
steps, all data are collated and star observations are matched and indexed
before the fit runs are started.  Depending on how many visits are being
processed, and where the data are located, this step can take quite a while.
(In the future, with a database backend for the butler, I believe this step
will be much faster).

The FGCM code can run on a list of visits specified by the `--id` parameter on
the command line, or it can search the input repository for all visits with
`src` catalogs generated from `calexps`.  Be aware that FGCM uses **entire
visits** even if only one CCD is specified.  Also note that input processing
will be much faster if you specify a single CCD with the visits if `--id` is
used.  E.g., `--id visit=13376^13450 ccd=13`.  For HSC, using a reference CCD
to scan for available `src` catalogs speeds things up by x100, which is
necessary.  You can also specify by field or other valid dataId. A [sample
config](fgcmBuildStarsHsc.py) for HSC is available.

### Running `fgcmBuildStars.py`

This is also a simple command-line task.  E.g.,

```bash
fgcmBuildStars.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/lut:private/${USER}/${COOKBOOKRERUN}/wide \
--configfile $FGCMCAL_DIR/cookbook/fgcmBuildStarsHsc.py --id field=SSP_WIDE ccd=13 \
filter=HSC-G^HSC-R^HSC-I^HSC-Z^HSC-Y
```

## Running a Fit Cycle

### The First Fit Cycle (Cycle 0)

See the [sample config](fgcmFitCycleHsc_cycle00_config.py) for a sample config
file.  In this config is a place to describe the observational epochs (divided
by MJD) which should correspond to changes in observing (years, camera
warm-ups, etc), preferably coordinated with when new flat-fields were
generated.  At the moment flats are regenerated for HSC every couple of weeks
which is too frequent for a good fit from FGCM.  This is being
investigated. There is also a place to describe when the mirror was washed (in
MJD), as the system throughput usually jumps at these times.

There are three important differences between the first fit cycle and subsequent
cycles.  The first is that a "bright observation" algorithm is employed in the
first cycle to choose approximately photometric observations, while for
subsequent iterations the fit parameters from the previous iterations are used
to find photometric observations.  The second is that it is recommended that
you freeze the atmosphere to the "standard parameters" for the first fit
cycle.  The third is that for HSC it can be useful to set
`precomputeSuperStarInitialCycle` which can build a "super-star flat" based on
the bright observation selection.

At the end of the first (and all subsequent) cycles a whole bunch of diagnostic
plots are made in a subdirectory of the current directory where you invoked the
FGCM fit.  The directory name and plot filenames are generated from the
`outfileBase` and the `cycleNumber`.  Finally, a new config file is output
for the next fit cycle that automatically increments the `cycleNumber` and
turns off `freezeStdAtmosphere`.

When running the output is logged to stdout, though it's useful to capture the
output into a file with `tee`.  Note that if you want to start a new fit with
different parameters with the same stars/LUT then you can simply specify a new
output rerun and go from there.

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/wide:private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile $FGCMCAL_DIR/cookbook/fgcmFitCycleHscCookbook_cycle00_config.py \
|& tee ${COOKBOOKRERUN}_cycle00.log
```

### Subsequent Fit Cycles

Before running any subsequent cycle, you should look at the diagnostic plots.
The most important plots are `cycleNUMBER_expgray_BAND.png` plots.  Use these
histograms to choose appropriate cuts for each band to select "photometric"
exposures in the next cycle with the `expGrayPhotometricCut` (negative side)
and `expGrayHighCut` (positive side) parameters.  You want to capture the core
and reject outliers, but not too tight that you lose a lot of visits and can't
constrain the fit.  You can also up the number of iterations per fit cycle.  In
my experience, the fit does not improve if you go beyond ~50 iterations.  The
best way to get the fit to improve is to remove non-photometric exposures.

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile fgcmFitCycleHscCookbook_cycle01_config.py |& tee \
${COOKBOOKRERUN}_cycle01.log
```

### Final fit cycle

After the user had concluded that the fit has converged to her satisfaction,
one last run should be made with zero iterations and with `outputStandards =
True`.  This will ensure the final products have consistently applied superstar
flats and internal aperture corrections.

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile fgcmFitCycleHscCookbook_cycle02_config.py \
--config maxIter=0 --config outputStandards=True |& tee \
${COOKBOOKRERUN}_cycle02.log
```

## Outputs

The FGCM fitter spits out a number of outputs.  The first set is the diagnostic
plots which are put in a subdirectory of the current directory where the
command-line task is run.  The second set is the data products with zeropoint
and atmosphere information that is stored in the output repository.

### Diagnostic Plots

There are a wide range of diagnostic plots that get output.  This is a brief
explanation of what's in each one.  For details on how quantities are computed,
please see the [FGCM paper](http://adsabs.harvard.edu/abs/2018AJ....155...41B).

* `_airmass_expgray_BAND.png`: Average exposure gray residual (EXP^gray) as a
  function of airmass.  Used to check for airmass residuals.
* `_UT_expgray_BAND.png`: Average EXP^gray as a function of UT.  Used to check
  for residuals through the night.
* `_expgray_BAND.png`: Histogram of all EXP^gray with Gaussian fit.  Used to
  determine photometric cuts for next fit cycle.
* `_initial_expgray_BAND.png`: Histogram of all EXP^gray estimated from
  bright-observation algorithm (cycle 0 only).  Used for initial rough
  determination of photoemtric cuts for 0th fit cycle.
* `_expgray_redblue_compare_BAND.png`: Difference in gray residual between
  reddest and bluest 25% of stars on each exposure.  Used to check for
  chromatic residuals.
* `_i1r1_BAND.png`: Chromatic slope factor retrieved from stars (R1) vs predicted
  from FGCM fit (I1), for each visit/ccd.  Used to confirm filter/ccd
  throughputs.
* `_I1_BAND.png`: Map of producted chromatic slope (I1) from input throughput.
* `_R1_BAND.png`: Map of retrieved chromatic slope (R1) from stars.
* `_R1-I1_BAND.png`: Residual between R1 - I1 over the field.  Used to check for
  outliers which are usually caused by incorrect system throughput predictions.
* `_mjd_deep_expgray.png`: Not currently used for FGCM in the stack.
* `_pwv_vs_rpwv.png`: Model PWV vs "retrieved" PWV from z-band colors
  (experimental).
* `_pwv_vs_rpwv_scaled.png`: Model PWV vs scaled retrieved PWV from z-band colors
  (experimental).
* `_rpwv_vs_rpwv_in.png`: Comparison of retrieved PWV from z-band colors from
  previous cycle to present cycle (experimental).
* `_rpwv_vs_rpwv_smooth.png`: Comparison of retrieved PWV to smoothed retrieved
  PWV (experimental).
* `_sigfgcm_all_BAND.png`: Four panels showing individual star repeatability for
  all colors; blue/middle/red 25%/50%/25% of stars.  Used to determine
  intrinsic calibration error.
* `_sigfgcm_reserved_BAND.png`: Four panels showing individual star repeatability
  for reserved stars not used in fit.
* `_sigfgcm_reseved_crunched_BAND.png`: Four panels showing individual star
  repeatability after "gray crunch" for reserved stars.  Used to estimate final
  output repeatability.
* `_superstar_BAND_EPOCH.png`: Two panels showing the overall super-star flat
  illumination correction (left), and residual from last fit cycle (right).
  Used to check for superstar convergence.
* `_nightly_PAR.png`: Average nightly parameter value for alpha, o3, pwv, tau.
* `_qesys_washes.png`: Overall model for systematic dimming.
* `_zeropoints.png`: All zeropoints as a function of MJD from first exposure of
  calibration run.

## Outputs and Data Products

FGCM has both internal data products that are produced for each fit cycle, as
well as the ability to output photometric calibration data in a way that can be
used by downstream processing in the stack.

### Making Stack Output Products

The final task can output atmosphere transmission curves derived from the model
fit for each visit (`doAtmosphereOutput'); output a reference catalog of
"standard stars" that may be used for single-ccd calibration
(`doRefcatOutput`); and `jointcal_photoCalib` files for each visit/ccd
(`doZeropointOutput`).  All of these options are on by default.

One key thing is that FGCM does not use any absolute reference sources, so the
overall throughput for each band is unconstrained.  Therefore, to get the
output reference catalog and zeropoints on a physical scale, we must use
calspec standards (preferred, and deferred to the future) or a reference
catalog.  At the moment, this uses the same PS1 catalog that is used as a
reference for the `PhotoCal` and `jointcal` tasks.  This is set with
`doReferenceCalibration`, and the default is on.  Please see the [sample
config](fgcmOutputProductsHsc.py) for example usage of setting up the reference
catalog.

To build the output products you use the `fgcmOutputProducts.py` command-line
task.  The only config variable is the cycle number that the user has deemed to
be converged and should be output.

```bash
fgcmOutputProducts.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/wide:private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile $FGCMCAL_DIR/cookbook/fgcmOutputProductsHsc.py \
--config cycleNumber=2 |& tee ${COOKBOOKRERUN}_output.log
```

In the output repo
`/datasets/hsc/repo/rerun/private/${USER}/${COOKBOOKRERUN}/fit1/transmission/` you will now
find a whole bunch of `atmosphere_VISIT.fits` files of the
`lsst.afw.image.TransmissionCurve` type.  Note that even if you have frozen the
atmosphere parameters to the "standard" values, these will be computed
specifically for the airmass of each individual observation.

In the output repo
`/datasets/hsc/repo/rerun/private/erykoff/rc2_w_2018_32/fit1/ref_cats/fgcm_stars/'
you will find a sharded reference catalog suitable to use for any observations
overlapping the survey images that have been calibrated.

And in the output repo
`/datasets/hsc/repo/rerun/private/erykoff/rc2_w_2018_32/fit1/jointcal-results`
you will find all of the `jointcal_photoCalib` spatially-variable zeropoint
files that can be used in coaddition.  At the moment, these combine information
from the WCS Jacobian distortions as well as the spatially variable
transmission (due to the illumination correction derived by FGCM).  In the
future, the Jacobian information will be factored out.

Note that all the files are set to tract `0000`, as FGCM has no concept of tracts
and patches.

### FGCM-Specific Data Products

The data products that may be of use to the end user, and grabbed from the butler,
e.g.:

```python
zeropoints = butler.get('fgcmZeropoints', dataId={'fgcmcycle': last_cycle_run})
```

* `fgcmZeropoints-%(fgcmcycle)02d.fits`: Table of zeropoints (`FGCM_ZPT`) and
  first-order chromatic corrections (`FGCM_I10`) per visit/ccd.  Quality of
  exposure is given by `FGCM_FLAG`, see FGCM paper for details.
* `fgcmAtmosphereParameters-%(fgcmcycle)02d.fits`: Atmosphere model parameters
  (per visit) that can be used as input to MODTRAN or `fgcm` to generate atmosphere
  transmission.
* Other data products are used as an input to the next fit cycle.
