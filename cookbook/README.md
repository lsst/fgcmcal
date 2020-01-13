FGCM Cookbook
=============

This is a simple cookbook/example for running FGCM on the HSC RC dataset on the
lsst-dev server.  Please see [the FGCM README](../README.md) for general installation
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

export RCRERUN=RC/w_2019_38/DM-21386-sfm
export COOKBOOKRERUN=fgcm_cookbook_w_2019_38
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
environment variable ${USER}).  This task should take around 3.5 hours on
`lsst-dev01`.  The recommended command is:

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
(In the future, with a database backend for the Gen3 butler, I believe this step
will be much faster).

The FGCM code runs on calexp source catalogs from visits constrained by the
`--id` parameter on the command line.  Best results are obtained when FGCM is
run with **full visits**.  Due to limitations in the Gen2 Butler (the only
Butler currently supported by `fgcmcal`), optimal performance is obtained by
specifying a single "reference" ccd on the command line (e.g. `ccd=13`) and
setting the config variable `checkAllCcds = True` (which is the default).  The
alternative is to specify all the desired CCDs and set `checkAllCcds=False`,
e.g., "--ccd 0..8^10..103".  However, this is slower than the first option, and
the improvement in speed in the first option is greater the more visits are
specified.  If instead you want to process all the visits in a rerun selected
by filter, field, or some other dataid field, then by using a reference ccd and
setting `checkAllCcds=True` you can speed things up by a factor of
approximately 100 relative to the alternative (naming CCDs specifically).  For
config settings, please see the [sample config](fgcmBuildStarsHsc.py).

If `doReferenceCalibration = True` in the configuration (the default), then
stars from a reference catalog (e.g. PS1) will be loaded and matched to the
internal stars.  The signal-to-noise cut specified here should be the minimum
one expects to use in the actual fit, and the fit may also be performed without
the use of reference stars if desired.

### Running `fgcmBuildStars.py`

This is also a simple command-line task.  This task will take around 2.5 hours
on `lsst-dev01`:

```bash
fgcmBuildStars.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/lut:private/${USER}/${COOKBOOKRERUN}/wide+deep \
--configfile $FGCMCAL_DIR/cookbook/fgcmBuildStarsHsc.py \
--id ccd=13 filter=HSC-G^HSC-R^HSC-I^HSC-Z^HSC-Y
```

## Running a Fit Cycle

### The First Fit Cycle (Cycle 0)

See the [sample config](fgcmFitCycleHsc_cycle00_config.py) for a sample config
file.  In this config is a place to describe the observational epochs (divided
by MJD) which should correspond to changes in observing (years, camera
warm-ups, etc), preferably coordinated with when new flat-fields were
generated.  At the moment flats are regenerated for HSC every couple of weeks
which is too frequent for a good fit from FGCM.  This is being
investigated. There is also a place to describe when the mirror was washed or
recoated (in MJD), as the system throughput usually jumps at these times.

There are three important differences between the first fit cycle ("cycle 0")
and subsequent cycles.  The first is that a "bright observation" algorithm is
employed in the first cycle to choose approximately photometric observations,
while for subsequent iterations the best fit parameters from the previous
iterations are used to determine photometric observations.  The second is that
it is recommended that you freeze the atmosphere to the "standard parameters"
for the first fit cycle (the default).  The third is that for HSC it can be
useful to set `precomputeSuperStarInitialCycle` which can build a "super-star
flat" based on the bright observation selection to smooth out the large
illumination correction term.

At the end of the first (and all subsequent) cycles a whole bunch of diagnostic
plots are made in a subdirectory of the current directory where you invoked the
FGCM fit.  Therefore, you should run in a clean directory.  The directory name
and plot filenames are generated from the `outfileBase` and the `cycleNumber`.
Finally, a new config file is output for the next fit cycle that automatically
increments the `cycleNumber` and turns off `freezeStdAtmosphere`.

When running the output is logged to stdout, though it's useful to capture the
output into a file with `tee`.  Note that if you want to start a new fit with
different parameters with the same stars/LUT then you can simply specify a new
output rerun and go from there.  This stage should take less than 5 minutes to
run on `lsst-dev01`.

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/wide+deep:private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile $FGCMCAL_DIR/cookbook/fgcmFitCycleHscCookbook_cycle00_config.py \
|& tee fgcmFitCycleHscCookbook_cycle00.log
```

### Subsequent Fit Cycles

Before running subsequent fit cycles, it is helpful to look at the diagnostic
plots.  The most important plots are `cycleNUMBER_expgray_BAND.png` plots.  As
of [DM-17305](https://jira.lsstcorp.org/browse/DM-17305), the width of these
histograms are used to set a default photometricity cut for the next fit
cycle.  However, it is worth inspecting these histograms to ensure that they
look "normal".  The expectation is that they have a Gaussian core with an
excess non-photometric tail out to the negative side.  Any significant
bi-modality is a sign that something has gone horribly awry with the fits,
perhaps indicating issues with the input data or the configuration parameters.

The relevant configuration parameters that you might want to adjust based on
these plots are `expGrayPhotometricCut` (on the negative side) and
`expGrayHighCut` (on the positive side).  You want to capture the core and
reject outliers, but not too tight that you lose a lot of visits and can't
constrain the fit.

You can also increase the number of iterations per fit cycle.  In my
experience, the fit does not improve if you go beyond ~50 iterations with small
datasets, but large datasets can use ~100 iterations per cycle.

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile fgcmFitCycleHscCookbook_cycle01_config.py |& tee \
fgcmFitCycleHscCookbook_cycle01.log
```

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile fgcmFitCycleHscCookbook_cycle02_config.py |& tee \
fgcmFitCycleHscCookbook_cycle02.log
```


### Final fit cycle

After the user has concluded that the fit has converged to her satisfaction,
one last run should be made with `isFinalCycle=True`.  This will tell `fgcm` to
run one final selection of stars and photometric exposures, and will output
zeropoints and standard stars for use in `fgcmOutputProducts.py`.  Running the
final cycle ensures that the final products have consistently applied
selections, superstar flats, and internal aperture corrections.  The final
cycle run should only take 2-3 minutes.

```bash
fgcmFitCycle.py /datasets/hsc/repo --rerun private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile fgcmFitCycleHscCookbook_cycle03_config.py \
--config isFinalCycle=True |& tee \
fgcmFitCycleHscCookbook_cycle03.log
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
(`doRefcatOutput`); and `fgcm_photoCalib` files for each visit/ccd
(`doZeropointOutput`).  All of these options are on by default.

The FGCM calibration can be run with or without reference stars as an anchor.
The new default, including in this cookbook, is to run with reference sources
from PS1.  However, if you run without reference sources you get a good
relative calibraiton but the overall throughput for each band is
unconstrained.  To get the output reference catalog and zeropoints on a
physical scale, we must use calspec standards (which is preferred, but not
currently supported) or a reference catalog such as PS1.  To turn on this final
computation of the system throughput, set `doReferenceCalibration=True` in the
`FgcmOutputProducts` task.  (The default is now `False`).  Please see the [sample
config](fgcmOutputProductsHsc.py) for example usage of setting up the reference
catalog.

To build the output products you use the `fgcmOutputProducts.py` command-line
task.  The only config variable is the cycle number that the user has deemed to
be converged and should be output.  This should take around 5 minutes on `lsst-dev01`.

```bash
fgcmOutputProducts.py /datasets/hsc/repo --rerun \
private/${USER}/${COOKBOOKRERUN}/wide:private/${USER}/${COOKBOOKRERUN}/fit1 \
--configfile $FGCMCAL_DIR/cookbook/fgcmOutputProductsHsc.py \
--config cycleNumber=3 |& tee fgcmFitCycleHscCookbook_output.log
```

In the output repo
`/datasets/hsc/repo/rerun/private/${USER}/${COOKBOOKRERUN}/fit1/transmission/` you will now
find a whole bunch of `atmosphere_VISIT.fits` files of the
`lsst.afw.image.TransmissionCurve` type.  Note that even if you have frozen the
atmosphere parameters to the "standard" values, these will be computed
specifically for the airmass of each individual observation.

In the output repo
`/datasets/hsc/repo/rerun/private/${USER}/${COOKBOOKRERUN}/fit1/ref_cats/fgcm_stars/'
you will find a sharded reference catalog suitable to use for any observations
overlapping the survey images that have been calibrated.

And in the output repo
`/datasets/hsc/repo/rerun/private/${USER}/${COOKBOOKRERUN}/fit1/fgcmcal-results`
you will find all of the `fgcmcal_photoCalib` spatially-variable zeropoint
files that can be used in coaddition.  At the moment, these combine information
from the WCS Jacobian distortions as well as the spatially variable
transmission (due to the illumination correction derived by FGCM).  In the
future, the Jacobian information will be factored out.

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
