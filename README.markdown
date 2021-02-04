FAUST Line Imaging README
=========================
This document accompanies the Python file `faust_imaging.py`, which is a
framework for imaging ALMA FAUST spectral line data. Because the script is
long, I recommend using an editor that can 'fold' code at the function level.
To run this script in CASA (v5), run in the CASA prompt:

```python
# using an absolute path
execfile('/PATH/TO/faust_imaging.py')
# or a relative path
execfile('../casa_scripts/faust_imaging.py')
```


Paths and directories
---------------------
This script currently requires an explicit directory structure to read the
measurement sets and where to write the images. Note the global paths:

```
DATA_DIR -- where the measurement sets are stored
PROD_DIR -- where CASA will be run and the output products stored
```

Directories in `PROD_DIR` will be automatically created to store images,
moment maps, and plots.

Visualizing the directory structure as a tree would have files organized as
such:


```
/PATH/TO/DATA_DIR/   <- location on user's system
  - CB68-Setup1      <- MS files per target per setup
  - CB68-Setup2
  - ...
/PATH/TO/PROD_DIR/   <- location on user's system
  - images/          <- these will be made automatically
    + CB68/          <- image cubes and products per target
    + L1527/
    + ...
  - moments/
    + CB68/          <- moment maps per target
    + L1527/
    + ...
  - plots/           <- diagnostic and QA plots
```

Some CASA tasks do not work with absolute paths, so please note that (as
currently written) this CASA script must be run from the directory specified by
`PROD_DIR`.  The MS files should be present in the directories:

```bash
$DATA_DIR/<TARGET>-Setup1/
$DATA_DIR/<TARGET>-Setup2/
$DATA_DIR/<TARGET>-Setup2/
```

where `<TARGET>` is the FAUST target "field" name, e.g. "CB68". The names may
be found in the `ALL_TARGETS` dictionary, I retrieved these values from the
proposal, so they may be inconsistent for some targets. The paths above may
modified by editing the format string attribute `DataSet.ms_fmt`.


Running the pipeline
--------------------
To run the imaging pipeline for all SPWs of a target, use the `run_pipeline`
function.

```python
# for just Setup 1
run_pipeline('CB68', setup=1)
# for all setups (default)
run_pipeline('CB68', setup=None)
```

A list of uv-weightings to use and whether to image the fullcube or a windowed
sub-cube may be passed.

```python
# image both Briggs robust=0.5 and Natural uv-weightings
# image a 20km/s window centered on the priamry line instead of all channels
run_pipeline('CB68', setup=1, weightings=(0.5, 'natural'), fullcube=False)
```

Please see the docstring for further configuration information. For
imaging individual SPWs or more fine-grained control of the pipeline,
please see the following section.


Making images
-------------
The `ImageConfig` class abstracts properties specific to a desired image
configuration. Please see the docstring (either in the script or with
`ImageConfig?` from the prompt) for more information. After executing the
script using `execfile`, an instance may created by:

```python
config = ImageConfig.from_name('CB68', '244.936GHz_CS', fullcube=True, weighting=0.5)
config.run_pipeline()
```

The full list of SPW labels may be found in the `ALL_SPW_LABELS` variable.
The steps of the pipeline may be called individually for further control.

```python
config.make_dirty_cube()
config.clean_line_nomask(sigma=4.5)
config.make_seed_mask(sigma=5.0)
config.clean_line(mask_method='seed+multithresh', ext='clean')
config.postprocess(ext='clean', make_fits=True)
```

If run with `mpicasa`, `tclean` will be run in parallel, otherwise it will run
normally (i.e., single-threaded). At the current the method of seeding the
masking with auto-multithresh (`mask_method='seed+multithresh'`) does not
currently work with parallel `tclean`.

