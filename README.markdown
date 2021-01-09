FAUST Line Imaging README
=========================
This document accompanies the Python file `faust_imaging.py`, which is a
framework for imaging ALMA FAUST spectral line data. Because the script is
long, I recommend using an editor that can 'fold' code at the function level.
To run this script in CASA (v5), run in the CASA prompt:

```python
execfile('/PATH/TO/faust_imaging.py')
```

Note that relative paths are also valid:

```python
execfile('../casa_scripts/faust_imaging.py')
```


Paths and directories
---------------------
This script currently requires an explicit directory structure to read the
measurement sets and where to write the images. Note the global paths:

```
ROOT_DIR -- absolute path to workspace
DATA_DIR -- where CASA should be run
IMAG_DIR -- where output images will be written
MOMN_DIR -- where output moment maps will be written
```

Visualizing the directory structure as a tree would have files organized as
such:

```
/path/to/root/       <- path arbitrary
  - data/            <- run CASA from here
    + CB68-Setup1    <- create and put MS files in here
    + CB68-Setup2
    + ...
    + images/        <- these will be made automatically
      ~ CB68/
      ~ L1527/
      ~ ...
    + moments/
      ~ CB68/
      ~ L1527/
      ~ ...
```

Some CASA tasks do not work with absolute paths, so please note that (as
currently written) this CASA script must be run from the directory specified by
`DATA_DIR`.  The MS files should be present in the directories:

```bash
$DATA_DIR/<TARGET>-Setup1/
$DATA_DIR/<TARGET>-Setup2/
$DATA_DIR/<TARGET>-Setup2/
```

where `<TARGET>` is the FAUST target "field" name, e.g. "CB68". The names may
be found in the `ALL_TARGETS` dictionary, I retrieved these values from the
proposal, so they may be inconsistent for some targets. The paths above may
modified by editing the format string attribute `DataSet.ms_fmt`.


Making images
-------------
The `DataSet` class abstracts properties specific to a spectral setup. Please
see the docstring for more information. After executing the script using
`execfile`, an instance may created by:

```python
dset = DataSet('CB68', setup=2, kind='joint')
```

In the current scheme, dirty cubes must be made first to get estimates of the
RMS for each window and uv-weighting. To generate cut-out dirty cubes for all
lines with natural and robust=0 weighting, run:

```python
image_all_lines_inspect_dirty(dset)
```

If run with `mpicasa`, `tclean` will be run in parallel, otherwise it will run
normally (i.e., single-threaded). After these are created, auto-masking using
auto-multithresh can be used to create cleaned cubes by running:

```python
clean_all_lines_with_masking(dset, fullcube=False)
```

For full cubes that are not cut-out around the source velocity, switch to
`fullcube=True`.

Individual lines can be imaged by retrieving an instance of the `Spw` class
that abstracts properties specific to an Spw.

```python
dset = DataSet('CB68', setup=2, kind='joint')
spw = ALL_SPW_SETS[2][29]  # Setup 2, CS (5-4)
clean_line_target(dset, spw, weighting='natural', ext='test1')
clean_line_target(dset, spw, weighting=0, ext='test2')
```

Note the default parameters are `dirty=True` and `fullcube=True`, implicit in
the above. Please see the docstring via `clean_line_target?` for more
information. This function is also where most things "actually happen" in terms
of program logic, but the details of calculating most values occurs in the
attributes and methods of the classes `Target`, `Spw`, and `DataSet`, in
addition to some stand-alone utility functions for book-keeping purposes.
Attributes and methods are terms in Python for variables and functions bound to
classes.

To image all combinations of SPWs and uv-weightings, it is sufficient to run:

```python
dset = DataSet('CB68', setup=2, kind='joint')
run_pipeline(dset)
```

The default weightings to use are 'natural' and briggs robust 0.5 . A list
of weightings may be passed to the `run_pipeline` call using the `weightings`
keyword arguement for other combinations.

