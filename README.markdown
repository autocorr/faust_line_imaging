FAUST Line Imaging README
=========================
This file accompanies `faust_imaging.py`, a framework for imaging ALMA FAUST
spectral line data. For viewing the file, I recommend using an editor that can
'fold' code at the function level. To run this script in CASA v5, execute in
the CASA prompt:

```python
execfile('/PATH/TO/faust_imaging.py')
```

Note that relative paths are also valid:

```python
execfile('../casa_scripts/faust_imaging.py')
```


Paths and directories
---------------------
This script currently requires an explicit directory structure for where to
retrieve the measurement sets and and where to write the images. Note the three
global paths:

```
ROOT_DIR -- absolute path to workspace
DAT_DIR -- directory where CASA should be run
IMG_DIR -- directory where output images will be written.
```

with files organized as a tree:

```
root/
  - data/   <- run CASA from here
    + CB68-Setup1
    + CB68-Setup2
    + ...
    + images/
      ~ CB68/
      ~ L1527/
      ~ ...
    + moments/
      ~ CB68/
      ~ L1527/
      ~ ...
```

Some tasks do not work with absolute paths, so please note that (as currently
written) CASA must be run from the directory specified by `DAT_DIR`. The MS files
should be present (by default) in the directories:

```bash
$DAT_DIR/<TARGET>-Setup1/
$DAT_DIR/<TARGET>-Setup2/
$DAT_DIR/<TARGET>-Setup2/
```

where `<TARGET>` is the FAUST target "field" name, e.g. "CB68" (names may be
found in the `ALL_TARGETS` dictionary, I retrieved these values from the
proposal, so they may be inconsistent for some targets). The paths above may
modified by editing the format string attribute `DataSet.ms_fmt`.


Making images
-------------
The `DataSet` class abstracts properties specific to a spectral setup. Please
see the docstring for more information. After execfile'ing the script, an
instance may created by:

```python
dset = DataSet('CB68', setup=2, kind='joint')
```

As currently written, dirty cubes must be made first to get estimates of the
RMS for each window and uv-weighting. To generate cut-out dirty cubes for
all lines with natural and robust=0 weighting, run:

```python
image_all_lines_inspect_dirty(dset)
```

If run with `mpicasa`, `tclean` will be run in parallel, otherwise it will run
normally. After these are created, auto-masking using auto-multithresh can be
used to create cleaned cubes by running:

```python
clean_all_lines_with_masking(dset, fullcube=False)
```

For full cubes that are not cut out around the line, switch to `fullcube=True`.

Individual lines can be imaged by retrieving an instance of the `Spw` class
that abstracts properties specific to an Spw.

```python
dset = DataSet('CB68', setup=2, kind='joint')
spw = ALL_SPW_SETS[2][29]  # Setup 2, CS (5-4)
clean_line_target(dset, spw, weighting='natural', ext='test1')
clean_line_target(dset, spw, weighting=0, ext='test2')
```

Note the default parameters are `dirty=True` and `fullcube=True`, implicit in
the above. Please see the docstring for `clean_line_target` for more
information. This is also where most things "actually happen" in terms of program
logic, but the details of calculating most values occurs in the attributes and
methods of the classes `Target`, `Spw`, and `DataSet`, in addition to some
stand-alone utility functions for book-keeping purposes. Attributes and methods
are terms in Python for variables and functions bound to classes.


