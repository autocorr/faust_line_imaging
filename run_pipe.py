#!/usr/bin/env python
"""
Running this script
-------------------
This script is meant to be run from the command line as:

    $ casa --nogui -c "execfile('run_pipe.py'); _run_subset(0)"

or from a torque shell script to be run with qsub as:

    xvfb-run -d casa --nogui -c "execfile('run_pipe.py'); _run_subset(0)" &

Lastly, it can be executed from a CASA session with `execfile`.

    $ casa
    > execfile('run_pipe.py')
    > full_config, chunked_configs = _get_config()
    > problematic_one = chunked_configs[42]
    > problematic_one.clean_line(ext='clean', interactive=True, restart=True)

To perform the postprocessing, run:

    $ casa
    > execfile('run_pipe.py')
    > _postprocess()

The post-processing cannot take advantage of parallelization, so it can be run
interactively from CASA after the previous jobs have finished.

How to configure
----------------
To modify this template, set the run properties using the global variables
`_NBATCHES` and `_RUN_SUFFIX`. Note that this script is to be run from the same
directory as defined by `PROD_DIR` in `faust_imaging.py`.

Global source and image configuration settings (field, line, perhaps the RMS or
auto-multithresh parameters, etc.) may be set in :func:`_get_config`.  Pipeline
procedures specific to a chunk may be added in :func:`_run_subset`.
"""

# NOTE To "import" the functions from `faust_imaging.py` under CASA v5,
# it must be first evaluated with `execfile` before running this script.
execfile('../casa_scripts/faust_imaging.py')

# The "_" are to help ensure we don't accidentally shadow variables set in
# `faust_imaging.py` that are now in scope.
# The number of batches should be defined in the torque shell script and
# should be either twice the number of chunks or half the number of CPUs.
_NBATCHES = int(os.getenv('NBATCHES', default=4))
_RUN_SUFFIX = 'clean'


def _get_config():
    # Edit values for desired target setup here.
    field = 'CB68'
    label = '244.936GHz_CS'
    nchunks = 4
    # Configure global settings used for all chunks here.
    full_config = ImageConfig.from_name(field, label)
    chunked_configs = full_config.duplicate_into_chunks(nchunks=nchunks)
    return full_config, chunked_configs


def _run_subset(batch_ix):
    """
    Run the pipeline for a subset of chunks. Chunks are processed in stride,
    e.g.: 0, 10, 20... or 1, 11, 21...

    Parameters
    ----------
    batch_ix : int
        Batch index number. For 100 chunks and 10 batches, ``batch_ix=0``
        would process chunks 0, 10, ..., 90.
    """
    nbatches = _NBATCHES
    batch_ix = int(batch_ix)
    assert batch_ix < nbatches
    log_post(':: Running batch index: {0}'.format(batch_ix))
    log_post('-- Batch: {0} / {1}'.format(batch_ix+1, nbatches))
    full_conifg, chunked_configs = _get_config()
    # Process each chunk. Starting at `batch_ix` step by every `_NBATCHES`.
    for config in chunked_configs[batch_ix::nbatches]:
        # Pipeline processes specific to a chunk may included here.
        config.run_pipeline(ext=_RUN_SUFFIX)


def _postprocess():
    full_config, chunked_configs = _get_config()
    chunked_configs.postprocess(ext=_RUN_SUFFIX)


