#!/usr/bin/env python
"""
For instructions on using this template script, please refer to the Cookbook
entry "Parallelized computation with multiple CASA processes" in the
documentation.

In this example:

    * "qsub_run_pipe.sh" requests 16 cores and spawns 8 CASA processes
    * The CB68 Setup 1 continuum cube (12m-array SPW ID number 49) is sharded
      into 40 chunks (the default heuristic) and distributed across the 8 CASA
      processes (5 chunks per process).
    * The emission in the continuum cube tends to be compact and over a smaller
      number of channels (channels ~1.2 km/s versus ~0.15 km/s in the narrow),
      so the RMS can be determined reliably per chunk.
    * Post-processing must be performed separately after the jobs finish
      running by calling the `_postprocess` function.
"""

# FIXME Please change the path below such that is correct for the location of
# the "faust_imaging.py" script on your system. The pipeline is intended to be
# run under CASA v5 and Python v2.7; once `execfile` is run, all functions and
# variables from the pipeline will be in scope.  The "_" prefixing the names
# below helps ensure we don't accidentally shadow variables set in
# `faust_imaging.py`.
execfile('../casa_scripts/faust_imaging.py')

# The number of batches should be defined in the torque shell script and
# should be either twice the number of chunks or half the number of CPUs.
_NBATCHES = int(os.getenv('NBATCHES', default=4))
_FIELD = 'CB68'
_RUN_EXT = 'clean'


def _get_config():
    # Edit values for desired target setup here.
    label = '233.796GHz_cont'
    # Configure global settings used for all chunks here.
    full_config = ImageConfig.from_name(_FIELD, label)
    chunked_configs = full_config.duplicate_into_chunks()
    return full_config, chunked_configs


def _preprocess():
    """
    Before starting the jobs, first create the image files required for
    determining the chunk starting frequencies (i.e., "_tinyimg.sumwt"). If
    this file already exists it will move on.
    """
    _get_config()


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
    full_config, chunked_configs = _get_config()
    # Process each chunk. Starting at `batch_ix` step by every `_NBATCHES`.
    for config in chunked_configs[batch_ix::nbatches]:
        # Pipeline processes specific to a chunk may included here.
        config.run_pipeline(ext=_RUN_EXT)


def _postprocess():
    full_config, chunked_configs = _get_config()
    chunked_configs.postprocess(ext=_RUN_EXT)
    imagebase = full_config.get_imagebase(ext=_RUN_EXT)
    imagepath = '{0}.image'.format(imagebase)
    make_qa_plots_from_image(imagepath, overwrite=True)
    make_moments_from_image(imagepath, overwrite=True)


