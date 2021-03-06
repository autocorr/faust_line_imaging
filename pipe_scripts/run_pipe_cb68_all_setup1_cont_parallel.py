#!/usr/bin/env python
"""
For instructions on using this template script, please refer to the Cookbook
entry "Parallelized computation with multiple CASA processes" in the
documentation.

In this example:

    * "qsub_run_pipe.sh" requests 16 cores and spawns 8 CASA processes.
      Ensure that the "_get_config" function is called before the parallel
      run is started and the "_postprocess" function is called afterward.
    * This script is similar to `run_pipe_cb68_all_setup1.py` except that
      the continuum window chunks are run in parallel and distributed among
      the CASA instances. Thus, each CASA instance should process one or two
      narrow-band SPWs and several chunks of the continuum SPW.
"""

# FIXME Please change the path below such that is correct for the location of
# the "faust_imaging.py" script on your system.  Once the code is executed, all
# functions and variables from the pipeline will be in scope.  The "_"
# prefixing the names below helps ensure we don't accidentally shadow variables
# set in `faust_imaging.py`.
_SCRIPT_NAME = '../casa_scripts/faust_imaging.py'
if sys.version_info >= (3, 0, 0):
    with open(_SCRIPT_NAME) as f:
        code = compile(f.read(), _SCRIPT_NAME, 'exec')
        exec(code)
else:
    execfile(_SCRIPT_NAME)

if os.getenv('USING_SHARED_MEM') == 'true':
    DATA_DIR = os.getenv('SHM_DIR')  # mutating definition from `faust_imaging.py`
    os.chdir(DATA_DIR)
    if not os.path.exists('images'):
        os.symlink(IMAG_DIR, 'images')
        os.symlink(PLOT_DIR, 'plots')
        os.symlink(MOMA_DIR, 'moments')

# The number of batches should be defined in the torque shell script and
# should be either twice the number of chunks or half the number of CPUs.
_NBATCHES = int(os.getenv('NBATCHES', default=4))
_FIELD = 'CB68'
_RUN_EXT = 'clean'
_SETUP = 1
_SPW_SET = SPWS_BY_SETUP[_SETUP]  # defined in `faust_imaging.py`


def _get_cont_chunks():
    # There should be only one "cont" SPW per Setup.
    label = [s for s in _SPW_SET if 'cont' in s][0]
    full_config = ImageConfig.from_name(_FIELD, label)
    chunked_configs = full_config.duplicate_into_chunks()
    return chunked_configs


def _get_config():
    """
    Merge the configurations for the narrow-band SPWs to be processed in serial
    by chunk and the continuum SPW to be processed in parallel by chunk.

    Returns
    -------
    [ImageConfig]
    """
    cont_configs = _get_cont_chunks()
    # Create narrow-band SPW configs and merge with the continuum configs
    narrow_configs = [
            ImageConfig.from_name(_FIELD, label)
            for label in _SPW_SET if 'cont' not in label
    ]
    return narrow_configs + list(cont_configs)


def _preprocess():
    """
    Before starting the jobs, first create the image files required for
    determining the chunk starting frequencies (i.e., "_tinyimg.sumwt"). If
    this file already exists it will move on.
    """
    _get_config()


def _run_subset(batch_ix):
    """
    Run the pipeline for a subset of image configurations. SPWs are processed
    in stride based on the number of batches. For 13 SPws in Setup 1 and 8
    batches, ``batch_ix=0`` will process configurations:

        (0) DCOp,
        (7) SO,
        (14) 233.796 GHz continuum chunk2 (indexed from 0),
        (21) 233.796 GHz continuum chunk9,
        (..) etc.

    Parameters
    ----------
    batch_ix : int
        Batch index number.
    """
    all_configs = _get_config()
    # Run batch.
    nbatches = _NBATCHES
    batch_ix = int(batch_ix)
    assert batch_ix < nbatches
    log_post(':: Running batch index: {0}'.format(batch_ix))
    log_post('-- Batch: {0} / {1}'.format(batch_ix+1, nbatches))
    # Process each SPW or chunk. Starting at `batch_ix` step by every `_NBATCHES`.
    for config in all_configs[batch_ix::nbatches]:
        # Pipeline processes specific to the config may be included here. The
        # `.run_pipeline` method will automatically perform the post-
        # processing for narrow-band SPWs but not for the continuum chunks.
        config.run_pipeline(ext=_RUN_EXT)


def _postprocess():
    chunked_configs = _get_cont_chunks()
    chunked_configs.postprocess(ext=_RUN_EXT)
    make_all_qa_plots(_FIELD, ext=_RUN_EXT, overwrite=True)
    make_all_moment_maps(_FIELD, ext=_RUN_EXT, overwrite=True)


