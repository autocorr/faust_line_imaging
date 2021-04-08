#!/usr/bin/env python
"""
For instructions on using this template script, please refer to the Cookbook
entry "Parallelized computation with multiple CASA processes" in the
documentation.

In this example:

    * "qsub_run_pipe.sh" requests 16 cores and spawns 8 CASA processes
    * There are 13 SPWs for Setup 1, so each CASA instance will process
      1-2 SPWs. Because the continuum takes longer to process than the
      others, and some jobs may get SPWs that are mostly noise, it is likely
      that the time it takes each job to finish will vary considerably.
    * Unlike the "parallelized single SPW" examples, this script does
      perform the post-processing for each SPW.
    * This is currently untested! Hopefully mixing the post-processing
      with the imaging among different instances of CASA does not cause
      memory issues.
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


def _run_subset(batch_ix):
    """
    Run the pipeline for a subset of SPWs. SPWs are processed in stride
    based on the number of batches. For 13 SPws in Setup 1 and 8 batches,
    ``batch_ix=0`` would process SPWs (0) DCOp and (7) SO.

    Parameters
    ----------
    batch_ix : int
        Batch index number.
    """
    # Configure parameters here.
    field = 'CB68'
    spw_set = SPW_S1  # defined in `faust_imaging.py`
    ext = 'clean'
    # Run batch.
    nbatches = _NBATCHES
    batch_ix = int(batch_ix)
    assert batch_ix < nbatches
    log_post(':: Running batch index: {0}'.format(batch_ix))
    log_post('-- Batch: {0} / {1}'.format(batch_ix+1, nbatches))
    spw_labels = [s for s in spw_set.keys()]
    # Process each SPW. Starting at `batch_ix` step by every `_NBATCHES`.
    for label in spw_labels[batch_ix::nbatches]:
        # Pipeline processes specific to SPW may be included here. The
        # `.run_pipeline` method will automatically perform the post-
        # processing.
        config = ImageConfig.from_name(field, label)
        config.run_pipeline(ext=ext)


