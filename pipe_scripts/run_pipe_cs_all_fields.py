#!/usr/bin/env python
"""
For instructions on using this template script, please refer to the Cookbook
entry "Parallelized computation with multiple CASA processes" in the
documentation.

In this example:

    * "qsub_run_pipe.sh" requests 16 cores and spawns 8 CASA processes
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
_LABEL = '244.936GHz_CS'


def _run_subset(batch_ix):
    """
    Run the pipeline for a given SPW for the full set of fields.

    Parameters
    ----------
    batch_ix : int
        Batch index number.
    """
    # Configure parameters here.
    ext = 'clean'
    # Run batch.
    nbatches = _NBATCHES
    batch_ix = int(batch_ix)
    assert batch_ix < nbatches
    log_post(':: Running batch index: {0}'.format(batch_ix))
    log_post('-- Batch: {0} / {1}'.format(batch_ix+1, nbatches))
    # Process each SPW. Starting at `batch_ix` step by every `_NBATCHES`.
    for field in ALL_FIELD_NAMES[batch_ix::nbatches]:
        # Pipeline processes specific to SPW may be included here. The
        # `.run_pipeline` method will automatically perform the post-
        # processing.
        config = ImageConfig.from_name(field, _LABEL)
        config.run_pipeline(ext=ext)


