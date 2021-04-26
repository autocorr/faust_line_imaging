#!/bin/sh

# Example name prefix. $JOBNAME affects the names of the terminal STDOUT log
# files generated by each CASA process. This should be set to a unique name
# in order that log files of different processes don't overwrite each other.
# Please ensure that the script name matches the intended Python target.
JOBNAME="pipe_log"
SCRIPTNAME=run_pipe.py
alias RUN_CASA="casa --nogui --nologger -c"

# The number of CASA instances or jobs/batches should be about half the number
# of CPUs requested (i.e., ~two CPUs per job). The shell environment variable
# will be retrieved in the script. Note that there must be more chunks than
# batches (or an equal number, with the minimum of one chunk per CASA job).
export NBATCHES=8


# Before starting the jobs, first create the image files required for
# determining the chunk starting frequencies. If this file already exists, it
# will use the existing and move on.
# If the script does not require them, simply comment out or remove the startup
# and postprocessing lines.
RUN_CASA "execfile('$SCRIPTNAME'); _get_config()" >& casa_imaging_${JOBNAME}_startup.out
# Start up the number of CASA instances asynchronously by appending "&"
for ((i=0; i<$NBATCHES; i++))
do
    RUN_CASA "execfile('$SCRIPTNAME'); _run_subset($i)" >& casa_imaging_${JOBNAME}_${i}.out &
    # CASA log files are indexed by the date and time to the nearest second.
    # Sleep for a couple seconds to ensure that logs do not clobber each other.
    sleep 2
done
wait
RUN_CASA "execfile('$SCRIPTNAME'); _postprocess()" >& casa_imaging_${JOBNAME}_postprocess.out


