#!/bin/sh
#PBS -V
#PBS -l mem=500gb
#PBS -w /lustre/aoc/users/bsvoboda/faust/faust_alma/data
#PBS -m abe
#PBS -l nodes=1:ppn=16
#PBS -N "run_pipe"

# The name of the job set with "-N" above (put into $PBS_JOBNAME) affects the
# name of the logs generated from CASA and the name of the process shown with
# `qstat`. While not strictly necessary, it is recommended to give this a
# unique label so separate processes do not overwrite one anothers log files.

cd $PBS_O_WORKDIR  # directory where `qsub` was executed.
CASAPATH=/home/casa/packages/RHEL7/release/current
PATH=$CASAPATH/bin:$PATH
SCRIPTNAME=$PBS_O_WORKDIR/run_pipe.py

# The number of CASA instances or jobs/batches should be about half the number
# of CPUs requested (i.e., ~two CPUs per job). This shell environment variable
# will be retrieved in the Python script. Note that there must be an equal or
# greater number of chunks compared to jobs (with the minimum of one chunk per
# CASA job).
export NBATCHES=8


# Before starting the jobs, first create the image files required for
# determining the chunk starting frequencies (i.e., "_tinyimg.sumwt"). If this
# file already exists, it will use the existing and move on. Note that this can
# be commented ouut for the "run one SPW per process" sorts of batch scripts.
casa --nogui --nologger -c "execfile('$SCRIPTNAME'); _get_config()" >& casa_imaging_${PBS_JOBNAME}_startup.out
# Start up the number of jobs asynchronously by appending "&"
for ((i=0; i<$NBATCHES; i++))
do
    xvfb-run -d $CASAPATH/bin/casa --nogui --nologger -c "execfile('$SCRIPTNAME'); _run_subset($i)" >& casa_imaging_${PBS_JOBNAME}_${i}.out &
done
wait


