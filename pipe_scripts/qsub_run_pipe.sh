#!/bin/sh
#PBS -V
#PBS -l mem=500gb
#PBS -w /lustre/aoc/users/bsvoboda/faust/faust_alma/data
#PBS -m abe
#PBS -l nodes=1:ppn=16
#PBS -N "run_pipe"
#       ^ The above is the job name (put into $PBS_JOBNAME by Torque) and
#         also used by this script as a suffix for log files containing the
#         STDOUT and STDERR of each CASA process. It is recommend to give the
#         job a unique name so that concurrant jobs do not clobber one
#         another's log files.

# NOTE In editing this script, please take care that:
#  * The PBS "-N" directive is given a unique name.
#  * SCRIPTNAME points to the correct corresponding Python script.
#  * NBATCHES is equal to or fewer than the number of chunks to be processed.
#  * The correct field name is used in the path to the MS files if copying the
#    data into RAM using /dev/shm.
#  * The lines containing "_get_config()" and "_postprocess()" are removed
#    or commented out if the Python script does not require them (such as
#    the "run one SPW in serial per CASA instance" sorts of batch jobs).
#  * That the number of batches is commensurate with the requested memory
#    and core count.


cd $PBS_O_WORKDIR  # directory where `qsub` was executed.
CASAPATH=/home/casa/packages/RHEL7/release/current
PATH=$CASAPATH/bin:$PATH
SCRIPTNAME=$PBS_O_WORKDIR/run_pipe.py
alias RUN_CASA="xvfb-run -d $CASAPATH/bin/casa --nogui --nologger -c"

# The number of CASA instances or jobs/batches should be about half the number
# of CPUs requested (i.e., ~two CPUs per job). This shell environment variable
# will be retrieved in the Python script. Note that there must be an equal or
# greater number of chunks compared to jobs (with the minimum of one chunk per
# CASA job).
export NBATCHES=8


# Move the visibility data into shared memory for fast IO. Appropriate for MSs
# totalling less than 250 GB and on nodes with 500 GB memory.
#DATA_DIR=/lustre/aoc/users/cchandle/FAUST/2018.1.01205.L/completed_SBs
#cp -r $DATA_DIR/CB68-Setup1-mosaic /dev/shm/CB68-Setup1
#export USING_SHM=True

# Before starting the jobs, first create the image files required for
# determining the chunk starting frequencies (i.e., "_tinyimg.sumwt"). If this
# file already exists, it will use the existing and move on.
RUN_CASA "execfile('$SCRIPTNAME'); _get_config()" >& casa_imaging_${PBS_JOBNAME}_startup.out
# Start up the number of jobs asynchronously by appending "&"
for ((i=0; i<$NBATCHES; i++))
do
    RUN_CASA "execfile('$SCRIPTNAME'); _run_subset($i)" >& casa_imaging_${PBS_JOBNAME}_${i}.out &
done
wait
RUN_CASA "execfile('$SCRIPTNAME'); _postprocess()" >& casa_imaging_${PBS_JOBNAME}_postprocess.out


