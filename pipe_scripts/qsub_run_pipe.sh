#!/bin/bash
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
#  * That the number of batches is commensurate with the requested memory
#    and core count.


cd $PBS_O_WORKDIR  # directory where `qsub` was executed.
CASAPATH=/home/casa/packages/RHEL7/release/current
PATH=$CASAPATH/bin:$PATH
SCRIPTNAME=$PBS_O_WORKDIR/run_pipe.py
function run_casa {
    xvfb-run -d $CASAPATH/bin/casa --nogui --nologger -c "$1"
}

# The number of CASA instances or jobs/batches should be about half the number
# of CPUs requested (i.e., ~two CPUs per job). This shell environment variable
# will be retrieved in the Python script. Note that there must be an equal or
# greater number of chunks compared to jobs (with the minimum of one chunk per
# CASA job).
export NBATCHES=8


# Move the visibility data into shared memory for fast IO. Appropriate for MSs
# totalling less than 250 GB and on nodes with 500 GB memory.
export USING_SHARED_MEM=false
if [ $USING_SHARED_MEM = true ] ; then
    export SHM_DIR=/dev/shm/faust_pipeline
    DATA_DIR=/lustre/aoc/users/USER/FAUST/2018.1.01205.L/completed_SBs
    VIS_DATA=$DATA_DIR/CB68-Setup1-mosaic
    SHM_DATA=$SHM_DIR/CB68-Setup1
    if [[ ! -d $SHM_DATA ]] ; then
        cp -r $VIS_DATA $SHM_DATA
    fi
fi


# Run the CASA processes for:
#   * synchronous image preprocessing with "_preprocess()"
#   * asynchronous subsets for chunked imaging with "_run_subset($i)"
#   * synchronous image postprocessing with "_postprocess()"
run_casa "execfile('$SCRIPTNAME'); _preprocess()" >& casa_imaging_${PBS_JOBNAME}_preprocess.out
for ((i=0; i<$NBATCHES; i++))
do
    run_casa "execfile('$SCRIPTNAME'); _run_subset($i)" >& casa_imaging_${PBS_JOBNAME}_${i}.out &
    # CASA log files are indexed by the date and time to the nearest second.
    # Sleep for a couple seconds to ensure that logs do not clobber each other.
    sleep 2
done
wait
run_casa "execfile('$SCRIPTNAME'); _postprocess()" >& casa_imaging_${PBS_JOBNAME}_postprocess.out


# Clean up the temporary files stored in memory and links. This may not be
# strictly necessary but could cause problems if the files persist.
if [ $USING_SHARED_MEM = true ] ; then
    rm -rf $SHM_DATA
    for PROD_LINK_NAME in images moments plots
    do
        LINK_NAME=$SHM_DIR/$PROD_LINK_NAME
        [[ -L $LINK_NAME ]] && unlink $LINK_NAME
    done
fi


