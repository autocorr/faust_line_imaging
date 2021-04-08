#!/bin/sh
#PBS -V
#PBS -l mem=500gb
#PBS -w /lustre/aoc/users/bsvoboda/faust/faust_alma/data
#PBS -m abe
#PBS -l nodes=1:ppn=16
#PBS -N "run_pipe"

# The name of the job set with "-N" above (put into $PBS_JOBNAME) will only
# affect the name of the logs generated from CASA and the name of the process
# shown with `qstat`. It may be helpful to label it by source and transition,
# but it is not strictly necessary.

cd $PBS_O_WORKDIR
CASAPATH=/home/casa/packages/RHEL7/release/current
PATH=$CASAPATH/bin:$PATH
SCRIPTNAME=$PBS_O_WORKDIR/run_pipe.py

# The number of CASA instances or jobs/batches should be about half the number
# of CPUs requested (i.e., ~two CPUs per job). This shell environment variable
# will be retrieved in the Python script. Note that there must be an equal or
# greater number of chunks compared to jobs (with the minimum of one chunk per
# CASA job).
NBATCHES=8


# Start up the number of jobs asynchronously by appending "&"
for ((i=0; i<$NBATCHES; i++))
do
    xvfb-run -d $CASAPATH/bin/casa --nogui --nologger -c "execfile('$SCRIPTNAME'); _run_subset($i)" >& casa_imaging_${PBS_JOBNAME}_${i}.out &
done
wait


