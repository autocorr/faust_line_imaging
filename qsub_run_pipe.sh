#!/bin/sh
#PBS -V
#PBS -l mem=500gb
#PBS -w /lustre/aoc/users/bsvoboda/faust/faust_alma/data
#PBS -m abe
#PBS -l nodes=1:ppn=16
#PBS -N "CB68_CS"

cd $PBS_O_WORKDIR

CASAPATH=/home/casa/packages/RHEL7/release/current
PATH=$CASAPATH/bin:$PATH

SCRIPTNAME=$PBS_O_WORKDIR/run_pipe.py
# The number of CASA instances or jobs/batches should be about half the number
# of CPUs requested (i.e., ~two CPUs per job). The shell environment variable
# will be retrieved in the script. Note that the number of batches must be
# equal to or fewer than the number of chunks defined in the script (minimum
# one chunk for one job/batch).
NBATCHES=8


# Start up the number of jobs asynchronously by appending "&"
for ((i=0; i<$NBATCHES; i++))
do
    xvfb-run -d $CASAPATH/bin/casa --nogui --nologger -c "execfile('$SCRIPTNAME'); _run_subset($i)" >& casa_imaging_${PBS_JOBNAME}_${i}.out &
done
wait


