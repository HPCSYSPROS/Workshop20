#! /bin/bash
#SBATCH -o zelmani_skx.o%j
#SBATCH -e zelmani_skx.e%j
#SBATCH -p skx-normal
#SBATCH -t 01:00:00
#SBATCH -N 100
#SBATCH --tasks-per-node 2
#SBATCH -A A-ccsc
##SBATCH -A TG-PHY100033 
#SBATCH -J zelmanibench_skx

source ./sourceme.intel18.sh
# restore stuff that was purged
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/apps/intel17/hdf5/1.8.16/x86_64/lib

# compute cores and threads and stuff
export N=${SLURM_NNODES}
export tpn=${SLURM_NTASKS_PER_NODE}
export n=$(( tpn * N ))
#### Note: 
export t=$(( 48 / tpn ))
export OMP_NUM_THREADS=${t}
export OMP_PROC_BIND=true

##
## Flop counting switch
##
export FLOPCOUNT=0

##
## run, tee'ing all output to a timing file
##
if [ ${FLOPCOUNT} -eq 1 ] ; then
  export timingfile=timing-skx-fl-${N}-${tpn}.txt
else
  export timingfile=timing-skx-${N}-${tpn}.txt
fi

echo "Job: ${SLURM_JOBID}" | tee ${timingfile}
echo "Run ${N} nodes, tnp=${tpn}, ${t} way thread parallel" | tee -a ${timingfile}

echo "Starting:"
export CACTUS_STARTTIME=$(date +%s)
echo "start at ${CACTUS_STARTTIME}"
if [ ${FLOPCOUNT} -eq 1 ] ; then
  ibrun ./zelmani.sh ./cactus_zelmani s15WH07_R0001.par
else
  ibrun ./cactus_zelmani s15WH07_R0001.par
fi
echo "Stopping:"
export CACTUS_STOPTIME=$(date +%s)
echo "start at ${CACTUS_STOPTIME}"
echo "RunTime: $(( CACTUS_STOPTIME - CACTUS_STARTTIME )) s" | tee -a ${timingfile}
