#!/bin/bash
pwd
date
export MKL_DEBUG_CPU_TYPE=7
export HOST=`hostname -s`
export MYRANK=$PMI_RANK

ml load boost
ml load fftw3
ml load hdf5

####################################################################
# Modify these to control the run setup
# Note that these values must allow the job to fit in the
#   reservation requested by the #PBS -l nodes=<num_nodes>:ppn=<ppn>
#   line above
# qmcexe - name of the executable
# NP - total number of MPI processes (aprun -n value)
# NPPNODE - number of MPI processes per node (aprun -N value)
# OMP_NUM_THREADS - number of OpenMP threads per MPI process (aprun -d value)
# ------------------------------------------------------------------
export qmcexe=qmcpack
export NP=1200
export NPPNODE=1
####################################################################

####################################################################
# Modify these to control the duration and accuracy of the run
# DMCBLOCKS - number of DMC blocks to run
# DMC_SAMPLES_PER_THREAD - number of DMC samples (walkers) run by each thread;
#       this needs to be high enough to saturate each thread computationally
# ------------------------------------------------------------------
export DMCBLOCKS=30
export DMC_SAMPLES_PER_THREAD=16
####################################################################

####################################################################
# More parameters set based on the above values; these shouldn't be modified 
# unless you know what you're doing
# VMCBLOCKS - number of VMC blocks to run
# VMCSTEPS - number of steps per VMC block
# VMCWALKERS - number of VMC walkers per thread
# SAMPLES - total number of DMC samples
####################################################################
export VMCBLOCKS=$DMC_SAMPLES_PER_THREAD
export VMCSTEPS=10
export VMCWALKERS=$OMP_NUM_THREADS
#let SAMPLES=$DMC_SAMPLES_PER_THREAD*$NP*$OMP_NUM_THREADS
let SAMPLES=2560000
####################################################################
export WFSXML=gr4x4x1.wfs.xml
export title=gr4x4x1.p${NP}x${OMP_NUM_THREADS}

export HUGETLB_MORECORE=yes
export COUNTER1=GET_TIME_OF_DAY

# create the .xml input file based on qmc.xml
cat qmc.xml\
| sed s/TITLE/$title/ \
| sed s/WFSXML/$WFSXML/ \
| sed s/VMCBLOCKS/$VMCBLOCKS/ \
| sed s/VMCSTEPS/$VMCSTEPS/ \
| sed s/WALKERS/$VMCWALKERS/ \
| sed s/SAMPLES/$SAMPLES/ \
| sed s/DMCBLOCKS/$DMCBLOCKS/ \
| sed s/SEED/13/ \
> ${title}.xml

# print some input parameters
date
echo "VMC blocks: " ${VMCBLOCKS}
echo "VMC walkers per MPI rank: " ${VMCWALKERS}
echo "DMC blocks: " ${DMCBLOCKS}
echo "DMC samples per thread: " ${DMC_SAMPLES_PER_THREAD}
echo "Total DMC samples: " ${SAMPLES}
echo ibrun -np ${NP} /scratch/03024/chenk/spp_frontiers/qmcpack_dev_20170512-1310/qmcpack/build/bin/qmcpack ${title}.xml

##########################################################################################################################################


executable=qmcpack

frequency=22
export I_MPI_STATS=0

data_dir=${executable}_${frequency}_${SLURM_JOBID}

if [ "$ANTIAFFINITY" -eq "1" ]; then
    data_dir=anti_$data_dir
fi

mkdir -p $data_dir

export I_MPI_STATS_FILE=$data_dir/mpi_stats.txt

export HOST=`hostname -s`
MY_TASKS_PER_NODE=$(( $PMI_SIZE$MPIRUN_SIZE / $SLURM_NNODES ))
MY_LOCAL_RANK=$(( $PMI_RANK$MPIRUN_RANK % $MY_TASKS_PER_NODE ))

if [ "$MY_LOCAL_RANK" = "0" ]; then
   if [ "$frequency" != "performance" ]; then
      /scratch/hpc_tools/TACC_HWP_set --freq $frequency
   fi
   unbuffer /scratch/hpc_tools/read_imc_skx > $data_dir/imc.$HOST.$SLURM_JOBID
fi

perf stat -x' ' -o $data_dir/freq.${HOST}.${PMI_RANK}${MPIRUN_RANK}.${SLURM_JOBID} /scratch/03024/chenk/spp_frontiers/qmcpack_dev_20170512-1310_amdlibm/skylake/qmcpack_-xCOMMON-AVX512/qmcpack/build//bin/${executable} ${title}.xml


if [ "$MY_LOCAL_RANK" = "0" ]; then
   unbuffer /scratch/hpc_tools/read_imc_skx >> $data_dir/imc.$HOST.$SLURM_JOBID
fi


#########################################################################################################################################


#perf stat -x' ' -o perf.out.${HOST}.${MYRANK}.${SLURM_JOBID} -e r01c7 -e r04c7 -e r10c7 -e r40c7 -e r02c7 -e r08c7 -e r20c7 -e r80c7 /scratch/03024/chenk/spp_frontiers/qmcpack_dev_20170512-1310_amdlibm/skylake/qmcpack_-xCOMMON-AVX512/qmcpack/build/bin/${qmcexe} ${titlf}.xml


