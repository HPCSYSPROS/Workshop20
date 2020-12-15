################################################
#!/bin/bash
executable=su3_rhmd_hisq
frequency=performance
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
    /scratch/hpc_tools/read_imc_skx > $data_dir/imc.$HOST.$SLURM_JOBID
fi

perf stat -x' ' -o $data_dir/freq.${HOST}.${PMI_RANK}${MPIRUN_RANK}.${SLURM_JOBID} ./$executable

if [ "$MY_LOCAL_RANK" = "0" ]; then
/scratch/hpc_tools/read_imc_skx >> $data_dir/imc.$HOST.$SLURM_JOBID
fi
##################################################

