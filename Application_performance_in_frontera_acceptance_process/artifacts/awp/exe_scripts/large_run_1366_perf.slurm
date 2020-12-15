#!/bin/bash
#SBATCH -J awp_large    
#SBATCH -o awp_large.%j.out
#SBATCH -e awp_large.%j.out
#SBATCH -p test2
#SBATCH -N 1366 
#SBATCH -n 65536            
#SBATCH -t 02:00:00     
#SBATCH -A A-ccsc       


mkdir -p large/$SLURM_JOBID
cd large/$SLURM_JOBID

bash ../../pre-run
source /scratch/hpc_tools/spp_ibrun/sourceme.sh
#source /scratch/hpc_tools/spp_ibrun/antisourceme.sh
cp ../../IN3D_large IN3D

export executable=pmcl3d
export frequency=performance
export data_dir=${executable}_${frequency}_${SLURM_JOBID}
if [ "$ANTIAFFINITY" -eq "1" ]; then
    export data_dir=anti_$data_dir
fi
export I_MPI_STATS=0
ln -s ../../../../../src-v1.1.2/pmcl3d

ls -al

start_time=$(date +%s)
ibrun tacc_affinity instrumented_run.sh ./pmcl3d # tacc_affinity
#ibrun tacc_affinity ./pmcl3d # tacc_affinity
#ibrun ../../spp_stats.sh # no affinity set
end_time=$(date +%s)
echo Completed in $(( end_time - start_time )) seconds

cd ${data_dir}
/scratch/hpc_tools/average_freq.py
/scratch/hpc_tools/average_imc.py
