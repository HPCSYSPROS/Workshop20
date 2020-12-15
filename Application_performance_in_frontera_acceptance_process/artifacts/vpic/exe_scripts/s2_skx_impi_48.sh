#!/bin/bash
#SBATCH -J vpic_73728_48     # job name
#SBATCH -o vpic_73728.%j     # output and error file name (%j expands to jobID)
#SBATCH -N 1536              # number of nodes requested
#SBATCH --ntasks-per-node 48 # total number of mpi tasks requested
#SBATCH -p test2             # queue (partition) -- normal, development, etc.
#SBATCH -t 03:00:00          # run time (hh:mm:ss)
#SBATCH -A A-ccsc

#SBATCH --mail-user=dmcdougall@tacc.utexas.edu
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

export vpicexe=/scratch/02664/damon/spp/vpic/skx/73728/test_73728.s2_skx_intel-17.0.4_impi-17.0.3_core_avx512
export NP=73728
export NPPNODE=48

module purge
module load xalt/1.7.7
module load intel/17.0.4
module load impi/17.0.3

mkdir ${SLURM_JOBID}
cd ${SLURM_JOBID}

date
export PATH=/scratch/hpc_tools/spp_ibrun:$PATH
time ibrun tacc_affinity instrumented_run.sh ${vpicexe} -tpp=1

cp ../s2_skx_impi_48.sh .
