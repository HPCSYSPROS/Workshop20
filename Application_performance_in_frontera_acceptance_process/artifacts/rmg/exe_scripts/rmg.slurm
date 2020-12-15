#!/bin/bash
#SBATCH -J large_rmg
#SBATCH -o rmg.%J
#SBATCH -e rmg.%J
#SBATCH -N 1536
#SBATCH -n 18432
#SBATCH -t 02:00:00
#SBATCH -p test2
#SBATCH -A A-ccsc

ml gcc
ml fftw3
ml boost

export base=/scratch/03078/cproctor/apps/spp/rmg/for_damon/code
export LD_LIBRARY_PATH=${base}/openbabel/install/gcc/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${base}/plplot/install/gcc/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${base}/rmg/install/gcc/lib:${LD_LIBRARY_PATH}
export PATH=${base}/rmg/install/gcc/bin:${PATH}

export OMP_NUM_THREADS=4
ibrun rmg-cpu in.dnv4096_3456xe_large > job.${SLURM_JOB_ID}.out
