#!/bin/bash

#SBATCH -J namd
#SBATCH -N 1600
#SBATCH -n 76800
#SBATCH -t 00:30:00
#SBATCH -p test2
#SBATCH -A A-ccsc

module load intel/16.0.3

ibrun /work/00410/huang/namd/build/2.12_skx/NAMD_2.12_Source/Linux-KNL-icc/namd2 chromat100-bench.namd &> log


