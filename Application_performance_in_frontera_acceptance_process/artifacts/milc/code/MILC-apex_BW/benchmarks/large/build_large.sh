#!/bin/bash
#SBATCH --partition=regular
#SBATCH --nodes=1728
#SBATCH --time=01:00:00
#SBATCH --job-name=milc-large

build_lattice=true
. run_large.sh

