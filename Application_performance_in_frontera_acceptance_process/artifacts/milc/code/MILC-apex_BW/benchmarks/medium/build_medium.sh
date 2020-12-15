#!/bin/bash
#SBATCH --partition=regular
#SBATCH --nodes=108
#SBATCH --time=01:00:00
#SBATCH --job-name=milc-medium

build_lattice=true
. run_medium.sh

