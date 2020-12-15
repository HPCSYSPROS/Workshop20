#!/bin/bash
#SBATCH --partition=regular
#SBATCH --nodes=3
#SBATCH --time=01:00:00
#SBATCH --job-name=milc-small

build_lattice=true
. run_small.sh

