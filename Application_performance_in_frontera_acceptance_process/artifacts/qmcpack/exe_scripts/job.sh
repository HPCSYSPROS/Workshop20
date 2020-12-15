#!/bin/bash
#SBATCH -J myjob_3200 # Job name
#SBATCH -o myjob.o%j # Name of stdout output file
#SBATCH -e myjob.e%j # Name of stderr error file
#SBATCH -p test2 # Queue name
#SBATCH -N 1200 # Total # of nodes (now required)
#SBATCH -n 1200 # Total # of mpi tasks
#SBATCH -t 8:00:00 # Run time (hh:mm:ss)
#SBATCH --mail-user=myname@myschool.edu
#SBATCH --mail-type=all # Send email at begin and end of job
#SBATCH -A A-ccsc # Allocation name (req'd if more than 1)
# Other commands must follow all #SBATCH directives...
pwd
date
# Launch MPI application...
# Use ibrun instead of mpirun or mpiexec


ml load boost
ml load fftw3
ml load hdf5
export MKL_DEBUG_CPU_TYPE=7
export OMP_NUM_THREADS=48
#export I_MPI_STATS=20
ibrun wrapper.sh
