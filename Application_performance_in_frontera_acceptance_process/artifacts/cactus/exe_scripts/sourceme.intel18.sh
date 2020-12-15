### source /scratch/projects/compilers/sourceme.intel18.sh

compdir="/scratch/projects/compilers"
topdir=$compdir/intel18
ver="2018.0.033"

module unload impi
module unload intel
module unload xalt    # xalt expects a certain linker; we can't guarantee this here

source $topdir/parallel_studio_xe_$ver/psxevars.sh

export TACC_MPI_GETMODE=impi_hydra
export MPICH_HOME=$I_MPI_ROOT

