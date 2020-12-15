#!/bin/bash
#SBATCH -A A-ccsc
#SBATCH -t 01:00:00
#SBATCH -o milc-large.%j.out
#SBATCH -p test2
#SBATCH -N 1296
#SBATCH --ntasks-per-node 48

#------------------------------------------------------------------
ml rm xalt
# Specify ranks per node and threads per rank to use in this run

ranks_per_node=48
threads_per_rank=1

# Specify OpenMP thread number, affinity, etc.
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
# Intel default is passive; use only if hyperthreading?
export OMP_WAIT_POLICY=passive
#  export OMP_WAIT_POLICY=active
#  export KMP_AFFINITY=disabled

NODES=$SLURM_JOB_NUM_NODES
CORES=$SLURM_CPUS_ON_NODE
run_type=ibrun
echo "$SLURM_JOB_NODELIST" > node_list.$SLURM_JOBID
echo "NODES in allocation is $NODES"

# load the run_milc() function
. milc_in.sh

#----- START USER DEFINED PARAMETERS -----#
#build_lattice=true        # Set only if you want to create a lattice file
# Comment out otherwise
# You really only want to do this once

total_nodes=$NODES             # match above allocation request!
KNL="`echo $LOADEDMODULES | grep 'knl'`"
# For Knights Landing
hyperthreads=4
numa_per_node=1
  
let CORES="$CORES / $hyperthreads"
cores_per_node=$CORES
let cpus_per_node="$cores_per_node*$hyperthreads"
echo "cores_per_node on each node is $cores_per_node"
echo "cpus_per_node on each node is $cpus_per_node"

# Parse command line parameters
debug="false"
while [ "$#" -ge 1 ]; do
    case "$1" in
	"--build" )           build_lattice=true; shift 1;;
	"--debug" )           debug="true"; shift 1;;
	"--ranks" | "-r" )    ranks_per_node=$2; shift 2;;
	"--threads" | "-t" )  threads_per_rank=$2; shift 2;;
	*)                    break;;
    esac
done
#----- END USER DEFINED PARAMETERS -----#

export OMP_NUM_THREADS=$threads_per_rank

# Set problem size
nx=72
ny=72
nz=72
nt=144

N=$(( $total_nodes*$ranks_per_node ))  # total ranks
S=$(( $ranks_per_node/$numa_per_node ))                 # ranks per socket
threads_per_node=$(( $ranks_per_node*$threads_per_rank ))

#S="ibrun numactl --preferred=1"
source /scratch/hpc_tools/spp_ibrun/sourceme.sh
S="ibrun tacc_affinity"

# Lattice size
F1=$nx
F2=$ny
F3=$nz
F4=$nt

# sanity check for parameters
if [ $(( $N*$threads_per_rank )) -gt $(( total_nodes*$cpus_per_node )) ]; then
    echo "$0: ERROR, number of threads exceeds total concurrency, aborting"
    exit;
else
    echo "$0: Using $N MPI rank(s), $S rank(s) per NUMA domain, $threads_per_rank thread(s) per rank"
fi

dimension=${nx}x${ny}x${nz}x${nt}
checklat=$dimension.chklat
if [ "$build_lattice" == "true" ]; then
    warms=40
    trajecs=2
    save_cmd="save_serial $checklat"
    reload_cmd="continue"
    echo "$0: Building lattice $dimension on $total_nodes nodes"
else
    warms=0
    trajecs=2
    save_cmd="forget"
    reload_cmd="reload_parallel $checklat"
    echo "$0: Running lattice $dimension on $total_nodes nodes"
  fi

# Run MILC
echo -n "$0: Begin time: "; date
echo "UNIX_TIME_START=`date +%s`"
# Write output to a file that we can see during execution
#  run_milc
run_milc > milc_${dimension}_on_${total_nodes}_${ranks_per_node}_${threads_per_rank}.$SLURM_JOBID
echo -n "$0: End time: "; date
echo "UNIX_TIME_END=`date +%s`"
