#!/bin/bash
#PBS -N milc-large
#PBS -l nodes=1296:ppn=32:xe
#PBS -l walltime=02:55:00
#PBS -j oe
#

#------------------------------------------------------------------

# Specify ranks per node and threads per rank to use in this run
  ranks_per_node=32
  threads_per_rank=1

# Shorten the benchmark drastically by stopping after 10000 iterations
#  export PROF_NUM_ITERS=10000

# Disable Craypat-lite collecting performance data at run time
#  export PAT_RT_RECORD=0

# Disable pat_report running at end of job
  export PAT_RT_REPORT_METHOD=0

#------------------------------------------------------------------

# Set up environment same as when exe was built.
  module swap PrgEnv-cray PrgEnv-intel

# Memory management -- larger chunks
  export MALLOC_MMAP_MAX_=0
  export MALLOC_TRIM_THRESHOLD_=536870912

  export MPICH_ENV_DISPLAY=1
  export MPICH_VERSION_DISPLAY=1
 
  export MPICH_OPTIMIZED_MEMCPY=2
  export MPICH_USE_DMAPP_COLL=1
  export MPICH_SHARED_MEM_COLL_OPT=1
#  export MPICH_GNI_DMAPP_INTEROP=disabled
 
#  export PGT_ENABLE_PROCMASK=-1
 
#  export MPICH_NEMESIS_ASYNC_PROGRESS=SC
#  export MPICH_MAX_THREAD_SAFETY=multiple

  JOBDIR=`pwd`
  JOBDIR=${PBS_O_WORKDIR:-$JOBDIR}
  echo "JOBDIR is $JOBDIR"
  
  JOBID=$$
  JOBID=${PBS_JOBID:-$JOBID}
  JOB_ID="`echo $JOBID | sed -e 's/\..*//'`"
  echo "JOB_ID is $JOB_ID"

  JOBDIR=${JOBDIR}"/results_${JOBID}"
  mkdir -p ${JOBDIR}

  cd $JOBDIR
  ln -s ../su3_rhmd_hisq .
  ln -s ../rationals_m001m05m5.test1 .
  ln -s ../milc_in.sh .
  ln -s ../72x72x72x144.chklat .

# Handle different workload managers/batch systems
  NODES=1296
  WLM="INTERACTIVE"
  NN=99999
  NN=${PBS_NUM_NODES:-$NN}
  if [ $NN -ne 99999 ]; then
    WLM="PBS"
    NODES=$PBS_NUM_NODES
    CORES=$PBS_NUM_PPN
    run_type=aprun
  fi
  echo "WLM is $WLM"
  echo "NODES in allocation is $NODES"

  # load the run_milc() function
  . milc_in.sh

  #----- START USER DEFINED PARAMETERS -----#
  #build_lattice=true        # Set only if you want to create a lattice file
                            # Comment out otherwise
                            # You really only want to do this once

  total_nodes=$NODES             # match above allocation request!

  # this is for Interlagos
  let CORES="$CORES / 2"
  hyperthreads=2
  numa_per_node=4
  module load craype-hugepages2M
  
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

  if [ $threads_per_node -gt $cores_per_node ]; then
# Compute the required j factor
    let threads_per_core="$threads_per_node/$cores_per_node"
    let check="$threads_per_core*$cores_per_node"
    if [ $check -lt $threads_per_node ]; then
      let threads_per_core="$threads_per_core + 1"
    fi
    S="$S -j $threads_per_core"
  else
# One thread per core
    S="$S -j 1"
  fi

  PPN=$ranks_per_node
  let WIDTH="$NODES * $PPN"

  echo "NODES to use = $NODES"
  echo "PPN (ranks per node) to use = $PPN"

  # sanity check for parameters
#  if [ $(( $N*$threads_per_rank )) -gt $(( total_nodes*$cpus_per_node )) ]; then
#    echo "$0: ERROR, number of threads exceeds total concurrency, aborting"
#    exit;
#  else
    echo "$0: Using $N MPI rank(s), $S rank(s) per NUMA domain, $threads_per_rank thread(s) per rank"
#  fi

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
  run_milc > milc_${dimension}_on_${total_nodes}_${ranks_per_node}_${threads_per_rank}.$JOB_ID
  echo -n "$0: End time: "; date
  echo "UNIX_TIME_END=`date +%s`"

  exit
