#!/bin/bash
#SBATCH --partition=regular
#SBATCH --nodes=108
#SBATCH --time=01:00:00
#SBATCH --job-name=milc-medium

  # load the run_milc() function
  . milc_in.sh

  #----- START USER DEFINED PARAMETERS -----#
  #build_lattice=true        # Set only if you want to create a lattice file
                            # Comment out otherwise
                            # You really only want to do this once

  total_nodes=108             # match above allocation request!
  cores_per_node=24
  numa_per_node=2
  threads_per_rank=6
  run_type=srun            # comment out to debug script

  # Parse command line parameters
  debug="false"
  while [ "$#" -ge 1 ] ; do
    case "$1" in
      "--build" )           build_lattice=true; shift 1;;
      "--debug" )           debug="true"; shift 1;;
      "--cores" | "-c" )    cores_per_node=$2; shift 2;;
      "--threads" | "-t" )  threads_per_rank=$2; shift 2;;
      *)                    break;;
    esac
  done

  #----- END USER DEFINED PARAMETERS -----#

  # Set problem size
  nx=36
  ny=36
  nz=36
  nt=72

  N=$(( $total_nodes*$cores_per_node/$threads_per_rank ))  # total ranks
  S=$(( $N/$total_nodes/$numa_per_node ))                 # ranks per socket

  # sanity check for parameters
  if [ $(( $N*$threads_per_rank )) -gt $(( total_nodes*$cores_per_node )) ]; then
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
  export OMP_NUM_THREADS=$threads_per_rank
  echo -n "$0: Begin time: "; date
  run_milc
  echo -n "$0: End time: "; date

