#!/bin/bash

function usage {
   echo "Usage:  run_tests.sh <# of MPI tasks> <# of OMP threads > "
}


#Set testing modes
#  denoted by queue
my_host=`hostname -s`
my_queue=$SLURM_JOB_PARTITION

#Possible configs
#  Flat-All2All        
#  Flat-Quadrant
#  Flat-SNC-4
#  Cache-All2All
#  Cache-Quadrant
#  Cache-SNC-4

echo "my_queue : $my_queue "
case $my_queue in
Flat-All2All)
  mcmode="flat_a2a"
  #run using mcdram only using numactl
  mcdram_enable="numactl --preferred=1 "
  ;;
#normal|Flat-Quadrant)
Flat-Quadrant)
  mcmode="flat_quad"
  mcdram_enable="numactl --preferred=1 "
  ;;
Flat-SNC-4)
  mcmode="flat_sn4"
  mcdram_enable="numactl --preferred=4,5,6,7 "
  ;;
Cache-All2All)
  mcmode="cache_a2a"
  mcdram_enable=""
  ;;
normal|Cache-Quadrant)
  mcmode="cache_quad"
  mcdram_enable=""
  ;;
Cache-SNC-4)
  mcmode="cache_sn4"
  mcdram_enable=""
  ;;
test1|test2)
  mcmode="ddr4"
  mcdram_enable="tacc_affinity"
  ;;
*)
  echo "Unknown configuration -- exiting "
  exit
  ;;
esac


affinity="spread"
#affinity="compact"
export OMP_PROC_BIND="${affinity}"
#export KMP_AFFINITY=verbose
export KMP_STACKSIZE=192m


if [ x$1 == "x" ]; then
  usage
  exit 1
else
  num_tasks=$1;
  if [ x$2 == "x" ]; then
    usage
    exit 1
  else
    num_threads=$2;
    if [ "$num_tasks" -lt 30 ]; then
        export KMP_STACKSIZE=512m
    fi
  fi
fi


for run in 1; do
#for run in 3  ; do

  # Create output dir and prefix
  out_prefix="spplarge_${mcmode}_${affinity}_${TACC_NODECT}_${num_tasks}p_${num_threads}t_${run}"
  echo " Running $out_prefix "
  mkdir -p ${out_prefix}

  #Set # of threads
  export OMP_NUM_THREADS=$num_threads
  #Save the env
  env > ${out_prefix}/env.out

  #Begin the run
  executable="wrf.exe"
  frequency=performance
  date  > ${out_prefix}/wrf_spp_${out_prefix}.out
  begin=`date +%s`

   echo " command:  ibrun -np $num_tasks $mcdram_enable ./instrumented_run.sh ./wrf.exe  >>  ${out_prefix}/wrf_spp_${out_prefix}.out " &>> ${out_prefix}/wrf_spp_${out_prefix}.out
   ibrun -np $num_tasks $mcdram_enable ./instrumented_run.sh ./wrf.exe  &>> ${out_prefix}/wrf_spp_${out_prefix}.out

  end=`date +%s`

  elapsed=$((end - begin ))
  
  echo "Elapsed:  ${elapsed} " &>> ${out_prefix}/wrf_spp_${out_prefix}.out 
# mv rsl.* ${out_prefix}

done

