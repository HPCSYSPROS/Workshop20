#!/bin/bash

#script to collect memory b/w data and ib traffic
# usage: 

executable=$1
frequency=performance
export I_MPI_STATS=0

data_dir=${executable}_${frequency}_${SLURM_JOBID}

if [ x"$ANTIAFFINITY" == "x1" ]; then
    data_dir=anti_$data_dir
fi

mkdir -p $data_dir

export I_MPI_STATS_FILE=$data_dir/mpi_stats.txt

export HOST=`hostname -s`

#Determine local MPI rank
my_rank=$(( ${PMI_RANK-0} + ${PMI_ID-0} + ${MPIRUN_RANK-0} + ${OMPI_COMM_WORLD_RANK-0} + ${OMPI_MCA_ns_nds_vpid-0} + ${SLURM_PROCID-0} + ${ibrun_o_option-0} ))

#set -x 
SCHEDULER=SLURM
if [ "$SCHEDULER" == "SGE" ];then
    myway=`builtin echo $PE | /bin/sed s/way//`
elif [ "$SCHEDULER" == "SLURM" ];then
  # Figure out number of MPI tasks per node.   Added by McCalpin 2012-12-05
  # If running under "ibrun", NODE_TASKS_PPN_INFO will already be set
  # else get info from SLURM_TASKS_PER_NODE (not propagated by "ibrun" due to illegal characters)
  if [ -z "$NODE_TASKS_PPN_INFO" ]
  then
    myway=`echo $SLURM_TASKS_PER_NODE | awk -F '(' '{print $1}'`
  else
    #Because slurm will spread tasks evenly across nodes, the wayness of each node 
    # may be different.  The env variable NODE_TASKS_PPN_INFO propagates the wayness 
    # per node cluster with this format:
    #"{# of tasks per node},{#initial task id}_[repeats if necessary]"
    NODE_TASKS_PPN_INFO=`echo $NODE_TASKS_PPN_INFO | sed -e's/_/ /g'`
    for cluster in $NODE_TASKS_PPN_INFO ; do 
        way=`echo $cluster | awk -F ',' '{print $1}'` ; 
        task_cutoff=`echo $cluster | awk -F ',' '{print $2}'`; 
	if [ $my_rank -ge $task_cutoff ] ; then
           myway=$way
           mytask_cutoff=$task_cutoff
        fi
    done 
  fi
 
else
    echo "ERROR: Unknown batch system"
    exit 1
fi

#Calculate local MPI rank 
local_rank=$(( ( $my_rank - $mytask_cutoff) % $myway ))

MY_LOCAL_RANK=$local_rank

if [ "$MY_LOCAL_RANK" = "0" ]; then
   if [ "$frequency" != "performance" ]; then
      /scratch/hpc_tools/TACC_HWP_set --freq $frequency
   fi
   unbuffer /scratch/hpc_tools/read_imc_skx > $data_dir/imc.$HOST.$SLURM_JOBID
  #Capture network info
  if [[ -e /sys/class/infiniband/hfi1_0/ports/1/hw_counters/RxWords ]] ; then
    IB_RECV_FILE=/sys/class/infiniband/hfi1_0/ports/1/hw_counters/RxWords
    IB_XMIT_FILE=/sys/class/infiniband/hfi1_0/ports/1/hw_counters/TxWords
    RX_WORDS_BEFORE=`cat $IB_RECV_FILE`
    TX_WORDS_BEFORE=`cat $IB_XMIT_FILE`
  fi
fi


perf stat -x' ' -o $data_dir/freq.${HOST}.${PMI_RANK}${MPIRUN_RANK}.${SLURM_JOBID} $*

if [ "$MY_LOCAL_RANK" = "0" ]; then
  # Capture memory data
  /scratch/hpc_tools/read_imc_skx >> $data_dir/imc.$HOST.$SLURM_JOBID

  # Capture network data
  if [[ -e /sys/class/infiniband/hfi1_0/ports/1/hw_counters/RxWords ]] ; then
    RX_WORDS_AFTER=`cat $IB_RECV_FILE`
    TX_WORDS_AFTER=`cat $IB_XMIT_FILE`
    echo $(( ($RX_WORDS_AFTER - $RX_WORDS_BEFORE) * 4 )) >  ${data_dir}/ib_bytes.$HOST.$SLURM_JOBID
    echo $(( ($TX_WORDS_AFTER - $TX_WORDS_BEFORE) * 4 )) >> ${data_dir}/ib_bytes.$HOST.$SLURM_JOBID
  fi
fi

##################################################


