#!/bin/bash

#Run only on a subset of cores avoiding core 6 and the last (68-64+1= 3 ) cores
#export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=6,65,66,67,74,133,134,135,142,201,202,203,210,269,270,271
#export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=4,5,34,35,72,73,104,105,140,141,172,173,208,209,240,241

mydir=`pwd`

## #Hybrid tasks
###1680 nodes 

for node_task_thread in 1680_3360_24; do
  echo $node_task_thread 
  node=`echo $node_task_thread | awk -F_ '{print $1}'`
  tpn=`echo $node_task_thread | awk -F_ '{print int($2/$1)}'`
  task=`echo $node_task_thread | awk -F_ '{print $2}'`
  thread=`echo $node_task_thread | awk -F_ '{print $3}'`
  echo " node = $node tpn = $tpn task = $task thread = $thread "

  #Set the tasks per node
  export TACC_TASKS_PER_NODE=$tpn
  export TACC_NODECT=$node
  
  ./run_tests_spp.sh $task $thread
done

