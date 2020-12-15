#!/bin/bash
#PBS -l nodes=1:ppn=1:xkhimem
#PBS -l walltime=02:00:00
#PBS -N caffe
#PBS -e $PBS_JOBID.err
#PBS -o $PBS_JOBID.out

cd /scratch/sciteam/tp2u2/caffe/caffe-intel/caffe-1.0.0

module load cmake/3.1.3
module load gcc/4.9.3
module load cudatoolkit
export PATH=/scratch/sciteam/tp2u2/caffe/protocol-buffer/bin:$PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/protocol-buffer/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/hdf5/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/sw/xe/boost/1.63.0/sles11.3_gnu4.9.3/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/glog/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/gflags/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/leveldb/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/lmdb/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/opencv/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/snappy/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/mlsl/intel64/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/cudnn/cudnn-5.1/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/sciteam/tp2u2/caffe/caffe-intel/caffe-1.0.0/build/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=/scratch/sciteam/tp2u2/caffe/cudnn/cudnn-5.1/include:$CPLUS_INCLUDE_PATH

export MLSL_HOME=/scratch/sciteam/tp2u2/caffe/mlsl/intel64

date +"%T"
aprun -n 1 -N 1 cp -r /scratch/sciteam/tp2u2/imagenet/ilsvrc12_*_small_lmdb /tmp/
date +"%T"
aprun -n 1 -N 1 cp /scratch/sciteam/tp2u2/imagenet/imagenet_mean_small.binaryproto /tmp/
date +"%T"

for i in 1 2 3;
do
  echo "Running ResNet-50 Case ${i}"
  time aprun -n 1 -N 1 build/tools/caffe train --solver=models/default_resnet_50/solver_small.prototxt > log-resn
et50-1node-1K20-5k-${i}.log 2>&1
done
