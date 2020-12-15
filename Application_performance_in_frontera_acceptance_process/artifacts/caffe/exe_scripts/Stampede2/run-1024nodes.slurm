#!/bin/bash
#SBATCH -J Caffe-1024nodes	   # job name
#SBATCH -o Caffe-1024nodes-%j.out         # output and error file name (%j expands to jobID)
#SBATCH -N 1024              # total number of nodes
#SBATCH -n 1024
#SBATCH -p test2	   # queue (partition) -- normal, development, etc.
#SBATCH -t 3:00:00        # run time (hh:mm:ss) - 4 hours
#SBATCH -A A-ccsc    #TG-CCR150011    # project name

export OMP_NUM_THREADS=48
export PATH=/scratch/00946/zzhang/skylake/caffe/protocol-buffer/bin:$PATH

module load hdf5/1.8.16
export MLSL_ROOT=/scratch/00946/zzhang/skylake/caffe/caffe-intel-1.0.3-AVX512/external/mlsl/l_mlsl_2017.1.016
export LD_LIBRARY_PATH=/scratch/hpc_tools/apps/intel18/boost-common-avx512/1.65/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/protocol-buffer/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/glog/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/gflags/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$TACC_HDF5_LIB:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/leveldb/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/lmdb/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/opencv/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/snappy/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/mkl-dnn/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/scratch/00946/zzhang/skylake/caffe/caffe-intel-1.0.3-AVX512/external/mlsl/l_mlsl_2017.1.016/intel64/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries/linux/mkl/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/apps/intel17/python/2.7.13/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/00946/zzhang/skylake/caffe/caffe-intel-1.0.3-AVX512/external/mkl/mklml_lnx_2017.0.2.20170110/lib


date +"%T"
cp -r /scratch/00946/zzhang/stampede2/caffe-knl/caffe-intel-mkldnn-mpi/examples/imagenet/ilsvrc12_*_small_lmdb /dev/shm/
cp /scratch/00946/zzhang/stampede2/caffe-knl/caffe-intel-mkldnn-mpi/data/ilsvrc12/imagenet_mean_small.binaryproto /dev/shm/

/scratch/00946/zzhang/stampede2/imagenet/broadcast-mpi.sh /dev/shm/ilsvrc12_train_small_lmdb/data.mdb 1024
sleep 5
/scratch/00946/zzhang/stampede2/imagenet/broadcast-mpi.sh /dev/shm/ilsvrc12_val_small_lmdb/data.mdb 1024
sleep 5
/scratch/00946/zzhang/stampede2/imagenet/broadcast-mpi.sh /dev/shm/imagenet_mean_small.binaryproto 1024
sleep 5

date +"%T"

for i in 1;
do
  echo "Running ResNet-50 Case ${i}"
  time ibrun -np 1024 /scratch/00946/zzhang/skylake/caffe/caffe-intel-1.0.3-AVX512/build/tools/caffe train -engine "MKL2017" --solver=models/default_resnet_50/solver_small.prototxt > log-resnet50-1024skx-mkl2017-${i}.log 2>&1
  sleep 5
done

