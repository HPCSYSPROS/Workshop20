#!/bin/bash 

module load intel/16.0.3
tar -xvf NAMD_2.12_Source.tar.gz
cd NAMD_2.12_Source
tar -xvf charm-6.7.1.tar
cd charm-6.7.1/
env MPICXX=mpicxx CC=icc CXX=icpc ./build charm++ mpi-linux-x86_64 mpicxx --incdir /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include --libdir /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib --no-build-shared --with-production -j16 -xCORE-AVX512
cd ../src
patch < ../../../patches/Output.patch
# This patch will turn off writing output files (coordinates and velocities). Only for bechmarking to eliminate I/O at the end. Not for real production runs. It does not affect reported performance at all.
cd ..
sed -i 's/xMIC-AVX512/xCORE-AVX512/' arch/Linux-KNL-icc.arch
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64-threaded.tar.gz
tar -xvf tcl8.5.9-linux-x86_64-threaded.tar.gz
./config Linux-KNL-icc --with-memopt --charm-arch mpi-linux-x86_64-mpicxx --with-mkl --mkl-prefix /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/ --with-tcl --tcl-prefix `pwd`/tcl8.5.9-linux-x86_64-threaded
cd Linux-KNL-icc
make -j 16

