#!/bin/sh

module swap PrgEnv-cray PrgEnv-intel
module load craype-hugepages2M
#module unload darshan xalt
#module load perftools-base
#module load perftools-lite
export CRAYPE_LINK_TYPE=static
cd ks_imp_rhmc/
make clean
rm ../libraries/*.[ao]
make -f Makefile_BW_noOMP su3_rhmd_hisq
cd ..
