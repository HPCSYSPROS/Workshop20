#!/bin/bash

rm -rf configs/*

module load intel/18.0.0

## old: gmake zelmani-config options=stampede2-knl-intel18-O3.cfg
## new:
gmake zelmani-config options=stampede2-victor.cfg

cp ThornList configs/zelmani
gmake -j20 zelmani
#gmake -j20 zelmani
cp exe/cactus_zelmani ../
