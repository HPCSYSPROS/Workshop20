#!/bin/sh

if (( $# != 2 )); then {
  echo "usage: $0 <output_dir> <num_replicas>"
  exit -1
}; fi

mkdir $1

for (( i=0; i<$2; ++i )); do mkdir $1/$i; done

