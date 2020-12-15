#!/bin/bash

FSNAME=fs0
FSPATH=/gpfs/${FSNAME}
WORKTMP=${FSPATH}/ddn/tmp

if [ $1 == "LIST" ]
then
    TEMPFILE=$(mktemp -p $WORKTMP --suffix=.chgprj.$3.$4.lst)
    awk -F/ {'print $7'} $2 >> $TEMPFILE
else
   echo $@
fi
