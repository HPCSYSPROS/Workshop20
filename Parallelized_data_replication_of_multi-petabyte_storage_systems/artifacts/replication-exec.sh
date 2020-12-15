#!/bin/bash

export FSNAME=fs0
export FSPATH=/gpfs/${FSNAME}
export WORKTMP=${FSPATH}/ddn/tmp
export PATH=${FSPATH}/ddn/parallel/bin:${FSPATH}/ddn/rsync-ss/bin:$PATH

#############################################################################################
##                                                                                         ##
## Get remote host name.                                                                   ##
## E.g. if local host is s1-ces-1, the remote host is s2-ddngs-1                           ##
## ddngs is a virtual hostname with virtual IP, used in NFS cluster for high availability. ##
##                                                                                         ## 
#############################################################################################

REMOTEHOST=$(echo $HOSTNAME | sed 's/s1-ces/s2-ddngs/g')

case $4 in
    "fs-1")
        remotedir=fs-1-dr;;
    *)
        echo Error
        exit -1;;
esac

if [ $1 == "LIST" ]
then
    [[ -z $HOME ]] && export HOME=/root
    TEMPFILE=$(mktemp -p ${WORKTMP} --suffix=.$4)
    cat $2 | sed "s/.*\/gpfs\/${FSNAME}\/\.snapshots\/$3\/$4\\(.*\)/\1/g" | sed 's/^\///g' > $TEMPFILE
    parallel -j 8 rsync -RltgoWd --gpfs-attrs --rsync-path='/gpfs/fs0/ddn/rsync-ss/bin/rsync-ss' ${FSPATH}/.snapshots/$3/$4/./{} $REMOTEHOST:${FSPATH}/$remotedir/ :::: $TEMPFILE
    rm $TEMPFILE
else
   echo $@
fi

#######################################################################
##                                                                   ##
## /gpfs/fs0/ddn/rsync-ss/bin/rsync-ss is a wrapper at remote site.  ##
##                                                                   ##
##    cat /gpfs/fs0/ddn/rsync-ss/bin/rsync-ss                        ##
##    #!/bin/bash                                                    ##
##    /gpfs/fs0/ddn/rsync-ss/bin/rsync --gpfs-attrs $@               ##
##                                                                   ##
#######################################################################
