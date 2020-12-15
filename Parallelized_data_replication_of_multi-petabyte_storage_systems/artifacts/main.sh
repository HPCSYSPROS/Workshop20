#!/bin/bash

echo
echo [$(date +%FT%H:%M:%S)] Setting execution environment
FSNAME=fs0
FSPATH=/gpfs/${FSNAME}
SCRIPTPATH=${FSPATH}/ddn/scripts/bin
WORKTMP=${FSPATH}/ddn/tmp
NFSMNT=/mnt
export PATH=/usr/lpp/mmfs/bin:${FSPATH}/ddn/parallel/bin:${FSPATH}/ddn/rsync-ss/bin:${FSPATH}/ddn/openmpi/bin:${FSPATH}/ddn/mpifileutils/bin:$PATH

[[ -f ${WORKTMP}/sync-lock ]] && { echo [$(date +%FT%H:%M:%S)] Previous process $(awk {'print $2'} ${WORKTMP}/sync-lock) is still running since $(awk {'print $1'} ${WORKTMP}/sync-lock).; exit -1; }
echo [$(date +%FT%H:%M:%S)] $$ > ${WORKTMP}/sync-lock

echo [$(date +%FT%H:%M:%S)] Creating hourly snapshot.
mmcrsnapshot ${FSNAME} hourly@$(date +%F_%H%M%S)
if [ $? != 0 ]
then
    echo [$(date +%FT%H:%M:%S)] "Hourly snapshot creation failed."
fi

snaplist=$(mmlssnapshot ${FSNAME} -s global | tail -n +3 | awk {'print $1'} | grep ^hourly@)

[ $(echo $snaplist | xargs -n1 | wc -l) -ge 2 ] || { echo [$(date +%FT%H:%M:%S)] This program cannot be executed with less than two snapshots; rm -f ${WORKTMP}/sync-lock; exit -1; }


echo [$(date +%FT%H:%M:%S)] Reviewing snapshots
SNAPS=($(echo $snaplist | xargs -n1 | tail -n 2 | xargs))
echo [$(date +%FT%H:%M:%S)] Comparing ${SNAPS[0]} and ${SNAPS[1]}

for FILESET in fs-1
do
    echo [$(date +%FT%H:%M:%S)] Processing $FILESET fileset.
    case $FILESET in
        "fs-1")
            REMOTEDIR="fs-1-dr"
            ;;
        *)
            echo [$(date +%FT%H:%M:%S)] Invalid fileset
            exit -1
            ;;
    esac

    for node in s1-ces-1 s1-ces-2 s1-ces-3
    do
        ping -c 1 $node > /dev/null 2>&1
        [[ $? == 0 ]] && DATAMOVERS="$DATAMOVERS $node"
    done

    DATAMOVERS=$(echo $DATAMOVERS | tr ' ' ',')

    echo [$(date +%FT%H:%M:%S)] Replicating list of new and updated files from local $FILESET to remote DR $REMOTEDIR. 
    echo [$(date +%FT%H:%M:%S)] Executing -- \'mmapplypolicy ${FSPATH}/.snapshots/${SNAPS[1]}/${FILESET}/ -P ${SCRIPTPATH}/replication.policy -M prevsnapid=\"SNAPID\(${SNAPS[0]}\)\" -M lastsnap=${SNAPS[1]} -M workfileset=$FILESET -N $DATAMOVERS -m 3 -B 1000\'
    mmapplypolicy ${FSPATH}/.snapshots/${SNAPS[1]}/${FILESET}/ -P ${SCRIPTPATH}/replication.policy -M prevsnapid="SNAPID(${SNAPS[0]})" -M lastsnap=${SNAPS[1]} -M workfileset=$FILESET -N $DATAMOVERS -m 3 -B 1000
    if [[ ! $(mount -t nfs4 | grep ${NFSMNT}/${FSNAME}) || ! $(mount -t nfs4 | grep ${NFSMNT}/${REMOTEDIR}) ]]
    then
        echo "[$(date +%FT%H:%M:%S)] Local or remote NFS mount point not available, unable to do final sync. "
    else
        echo "[$(date +%FT%H:%M:%S)] Syncing again to update soft links and deleted files at remote site. "
        if [[ $(ls ${WORKTMP}/tmp.*.chgprj.${SNAPS[1]}.$FILESET.lst 2>/dev/null) ]]
        then
            PROJECTS=$(sort -u ${WORKTMP}/tmp.*.chgprj.${SNAPS[1]}.$FILESET.lst)
        else
            PROJECTS=$(ls ${FSPATH}/.snapshots/${SNAPS[1]}/$FILESET)
        fi
        if [ -n "$PROJECTS" ]
        then
            for PROJECT in $PROJECTS
            do
                if [[ -d ${NFSMNT}/${FSNAME}/.snapshots/${SNAPS[1]}/${FILESET}/${PROJECT} && -d ${NFSMNT}/${REMOTEDIR}/${PROJECT} ]]
                then
                    echo [$(date +%FT%H:%M:%S)] Syncing $PROJECT
                    DATAMOVERS=$(parallel --timeout 2 -j3 'ping -c 1 {} >/dev/null 2>&1 && echo {}' ::: s1-ces-1 s1-ces-2 s1-ces-3 | xargs)
                    MPIPROCS=$(($(echo $DATAMOVERS | wc -w) * 8))
                    MPIHOSTS=$(echo $DATAMOVERS | tr ' ' ',' | sed 's/\(ces-[1-3]\)/\1:8/g')
                    echo [$(date +%FT%H:%M:%S)] Executing -- \'mpirun -allow-run-as-root -np $MPIPROCS -N 8 -host $MPIHOSTS -map-by socket -bind-to core -mca btl ^openib -mca pml ucx -x UCX_MAX_EAGER_RAILS=2 -x UCX_MAX_RNDV_RAILS=2 dsync -s -D ${NFSMNT}/${FSNAME}/.snapshots/${SNAPS[1]}/${FILESET}/${PROJECT} ${NFSMNT}/${REMOTEDIR}/${PROJECT}\'
                    mpirun -allow-run-as-root -np $MPIPROCS -N 8 -host $MPIHOSTS -map-by socket -bind-to core -mca btl ^openib -mca pml ucx -x UCX_MAX_EAGER_RAILS=2 -x UCX_MAX_RNDV_RAILS=2 dsync -s -D ${NFSMNT}/{FSNAME}/.snapshots/${SNAPS[1]}/${FILESET}/${PROJECT} ${NFSMNT}/${REMOTEDIR}/${PROJECT}
                fi
            done
        fi
        [[ $(ls ${WORKTMP}/tmp.*.chgprj.${SNAPS[1]}.$FILESET.lst 2>/dev/null) ]] && rm -f ${WORKTMP}/tmp.*.chgprj.${SNAPS[1]}.$FILESET.lst
    fi
done

echo [$(date +%FT%H:%M:%S)] Clearing system cache.
clush -w s1-ces-[1-3] "sync; echo 3 > /proc/sys/vm/drop_caches"
echo [$(date +%FT%H:%M:%S)] System cache cleared.

[[ -f ${WORKTMP}/sync-lock ]] && rm -f ${WORKTMP}/sync-lock

num_snaps=$(echo $snaplist | xargs -n1 | wc -l)
if [ ${num_snaps} -gt 12 ]
then
    echo [$(date +%FT%H:%M:%S)] Deleting old hourly snapshots.
    delsnaplist=$(echo $snaplist|xargs -n1 | head -n -12)
    for snap in $delsnaplist
    do
        echo [$(date +%FT%H:%M:%S)] Deleting $snap.
        DATAMOVERS=$(parallel --timeout 2 -j3 'ping -c 1 {} >/dev/null 2>&1 && echo {}' ::: s1-ces-1 s1-ces-2 s1-ces-3 | xargs | tr ' ' ',')
        mmdelsnapshot ${FSNAME} $snap -N $DATAMOVERS
        if [ $? != 0 ]
        then
            echo [$(date +%FT%H:%M:%S)] Deletion of $snap failed.
        else
            echo [$(date +%FT%H:%M:%S)] Deletion of $snap completed.
        fi
    done
fi
