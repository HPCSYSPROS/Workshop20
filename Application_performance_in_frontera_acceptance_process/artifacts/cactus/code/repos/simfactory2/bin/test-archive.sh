#! /bin/bash

# Archive Cactus test results

set -e
set -x
set -u

# Read arguments
machine=$1         # machine name
source=$2          # source name
configuration=$3   # configuration name (identifying the source tree)
simulation=$4      # simulation name (identifying the run)
sourcedir=$5       # Cactus directory (sourcebasedir/source)
simdir=$6          # simulation directory (basedir)
sim1=$7            # first simulation name
nprocs1=$8         # first number of processes
sim2=$9            # second simulation name
nprocs2=${10}      # second number of processes



# Read configuration options
# TODO: Use Simfactory's MDB for this
if [ -r $sourcedir/cactus.config ]; then
    export CCTK_HOME=$sourcedir
    source $sourcedir/cactus.config
fi

if [ -z "${CACTUS_TESTSUITE_GIT_REPO:-}" -o     \
    -z "${CACTUS_TESTSUITE_SSH_HOST_ID:-}" -o   \
    -z "${CACTUS_TESTSUITE_SSH_ID:-}" ]
then
    echo "Configuration options are not set -- not pushing test results to repository"
    exit 0
fi



GIT_CMD=$(env PATH="$HOME/git/bin:$PATH" which git)

# Prepare git repository
resultsdir=$sourcedir/configs/$configuration/build-and-test.git
if [ ! -e $resultsdir/.git ]; then
    rm -rf $resultsdir
    mkdir $resultsdir
    cd $resultsdir
    $GIT_CMD init-db
    $GIT_CMD config receive.denyCurrentBranch false
fi
cd $resultsdir
rm -rf *



# Add source tree revision information
cp $sourcedir/REVISION_TIMESTAMP . || true
cp $sourcedir/REVISION_MANIFEST . || true

# TODO: Add real build log
echo '[Build log is not yet captured]' > build.log



# Copy test results
rsync -av                                               \
    --include "/*"                                      \
    --include "/*/SIMFACTORY"                           \
    --exclude "/*/SIMFACTORY/exe"                       \
    --exclude "/*/SIMFACTORY/exe/*"                     \
    --include "/*/SIMFACTORY/*"                         \
    --include "/*/SIMFACTORY/*/*"                       \
    --include "/*/output-0000"                          \
    --include "/*/output-0000/TEST"                     \
    --include "/*/output-0000/TEST/*"                   \
    --include "/*/output-0000/TEST/*/*"                 \
    --include "/*/output-0000/TEST/*/*/*.log"           \
    --include "/*/output-0000/TEST/*/*/*.diffs"         \
    --exclude "*"                                       \
    $simdir/$sim1 .

rsync -av                                               \
    --include "/*"                                      \
    --include "/*/SIMFACTORY"                           \
    --exclude "/*/SIMFACTORY/exe"                       \
    --exclude "/*/SIMFACTORY/exe/*"                     \
    --include "/*/SIMFACTORY/*"                         \
    --include "/*/SIMFACTORY/*/*"                       \
    --include "/*/output-0000"                          \
    --include "/*/output-0000/TEST"                     \
    --include "/*/output-0000/TEST/*"                   \
    --include "/*/output-0000/TEST/*/*"                 \
    --include "/*/output-0000/TEST/*/*/*.log"           \
    --include "/*/output-0000/TEST/*/*/*.diffs"         \
    --exclude "*"                                       \
    $simdir/$sim2 .

# Create summary log files if missing, ensuring that the directory
# hierarchy is not empty (and will thus not be ignored by git)
mkdir -p $sim1/output-0000/TEST/$configuration
mkdir -p $sim2/output-0000/TEST/$configuration
: >> $sim1/output-0000/TEST/$configuration/summary.log
: >> $sim2/output-0000/TEST/$configuration/summary.log

mv $sim1 build-and-test_$nprocs1
mv $sim2 build-and-test_$nprocs2




# Create release-info files
mkdir release-info
nthreads1=$(grep '^ *numthreads *=' $simdir/$sim1/output-0000/SIMFACTORY/properties.ini | sed -e 's/^[^=]*= *//')
nthreads2=$(grep '^ *numthreads *=' $simdir/$sim2/output-0000/SIMFACTORY/properties.ini | sed -e 's/^[^=]*= *//')
cp "build-and-test_$nprocs1/output-0000/TEST/$configuration/summary.log" \
    "release-info/${machine}__${nprocs1}_${nthreads1}.log"
cp "build-and-test_$nprocs2/output-0000/TEST/$configuration/summary.log" \
    "release-info/${machine}__${nprocs2}_${nthreads2}.log"



# Shorten all overly long files
for file in $(find . -size +1000k -type f -print); do
    echo "File '$file' is very large; shortening it"
    ls -l $file
    {
        head -n 1000 $file
        echo '[FILE TOO LONG -- TRUNCATED]'
        tail -n 1000 $file
    } > $file.tmp
    mv $file.tmp $file
done



# Wrap test results
$GIT_CMD add -A .
$GIT_CMD commit -m "Update build-and-test results for $machine-$source-$configuration"
$GIT_CMD tag $simulation
branchname="$machine-$source-$configuration"
branchname=$(echo $branchname | sed -e 's/-Cvanilla-/-/')
branchname=$(echo $branchname | sed -e 's/-sim-/-/')
branchname=$(echo $branchname | sed -e 's/-sim$//')
$GIT_CMD branch -f $branchname



# Create git ssh wrapper
export SSH_KNOWN_HOSTS=$resultsdir/known_hosts
echo $CACTUS_TESTSUITE_SSH_HOST_ID > $SSH_KNOWN_HOSTS
export GIT_SSH=$resultsdir/git-ssh.sh
# IdentityFile: use this identity
# IdentitiesOnly: don't try other identities; those might be accepted
#    by the remote ssh, but then rejected by git
# PasswordAuthentication: don't output a password prompt, we are not
#    interactive
# UserKnownHostsFile: allow connecting to the git server
echo ssh                                                \
    -o IdentityFile="'$CACTUS_TESTSUITE_SSH_ID'"        \
    -o IdentitiesOnly=yes                               \
    -o PasswordAuthentication=no                        \
    -o UserKnownHostsFile="'$SSH_KNOWN_HOSTS'"          \
    '"$@"' > $GIT_SSH
chmod a+x $GIT_SSH



# Push results
# TODO: push only the current branch, not all branches
#$GIT_CMD push -v -f --all $CACTUS_TESTSUITE_GIT_REPO
#$GIT_CMD push -v -f --tags $CACTUS_TESTSUITE_GIT_REPO
$GIT_CMD push -v -f --tags $CACTUS_TESTSUITE_GIT_REPO $branchname
