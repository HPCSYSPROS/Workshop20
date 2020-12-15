#! /bin/bash

################################################################################
# Build
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



# Set locations
THORN=GSL
NAME=gsl-1.16
SRCDIR="$(dirname $0)"
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${GSL_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "BEGIN MESSAGE"
    echo "Installing GSL into ${GSL_INSTALL_DIR} "
    echo "END MESSAGE"
    INSTALL_DIR=${GSL_INSTALL_DIR}
fi
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
GSL_DIR=${INSTALL_DIR}

# Set up environment
unset LIBS
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi

echo "GSL: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "GSL: Unpacking archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/${NAME}.tar.gz
pushd ${NAME}
${PATCH?} -p1 < ${SRCDIR}/../dist/stdarg.patch
# Some (ancient but still used) versions of patch don't support the
# patch format used here but also don't report an error using the exit
# code. So we use this patch to test for this
${PATCH?} -p1 < ${SRCDIR}/../dist/patchtest.patch
if [ ! -e .patch_tmp ]; then
    echo 'BEGIN ERROR'
    echo 'The version of patch is too old to understand this patch format.'
    echo 'Please set the PATCH environment variable to a more recent '
    echo 'version of the patch command.'
    echo 'END ERROR'
    exit 1
fi
rm -f .patch_tmp
popd

echo "GSL: Configuring..."
cd ${NAME}
./configure --prefix=${GSL_DIR} --enable-shared=no

echo "GSL: Building..."
${MAKE}

echo "GSL: Installing..."
${MAKE} install
popd

echo "GSL: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "GSL: Done."
