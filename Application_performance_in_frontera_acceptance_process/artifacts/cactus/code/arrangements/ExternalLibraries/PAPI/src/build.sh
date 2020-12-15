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
THORN=PAPI
NAME=papi-5.3.0
TARNAME=papi-5.3.0
SRCDIR="$(dirname $0)"
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${PAPI_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "BEGIN MESSAGE"
    echo "Installing PAPI into ${PAPI_INSTALL_DIR} "
    echo "END MESSAGE"
    INSTALL_DIR=${PAPI_INSTALL_DIR}
fi
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
PAPI_DIR=${INSTALL_DIR}
    
# Set up environment
unset CPP
unset LIBS
export MPICC=:          # disable MPI tests
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi

echo "PAPI: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "PAPI: Unpacking archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/${TARNAME}.tar.gz

echo "PAPI: Applying patches..."
pushd ${NAME}
# Replace <malloc.h> by <stdlib.h>
find . -type f -print | xargs perl -pi -e 's/malloc.h/stdlib.h/'
# disable Werror since new compiler warns about PAPI
find . -name config.mk -print | xargs perl -pi -e 's/-Werror//g'
popd

echo "PAPI: Configuring..."
cd ${NAME}/src
# CC          C compiler command
# CFLAGS      C compiler flags
# LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
#             nonstandard directory <lib dir>
# CPPFLAGS    C/C++ preprocessor flags, e.g. -I<include dir> if you have
#             headers in a nonstandard directory <include dir>
# F77         Fortran 77 compiler command
# FFLAGS      Fortran 77 compiler flags
# CPP         C preprocessor
#unset CC
#unset CFLAGS
#unset LDFLAGS
#unset CPPFLAGS
#unset F77
#unset FFLAGS
#unset CPP
#export MPICC=

./configure --prefix=${PAPI_DIR} --with-shared-lib=no

echo "PAPI: Building..."
${MAKE}

echo "PAPI: Installing..."
${MAKE} install
popd

echo "PAPI: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "PAPI: Done."
