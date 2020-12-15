#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



################################################################################
# Search
################################################################################

# TODO: Look not just for any BLAS, but for OpenBLAS specifically to
# ensure we are using an efficient BLAS

# if [ -z "${OPENBLAS_DIR}" ]; then
#     echo "BEGIN MESSAGE"
#     echo "OpenBLAS selected, but OPENBLAS_DIR not set. Checking some places..."
#     echo "END MESSAGE"
#     
#     FILES="include/gsl/gsl_math.h"
#     DIRS="/usr /usr/local /usr/local/gsl /usr/local/packages/gsl /usr/local/apps/gsl ${HOME} c:/packages/gsl"
#     for dir in $DIRS; do
#         OPENBLAS_DIR="$dir"
#         for file in $FILES; do
#             if [ ! -r "$dir/$file" ]; then
#                 unset OPENBLAS_DIR
#                 break
#             fi
#         done
#         if [ -n "$OPENBLAS_DIR" ]; then
#             break
#         fi
#     done
#     
#     if [ -z "$OPENBLAS_DIR" ]; then
#         echo "BEGIN MESSAGE"
#         echo "OpenBLAS not found"
#         echo "END MESSAGE"
#     else
#         echo "BEGIN MESSAGE"
#         echo "Found OpenBLAS in ${OPENBLAS_DIR}"
#         echo "END MESSAGE"
#     fi
# fi

THORN=OpenBLAS



################################################################################
# Configure
################################################################################

: ${OPENBLAS_INT8:=0}
export OPENBLAS_INT8
if [ "${OPENBLAS_INT8}" != 0 -a "${OPENBLAS_INT8}" != 1 ]; then
    echo "BEGIN ERROR"
    echo "OPENBLAS_INT8 must be either 0 or 1"
    echo "END ERROR"
fi

################################################################################
# Build
################################################################################

if [ -z "${OPENBLAS_DIR}"                                            \
     -o "$(echo "${OPENBLAS_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Using bundled OpenBLAS..."
    echo "END MESSAGE"
    
    # Check for required tools. Do this here so that we don't require
    # them when using the system library.
    if [ "x$TAR" = x ] ; then
        echo 'BEGIN ERROR'
        echo 'Could not find tar command.'
        echo 'Please make sure that the (GNU) tar command is present,'
        echo 'and that the TAR variable is set to its location.'
        echo 'END ERROR'
        exit 1
    fi
    if [ "x$PATCH" = x ] ; then
        echo 'BEGIN ERROR'
        echo 'Could not find patch command.'
        echo 'Please make sure that the patch command is present,'
        echo 'and that the PATCH variable is set to its location.'
        echo 'END ERROR'
        exit 1
    fi

    # Set locations
    SRCDIR="$(dirname $0)"
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${OPENBLAS_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing OpenBLAS into ${OPENBLAS_INSTALL_DIR}"
        echo "END MESSAGE"
        INSTALL_DIR=${OPENBLAS_INSTALL_DIR}
    fi
    OPENBLAS_BUILD=1
    OPENBLAS_DIR=${INSTALL_DIR}
else
    OPENBLAS_BUILD=
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    if [ ! -e ${DONE_FILE} ]; then
        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
        date > ${DONE_FILE}
    fi
fi



################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "OPENBLAS_BUILD       = ${OPENBLAS_BUILD}"
echo "OPENBLAS_INSTALL_DIR = ${OPENBLAS_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
if [ "${OPENBLAS_DIR}" != 'NO_BUILD' ]; then
    : ${OPENBLAS_INC_DIRS="${OPENBLAS_DIR}/include"}
    : ${OPENBLAS_LIB_DIRS="${OPENBLAS_DIR}/lib"}
fi
: ${OPENBLAS_LIBS='openblas'}

OPENBLAS_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${OPENBLAS_INC_DIRS})"
OPENBLAS_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${OPENBLAS_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN DEFINE"
echo "CCTK_BLAS_INT8 ${OPENBLAS_INT8}"
echo "END DEFINE"
echo "BEGIN MAKE_DEFINITION"
echo "OPENBLAS_DIR      = ${OPENBLAS_DIR}"
echo "OPENBLAS_INC_DIRS = ${OPENBLAS_INC_DIRS}"
echo "OPENBLAS_LIB_DIRS = ${OPENBLAS_LIB_DIRS}"
echo "OPENBLAS_LIBS     = ${OPENBLAS_LIBS}"
echo "CCTK_BLAS_INT8    = ${OPENBLAS_INT8}"
echo "BLAS_DIR          = ${OPENBLAS_DIR}"
echo "BLAS_INC_DIRS     = ${OPENBLAS_INC_DIRS}"
echo "BLAS_LIB_DIRS     = ${OPENBLAS_LIB_DIRS}"
echo "BLAS_LIBS         = ${OPENBLAS_LIBS}"
echo "LAPACK_DIR        = ${OPENBLAS_DIR}"
echo "LAPACK_INC_DIRS   = ${OPENBLAS_INC_DIRS}"
echo "LAPACK_LIB_DIRS   = ${OPENBLAS_LIB_DIRS}"
echo "LAPACK_LIBS       = ${OPENBLAS_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(OPENBLAS_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(OPENBLAS_LIB_DIRS)'
echo 'LIBRARY           $(OPENBLAS_LIBS)'
