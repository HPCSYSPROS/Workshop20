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

if [ -z "${PAPI_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "PAPI selected, but PAPI_DIR not set. Checking some places..."
    echo "END MESSAGE"
    
    DIRS="/usr /usr/local /opt/local ${HOME}"
    # look into each directory
    for dir in $DIRS; do
        # libraries might have different file extensions
        for libext in a dll dll.a dylib lib so; do
            # libraries can be in lib or lib64 (or libx32?)
            for libdir in lib64 lib; do
                # These files must exist
                FILES="include/papi.h ${libdir}/libpapi.${libext}"
                # assume this is the one and check all needed files
                PAPI_DIR="$dir"
                for file in $FILES; do
                    # discard this directory if one file was not found
                    if [ ! -r "$dir/$file" ]; then
                        unset PAPI_DIR
                        break
                    fi
                done
                # don't look further if all files have been found
                if [ -n "$PAPI_DIR" ]; then
                    break
                fi
            done
            # don't look further if all files have been found
            if [ -n "$PAPI_DIR" ]; then
                break
            fi
        done
        # don't look further if all files have been found
        if [ -n "$PAPI_DIR" ]; then
            break
        fi
    done
    
    if [ -z "$PAPI_DIR" ]; then
        echo "BEGIN MESSAGE"
        echo "PAPI not found"
        echo "END MESSAGE"
    else
        echo "BEGIN MESSAGE"
        echo "Found PAPI in ${PAPI_DIR}"
        echo "END MESSAGE"
    fi
fi

THORN=PAPI



################################################################################
# Build
################################################################################

if [ -z "${PAPI_DIR}"                                            \
     -o "$(echo "${PAPI_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Using bundled PAPI..."
    echo "END MESSAGE"
    
    # Set locations
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
    PAPI_BUILD=1
    PAPI_DIR=${INSTALL_DIR}
else    
    PAPI_BUILD=
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
echo "PAPI_BUILD       = ${PAPI_BUILD}"
echo "PAPI_INSTALL_DIR = ${PAPI_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
PAPI_INC_DIRS="${PAPI_DIR}/include"
PAPI_LIB_DIRS="${PAPI_DIR}/lib"
PAPI_LIBS="papi"

if nm ${PAPI_LIB_DIRS}/libpapi.a 2>/dev/null | grep -q pm_initialize; then
    PAPI_LIBS="${PAPI_LIBS} pmapi"
fi

PAPI_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${PAPI_INC_DIRS})"
PAPI_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${PAPI_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "PAPI_DIR      = ${PAPI_DIR}"
echo "PAPI_INC_DIRS = ${PAPI_INC_DIRS}"
echo "PAPI_LIB_DIRS = ${PAPI_LIB_DIRS}"
echo "PAPI_LIBS     = ${PAPI_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(PAPI_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(PAPI_LIB_DIRS)'
echo 'LIBRARY           $(PAPI_LIBS)'
