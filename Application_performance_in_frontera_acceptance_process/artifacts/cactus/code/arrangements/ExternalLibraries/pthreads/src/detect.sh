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
# Check for old mechanism
################################################################################

if [ -n "${PTHREADS}" ]; then
    echo 'BEGIN ERROR'
    echo "Setting the option \"PTHREADS\" is incompatible with the Pthreads thorn. Please remove the option PTHREADS=${PTHREADS}."
    echo 'END ERROR'
    exit 1
fi



################################################################################
# Set default libraries to link if nothing provided by user
################################################################################
: ${PTHREADS_LIBS=pthread}

################################################################################
# Search
################################################################################

if [ -z "${PTHREADS_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "PTHREADS selected, but PTHREADS_DIR not set. Checking some places..."
    echo "END MESSAGE"
    
    # We look in these directories
    DIRS="/usr /usr/local /opt/local ${HOME} c:/packages/PTHREADS"
    # look into each directory
    for dir in $DIRS; do
        # libraries might have different file extensions
        for libext in a so dylib; do
            # libraries can be in /lib or /lib64
            for libdir in lib64 lib lib/x86_64-linux-gnu lib/i386-linux-gnu lib/arm-linux-gnueabihf; do
                # These files must exist
                FILES="include/pthread.h $(for lib in ${PTHREADS_LIBS}; do echo ${libdir}/lib${lib}.${libext}; done)"
                # assume this is the one and check all needed files
                PTHREADS_DIR="$dir"
                for file in $FILES; do
                    # discard this directory if one file was not found
                    if [ ! -r "$dir/$file" ]; then
                        unset PTHREADS_DIR
                        break
                    fi
                done
                # don't look further if all files have been found
                if [ -n "$PTHREADS_DIR" ]; then
                    break
                fi
           done
           # don't look further if all files have been found
           if [ -n "$PTHREADS_DIR" ]; then
               break
           fi
        done
        # don't look further if all files have been found
        if [ -n "$PTHREADS_DIR" ]; then
            break
        fi
    done
    
    if [ -z "$PTHREADS_DIR" ]; then
        echo "BEGIN ERROR"
        echo "Did not find PTHREADS"
        echo "END ERROR"
        exit 1
    else
        echo "BEGIN MESSAGE"
        echo "Found PTHREADS in ${PTHREADS_DIR}"
        echo "END MESSAGE"
    fi
fi



################################################################################
# Configure Cactus
################################################################################

# PTHREADS_LIBS and PTHREADS_DIR are already set; set remaining unset
# options based on these
if [ "x$PTHREADS_DIR" != xNO_BUILD ] ; then
    : ${PTHREADS_INC_DIRS="${PTHREADS_DIR}/include"}
    : ${PTHREADS_LIB_DIRS="${PTHREADS_DIR}/lib"}
fi
PTHREADS_INC_DIRS=$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${PTHREADS_INC_DIRS})
PTHREADS_LIB_DIRS=$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${PTHREADS_LIB_DIRS})

# Pass options to Cactus
echo "BEGIN DEFINE"
echo "CCTK_PTHREADS 1"
echo "END DEFINE"

echo "BEGIN MAKE_DEFINITION"
echo "PTHREADS_DIR            = ${PTHREADS_DIR}"
echo "PTHREADS_INC_DIRS       = ${PTHREADS_INC_DIRS}"
echo "PTHREADS_LIB_DIRS       = ${PTHREADS_LIB_DIRS}"
echo "PTHREADS_LIBS           = ${PTHREADS_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(PTHREADS_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(PTHREADS_LIB_DIRS)'
echo 'LIBRARY           $(PTHREADS_LIBS)'
