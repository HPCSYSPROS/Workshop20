#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors

. $CCTK_HOME/lib/make/bash_utils.sh

# Take care of requests to build the library in any case
GSL_DIR_INPUT=$GSL_DIR
if [ "$(echo "${GSL_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]; then
    GSL_BUILD=yes
    GSL_DIR=
else
    GSL_BUILD=
fi

# Try to find the library if build isn't explicitly requested
if [ -z "${GSL_BUILD}" ]; then
    find_lib GSL gsl 1 1.0 "gsl" "gsl/gsl_version.h" "$GSL_DIR"
fi

THORN=GSL

################################################################################
# Build
################################################################################

if [ -n "$GSL_BUILD" -o -z "${GSL_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "Using bundled GSL..."
    echo "END MESSAGE"
    
    check_tools "tar patch"

    # Set locations
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${GSL_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing GSL into ${GSL_INSTALL_DIR} "
        echo "END MESSAGE"
        INSTALL_DIR=${GSL_INSTALL_DIR}
    fi
    GSL_BUILD=1
    GSL_DIR=${INSTALL_DIR}
    GSL_INC_DIRS="$GSL_DIR/include"
    GSL_LIB_DIRS="$GSL_DIR/lib"
    GSL_LIBS="gsl gslcblas"
else
    GSL_BUILD=
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
echo "GSL_BUILD       = ${GSL_BUILD}"
echo "GSL_INSTALL_DIR = ${GSL_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
if [ -x ${GSL_DIR}/bin/gsl-config ]; then
    inc_dirs="$(${GSL_DIR}/bin/gsl-config --cflags)"
    lib_dirs="$(${GSL_DIR}/bin/gsl-config --libs)"
    libs="$(${GSL_DIR}/bin/gsl-config --libs)"
    # Translate option flags into Cactus options:
    # - for INC_DIRS, remove -I prefix from flags
    # - for LIB_DIRS, remove all -l flags, and remove -L prefix from flags
    # - for LIBS, keep only -l flags, and remove -l prefix from flags
    GSL_INC_DIRS="$(echo '' $(for flag in $inc_dirs; do echo '' $flag; done | sed -e 's/^ -I//'))"
    GSL_LIB_DIRS="$(echo '' $(for flag in $lib_dirs; do echo '' $flag; done | grep -v '^ -l' | sed -e 's/^ -L//'))"
    GSL_LIBS="$(echo '' $(for flag in $libs; do echo '' $flag; done | grep '^ -l' | sed -e 's/^ -l//'))"
fi

set_make_vars "GSL" "$GSL_LIBS" "$GSL_LIB_DIRS" "$GSL_INC_DIRS"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "GSL_DIR      = ${GSL_DIR}"
echo "GSL_INC_DIRS = ${GSL_INC_DIRS}"
echo "GSL_LIB_DIRS = ${GSL_LIB_DIRS}"
echo "GSL_LIBS     = ${GSL_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(GSL_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(GSL_LIB_DIRS)'
echo 'LIBRARY           $(GSL_LIBS)'
