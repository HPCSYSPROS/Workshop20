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

if [ -z "${HWLOC_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "hwloc selected, but HWLOC_DIR not set. Checking some places..."
    echo "END MESSAGE"
    
    DIRS="/usr /usr/local /usr/local/packages /usr/local/apps /opt/local ${HOME} c:/packages"
    for dir in $DIRS; do
        # libraries might have different file extensions
        for libext in a dll dll.a dylib lib so; do
            # libraries can be in lib or lib64 (or libx32?)
            for libdir in lib64 lib/x86_64-linux-gnu lib lib/i386-linux-gnu; do
                FILES="include/hwloc.h $libdir/libhwloc.$libext"
                # assume this is the one and check all needed files
                HWLOC_DIR="$dir"
                for file in $FILES; do
                    # discard this directory if one file was not found
                    if [ ! -r "$dir/$file" ]; then
                        unset HWLOC_DIR
                        break
                    fi
                done
                # don't look further if all files have been found
                if [ -n "$HWLOC_DIR" ]; then
                    hwloc_lib_dir="$HWLOC_DIR/$libdir"
                    break
                fi
            done
            # don't look further if all files have been found
            if [ -n "$HWLOC_DIR" ]; then
                break
            fi
        done
        # don't look further if all files have been found
        if [ -n "$HWLOC_DIR" ]; then
            break
        fi
    done
    
    if [ -z "$HWLOC_DIR" ]; then
        echo "BEGIN MESSAGE"
        echo "hwloc not found"
        echo "END MESSAGE"
    else
        echo "BEGIN MESSAGE"
        echo "Found hwloc in ${HWLOC_DIR}"
        echo "END MESSAGE"
        # Check that version is sufficient
        export PKG_CONFIG_PATH=${hwloc_lib_dir}/pkgconfig:${PKG_CONFIG_PATH}
        # we negate the return code since in perl true == 1 but for the shell true == 0
        if ( pkg-config hwloc && pkg-config --modversion hwloc | 
               perl -ne 'm/^0*(\d+)[.]0*(\d+)/; exit !($1 < 1 or $2 < 6)' ) || \
           ( [ -r ${HWLOC_DIR}/include/hwloc.h ] &&
             perl -ne 'exit !($1 lt "0x00010600") if m/^#define HWLOC_API_VERSION (.*)/' \
               ${HWLOC_DIR}/include/hwloc.h ) ; then
            echo "BEGIN MESSAGE"
            echo "hwloc too old (require at least version 1.6)"
            echo "END MESSAGE"
            HWLOC_DIR='BUILD'
        fi
    fi
fi

THORN=hwloc

################################################################################
# Build
################################################################################

if [ -z "${HWLOC_DIR}"                                                  \
     -o "$(echo "${HWLOC_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Using bundled hwloc..."
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
    NAME=hwloc-1.10.1
    SRCDIR="$(dirname $0)"
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${HWLOC_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing hwloc into ${HWLOC_INSTALL_DIR}"
        echo "END MESSAGE"
        INSTALL_DIR=${HWLOC_INSTALL_DIR}
    fi
    HWLOC_BUILD=1
    HWLOC_DIR=${INSTALL_DIR}
    HWLOC_INC_DIRS="${HWLOC_DIR}/include"
    HWLOC_LIB_DIRS="${HWLOC_DIR}/lib"
    HWLOC_LIBS='hwloc'
else
    HWLOC_BUILD=
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    if [ ! -e ${DONE_FILE} ]; then
        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
        date > ${DONE_FILE}
    fi
    
    if [ -z "${hwloc_lib_dir}" ]; then
        hwloc_lib_dir="${HWLOC_DIR}/lib"
    fi
    
    # Check whether pkg-config works
    export PKG_CONFIG_PATH=${hwloc_lib_dir}/pkgconfig:${PCIUTILS_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH}
    if ! pkg-config hwloc; then
        echo "BEGIN MESSAGE"
        echo "pkg-config not found; attempting to use reasonable defaults"
        echo "END MESSAGE"
        
        HWLOC_INC_DIRS="${HWLOC_DIR}/include"
        HWLOC_LIB_DIRS="${HWLOC_DIR}/lib"
        HWLOC_LIBS='hwloc'
    else
        inc_dirs="$(pkg-config hwloc --static --cflags 2>/dev/null || pkg-config hwloc --cflags)"
        lib_dirs="$(pkg-config hwloc --static --libs 2>/dev/null || pkg-config hwloc --libs)"
        libs="$(pkg-config hwloc --static --libs 2>/dev/null || pkg-config hwloc --libs)"
        # Translate option flags into Cactus options:
        # - for INC_DIRS, remove -I prefix from flags
        # - for LIB_DIRS, remove all -l flags, and remove -L prefix from flags
        # - for LIBS, keep only -l flags, and remove -l prefix from flags
        HWLOC_INC_DIRS="$(echo '' $(for flag in $inc_dirs; do echo '' $flag; done | sed -e 's/^ -I//'))"
        HWLOC_LIB_DIRS="$(echo '' $(for flag in $lib_dirs; do echo '' $flag; done | grep -v '^ -l' | sed -e 's/^ -L//'))"
        HWLOC_LIBS="$(echo '' $(for flag in $libs; do echo '' $flag; done | grep '^ -l' | sed -e 's/^ -l//'))"
    fi
    
    # Add libnuma manually, if necessary
    if grep -q '[-]lnuma' ${hwloc_lib_dir}/libhwloc.la 2>/dev/null; then
        if ! echo '' ${HWLOC_LIBS} '' | grep -q ' numa '; then
            HWLOC_LIBS="${HWLOC_LIBS} numa"
        fi
    fi

fi



################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "HWLOC_BUILD       = ${HWLOC_BUILD}"
echo "HWLOC_INSTALL_DIR = ${HWLOC_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

HWLOC_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${HWLOC_INC_DIRS})"
HWLOC_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${HWLOC_LIB_DIRS})"
HWLOC_LIBS="${HWLOC_LIBS} ${HWLOC_EXTRA_LIBS}"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "HWLOC_DIR      = ${HWLOC_DIR}"
echo "HWLOC_INC_DIRS = ${HWLOC_INC_DIRS}"
echo "HWLOC_LIB_DIRS = ${HWLOC_LIB_DIRS}"
echo "HWLOC_LIBS     = ${HWLOC_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(HWLOC_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(HWLOC_LIB_DIRS)'
echo 'LIBRARY           $(HWLOC_LIBS)'
