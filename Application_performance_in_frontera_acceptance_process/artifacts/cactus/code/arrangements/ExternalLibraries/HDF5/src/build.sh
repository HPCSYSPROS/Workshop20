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
THORN=HDF5
NAME=hdf5-1.8.17
SRCDIR="$(dirname $0)"
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${HDF5_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "BEGIN MESSAGE"
    echo "Installing HDF5 into ${HDF5_INSTALL_DIR}"
    echo "END MESSAGE"
    INSTALL_DIR=${HDF5_INSTALL_DIR}
fi
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
HDF5_DIR=${INSTALL_DIR}

# Set up environment
if [[ ${F90} == none ]]; then
    echo 'BEGIN MESSAGE'
    echo 'No Fortran 90 compiler available. Building HDF5 library without Fortran support.'
    echo 'END MESSAGE'
    unset FC
    unset FCFLAGS
else
    export FC="${F90}"
    export FCFLAGS="${F90FLAGS}"
fi
for dir in $SYS_INC_DIRS; do
    CPPFLAGS="$CPPFLAGS -I$dir"
done
export CPPFLAGS
export LDFLAGS
unset CPP
unset LIBS
unset RPATH
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi

echo "HDF5: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

# Build core library
echo "HDF5: Unpacking archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/${NAME}.tar.gz

echo "HDF5: Configuring..."
cd ${NAME}
# Do not build C++ API if it has been disabled.
: ${HDF5_ENABLE_CXX:=yes}
# Do not build Fortran API if it has been disabled, or if there is no
# Fortran 90 compiler.
if [[ -n ${FC} ]]; then
    : ${HDF5_ENABLE_FORTRAN:=yes}
else
    HDF5_ENABLE_FORTRAN=no
fi
./configure --prefix=${HDF5_DIR} --with-zlib=${ZLIB_DIR} --enable-cxx=${HDF5_ENABLE_CXX} --enable-fortran=${HDF5_ENABLE_FORTRAN} --enable-fortran2003=${HDF5_ENABLE_FORTRAN} --disable-shared --enable-static-exec

echo "HDF5: Building..."
${MAKE}

echo "HDF5: Installing..."
${MAKE} install
popd

# Build checker
echo "HDF5: Unpacking checker archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/h5check_2_0.tar.gz

echo "HDF5: Configuring checker..."
cd h5check_2_0
# Point the checker to the just-installed library
export CPPFLAGS="${CPPFLAGS} -I${HDF5_DIR}/include"
export LDFLAGS="${LDFLAGS} ${LIBDIR_PREFIX}${HDF5_DIR}/lib ${RUNDIR_PREFIX}${HDF5_DIR}/lib"
export H5CC="${CC}"
export H5CC_PP="${CPP}"
export H5FC="${FC}"
export H5FC_PP="${FPP}"
export H5CPP="${CXX}"
./configure --prefix=${HDF5_DIR} --with-zlib=${ZLIB_DIR}

echo "HDF5: Building checker..."
#${MAKE}
(cd src && ${MAKE})
(cd tool && ${MAKE})

echo "HDF5: Installing checker..."
# The build fails in the "test" subdirectory, because
# /usr/include/hdf5.h (if it exists) is used instead of the the one we
# just installed. We therefore skip the build in the "test"
# subdirectory.
#${MAKE} install
(cd src && ${MAKE} install)
(cd tool && ${MAKE} install)
popd

echo "HDF5: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "HDF5: Done."
