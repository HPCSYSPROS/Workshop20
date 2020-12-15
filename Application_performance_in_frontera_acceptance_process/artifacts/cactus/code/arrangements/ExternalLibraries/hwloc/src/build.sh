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
THORN=hwloc
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
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
HWLOC_DIR=${INSTALL_DIR}

# Set up environment
export CPPFLAGS="${CPPFLAGS} $(echo $(for dir in ${SYS_INC_DIRS}; do echo '' -I${dir}; done))"
export LDFLAGS="${LDFLAGS} $(echo $(for dir in ${LIBDIRS}; do echo '' -L${dir} -Wl,-rpath,${dir}; done))"
unset CPP
unset LIBS
unset LD
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi
export HWLOC_PCIUTILS_CFLAGS="$(echo '' $(for dir in ${PCIUTILS_INC_DIRS} ${ZLIB_INC_DIRS}; do echo '' $dir; done | sed -e 's/^ /-I/'))"
export HWLOC_PCIUTILS_LIBS="$(echo '' $(for dir in ${PCIUTILS_LIB_DIRS} ${ZLIB_LIB_DIRS}; do echo '' $dir; done | sed -e 's/^ /-L/') $(for dir in ${PCIUTILS_LIBS} ${ZLIB_LIBS}; do echo '' $dir; done | sed -e 's/^ /-l/'))"
echo "hwloc: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "hwloc: Unpacking archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/${NAME}.tar.gz

echo "hwloc: Configuring..."
cd ${NAME}
# Provide a special option for Blue Gene/Q; this is a cross-compile,
# so we can't easily detect this automatically
if echo ${CC} | grep -q bgxlc; then
    bgq='--host=powerpc64-bgq-linux'
else
    bgq=''
fi
# Disable Cairo and XML explicitly, since configure may pick it up if
# it is installed on the system, but our final link line may not link
# against these libraries. (We could use our own libxml2 library if we
# want.)
if test -n "${HAVE_CAPABILITY_PCIUTILS}"; then
    handle_pci='--enable-libpci'
else
    handle_pci='--disable-pci'
fi
## Disable pciaccess by forcing compiler errors
#export HWLOC_PCIACCESS_CFLAGS=DISABLE-PCIACCESS
./configure --prefix=${HWLOC_DIR} ${bgq} ${handle_pci} --disable-cairo --disable-libxml2 --disable-cuda --disable-opencl --enable-shared=no --enable-static=yes

echo "hwloc: Building..."
${MAKE}

echo "hwloc: Installing..."
${MAKE} install
popd

echo "hwloc: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "hwloc: Done."
