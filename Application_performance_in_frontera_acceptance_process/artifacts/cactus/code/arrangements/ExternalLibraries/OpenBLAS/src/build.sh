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
THORN=OpenBLAS
NAME=OpenBLAS-0.2.19
TARNAME=v0.2.19
SRCDIR="$(dirname $0)"
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${OPENBLAS_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "BEGIN MESSAGE"
    echo "Installing OpenBLAS into ${OPENBLAS_INSTALL_DIR} "
    echo "END MESSAGE"
    INSTALL_DIR=${OPENBLAS_INSTALL_DIR}
fi
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
OPENBLAS_DIR=${INSTALL_DIR}

unset LIBS
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi

echo "OpenBLAS: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "OpenBLAS: Unpacking archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/${TARNAME}.tar.gz

echo "OpenBLAS: Configuring..."
cd ${NAME}
# no configuration necessary

echo "OpenBLAS: Building..."
${MAKE} libs netlib shared MAKE="$MAKE" CC="$CC" FC="$F90" CFLAGS="$CFLAGS" FFLAGS="$F90FLAGS" LDFLAGS="$LDFLAGS" LIBS="$LIBS" INTERFACE64=${OPENBLAS_INT8}

echo "OpenBLAS: Installing..."
if [ "$(uname)" = "Darwin" ]; then
    # Create a script "install" that can handle the "-D" option that
    # OpenBLAS is using
    bindir="${BUILD_DIR}/bin"
    mkdir -p $bindir
    echo "exec ginstall" '"$@"' >$bindir/install
    chmod a+x $bindir/install
    export PATH="$bindir:$PATH"
fi
${MAKE} install MAKE="$MAKE" PREFIX="${INSTALL_DIR}" INTERFACE64=1
{
    cat <<EOT
#ifndef LAPACK_ILP64
#  define LAPACK_ILP64
#endif
EOT
    cat ${INSTALL_DIR}/include/lapacke_config.h
} > ${INSTALL_DIR}/include/lapacke_config.h.tmp
mv -f ${INSTALL_DIR}/include/lapacke_config.h.tmp ${INSTALL_DIR}/include/lapacke_config.h
{
    cat <<EOT
#ifndef HAVE_LAPACK_CONFIG_H
#  define HAVE_LAPACK_CONFIG_H
#endif
EOT
    cat ${INSTALL_DIR}/include/lapacke.h
} > ${INSTALL_DIR}/include/lapacke.h.tmp
mv -f ${INSTALL_DIR}/include/lapacke.h.tmp ${INSTALL_DIR}/include/lapacke.h
popd

echo "OpenBLAS: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "OpenBLAS: Done."
