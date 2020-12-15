#! /bin/sh
# /*@@
#   @file      setup.sh
#   @date      Sat May 31 17:56:34 CEST 2003
#   @author    Jonathan Thornburg, borrowing ++heavily from
#              Thomas Radke's LAPACK setup.sh
#   @desc
#              Setup for compilation with the BLAS library system installation
#              See http://www.netlib.org/blas/ for more info on the BLAS.
#   @enddesc
#   @version   $Header$
# @@*/

choose_blas=`echo $BLAS | tr '[:upper:]' '[:lower:]'`

if test "X$choose_blas" = 'Xyes' ; then

  echo 'Configuring with BLAS'

  # Search for BLAS installation
  if test -z "$BLAS_DIR"; then
    echo 'BLAS selected but no BLAS_DIR set... Checking some places'
    CCTK_Search BLAS_DIR '/usr/lib /usr/local/lib' libblas.a
    if test -z "$BLAS_DIR"; then
      CCTK_Search BLAS_DIR '/usr/lib /usr/local/lib' libblas.so
    fi
    if test -z "$BLAS_DIR"; then
      echo 'Unable to locate the BLAS library - please set BLAS_DIR'
      exit 2
    fi
    echo "Found a BLAS package in $BLAS_DIR"
    # don't explicitely add standard include and library search paths
    if [ "$BLAS_DIR" = '/usr/lib' -o "$BLAS_DIR" = '/usr/local/lib' ]; then
      BLAS_DIR=''
    fi
  elif test "$BLAS_DIR" = 'none'; then
    # user doesn't want the library path added
    BLAS_DIR=''
  fi

  if test -z "$BLAS_LIBS"; then
    BLAS_LIBS='blas'
  fi

  # write the variables out to the header and makefiles
  CCTK_WriteLine cctk_Extradefs.h '#define CCTK_BLAS 1'

  CCTK_WriteLine make.extra.defn "HAVE_BLAS     = 1"
  CCTK_WriteLine make.extra.defn "BLAS_LIBS     = $BLAS_LIBS $BLAS_EXTRA_LIBS m"
  CCTK_WriteLine make.extra.defn "BLAS_LIB_DIRS = $BLAS_DIR $BLAS_EXTRA_LIB_DIRS"
  CCTK_WriteLine make.extra.defn ''
  CCTK_WriteLine make.extra.defn 'LIBS           += $(BLAS_LIBS)'
  CCTK_WriteLine make.extra.defn 'LIBDIRS        += $(BLAS_LIB_DIRS)'

elif test "X$choose_blas" != 'Xno' -a "X$choose_blas" != 'X'; then

  echo "Don't understand the setting \"BLAS=$BLAS\" !"
  echo 'Please set it to either "yes" or "no", or leave it blank (same as "no") !'
  exit 1

fi
