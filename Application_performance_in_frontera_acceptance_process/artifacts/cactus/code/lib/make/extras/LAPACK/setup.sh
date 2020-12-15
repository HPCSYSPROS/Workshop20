#! /bin/sh
# /*@@
#   @file      setup.sh
#   @date      Wed 8 January 2003
#   @author    Thomas Radke
#   @desc
#              Setup for compilation with the lapack library system installation
#   @enddesc
# @@*/

choose_lapack=`echo $LAPACK | tr '[:upper:]' '[:lower:]'`

if test "X$choose_lapack" = "Xyes" ; then

  echo "Configuring with LAPACK"

  # Search for LAPACK installation
  if test -z "$LAPACK_DIR"; then
    echo '  LAPACK selected but no LAPACK_DIR set... Checking some places'
    CCTK_Search LAPACK_DIR '/usr/lib /usr/local/lib' liblapack.a
    if test -z "$LAPACK_DIR"; then
      CCTK_Search LAPACK_DIR '/usr/lib /usr/local/lib' liblapack.so
    fi
    if test -z "$LAPACK_DIR"; then
      echo '  Unable to locate the LAPACK library - please set LAPACK_DIR'
      exit 2
    fi
    echo "  Found a LAPACK package in $LAPACK_DIR"
    # don't explicitely add standard include and library search paths
    if [ "$LAPACK_DIR" = '/usr/lib' -o "$LAPACK_DIR" = '/usr/local/lib' ]; then
      LAPACK_DIR=''
    fi
  elif test "$LAPACK_DIR" = 'none'; then
    # user doesn't want the library path added
    LAPACK_DIR=''
  fi

  if test -z "$LAPACK_LIBS"; then
    LAPACK_LIBS='lapack blas'
  fi

  # write the variables out to the header and makefiles
  CCTK_WriteLine cctk_Extradefs.h '#define CCTK_LAPACK 1'

  CCTK_WriteLine make.extra.defn "HAVE_LAPACK     = 1"
  CCTK_WriteLine make.extra.defn "LAPACK_LIBS     = $LAPACK_LIBS $LAPACK_EXTRA_LIBS m"
  CCTK_WriteLine make.extra.defn "LAPACK_LIB_DIRS = $LAPACK_DIR $LAPACK_EXTRA_LIB_DIRS"
  CCTK_WriteLine make.extra.defn ''
  CCTK_WriteLine make.extra.defn 'LIBS           += $(LAPACK_LIBS)'
  CCTK_WriteLine make.extra.defn 'LIBDIRS        += $(LAPACK_LIB_DIRS)'

elif test "X$choose_lapack" != "Xno" -a "X$choose_lapack" != "X"; then

  echo "  Don't understand the setting \"LAPACK=$LAPACK\" !"
  echo '  Please set it to either "yes" or "no", or leave it blank (same as "no") !'
  exit 1

fi
