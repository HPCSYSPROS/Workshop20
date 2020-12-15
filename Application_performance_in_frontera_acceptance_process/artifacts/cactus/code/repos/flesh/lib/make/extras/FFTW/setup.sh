#! /bin/sh
#/*@@
#  @file   setup.sh
#  @date   Wed 03 Mar 2004
#  @author Erik Schnetter
#  @desc
#          Setup for an external FFTW installation
#  @enddesc
# @@*/

choose_fftw=`echo $FFTW | tr '[:upper:]' '[:lower:]'`

if [ "X$choose_fftw" = 'Xyes' ]; then

echo 'Configuring with FFTW'

# Work out FFTW's installation directory
if [ -z "$FFTW_DIR" ]; then
  echo '  FFTW selected but no FFTW_DIR set.  Checking some places...'
  CCTK_Search FFTW_DIR '/usr /usr/local /usr/local/fftw /usr/local/packages/fftw /usr/local/apps/fftw' include/rfftw.h
  if [ -z "$FFTW_DIR" ] ; then
     echo '  Unable to locate the FFTW directory - please set FFTW_DIR'
     exit 2
  fi
  echo "  Found a FFTW package in $FFTW_DIR"
else
  echo "  Using FFTW package in $FFTW_DIR"
fi


# Set platform-specific libraries
if [ -z "$FFTW_LIBS" ]; then
  if test -n "$MPI" -a "$MPI" != 'none' ; then
    FFTW_LIBS='rfftw_mpi fftw_mpi $(MPI_LIBS) rfftw fftw m'
  else
    FFTW_LIBS='rfftw fftw m'
  fi
else
  echo "  Using FFTW libraries '$FFTW_LIBS'"
fi

# Set the FFTW libs, libdirs and includedirs

# don't explicitely add standard include and library search paths
if [ "$FFTW_DIR" != '/usr' -a "$FFTW_DIR" != '/usr/local' ]; then
  FFTW_LIB_DIRS="$FFTW_DIR/lib"
  FFTW_INC_DIRS="$FFTW_DIR/include"
fi


# Write the data out to the header and make files.
CCTK_WriteLine cctk_Extradefs.h "#define CCTK_FFTW 1"

CCTK_WriteLine make.extra.defn "HAVE_FFTW     = 1"
CCTK_WriteLine make.extra.defn "FFTW_DIR      = $FFTW_DIR"
CCTK_WriteLine make.extra.defn "FFTW_LIBS     = $FFTW_LIBS"
CCTK_WriteLine make.extra.defn "FFTW_LIB_DIRS = $FFTW_LIB_DIRS"
CCTK_WriteLine make.extra.defn "FFTW_INC_DIRS = $FFTW_INC_DIRS"
CCTK_WriteLine make.extra.defn ''
CCTK_WriteLine make.extra.defn 'LIBS         += $(FFTW_LIBS)'
CCTK_WriteLine make.extra.defn 'LIBDIRS      += $(FFTW_LIB_DIRS)'
CCTK_WriteLine make.extra.defn 'SYS_INC_DIRS += $(FFTW_INC_DIRS)'

elif [ "X$choose_fftw" != 'Xno' -a "X$choose_fftw" != 'X' ]; then

  echo "  Don't understand the setting \"FFTW=$FFTW\" !"
  echo '  Please set it to either "yes" or "no", or leave it blank (same as "no") !'
  exit 1

fi
