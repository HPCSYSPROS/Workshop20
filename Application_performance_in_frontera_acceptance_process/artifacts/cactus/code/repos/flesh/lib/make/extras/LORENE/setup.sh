#! /bin/sh
#/*@@
#  @file   setup.sh
#  @date   Wed 03 Mar 2004
#  @author Erik Schnetter
#  @desc
#          Setup for an external Lorene installation
#  @enddesc
# @@*/

choose_lorene=`echo $LORENE | tr '[:upper:]' '[:lower:]'`

if [ "X$choose_lorene" = 'Xyes' ]; then

echo 'Configuring with Lorene'

# Work out Lorene's installation directory
if [ -z "$LORENE_DIR" ]; then
  echo '  Lorene selected but no LORENE_DIR set.  Checking some places...'
  CCTK_Search LORENE_DIR '/usr /usr/local /usr/local/lorene /usr/local/packages/lorene /usr/local/apps/lorene' Export/C++/Include/bin_bh.h
  if [ -z "$LORENE_DIR" ] ; then
     echo '  Unable to locate the Lorene directory - please set LORENE_DIR'
     exit 2
  fi
  echo "  Found a Lorene package in $LORENE_DIR"
else
  echo "  Using Lorene package in $LORENE_DIR"
fi


# Set platform-specific libraries
if [ -z "$LORENE_LIBS" ]; then
  LORENE_LIBS='lorene_export lorene lorenef77'
else
  echo "  Using Lorene libraries '$LORENE_LIBS'"
fi

# Set the Lorene libs, libdirs and includedirs
LORENE_LIB_DIRS='$(LORENE_DIR)/Lib $(LORENE_DIR)/Export/BinBH'
LORENE_INC_DIRS='$(LORENE_DIR)/Export/C++/Include'


# Write the data out to the header and make files.
CCTK_WriteLine cctk_Extradefs.h "#define CCTK_LORENE 1"

CCTK_WriteLine make.extra.defn "HAVE_LORENE     = 1"
CCTK_WriteLine make.extra.defn "LORENE_DIR      = $LORENE_DIR"
CCTK_WriteLine make.extra.defn "LORENE_LIBS     = $LORENE_LIBS"
CCTK_WriteLine make.extra.defn "LORENE_LIB_DIRS = $LORENE_LIB_DIRS"
CCTK_WriteLine make.extra.defn "LORENE_INC_DIRS = $LORENE_INC_DIRS"
CCTK_WriteLine make.extra.defn ''
CCTK_WriteLine make.extra.defn 'LIBS         += $(LORENE_LIBS)'
CCTK_WriteLine make.extra.defn 'LIBDIRS      += $(LORENE_LIB_DIRS)'
CCTK_WriteLine make.extra.defn 'SYS_INC_DIRS += $(LORENE_INC_DIRS)'

elif [ "X$choose_lorene" != 'Xno' -a "X$choose_lorene" != 'X' ]; then

  echo "  Don't understand the setting \"LORENE=$LORENE\" !"
  echo '  Please set it to either "yes" or "no", or leave it blank (same as "no") !'
  exit 1

fi
