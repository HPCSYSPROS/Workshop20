#! /bin/sh
# /*@@
#  @file   setup.sh
#  @date   Fri Jul 30 1999
#  @author Thomas Radke
#  @desc
#          Setup HDF5
#  @enddesc
# @@*/

choose_hdf5=`echo $HDF5 | tr '[:upper:]' '[:lower:]'`

if [ "X$choose_hdf5" != 'Xyes' -a "X$choose_hdf5" != "Xno" -a "X$choose_hdf5" != "X" ]; then
  echo "  Don't understand the setting \"HDF5=$HDF5\" !"
  echo '  Please set it to either "yes" or "no", or leave it blank (same as "no") !'
  exit 1
fi

if test  "X$choose_hdf5" = 'Xyes' ; then 

  echo 'Configuring with HDF5.'

# Work out which variation of HDF5 is installed
# and set the HDF5 libs, libdirs and includedirs
  if [ -z "$HDF5_DIR" ]; then
    echo '  HDF5 selected but no HDF5_DIR set. Checking some places...'
    CCTK_Search HDF5_DIR '/usr /usr/local /usr/local/hdf5 /usr/local/packages/hdf5 /usr/local/apps/hdf5 /usr/local/hdf5/serial c:/packages/hdf5' include/hdf5.h
    if [ -z "$HDF5_DIR" ]; then
       echo '  Unable to locate the HDF5 directory - please set HDF5_DIR'
       exit 2
    fi
    echo "  Found an HDF5 package in $HDF5_DIR"

    # don't explicitely add standard include and library search paths
    if [ "$HDF5_DIR" != '/usr' -a "$HDF5_DIR" != '/usr/local' ]; then
      HDF5_LIB_DIRS="$HDF5_DIR/lib"
      HDF5_INC_DIRS="$HDF5_DIR/include"
    fi
  else
    echo "  Using HDF5 package in $HDF5_DIR"

    HDF5_LIB_DIRS="$HDF5_DIR/lib"
    HDF5_INC_DIRS="$HDF5_DIR/include"
  fi
  HDF5_LIBS=hdf5

# Check if HDF5 was built with Linux Large File Support (LFS)
  grep -qe _LARGEFILE_SOURCE ${HDF5_LIB_DIRS}/libhdf5.settings 2> /dev/null
  test_LFS=$?

# Check if the Stream VFD was compiled in
  grep -qe '#define H5_HAVE_STREAM 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
  test_stream_vfd=$?

# Check if the MPI I/O VFD was compiled in
  grep -qe '#define H5_HAVE_PARALLEL 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
  test_phdf5=$?

  if [ -n "$MPI" ]; then
    if [ $test_phdf5 -eq 0 ]; then
      echo '  Found parallel HDF5 library, so Cactus will potentially make use of parallel HDF5 support.'
#    else
#      echo '  Found serial HDF5 library, so Cactus can't make use of parallel HDF5 support.'
    fi
  else
    if [ $test_phdf5 -eq 0 ]; then
      echo "  Found parallel HDF5 library, but Cactus wasn't configured with MPI."
      echo '  Please set HDF5_DIR to point to a serial HDF5 package, or configure Cactus with MPI.'
      exit 2
    fi
  fi


# check that we have the right version of HDF5 under 32/64 bit IRIX
# This should better be checked by some autoconf script.
  if [ -n "$IRIX_BITS" ]; then
    if [ -r "$HDF5_LIB_DIRS/libhdf5.a" ]; then
      hdf5_lib="$HDF5_LIB_DIRS/libhdf5.a"
    elif [ -r "$HDF5_LIB_DIRS/libhdf5.so" ]; then
      hdf5_lib="$HDF5_LIB_DIRS/libhdf5.so"
    else
      hdf5_lib=
    fi

    if [ -n "$hdf5_lib" ]; then
      file $hdf5_lib | grep -qe $IRIX_BITS 2> /dev/null
      if [ $? -ne 0 ]; then
        echo "  The HDF5 library found in \"$HDF5_LIB_DIRS\" was not compiled as $IRIX_BITS bits !"
        echo '  Please reconfigure Cactus with the correct setting for HDF5_DIR !'
        exit 1
      fi
    fi
  fi


# check whether we run Windows or not
  $PERL -we 'exit (`uname` =~ /^CYGWIN/)'
  is_windows=$?


# check whether we run MacOS or not
  $PERL -we 'exit (`uname` =~ /^Darwin/)' 
  is_macos=$?



# Check whether we have to link with libsz.a
  grep -qe '#define H5_HAVE_LIBSZ 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
  test_szlib=$?
  if [ $test_szlib -eq 0 ]; then
    if [ $is_windows -ne 0 ]; then
      libsz='szlib.lib'
    elif [ $is_macos -ne 0 ]; then
      libsz='libsz.dylib'
    else
      libsz='libsz.a'
    fi

    if [ -z "$LIBSZ_DIR" -a ! -r /usr/lib/$libsz ]; then
      echo "  HDF5 library was built with external szlib I/O filter, searching for library $libsz ..."
      CCTK_Search LIBSZ_DIR "$HDF5_LIB_DIRS /usr/local/lib c:/packages/libsz/lib c:/packages/hdf5/lib" $libsz
      if [ -z "$LIBSZ_DIR" ]; then
        echo "  Unable to locate the library $libsz - please set LIBSZ_DIR"
        exit 2
      fi
      echo "  Found library $libsz in $LIBSZ_DIR"
    fi
    if [ $is_windows -eq 0 ]; then
      HDF5_LIBS="$HDF5_LIBS sz"
    else
      HDF5_LIBS="$HDF5_LIBS szlib"
    fi
    if [ ! -z "$LIBSZ_DIR" -a "$LIBSZ_DIR" != "$HDF5_LIB_DIRS" ]; then
      HDF5_LIB_DIRS="$HDF5_LIB_DIRS $LIBSZ_DIR"
    fi
  fi


# Check whether we have to link with libz.a
# this is for current versions of HDF5 (starting from 1.4.x)
  grep -qe '#define H5_HAVE_LIBZ 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
  test_zlib=$?

# this is for old versions of HDF5 (before 1.4.x)
  if [ $test_zlib -ne 0 ]; then
    grep -qe '#define HAVE_LIBZ 1' ${HDF5_DIR}/include/H5config.h 2> /dev/null
    test_zlib=$?
  fi

  if [ $test_zlib -eq 0 ]; then
    if [ $is_windows -ne 0 ]; then
      libz='zlib.lib'
    elif [ $is_macos -ne 0 ]; then
      libz='libz.dylib'
    else
      libz='libz.a'
    fi
    if [ -z "$LIBZ_DIR" -a ! -r /usr/lib/$libz -a ! -r /usr/lib64/$libz ]; then
      echo "  HDF5 library was built with external deflate I/O filter, searching for library $libz ..."
      CCTK_Search LIBZ_DIR '/usr/local/lib c:/packages/libz/lib c:/packages/hdf5/lib' $libz
      if [ -z "$LIBZ_DIR" ]; then
         echo "  Unable to locate the library $libz - please set LIBZ_DIR"
         exit 2
      fi
      echo "  Found library $libz in $LIBZ_DIR"
    fi
    if [ $is_windows -eq 0 ]; then
      HDF5_LIBS="$HDF5_LIBS z"
    else
      HDF5_LIBS="$HDF5_LIBS zlib"
    fi
    HDF5_LIB_DIRS="$HDF5_LIB_DIRS $LIBZ_DIR"
  fi

# Add MPI libraries for parallel HDF5
  if [ $test_phdf5 -eq 0 ]; then
    HDF5_LIBS="$HDF5_LIBS \$(MPI_LIBS)" 
    HDF5_INC_DIRS="$HDF5_INC_DIRS \$(MPI_INC_DIRS)"
    HDF5_LIB_DIRS="$HDF5_LIB_DIRS \$(MPI_LIB_DIRS)"
  fi

# Finally, add the math lib which might not be linked against by default
  if [ $is_windows -eq 0 ]; then
    HDF5_LIBS="$HDF5_LIBS m"
  fi

# Write the data out to the header and make files.
  CCTK_WriteLine cctk_Extradefs.h '#define CCTK_HDF5'
  CCTK_WriteLine make.extra.defn 'HAVE_HDF5            = 1'
  if [ $test_stream_vfd -eq 0 ]; then
    CCTK_WriteLine make.extra.defn 'HAVE_HDF5_STREAM_VFD = 1'
  fi
  if [ $test_LFS -eq 0 ]; then
    CCTK_WriteLine make.extra.defn 'HDF5_LFS_FLAGS       = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64'
  fi
  CCTK_WriteLine make.extra.defn "HDF5_LIBS            = $HDF5_LIBS"
  CCTK_WriteLine make.extra.defn "HDF5_LIB_DIRS        = $HDF5_LIB_DIRS"
  CCTK_WriteLine make.extra.defn "HDF5_INC_DIRS        = $HDF5_INC_DIRS"

  CCTK_WriteLine make.extra.defn ''

  CCTK_WriteLine make.extra.defn 'LIBS         += $(HDF5_LIBS)'
  CCTK_WriteLine make.extra.defn 'LIBDIRS      += $(HDF5_LIB_DIRS)'
  CCTK_WriteLine make.extra.defn 'SYS_INC_DIRS += $(HDF5_INC_DIRS)'

fi
