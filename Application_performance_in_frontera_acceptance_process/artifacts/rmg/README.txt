RMG
===

Description of the benchmark may be found here:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Note that not all job output is included here because it is over 3TB.

Code modifications
==================

Openbabel:
  1.) Disable override of CMAKE_CXX_FLAGS
  2.) Disable override of LD_LIBRARY_PATH
  3.) Disable broken check to see if GCC is newer than 4.0

RMG:
  1.) Enable the use of an old CMAKE policy (CMP0046) to continue building
  2.) Enable the use of ${MKLROOT} environment variable
  3.) Updated hard-coded Intel MKL link line
  4.) Updated hard-coded SCALAPACK link line
  5.) Added macro _BSD_SOURCE
  6.) Fixed RMG library install location
  7.) Disabled unnecessary poisson_pbc executable build

Patches are generated at build time from the build scripts themselves.  The
source code is downloaded from the internet, modified and a patch file is
generated. Currently, both openbabel and rmg require code/build harness
modifications to build properly. They are summarized by their respective patch
files.

For build instructions, see code/README.

Code optimizations
==================

The following GCC compiler flags were used:
-O3 -march=skylake-avx512 -mtune=skylake-avx512 -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi
With the addition of "-std=c++03" for RMG C++ specifically.

Third-party library dependencies
================================

Third party libraries for this version of RMG include:
1.) eigen/3.3.4.5a0156e40feb
  a.) fftw3/3.3.6
  b.) petsc/3.7
  c.) boost/1.64
  d.) superlu/5.2.1
2.) openbabel/2.3.2
  a.) eigen/3.3.4.5a0156e40feb
  b.) boost/1.64
3.) plplot/5.12.0
4.) rmg/2.2.1
  a.) boost/1.64
  b.) fftw3/3.3.6
  c.) intel mkl/17.0.4
  d.) openbabel/2.3.2
  e.) plplot/5.12.0

Notes
=====

From https://sourceforge.net/p/rmgdft/wiki/Building%20RMG/
"""
  Supported OS platforms
  Opensuse Linux 13.1.
  Cray XK/XE using the PrgEnv-gnu,boost,cudatoolkit,cmake and optionally magma modules.
  Microsoft windows using the Cygwin Posix emulation layer.
  Native Microsoft Windows using Visual Studio 2013, Intel Parallel Studio XE2015 and MSMPI.
  Mac with Macports
"""
The version provided by the SPP benchmark website is dependent on
libraries and headers found in particular on the BW Cray system. A more
modern official version of the code was used and modified to to run
on the CentOS Linux 7.3 distribution.

Timing
======

For timing information see output/in.dnv4096_3456xe_large.00.log
