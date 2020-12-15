WRF
===

The Weather Research and Forecasting (WRF) Model is a next-generation mesoscale
numerical weather prediction system designed for both atmospheric research and
operational forecasting applications. It features two dynamical cores, a data
assimilation system, and a software architecture supporting parallel computation
and system extensibility. The model serves a wide range of meteorological
applications across scales from tens of meters to thousands of kilometers.

For more detailed information, see:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Inputs
======

The WRF inputs are too large to include here.  They were used unmodified from
the SPP benchmark website: https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Code modification
=================

1. Nothing has been changed in the WRF source code, but we used a newer version
than the one provided on the SPP website (3.3.1).  This is because the version
on the SPP website does not scale with multithreading.

2. The namelist is modified to work with WRF 3.6.1 (in the Patches directory).

3. Extra vectorization flags are added to compile WRF 3.6.1 (in the Patches
   directory).

Optimizations
=============

No additional optimization is implemented except for the AVX512 instruction.

Third-party library dependencies
================================

netcdf 3.6.3

Timing
======

Timing for all steps can be found in output/rsl.error.0000.  The total time can
be found in output/time_N1680_n3360_T24.
