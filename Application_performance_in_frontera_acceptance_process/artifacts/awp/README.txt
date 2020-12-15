AWP-ODC
=======

For a description of the benchmark see here:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Inputs
======

Inputs are used unmodified (with one exception, described below) from the SPP
website.

Code modifications
==================

The only modification made to AWP-ODC was to the IN3D input file.  We modified
lines 19--21:

32 NPX           = number of procs in the x direction
32 NPY           = number of procs in the Y direction
64 NPZ           = number of procs in the Z direction

These numbers reflect the number of MPI processes used during the run, and we
modified them to perform scaling studies.

Optimisations
=============

The only optimisations done were the use of the following compiler flags:

-O3 -ipo -ip -inline-forceinline -xCORE-AVX512 -mp1 -c -extend_source -assume byterecl

Third-party library dependencies
================================

There are no third-party library dependencies.

Timing
======

Timing output is in output/awp_large.394085.out
