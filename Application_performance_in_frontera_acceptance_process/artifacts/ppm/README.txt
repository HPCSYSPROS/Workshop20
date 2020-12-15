PPM
===

Piecewise Parabolic Method (PPM) for gas dynamical simulations. Website:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Code modifications
====================

The number of MPI tasks was reduced by a factor of 8. For this the number of
teams in x, y, and z direction were reduced by a factor of 2. The total work
load was kept constant by increasing the number of briquettes per task from
40x40x40 to 80x80x80.

This should keep the total grid size in briquettes the same as the original
problem, thus producing the same workload as the original problem.

Optimizations
=============

No optimizations were performed other than a general optimization flag
and a archticture flag: -xCORE-AVX512 -O3

Third-party library dependencies
================================

There are no third-party library dependencies.

Timing
======

Timing is at the end of output/run-JK_ompi.o393000
