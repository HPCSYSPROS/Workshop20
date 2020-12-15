NAMD
====

Scalable Molecular Dynamics software.

SPP benchmark of NAMD large case on Stampede2 Skylake nodes. Find the detailed
description of the case at: https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Inputs
======

Three large files:

1. chromat100-in.vel
2. chromat100-in.coor
3. final.js.inter.bin

are not provided but are used unmodified in the namd-spp-1.1.tar.gz file from
the namd directory in the shared benchmarks directory on the SPP website:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Code modifications
==================

The provided patch in the patches/ directory will turn off writing output files
(coordinates and velocities) at the end of the benchmark run.  Reported
SPP benchmark results do not include the final I/O at the end of a run, so this
patch does not affect the reported results.

Optimisations
=============

Architecture-specific compiler flags were used to enable use of the AVX512
instruction set on the Xeon Scalable processors.

Third-party library dependency
==============================

mkl 2016.3.210

Timing
======

Timing output is in, for example, output/log_2.12_mpi_long_1600_48_1_0.  Subtract
the timings at timesteps 69800 and 79800:

$ grep "TIMING: 69800" log_2.12_mpi_long_1600_48_1_0
TIMING: 69800  CPU: 657.666, 0.00960951/step  Wall: 657.666, 0.00960952/step, 0.027227 hours remaining, 1845.546875 MB of memory in use.
$ grep "TIMING: 79800" log_2.12_mpi_long_1600_48_1_0
TIMING: 79800  CPU: 753.787, 0.00969131/step  Wall: 753.787, 0.00969131/step, 0.000538406 hours remaining, 1845.546875 MB of memory in use.
