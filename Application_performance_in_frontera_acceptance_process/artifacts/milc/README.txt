MILC
====

Code Modifications
==================

We used the MILC_2016DEC_BW_versD.tgz version of the benchmark although the
MILC_2017_BW_versE.tgz version is currently on the SPP website.  This is because
the version was updated after we had started working with versD and the source
code and problems are not different.  The only differences are in the job
submission scripts and flop counting tools, neither of which we use to run or
determine performance. We have made two modifications

1) We have modified the build system to use architecture specific flags
2) We have modifed the milc_in.sh file to be compatible with TACC's
MPI launcher ibrun.

Optimisations
=============

We have made no heavy-weight optimisations, but have added architecure-specific
compiler flags for the CORE-AVX2 instructions set and the ipo flag to enable
interprocedural optimisations.

Third-Party library dependencies
================================

The MILC benchmark requires only the default module set for Stampede2:

1. intel/17.0.4
2. impi/17.0.3
3. git/2.9.0
4. autotools/1.1
5. python/2.7.13
6. xalt/1.7.7
7. TACC

Timing
======

Timing output is in output/milc-large.359035.out.  Subtract UNIX_TIME_END and
UNIX_TIME_START for the time.
