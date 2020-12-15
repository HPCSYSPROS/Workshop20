QMCPACK
=======

QMCPACK is a Quantum Monte Carlo (QMC) code used for many-body ab initio
simulation of atoms, molecules, and solid materials.  Open-source and written in
C++, QMCPACK solves the many-body Schrödinger equation by stochastically
sampling the configuration domain.  A Variational Monte Carlo (VMC) algorithm
can be used either to obtain the results directly or to quickly find a
“ballpark” estimate of the solution, which is then refined by a Diffusion Monte
Carlo (DMC) algorithm.  In both configurations, a number of VMC walkers (each a
complete state representation) randomly move through the energy domain each time
step.  If the VMC-DMC setup is used, the VMC walkers are sampled to create
walkers for the DMC phase.  The output is the lowest energy state within a
statistical uncertainty, which can be reduced by taking more samples (i.e.,
using more walkers).

The 4x4x1 graphite problem consists of 4x4x1 blocks of carbon atom supercells,
each containing four carbon atoms (64 total per 4x4x1 block), which are repeated
in three dimensions (i.e., the boundary conditions are periodic).  With a total
of 256 valence electrons per 4x4x1 block (four per atom), this problem
represents two graphene layers stacked in the third dimension.  The goal is to
find the lowest energy state that describes the system.  Both VMC and DMC
algorithms are used to solve this problem.

Inputs
======

The input file lda.pwscf.h5 is too large to include here.  It is used
unmodified from the qmcpack_dev_20170512-1310.tgz tarball on the SPP website
here: https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Code modifications
==================

None.

Optimisations
=============

-O3 -fopenmp -xCOMMON-AVX512

Third-party library dependencies
================================

1. boost 1.64
2. hdf5 1.8.16
3. mkl 2017.4.196
4. libxml2 2.9.3
5. zlib

Timing
======

See the end of output/myjob.o384573
