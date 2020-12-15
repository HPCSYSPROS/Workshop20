VPIC
====

To simulate plasma, the Vector Particle-In-Cell (VPIC) code follows the movement
of charged particles in simulated electric and magnetic fields that interact
with the particles.  VPIC integrates the relativistic Maxwell-Boltzmann system
in a linear background medium for multiple particle species, in time with an
explicit-implicit mixture of velocity Verlet, leapfrog, Boris rotation and
exponential differencing based on a reversible phase-space volume conserving 2nd
order Trotter factorization.  VPIC can be used in studies of magnetic
reconnection of high temperature plasmas (H^+ and e^-).

For a detailed description of the VPIC benchmark, see here:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Output from VPIC is about 737k files and is too large to include here.  We have
provided output from stdout and stderr and a list of energies every 500
timesteps.

Code modifications
==================

We have made two modifications to VPIC:

1.  We changed all of the calls to MPI_Issend to MPI_Isend to work around a
potential bug in Intel's MPI stack.  This change does not change the computed
result;

2.  We modified the input file for the large case SPP benchmark to solve the
same problem on a different number of processes for a scaling study.  For
example, for the input file provided, the number of processes in the z-direction
is changed from 64 to 32.

Optimisations
=============

The only optimisations we make to the code are to use architecture-specific
compiler flags that generate code for the AVX512 instruction set.

Third-party library dependencies
================================

VPIC has no third-party library dependencies.

Timing
======

See the end of output/vpic_73728.401507
