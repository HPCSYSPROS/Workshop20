PSDNS
=====

PSDNS is a highly parallelized application code used for performing direct
numerical simulations (DNS) of three-dimensional unsteady turbulent fluid flows,
under the assumption of statistical homogeneity in space. A system of partial
differential equations expressing the fundamental laws of conservation of mass
and momentum is solved using Fourier pseudo-spectral methods in space and
explicit Runge-Kutta integration in time. The pseudo-spectralapproach requires
multiple FFTs to be taken in three directions per time step, resulting in
communication-intensive operations due to transposes involving collective
communication among processors. However, remote-memory addressing techniques
such as Co-Array Fortran on Blue Waters have been found to be helpful.

Simulation outputs can include flow field information stored at a large number
of grid points, as well as the trajectories of infinitesimal fluid elements
which are tracked over sustained periods of time for the study of turbulent
dispersion. The current code has enabled a production simulation of turbulent
flow using 4 trillion grid points on 262,144 CPU cores.

For the SPP benchmark, PSDNS is initialized with a simple sinusoidal velocity
field, and solves for the velocity field variables using the 4th order
Runge-Kutta method. Its IO and checkpoint are performed at the frequency
specified in the input file.

Code modifications
==================

None.

Optimisations
=============

Built with "-O3 -CORE-AVX2".

Third-party library dependencies
================================

Uses fftw3/3.3.6 and phdf5/1.8.16

Timing
======

We use the wallclock time reported at the end of the output from the batch
output file.  See output/psdns_core_inst_1640nodes_32x2048_1t.o392095
