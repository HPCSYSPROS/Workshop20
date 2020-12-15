NWChem
======

NWChem provides scalable computational chemistry tools that are able to treat
large scientific computational chemistry problems, and efficiently use available
parallel computing resources from conventional workstation clusters to
high-performance parallel supercomputers. NWChem development community
represents a consortium of developers that is lead by the EMSL team located at
the Pacific Northwest National Laboratory (PNNL) in Washington State. The NWChem
development strategy focuses on delivering essential scientific capabilities to
its users in the areas of kinetics and dynamics of chemical transformations,
which undergo in gas phase, at interfaces, and in condensed phase.

The benchmark test involves an electronic structure calculation of Guanine dimer
at Coupled Cluster Singles and Doubles excitation including perturbative Triple
excitations known as CCSD(T) level of theory. The computation employs
correlation-consitent aug-cc-pvtz basis set.  The outcome of the computation is
the total energy of the molecular system at the fixed geometry representing a
stationary solution of the Schroedinger equation.

Find the detailed description of the case at:
https://bluewaters.ncsa.illinois.edu/spp-benchmarks

Input
=====

The benchmark includes two input files flops.nw.ccsd and flops.nw.t, which
constitute a composite job. Each of the input files utilizes 1024 SKX nodes on
Stampede2.  These jobs use a different number of cores per node, therefore they
have to be executed one after another within a single Slurm script. The total
run takes about 9263+7324 secs.  The first input file flops.nw.ccsd performs
iterative computation of single and double excitations, and saves the fully
converged amplitues in a binary restart file.  The second input file flops.nw.t
reads the previously saved amplitudes from the file and performs perturbative
calculation of triple excitations.

Output
======

NWChem produces some intermediate and restart files during a job.  These files
are not included here.

Code modifications
==================

Modifications to the source code involve only updating compiler flags to include
optimisations for the AVX512 instuction set.

A patch is provided that will modify the appropriate header file to make this
change for your perusal.  The build_common_512 script does everything needed to
build NWChem automatically and is provided for convenience.

Optimisations
=============

The only optimisations made to the code are during the build process to include
architecture-specific compiler flags for the AVX512 instruction set.

Third-party library dependencies
================================

There are no third-party library dependencies.

Timing
======

Timings are at the end of output/job.ccsd.406096.out and output/job.t.406096.out
