==============================================================
  Description
==============================================================
The benchmark code MILC represents part of a set of codes written by the MIMD 
Lattice Computation (MILC) collaboration used to study quantum chromodynamics 
(QCD), the theory of the strong interactions of subatomic physics. It performs 
simulations of four dimensional SU(3) lattice gauge theory on MIMD parallel 
machines. "Strong interactions" are responsible for binding quarks into protons 
and neutrons and holding them all together in the atomic nucleus. QCD discretizes 
space and evaluates field variables on sites and links of a regular hypercube 
lattice in four-dimensional space time. Each link between nearest neighbors in 
this lattice is associated with a 3-dimensional SU(3) complex matrix for a given field.

The MILC collaboration has produced application codes to study several different 
QCD research areas, only one of which is used here. This code generates lattices 
using rational function approximations for the fermion determinants using the 
Rational Hybrid Monte Carlo (RHMC) algorithm and implementing the HISQ action. 

==============================================================
  How to compile
==============================================================
The directories with source code include:

  ./ks_imp_rhmc
  ./generic
  ./generic_ks
  ./libraries

In the top level directory a script is provided, "build.sh". This should be used, 
or followed, for a complete clean and build of the MILC code.    
 
The libraries subdirectory contains the makefile "Make_vanilla". This includes 
compiler options that affect the libraries directory only.  The compiler command 
used to create objects in this subdirectory can be a serial (i.e., not MPI) 
compiler. Modify Make_vanilla for your particular environment. The file 
Make_template should not need to be changed.

In the ks_imp_rhmc subdirectory the top level makefile is "Makefile".  
You can edit this file to change compiler options (e.g. modify the variables 
"CC" and "OPT") for your particular environment. This builds source in the 
directories ks_imp_rhmc, generic and generic_ks. The Make_template files in 
these directories should not need to be changed.

The PRECISION variable SHALL remain set to 2 for double precision. 

==============================================================
  Problem sizes and how to run
==============================================================
Input decks and sample batch scripts are provided in the "benchmarks" directory. 
The provided problem sizes include "small", "medium", "large", and "grand_challenge", 
each having their respective directory.

For each problem, a run script is provided. The scripts are set up for Slurm and 
NERSC's Edison platform. You will need to modify this for the scheduler used in 
your environment. The script sets up all the parameters required by "milc_in.sh". 
This includes the problem size and the parallel decomposition (e.g. node geometry 
and OpenMP threads). In general, these files can be modified to suit your environment 
and needs, EXCEPT any parameters associated with the problem definition, which 
includes the parameters nx, ny, nz, nt, warms and trajecs.

In order to ensure that a given lattice is in a reasonably equilibrated state it 
is necessary to run several warmup phases, which can take 10's of minutes to hours. 
This warmup is of no interest in benchmarking. Hence, lattices have been pre-generated 
for all problems sizes and can be found in the "lattices" directory. The provided run 
scripts are set up to load the appropriate lattice file.

The benchmarks directory also contains the input files to be used by all problems, 
"milc_in.sh" and "rationals_m001m05m5.test1". The milc_in.sh file SHALL NOT be modified 
EXCEPT for any logic setting up the parameter "command" used to launch a job. 
The rationals_m001m05m5.test1 file SHALL NOT be modified.

==============================================================
  How to create your own problem size (lattice)
==============================================================
The provided lattices were built with the respecitve "build_xxx.sh". 
If you would like to create a lattice of different dimensions for
you own purposes, simply copy one of the supplied directories and modify the the 
provided scripts. 

Of course, this is for development use only and not to be used in place of 
the problem sizes stated above. 

==============================================================
  How to report results
==============================================================
The entire standard output reporting for each problem size defined above
is to be provided.  The metric to be used in the calculation of SSI  is 
"TOTAL_TIME", in seconds, as reported towards the end of the output. 

==============================================================
  How to verify results
==============================================================
Verification is difficult, and there is no one value in the output that you can check 
for correctness of the result. If a variant of the code is made, or the parallel 
decomposition is changed for a given lattice size, there may be differences due to 
different orders of floating-point operations that will contribute to variances in the 
results. However, the final values of "PLAQUETTE ACTION" and "CHECK: delta S" should 
be close. You can check those against results found in the "benchmarks/sample_outputs"
directory.

The run rules require you submit the generated output file for benchmark submission. 
The Crossroads/NERSC-9 benchmark evaluation team will work the MILC developers if 
necessary to ensure a valid result.

==============================================================
  Performanc Analysis, Optimization  and Profiling 
==============================================================
The run time, and hence turn around time,  for a given problem may be excessive 
for performance anaylsis and code optimization tasks. A feature added to the 
production MILC code to support development is the ability to exit the simulation after a given 
number of iterations. This feature is enabled by compiling with PROF_MARK defined,
which is done by modifying ks_rhms_hisq/Makefile and adding -DPROF_MARK to the 
OCFLAGS variable. The early exit feature is then enabled by setting the environment
variable PROF_NUM_ITERS to a given value. For example:

  $ export PROF_NUM_ITERS=10000
  $ ./run_small.sh
  .
  .
  .
  ACTION: g,h,f = 3.31881361928098e+06  9.32699235594819e+02  1.57766672068415e+06  4.89741303920073e+06
  DG = 3.318814e+06, DH = 9.326992e+02, DF = 1.577667e+06, D = 4.897413e+06
  ------> PROFILING STARTED: start iters = 20848, stop iters >= 30848 iters
  step iters=9583: t_cg=3.372173e+01, t_force=9.760920e+00, t_link=1.039341e+00
  step iters=9479: t_cg=6.722683e+01, t_force=1.953130e+01, t_link=2.076335e+00
  ------> PROFILING STOPPED: stop iters = 39984 profiled iters = 19136
  TOTAL_TIME    171.639 secs
  PROFILED_TIME 94.137 secs
  CG_TIME       67.227 secs
  FORCE_TIME    19.531 secs
  LINK_TIME     2.076 secs
  Approximate memory usage = 3408.141 MiB
  Exiting due to profiling iteration count exceeded
  .
  .
  .

If PROF_NUM_ITERS is not defined, then the executable will run normally.

