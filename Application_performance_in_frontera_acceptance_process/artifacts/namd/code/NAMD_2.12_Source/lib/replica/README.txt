
This directory contains Tcl scripts that implement replica exchange.
This replaces the old Tcl server and socket connections driving a
separate NAMD process for every replica used in the simulation.

*** NET, IBVERBS, AND MULTICORE BUILDS ARE NOT SUPPORTED. ***

*** See EXPERIMENTAL MULTIPLEXED REPLICAS below for new option! ***

Charm++ 6.5.0 or newer for netlrts, verbs, mpi, gni, pamilrts, or other
machine layers based on the lrts low-level runtime implementation
are required for replica functionality based on Charm++ partitions.
A patched version of Charm++ is no longer required.

Replica exchanges and energies are recorded in the .history files
written in the output directories.  These can be viewed with, e.g.,
"xmgrace output/*/*.history" and processed via awk or other tools.
There is also a script to load the output into VMD and color each
frame according to replica index.  An example simulation folds
a 66-atom model of a deca-alanine helix in about 10 ns.

replica.namd - master script for replica exchange simulations
  to run: cd example; mkdir output
          (cd output; mkdir 0 1 2 3 4 5 6 7)
          charmrun +p8 namd2 +replicas 8 job0.conf +stdout output/%d/job0.%d.log
          charmrun +p8 namd2 +replicas 8 job1.conf +stdout output/%d/job1.%d.log

The number of NAMD processes must be a multiple of the number of replicas
(+replicas).  Be sure to increment jobX for +stdout option on command line.

replica_util.namd - functions for reduction, broadcast, and timing,
    useful for user-developed scripts and performance analysis

show_replicas.vmd - script for loading replicas into VMD, first source
    the replica exchange conf file and then this script, repeat for
    restart conf file or for example just do vmd -e load_all.vmd

clone_reps.vmd - provides "clone_reps" commmand to copy graphical
  representation from top molecule to all others

sortreplicas - found in namd2 binary directory, program to un-shuffle
  replica trajectories to place same-temperature frames in the same file.
  Usage: "sortreplicas <job_output_root> <num_replicas> <runs_per_frame> [final_step]"
  where <job_output_root> the job specific output base path, including
  %s or %d for separate directories as in output/%s/fold_alanin.job1
  Will be extended with .%d.dcd .%d.history for input files and
  .%d.sort.dcd .%d.sort.history for output files.  The optional final_step
  parameter will truncate all output files after the specified step,
  which is useful in dealing with restarts from runs that did not complete.
  Colvars trajectory files are similarly processed if they are found.

example subdirectory:

alanin_base.namd - basic config options for NAMD
alanin.params - parameters
alanin.psf - structure
unfolded.pdb - initial coordinates
alanin.pdb - folded structure for fitting in show_replicas.vmd

fold_alanin.conf - config file for replica_exchange.tcl script
job0.conf - config file to start alanin folding for 10 ns
job1.conf - config file to continue alanin folding another 10 ns

load_all.vmd - load all output into VMD and color by replica index


*** EXPERIMENTAL MULTIPLEXED REPLICAS - SUPPORTS MULTICORE ***

NAMD 2.11 provides new asynchronous inter-replica communication commands
that enable an arbitrary number of "logical" replicas to be multiplexed
onto a smaller number of "physical" replical partitions.  In short, the
+replicas argument no longer needs to match the num_replicas variable.
Multiplexing reduces the impact of performance differences between replicas
and avoids synchronizing writes to output files.

multiplex.namd - master script for multiplexed replica exchange simulations
  to run: cd example; mkdir output
          (cd output; mkdir 0 1 2 3 4 5 6 7)
          charmrun +p6 namd2 +replicas 6 mplexjob0.conf +stdout output/%d/job0.%d.log
          charmrun +p6 namd2 +replicas 6 mplexjob1.conf +stdout output/%d/job1.%d.log

Note that the stdout files will correspond to a single NAMD image onto
which multiple virtual replicas have been multiplexed, so they are not
as useful as for non-multiplexed replicas.  Dynamic stdout redirect
currently requires a file open for every run and does not work for
builds based on charmrun.  If you are using multicore or another build
that does not use charmrun and want to try it search for "stdout redirect"
in multiplex.namd and uncomment the two lines of code in the comment.


