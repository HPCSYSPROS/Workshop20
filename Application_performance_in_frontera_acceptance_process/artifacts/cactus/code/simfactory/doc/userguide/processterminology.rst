
Terminology for nodes, cores, processes, and threads
===============================================================

SimFactory uses some terminology for nodes, cores, processes, and
threads. These terms are used somewhat inconsistently in the wild, and
marketing changes their definitions every few years. Below we define
the terms unambiguously. It is unfortunate that SimFactory's variable
names and command line options are somewhat outdated and don't
correspond to modern terminology any more.

Definitions
-----------

- A *machine* consists of a certain number of *nodes*, each of which
  consists of a certain number of *cores*.

- A job requests (from the queuing system) a certain number of nodes,
  and requests a certain number of cores on each node.

- The MDB entries [maxnodes], [minppn], and [maxppn] define how many
  nodes a job can request at most, and thus define the range of
  allowed values for cores-per-node.

- SimFactory starts a number of MPI *processes*, choosing how many MPI
  processes should be placed on every node. Each MPI process starts a
  certain number of OpenMP *threads*. The distribution of threads onto
  cores is performed automatically by the operating system and usually
  cannot be influenced.

Note that nodes and cores are requested from the queuing system, while
processes and threads are started by SimFactory. These numbers may
differ, allowing under- and over-subscription.

Variables and command-line options
----------------------------------

SimFactory expands variables in files and uses certain command line
options. These variables and options have the following definitions:

===============   =============   ===================
Variable          Option          Definition
===============   =============   ===================
NODES                             nodes
PROCS_REQUESTED                   cores
PPN               --ppn           cores per node
NUM_PROCS                         processes
NODE_PROCS                        processes per node
PROCS             --procs         threads
NUM_THREADS       --num-threads   threads per process
PPN_USED          --ppn-used      threads per node
NUM_SMT           --num-smt       threads per core
===============   =============   ===================

Choice process
--------------

The user chooses the total number of threads (--procs). The user can
also choose the number of threads per process (--num-threads) and the
number of threads per core (--num-smt). Additionally, the user can
also specify the number of cores per node (--ppn) and the number of
threads per node (--ppn-used), allowing for under- or over-subscribing
or cores. The number of nodes is always chosen automatically. Values
that are not specified are taken from a previous restart (if one
exists), or from the system's MDB entry.

The number of cores per node that can be requested from the queuing
system define a granularity that may be inconsistent with the total
number of threads. In this case, the last node may be used only
partially. Similarly, the number of threads per process defines a
granularity that may be inconsistent with the total number of threads.
In this case, the total number of threads is rounded up, so that the
job will have more threads running.

Definitions and constraints
---------------------------

Number of MPI processes:
   NUM_PROCS := PROCS / NUM_THREADS
constraint:
   PROCS % NUM_THREADS = 0

Number of nodes:
   NODES := ceil(PROCS / PPNUSED)

Number of requested cores:
   PROCS_REQUESTED := NODES * PPN

Number of MPI processes per node:
   NODE_PROCS := PPNUSED * NUM_SMT/ NUM_THREADS
constraint:
   PPNUSED * NUM_SMT % NUM_THREADS = 0
