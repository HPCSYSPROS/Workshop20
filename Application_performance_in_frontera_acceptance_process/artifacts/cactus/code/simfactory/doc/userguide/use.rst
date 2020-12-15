
Using SimFactory
================

.. toctree::
   :maxdepth: 2

Invoking SimFactory
-------------------

All SimFactory commands should be executed from your Cactus directory.
SimFactory commands are executed as follows:

  sim *command* [ *options* ]

The order of the *options* is not important.  The *Using SimFactory*
chapter describes how to use the different commands, and a full
listing of commands and their options is available in the *Reference*
chapter.

Building a configuration
------------------------

SimFactory provides the :ref:`command_build` command for building Cactus
configurations.  The most important argument to this command is the
:option:`--thornlist` option, as this tells Cactus which thorns from your source
tree you want to include in the configuration.  If you do not specify
a configuration name, a default name is used.  This default name is
"sim", though this can be customised in defs.local.ini.

Creating a simulation
---------------------

A simulation can be created using the :ref:`command_create` command.  SimFactory
needs to know a name for the simulation as well as what parameter file
to use.  You can either specify the name on the command line and give
the parameter file with the :option:`--parfile` option, or you can give the path
to the parameter file directly, and SimFactory will use the parameter
file name (without an extension) as the simulation name.  Instead of
specifying a parameter file, you can also specify :option:`--testsuite`, and
this instructs SimFactory to create a simulation which runs the Cactus
Test Suite for the configuration (see Running the Test Suite).

Running a simulation directly
-----------------------------

If you are working on a laptop or workstation, you can run SimFactory
simulations directly in your terminal without going via a queuing
system.  You can use the :ref:`command_run` command for this.  This command takes
the name of the simulation as its argument.  You probably also want to
specify the number of processors to use using :option:`--procs`.

Submitting a simulation
-----------------------

If you are running SimFactory on a supercomputer with a queuing
system, you cannot run simulations directly using the run command -
they must instead be *submitted* to the queuing system.  The :ref:`command_submit`
command takes the name of a simulation as an argument and submits a
job to the queuing system.  When the job starts, Cactus will run.  You
should also specify the number of processors using :option:`--procs` and the
amount of time that the simulation should run for using :option:`--walltime`.

Listing simulations
-------------------

You can use the :ref:`command_list-simulations` command to obtain a list of all the
simulations on the current machine.  The simulation will be listed as
either *active* or *inactive*, which indicates whether the simulation
is present in the queuing system (either in the queued or running
state).  If the simulation is currently running, the output of
list-simulations will also show this.

Showing output
--------------

The :ref:`command_show-output` command will show the output of a simulation if it has
started running.

Compound commands
-----------------

Simulations must always be created before they can be submitted or
run.  Since it is very common to want to create a simulation and
immediately submit or run it, SimFactory provides the :ref:`command_create-run` and
:ref:`command_create-submit` commands.  These commands create the simulation and then
either run or submit it immediately.

Remote commands
---------------

SimFactory is designed to be used both locally and remotely.  What
this means is that you can have a central source tree and installation
of SimFactory on your laptop or workstation, and *sync* to each remote
machine that you want to work on.  Development (and backup) would
happen only with the central source tree. This helps to ensure that
you don't make different modifications on different machines and get
confused about what is where.

SimFactory provides the :ref:`command_sync` command which takes a machine name as
an argument.  This command copies your source tree to the remote
machine.  It uses rsync internally which means that only files which
have changed are copied.

You can run SimFactory locally and have the commands executed on the
remote machine using the :option:`--remote` option.  All SimFactory commands
accept the :option:`--remote` option, which takes a machine name as argument.
The advantage of this is that you don't have to log in to each remote
machine, and you can write scripts on your central machine which
manage simulations on several machines.

Running the Test Suite
----------------------

Cactus thorns often contain regression tests consisting of test
parameter files and the resulting output, and Cactus has a mechanism
for verifying that the output of the parameter files with the current
version of the code matches the output stored in the thorn. SimFactory
can run these tests, using a queuing system if necessary. To run the
tests, you use the usual SimFactory commands for creating, submitting
or running a simulation, but you do not need to specify a parameter
file. Instead, you include the :option:`--testsuite` option, and if you want to
run specific tests, the :option:`--select-tests` option.

To run all tests immediately on two processors (cores):

 sim :ref:`command_create-run` mytests :option:`--testsuite` :option:`--procs` 2

where "mytests" is the name of the simulation that will be created.

To run the tests using a queuing system:

 sim :ref:`command_create-submit` mytests :option:`--testsuite` :option:`--procs` 2

You can use all the usual SimFactory commands and option for creating,
running and submitting simulations.  By default, the entire test suite
is run. If you want to run only specific tests, you can additionally
use the :option:`--select-tests` option. You can give this option a test name
(ending in .par), an arrangement name or a thorn specification in the
form <arrangement>/<thorn>.

 sim :ref:`command_create-run` mytests :option:`--testsuite` :option:`--procs` 2 :option:`--select-tests` McLachlan

 sim :ref:`command_create-run` mytests :option:`--testsuite` :option:`--procs` 2 :option:`--select-tests` McLachlan/ML_BSSN

 sim :ref:`command_create-run` mytests :option:`--testsuite` :option:`--procs` 2 :option:`--select-tests` ML_BSSN_sgw3d.par

Whether run using :ref:`command_create-run` or :ref:`command_create-submit`, a summary.log file will
be created in mytests/TEST/<config>/summary.log.  Note that for some
machines, you may need to use :option:`--ppn-used` to run on the correct number
of processors. Note that many of the tests will only run on 1 or 2
processes.

Since it is necessary to have the test data and Cactus flesh scripts
available when the job starts, the required data is copied into the
simulation restart directory on job submission (or interactive
running). This ensures that the test data and scripts are available
when the home directory is not mounted on the compute nodes, and that
the test data is not modified between job submission and job running.
On some machines, this copying process can take a long time.

When you use the :option:`--testsuite` option, it is not necessary to specify a
parameter file. The positional arguments syntax (parfile, cores,
walltime) is not supported for running the test suite.

Parameter file scripts
----------------------

It is often useful to specify simulations by higher-level descriptions
than Cactus parameter files. For example, when performing a
convergence test, many parameters might change between simulations at
different resolutions, and changing them all manually is tedious and
error-prone. Similarly, it can be very useful to set parameters to
values computed from simple expressions. For this reason, it is useful
to use a *parameter file script*, rather than a parameter file, as a
basic description of a simulation.

A parameter file script is a file with a ".rpar" extension which, when
executed, generates a file in the same place but with a ".par"
extension. The resulting file should be a valid SimFactory parameter
file. SimFactory supports such scripts directly.  You can use a script
in place of a parameter file when invoking SimFactory. When a
simulation is run, the script will be executed and the resulting
parameter file will be used by Cactus. NB: Remember to use the full path
of the original script when determining the output filename, not just the
base name.

For example,

  sim :ref:`command_create-submit` bbh :option:`--parfile` bbh.rpar :option:`--procs` 32 :option:`--walltime` 12:00:00

  sim :ref:`command_create-submit` bbh.rpar 32 12:00:00

The script is stored in the SIMFACTORY directory of the simulation, as
well as in the SIMFACTORY directories in each individual restart.

You can write a parameter file script in any language.  We provide
examples written in Python and in Perl:

Parameter file script written in Python::

  #!/usr/bin/env python

  import sys
  import re
  from string import Template

  ######################################################################

  dtfac = 0.5

  ######################################################################

  lines = """

  ActiveThorns = "Carpet CarpetIOBasic CarpetIOHDF5 CarpetLib CartGrid3D
                  CoordBase IOUtil SymBase Time"

  Cactus::cctk_itlast           = 10000000000000
  Cactus::max_runtime           = 0.1 # 6 seconds
  Cactus::terminate             = "runtime"

  Time::dtfac                   = $dtfac

  CarpetIOHDF5::checkpoint      = yes

  IO::checkpoint_on_terminate   = yes
  IO::recover                   = autoprobe

  IOBasic::outinfo_every        = 2000

  """

  data = open(re.sub(r'(.*)\.rpar$', r'\1.par', sys.argv[0]), 'w')
  data.write(Template(lines).substitute(locals()))

Parameter file script written in Perl::

  #!/usr/bin/perl -W
  
  $dtfac = 0.5;
  
  ######################################################################
  
  $lines = <<EOF;
  
  ActiveThorns = "Carpet CarpetIOBasic CarpetIOHDF5 CarpetLib CartGrid3D
                  CoordBase IOUtil SymBase Time"
  
  Cactus::cctk_itlast           = 10000000000000
  Cactus::max_runtime           = 0.1 # 6 seconds
  Cactus::terminate             = "runtime"
  
  Time::dtfac                   = $dtfac
  
  CarpetIOHDF5::checkpoint      = yes
  
  IO::checkpoint_on_terminate   = yes
  IO::recover                   = autoprobe
  
  IOBasic::outinfo_every        = 2000
  
  EOF
  
  $filename = "$0";
  $filename =~ s/\.rpar/.par/g; 
  
  open(OUT,">$filename");
  print OUT "$lines";
  close(OUT);

These parameter file scripts look like standard Cactus parameter files
but with $var variable replacements (in this case for the *dtfac*
variable), and Python/Perl headers and footers.  You can define new
variables and do calculations in the header and use the variables in
the main body.

Notes:

  * If you want to use the Cactus $parfile syntax, you need to escape the
    $.  In Python, this is done by repeating it:

      IO::out_dir = $$parfile

    and in Perl by using a backslash:

      IO::out_dir = \\$parfile

  * If you want to use SimFactory @...@ replacements, you can use these
    directly in Python but in Perl the @ signs must be escaped with a
    backslash:

      ManualTermination::max_walltime = \\@WALLTIME_HOURS\\@
