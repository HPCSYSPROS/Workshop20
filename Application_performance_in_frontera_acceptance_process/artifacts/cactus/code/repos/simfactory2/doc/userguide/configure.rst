
Configuring SimFactory
======================

Shell Alias
-----------

The main SimFactory binary is called "sim" and is located in
simfactory/bin.  You can execute SimFactory explicitly as
simfactory/bin/sim, but we recommend that you set up a shell alias in
your shell startup file so that you can just use the command "sim".
For bash users this file is .bashrc on Linux and .profile on Mac OS.
Add the following to the shell startup file:

  alias sim=simfactory/bin/sim

Setup
-----

SimFactory needs to be configured before it can be used.  The
recommended way to do this is to use the *setup* command which prompts
the user for the required information.

  sim setup

Follow the on-screen prompts.  This will output your choices in the
configuration file simfactory/etc/defs.local.ini.  It is likely that
you will have to further customise this file.  You can see some
possible option settings in simfactory/etc/defs.local.ini.example.

SimFactory needs to have a *machine definition* for every machine that
it is run on.  If you are using a machine that SimFactory already has
a definition for, such as a well-known supercomputer used by others in
the Cactus community, then no additional setup is required.  If,
however, you are running SimFactory on an individual laptop or
workstation, or on an unsupported supercomputer, the setup command
will also create a new machine definition for the local machine in
simfactory/mdb/machines.  You may also have to add extra information
to this file.
