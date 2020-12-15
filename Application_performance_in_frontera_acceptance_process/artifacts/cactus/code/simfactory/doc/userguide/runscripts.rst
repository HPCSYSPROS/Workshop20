
Run Scripts
===========

A run script is a shell script template which determines how Cactus
should be run on a particular machine.  It includes the necessary
mpirun variant for that machine, and the options that need to be
passed to it.  Run scripts make use of several variables defined by
SimFactory, such as the number of processors to run on, and the name
of the parameter file to run.  These definitions are refered to in the
run script using the @DEF@ syntax, where DEF is the name of the
definition.  SimFactory converts the templated run script into a real
shell script by expanding these quantities.

Run Script Definitions
----------------------

The definitions available to a run script are:

.. todo:: Include automatically generated list of definitions and their documentation from the code
