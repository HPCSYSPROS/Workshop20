
Submit Scripts
==============

A submit script is a template script which determines how jobs should
be submitted on a particular machine.  For example, if the machine was
running the PBS queuing system, the submit script would be a "qsub"
script.

Submit scripts make use of several variables defined by SimFactory,
such as the number of nodes required for a job, and the walltime to
allocate.  These definitions are refered to in the submit script using
the @DEF@ syntax, where DEF is the name of the definition.  SimFactory
converts the templated submit script into a real submit script by
expanding these quantities.

Submit Script Definitions
-------------------------

The definitions available to a submit script are:

.. todo:: Include automatically generated list of definitions and their documentation from the code
