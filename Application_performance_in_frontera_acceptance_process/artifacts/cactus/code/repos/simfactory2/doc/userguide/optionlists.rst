
Option Lists
============

An option list is a file containing a set of variable definitions used
to configure a Cactus build on a particular machine.  For example, it
might contain the path to the C compiler, the desired optimisation
settings, and the location of the HDF5 library on that machine.

.. todo:: Include link to description of option lists on cactuscode.org

SimFactory allows Cactus option lists to contain variable definitions
which it substitutes.  These definitions are refered to in the
optionlist using the @DEF@ syntax, where DEF is the name of the
definition.  SimFactory converts the templated option list into a real
option list by expanding these quantities.

Option List Definitions
-----------------------

The definitions available to an option list are:

.. todo:: Include automatically generated list of definitions and their documentation from the code

.. todo:: Include guidelines for writing option lists for use with SimFactory
