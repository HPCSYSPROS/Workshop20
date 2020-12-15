/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Molecule.C is compiled twice!
// MOLECULE2_C undefined only compiles first half of file
// MOLECULE2_C defined only compiles second half of file
// This is shameful but it works.  Molecule needs refactoring badly.

#define MOLECULE2_C
#include "Molecule.C"

