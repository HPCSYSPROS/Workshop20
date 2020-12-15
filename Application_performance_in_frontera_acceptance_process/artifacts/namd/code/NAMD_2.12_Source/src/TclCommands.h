/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef TCLCOMMANDS_H
#define TCLCOMMANDS_H

#ifdef NAMD_TCL

#define USE_COMPAT_CONST
#include <tcl.h>

int tcl_vector_math_init(Tcl_Interp *);

#endif
#endif
