/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Measurements on coordinates of entire system during run.
*/

#ifndef MEASURE_H
#define MEASURE_H

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>

class Measure {
public:
  static void createCommands(Tcl_Interp *);
  static void deleteCommands(Tcl_Interp *);
private:
  static int wrapCommand(ClientData, Tcl_Interp*, int, char**);
};

#endif

#endif

