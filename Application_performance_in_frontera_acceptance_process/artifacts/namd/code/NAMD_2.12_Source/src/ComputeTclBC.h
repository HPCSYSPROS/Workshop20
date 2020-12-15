/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTETCLBC_H
#define COMPUTETCLBC_H

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif

#include "ComputeHomePatches.h"
#include "ReductionMgr.h"
#include "Tensor.h"
#ifndef WIN32
#include <strings.h>
#endif

class ComputeMgr;

class ComputeTclBC : public ComputeHomePatches {

public:
  ComputeTclBC(ComputeID c);
  virtual ~ComputeTclBC();
  void doWork();

private:
  int wrapmode;
  Lattice *lattice;

  ResizeArray<char> drops;
  void cleardrops() {
    memset((void*)drops.begin(), 0, drops.size()*sizeof(char));
  }

  ResizeArrayIter<PatchElem> ap;
  int i_atom, n_atom;
  CompAtom *atoms;
  FullAtom *fullatoms;
  Force *forces;
  BigReal energy;
  SubmitReduction *reduction;

#ifdef NAMD_TCL
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_wrapmode(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_cleardrops(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_dropatom(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_nextatom(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getcoord(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getcell(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getmass(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getcharge(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getid(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_addforce(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_addenergy(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
#endif

};

#endif







