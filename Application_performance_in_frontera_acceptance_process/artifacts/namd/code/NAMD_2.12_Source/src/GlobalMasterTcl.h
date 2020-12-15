/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef COMPUTETCL_H
#define COMPUTETCL_H

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif

class GlobalMasterTcl : public GlobalMaster {
 public:
  GlobalMasterTcl();
  ~GlobalMasterTcl();
 protected:
  virtual void calculate();
 private:
  SubmitReduction *reduction;
  /* sets up the initial list of requested atoms */
  void initialize();
#ifdef NAMD_TCL
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_atomid(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_getstep(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addatom(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addgroup(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reconfig(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_clearconfig(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_loadcoords(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_loadmasses(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_loadforces(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_enabletotalforces(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_disabletotalforces(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_loadtotalforces(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_addforce(ClientData, Tcl_Interp *, int, Tcl_Obj * const []); 
  static int Tcl_addenergy(ClientData, Tcl_Interp *, int, char **);
#endif
};

#endif
