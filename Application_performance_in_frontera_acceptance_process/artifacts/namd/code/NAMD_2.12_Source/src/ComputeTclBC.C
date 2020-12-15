/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeTclBC.h"
#include "Node.h"
#include "SimParameters.h"
#include "Patch.h"
#include "Molecule.h"

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif
#include "TclCommands.h"

#define WRAPMODE_PATCH 0
#define WRAPMODE_INPUT 1
#define WRAPMODE_CELL 2
#define WRAPMODE_NEAREST 3

ComputeTclBC::ComputeTclBC(ComputeID c)
  : ComputeHomePatches(c), ap(patchList) {

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  SimParameters *simParams = Node::Object()->simParameters;

  wrapmode = WRAPMODE_PATCH;

  drops.resize(Node::Object()->molecule->numAtoms);
  cleardrops();

#ifdef NAMD_TCL
  if ( CkMyRank() && ! simParams->tclIsThreaded ) {
    NAMD_die("Sorry, tclBC requires TCL to be built with --enable-threads to use multiple threads per process.");
  }

  interp = Tcl_CreateInterp();
  tcl_vector_math_init(interp);

  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "wrapmode", Tcl_wrapmode,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "cleardrops", Tcl_cleardrops,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

  // run script to define calcforces, etc.
  if ( simParams->tclBCScript ) {
    int code = Tcl_Eval(interp,simParams->tclBCScript);
    const char *result = Tcl_GetStringResult(interp);
    if (result && *result != 0) CkPrintf("TCL: %s\n",result);
    if (code != TCL_OK) {
      const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
      NAMD_die(errorInfo ? errorInfo : "Unknown Tcl error");
    }
  } else NAMD_bug("tclBCScript pointer was NULL");

  // don't want these available until calcforces call
  Tcl_CreateObjCommand(interp, "dropatom", Tcl_dropatom,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "nextatom", Tcl_nextatom,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getcoord", Tcl_getcoord,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getcell", Tcl_getcell,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getmass", Tcl_getmass,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getcharge", Tcl_getcharge,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getid", Tcl_getid,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "addforce", Tcl_addforce,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "addenergy", Tcl_addenergy,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

#else

  NAMD_die("Sorry, tclBC is not available; built without TCL.");

#endif

}


ComputeTclBC::~ComputeTclBC() {
#ifdef NAMD_TCL
  Tcl_DeleteInterp(interp);
#endif
  delete reduction;
}


void ComputeTclBC::doWork() {

  SimParameters *simParams = Node::Object()->simParameters;
  lattice = &(patchList[0].p->lattice);
  const int step = patchList[0].p->flags.step;
  char cmd[128];

  energy = 0;
  n_atom = -1;  // set initial flags for iteration by nextatom

#ifdef NAMD_TCL
  sprintf(cmd,"calcforces %d %d %s",step,hasPatchZero,simParams->tclBCArgs);
  int code = Tcl_Eval(interp,cmd);
  if (code != TCL_OK) {
    const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo ? errorInfo : "Unknown Tcl error");
  }
  if (n_atom != -2) {
    NAMD_die("tclBCScript failed to call nextatom until failure");
  }
#endif

  reduction->item(REDUCTION_BC_ENERGY) += energy;
  reduction->submit();

}

#ifdef NAMD_TCL

int ComputeTclBC::Tcl_print(ClientData,
        Tcl_Interp *, int argc, char *argv[]) {
  Tcl_DString msg;
  Tcl_DStringInit(&msg);
  for ( int i = 1; i < argc; ++i ) {
    Tcl_DStringAppend(&msg," ",-1);
    Tcl_DStringAppend(&msg,argv[i],-1);
  }
  CkPrintf("TCL:%s\n",Tcl_DStringValue(&msg));
  Tcl_DStringFree(&msg);
  return TCL_OK;
}

int ComputeTclBC::Tcl_wrapmode(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp,"usage: wrapmode patch|input|cell|nearest",
								TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;

  if ( ! strcmp(argv[1],"patch") ) self->wrapmode = WRAPMODE_PATCH;
  else if ( ! strcmp(argv[1],"input") ) self->wrapmode = WRAPMODE_INPUT;
  else if ( ! strcmp(argv[1],"cell") ) self->wrapmode = WRAPMODE_CELL;
  else if ( ! strcmp(argv[1],"nearest") ) self->wrapmode = WRAPMODE_NEAREST;
  else {
    Tcl_SetResult(interp,"usage: wrapmode patch|input|cell|nearest",
								TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int ComputeTclBC::Tcl_cleardrops(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;

  self->cleardrops();

  return TCL_OK;
}

int ComputeTclBC::Tcl_dropatom(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  self->drops[self->fullatoms[self->i_atom].id] = 1;

  return TCL_OK;
}

int ComputeTclBC::Tcl_nextatom(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;

  // n_atom = -2 after all atoms processed
  if (self->n_atom < -1) {
    Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(0)));
    return TCL_OK;
  }

  // assume n_atom = -1 before first call
  do while ( self->n_atom < 0 || ++self->i_atom >= self->n_atom ) {
    if ( self->n_atom < 0 ) {  // first call
      self->ap = self->ap.begin();
    } else {
      (*(self->ap)).positionBox->close(&(self->atoms));
      (*(self->ap)).forceBox->close(&((*(self->ap)).r));
      self->ap++;
    }
    if ( self->ap == self->ap.end() ) {
      self->n_atom = -2;  // set error condition
      Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(0)));
      return TCL_OK;
    }
    self->i_atom = -1;
    self->n_atom = (*(self->ap)).p->getNumAtoms();
    self->fullatoms = (*(self->ap)).p->getAtomList().begin();
    self->atoms = (*(self->ap)).positionBox->open();
    (*(self->ap)).r = (*(self->ap)).forceBox->open();
    self->forces = (*(self->ap)).r->f[Results::normal];
  } while ( self->drops[self->fullatoms[self->i_atom].id] );

  Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(1)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_getcoord(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Position pos = self->atoms[i].position;
  switch ( self->wrapmode ) {
  case WRAPMODE_PATCH:
    break;
  case WRAPMODE_INPUT:
    pos = self->lattice->reverse_transform(pos,self->fullatoms[i].transform);
    break;
  case WRAPMODE_CELL:
    pos += self->lattice->wrap_delta(pos);
    break;
  case WRAPMODE_NEAREST:
    pos += self->lattice->wrap_nearest_delta(pos);
    break;
  }

  Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
  Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(pos.x)));
  Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(pos.y)));
  Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(pos.z)));
  Tcl_SetObjResult(interp, newlist);
  return TCL_OK;
}

int ComputeTclBC::Tcl_getcell(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;

  Tcl_Obj *newcell = Tcl_NewListObj(0, NULL);
  Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
  Vector o(self->lattice->origin());  // always have origin
  Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(o.x)));
  Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(o.y)));
  Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(o.z)));
  Tcl_ListObjAppendElement(interp, newcell, newlist);
  if (self->lattice->a_p()) {  // only if periodic
    newlist = Tcl_NewListObj(0, NULL);
    Vector a(self->lattice->a());
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(a.x)));
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(a.y)));
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(a.z)));
    Tcl_ListObjAppendElement(interp, newcell, newlist);
  }
  if (self->lattice->b_p()) {  // only if periodic
    newlist = Tcl_NewListObj(0, NULL);
    Vector b(self->lattice->b());
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(b.x)));
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(b.y)));
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(b.z)));
    Tcl_ListObjAppendElement(interp, newcell, newlist);
  }
  if (self->lattice->c_p()) {  // only if periodic
    newlist = Tcl_NewListObj(0, NULL);
    Vector c(self->lattice->c());
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(c.x)));
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(c.y)));
    Tcl_ListObjAppendElement(interp, newlist, Tcl_NewDoubleObj((double)(c.z)));
    Tcl_ListObjAppendElement(interp, newcell, newlist);
  }
  Tcl_SetObjResult(interp, newcell);
  return TCL_OK;
}

int ComputeTclBC::Tcl_getmass(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((double)(self->fullatoms[i].mass)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_getcharge(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((double)(self->atoms[i].charge)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_getid(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(self->fullatoms[i].id + 1)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_addforce(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  Tcl_Obj **force;  int fnum;  double x,y,z;
  if (Tcl_ListObjGetElements(interp, objv[1], &fnum, &force) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDoubleFromObj(interp, force[0],&x) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[1],&y) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"force not a vector",TCL_VOLATILE);
    return TCL_ERROR;
  }

  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int i = self->i_atom;
  self->forces[i].x += x;
  self->forces[i].y += y;
  self->forces[i].z += z;

  return TCL_OK;
}

int ComputeTclBC::Tcl_addenergy(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  double energy;
  if ( Tcl_GetDoubleFromObj(interp, objv[1], &energy) != TCL_OK ) {
    Tcl_SetResult(interp,"energy not a number",TCL_VOLATILE);
    return TCL_ERROR;
  }

  ComputeTclBC *self = (ComputeTclBC *)clientData;
  self->energy += energy;

  return TCL_OK;
}

#endif

