/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Measure.h"
#include "Node.h"
#include "Parameters.h"
#include "Molecule.h"

#ifdef NAMD_TCL

int Measure::wrapCommand(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {

/*
  double (wraped *)(Vector*, Molecule*, Parameters*);
  wraped = 
  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;

  double result = *wraped(coordinates,molecule,parameters);

*/
  return TCL_OK;
}

static int Tcl_centerOfNumber(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {

  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  // Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;
  int numAtoms = molecule->numAtoms;

  int number = 0;
  Vector center = 0;
  for( int i = 0; i < numAtoms; ++i ) {
    number += 1;
    center += coordinates[i];
  }
  center /= number;

  char s[1024];
  sprintf(s,"%g %g %g", center.x, center.y, center.z);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  return TCL_OK;
}

static int Tcl_centerOfMass(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {

  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  // Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;
  int numAtoms = molecule->numAtoms;

  Vector center = 0;
  BigReal totalMass = 0;
  for( int i = 0; i < numAtoms; ++i ) {
    BigReal mass = molecule->atommass(i);
    totalMass += mass;
    center += mass * coordinates[i];
  }
  center /= totalMass;

  char s[1024];
  sprintf(s,"%g %g %g", center.x, center.y, center.z);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  return TCL_OK;
}

static int Tcl_radiusOfGyration(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {

  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  // Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;
  int numAtoms = molecule->numAtoms;

  Vector center = 0;
  BigReal totalMass = 0;
  int i;
  for( i = 0; i < numAtoms; ++i ) {
    BigReal mass = molecule->atommass(i);
    totalMass += mass;
    center += mass * coordinates[i];
  }
  center /= totalMass;

  BigReal moment = 0;
  for( i = 0; i < numAtoms; ++i ) {
    BigReal mass = molecule->atommass(i);
    moment += mass * (coordinates[i] - center).length2();
  }
  BigReal radius = sqrt(moment/totalMass);

  char s[1024];
  sprintf(s,"%g", radius);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  return TCL_OK;
}

static int Tcl_loadCoords(ClientData, Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  
  if (objc < 2 || objc > 3) {
    Tcl_SetResult(interp,"loadCoords: wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const vname = objv[1];
  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  // Parameters *parameters = node->parameters;
  const Vector *coords = node->coords;
  int numAtoms = molecule->numAtoms;
 
  if (objc == 2) {
    // get all the coordinates
    for (int i=0; i<numAtoms; i++) {
      Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
      Tcl_Obj *arrkey = Tcl_NewIntObj(i+1);
      Tcl_IncrRefCount(arrkey);
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj(coords[i].x));
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj(coords[i].y));
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj(coords[i].z));
      if (!Tcl_ObjSetVar2(interp, vname, arrkey, newlist, 0))
        return TCL_ERROR;
      Tcl_DecrRefCount(arrkey);
    }
  } else {
    // third argument must be a list of indices
    int nelems;
    Tcl_Obj **elems;
    if (Tcl_ListObjGetElements(interp, objv[2], &nelems, &elems) != TCL_OK) {
      return TCL_ERROR;
    }
    for (int i=0; i<nelems; i++) {
      int ind; 
      if (Tcl_GetIntFromObj(interp, elems[i], &ind) != TCL_OK)
        return TCL_ERROR;
      Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
      Tcl_Obj *arrkey = Tcl_NewIntObj(ind);
      Tcl_IncrRefCount(arrkey);
      --ind;
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj(coords[ind].x));
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj(coords[ind].y));
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj(coords[ind].z));
      if (!Tcl_ObjSetVar2(interp, vname, arrkey, newlist, 0))
        return TCL_ERROR;
      Tcl_DecrRefCount(arrkey);
    }
  }
  return TCL_OK;
}

void Measure::createCommands(Tcl_Interp *interp) {
  Tcl_CreateCommand(interp, "centerOfNumber", Tcl_centerOfNumber,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "centerOfMass", Tcl_centerOfMass,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "radiusOfGyration", Tcl_radiusOfGyration,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "loadCoords", Tcl_loadCoords, 
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
}

void Measure::deleteCommands(Tcl_Interp *interp) {
  Tcl_DeleteCommand(interp, "centerOfNumber");
  Tcl_DeleteCommand(interp, "centerOfMass");
  Tcl_DeleteCommand(interp, "radiusOfGyration");
  Tcl_DeleteCommand(interp, "loadCoords");
}

#endif

