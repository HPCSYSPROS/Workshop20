/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ScriptTcl.h"
#include <stdio.h>

#include "GlobalMaster.h"
#include "GlobalMasterTcl.h"

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"


#ifdef NAMD_TCL
int GlobalMasterTcl::Tcl_print(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  int arglen = 1;  int ai;
  for (ai=1; ai<argc; ++ai) { arglen += strlen(argv[ai]) + 1; }
  char *buf = new char[arglen];  *buf = 0;
  for (ai=1; ai<argc; ++ai) { strcat(buf,argv[ai]); strcat(buf," "); }
  ai = strlen(buf);  if ( ai ) buf[ai-1] = 0;
  CkPrintf("TCL: %s\n",buf);
  delete [] buf;
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_atomid(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 4) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  char *segid = argv[1];
  int resid;
  if (Tcl_GetInt(interp,argv[2],&resid) != TCL_OK) {
    return TCL_ERROR;
  }
  char *aname = argv[3];

  Molecule *mol = (Molecule *)clientData;
  int atomid = mol->get_atom_from_name(segid,resid,aname);

  if (atomid < 0) {
    Tcl_SetResult(interp,"atom not found",TCL_VOLATILE);
    return TCL_ERROR;
  }
  atomid += 1;

  char s[10];  sprintf(s,"%d",atomid);
  Tcl_SetResult(interp,s,TCL_VOLATILE);
  DebugM(4,"Atom ID " << atomid << " identified by name\n");
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_addatom(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  DebugM(2,"Tcl_addatom called\n");
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int atomid;
  if (Tcl_GetInt(interp,argv[1],&atomid) != TCL_OK) {
    return TCL_ERROR;
  }
  Molecule *mol = Node::Object()->molecule;
  int numAtoms = mol->numAtoms;
  if ( (atomid-1) < 0 || (atomid-1) >= numAtoms ) {
    char errmsg[128];
    sprintf(errmsg,"illegal atomid %d",atomid);
    Tcl_SetResult(interp,errmsg,TCL_VOLATILE);
    return TCL_ERROR;
  }
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  self->modifyRequestedAtoms().add(atomid-1);
  DebugM(4,"Atom ID " << atomid << " added to config list\n");
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_addgroup(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  DebugM(2,"Tcl_addgroup called\n");
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  ResizeArray<AtomIDList> &group_list = self->modifyRequestedGroups();

  /* set gcount to the number of groups after we add one, and add it! */
  int gcount = 1 + group_list.size();
  group_list.resize(gcount);

  /* get the list of atoms that go in the group */
  int listc, i;  char **listv;
  if (Tcl_SplitList(interp,argv[1],&listc,&listv) != TCL_OK) {
    return TCL_ERROR;
  }

  Molecule *mol = Node::Object()->molecule;
  int numAtoms = mol->numAtoms;
  /* add each atom to the new group */
  for ( i = 0; i < listc; ++i ) {
    int atomid;
    if (Tcl_GetInt(interp,listv[i],&atomid) != TCL_OK) { // error getting int
      group_list.resize(gcount-1); // remove the group we made
      Tcl_Free((char*) listv);
      return TCL_ERROR;
    }
    if ( (atomid-1) < 0 || (atomid-1) >= numAtoms ) {
      char errmsg[128];
      sprintf(errmsg,"illegal atomid %d",atomid);
      Tcl_SetResult(interp,errmsg,TCL_VOLATILE);
      return TCL_ERROR;
    }
    group_list[gcount-1].add(atomid-1); // add the atom to the group
  }
  Tcl_Free((char*) listv);

  /* return the group number to TCL */
  char s[10];  sprintf(s,"g%d",gcount);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  DebugM(4,"Group " << s << " added to config list\n");
  return TCL_OK;
}

/* this function is useless - it reconfigures whenever you add atoms! */
int GlobalMasterTcl::Tcl_reconfig(ClientData clientData,
	Tcl_Interp *interp, int argc, char **) {
  DebugM(2,"Tcl_reconfig called\n");
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  iout << iWARN << "'reconfig' is obsolete - reconfiguration is now automatic." << endi;
  iout << iWARN << "Use 'clearconfig' to clear the list of atoms and groups." << endi;
  DebugM(4,"Reconfiguration turned on\n");
  return TCL_OK;
}

int GlobalMasterTcl::Tcl_clearconfig(ClientData clientData,
	Tcl_Interp *interp, int argc, char **) {
  DebugM(2,"Tcl_reconfig called\n");
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  self->modifyRequestedGroups().resize(0);
  self->modifyRequestedAtoms().resize(0);
  return TCL_OK;
}

int GlobalMasterTcl::Tcl_getstep(ClientData clientData,
	Tcl_Interp *interp, int argc, char **) {
  DebugM(2,"Tcl_reconfig called\n");
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  
  char s[16];  sprintf(s,"%d",self->step);
  Tcl_SetResult(interp,s,TCL_VOLATILE);
  return TCL_OK;
}

int GlobalMasterTcl::Tcl_loadforces(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  DebugM(1,"Making tcl force array\n");
  if(objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const force_array_name = objv[1];
  
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  AtomIDList::const_iterator forced_ids_i = self->getLastAtomsForcedBegin();
  AtomIDList::const_iterator forced_ids_e = self->getLastAtomsForcedEnd();
  ForceList::const_iterator forces_i = self->getLastForcesBegin();

  // plf -- changed 06/12/2008 to check for more than one force on each atom

  // now make a Tcl array containing all of the requested atoms and
  // their forces
  DebugM(1,"Making Tcl array\n");
  while(forced_ids_i != forced_ids_e) {
    Tcl_Obj *array_key = Tcl_NewIntObj((int)((*forced_ids_i)+1)); // the id
    Tcl_IncrRefCount(array_key);

    // Check if the element is already defined, and if so, add to it
    Tcl_Obj *oldlist = Tcl_ObjGetVar2(interp, force_array_name, array_key, 0);
    Tcl_Obj *newlist = Tcl_NewListObj(0,NULL); // the list <fx,fy,fz>
    if (oldlist == NULL) {
      Tcl_ListObjAppendElement(interp, newlist,
        Tcl_NewDoubleObj((double)((*forces_i).x)));
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj((double)((*forces_i).y)));
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj((double)((*forces_i).z)));
    } else {
      Tcl_Obj** old_elems;
      int num_old_elems;
      double currval = 0.0;
      Tcl_ListObjGetElements(interp, oldlist, &num_old_elems, &old_elems);
      if (num_old_elems != 3) {
        NAMD_die("TCL error in loadforces! Force list doesn't have 3 elements!");
      }
      Tcl_GetDoubleFromObj(interp, old_elems[0], &currval);
      Tcl_ListObjAppendElement(interp, newlist,
        Tcl_NewDoubleObj((double)((*forces_i).x) + currval));
      Tcl_GetDoubleFromObj(interp, old_elems[1], &currval);
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj((double)((*forces_i).y + currval)));
      Tcl_GetDoubleFromObj(interp, old_elems[2], &currval);
      Tcl_ListObjAppendElement(interp, newlist, 
        Tcl_NewDoubleObj((double)((*forces_i).z + currval)));
    }

    // add the pair (id,F) to the array
    if (!Tcl_ObjSetVar2(interp, force_array_name, array_key, newlist, 0)) {
      NAMD_die("TCL error in loadforces!");
      return TCL_ERROR;
    }

    Tcl_DecrRefCount(array_key);

    // go to the next atom
    forced_ids_i++;
    forces_i++;
  }

  DebugM(1,"Done making tcl force array\n");
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_enabletotalforces(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[])
{
  DebugM(2,"Tcl_enabletotalforces called\n");
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  self->requestTotalForce(true);
  return TCL_OK;
}

int GlobalMasterTcl::Tcl_disabletotalforces(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[])
{
  DebugM(2,"Tcl_disabletotalforces called\n");
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  self->requestTotalForce(false);
  return TCL_OK;
}


// Here I simply copied the code from "Tcl_loadforces" above. The
// only difference is the source data.
int GlobalMasterTcl::Tcl_loadtotalforces(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[])
{
  if(objc != 2)
  { Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const force_array_name = objv[1];
  
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  if ( ! self->requestedTotalForces() ) {
    Tcl_SetResult(interp,"must call enabletotalforces before loadtotalforces",TCL_VOLATILE);
    return TCL_ERROR;
  }

  AtomIDList::const_iterator forced_ids_i = self->getForceIdBegin();
  AtomIDList::const_iterator forced_ids_e = self->getForceIdEnd();
  ForceList::const_iterator forces_i = self->getTotalForce();
  
  // now make a Tcl array containing all of the requested atoms and
  // their forces
  while(forced_ids_i != forced_ids_e) {
    Tcl_Obj *array_key = Tcl_NewIntObj((int)((*forced_ids_i)+1)); // the id
    Tcl_IncrRefCount(array_key);
    Tcl_Obj *newlist = Tcl_NewListObj(0,NULL); // the list <fx,fy,fz>
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*forces_i).x)));
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*forces_i).y)));
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*forces_i).z)));

    // add the pair (id,F) to the array
    if (!Tcl_ObjSetVar2(interp, force_array_name, array_key, newlist, 0)) {
      NAMD_die("TCL error in loadtotalforces!");
      return TCL_ERROR;
    }

    Tcl_DecrRefCount(array_key);

    // go to the next atom
    forced_ids_i++;
    forces_i++;
  }

  /* do the group stuff */
  ForceList::const_iterator tf_i = self->getGroupTotalForceBegin();
  ForceList::const_iterator tf_e = self->getGroupTotalForceEnd();
  int gcount = 1;
  for ( ; tf_i != tf_e; ++tf_i, ++gcount ) {
    Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
    char buf[10];
    sprintf(buf, "g%d", gcount);
    Tcl_Obj *arrkey = Tcl_NewStringObj(buf, -1);
    Tcl_IncrRefCount(arrkey);
 
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*tf_i).x)));
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*tf_i).y)));
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*tf_i).z)));
   
    if (!Tcl_ObjSetVar2(interp, force_array_name, arrkey, newlist, 0)) {
      NAMD_die("TCL error in loadtotalforces for groups!");
      return TCL_ERROR;
    }
    Tcl_DecrRefCount(arrkey);
  }
  return TCL_OK;
}

  
int GlobalMasterTcl::Tcl_loadcoords(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const vname = objv[1];
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  AtomIDList::const_iterator a_i = self->getAtomIdBegin();
  AtomIDList::const_iterator a_e = self->getAtomIdEnd();
  PositionList::const_iterator p_i = self->getAtomPositionBegin();
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
    Tcl_Obj *arrkey = Tcl_NewIntObj((int)((*a_i)+1));
    Tcl_IncrRefCount(arrkey);
    
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*p_i).x)));
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*p_i).y)));
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*p_i).z)));
   
    if (!Tcl_ObjSetVar2(interp, vname, arrkey, newlist, 0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
    Tcl_DecrRefCount(arrkey);
  }

  /* do the group stuff */
  PositionList::const_iterator c_i = self->getGroupPositionBegin();
  PositionList::const_iterator c_e = self->getGroupPositionEnd();
  int gcount = 1;
  for ( ; c_i != c_e; ++c_i, ++gcount ) {
    Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
    char buf[10];
    sprintf(buf, "g%d", gcount);
    Tcl_Obj *arrkey = Tcl_NewStringObj(buf, -1);
    Tcl_IncrRefCount(arrkey);
 
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*c_i).x)));
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*c_i).y)));
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*c_i).z)));
   
    if (!Tcl_ObjSetVar2(interp, vname, arrkey, newlist, 0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
    Tcl_DecrRefCount(arrkey);
  }
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_loadmasses(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const vname = objv[1];
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::const_iterator a_i = self->getAtomIdBegin();
  AtomIDList::const_iterator a_e = self->getAtomIdEnd();
  for ( ; a_i != a_e; ++a_i) {
    Tcl_Obj *arrkey = Tcl_NewIntObj((int)((*a_i)+1));
    Tcl_IncrRefCount(arrkey);
    if (!Tcl_ObjSetVar2(interp, vname, arrkey,
                        Tcl_NewDoubleObj((double)(mol->atommass(*a_i))),
                        0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
    Tcl_DecrRefCount(arrkey);
  }

  const BigReal *g_i, *g_e;
  g_i = self->getGroupMassBegin();
  g_e = self->getGroupMassEnd();
  int gcount = 1;
  for ( ; g_i != g_e; ++g_i, ++gcount) {
    char buf[10];
    sprintf(buf, "g%d", gcount);
    Tcl_Obj *arrkey = Tcl_NewStringObj(buf, -1);
    Tcl_IncrRefCount(arrkey);
    if (!Tcl_ObjSetVar2(interp, vname, arrkey,
                        Tcl_NewDoubleObj((double)(*g_i)),
                        0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
    Tcl_DecrRefCount(arrkey);
  }
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_addforce(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  DebugM(2,"Tcl_addforce called\n");
  if (objc != 3) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj **force;  int fnum;  int atomid;  double x, y, z;
  int isgroup = 0;
  char *id = Tcl_GetStringFromObj(objv[1], NULL); 
  if ( id[0] == 'g' ) {
    isgroup = 1;
    if ( Tcl_GetInt(interp,id+1,&atomid) != TCL_OK ) return TCL_ERROR;
  } else {
    if ( Tcl_GetInt(interp,id,&atomid) != TCL_OK ) return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[2], &fnum, &force) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDoubleFromObj(interp, force[0],&x) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[1],&y) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"force not a vector",TCL_VOLATILE);
    return TCL_ERROR;
  }

  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  if ( isgroup ) {
    int ngrps = self->getGroupMassEnd() - self->getGroupMassBegin();
    if ( atomid < 1 || atomid > ngrps ) {
      Tcl_SetResult(interp,"requested group not available",TCL_VOLATILE);
      return TCL_ERROR;
    }
    self->modifyGroupForces().item(atomid-1) += Vector(x,y,z);
  } else {
    self->modifyForcedAtoms().add(atomid-1);
    self->modifyAppliedForces().add(Vector(x,y,z));
  }
  DebugM(4,"Atom ID " << atomid << " added to force list\n");
  return TCL_OK;
}


int GlobalMasterTcl::Tcl_addenergy(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[])
{
  double energy;
  
  if (argc != 2)
    return TCL_ERROR;
  if (Tcl_GetDouble(interp,argv[1],&energy) != TCL_OK)
    return TCL_ERROR;
  
  GlobalMasterTcl *self = (GlobalMasterTcl *)clientData;
  self->reduction->item(REDUCTION_MISC_ENERGY) += energy;
  
  return TCL_OK;
}
#endif


GlobalMasterTcl::GlobalMasterTcl() {
  DebugM(3,"Constructing GlobalMasterTcl\n");
#ifdef NAMD_TCL
  interp = 0;
#endif
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  initialize();
  DebugM(2,"Done constructing ("<<requestedGroups().size()<<" initial groups)\n");
}

GlobalMasterTcl::~GlobalMasterTcl() {
  DebugM(3,"Destructing GlobalMasterTcl\n");
#ifdef NAMD_TCL
/*
  if ( interp ) Tcl_DeleteInterp(interp);
*/
#endif
  delete reduction;
}


void GlobalMasterTcl::initialize() {
  DebugM(4,"Initializing master\n");
#ifdef NAMD_TCL
  DebugM(1,"here\n");
  if(Node::Object() == NULL) NAMD_die("Node::Object() == NULL");
  if(Node::Object()->getScript() == NULL)
    NAMD_die("Node::Object()->getScript() == NULL");

  interp = Node::Object()->getScript()->interp;
  DebugM(1,"here\n");

  Tcl_CreateCommand(interp, "atomid", Tcl_atomid,
    (ClientData) (Node::Object()->molecule), (Tcl_CmdDeleteProc *) NULL);

  DebugM(1,"here\n");
  // Call interpreter to determine requested atoms
  Tcl_CreateCommand(interp, "addatom", Tcl_addatom,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "addgroup", Tcl_addgroup,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"enabletotalforces", Tcl_enabletotalforces,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"disabletotalforces", Tcl_disabletotalforces,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

  DebugM(1,"here\n");
  // Get the script
  StringList *script = Node::Object()->configList->find("tclForcesScript");

  DebugM(1,"here\n");
  for ( ; script; script = script->next ) {
    int code;
    DebugM(1,"here "<<script->data<<"\n");
    if ( strstr(script->data,"\n") ) {
       code = Tcl_Eval(interp,script->data);
    }
    else code = Tcl_EvalFile(interp,script->data);
    DebugM(1,"here\n");
    const char *result = Tcl_GetStringResult(interp);
    DebugM(1,"here\n");
    if (*result != 0) CkPrintf("TCL: %s\n",result);
    DebugM(1,"here\n");
    if (code != TCL_OK) {
      const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
      NAMD_die(errorInfo ? errorInfo : "Unknown Tcl error");
    }
  }

  DebugM(1,"here\n");
  Tcl_CreateObjCommand(interp, (char *)"loadforces", Tcl_loadforces,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"loadtotalforces", Tcl_loadtotalforces,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"loadcoords", Tcl_loadcoords,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"loadmasses", Tcl_loadmasses,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"addforce", Tcl_addforce,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"addenergy", Tcl_addenergy,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"reconfig", Tcl_reconfig,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"clearconfig", Tcl_clearconfig,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"getstep", Tcl_getstep,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
#else

  NAMD_die("Sorry, tclForces is not available; built without TCL.");

#endif
  DebugM(2,"done initializing master\n");
}


void GlobalMasterTcl::calculate() {
  DebugM(4,"Calculating forces on master\n");

  /* clear out the requested forces first! */
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);
  modifyGroupForces().resize(getGroupMassEnd() - getGroupMassBegin());
  modifyGroupForces().setall(Vector(0,0,0));

#ifdef NAMD_TCL
  // Call interpreter to calculate forces

  char cmd[129];  int code;
  strcpy(cmd,"calcforces");  code = Tcl_Eval(interp,cmd);
  const char *result = Tcl_GetStringResult(interp);
  if (*result != 0) CkPrintf("TCL: %s\n",result);
  if (code != TCL_OK) {
    const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo ? errorInfo : "Unknown Tcl error");
  }
#endif

  reduction->submit();

}
