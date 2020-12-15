/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeExt.h"
#include "ComputeExtMgr.decl.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif

struct ComputeExtAtom {
    Position position;
    float charge;
    int id;
};

class ExtCoordMsg : public CMessage_ExtCoordMsg {
public:
  int sourceNode;
  int numAtoms;
  Lattice lattice;
  ComputeExtAtom *coord;
};

class ExtForceMsg : public CMessage_ExtForceMsg {
public:
  BigReal energy;
  BigReal virial[3][3];
  ExtForce *force;
};

class ComputeExtMgr : public CBase_ComputeExtMgr {
public:
  ComputeExtMgr();
  ~ComputeExtMgr();

  void setCompute(ComputeExt *c) { extCompute = c; }

  void recvCoord(ExtCoordMsg *);
  void recvForce(ExtForceMsg *);

private:
  CProxy_ComputeExtMgr extProxy;
  ComputeExt *extCompute;

  int numSources;
  int numArrived;
  ExtCoordMsg **coordMsgs;
  int numAtoms;
  ComputeExtAtom *coord;
  ExtForce *force;
  ExtForceMsg *oldmsg;
};

ComputeExtMgr::ComputeExtMgr() :
  extProxy(thisgroup), extCompute(0), numSources(0), numArrived(0),
  coordMsgs(0), coord(0), force(0), oldmsg(0), numAtoms(0) {
  CkpvAccess(BOCclass_group).computeExtMgr = thisgroup;
}

ComputeExtMgr::~ComputeExtMgr() {
  for ( int i=0; i<numSources; ++i ) { delete coordMsgs[i]; }
  delete [] coordMsgs;
  delete [] coord;
  delete [] force;
  delete oldmsg;
}

ComputeExt::ComputeExt(ComputeID c) :
  ComputeHomePatches(c)
{
  CProxy_ComputeExtMgr::ckLocalBranch(
	CkpvAccess(BOCclass_group).computeExtMgr)->setCompute(this);

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}

ComputeExt::~ComputeExt()
{
}

void ComputeExt::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);

#if 0
  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics )
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
    }
    reduction->submit();
    return;
  }
#endif

  // allocate message
  int numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  ExtCoordMsg *msg = new (numLocalAtoms, 0) ExtCoordMsg;
  msg->sourceNode = CkMyPe();
  msg->numAtoms = numLocalAtoms;
  msg->lattice = patchList[0].p->flags.lattice;
  ComputeExtAtom *data_ptr = msg->coord;

  // get positions
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    CompAtom *x = (*ap).positionBox->open();
    CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
#if 0
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
#endif
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      data_ptr->position = x[i].position;
      data_ptr->charge = x[i].charge;
      data_ptr->id = xExt[i].id;
      ++data_ptr;
    }

#if 0
    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
#endif
    (*ap).positionBox->close(&x);
  }

  CProxy_ComputeExtMgr extProxy(CkpvAccess(BOCclass_group).computeExtMgr);
  extProxy[0].recvCoord(msg);

}

void ComputeExtMgr::recvCoord(ExtCoordMsg *msg) {
  if ( ! numSources ) {
    numSources = (PatchMap::Object())->numNodesWithPatches();
    coordMsgs = new ExtCoordMsg*[numSources];
    for ( int i=0; i<numSources; ++i ) { coordMsgs[i] = 0; }
    numArrived = 0;
    numAtoms = Node::Object()->molecule->numAtoms;
    coord = new ComputeExtAtom[numAtoms];
    force = new ExtForce[numAtoms];
  }

  int i;
  for ( i=0; i < msg->numAtoms; ++i ) {
    coord[msg->coord[i].id] = msg->coord[i];
  }

  coordMsgs[numArrived] = msg;
  ++numArrived;

  if ( numArrived < numSources ) return;
  numArrived = 0;

  // ALL DATA ARRIVED --- CALCULATE FORCES
  Lattice lattice = msg->lattice;
  SimParameters *simParams = Node::Object()->simParameters;
  FILE *file;
  int iret;

  // write coordinates to file
  //iout << "writing to file " << simParams->extCoordFilename << "\n" << endi;
  file = fopen(simParams->extCoordFilename,"w");
  if ( ! file ) { NAMD_die(strerror(errno)); }
  for ( i=0; i<numAtoms; ++i ) {
    int id = coord[i].id + 1;
    double charge = coord[i].charge;
    double x = coord[i].position.x;
    double y = coord[i].position.y;
    double z = coord[i].position.z;
    iret = fprintf(file,"%d %f %f %f %f\n",id,charge,x,y,z);
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
  }
  // write periodic cell lattice (0 0 0 if non-periodic)
  Vector a = lattice.a();  if ( ! lattice.a_p() ) a = Vector(0,0,0);
  Vector b = lattice.b();  if ( ! lattice.b_p() ) b = Vector(0,0,0);
  Vector c = lattice.c();  if ( ! lattice.c_p() ) c = Vector(0,0,0);
  iret = fprintf(file,"%f %f %f\n%f %f %f\n%f %f %f\n",
                 a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
  if ( iret < 0 ) { NAMD_die(strerror(errno)); }
  fclose(file);

  // run user-specified command
  //iout << "running command " << simParams->extForcesCommand << "\n" << endi;
  iret = system(simParams->extForcesCommand);
  if ( iret == -1 ) { NAMD_die(strerror(errno)); }
  if ( iret ) { NAMD_die("Error running command for external forces."); }

  // remove coordinate file
  iret = remove(simParams->extCoordFilename);
  if ( iret ) { NAMD_die(strerror(errno)); }

  // read forces from file (overwrite positions)
  //iout << "reading from file " << simParams->extForceFilename << "\n" << endi;
  file = fopen(simParams->extForceFilename,"r");
  if ( ! file ) { NAMD_die(strerror(errno)); }
  for ( i=0; i<numAtoms; ++i ) {
    int id, replace;
    double x, y, z;
    iret = fscanf(file,"%d %d %lf %lf %lf\n", &id, &replace, &x, &y, &z);
    if ( iret != 5 ) { NAMD_die("Error reading external forces file."); }
    if ( id != i + 1 ) { NAMD_die("Atom ID error in external forces file."); }
    force[i].force.x = x; force[i].force.y = y; force[i].force.z = z;
    force[i].replace = replace;
  }
  // read energy and virial if they are present
  // virial used by NAMD is -'ve of normal convention, so reverse it!
  // virial[i][j] in file should be sum of -1 * f_i * r_j
  double energy;
  double virial[3][3];
  iret = fscanf(file,"%lf\n", &energy);
  if ( iret != 1 ) {
    energy = 0;
    for ( int k=0; k<3; ++k ) for ( int l=0; l<3; ++l ) virial[k][l] = 0;
  } else {
    iret = fscanf(file,"%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n",
      &virial[0][0], &virial[0][1], &virial[0][2],
      &virial[1][0], &virial[1][1], &virial[1][2],
      &virial[2][0], &virial[2][1], &virial[2][2]);
    if ( iret != 9 ) {
      for ( int k=0; k<3; ++k ) for ( int l=0; l<3; ++l ) virial[k][l] = 0;
    } else {
      // virial used by NAMD is -'ve of normal convention, so reverse it!
      for ( int k=0; k<3; ++k ) for ( int l=0; l<3; ++l ) virial[k][l] *= -1.0;
    }
  }
  fclose(file);

  // remove force file
  iret = remove(simParams->extForceFilename);
  if ( iret ) { NAMD_die(strerror(errno)); }

  // distribute forces

  for ( int j=0; j < numSources; ++j ) {
    ExtCoordMsg *cmsg = coordMsgs[j];
    coordMsgs[j] = 0;
    ExtForceMsg *fmsg = new (cmsg->numAtoms, 0) ExtForceMsg;
    for ( int i=0; i < cmsg->numAtoms; ++i ) {
      fmsg->force[i] = force[cmsg->coord[i].id];
    }
    if ( ! j ) {
      fmsg->energy = energy;
      for ( int k=0; k<3; ++k ) for ( int l=0; l<3; ++l )
        fmsg->virial[k][l] = virial[k][l];
    } else {
      fmsg->energy = 0;
      for ( int k=0; k<3; ++k ) for ( int l=0; l<3; ++l )
        fmsg->virial[k][l] = 0;
    }
    extProxy[cmsg->sourceNode].recvForce(fmsg);
    delete cmsg;
  }

}

void ComputeExtMgr::recvForce(ExtForceMsg *msg) {
  extCompute->saveResults(msg);
  delete oldmsg;
  oldmsg = msg;
}

void ComputeExt::saveResults(ExtForceMsg *msg)
{
  ResizeArrayIter<PatchElem> ap(patchList);

  ExtForce *results_ptr = msg->force;

  // add in forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::normal];
    int numAtoms = (*ap).p->getNumAtoms();

    int replace = 0;
    ExtForce *replacementForces = results_ptr;
    for(int i=0; i<numAtoms; ++i) {
      if ( results_ptr->replace ) replace = 1;
      else f[i] += results_ptr->force;
      ++results_ptr;
    }
    if ( replace ) (*ap).p->replaceForces(replacementForces);
  
      (*ap).forceBox->close(&r);
    }

    reduction->item(REDUCTION_MISC_ENERGY) += msg->energy;
    reduction->item(REDUCTION_VIRIAL_NORMAL_XX) += msg->virial[0][0];
    reduction->item(REDUCTION_VIRIAL_NORMAL_XY) += msg->virial[0][1];
    reduction->item(REDUCTION_VIRIAL_NORMAL_XZ) += msg->virial[0][2];
    reduction->item(REDUCTION_VIRIAL_NORMAL_YX) += msg->virial[1][0];
    reduction->item(REDUCTION_VIRIAL_NORMAL_YY) += msg->virial[1][1];
    reduction->item(REDUCTION_VIRIAL_NORMAL_YZ) += msg->virial[1][2];
    reduction->item(REDUCTION_VIRIAL_NORMAL_ZX) += msg->virial[2][0];
    reduction->item(REDUCTION_VIRIAL_NORMAL_ZY) += msg->virial[2][1];
    reduction->item(REDUCTION_VIRIAL_NORMAL_ZZ) += msg->virial[2][2];
    reduction->submit();
}

#include "ComputeExtMgr.def.h"

