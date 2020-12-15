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
#include "ComputeFmmSerial.h"
#include "ComputeFmmSerialMgr.decl.h"
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


#ifdef FMM_SOLVER

// Calculate FMM:
//
// ufmmlap library from J. Huang - http://fastmultipole.org/Main/FMMSuite/
//   (ufmmlap = Uniform FMM Laplace Solver)
//
// The positions must translated and scaled to be within the unit
// cube centered at the origin, i.e., [-1/2,1/2]^3.  The potentials 
// and forces must then also be scaled.
//
extern "C" void fmmlap_uni_(
    int *natoms,        // input: number of particles
    double zat[][3],    // input: particle positions, length natoms
    double charge[],    // input: particle charges, length natoms
    double pot[],       // output: particle potentials, length natoms
    double field[][3],  // output: particle forces, length natoms
    int *nlev,          // input: number of oct-tree levels
    int *ier            // output: error code
    );

#endif // FMM_SOLVER


struct ComputeFmmSerialAtom {
  Position position;
  float charge;
  int id;
};

typedef Force FmmSerialForce;

class FmmSerialCoordMsg : public CMessage_FmmSerialCoordMsg {
public:
  int sourceNode;
  int numAtoms;
  Lattice lattice;
  ComputeFmmSerialAtom *coord;
};

class FmmSerialForceMsg : public CMessage_FmmSerialForceMsg {
public:
  BigReal energy;
  BigReal virial[3][3];
  FmmSerialForce *force;
};

class ComputeFmmSerialMgr : public CBase_ComputeFmmSerialMgr {
public:
  ComputeFmmSerialMgr();
  ~ComputeFmmSerialMgr();

  void setCompute(ComputeFmmSerial *c) { fmmCompute = c; }

  void recvCoord(FmmSerialCoordMsg *);
  void recvForce(FmmSerialForceMsg *);

private:
  CProxy_ComputeFmmSerialMgr fmmProxy;
  ComputeFmmSerial *fmmCompute;

  int numSources;
  int numArrived;
  FmmSerialCoordMsg **coordMsgs;
  int numAtoms;
  ComputeFmmSerialAtom *coord;
  FmmSerialForce *force;
  FmmSerialForceMsg *oldmsg;

  int fmmsolver;
  int nlevels;
  double scaling;
  Vector center;
  double *chargebuffer;
  double *epotbuffer;
  double (*posbuffer)[3];
  double (*forcebuffer)[3];
};

ComputeFmmSerialMgr::ComputeFmmSerialMgr() :
  fmmProxy(thisgroup), fmmCompute(0), numSources(0), numArrived(0),
  coordMsgs(0), coord(0), force(0), oldmsg(0), numAtoms(0),
  fmmsolver(0), nlevels(0), scaling(0), center(0),
  chargebuffer(0), epotbuffer(0), posbuffer(0), forcebuffer(0)
{
  CkpvAccess(BOCclass_group).computeFmmSerialMgr = thisgroup;
}

ComputeFmmSerialMgr::~ComputeFmmSerialMgr()
{
  for (int i=0;  i < numSources;  i++)  delete coordMsgs[i];
  delete [] coordMsgs;
  delete [] coord;
  delete [] force;
  delete oldmsg;
  if (chargebuffer) delete[] chargebuffer;
  if (epotbuffer)   delete[] epotbuffer;
  if (posbuffer)    delete[] posbuffer;
  if (forcebuffer)  delete[] forcebuffer;
}

ComputeFmmSerial::ComputeFmmSerial(ComputeID c) :
  ComputeHomePatches(c)
{
  CProxy_ComputeFmmSerialMgr::ckLocalBranch(
	CkpvAccess(BOCclass_group).computeFmmSerialMgr)->setCompute(this);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeFmmSerial::~ComputeFmmSerial()
{
}

void ComputeFmmSerial::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);

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

  // allocate message
  int numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  FmmSerialCoordMsg *msg = new (numLocalAtoms, 0) FmmSerialCoordMsg;
  msg->sourceNode = CkMyPe();
  msg->numAtoms = numLocalAtoms;
  msg->lattice = patchList[0].p->flags.lattice;
  ComputeFmmSerialAtom *data_ptr = msg->coord;

  // get positions
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    CompAtom *x = (*ap).positionBox->open();
    CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0;  i < numAtoms;  i++)
    {
      data_ptr->position = x[i].position;
      data_ptr->charge = x[i].charge;
      data_ptr->id = xExt[i].id;
      ++data_ptr;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  }

  CProxy_ComputeFmmSerialMgr fmmProxy(
      CkpvAccess(BOCclass_group).computeFmmSerialMgr);
  fmmProxy[0].recvCoord(msg);
}


void ComputeFmmSerialMgr::recvCoord(FmmSerialCoordMsg *msg) {
  int i;
  if ( ! numSources ) {
    numSources = (PatchMap::Object())->numNodesWithPatches();
    coordMsgs = new FmmSerialCoordMsg*[numSources];
    for ( i=0; i<numSources; ++i ) { coordMsgs[i] = 0; }
    numArrived = 0;
    numAtoms = Node::Object()->molecule->numAtoms;
    coord = new ComputeFmmSerialAtom[numAtoms];
    force = new FmmSerialForce[numAtoms];
  }

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

  double energy = 0;
  double virial[3][3] = { 0 };  // FMM does not calculate virial

  int rc = 0;  // return code

  if ( ! fmmsolver ) {
    //
    // setup FMM solver
    //

    // check boundary conditions
    if (lattice.a_p() || lattice.b_p() || lattice.c_p()) {
      NAMD_die("FMM solver requires non-periodic boundaries");
    }

    // setup number of levels
    if (simParams->FMMLevels > 0) {
      // explicitly set in config file
      nlevels = simParams->FMMLevels;
    }
    else {
      // otherwise estimate number of levels as round(log_8(numAtoms))
      nlevels = (int) floor(log((double)numAtoms) / log(8.) + 0.5);
    }

    // find bounding cube length
    if (numAtoms <= 0) {
      NAMD_die("setting up FMM with no atoms");
    }
    Vector min, max;
    min = max = coord[0].position;
    for (i = 1;  i < numAtoms;  i++) {
      Vector r = coord[i].position;
      if      (r.x < min.x)  min.x = r.x;
      else if (r.x > max.x)  max.x = r.x;
      if      (r.y < min.y)  min.y = r.y;
      else if (r.y > max.y)  max.y = r.y;
      if      (r.z < min.z)  min.z = r.z;
      else if (r.z > max.z)  max.z = r.z;
    }
    double padding = simParams->FMMPadding;
    if (padding <= 0) padding = 0.01;  // pad by at least a delta margin
    min -= padding;
    max += padding;
    double len = max.x - min.x;
    if (len < max.y - min.y)  len = max.y - min.y;
    if (len < max.z - min.z)  len = max.z - min.z;
    scaling = 1.0 / len;  // scale coordinates by length of the cube
    center = 0.5*(min + max);  // center of cube
    iout << iINFO << "FMM scaling length set to " << len << " A\n" << endi;
    iout << iINFO << "FMM center set to " << center << "\n" << endi;

    // allocate buffer space
    chargebuffer = new double[numAtoms];     // double *
    epotbuffer   = new double[numAtoms];     // double *
    posbuffer    = new double[numAtoms][3];  // double(*)[3]
    forcebuffer  = new double[numAtoms][3];  // double(*)[3]
    if (chargebuffer == 0 || epotbuffer == 0 ||
        posbuffer == 0 || forcebuffer == 0) {
      NAMD_die("can't allocate buffer space for FMM data");
    }

    // scale the charges - these won't change
    double celec = sqrt(COULOMB / simParams->dielectric);
    for (i = 0;  i < numAtoms;  i++) {
      chargebuffer[i] = celec * coord[i].charge;
    }
    fmmsolver = 1;  // is setup
  }

  // translate and scale positions into [-1/2,1/2]^3
  for (i = 0;  i < numAtoms;  i++) {
    posbuffer[i][0] = scaling*(coord[i].position.x - center.x);
    posbuffer[i][1] = scaling*(coord[i].position.y - center.y);
    posbuffer[i][2] = scaling*(coord[i].position.z - center.z);
  }

  // call the FMM solver
  int errcode = 0;
#ifdef FMM_SOLVER
  fmmlap_uni_(&numAtoms, posbuffer, chargebuffer, epotbuffer, forcebuffer,
      &nlevels, &errcode);
#else
  NAMD_die("Must link to FMM library to use FMM\n");
#endif

  // scale force and potentials
  for (i = 0;  i < numAtoms;  i++) {
    double qs = chargebuffer[i] * scaling;
    force[i].x = qs * scaling * forcebuffer[i][0];
    force[i].y = qs * scaling * forcebuffer[i][1];
    force[i].z = qs * scaling * forcebuffer[i][2];
    energy += qs * epotbuffer[i];
  }
  energy *= 0.5;

  // distribute forces
  for (int j=0;  j < numSources;  j++) {
    FmmSerialCoordMsg *cmsg = coordMsgs[j];
    coordMsgs[j] = 0;
    FmmSerialForceMsg *fmsg = new (cmsg->numAtoms, 0) FmmSerialForceMsg;

    for (int i=0;  i < cmsg->numAtoms;  i++) {
      fmsg->force[i] = force[cmsg->coord[i].id];
    }

    if ( ! j ) {  // set virial and energy only for first message
      fmsg->energy = energy;
      for (int k=0;  k < 3;  k++) {
        for (int l=0;  l < 3;  l++) {
          fmsg->virial[k][l] = virial[k][l];
        }
      }
    }
    else {  // set other messages to zero, add into reduction only once
      fmsg->energy = 0;
      for (int k=0;  k < 3;  k++) {
        for (int l=0;  l < 3;  l++) {
          fmsg->virial[k][l] = 0;
        }
      }
    }

    fmmProxy[cmsg->sourceNode].recvForce(fmsg);
    delete cmsg;
  }
}

void ComputeFmmSerialMgr::recvForce(FmmSerialForceMsg *msg) {
  fmmCompute->saveResults(msg);
  delete oldmsg;
  oldmsg = msg;
}

void ComputeFmmSerial::saveResults(FmmSerialForceMsg *msg)
{
  ResizeArrayIter<PatchElem> ap(patchList);

  FmmSerialForce *results_ptr = msg->force;

  // add in forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i) {
      f[i].x += results_ptr->x;
      f[i].y += results_ptr->y;
      f[i].z += results_ptr->z;
      ++results_ptr;
    }
  
    (*ap).forceBox->close(&r);
  }

  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += msg->energy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += msg->virial[0][0];
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += msg->virial[0][1];
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += msg->virial[0][2];
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += msg->virial[1][0];
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += msg->virial[1][1];
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += msg->virial[1][2];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += msg->virial[2][0];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += msg->virial[2][1];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += msg->virial[2][2];
  reduction->submit();
}

#include "ComputeFmmSerialMgr.def.h"
