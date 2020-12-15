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
#include "ComputeMsmSerial.h"
#include "ComputeMsmSerialMgr.decl.h"
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
#include "msm.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>


struct ComputeMsmSerialAtom {
  Position position;
  float charge;
  int id;
};

typedef Force MsmSerialForce;

class MsmSerialCoordMsg : public CMessage_MsmSerialCoordMsg {
public:
  int sourceNode;
  int numAtoms;
  Lattice lattice;
  ComputeMsmSerialAtom *coord;
};

class MsmSerialForceMsg : public CMessage_MsmSerialForceMsg {
public:
  BigReal energy;
  BigReal virial[3][3];
  MsmSerialForce *force;
};

class ComputeMsmSerialMgr : public CBase_ComputeMsmSerialMgr {
public:
  ComputeMsmSerialMgr();
  ~ComputeMsmSerialMgr();

  void setCompute(ComputeMsmSerial *c) { msmCompute = c; }

  void recvCoord(MsmSerialCoordMsg *);
  void recvForce(MsmSerialForceMsg *);

private:
  CProxy_ComputeMsmSerialMgr msmProxy;
  ComputeMsmSerial *msmCompute;

  int numSources;
  int numArrived;
  MsmSerialCoordMsg **coordMsgs;
  int numAtoms;
  ComputeMsmSerialAtom *coord;
  MsmSerialForce *force;
  MsmSerialForceMsg *oldmsg;

  NL_Msm *msmsolver;
  double *msmcoord;
  double *msmforce;
};

ComputeMsmSerialMgr::ComputeMsmSerialMgr() :
  msmProxy(thisgroup), msmCompute(0), numSources(0), numArrived(0),
  coordMsgs(0), coord(0), force(0), oldmsg(0), numAtoms(0),
  msmsolver(0), msmcoord(0), msmforce(0)
{
  CkpvAccess(BOCclass_group).computeMsmSerialMgr = thisgroup;
}

ComputeMsmSerialMgr::~ComputeMsmSerialMgr()
{
  for (int i=0;  i < numSources;  i++)  delete coordMsgs[i];
  delete [] coordMsgs;
  delete [] coord;
  delete [] force;
  delete oldmsg;
  if (msmsolver) NL_msm_destroy(msmsolver);
  if (msmcoord) delete[] msmcoord;
  if (msmforce) delete[] msmforce;
}

ComputeMsmSerial::ComputeMsmSerial(ComputeID c) :
  ComputeHomePatches(c)
{
  CProxy_ComputeMsmSerialMgr::ckLocalBranch(
	CkpvAccess(BOCclass_group).computeMsmSerialMgr)->setCompute(this);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeMsmSerial::~ComputeMsmSerial()
{
}

void ComputeMsmSerial::doWork()
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

  MsmSerialCoordMsg *msg = new (numLocalAtoms, 0) MsmSerialCoordMsg;
  msg->sourceNode = CkMyPe();
  msg->numAtoms = numLocalAtoms;
  msg->lattice = patchList[0].p->flags.lattice;
  ComputeMsmSerialAtom *data_ptr = msg->coord;

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

  CProxy_ComputeMsmSerialMgr msmProxy(
      CkpvAccess(BOCclass_group).computeMsmSerialMgr);
  msmProxy[0].recvCoord(msg);
}


static void rescale_nonperiodic_cell(
    Vector &u, Vector &v, Vector &w, Vector &c,
    Vector &ru, Vector &rv, Vector &rw,
    int isperiodic, int numatoms, const ComputeMsmSerialAtom *coord,
    double padding, double gridspacing) {
  const double NL_ZERO_TOLERANCE = 1e-4;
  double xlen = 1, inv_xlen = 1;  /* don't rescale periodic directions */
  double ylen = 1, inv_ylen = 1;
  double zlen = 1, inv_zlen = 1;
  double ulen_1 = u.rlength();
  double vlen_1 = v.rlength();
  double wlen_1 = w.rlength();
  double rpadx = padding * ulen_1;  /* padding distance in recip space */
  double rpady = padding * vlen_1;
  double rpadz = padding * wlen_1;
  double rhx = gridspacing * ulen_1;  /* grid spacing in recip space */
  double rhy = gridspacing * vlen_1;
  double rhz = gridspacing * wlen_1;
  Vector s, d;
  Vector rc = 0;  // don't move center along periodic directions
  double xmin, xmax, ymin, ymax, zmin, zmax;
  int is_periodic_x = (isperiodic & NL_MSM_PERIODIC_VEC1);
  int is_periodic_y = (isperiodic & NL_MSM_PERIODIC_VEC2);
  int is_periodic_z = (isperiodic & NL_MSM_PERIODIC_VEC3);
  int i;

  /* affine linear transformation of coordinates to reciprocal space,
   * where periodic vector directions map into [-0.5, 0.5) */
  //printf("*** center=%.4f %.4f %.4f\n", c.x, c.y, c.z);
  d = coord[0].position - c;
  s.x = ru * d;  // dot product
  s.y = rv * d;
  s.z = rw * d;
  xmin = xmax = s.x;
  ymin = ymax = s.y;
  zmin = zmax = s.z;
  for (i = 1;  i < numatoms;  i++) {
    d = coord[i].position - c;
    s.x = ru * d;  // dot product
    s.y = rv * d;
    s.z = rw * d;
    if      (s.x < xmin)  xmin = s.x;
    else if (s.x > xmax)  xmax = s.x;
    if      (s.y < ymin)  ymin = s.y;
    else if (s.y > ymax)  ymax = s.y;
    if      (s.z < zmin)  zmin = s.z;
    else if (s.z > zmax)  zmax = s.z;
  }
#if 0
  printf("*** xmin=%.4f  xmax=%.4f\n", xmin, xmax);
  printf("*** ymin=%.4f  ymax=%.4f\n", ymin, ymax);
  printf("*** zmin=%.4f  zmax=%.4f\n", zmin, zmax);
#endif

  if ( ! is_periodic_x) {
    xmax += rpadx;  /* pad the edges */
    xmin -= rpadx;
    if (rhx > 0) {  /* restrict center to rhx lattice points */
      double mupper = ceil(xmax / (2*rhx));
      double mlower = floor(xmin / (2*rhx));
      xmax = 2*rhx*mupper;
      xmin = 2*rhx*mlower;
    }
    rc.x = 0.5*(xmin + xmax);
    xlen = xmax - xmin;
    if (xlen < NL_ZERO_TOLERANCE) {
      xlen = 1;  /* leave scaling unchanged */
    }
    else {
      inv_xlen = 1/xlen;
    }
  }

  if ( ! is_periodic_y) {
    ymax += rpady;  /* pad the edges */
    ymin -= rpady;
    if (rhy > 0) {  /* restrict center to rhy lattice points */
      double mupper = ceil(ymax / (2*rhy));
      double mlower = floor(ymin / (2*rhy));
      ymax = 2*rhy*mupper;
      ymin = 2*rhy*mlower;
    }
    rc.y = 0.5*(ymin + ymax);
    ylen = ymax - ymin;
    if (ylen < NL_ZERO_TOLERANCE) {
      ylen = 1;  /* leave scaling unchanged */
    }
    else {
      inv_ylen = 1/ylen;
    }
  }

  if ( ! is_periodic_z) {
    zmax += rpadz;  /* pad the edges */
    zmin -= rpadz;
    if (rhz > 0) {  /* restrict center to rhz lattice points */
      double mupper = ceil(zmax / (2*rhz));
      double mlower = floor(zmin / (2*rhz));
      zmax = 2*rhz*mupper;
      zmin = 2*rhz*mlower;
    }
    rc.z = 0.5*(zmin + zmax);
    zlen = zmax - zmin;
    if (zlen < NL_ZERO_TOLERANCE) {
      zlen = 1;  /* leave scaling unchanged */
    }
    else {
      inv_zlen = 1/zlen;
    }
  }

#if 0
  printf("xmin=%g  xmax=%g\n", xmin, xmax);
  printf("ymin=%g  ymax=%g\n", ymin, ymax);
  printf("zmin=%g  zmax=%g\n", zmin, zmax);
  printf("xlen=%g  ylen=%g  zlen=%g\n", xlen, ylen, zlen);
  printf("rc_x=%g  rc_y=%g  rc_z=%g\n", rc.x, rc.y, rc.z);
#endif

  /* transform new center back to real space, then rescale basis vectors */
  c.x = (u.x*rc.x + v.x*rc.y + w.x*rc.z) + c.x;
  c.y = (u.y*rc.x + v.y*rc.y + w.y*rc.z) + c.y;
  c.z = (u.z*rc.x + v.z*rc.y + w.z*rc.z) + c.z;

#if 0
  printf("c_x=%g  c_y=%g  c_z=%g\n", c.x, c.y, c.z);
#endif

  u *= xlen;
  v *= ylen;
  w *= zlen;

  ru *= inv_xlen;
  rv *= inv_ylen;
  rw *= inv_zlen;
}


void ComputeMsmSerialMgr::recvCoord(MsmSerialCoordMsg *msg) {
  if ( ! numSources ) {
    numSources = (PatchMap::Object())->numNodesWithPatches();
    coordMsgs = new MsmSerialCoordMsg*[numSources];
    for ( int i=0; i<numSources; ++i ) { coordMsgs[i] = 0; }
    numArrived = 0;
    numAtoms = Node::Object()->molecule->numAtoms;
    coord = new ComputeMsmSerialAtom[numAtoms];
    force = new MsmSerialForce[numAtoms];
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

  double energy = 0;
  double virial[3][3];

  int rc = 0;  // return code

  if ( ! msmsolver ) {
    //
    // setup MSM solver
    //
    msmsolver = NL_msm_create();
    if ( ! msmsolver ) NAMD_die("unable to create MSM solver");
    double dielectric = simParams->dielectric;
    double cutoff = simParams->cutoff;
    double gridspacing = simParams->MSMGridSpacing;
    double padding = simParams->MSMPadding;
    int approx = simParams->MSMApprox;
    int split = simParams->MSMSplit;
    int nlevels = simParams->MSMLevels;
    int msmflags = 0;
    msmflags |= (lattice.a_p() ? NL_MSM_PERIODIC_VEC1 : 0);
    msmflags |= (lattice.b_p() ? NL_MSM_PERIODIC_VEC2 : 0);
    msmflags |= (lattice.c_p() ? NL_MSM_PERIODIC_VEC3 : 0);
    msmflags |= NL_MSM_COMPUTE_LONG_RANGE;  // compute only long-range part
    //msmflags |= NL_MSM_COMPUTE_ALL;
    //printf("msmflags = %x\n", msmflags);
    rc = NL_msm_configure(msmsolver, gridspacing, approx, split, nlevels);
    if (rc) NAMD_die("unable to configure MSM solver");
    Vector u=lattice.a(), v=lattice.b(), w=lattice.c(), c=lattice.origin();
    Vector ru=lattice.a_r(), rv=lattice.b_r(), rw=lattice.c_r();
    if ((msmflags & NL_MSM_PERIODIC_ALL) != NL_MSM_PERIODIC_ALL) {
      // called only if there is some non-periodic boundary
      int isperiodic = (msmflags & NL_MSM_PERIODIC_ALL);
      //printf("calling rescale\n");
      rescale_nonperiodic_cell(u, v, w, c, ru, rv, rw,
          isperiodic, numAtoms, coord, padding, gridspacing);
    }
    double vec1[3], vec2[3], vec3[3], center[3];
    vec1[0] = u.x;
    vec1[1] = u.y;
    vec1[2] = u.z;
    vec2[0] = v.x;
    vec2[1] = v.y;
    vec2[2] = v.z;
    vec3[0] = w.x;
    vec3[1] = w.y;
    vec3[2] = w.z;
    center[0] = c.x;
    center[1] = c.y;
    center[2] = c.z;
#if 0
    printf("dielectric = %g\n", dielectric);
    printf("vec1 = %g %g %g\n", vec1[0], vec1[1], vec1[2]);
    printf("vec2 = %g %g %g\n", vec2[0], vec2[1], vec2[2]);
    printf("vec3 = %g %g %g\n", vec3[0], vec3[1], vec3[2]);
    printf("center = %g %g %g\n", center[0], center[1], center[2]);
    printf("cutoff = %g\n", cutoff);
    printf("numatoms = %d\n", numAtoms);
#endif
    rc = NL_msm_setup(msmsolver, cutoff, vec1, vec2, vec3, center, msmflags);
    if (rc) NAMD_die("unable to set up MSM solver");
    msmcoord = new double[4*numAtoms];
    msmforce = new double[3*numAtoms];
    if (msmcoord==0 || msmforce==0) NAMD_die("can't allocate MSM atom buffers");
    // scale charges - these won't change
    double celec = sqrt(COULOMB / simParams->dielectric);
    for (i = 0;  i < numAtoms;  i++) {
      msmcoord[4*i+3] = celec * coord[i].charge;
    }
  }

  // evaluate long-range MSM forces
  for (i = 0;  i < numAtoms;  i++) {
    msmcoord[4*i  ] = coord[i].position.x;
    msmcoord[4*i+1] = coord[i].position.y;
    msmcoord[4*i+2] = coord[i].position.z;
  }
  for (i = 0;  i < numAtoms;  i++) {
    msmforce[3*i  ] = 0;
    msmforce[3*i+1] = 0;
    msmforce[3*i+2] = 0;
  }
  rc = NL_msm_compute_force(msmsolver, msmforce, &energy, msmcoord, numAtoms);
  if (rc) NAMD_die("error evaluating MSM forces");
  for (i = 0;  i < numAtoms;  i++) {
    force[i].x = msmforce[3*i  ];
    force[i].y = msmforce[3*i+1];
    force[i].z = msmforce[3*i+2];
  }
  // MSM does not yet calculate virial
  for (int k=0;  k < 3;  k++) {
    for (int l=0;  l < 3;  l++) {
      virial[k][l] = 0;
    }
  }

  // distribute forces
  for (int j=0;  j < numSources;  j++) {
    MsmSerialCoordMsg *cmsg = coordMsgs[j];
    coordMsgs[j] = 0;
    MsmSerialForceMsg *fmsg = new (cmsg->numAtoms, 0) MsmSerialForceMsg;

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

    msmProxy[cmsg->sourceNode].recvForce(fmsg);
    delete cmsg;
  }
}

void ComputeMsmSerialMgr::recvForce(MsmSerialForceMsg *msg) {
  msmCompute->saveResults(msg);
  delete oldmsg;
  oldmsg = msg;
}

void ComputeMsmSerial::saveResults(MsmSerialForceMsg *msg)
{
  ResizeArrayIter<PatchElem> ap(patchList);

  MsmSerialForce *results_ptr = msg->force;

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

#include "ComputeMsmSerialMgr.def.h"

