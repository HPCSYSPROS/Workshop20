/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifdef CHARM_HAS_MSA

#include "InfoStream.h"
#include "Node.h"
#include "PDB.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeMsmMsa.h"
#include "ComputeMsmMsaMgr.decl.h"
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
#include "pup_stl.h"
#include "MsmMacros.h"


void MsmMsaData::pup(PUP::er &p)
{
  p|ispx, p|ispy, p|ispz;
  p|hx_1, p|hy_1, p|hz_1;
  p|a;
  p|origin_x, p|origin_y, p|origin_z;
  p|nlevels, p|maxlevels, p|toplevel;
  p|approx, p|split;
  p|qh, p|eh;
  p|grid_len;
  p|grid_idstart;
  p|scaling;
  p|gc_len;
  p|gc_idstart;
  p|gc;
  p|gctop_len;
  p|gctop_idstart;
  p|gctop;
  p|num_clients_qh, p|num_clients_eh;
  p|num_anterpolation_chares;
  p|num_interpolation_chares;
  p|num_restriction_chares;
  p|num_prolongation_chares;
  p|num_gridcutoff_chares;
  p|dim_gridcutoff_chares;
  p|dim_gridtransfer_chares;
  p|num_total_restriction_chares;
  p|num_total_prolongation_chares;
  p|num_total_gridcutoff_chares;
  p|num_energy_chares;
  p|num_points_per_chare;
  p|self_energy_const;
}

void MsmMsaData::print()
{
#if 0
  printf("MSM data:\n");
  printf("ispx=%d ispy=%d ispz=%d\n", ispx, ispy, ispz);
  printf("hx=%g  hy=%g  hz=%g\n", 1/hx_1, 1/hy_1, 1/hz_1);
  printf("a=%g\n", a);
  printf("origin = %g %g %g\n", origin_x, origin_y, origin_z);
  printf("nlevels=%d maxlevels=%d toplevel=%d\n", nlevels, maxlevels, toplevel);
  printf("approx=%d split=%d\n", approx, split);
#endif
}

struct ComputeMsmMsaAtom {
  Position position;
  float msmcharge;    // scaled for MSM
  int id;
};

class MsmMsaCoordMsg : public CMessage_MsmMsaCoordMsg {
public:
  int numAtoms;
  Lattice lattice;
  ComputeMsmMsaAtom *coord;
};

class ComputeMsmMsaMgr : public CBase_ComputeMsmMsaMgr {
public:
  ComputeMsmMsaMgr();                    // entry
  ~ComputeMsmMsaMgr();

  void initialize(CkQdMsg *);         // entry with message

  void setCompute(ComputeMsmMsa *c) { msmCompute = c; }  // local

  void recvMsmMsaData(const MsmMsaData &);  // entry with parameter marshalling

  void initWorkers(CkQdMsg *);        // entry with message
  void startWorkers(CkQdMsg *);       // entry with message

  void anterpolate(MsmMsaCoordMsg *);    // entry with coordinates message
  void interpolate(CkQdMsg *);        // entry with message

  const MsmMsaData &getMsmMsaData() const { return msmData; }  // local

private:
  CProxy_ComputeMsmMsaMgr msmProxy;
  ComputeMsmMsa *msmCompute;

  MsmMsaCoordMsg *coordMsg;
  int numAtoms;
  MsmMsaForce *force;
  int numForces;

  MsmMsaData msmData;

  CProxy_MsmMsaLevel msmLevelProxy;   // 1D chare array (number of grid levels)
  CProxy_MsmMsaEnergy msmEnergyProxy; // 1D chare array

  MsmMsaGrid::Accum qhacc;  // accumulate charge grid for anterpolation
  MsmMsaGrid::Read qhread;  // read charge grid after anterpolation
  MsmMsaGrid::Accum ehacc;  // accumulate potential grid before interpolation
  MsmMsaGrid::Read ehread;  // read potential grid for interpolation
  int msa_setup_anterpolation;   // have MSAs been initially set up?
  int msa_setup_interpolation;   // have MSAs been initially set up?
};

class MsmMsaLevel : public CBase_MsmMsaLevel {
public:
  MsmMsaLevel(MsmMsaGrid &qh, MsmMsaGrid &eh, MsmMsaGrid &q2h, MsmMsaGrid &e2h);
  MsmMsaLevel(MsmMsaGrid &qh, MsmMsaGrid &eh);
  MsmMsaLevel(CkMigrateMessage *m) { }
  void compute();
private:
  int lastlevel;
  CProxy_MsmMsaGridCutoff gridcutoff;
  CProxy_MsmMsaRestriction restriction;
  CProxy_MsmMsaProlongation prolongation;
};

class MsmMsaGridCutoff : public CBase_MsmMsaGridCutoff {
public:
  MsmMsaGridCutoff(int level, MsmMsaGrid &qh, MsmMsaGrid &eh);
  MsmMsaGridCutoff(CkMigrateMessage *m) { }
  void compute();
private:
  MsmMsaGrid qh, eh;    // MSA handles to this level charge and potential grids
  MsmMsaGrid::Accum qhacc, ehacc;
  MsmMsaGrid::Read qhread, ehread;
  int msa_setup;
  int mylevel;       // which level am I on?
  int mia, mja, mka; // my lowest grid point index of my sub-cube
  int mib, mjb, mkb; // my highest grid point index of my sub-cube
};

class MsmMsaRestriction : public CBase_MsmMsaRestriction {
public:
  MsmMsaRestriction(int level, MsmMsaGrid &qh, MsmMsaGrid &q2h);
  MsmMsaRestriction(CkMigrateMessage *m) { }
  void compute();
private:
  MsmMsaGrid qh, q2h;   // MSA handles to charge grids mylevel, (mylevel+1)
  MsmMsaGrid::Accum qhacc, q2hacc;
  MsmMsaGrid::Read qhread, q2hread;
  int msa_setup;
  int mylevel;       // which level am I on?
  int mia, mja, mka; // my lowest grid point index of (mylevel+1) sub-cube
  int mib, mjb, mkb; // my highest grid point index of (mylevel+1) sub-cube
};

class MsmMsaProlongation : public CBase_MsmMsaProlongation {
public:
  MsmMsaProlongation(int level, MsmMsaGrid &eh, MsmMsaGrid &e2h);
  MsmMsaProlongation(CkMigrateMessage *m) { }
  void compute();
private:
  MsmMsaGrid eh, e2h;   // MSA handles to potential grids mylevel, (mylevel+1)
  MsmMsaGrid::Accum ehacc, e2hacc;
  MsmMsaGrid::Read ehread, e2hread;
  int msa_setup;
  int mylevel;       // which level am I on?
  int mia, mja, mka; // my lowest grid point index of (mylevel+1) sub-cube
  int mib, mjb, mkb; // my highest grid point index of (mylevel+1) sub-cube
};

class MsmMsaEnergy : public CBase_MsmMsaEnergy {
public:
  MsmMsaEnergy(MsmMsaGrid &qh, MsmMsaGrid &eh);
  MsmMsaEnergy(CkMigrateMessage *m) { }
  ~MsmMsaEnergy();
  void compute();
private:
  MsmMsaGrid qh, eh;
  MsmMsaGrid::Accum qhacc, ehacc;
  MsmMsaGrid::Read qhread, ehread;
  float *qhbuffer;
  int msa_setup;
  int mli, mlj, mlk; // my 3D index within level 0
  int mia, mja, mka; // my lowest grid point index of my sub-cube
  int mib, mjb, mkb; // my highest grid point index of my sub-cube
  SubmitReduction *reduction;
};


ComputeMsmMsa::ComputeMsmMsa(ComputeID c) :
  ComputeHomePatches(c)
{
  CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->setCompute(this);
  SimParameters *simParams = Node::Object()->simParameters;
  qscaling = sqrt(COULOMB / simParams->dielectric);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeMsmMsa::~ComputeMsmMsa()
{
}

void ComputeMsmMsa::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics ) {
    for (ap = ap.begin();  ap != ap.end();  ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
      reduction->submit();
    }
    return;
  }

  // allocate message
  int numLocalAtoms = 0;
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  MsmMsaCoordMsg *msg = new (numLocalAtoms, 0) MsmMsaCoordMsg;
  msg->numAtoms = numLocalAtoms;
  msg->lattice = patchList[0].p->flags.lattice;
  ComputeMsmMsaAtom *data_ptr = msg->coord;

  // get positions
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    CompAtom *x = (*ap).positionBox->open();
    CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();

    for (int i=0;  i < numAtoms;  i++) {
      data_ptr->position = x[i].position;
      data_ptr->msmcharge = qscaling * x[i].charge;
      data_ptr->id = xExt[i].id;
      ++data_ptr;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  }

  CProxy_ComputeMsmMsaMgr msmProxy(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr);
  msmProxy[CkMyPe()].anterpolate(msg);
  msmProxy[CkMyPe()].interpolate(new CkQdMsg);
}


void ComputeMsmMsa::saveResults(int n, const MsmMsaForce force[], double self_energy)
{
  ResizeArrayIter<PatchElem> ap(patchList);

  // add in forces
  int j = 0;
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();
    for (int i=0;  i < numAtoms;  i++, j++) {
      f[i].x += force[j].x;
      f[i].y += force[j].y;
      f[i].z += force[j].z;
    }
    (*ap).forceBox->close(&r);
  }
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) -= self_energy;
  reduction->submit();
}


ComputeMsmMsaMgr::ComputeMsmMsaMgr() :
  msmProxy(thisgroup), msmCompute(0),
  coordMsg(0), numAtoms(0), force(0), numForces(0),
  msa_setup_anterpolation(0), msa_setup_interpolation(0)
{
  CkpvAccess(BOCclass_group).computeMsmMsaMgr = thisgroup;
}

ComputeMsmMsaMgr::~ComputeMsmMsaMgr()
{
  delete [] force;
}


static void rescale_nonperiodic_cell(
    Vector &u, Vector &v, Vector &w, Vector &c,
    Vector &ru, Vector &rv, Vector &rw,
    const Vector &smin, const Vector &smax,
    double padding, double gridspacing, int isperiodic)
{
  const double ZERO_TOLERANCE = 1e-4;
  double xlen = 1, inv_xlen = 1;  // don't rescale periodic directions
  double ylen = 1, inv_ylen = 1;
  double zlen = 1, inv_zlen = 1;
  double ulen_1 = u.rlength();
  double vlen_1 = v.rlength();
  double wlen_1 = w.rlength();
  double rpadx = padding * ulen_1;  // padding distance in recip space
  double rpady = padding * vlen_1;
  double rpadz = padding * wlen_1;
  double rhx = gridspacing * ulen_1;  // grid spacing in recip space
  double rhy = gridspacing * vlen_1;
  double rhz = gridspacing * wlen_1;
  double xmin = smin.x, xmax = smax.x;  // scaled coordinate extremes
  double ymin = smin.y, ymax = smax.y;
  double zmin = smin.z, zmax = smax.z;
  Vector rc = 0;  // don't move center along periodic directions
  int is_periodic_x = (isperiodic & PERIODIC_VEC1);
  int is_periodic_y = (isperiodic & PERIODIC_VEC2);
  int is_periodic_z = (isperiodic & PERIODIC_VEC3);

#if 0
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
#endif

  if ( ! is_periodic_x) {
    xmax += rpadx;  // pad the edges
    xmin -= rpadx;
    if (rhx > 0) {  // restrict center to rhx lattice points
      double mupper = ceil(xmax / (2*rhx));
      double mlower = floor(xmin / (2*rhx));
      xmax = 2*rhx*mupper;
      xmin = 2*rhx*mlower;
    }
    rc.x = 0.5*(xmin + xmax);
    xlen = xmax - xmin;
    if (xlen < ZERO_TOLERANCE) {
      xlen = 1;  // leave scaling unchanged
    }
    else {
      inv_xlen = 1/xlen;
    }
  }

  if ( ! is_periodic_y) {
    ymax += rpady;  // pad the edges
    ymin -= rpady;
    if (rhy > 0) {  // restrict center to rhy lattice points
      double mupper = ceil(ymax / (2*rhy));
      double mlower = floor(ymin / (2*rhy));
      ymax = 2*rhy*mupper;
      ymin = 2*rhy*mlower;
    }
    rc.y = 0.5*(ymin + ymax);
    ylen = ymax - ymin;
    if (ylen < ZERO_TOLERANCE) {
      ylen = 1;  // leave scaling unchanged
    }
    else {
      inv_ylen = 1/ylen;
    }
  }

  if ( ! is_periodic_z) {
    zmax += rpadz;  // pad the edges
    zmin -= rpadz;
    if (rhz > 0) {  // restrict center to rhz lattice points
      double mupper = ceil(zmax / (2*rhz));
      double mlower = floor(zmin / (2*rhz));
      zmax = 2*rhz*mupper;
      zmin = 2*rhz*mlower;
    }
    rc.z = 0.5*(zmin + zmax);
    zlen = zmax - zmin;
    if (zlen < ZERO_TOLERANCE) {
      zlen = 1;  // leave scaling unchanged
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

  // transform new center back to real space, then rescale basis vectors
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


struct MsmMsaInterpParams {
  int nu;
  int stencil;
  int omega;
};

// ordering must be identical to APPROX enum constants
static const MsmMsaInterpParams InterpParams[] = {
  { 1, 4, 6 },    // cubic
  { 2, 6, 10 },   // quintic
  { 2, 6, 10 },   // quintic, C2
  { 3, 8, 14 },   // septic
  { 3, 8, 14 },   // septic, C3
  { 4, 10, 18 },  // nonic
  { 4, 10, 18 },  // nonic, C4
  { 1, 4, 6 },    // B-spline
};


static int setup_hgrid_1d(
    double len,           // cell length
    double gridspacing,   // desired grid spacing
    double &hh,           // determine h to fit cell length
    int &nn,              // determine number of h's covering cell
    int &aindex,          // determine smallest grid index
    int &bindex,          // determine largest grid index
    int isperiodic,       // is this dimension periodic?
    int approx            // which approximation
    ) {

  const int nu = InterpParams[approx].nu;  // interp stencil radius

  // ASSERT(hmax > 0);
  if (isperiodic) {
    const double hmin = (4./5) * gridspacing;  // minimum bound on h
    const double hmax = 1.5 * hmin;
    double h = len;
    int n = 1;    // start with one grid point across domain
    while (h >= hmax) {
      h *= 0.5;   // halve h
      n <<= 1;    // double grid points
    }
    if (h < hmin) {
      if (n < 4) {  // either len is too small or hmin is too large
        return -1;
      }
      h *= (4./3);  // scale h by 4/3
      n >>= 2;      // scale n by 3/4
      n *= 3;
    }
    // now we have:  hmin <= h < hmax
    // now we have:  n is power of two times no more than one power of 3
    hh = h;
    nn = n;
    aindex = 0;
    bindex = n-1;
  }
  else {  // non-periodic
    double h = gridspacing;
    int n = (int) floorf(len / h) + 1;
    hh = h;
    nn = n;
    aindex = -nu;
    bindex = n + nu;
  }
  return 0;
}


void ComputeMsmMsaMgr::initialize(CkQdMsg *msg)
{
  delete msg;
  //printf("MSM initialize PE=%d\n", CkMyPe());
  if (CkMyPe() != 0) return;  // initialize only on PE 0, broadcast MsmMsaData

  // initialize MSM here
  SimParameters *simParams = Node::Object()->simParameters;
  const double a = simParams->cutoff;
  const double gridspacing = simParams->MSMGridSpacing;
  const double padding = simParams->MSMPadding;
  const int approx = simParams->MSMApprox;
  const int split = simParams->MSMSplit;
  int nlevels = simParams->MSMLevels;
  const Lattice &lattice = simParams->lattice;
  int msmflags = 0;
  msmflags |= (lattice.a_p() ? PERIODIC_VEC1 : 0);
  msmflags |= (lattice.b_p() ? PERIODIC_VEC2 : 0);
  msmflags |= (lattice.c_p() ? PERIODIC_VEC3 : 0);
  int allperiodic = (PERIODIC_VEC1 | PERIODIC_VEC2 | PERIODIC_VEC3);
  Vector u(lattice.a()), v(lattice.b()), w(lattice.c());
  Vector c(lattice.origin());
  Vector ru(lattice.a_r()), rv(lattice.b_r()), rw(lattice.c_r());
  if ( (msmflags & allperiodic) != allperiodic ) {
    ScaledPosition smin, smax;
    Node::Object()->pdb->get_extremes(smin, smax);
    //printf("smin = %g %g %g  smax = %g %g %g\n",
    //    smin.x, smin.y, smin.z, smax.x, smax.y, smax.z);
    rescale_nonperiodic_cell(u, v, w, c, ru, rv, rw, smin, smax,
        padding, gridspacing, msmflags);
  }

  if (u.x <= 0 || u.y != 0 || u.z != 0 ||
      v.x != 0 || v.y <= 0 || v.z != 0 ||
      w.x != 0 || w.y != 0 || w.z <= 0) {
    NAMD_die("MSM requires cell basis be along x, y, and z coordinate axes.");
  }

  // setup grids
  MsmMsaData &p = msmData;  // the MSM data will be broadcast to all PEs

  const int nu = InterpParams[approx].nu;
  const int omega = InterpParams[approx].omega;
  const int ispx = ((msmflags & PERIODIC_VEC1) != 0);
  const int ispy = ((msmflags & PERIODIC_VEC2) != 0);
  const int ispz = ((msmflags & PERIODIC_VEC3) != 0);
  const int ispany = ((msmflags & PERIODIC_ALL) != 0);

  const double xlen = u.x;  // XXX want to extend to non-orthogonal cells
  const double ylen = v.y;
  const double zlen = w.z;

  double scaling;
  double d;  // temporary for SPOLY derivative

  int ia, ib, ja, jb, ka, kb, ni, nj, nk;
  int nx, ny, nz;  // counts the grid points that span just the domain

  double hx, hy, hz;
  double gx, gy, gz;

  int i, j, k, n;
  int index;
  int level, toplevel, maxlevels;
  int lastnelems = 1;
  int isclamped = 0;
  int done, alldone;

  int rc = 0;  // return code

  //CkPrintf("ispx=%d  ispy=%d  ispz=%d\n", ispx, ispy, ispz);
  p.ispx = ispx;
  p.ispy = ispy;
  p.ispz = ispz;

  rc = setup_hgrid_1d(xlen, gridspacing, hx, nx, ia, ib, ispx, approx);
  if (rc)  NAMD_die("MSM failed to setup grid along x dimension.");

  rc = setup_hgrid_1d(ylen, gridspacing, hy, ny, ja, jb, ispy, approx);
  if (rc)  NAMD_die("MSM failed to setup grid along y dimension.");

  rc = setup_hgrid_1d(zlen, gridspacing, hz, nz, ka, kb, ispz, approx);
  if (rc)  NAMD_die("MSM failed to setup grid along z dimension.");

  p.a = (float) a;  // cutoff distance

  p.hx_1 = (float) (1/hx);  // inverse of grid spacing
  p.hy_1 = (float) (1/hy);
  p.hz_1 = (float) (1/hz);

  // XXX set coordinate for h-grid (0,0,0) point
  gx = c.x - ((nx >> 1) * hx);
  gy = c.y - ((ny >> 1) * hy);
  gz = c.z - ((nz >> 1) * hz);

  p.origin_x = (float) gx;  // grid(0,0,0) location
  p.origin_y = (float) gy;
  p.origin_z = (float) gz;

  ni = ib - ia + 1;
  nj = jb - ja + 1;
  nk = kb - ka + 1;

  p.approx = approx;
  p.split = split;

  if (nlevels <= 0) {
    // automatically set number of levels
    n = ni;
    if (n < nj) n = nj;
    if (n < nk) n = nk;
    for (maxlevels = 1;  n > 0;  n >>= 1)  maxlevels++;
    nlevels = maxlevels;
    if (ispany == 0) {  // no periodicity
      int omega3 = omega * omega * omega;
      int nhalf = (int) sqrtf(ni*nj*nk);  // scale down for performance?
      lastnelems = (nhalf > omega3 ? nhalf : omega3);
      isclamped = 1;
    }
  }
  else {
    // user-defined number of levels
    maxlevels = nlevels;
  }

  p.maxlevels = maxlevels;
  p.grid_len.resize(maxlevels);
  p.grid_idstart.resize(maxlevels);
  p.scaling.resize(maxlevels);

#if 0
  /* allocate any additional levels that may be needed */
  if (pm->maxlevels < maxlevels) {
    void *vqh, *veh, *vgc;
    if (issprec) {
      vqh = realloc(pm->qh, maxlevels * sizeof(NL_MsmMsagrid_float));
      if (NULL == vqh) return NL_MSM_ERROR_MALLOC;
      veh = realloc(pm->eh, maxlevels * sizeof(NL_MsmMsagrid_float));
      if (NULL == veh) return NL_MSM_ERROR_MALLOC;
      vgc = realloc(pm->gc, maxlevels * sizeof(NL_MsmMsagrid_float));
      if (NULL == vgc) return NL_MSM_ERROR_MALLOC;
      pm->qh_f = (NL_MsmMsagrid_float *) vqh;
      pm->eh_f = (NL_MsmMsagrid_float *) veh;
      pm->gc_f = (NL_MsmMsagrid_float *) vgc;
      /* initialize the newest grids appended to array */
      for (level = pm->maxlevels;  level < maxlevels;  level++) {
        GRID_INIT( &(pm->qh_f[level]) );
        GRID_INIT( &(pm->eh_f[level]) );
        GRID_INIT( &(pm->gc_f[level]) );
      }
    }
    else {
      vqh = realloc(pm->qh, maxlevels * sizeof(NL_MsmMsagrid_double));
      if (NULL == vqh) return NL_MSM_ERROR_MALLOC;
      veh = realloc(pm->eh, maxlevels * sizeof(NL_MsmMsagrid_double));
      if (NULL == veh) return NL_MSM_ERROR_MALLOC;
      vgc = realloc(pm->gc, maxlevels * sizeof(NL_MsmMsagrid_double));
      if (NULL == vgc) return NL_MSM_ERROR_MALLOC;
      pm->qh = (NL_MsmMsagrid_double *) vqh;
      pm->eh = (NL_MsmMsagrid_double *) veh;
      pm->gc = (NL_MsmMsagrid_double *) vgc;
      /* initialize the newest grids appended to array */
      for (level = pm->maxlevels;  level < maxlevels;  level++) {
        GRID_INIT( &(pm->qh[level]) );
        GRID_INIT( &(pm->eh[level]) );
        GRID_INIT( &(pm->gc[level]) );
      }
    }
    pm->maxlevels = maxlevels;
  }
#endif

  level = 0;
  done = 0;
  alldone = 0;
  scaling = 1.0;
  do {
#if 0
    if (issprec) {
      GRID_RESIZE( &(pm->qh_f[level]), float, ia, ni, ja, nj, ka, nk);
      GRID_RESIZE( &(pm->eh_f[level]), float, ia, ni, ja, nj, ka, nk);
    }
    else {
      GRID_RESIZE( &(pm->qh[level]), double, ia, ni, ja, nj, ka, nk);
      GRID_RESIZE( &(pm->eh[level]), double, ia, ni, ja, nj, ka, nk);
    }
#endif
    p.grid_len[level] = Int3(ni,nj,nk);
    p.grid_idstart[level] = Int3(ia,ja,ka);
    p.scaling[level] = (float) scaling;
    scaling *= 0.5;

    if (++level == nlevels)    done |= 0x07;  // user limit on levels

    alldone = (done == 0x07);  // make sure all dimensions are done

    if (isclamped) {
      int nelems = ni * nj * nk;
      if (nelems <= lastnelems)  done |= 0x07;
    }

    if (ispx) {
      ni >>= 1;
      ib = ni-1;
      if (ni & 1)              done |= 0x07;  // == 3 or 1
      else if (ni == 2)        done |= 0x01;  // can do one more
    }
    else {
      ia = -((-ia+1)/2) - nu;
      ib = (ib+1)/2 + nu;
      ni = ib - ia + 1;
      if (ni <= omega)         done |= 0x01;  // can do more restrictions
    }

    if (ispy) {
      nj >>= 1;
      jb = nj-1;
      if (nj & 1)              done |= 0x07;  // == 3 or 1
      else if (nj == 2)        done |= 0x02;  // can do one more
    }
    else {
      ja = -((-ja+1)/2) - nu;
      jb = (jb+1)/2 + nu;
      nj = jb - ja + 1;
      if (nj <= omega)         done |= 0x02;  // can do more restrictions
    }

    if (ispz) {
      nk >>= 1;
      kb = nk-1;
      if (nk & 1)              done |= 0x07;  // == 3 or 1
      else if (nk == 2)        done |= 0x04;  // can do one more
    }
    else {
      ka = -((-ka+1)/2) - nu;
      kb = (kb+1)/2 + nu;
      nk = kb - ka + 1;
      if (nk <= omega)         done |= 0x04;  // can do more restrictions
    }

  } while ( ! alldone );
  nlevels = level;
  p.nlevels = nlevels;

  toplevel = (ispany ? nlevels : nlevels - 1);
  p.toplevel = -1;  // set only for non-periodic system

  // ellipsoid axes for grid cutoff weights
  ni = (int) ceil(2*a/hx) - 1;
  nj = (int) ceil(2*a/hy) - 1;
  nk = (int) ceil(2*a/hz) - 1;

  p.gc_len = Int3(2*ni+1, 2*nj+1, 2*nk+1);
  p.gc_idstart = Int3(-ni, -nj, -nk);
  p.gc.resize( (2*ni+1)*(2*nj+1)*(2*nk+1) );

  index = 0;
  for (k = -nk;  k <= nk;  k++) {
    for (j = -nj;  j <= nj;  j++) {
      for (i = -ni;  i <= ni;  i++) {
        double s, t, gs, gt, g;
        s = sqrt((i*hx)*(i*hx) + (j*hy)*(j*hy) + (k*hz)*(k*hz)) / a;
        t = 0.5 * s;
        if (t >= 1) {
          g = 0;
        }
        else if (s >= 1) {
          gs = 1/s;
          SPOLY(&gt, &d, t, split);
          g = (gs - 0.5 * gt) / a;
        }
        else {
          SPOLY(&gs, &d, s, split);
          SPOLY(&gt, &d, t, split);
          g = (gs - 0.5 * gt) / a;
        }
        p.gc[index] = (float) g;
        index++;
      }
    }
  }

  if (toplevel < nlevels) {
    // nonperiodic in all dimensions,
    // calculate top level weights, ellipsoid axes are length of grid
#ifdef MSM_DEBUG
    CkPrintf("MSM setting top level weights\n");
#endif
    ni = p.grid_len[toplevel].nx;
    nj = p.grid_len[toplevel].ny;
    nk = p.grid_len[toplevel].nz;

    p.toplevel = toplevel;
    p.gctop_len = Int3(2*ni+1, 2*nj+1, 2*nk+1);
    p.gctop_idstart = Int3(-ni, -nj, -nk);
    p.gctop.resize( (2*ni+1)*(2*nj+1)*(2*nk+1) );

    index = 0;
    for (k = -nk;  k <= nk;  k++) {
      for (j = -nj;  j <= nj;  j++) {
        for (i = -ni;  i <= ni;  i++) {
          double s, gs;
          s = sqrt((i*hx)*(i*hx) + (j*hy)*(j*hy) + (k*hz)*(k*hz)) / a;
          if (s >= 1) {
            gs = 1/s;
          }
          else {
            SPOLY(&gs, &d, s, split);
          }
          p.gctop[index] = (float) (gs / a);
          index++;
        }
      }
    } // end loops over k-j-i for coarsest level weights
  }

  // calculate self energy factor for splitting
  double s, gs;
  s = 0;
  SPOLY(&gs, &d, s, split);
  p.self_energy_const = gs / a;
#if 0
  s = 0;
  for (n = 0;  n < natoms;  n++) {
    double q = atom[4*n + 3];
    s += q * q;
  }
  self_energy *= 0.5*s;
#endif

  // count number of client chares to the MSAs
  p.num_points_per_chare.nx = simParams->MSMBlockSizeX;
  p.num_points_per_chare.ny = simParams->MSMBlockSizeY;
  p.num_points_per_chare.nz = simParams->MSMBlockSizeZ;
  if (nlevels > 1) {
    p.num_restriction_chares.resize(nlevels-1);
    p.num_prolongation_chares.resize(nlevels-1);
    p.dim_gridtransfer_chares.resize(nlevels-1);
  }
  p.num_gridcutoff_chares.resize(nlevels);
  p.dim_gridcutoff_chares.resize(nlevels);
  p.num_clients_qh.resize(nlevels);
  p.num_clients_eh.resize(nlevels);

  p.num_anterpolation_chares = (PatchMap::Object())->numNodesWithPatches();
  p.num_interpolation_chares = p.num_anterpolation_chares;
  p.num_total_gridcutoff_chares = 0;
  p.num_total_restriction_chares = 0;
  p.num_total_prolongation_chares = 0;
  for (n = 0;  n < nlevels;  n++) {
    int ni = p.grid_len[n].nx;
    int nj = p.grid_len[n].ny;
    int nk = p.grid_len[n].nz;
    int nci = ROUNDUP_QUOTIENT(ni, p.num_points_per_chare.nx);
    int ncj = ROUNDUP_QUOTIENT(nj, p.num_points_per_chare.ny);
    int nck = ROUNDUP_QUOTIENT(nk, p.num_points_per_chare.nz);
    int nc = nci * ncj * nck;
#ifdef MSM_DEBUG
    CkPrintf("n=%d  nc=%d  nci=%d ncj=%d nck=%d\n", n, nc, nci, ncj, nck);
#endif
    p.num_gridcutoff_chares[n] = nc;
    p.dim_gridcutoff_chares[n].nx = nci;
    p.dim_gridcutoff_chares[n].ny = ncj;
    p.dim_gridcutoff_chares[n].nz = nck;
    p.num_total_gridcutoff_chares += nc;
    if (n > 0) {
      p.num_restriction_chares[n-1] = nc;
      p.num_prolongation_chares[n-1] = nc;
      p.dim_gridtransfer_chares[n-1].nx = nci;
      p.dim_gridtransfer_chares[n-1].ny = ncj;
      p.dim_gridtransfer_chares[n-1].nz = nck;
      p.num_total_restriction_chares += nc;
      p.num_total_prolongation_chares += nc;
    }
    else {
      p.num_energy_chares = nc;
    }
  }

  p.num_clients_qh[0] = p.num_anterpolation_chares +
    p.num_gridcutoff_chares[0] + p.num_energy_chares;
  p.num_clients_eh[0] = p.num_interpolation_chares +
    p.num_gridcutoff_chares[0] + p.num_energy_chares;
  if (nlevels > 1) {
    p.num_clients_qh[0] += p.num_restriction_chares[0];
    p.num_clients_eh[0] += p.num_prolongation_chares[0];
  }
  for (n = 1;  n < nlevels;  n++) {
    p.num_clients_qh[n] = p.num_gridcutoff_chares[n] +
      p.num_restriction_chares[n-1];
    p.num_clients_eh[n] = p.num_gridcutoff_chares[n] +
      p.num_prolongation_chares[n-1];
    if (n < nlevels-1) {
      p.num_clients_qh[n] += p.num_restriction_chares[n];
      p.num_clients_eh[n] += p.num_prolongation_chares[n];
    }
  }

#ifdef MSM_DEBUG
  CkPrintf("nlevels = %d\n", nlevels);
  for (n = 0;  n < nlevels;  n++) {
    CkPrintf("num clients qh[%d] = %d\n", n, p.num_clients_qh[n]);
    CkPrintf("num clients eh[%d] = %d\n", n, p.num_clients_eh[n]);
  }
  CkPrintf("num anterpolation chares = %d\n", p.num_anterpolation_chares);
  CkPrintf("num interpolation chares = %d\n", p.num_interpolation_chares);
  for (n = 0;  n < nlevels;  n++) {
    CkPrintf("num grid cutoff chares[%d] = %d\n",
        n, p.num_gridcutoff_chares[n]);
  }
  for (n = 0;  n < nlevels-1;  n++) {
    CkPrintf("num restriction chares[%d] = %d\n",
        n, p.num_restriction_chares[n]);
  }
  for (n = 0;  n < nlevels-1;  n++) {
    CkPrintf("num prolongation chares[%d] = %d\n",
        n, p.num_prolongation_chares[n]);
  }
#endif

#if 0
  // allocate MSAs
  if (qh || eh) {
    CkAbort("attempted to reallocate MSAs");
  }
  char *qh_memory = new char[nlevels * sizeof(MsmMsaGrid)];
  qh = static_cast<MsmMsaGrid *>((void *)qh_memory);
  char *eh_memory = new char[nlevels * sizeof(MsmMsaGrid)];
  eh = static_cast<MsmMsaGrid *>((void *)eh_memory);
#else
  p.qh.resize(nlevels);
  p.eh.resize(nlevels);
#endif
  for (n = 0;  n < nlevels;  n++) {
    ia = p.grid_idstart[n].nx;
    ib = p.grid_len[n].nx + ia - 1;
    ja = p.grid_idstart[n].ny;
    jb = p.grid_len[n].ny + ja - 1;
    ka = p.grid_idstart[n].nz;
    kb = p.grid_len[n].nz + ka - 1;
#if 0
    // using placement new to call non-default constructor
    // on each MsmMsaGrid element
    new(qh_memory + n*sizeof(MsmMsaGrid))MsmMsaGrid(ia, ib, ja, jb, ka, kb,
        p.num_clients_qh[n]);
    new(eh_memory + n*sizeof(MsmMsaGrid))MsmMsaGrid(ia, ib, ja, jb, ka, kb,
        p.num_clients_eh[n]);
#else
    p.qh[n] = MsmMsaGrid(ia, ib, ja, jb, ka, kb, p.num_clients_qh[n]);
    p.eh[n] = MsmMsaGrid(ia, ib, ja, jb, ka, kb, p.num_clients_eh[n]);
#endif
  }
  msmProxy.recvMsmMsaData(msmData);  // broadcast MsmMsaData to chare group
}


void ComputeMsmMsaMgr::recvMsmMsaData(const MsmMsaData &m)
{
  //printf("MSM recvMsmMsaData PE=%d\n", CkMyPe());
  if (CkMyPe() != 0) {
    msmData = m;
  }
  else {
    msmData.print();  // PE 0
  }
}


void ComputeMsmMsaMgr::initWorkers(CkQdMsg *msg)
{
  delete msg;
  //printf("MSM initWorkers PE=%d\n", CkMyPe());
  if (CkMyPe() != 0) return;
  // PE 0 creates the compute chare arrays

  MsmMsaData &p = msmData;
  int n;
  msmLevelProxy = CProxy_MsmMsaLevel::ckNew();  // create empty chare array
  for (n = 0;  n < p.nlevels-1;  n++) {
    msmLevelProxy[n].insert(p.qh[n], p.eh[n], p.qh[n+1], p.eh[n+1]);
  }
  msmLevelProxy[n].insert(p.qh[n], p.eh[n]);  // top level
  msmLevelProxy.doneInserting();
#ifdef MSM_DEBUG
  CkPrintf("Created %d MsmMsaLevel chares\n", p.nlevels);
#endif

  msmEnergyProxy = CProxy_MsmMsaEnergy::ckNew(p.qh[0], p.eh[0], p.num_energy_chares);
#ifdef MSM_DEBUG
  CkPrintf("Created %d MsmMsaEnergy chares\n", p.num_energy_chares);
#endif
}


void ComputeMsmMsaMgr::startWorkers(CkQdMsg *msg)
{
  delete msg;
  //printf("MSM startWorkers PE=%d\n", CkMyPe());
  if (CkMyPe() != 0) return;
  // PE 0 starts the workers;
  // they loop forever, waiting for compute work from the MSAs

  msmLevelProxy.compute();
  msmEnergyProxy.compute();
}


/** Approximation formulaes are up to degree 9 polynomials. */
enum { MAX_POLY_DEGREE = 9 };

/** Degree of polynomial basis function Phi.
 * Must be listed in same order as APPROX enum from msm.h */
static const int PolyDegree[APPROX_END] = {
  3, 5, 5, 7, 7, 9, 9, 3,
};


void ComputeMsmMsaMgr::anterpolate(MsmMsaCoordMsg *msg)
{
#ifdef MSM_DEBUG
  CkPrintf("Anterpolation compute %d starting, PE %d\n", thisIndex, CkMyPe());
#endif
  coordMsg = msg;  // save message for interpolation
  ComputeMsmMsaAtom *coord = coordMsg->coord;
  int numAtoms = coordMsg->numAtoms;
  MsmMsaData &p = msmData;
  if ( ! msa_setup_anterpolation) {
    p.qh[0].enroll(p.num_clients_qh[0]);
    qhacc = p.qh[0].getInitialAccum();
    qhread = qhacc.syncToRead();  // immediately sync to read phase
    msa_setup_anterpolation = 1;
  }
  int ni = p.grid_len[0].nx;
  int nj = p.grid_len[0].ny;
  int nk = p.grid_len[0].nz;
  int ia = p.grid_idstart[0].nx;
  int ja = p.grid_idstart[0].ny;
  int ka = p.grid_idstart[0].nz;
  int ib = ia + ni - 1;
  int jb = ja + nj - 1;
  int kb = ka + nk - 1;
  qhacc = qhread.syncToEAccum();  // we accumulate to it
#ifdef MSM_DEBUG
  CkPrintf("Anterpolation compute %d doing work, PE %d\n", thisIndex, CkMyPe());
#endif
  //
  // here is the compute work
  //

  float xphi[MAX_POLY_DEGREE+1];  // Phi stencil along x-dimension
  float yphi[MAX_POLY_DEGREE+1];  // Phi stencil along y-dimension
  float zphi[MAX_POLY_DEGREE+1];  // Phi stencil along z-dimension

  int n;
  const int s_size = PolyDegree[p.approx] + 1;         // stencil size
  const int s_edge = (PolyDegree[p.approx] - 1) >> 1;  // stencil "edge" size

  for (n = 0;  n < numAtoms;  n++) {
    // atomic charge
    float q = coord[n].msmcharge;
    if (0==q) continue;

    // distance between atom and origin measured in grid points
    // (does subtraction and division in double, then converts to float)
    float rx_hx = (coord[n].position.x - p.origin_x) * p.hx_1;
    float ry_hy = (coord[n].position.y - p.origin_y) * p.hy_1;
    float rz_hz = (coord[n].position.z - p.origin_z) * p.hz_1;

    // find smallest numbered grid point in stencil
    int ilo = (int) floorf(rx_hx) - s_edge;
    int jlo = (int) floorf(ry_hy) - s_edge;
    int klo = (int) floorf(rz_hz) - s_edge;

    // calculate Phi stencils along each dimension
    float delta;
    delta = rx_hx - (float) ilo;
    STENCIL_1D(xphi, delta, p.approx);
    delta = ry_hy - (float) jlo;
    STENCIL_1D(yphi, delta, p.approx);
    delta = rz_hz - (float) klo;
    STENCIL_1D(zphi, delta, p.approx);

    // test to see if stencil is within edges of grid
    int iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                     ja <= jlo && (jlo+(s_size-1)) <= jb &&
                     ka <= klo && (klo+(s_size-1)) <= kb );

    if ( iswithin ) {  // no wrapping needed
      // determine charge on cube of grid points around atom
      int i, j, k;
      for (k = 0;  k < s_size;  k++) {
        float ck = zphi[k] * q;
        for (j = 0;  j < s_size;  j++) {
          float cjk = yphi[j] * ck;
          for (i = 0;  i < s_size;  i++) {
            qhacc.accumulate(i+ilo,j+jlo,k+klo) += xphi[i] * cjk;
          }
        }
      }
    } // if is within
    else {  // requires wrapping around grid
      int ip, jp, kp;

      // adjust ilo, jlo, klo so they are within grid indexing
      if (p.ispx) {
        if      (ilo < ia) do { ilo += ni; } while (ilo < ia);
        else if (ilo > ib) do { ilo -= ni; } while (ilo > ib);
      }
      else if (ilo < ia || (ilo+(s_size-1)) > ib) {
        CkPrintf("Error anterpolation:  ilo=%d  out-of-range [%d,%d]\n",
            ilo, ia, (ib-s_size+1));
        continue /* XXX NL_MSM_ERROR_RANGE */;
      }

      if (p.ispy) {
        if      (jlo < ja) do { jlo += nj; } while (jlo < ja);
        else if (jlo > jb) do { jlo -= nj; } while (jlo > jb);
      }
      else if (jlo < ja || (jlo+(s_size-1)) > jb) {
        CkPrintf("Error anterpolation:  jlo=%d  out-of-range [%d,%d]\n",
            jlo, ja, (jb-s_size+1));
        continue /* XXX NL_MSM_ERROR_RANGE */;
      }

      if (p.ispz) {
        if      (klo < ka) do { klo += nk; } while (klo < ka);
        else if (klo > kb) do { klo -= nk; } while (klo > kb);
      }
      else if (klo < ka || (klo+(s_size-1)) > kb) {
        CkPrintf("Error anterpolation:  klo=%d  out-of-range [%d,%d]\n",
            klo, ka, (kb-s_size+1));
        continue /* XXX NL_MSM_ERROR_RANGE */;
      }

      // determine charge on cube of grid points around atom, with wrapping
      int i, j, k;
      for (k = 0, kp = klo;  k < s_size;  k++, kp++) {
        if (kp > kb) kp = ka;  /* wrap stencil around grid */
        float ck = zphi[k] * q;
        for (j = 0, jp = jlo;  j < s_size;  j++, jp++) {
          if (jp > jb) jp = ja;  /* wrap stencil around grid */
          float cjk = yphi[j] * ck;
          for (i = 0, ip = ilo;  i < s_size;  i++, ip++) {
            if (ip > ib) ip = ia;  /* wrap stencil around grid */
            qhacc.accumulate(ip,jp,kp) += xphi[i] * cjk;
          }
        }
      }
    } // else
  } // end loop over atoms

  //
  // end of compute work
  //
  qhread = qhacc.syncToRead();  // let consumers read
#ifdef MSM_DEBUG
  CkPrintf("Anterpolation compute %d exiting, PE %d\n", thisIndex, CkMyPe());
#endif
}


void ComputeMsmMsaMgr::interpolate(CkQdMsg *msg)
{
  delete msg;
#ifdef MSM_DEBUG
  CkPrintf("Interpolation compute %d starting, PE %d\n", thisIndex, CkMyPe());
#endif
  ComputeMsmMsaAtom *coord = coordMsg->coord;
  int numAtoms = coordMsg->numAtoms;
  if (numForces < numAtoms) {
    delete [] force;
    force = new MsmMsaForce[numAtoms];
    numForces = numAtoms;
  }
  MsmMsaData &p = msmData;
  if ( ! msa_setup_interpolation) {
    p.eh[0].enroll(p.num_clients_eh[0]);
    ehacc = p.eh[0].getInitialAccum();
    ehread = ehacc.syncToRead();  // immediately sync to read phase
    msa_setup_interpolation = 1;
  }
  int ni = p.grid_len[0].nx;
  int nj = p.grid_len[0].ny;
  int nk = p.grid_len[0].nz;
  int ia = p.grid_idstart[0].nx;
  int ja = p.grid_idstart[0].ny;
  int ka = p.grid_idstart[0].nz;
  int ib = ia + ni - 1;
  int jb = ja + nj - 1;
  int kb = ka + nk - 1;
  ehacc = ehread.syncToEAccum();  // let producers accumulate
  ehread = ehacc.syncToRead();    // then we read it
#ifdef MSM_DEBUG
  CkPrintf("Interpolation compute %d doing work, PE %d\n", thisIndex, CkMyPe());
#endif
#ifdef MSM_DEBUG
  if (thisIndex==0) {
#if 1
    CkPrintf("+++ eh[0,0,0] = %g\n", ehread.get(0,0,0));
#else
    int i, j, k;
    for (k = ka;  k <= kb;  k++) {
      for (j = ja;  j <= jb;  j++) {
        for (i = ia;  i <= ib;  i++) {
          CkPrintf("+++ eh[%d,%d,%d] = %g\n", i, j, k, ehread.get(i,j,k));
        }
      }
    }
#endif
  }
#endif
  //
  // here is the compute work
  //
  double qself = 0;  // accumulate charge for self-energy

  float xphi[MAX_POLY_DEGREE+1];  // Phi stencil along x-dimension
  float yphi[MAX_POLY_DEGREE+1];  // Phi stencil along y-dimension
  float zphi[MAX_POLY_DEGREE+1];  // Phi stencil along z-dimension
  float dxphi[MAX_POLY_DEGREE+1]; // derivative of Phi along x-dimension
  float dyphi[MAX_POLY_DEGREE+1]; // derivative of Phi along y-dimension
  float dzphi[MAX_POLY_DEGREE+1]; // derivative of Phi along z-dimension
  float fx, fy, fz;  // accumulate each force

  int n;
  const int s_size = PolyDegree[p.approx] + 1;         // stencil size
  const int s_edge = (PolyDegree[p.approx] - 1) >> 1;  // stencil "edge" size

  for (n = 0;  n < numAtoms;  n++) {
    // atomic charge
    float q = coord[n].msmcharge;
    if (0==q) {
      force[n].x = 0;
      force[n].y = 0;
      force[n].z = 0;
      continue;
    }
    qself += q*q;

    // distance between atom and origin measured in grid points
    // (does subtraction and division in double, then converts to float)
    float rx_hx = (coord[n].position.x - p.origin_x) * p.hx_1;
    float ry_hy = (coord[n].position.y - p.origin_y) * p.hy_1;
    float rz_hz = (coord[n].position.z - p.origin_z) * p.hz_1;

    // find smallest numbered grid point in stencil
    int ilo = (int) floorf(rx_hx) - s_edge;
    int jlo = (int) floorf(ry_hy) - s_edge;
    int klo = (int) floorf(rz_hz) - s_edge;

    // calculate Phi stencils along each dimension
    float delta;
    delta = rx_hx - (float) ilo;
    D_STENCIL_1D(dxphi, xphi, p.hx_1, delta, p.approx);
    delta = ry_hy - (float) jlo;
    D_STENCIL_1D(dyphi, yphi, p.hy_1, delta, p.approx);
    delta = rz_hz - (float) klo;
    D_STENCIL_1D(dzphi, zphi, p.hz_1, delta, p.approx);

    // test to see if stencil is within edges of grid
    int iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                     ja <= jlo && (jlo+(s_size-1)) <= jb &&
                     ka <= klo && (klo+(s_size-1)) <= kb );

    if ( iswithin ) {  // no wrapping needed
      // determine charge on cube of grid points around atom
      int i, j, k;
      fx = fy = fz = 0;
      for (k = 0;  k < s_size;  k++) {
        for (j = 0;  j < s_size;  j++) {
          float cx = yphi[j] * zphi[k];
          float cy = dyphi[j] * zphi[k];
          float cz = yphi[j] * dzphi[k];
          for (i = 0;  i < s_size;  i++) {
            float e = ehread.get(i+ilo,j+jlo,k+klo);
            fx += e * dxphi[i] * cx;
            fy += e * xphi[i] * cy;
            fz += e * xphi[i] * cz;
          }
        }
      }
    } // if is within
    else {  // requires wrapping around grid
      int ip, jp, kp;

      // adjust ilo, jlo, klo so they are within grid indexing
      if (p.ispx) {
        if      (ilo < ia) do { ilo += ni; } while (ilo < ia);
        else if (ilo > ib) do { ilo -= ni; } while (ilo > ib);
      }
      else if (ilo < ia || (ilo+(s_size-1)) > ib) {
        CkPrintf("Error interpolation:  ilo=%d  out-of-range [%d,%d]\n",
            ilo, ia, (ib-s_size+1));
        continue /* XXX NL_MSM_ERROR_RANGE */;
      }

      if (p.ispy) {
        if      (jlo < ja) do { jlo += nj; } while (jlo < ja);
        else if (jlo > jb) do { jlo -= nj; } while (jlo > jb);
      }
      else if (jlo < ja || (jlo+(s_size-1)) > jb) {
        CkPrintf("Error interpolation:  jlo=%d  out-of-range [%d,%d]\n",
            jlo, ja, (jb-s_size+1));
        continue /* XXX NL_MSM_ERROR_RANGE */;
      }

      if (p.ispz) {
        if      (klo < ka) do { klo += nk; } while (klo < ka);
        else if (klo > kb) do { klo -= nk; } while (klo > kb);
      }
      else if (klo < ka || (klo+(s_size-1)) > kb) {
        CkPrintf("Error interpolation:  klo=%d  out-of-range [%d,%d]\n",
            klo, ka, (kb-s_size+1));
        continue /* XXX NL_MSM_ERROR_RANGE */;
      }

      // determine charge on cube of grid points around atom, with wrapping
      int i, j, k;
      fx = fy = fz = 0;
      for (k = 0, kp = klo;  k < s_size;  k++, kp++) {
        if (kp > kb) kp = ka;  /* wrap stencil around grid */
        for (j = 0, jp = jlo;  j < s_size;  j++, jp++) {
          if (jp > jb) jp = ja;  /* wrap stencil around grid */
          float cx = yphi[j] * zphi[k];
          float cy = dyphi[j] * zphi[k];
          float cz = yphi[j] * dzphi[k];
          for (i = 0, ip = ilo;  i < s_size;  i++, ip++) {
            if (ip > ib) ip = ia;  /* wrap stencil around grid */
            float e = ehread.get(ip,jp,kp);
            fx += e * dxphi[i] * cx;
            fy += e * xphi[i] * cy;
            fz += e * xphi[i] * cz;
          }
        }
      }
    } // else

    // update force
    force[n].x = -q * fx;
    force[n].y = -q * fy;
    force[n].z = -q * fz;

  } // end loop over atoms
  double self_energy = 0.5*p.self_energy_const*qself;

  //
  // end of compute work
  //
#ifdef MSM_DEBUG
  CkPrintf("Interpolation compute %d read sync done, PE %d\n",
      thisIndex, CkMyPe());
#endif
  delete coordMsg;  // get rid of earlier coordMsg
  msmCompute->saveResults(numAtoms, force, self_energy);
#if 0
  int startAtomID = thisIndex * MsmMsa::NUM_ATOMS_PER_CHARE;
  msmProxy.doneForces(startAtomID, nforces, forcebuffer);
#endif
#ifdef MSM_DEBUG
  CkPrintf("Interpolation compute %d exiting, PE %d\n", thisIndex, CkMyPe());
#endif
}


MsmMsaLevel::MsmMsaLevel(MsmMsaGrid &qh, MsmMsaGrid &eh, MsmMsaGrid &q2h, MsmMsaGrid &e2h)
{
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  int level = thisIndex;
  lastlevel = -1;  // this is an intermediate level
  const Int3 &dg = p.dim_gridcutoff_chares[level];
  gridcutoff = CProxy_MsmMsaGridCutoff::ckNew(level, qh, eh, dg.nx, dg.ny, dg.nz);
  const Int3 &dt = p.dim_gridtransfer_chares[level];
  restriction = CProxy_MsmMsaRestriction::ckNew(level, qh, q2h, dt.nx, dt.ny, dt.nz);
  prolongation =CProxy_MsmMsaProlongation::ckNew(level, eh, e2h, dt.nx, dt.ny, dt.nz);
#ifdef MSM_DEBUG
  CkPrintf("Created %d grid cutoff chares\n", p.num_gridcutoff_chares[level]);
#endif
}


MsmMsaLevel::MsmMsaLevel(MsmMsaGrid &qh, MsmMsaGrid &eh)
{
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  int level = thisIndex;
  lastlevel = level;
  const Int3 &dg = p.dim_gridcutoff_chares[level];
  gridcutoff = CProxy_MsmMsaGridCutoff::ckNew(level, qh, eh, dg.nx, dg.ny, dg.nz);
#ifdef MSM_DEBUG
  CkPrintf("Created %d grid cutoff chares\n", p.num_gridcutoff_chares[level]);
#endif
}


void MsmMsaLevel::compute()
{
  if (thisIndex != lastlevel) {
    restriction.compute();
    prolongation.compute();
  }
  gridcutoff.compute();
}


MsmMsaGridCutoff::MsmMsaGridCutoff(int level, MsmMsaGrid &qh_, MsmMsaGrid &eh_)
  : qh(qh_), eh(eh_)
{
  mylevel = level;
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  // find the points on my part of my level's qh and eh grids
  const Int3 &len = p.grid_len[mylevel];
  const Int3 &idstart = p.grid_idstart[mylevel];
  mia = idstart.nx + thisIndex.x * p.num_points_per_chare.nx;
  mja = idstart.ny + thisIndex.y * p.num_points_per_chare.ny;
  mka = idstart.nz + thisIndex.z * p.num_points_per_chare.nz;
  mib = mia + p.num_points_per_chare.nx - 1;
  mjb = mja + p.num_points_per_chare.ny - 1;
  mkb = mka + p.num_points_per_chare.nz - 1;
  if (mib > idstart.nx + len.nx - 1) {
    mib = idstart.nx + len.nx - 1;
  }
  if (mjb > idstart.ny + len.ny - 1) {
    mjb = idstart.ny + len.ny - 1;
  }
  if (mkb > idstart.nz + len.nz - 1) {
    mkb = idstart.nz + len.nz - 1;
  }
  msa_setup = 0;
}


void MsmMsaGridCutoff::compute()
{
#ifdef MSM_DEBUG
  CkPrintf("Grid cutoff compute (%d,%d,%d) PE %d\n",
      thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
#endif
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  if ( ! msa_setup) {
    qh.enroll(p.num_clients_qh[mylevel]);
    eh.enroll(p.num_clients_eh[mylevel]);
    qhacc = qh.getInitialAccum();
    qhread = qhacc.syncToRead();  // immediately sync to read phase
    ehacc = eh.getInitialAccum();
    ehread = ehacc.syncToRead();  // immediately sync to read phase
    msa_setup = 1;
  }
for (;;) {  // loop forever
  qhacc = qhread.syncToEAccum();  // let producers accumulate
  qhread = qhacc.syncToRead();    // then we read charge
  ehacc = ehread.syncToEAccum();  // and we accumulate potential
  //
  // here is the compute work
  //

  int ni = p.grid_len[mylevel].nx;
  int nj = p.grid_len[mylevel].ny;
  int nk = p.grid_len[mylevel].nz;
  int ia = p.grid_idstart[mylevel].nx;
  int ja = p.grid_idstart[mylevel].ny;
  int ka = p.grid_idstart[mylevel].nz;
  int ib = ia + ni - 1;
  int jb = ja + nj - 1;
  int kb = ka + nk - 1;
  int ispnone = ! ( p.ispx || p.ispy || p.ispz );
  int gia, gib, gja, gjb, gka, gkb, gni, gnj;
  const float *gc = 0;
  if (mylevel != p.toplevel) {
    gc = &(p.gc[0]);
    gia = p.gc_idstart.nx;
    gib = gia + p.gc_len.nx - 1;
    gja = p.gc_idstart.ny;
    gjb = gja + p.gc_len.ny - 1;
    gka = p.gc_idstart.nz;
    gkb = gka + p.gc_len.nz - 1;
    gni = p.gc_len.nx;
    gnj = p.gc_len.ny;
  }
  else {
    gc = &(p.gctop[0]);
    gia = p.gctop_idstart.nx;
    gib = gia + p.gctop_len.nx - 1;
    gja = p.gctop_idstart.ny;
    gjb = gja + p.gctop_len.ny - 1;
    gka = p.gctop_idstart.nz;
    gkb = gka + p.gctop_len.nz - 1;
    gni = p.gctop_len.nx;
    gnj = p.gctop_len.ny;
  }
  int i, j, k;
  int gia_clip, gib_clip;
  int gja_clip, gjb_clip;
  int gka_clip, gkb_clip;
  int id, jd, kd;
  int kgoff, jkgoff, ngindex;

  float scaling = p.scaling[mylevel];
  float eh_sum = 0;

#ifdef MSM_DEBUG
  if (mylevel==0 && thisIndex.x==0 && thisIndex.y==0 && thisIndex.z==0) {
    int index = 0;
    CkPrintf("+++ qh[0,0,0] = %g\n", qhread.get(0,0,0));
    CkPrintf("+++ p.toplevel = %d\n", p.toplevel);
    CkPrintf("+++ scaling = %g\n", scaling);
#if 0
    for (k = gka;  k <= gkb;  k++) {
      for (j = gja;  j <= gjb;  j++) {
        for (i = gia;  i <= gib;  i++, index++) {
          printf("+++ gc[%d,%d,%d] = %g\n", i, j, k, gc[index]);
        }
      }
    }
#endif
  }
#endif

  if ( ispnone ) {  // non-periodic boundaries

    // loop over my sub-grid points
    for (k = mka;  k <= mkb;  k++) {

      // clip gc ranges to keep offset for k index within grid
      gka_clip = (k + gka < ka ? ka - k : gka);
      gkb_clip = (k + gkb > kb ? kb - k : gkb);

      for (j = mja;  j <= mjb;  j++) {

        // clip gc ranges to keep offset for j index within grid
        gja_clip = (j + gja < ja ? ja - j : gja);
        gjb_clip = (j + gjb > jb ? jb - j : gjb);

        for (i = mia;  i <= mib;  i++) {

          // clip gc ranges to keep offset for i index within grid
          gia_clip = (i + gia < ia ? ia - i : gia);
          gib_clip = (i + gib > ib ? ib - i : gib);

          // sum over "sphere" of weighted charge
          eh_sum = 0;
          for (kd = gka_clip;  kd <= gkb_clip;  kd++) {
            kgoff = (kd-gka) * gnj;       // find gc flat index

            for (jd = gja_clip;  jd <= gjb_clip;  jd++) {
              jkgoff = (kgoff + (jd-gja)) * gni;       // find gc flat index

              for (id = gia_clip;  id <= gib_clip;  id++) {
                ngindex = jkgoff + (id-gia);       // gc flat index
                // sum weighted charge
                eh_sum += qhread.get(i+id,j+jd,k+kd) * gc[ngindex];
              }
            }
          } // end loop over "sphere" of charge

          // add contribution into MSA
          ehacc.accumulate(i,j,k) += scaling * eh_sum;
        }
      }
    } // end loop over my sub-grid points

  } // if non-periodic boundaries
  else {
    // some boundary is periodic
    int ilo, jlo, klo;
    int ip, jp, kp;

    // loop over my sub-grid points
    for (k = mka;  k <= mkb;  k++) {
      klo = k + gka;
      if ( ! p.ispz ) {  // non-periodic z
        // clip gc ranges to keep offset for k index within grid
        gka_clip = (k + gka < ka ? ka - k : gka);
        gkb_clip = (k + gkb > kb ? kb - k : gkb);
        if (klo < ka) klo = ka;  // keep lowest qh index within grid
      }
      else {  // periodic z
        gka_clip = gka;
        gkb_clip = gkb;
        if (klo < ka) do { klo += nk; } while (klo < ka);
      }
      // ASSERT(klo <= kb);

      for (j = mja;  j <= mjb;  j++) {
        jlo = j + gja;
        if ( ! p.ispy ) {  // non-periodic y
          // clip gc ranges to keep offset for j index within grid
          gja_clip = (j + gja < ja ? ja - j : gja);
          gjb_clip = (j + gjb > jb ? jb - j : gjb);
          if (jlo < ja) jlo = ja;  // keep lowest qh index within grid
        }
        else {  // periodic y
          gja_clip = gja;
          gjb_clip = gjb;
          if (jlo < ja) do { jlo += nj; } while (jlo < ja);
        }
        // ASSERT(jlo <= jb);

        for (i = mia;  i <= mib;  i++) {
          ilo = i + gia;
          if ( ! p.ispx ) {  // nonperiodic x
            // clip gc ranges to keep offset for i index within grid
            gia_clip = (i + gia < ia ? ia - i : gia);
            gib_clip = (i + gib > ib ? ib - i : gib);
            if (ilo < ia) ilo = ia;  // keep lowest qh index within grid
          }
          else {  // periodic x
            gia_clip = gia;
            gib_clip = gib;
            if (ilo < ia) do { ilo += ni; } while (ilo < ia);
          }
          // ASSERT(ilo <= ib);

          // sum over "sphere" of weighted charge
          eh_sum = 0;
          for (kd = gka_clip, kp = klo;  kd <= gkb_clip;  kd++, kp++) {
            // clipping makes conditional always fail for nonperiodic
            if (kp > kb) kp = ka;  // wrap z direction
            kgoff = (kd-gka) * gnj;      // find gc flat index

            for (jd = gja_clip, jp = jlo;  jd <= gjb_clip;  jd++, jp++) {
              // clipping makes conditional always fail for nonperiodic
              if (jp > jb) jp = ja;              // wrap y direction
              jkgoff = (kgoff + (jd-gja)) * gni;       // find gc flat index */

              for (id = gia_clip, ip = ilo;  id <= gib_clip;  id++, ip++) {
                // clipping makes conditional always fail for nonperiodic
                if (ip > ib) ip = ia;   // wrap x direction
                ngindex = jkgoff + (id-gia);  // gc flat index
                // sum weighted charge
                eh_sum += qhread.get(ip,jp,kp) * gc[ngindex];
              }
            }
          } // end loop over "sphere" of charge

          // add contribution into MSA
          ehacc.accumulate(i,j,k) += scaling * eh_sum;
        }
      }
    } // end loop over my sub-grid points

  } // else some boundary is periodic

  //
  // end of compute work
  //
  ehread = ehacc.syncToRead();  // let consumers read potential
#ifdef MSM_DEBUG
  CkPrintf("Grid cutoff compute %d exiting, PE %d\n", thisIndex, CkMyPe());
#endif
} // loop forever
}


/** Max stencil length is basically PolyDegree+2 for those approximations
 * that interpolate.  (We skip over zero in the complete stencils above.) */
enum { MAX_NSTENCIL = 11 };

/** The stencil array lengths below. */
static const int Nstencil[APPROX_END] = {
  5, 7, 7, 9, 9, 11, 11, 7
};

/** Index offsets from the stencil-centered grid element, to get
 * to the correct contributing grid element. */
static const int IndexOffset[APPROX_END][MAX_NSTENCIL] = {
  /* cubic */
  {-3, -1, 0, 1, 3},

  /* quintic C1 */
  {-5, -3, -1, 0, 1, 3, 5},

  /* quintic C2  (same as quintic C1) */
  {-5, -3, -1, 0, 1, 3, 5},

  /* septic C1 */
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  /* septic C3  (same as septic C3) */
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  /* nonic C1 */
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  /* nonic C4  (same as nonic C1) */
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  /* bspline */
  {-3, -2, -1, 0, 1, 2, 3},
};

/** The grid transfer stencils for the non-factored restriction and
 * prolongation procedures. */
static const float PhiStencil[APPROX_END][MAX_NSTENCIL] = {
  /* cubic */
  {-1.f/16, 9.f/16, 1, 9.f/16, -1.f/16},

  /* quintic C1 */
  {3.f/256, -25.f/256, 75.f/128, 1, 75.f/128, -25.f/256, 3.f/256},

  /* quintic C2  (same as quintic C1) */
  {3.f/256, -25.f/256, 75.f/128, 1, 75.f/128, -25.f/256, 3.f/256},

  /* septic C1 */
  { -5.f/2048, 49.f/2048, -245.f/2048, 1225.f/2048, 1, 1225.f/2048,
    -245.f/2048, 49.f/2048, -5.f/2048 },

  /* septic C3  (same as septic C3) */
  { -5.f/2048, 49.f/2048, -245.f/2048, 1225.f/2048, 1, 1225.f/2048,
    -245.f/2048, 49.f/2048, -5.f/2048 },

  /* nonic C1 */
  { 35.f/65536, -405.f/65536, 567.f/16384, -2205.f/16384, 
    19845.f/32768, 1, 19845.f/32768, -2205.f/16384, 567.f/16384, 
    -405.f/65536, 35.f/65536 },

  /* nonic C4  (same as nonic C1) */
  { 35.f/65536, -405.f/65536, 567.f/16384, -2205.f/16384, 
    19845.f/32768, 1, 19845.f/32768, -2205.f/16384, 567.f/16384, 
    -405.f/65536, 35.f/65536 },

  /* bspline */
  { 1.f/48, 1.f/6, 23.f/48, 2.f/3, 23.f/48, 1.f/6, 1.f/48 },
};


MsmMsaRestriction::MsmMsaRestriction(int level, MsmMsaGrid &qh_, MsmMsaGrid &q2h_)
  : qh(qh_), q2h(q2h_)
{
  mylevel = level;
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  // find the points on my sub-grid of (mylevel+1) grid
  const Int3 &len = p.grid_len[mylevel+1];
  const Int3 &idstart = p.grid_idstart[mylevel+1];
  mia = idstart.nx + thisIndex.x * p.num_points_per_chare.nx;
  mja = idstart.ny + thisIndex.y * p.num_points_per_chare.ny;
  mka = idstart.nz + thisIndex.z * p.num_points_per_chare.nz;
  mib = mia + p.num_points_per_chare.nx - 1;
  mjb = mja + p.num_points_per_chare.ny - 1;
  mkb = mka + p.num_points_per_chare.nz - 1;
  if (mib > idstart.nx + len.nx - 1) {
    mib = idstart.nx + len.nx - 1;
  }
  if (mjb > idstart.ny + len.ny - 1) {
    mjb = idstart.ny + len.ny - 1;
  }
  if (mkb > idstart.nz + len.nz - 1) {
    mkb = idstart.nz + len.nz - 1;
  }
  msa_setup = 0;
}


void MsmMsaRestriction::compute()
{
#ifdef MSM_DEBUG
  CkPrintf("Restriction compute %d, PE %d\n", thisIndex, CkMyPe());
#endif
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  if ( ! msa_setup) {
    qh.enroll(p.num_clients_qh[mylevel]);
    q2h.enroll(p.num_clients_qh[mylevel+1]);
    qhacc = qh.getInitialAccum();
    qhread = qhacc.syncToRead();  // immediately sync to read phase
    q2hacc = q2h.getInitialAccum();
    q2hread = q2hacc.syncToRead();  // immediately sync to read phase
    msa_setup = 1;
  }
for (;;) {  // loop forever
  qhacc = qhread.syncToEAccum();  // let producers accumulate h-level charge
  qhread = qhacc.syncToRead();  // then we read h-level charge
  q2hacc = q2hread.syncToEAccum();  // and we accumulate 2h-level charge
#ifdef MSM_DEBUG
  CkPrintf("Restriction compute %d doing work, PE %d\n", thisIndex, CkMyPe());
#endif
  //
  // here is the compute work
  //

  const int nstencil = Nstencil[p.approx];
  const int *offset = IndexOffset[p.approx];
  const float *phi = PhiStencil[p.approx];

  int ni1 = p.grid_len[mylevel].nx;
  int nj1 = p.grid_len[mylevel].ny;
  int nk1 = p.grid_len[mylevel].nz;
  int ia1 = p.grid_idstart[mylevel].nx;
  int ja1 = p.grid_idstart[mylevel].ny;
  int ka1 = p.grid_idstart[mylevel].nz;
  int ib1 = ia1 + ni1 - 1;
  int jb1 = ja1 + nj1 - 1;
  int kb1 = ka1 + nk1 - 1;
  int i, j, k, i1, j1, k1, i2, j2, k2;
  int i1off, j1off, k1off;

  float q2h_sum, cjk;

  for (k2 = mka;  k2 <= mkb;  k2++) {
    k1 = k2 * 2;
    for (j2 = mja;  j2 <= mjb;  j2++) {
      j1 = j2 * 2;
      for (i2 = mia;  i2 <= mib;  i2++) {
        i1 = i2 * 2;

        q2h_sum = 0;
        for (k = 0;  k < nstencil;  k++) {
          k1off = k1 + offset[k];
          if (k1off < ka1) {
            if (p.ispz) do { k1off += nk1; } while (k1off < ka1);
            else continue;
          }
          else if (k1off > kb1) {
            if (p.ispz) do { k1off -= nk1; } while (k1off > kb1);
            else break;
          }
          for (j = 0;  j < nstencil;  j++) {
            j1off = j1 + offset[j];
            if (j1off < ja1) {
              if (p.ispy) do { j1off += nj1; } while (j1off < ja1);
              else continue;
            }
            else if (j1off > jb1) {
              if (p.ispy) do { j1off -= nj1; } while (j1off > jb1);
              else break;
            }
            cjk = phi[j] * phi[k];
            for (i = 0;  i < nstencil;  i++) {
              i1off = i1 + offset[i];
              if (i1off < ia1) {
                if (p.ispx) do { i1off += ni1; } while (i1off < ia1);
                else continue;
              }
              else if (i1off > ib1) {
                if (p.ispx) do { i1off -= ni1; } while (i1off > ib1);
                else break;
              }
              q2h_sum += qhread.get(i1off,j1off,k1off) * phi[i] * cjk;
            }
          }
        } // end loop over finer grid stencil

        q2hacc.accumulate(i2,j2,k2) += q2h_sum;
      }
    }
  } // end loop over each coarser grid point

  //
  // end of compute work
  //
  q2hread = q2hacc.syncToRead();  // let consumers read 2h-level charge
#ifdef MSM_DEBUG
  CkPrintf("Restriction compute %d exiting, PE %d\n", thisIndex, CkMyPe());
#endif
} // loop forever
}


MsmMsaProlongation::MsmMsaProlongation(int level, MsmMsaGrid &eh_, MsmMsaGrid &e2h_) 
  : eh(eh_), e2h(e2h_)
{
  mylevel = level;
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  // find the points on my sub-grid of (mylevel+1) grid
  const Int3 &len = p.grid_len[mylevel+1];
  const Int3 &idstart = p.grid_idstart[mylevel+1];
  mia = idstart.nx + thisIndex.x * p.num_points_per_chare.nx;
  mja = idstart.ny + thisIndex.y * p.num_points_per_chare.ny;
  mka = idstart.nz + thisIndex.z * p.num_points_per_chare.nz;
  mib = mia + p.num_points_per_chare.nx - 1;
  mjb = mja + p.num_points_per_chare.ny - 1;
  mkb = mka + p.num_points_per_chare.nz - 1;
  if (mib > idstart.nx + len.nx - 1) {
    mib = idstart.nx + len.nx - 1;
  }
  if (mjb > idstart.ny + len.ny - 1) {
    mjb = idstart.ny + len.ny - 1;
  }
  if (mkb > idstart.nz + len.nz - 1) {
    mkb = idstart.nz + len.nz - 1;
  }
  msa_setup = 0;
}


void MsmMsaProlongation::compute()
{
#ifdef MSM_DEBUG
  CkPrintf("Prolongation compute %d, PE %d\n", thisIndex, CkMyPe());
#endif
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  if ( ! msa_setup) {
    eh.enroll(p.num_clients_eh[mylevel]);
    e2h.enroll(p.num_clients_eh[mylevel+1]);
    ehacc = eh.getInitialAccum();
    ehread = ehacc.syncToRead();  // immediately sync to read phase
    e2hacc = e2h.getInitialAccum();
    e2hread = e2hacc.syncToRead();  // immediately sync to read phase
    msa_setup = 1;
  }
for (;;) {  // loop forever
  ehacc = ehread.syncToEAccum();  // we accumulate h-level potential
  e2hacc = e2hread.syncToEAccum();  // let producers accumulate 2h-level
  e2hread = e2hacc.syncToRead();  // then we read 2h-level potentials
#ifdef MSM_DEBUG
  CkPrintf("Prolongation compute %d doing work, PE %d\n", thisIndex, CkMyPe());
#endif
  //
  // here is the compute work
  //

  const int nstencil = Nstencil[p.approx];
  const int *offset = IndexOffset[p.approx];
  const float *phi = PhiStencil[p.approx];

  int ni1 = p.grid_len[mylevel].nx;
  int nj1 = p.grid_len[mylevel].ny;
  int nk1 = p.grid_len[mylevel].nz;
  int ia1 = p.grid_idstart[mylevel].nx;
  int ja1 = p.grid_idstart[mylevel].ny;
  int ka1 = p.grid_idstart[mylevel].nz;
  int ib1 = ia1 + ni1 - 1;
  int jb1 = ja1 + nj1 - 1;
  int kb1 = ka1 + nk1 - 1;
  int i, j, k, i1, j1, k1, i2, j2, k2;
  int i1off, j1off, k1off;

  float cjk;

  for (k2 = mka;  k2 <= mkb;  k2++) {
    k1 = k2 * 2;
    for (j2 = mja;  j2 <= mjb;  j2++) {
      j1 = j2 * 2;
      for (i2 = mia;  i2 <= mib;  i2++) {
        i1 = i2 * 2;

        for (k = 0;  k < nstencil;  k++) {
          k1off = k1 + offset[k];
          if (k1off < ka1) {
            if (p.ispz) do { k1off += nk1; } while (k1off < ka1);
            else continue;
          }
          else if (k1off > kb1) {
            if (p.ispz) do { k1off -= nk1; } while (k1off > kb1);
            else break;
          }
          for (j = 0;  j < nstencil;  j++) {
            j1off = j1 + offset[j];
            if (j1off < ja1) {
              if (p.ispy) do { j1off += nj1; } while (j1off < ja1);
              else continue;
            }
            else if (j1off > jb1) {
              if (p.ispy) do { j1off -= nj1; } while (j1off > jb1);
              else break;
            }
            cjk = phi[j] * phi[k];
            for (i = 0;  i < nstencil;  i++) {
              i1off = i1 + offset[i];
              if (i1off < ia1) {
                if (p.ispx) do { i1off += ni1; } while (i1off < ia1);
                else continue;
              }
              else if (i1off > ib1) {
                if (p.ispx) do { i1off -= ni1; } while (i1off > ib1);
                else break;
              }
              ehacc.accumulate(i1off,j1off,k1off) +=
                e2hread(i2,j2,k2) * phi[i] * cjk;
            }
          }
        } // end loop over finer grid stencil

      }
    }
  } // end loop over each coarser grid point

  //
  // end of compute work
  //
  ehread = ehacc.syncToRead();  // let consumers read h-level potential
#ifdef MSM_DEBUG
  CkPrintf("Prolongation compute %d exiting, PE %d\n", thisIndex, CkMyPe());
#endif
} // loop forever
}


MsmMsaEnergy::MsmMsaEnergy(MsmMsaGrid &qh_, MsmMsaGrid &eh_) : qh(qh_), eh(eh_)
{
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  int npoints = p.num_points_per_chare.nx * 
    p.num_points_per_chare.ny * p.num_points_per_chare.nz;
  qhbuffer = new float[npoints];
  int ni = p.grid_len[0].nx;
  int nj = p.grid_len[0].ny;
  int nk = p.grid_len[0].nz;
  int nci = ROUNDUP_QUOTIENT(ni, p.num_points_per_chare.nx);
  int ncj = ROUNDUP_QUOTIENT(nj, p.num_points_per_chare.ny);
  int nck = ROUNDUP_QUOTIENT(nk, p.num_points_per_chare.nz);
  mlk = thisIndex / (ncj * nci);
  int krem = thisIndex % (ncj * nci);
  mlj = krem / nci;
  mli = krem % nci;
  // find the points on my part of my level's qh and eh grids
  mia = p.grid_idstart[0].nx + mli * p.num_points_per_chare.nx;
  mja = p.grid_idstart[0].ny + mlj * p.num_points_per_chare.ny;
  mka = p.grid_idstart[0].nz + mlk * p.num_points_per_chare.nz;
  mib = mia + p.num_points_per_chare.nx - 1;
  mjb = mja + p.num_points_per_chare.ny - 1;
  mkb = mka + p.num_points_per_chare.nz - 1;
  if (mib > p.grid_idstart[0].nx + ni - 1) {
    mib = p.grid_idstart[0].nx + ni - 1;
  }
  if (mjb > p.grid_idstart[0].ny + nj - 1) {
    mjb = p.grid_idstart[0].ny + nj - 1;
  }
  if (mkb > p.grid_idstart[0].nz + nk - 1) {
    mkb = p.grid_idstart[0].nz + nk - 1;
  }
  msa_setup = 0;
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

MsmMsaEnergy::~MsmMsaEnergy()
{
  delete[] qhbuffer;
}


void MsmMsaEnergy::compute()
{
#ifdef MSM_DEBUG
  CkPrintf("Energy compute %d, PE %d\n", thisIndex, CkMyPe());
#endif
  const MsmMsaData &p = CProxy_ComputeMsmMsaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMsaMgr)->getMsmMsaData();
  if ( ! msa_setup) {
    qh.enroll(p.num_clients_qh[0]);
    eh.enroll(p.num_clients_eh[0]);
    qhacc = qh.getInitialAccum();
    qhread = qhacc.syncToRead();  // immediately sync to read phase
    ehacc = eh.getInitialAccum();
    ehread = ehacc.syncToRead();
    msa_setup = 1;
  }
for (;;) {  // loop forever
  double sum = 0;
  qhacc = qhread.syncToEAccum();  // let producers accumulate charge
  qhread = qhacc.syncToRead();  // then we read it
  //
  // Must be very careful of phase sync order between eh and qh
  // qhread sync done must be after enrollment in eh
  // (otherwise block gridcutoff)
  // but before ehacc sync to read
  // (otherwise block prolongation)
  //
  // Since we can't have both qh and eh open for reading at once,
  // we read qh values into local storage
  //
  int i, j, k, index;
  for (index = 0, k = mka;  k <= mkb;  k++) {
    for (j = mja;  j <= mjb;  j++) {
      for (i = mia;  i <= mib;  i++, index++) {
        qhbuffer[index] = qhread.get(i,j,k);
      }
    }
  }
  //
  // Must do qhread sync done BEFORE ehacc sync to read
  //
  ehacc = ehread.syncToEAccum();  // let producers accumulate potential
  ehread = ehacc.syncToRead();  // then we read it
#ifdef MSM_DEBUG
  CkPrintf("Energy compute %d doing work, PE %d\n", thisIndex, CkMyPe());
#endif
  for (index = 0, k = mka;  k <= mkb;  k++) {
    for (j = mja;  j <= mjb;  j++) {
      for (i = mia;  i <= mib;  i++, index++) {
        sum += qhbuffer[index] * ehread.get(i,j,k);
      }
    }
  }
  //contribute(sizeof(double), &sum, CkReduction::sum_double);
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += 0.5*sum;
  reduction->submit();
#ifdef MSM_DEBUG
  CkPrintf("Energy compute %d exiting, PE %d\n", thisIndex, CkMyPe());
#endif
} // loop forever
}


#include "ComputeMsmMsaMgr.def.h"

#endif // CHARM_HAS_MSA

