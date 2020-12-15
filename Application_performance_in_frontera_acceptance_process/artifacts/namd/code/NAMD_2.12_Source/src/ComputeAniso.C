/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeAniso.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"

#define CALCULATE_ANISO

// static initialization
int AnisoElem::pressureProfileSlabs = 0;
int AnisoElem::pressureProfileAtomTypes = 1;
BigReal AnisoElem::pressureProfileThickness = 0;
BigReal AnisoElem::pressureProfileMin = 0;

void AnisoElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Aniso** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in AnisoElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numAnisos;
  *byatom = mol->anisosByAtom;
  *structarray = mol->anisos;
#endif
}

void AnisoElem::getParameterPointers(Parameters *p, const AnisoValue **v) {
  *v = NULL;  // parameters are stored in the structure
}

void AnisoElem::computeForce(AnisoElem *tuples, int ntuple, BigReal *reduction, 
                                BigReal *pressureProfileData)
{
 const Lattice & lattice = tuples[0].p[0]->p->lattice;

 for ( int ituple=0; ituple<ntuple; ++ituple ) {
  const AnisoElem &tup = tuples[ituple];
  enum { size = 4 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const AnisoValue * const(&value)(tup.value);

  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << " "
               << localIndex[3] << std::endl);

#ifdef CALCULATE_ANISO
  // used some comments from Ed Harder's implementation in CHARMM

  const BigReal kpar0  = 2*value->k11;  // force constants
  const BigReal kperp0 = 2*value->k22;
  const BigReal kiso0  = 2*value->k33;

  const Position & ri = p[0]->x[localIndex[0]].position;    // atom I
  const Position & rj = p[0]->x[localIndex[0]+1].position;  // atom I's Drude
  const Position & rl = p[1]->x[localIndex[1]].position;    // atom L
  const Position & rm = p[2]->x[localIndex[2]].position;    // atom M
  const Position & rn = p[3]->x[localIndex[3]].position;    // atom N

  // calculate parallel and perpendicular displacement vectors
  Vector r_il = lattice.delta(ri,rl);  // shortest vector image:  ri - rl
  Vector r_mn = lattice.delta(rm,rn);  // shortest vector image:  rm - rn

  BigReal r_il_invlen = r_il.rlength();  // need recip lengths of r_il, r_mn
  BigReal r_mn_invlen = r_mn.rlength();

  Vector u1 = r_il * r_il_invlen;  // normalize r_il, r_mn
  Vector u2 = r_mn * r_mn_invlen;

  Vector dr = rj - ri;  // Drude displacement vector (ri, rj are in same patch)

  BigReal dpar  = dr * u1;  // parallel displacement
  BigReal dperp = dr * u2;  // perpendicular displacement

  // aniso spring energy
  // kpar reduces response along carbonyl vector
  // kperp reduces response perp to bond vector
  //   (reg in and out of plane response)
  BigReal eaniso;
  eaniso = 0.5*kpar0*dpar*dpar + 0.5*kperp0*dperp*dperp + 0.5*kiso0*(dr*dr);

  // calculate force vectors in one direction only:
  // fi = -(fj + fl),  fn = -fm

  // force on atom j
  Vector fj = -kiso0 * dr;
  fj -= kpar0 * dpar * u1;
  fj -= kperp0 * dperp * u2;

  // force on atom l
  Vector fl = kpar0 * dpar * r_il_invlen * dr;
  fl -= kpar0 * dpar * dpar * r_il_invlen * u1;

  // force on atom m
  Vector fm = kperp0 * dperp * dperp * r_mn_invlen * u2;
  fm -= kperp0 * dperp * r_mn_invlen * dr;

  // accumulate forces
  p[0]->f[localIndex[0]] -= (fj + fl);
  p[0]->f[localIndex[0]+1] += fj;
  p[1]->f[localIndex[1]] += fl;
  p[2]->f[localIndex[2]] += fm;
  p[3]->f[localIndex[3]] -= fm;

  // update potential
  reduction[anisoEnergyIndex] += eaniso;

  // update virial
  reduction[virialIndex_XX] += fj.x * dr.x - fl.x * r_il.x + fm.x * r_mn.x;
  reduction[virialIndex_XY] += fj.x * dr.y - fl.x * r_il.y + fm.x * r_mn.y;
  reduction[virialIndex_XZ] += fj.x * dr.z - fl.x * r_il.z + fm.x * r_mn.z;
  reduction[virialIndex_YX] += fj.y * dr.x - fl.y * r_il.x + fm.y * r_mn.x;
  reduction[virialIndex_YY] += fj.y * dr.y - fl.y * r_il.y + fm.y * r_mn.y;
  reduction[virialIndex_YZ] += fj.y * dr.z - fl.y * r_il.z + fm.y * r_mn.z;
  reduction[virialIndex_ZX] += fj.z * dr.x - fl.z * r_il.x + fm.z * r_mn.x;
  reduction[virialIndex_ZY] += fj.z * dr.y - fl.z * r_il.y + fm.z * r_mn.y;
  reduction[virialIndex_ZZ] += fj.z * dr.z - fl.z * r_il.z + fm.z * r_mn.z;

  // update pressure profile data
  if (pressureProfileData) {
    BigReal zi = p[0]->x[localIndex[0]].position.z;
    BigReal zj = p[0]->x[localIndex[0]+1].position.z;
    BigReal zl = p[1]->x[localIndex[1]].position.z;
    BigReal zm = p[2]->x[localIndex[2]].position.z;
    BigReal zn = p[3]->x[localIndex[3]].position.z;
    int ni = (int)floor((zi-pressureProfileMin)/pressureProfileThickness);
    int nj = (int)floor((zj-pressureProfileMin)/pressureProfileThickness);
    int nl = (int)floor((zl-pressureProfileMin)/pressureProfileThickness);
    int nm = (int)floor((zm-pressureProfileMin)/pressureProfileThickness);
    int nn = (int)floor((zn-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(ni, pressureProfileSlabs);
    pp_clamp(nj, pressureProfileSlabs);
    pp_clamp(nl, pressureProfileSlabs);
    pp_clamp(nm, pressureProfileSlabs);
    pp_clamp(nn, pressureProfileSlabs);
    int pi = p[0]->x[localIndex[0]].partition;
    int pj = p[0]->x[localIndex[0]+1].partition;
    int pl = p[1]->x[localIndex[1]].partition;
    int pm = p[2]->x[localIndex[2]].partition;
    int pn = p[3]->x[localIndex[3]].partition;
    int pt = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs, nj, ni,
        pj, pi, pt, fj.x * dr.x, fj.y * dr.y, fj.z * dr.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, ni, nl,
        pi, pl, pt, -fl.x * r_il.x, -fl.y * r_il.y, -fl.z * r_il.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, nm, nn,
        pm, pn, pt, fm.x * r_mn.x, fm.y * r_mn.y, fm.z * r_mn.z,
        pressureProfileData);
  }
#endif

 }
}


void AnisoElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_BOND_ENERGY) += data[anisoEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

