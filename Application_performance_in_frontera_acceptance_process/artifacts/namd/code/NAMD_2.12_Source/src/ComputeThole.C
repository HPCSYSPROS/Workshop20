/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeThole.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"

#define CALCULATE_THOLE_CORRECTION

// static initialization
int TholeElem::pressureProfileSlabs = 0;
int TholeElem::pressureProfileAtomTypes = 1;
BigReal TholeElem::pressureProfileThickness = 0;
BigReal TholeElem::pressureProfileMin = 0;

void TholeElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Thole** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in TholeElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numTholes;
  *byatom = mol->tholesByAtom;
  *structarray = mol->tholes;
#endif
}

void TholeElem::getParameterPointers(Parameters *p, const TholeValue **v) {
  *v = NULL;  // parameters are stored in the structure
}

void TholeElem::computeForce(TholeElem *tuples, int ntuple, BigReal *reduction, 
                                BigReal *pressureProfileData)
{
 const Lattice & lattice = tuples[0].p[0]->p->lattice;

 for ( int ituple=0; ituple<ntuple; ++ituple ) {
  const TholeElem &tup = tuples[ituple];
  enum { size = 4 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const TholeValue * const(&value)(tup.value);

  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << " "
               << localIndex[3] << std::endl);

#ifdef CALCULATE_THOLE_CORRECTION
  const BigReal aa = value->aa;
  const BigReal qq = value->qq;

  //  Calculate the vectors between atoms
  const Position & rai = p[0]->x[localIndex[0]].position;  // atom i
  const Position & rdi = p[1]->x[localIndex[1]].position;  // atom i's Drude
  const Position & raj = p[2]->x[localIndex[2]].position;  // atom j
  const Position & rdj = p[3]->x[localIndex[3]].position;  // atom j's Drude

  // r_ij = r_i - r_j
  Vector raa = lattice.delta(rai,raj);  // shortest vector image:  rai - raj
  Vector rad = lattice.delta(rai,rdj);  // shortest vector image:  rai - rdj
  Vector rda = lattice.delta(rdi,raj);  // shortest vector image:  rdi - raj
  Vector rdd = lattice.delta(rdi,rdj);  // shortest vector image:  rdi - rdj

  // 1/r, r = |r_ij|
  BigReal raa_invlen = raa.rlength();  // reciprocal of length
  BigReal rad_invlen = rad.rlength();
  BigReal rda_invlen = rda.rlength();
  BigReal rdd_invlen = rdd.rlength();

  // ar
  BigReal auaa = aa / raa_invlen;
  BigReal auad = aa / rad_invlen;
  BigReal auda = aa / rda_invlen;
  BigReal audd = aa / rdd_invlen;

  // exp(-ar)
  BigReal expauaa = exp(-auaa);
  BigReal expauad = exp(-auad);
  BigReal expauda = exp(-auda);
  BigReal expaudd = exp(-audd);

  // (1 + ar/2)
  BigReal polyauaa = 1 + 0.5*auaa;
  BigReal polyauad = 1 + 0.5*auad;
  BigReal polyauda = 1 + 0.5*auda;
  BigReal polyaudd = 1 + 0.5*audd;

  // U(r) = qq/r (1 - (1 + ar/2) exp(-ar))
  BigReal ethole = 0;
  ethole += qq * raa_invlen * (1 - polyauaa * expauaa);
  ethole += -qq * rad_invlen * (1 - polyauad * expauad);
  ethole += -qq * rda_invlen * (1 - polyauda * expauda);
  ethole += qq * rdd_invlen * (1 - polyaudd * expaudd);

  polyauaa = 1 + auaa*polyauaa;
  polyauad = 1 + auad*polyauad;
  polyauda = 1 + auda*polyauda;
  polyaudd = 1 + audd*polyaudd;

  BigReal raa_invlen3 = raa_invlen * raa_invlen * raa_invlen;
  BigReal rad_invlen3 = rad_invlen * rad_invlen * rad_invlen;
  BigReal rda_invlen3 = rda_invlen * rda_invlen * rda_invlen;
  BigReal rdd_invlen3 = rdd_invlen * rdd_invlen * rdd_invlen;

  // df = (1/r) (dU/dr)
  BigReal dfaa = qq * raa_invlen3 * (polyauaa*expauaa - 1);
  BigReal dfad = -qq * rad_invlen3 * (polyauad*expauad - 1);
  BigReal dfda = -qq * rda_invlen3 * (polyauda*expauda - 1);
  BigReal dfdd = qq * rdd_invlen3 * (polyaudd*expaudd - 1);

  Vector faa = -dfaa * raa;
  Vector fad = -dfad * rad;
  Vector fda = -dfda * rda;
  Vector fdd = -dfdd * rdd;

  p[0]->f[localIndex[0]] += faa + fad;
  p[1]->f[localIndex[1]] += fda + fdd;
  p[2]->f[localIndex[2]] -= faa + fda;
  p[3]->f[localIndex[3]] -= fad + fdd;

  DebugM(3, "::computeForce() -- ending with delta energy " << ethole
      << std::endl);
  reduction[tholeEnergyIndex] += ethole;

  reduction[virialIndex_XX] += faa.x * raa.x + fad.x * rad.x
    + fda.x * rda.x + fdd.x * rdd.x;
  reduction[virialIndex_XY] += faa.x * raa.y + fad.x * rad.y
    + fda.x * rda.y + fdd.x * rdd.y;
  reduction[virialIndex_XZ] += faa.x * raa.z + fad.x * rad.z
    + fda.x * rda.z + fdd.x * rdd.z;
  reduction[virialIndex_YX] += faa.y * raa.x + fad.y * rad.x
    + fda.y * rda.x + fdd.y * rdd.x;
  reduction[virialIndex_YY] += faa.y * raa.y + fad.y * rad.y
    + fda.y * rda.y + fdd.y * rdd.y;
  reduction[virialIndex_YZ] += faa.y * raa.z + fad.y * rad.z
    + fda.y * rda.z + fdd.y * rdd.z;
  reduction[virialIndex_ZX] += faa.z * raa.x + fad.z * rad.x
    + fda.z * rda.x + fdd.z * rdd.x;
  reduction[virialIndex_ZY] += faa.z * raa.y + fad.z * rad.y
    + fda.z * rda.y + fdd.z * rdd.y;
  reduction[virialIndex_ZZ] += faa.z * raa.z + fad.z * rad.z
    + fda.z * rda.z + fdd.z * rdd.z;

  if (pressureProfileData) {
    BigReal zai = p[0]->x[localIndex[0]].position.z;
    BigReal zdi = p[1]->x[localIndex[1]].position.z;
    BigReal zaj = p[2]->x[localIndex[2]].position.z;
    BigReal zdj = p[3]->x[localIndex[3]].position.z;
    int nai = (int)floor((zai-pressureProfileMin)/pressureProfileThickness);
    int ndi = (int)floor((zdi-pressureProfileMin)/pressureProfileThickness);
    int naj = (int)floor((zaj-pressureProfileMin)/pressureProfileThickness);
    int ndj = (int)floor((zdj-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(nai, pressureProfileSlabs);
    pp_clamp(ndi, pressureProfileSlabs);
    pp_clamp(naj, pressureProfileSlabs);
    pp_clamp(ndj, pressureProfileSlabs);
    int pai = p[0]->x[localIndex[0]].partition;
    int pdi = p[1]->x[localIndex[1]].partition;
    int paj = p[2]->x[localIndex[2]].partition;
    int pdj = p[3]->x[localIndex[3]].partition;
    int pn = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs, nai, naj,
        pai, paj, pn, faa.x * raa.x, faa.y * raa.y, faa.z * raa.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, nai, ndj,
        pai, pdj, pn, fad.x * rad.x, fad.y * rad.y, fad.z * rad.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, ndi, naj,
        pdi, paj, pn, fda.x * rda.x, fda.y * rda.y, fda.z * rda.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, ndi, ndj,
        pdi, pdj, pn, fdd.x * rdd.x, fdd.y * rdd.y, fdd.z * rdd.z,
        pressureProfileData);
  }
#endif

 }
}


// The energy from the screened Coulomb correction of Thole is 
// accumulated into the electrostatic potential energy.
void TholeElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ELECT_ENERGY) += data[tholeEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

