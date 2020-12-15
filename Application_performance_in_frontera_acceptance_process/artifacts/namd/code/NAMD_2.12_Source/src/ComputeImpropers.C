/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeImpropers.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"


// static initialization
int ImproperElem::pressureProfileSlabs = 0;
int ImproperElem::pressureProfileAtomTypes = 1;
BigReal ImproperElem::pressureProfileThickness = 0;
BigReal ImproperElem::pressureProfileMin = 0;

void ImproperElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Improper** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in ImproperElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numImpropers;
  *byatom = mol->impropersByAtom;
  *structarray = mol->impropers;
#endif
}

void ImproperElem::getParameterPointers(Parameters *p, const ImproperValue **v) {
  *v = p->improper_array;
}

void ImproperElem::computeForce(ImproperElem *tuples, int ntuple, BigReal *reduction,
                                BigReal *pressureProfileData)
{
 const Lattice & lattice = tuples[0].p[0]->p->lattice;

 //fepb BKR
 SimParameters *const simParams = Node::Object()->simParameters;
 const int step = tuples[0].p[0]->p->flags.step;
 const BigReal alchLambda = simParams->getCurrentLambda(step);
 const BigReal alchLambda2 = simParams->alchLambda2;
 const BigReal bond_lambda_1 = simParams->getBondLambda(alchLambda);
 const BigReal bond_lambda_2 = simParams->getBondLambda(1-alchLambda);
 const BigReal bond_lambda_12 = simParams->getBondLambda(alchLambda2);
 const BigReal bond_lambda_22 = simParams->getBondLambda(1-alchLambda2);
 Molecule *const mol = Node::Object()->molecule;
 //fepe

 for ( int ituple=0; ituple<ntuple; ++ituple ) {
  const ImproperElem &tup = tuples[ituple];
  enum { size = 4 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const ImproperValue * const(&value)(tup.value);

  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << " " <<
	       localIndex[3] << std::endl);

  // Vector r12, r23, r34;	// vector between atoms
  Vector A,B,C;		// cross products
  BigReal rA, rB, rC;	// length of vectors A, B, and C
  BigReal energy=0;	// energy from the angle
  BigReal phi;		// angle between the plans
  double cos_phi;	// cos(phi)
  double sin_phi;	// sin(phi)
  Vector dcosdA;	// Derivative d(cos(phi))/dA
  Vector dcosdB;	// Derivative d(cos(phi))/dB
  Vector dsindC;	// Derivative d(sin(phi))/dC
  Vector dsindB;	// Derivative d(sin(phi))/dB
  BigReal K,K1;		// energy constants
  BigReal diff;		// for periodicity
  Force f1(0,0,0),f2(0,0,0),f3(0,0,0);	// force components

  //DebugM(3, "::computeForce() -- starting with improper type " << improperType << std::endl);

  // get the improper information
  int multiplicity = value->multiplicity;

  //  Calculate the vectors between atoms
  const Position & pos0 = p[0]->x[localIndex[0]].position;
  const Position & pos1 = p[1]->x[localIndex[1]].position;
  const Position & pos2 = p[2]->x[localIndex[2]].position;
  const Position & pos3 = p[3]->x[localIndex[3]].position;
  const Vector r12 = lattice.delta(pos0,pos1);
  const Vector r23 = lattice.delta(pos1,pos2);
  const Vector r34 = lattice.delta(pos2,pos3);

  //  Calculate the cross products
  A = cross(r12,r23);
  B = cross(r23,r34);
  C = cross(r23,A);

  //  Calculate the distances
  rA = A.length();
  rB = B.length();
  rC = C.length();

  //  Calculate the sin and cos
  cos_phi = A*B/(rA*rB);
  sin_phi = C*B/(rC*rB);

  //  Normalize B
  rB = 1.0/rB;
  B *= rB;

  phi= -atan2(sin_phi,cos_phi);

  if (fabs(sin_phi) > 0.1)
  {
    //  Normalize A
    rA = 1.0/rA;
    A *= rA;
    dcosdA = rA*(cos_phi*A-B);
    dcosdB = rB*(cos_phi*B-A);
  }
  else
  {
    //  Normalize C
    rC = 1.0/rC;
    C *= rC;
    dsindC = rC*(sin_phi*C-B);
    dsindB = rB*(sin_phi*B-C);
  }

  //  Loop through the multiple parameter sets for this
  //  bond.  We will only loop more than once if this
  //  has multiple parameter sets from Charmm22
  for (int mult_num=0; mult_num<multiplicity; mult_num++)
  {
    /* get angle information */
    Real k = value->values[mult_num].k * scale;
    Real delta = value->values[mult_num].delta;
    int n = value->values[mult_num].n;

    //  Calculate the energy
    if (n)
    {
      //  Periodicity is greater than 0, so use cos form
      K = k*(1+cos(n*phi - delta));
      K1 = -n*k*sin(n*phi - delta);
    }
    else
    {
      //  Periodicity is 0, so just use the harmonic form
      diff = phi-delta;
      if (diff < -PI)           diff += TWOPI;
      else if (diff > PI)       diff -= TWOPI;

      K = k*diff*diff;
      K1 = 2.0*k*diff;
    }

    //  Add the energy from this improper to the total energy
    energy += K;

    //  Next, we want to calculate the forces.  In order
    //  to do that, we first need to figure out whether the
    //  sin or cos form will be more stable.  For this,
    //  just look at the value of phi
    if (fabs(sin_phi) > 0.1)
    {
      //  use the sin version to avoid 1/cos terms
      K1 = K1/sin_phi;

      f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

      f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

      f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
               + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
               + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
             + r34.x*dcosdB.y - r34.y*dcosdB.x);
    }
    else
    {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      K1 = -K1/cos_phi;

      f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
      f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
      f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

      f3 += cross(K1,dsindB,r23);

      f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
             +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
             +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
             +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
             +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
             +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
             +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
             +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
             +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
             +dsindB.y*r34.x - dsindB.x*r34.y);
    }
  } /* for multiplicity */

  //fepb - BKR scaling of alchemical bonded terms
  //       NB: TI derivative is the _unscaled_ energy.
  if ( simParams->alchOn ) {
    switch ( mol->get_fep_bonded_type(atomID, 4) ) {
    case 1:
      reduction[improperEnergyIndex_ti_1] += energy;
      reduction[improperEnergyIndex_f] += (bond_lambda_12 - bond_lambda_1) * 
                                           energy;
      energy *= bond_lambda_1;
      f1 *= bond_lambda_1;
      f2 *= bond_lambda_1;
      f3 *= bond_lambda_1;
      break;
    case 2:
      reduction[improperEnergyIndex_ti_2] += energy;
      reduction[improperEnergyIndex_f] += (bond_lambda_22 - bond_lambda_2) *
                                           energy;
      energy *= bond_lambda_2;
      f1 *= bond_lambda_2;
      f2 *= bond_lambda_2;
      f3 *= bond_lambda_2;
      break;
    }
  }
  //fepe

  /* store the forces */
  p[0]->f[localIndex[0]] += f1;
  p[1]->f[localIndex[1]] += f2 - f1;
  p[2]->f[localIndex[2]] += f3 - f2;
  p[3]->f[localIndex[3]] += -f3;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << std::endl);
  reduction[improperEnergyIndex] += energy;
  reduction[virialIndex_XX] += ( f1.x * r12.x + f2.x * r23.x + f3.x * r34.x );
  reduction[virialIndex_XY] += ( f1.x * r12.y + f2.x * r23.y + f3.x * r34.y );
  reduction[virialIndex_XZ] += ( f1.x * r12.z + f2.x * r23.z + f3.x * r34.z );
  reduction[virialIndex_YX] += ( f1.y * r12.x + f2.y * r23.x + f3.y * r34.x );
  reduction[virialIndex_YY] += ( f1.y * r12.y + f2.y * r23.y + f3.y * r34.y );
  reduction[virialIndex_YZ] += ( f1.y * r12.z + f2.y * r23.z + f3.y * r34.z );
  reduction[virialIndex_ZX] += ( f1.z * r12.x + f2.z * r23.x + f3.z * r34.x );
  reduction[virialIndex_ZY] += ( f1.z * r12.y + f2.z * r23.y + f3.z * r34.y );
  reduction[virialIndex_ZZ] += ( f1.z * r12.z + f2.z * r23.z + f3.z * r34.z );

  if (pressureProfileData) {
    BigReal z1 = p[0]->x[localIndex[0]].position.z;
    BigReal z2 = p[1]->x[localIndex[1]].position.z;
    BigReal z3 = p[2]->x[localIndex[2]].position.z;
    BigReal z4 = p[3]->x[localIndex[3]].position.z;
    int n1 = (int)floor((z1-pressureProfileMin)/pressureProfileThickness);
    int n2 = (int)floor((z2-pressureProfileMin)/pressureProfileThickness);
    int n3 = (int)floor((z3-pressureProfileMin)/pressureProfileThickness);
    int n4 = (int)floor((z4-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(n1, pressureProfileSlabs);
    pp_clamp(n2, pressureProfileSlabs);
    pp_clamp(n3, pressureProfileSlabs);
    pp_clamp(n4, pressureProfileSlabs);
    int p1 = p[0]->x[localIndex[0]].partition;
    int p2 = p[1]->x[localIndex[1]].partition;
    int p3 = p[2]->x[localIndex[2]].partition;
    int p4 = p[3]->x[localIndex[3]].partition;
    int pn = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs, n1, n2,
                p1, p2, pn,
                f1.x * r12.x, f1.y * r12.y, f1.z * r12.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n2, n3,
                p2, p3, pn,
                f2.x * r23.x, f2.y * r23.y, f2.z * r23.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n3, n4,
                p3, p4, pn,
                f3.x * r34.x, f3.y * r34.y, f3.z * r34.z,
                pressureProfileData);
  }

 }
}

void ImproperElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_IMPROPER_ENERGY) += data[improperEnergyIndex];
  reduction->item(REDUCTION_BONDED_ENERGY_F) += data[improperEnergyIndex_f];
  reduction->item(REDUCTION_BONDED_ENERGY_TI_1) += data[improperEnergyIndex_ti_1];
  reduction->item(REDUCTION_BONDED_ENERGY_TI_2) += data[improperEnergyIndex_ti_2];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

