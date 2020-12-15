/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeDihedrals.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"


// static initialization
int DihedralElem::pressureProfileSlabs = 0;
int DihedralElem::pressureProfileAtomTypes = 1;
BigReal DihedralElem::pressureProfileThickness = 0;
BigReal DihedralElem::pressureProfileMin = 0;

void DihedralElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Dihedral** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in DihedralElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numDihedrals;
  *byatom = mol->dihedralsByAtom;
  *structarray = mol->dihedrals;
#endif
}

void DihedralElem::getParameterPointers(Parameters *p, const DihedralValue **v) {
  *v = p->dihedral_array;
}

void DihedralElem::computeForce(DihedralElem *tuples, int ntuple, BigReal *reduction, 
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
  const DihedralElem &tup = tuples[ituple];
  enum { size = 4 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const DihedralValue * const(&value)(tup.value);

  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << std::endl);

  //  Calculate the vectors between atoms
  const Position & pos0 = p[0]->x[localIndex[0]].position;
  const Position & pos1 = p[1]->x[localIndex[1]].position;
  const Vector r12 = lattice.delta(pos0,pos1);
  const Position & pos2 = p[2]->x[localIndex[2]].position;
  const Vector r23 = lattice.delta(pos1,pos2);
  const Position & pos3 = p[3]->x[localIndex[3]].position;
  const Vector r34 = lattice.delta(pos2,pos3);

  //  Calculate the cross products and distances
  Vector A = cross(r12,r23);
  register  BigReal rAinv = A.rlength();
  Vector B = cross(r23,r34);
  register  BigReal rBinv = B.rlength();
  Vector C = cross(r23,A);
  register  BigReal rCinv = C.rlength();

  //  Calculate the sin and cos
  BigReal cos_phi = (A*B)*(rAinv*rBinv);
  BigReal sin_phi = (C*B)*(rCinv*rBinv);

  BigReal phi= -atan2(sin_phi,cos_phi);

  BigReal K=0;		// energy
  BigReal K1=0;		// force

  // get the dihedral information
  int multiplicity = value->multiplicity;

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
      K += k*(1+cos(n*phi - delta));
      K1 += -n*k*sin(n*phi - delta);
    }
    else
    {
      //  Periodicity is 0, so just use the harmonic form
      BigReal diff = phi-delta;
      if (diff < -PI)           diff += TWOPI;
      else if (diff > PI)       diff -= TWOPI;

      K += k*diff*diff;
      K1 += 2.0*k*diff;
    }
  } /* for multiplicity */

  //fepb - BKR scaling of alchemical bonded terms
  //       NB: TI derivative is the _unscaled_ energy.
  if ( simParams->alchOn ) {
    switch ( mol->get_fep_bonded_type(atomID, 4) ) {
    case 1:
      reduction[dihedralEnergyIndex_ti_1] += K;
      reduction[dihedralEnergyIndex_f] += (bond_lambda_12 - bond_lambda_1)*K;
      K *= bond_lambda_1;
      K1 *= bond_lambda_1;
      break;
    case 2:
      reduction[dihedralEnergyIndex_ti_2] += K;
      reduction[dihedralEnergyIndex_f] += (bond_lambda_22 - bond_lambda_2)*K;
      K *= bond_lambda_2;
      K1 *= bond_lambda_2;
      break;
    }
  }
  //fepe

  Force f1,f2,f3;

  //  Normalize B
  //rB = 1.0/rB;
  B *= rBinv;

    //  Next, we want to calculate the forces.  In order
    //  to do that, we first need to figure out whether the
    //  sin or cos form will be more stable.  For this,
    //  just look at the value of phi
    if (fabs(sin_phi) > 0.1)
    {
      //  use the sin version to avoid 1/cos terms

      //  Normalize A
      A *= rAinv;
      Vector dcosdA;
      Vector dcosdB;

      dcosdA.x = rAinv*(cos_phi*A.x-B.x);
      dcosdA.y = rAinv*(cos_phi*A.y-B.y);
      dcosdA.z = rAinv*(cos_phi*A.z-B.z);
	    
      dcosdB.x = rBinv*(cos_phi*B.x-A.x);
      dcosdB.y = rBinv*(cos_phi*B.y-A.y);
      dcosdB.z = rBinv*(cos_phi*B.z-A.z);

      K1 = K1/sin_phi;

      f1.x = K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y = K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z = K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);
			     		      
      f3.x = K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y = K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z = K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);
			     		      
      f2.x = K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
               + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y = K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
               + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z = K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
	       + r34.x*dcosdB.y - r34.y*dcosdB.x);
    }
    else
    {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms

      //  Normalize C
      //      rC = 1.0/rC;
      C *= rCinv;
      
      Vector dsindC;
      Vector dsindB;

      dsindC.x = rCinv*(sin_phi*C.x-B.x);
      dsindC.y = rCinv*(sin_phi*C.y-B.y);
      dsindC.z = rCinv*(sin_phi*C.z-B.z);

      dsindB.x = rBinv*(sin_phi*B.x-C.x);
      dsindB.y = rBinv*(sin_phi*B.y-C.y);
      dsindB.z = rBinv*(sin_phi*B.z-C.z);

      K1 = -K1/cos_phi;

      f1.x = K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
      f1.y = K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
      f1.z = K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

      f3 = cross(K1,dsindB,r23);

      f2.x = K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
             +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
             +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
             +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y = K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
             +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
             +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
             +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z = K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
             +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
             +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
             +dsindB.y*r34.x - dsindB.x*r34.y);
    }

  /* store the forces */
  //  p[0]->f[localIndex[0]] += f1;
  //  p[1]->f[localIndex[1]] += f2 - f1;
  //  p[2]->f[localIndex[2]] += f3 - f2;
  //  p[3]->f[localIndex[3]] += -f3;

  p[0]->f[localIndex[0]].x += f1.x;
  p[0]->f[localIndex[0]].y += f1.y;
  p[0]->f[localIndex[0]].z += f1.z;

  p[1]->f[localIndex[1]].x += f2.x - f1.x;
  p[1]->f[localIndex[1]].y += f2.y - f1.y;
  p[1]->f[localIndex[1]].z += f2.z - f1.z;

  p[2]->f[localIndex[2]].x += f3.x - f2.x;
  p[2]->f[localIndex[2]].y += f3.y - f2.y;
  p[2]->f[localIndex[2]].z += f3.z - f2.z;

  p[3]->f[localIndex[3]].x += -f3.x;
  p[3]->f[localIndex[3]].y += -f3.y;
  p[3]->f[localIndex[3]].z += -f3.z;  

    /* store the force for dihedral-only accelMD */
  if ( p[0]->af ) {
    p[0]->af[localIndex[0]].x += f1.x;
    p[0]->af[localIndex[0]].y += f1.y;
    p[0]->af[localIndex[0]].z += f1.z;

    p[1]->af[localIndex[1]].x += f2.x - f1.x;
    p[1]->af[localIndex[1]].y += f2.y - f1.y;
    p[1]->af[localIndex[1]].z += f2.z - f1.z;

    p[2]->af[localIndex[2]].x += f3.x - f2.x;
    p[2]->af[localIndex[2]].y += f3.y - f2.y;
    p[2]->af[localIndex[2]].z += f3.z - f2.z;

    p[3]->af[localIndex[3]].x += -f3.x;
    p[3]->af[localIndex[3]].y += -f3.y;
    p[3]->af[localIndex[3]].z += -f3.z;
  }

  DebugM(3, "::computeForce() -- ending with delta energy " << K << std::endl);
  reduction[dihedralEnergyIndex] += K;
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


void DihedralElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_DIHEDRAL_ENERGY) += data[dihedralEnergyIndex];
  reduction->item(REDUCTION_BONDED_ENERGY_F) += data[dihedralEnergyIndex_f];
  reduction->item(REDUCTION_BONDED_ENERGY_TI_1) += data[dihedralEnergyIndex_ti_1];
  reduction->item(REDUCTION_BONDED_ENERGY_TI_2) += data[dihedralEnergyIndex_ti_2];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_AMD_DIHE,data,virialIndex);
}

