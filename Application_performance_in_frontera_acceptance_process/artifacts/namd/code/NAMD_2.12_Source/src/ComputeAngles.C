/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Methods for ComputeAngles.  Main code is for
   loading in the AngleElem information and
   for computing forces and energies for all angles on node's.
   HomePatch(es)
*/

#include "InfoStream.h"
#include "ComputeAngles.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"


// static initialization
int AngleElem::pressureProfileSlabs = 0;
int AngleElem::pressureProfileAtomTypes = 1;
BigReal AngleElem::pressureProfileThickness = 0;
BigReal AngleElem::pressureProfileMin = 0;

void AngleElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Angle** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in AngleElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numAngles;
  *byatom = mol->anglesByAtom;
  *structarray = mol->angles;
#endif
}

void AngleElem::getParameterPointers(Parameters *p, const AngleValue **v) {
  *v = p->angle_array;
}

void AngleElem::computeForce(AngleElem *tuples, int ntuple, BigReal *reduction, BigReal *pressureProfileData)
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
  const AngleElem &tup = tuples[ituple];
  enum { size = 3 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const AngleValue * const(&value)(tup.value);

  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << std::endl);

  const Position & pos1 = p[0]->x[localIndex[0]].position;
  const Position & pos2 = p[1]->x[localIndex[1]].position;
  Vector r12 = lattice.delta(pos1,pos2);
  register BigReal d12inv = r12.rlength();
  const Position & pos3 = p[2]->x[localIndex[2]].position;
  Vector r32 = lattice.delta(pos3,pos2);
  register BigReal d32inv = r32.rlength();
  int normal = value->normal;

  BigReal cos_theta = (r12*r32)*(d12inv*d32inv);
  //  Make sure that the cosine value is acceptable.  With roundoff, you
  //  can get values like 1.0+2e-16, which makes acos puke.  So instead,
  //  just set these kinds of values to exactly 1.0
  if (cos_theta > 1.0) cos_theta = 1.0;
  else if (cos_theta < -1.0) cos_theta = -1.0;

  BigReal k = value->k * scale;
  BigReal theta0 = value->theta0;

  //  Get theta
  BigReal theta = acos(cos_theta);

  //  Compare it to the rest angle
  BigReal diff;

  if (normal == 1) {
    diff = theta - theta0;
  } else {
    diff = cos_theta - cos(theta0);
  }
  
  //  Add the energy from this angle to the total energy
  BigReal energy = k *diff*diff;

  //  Normalize vector r12 and r32
  //BigReal d12inv = 1. / d12;
  //BigReal d32inv = 1. / d32;


  //  Calculate constant factor 2k(theta-theta0)/sin(theta)
  BigReal sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  if (normal != 1) {
    diff *= (2.0* k);
  } else if ( sin_theta < 1.e-6 ) {
    // catch case where bonds are parallel
    // this gets the force approximately right for theta0 of 0 or pi
    // and at least avoids small division for other values
    if ( diff < 0. ) diff = 2.0 * k;
    else diff = -2.0 * k;
  } else { 
    diff *= (-2.0* k) / sin_theta; 
  }
  BigReal c1 = diff * d12inv;
  BigReal c2 = diff * d32inv;

  //  Calculate the actual forces
  Force force1 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
  Force force2 = force1;
  Force force3 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
  force2 += force3;  force2 *= -1;

  //  Check to see if we need to do the Urey-Bradley term
  if (value->k_ub)
  {
	//  Non-zero k_ub value, so calculate the harmonic
	//  potential between the 1-3 atoms

  if (normal != 1) {
    NAMD_die("ERROR: Can't use cosAngles with Urey-Bradley angles");
  }
	BigReal k_ub = value->k_ub;
	BigReal r_ub = value->r_ub;
	Vector r13 = r12 - r32;
	BigReal d13 = r13.length();
	diff = d13- r_ub;

	energy += k_ub *diff*diff;

	diff *= -2.0*k_ub / d13;
	r13 *= diff;

	force1 += r13;
	force3 -= r13;
  }

  //fepb - BKR scaling of alchemical bonded terms
  //       NB: TI derivative is the _unscaled_ energy.
  if ( simParams->alchOn ) {
    switch ( mol->get_fep_bonded_type(atomID, 3) ) {
    case 1:
      reduction[angleEnergyIndex_ti_1] += energy;
      reduction[angleEnergyIndex_f] += (bond_lambda_12 - bond_lambda_1)*energy; 
      energy *= bond_lambda_1;
      force1 *= bond_lambda_1;
      force2 *= bond_lambda_1;
      force3 *= bond_lambda_1;
      break;
    case 2:
      reduction[angleEnergyIndex_ti_2] += energy;
      reduction[angleEnergyIndex_f] += (bond_lambda_22 - bond_lambda_2)*energy;
      energy *= bond_lambda_2;
      force1 *= bond_lambda_2;
      force2 *= bond_lambda_2;
      force3 *= bond_lambda_2;
      break;
    }
  }
  //fepe

  p[0]->f[localIndex[0]].x += force1.x;
  p[0]->f[localIndex[0]].y += force1.y;
  p[0]->f[localIndex[0]].z += force1.z;

  p[1]->f[localIndex[1]].x += force2.x;
  p[1]->f[localIndex[1]].y += force2.y;
  p[1]->f[localIndex[1]].z += force2.z;

  p[2]->f[localIndex[2]].x += force3.x;
  p[2]->f[localIndex[2]].y += force3.y;
  p[2]->f[localIndex[2]].z += force3.z;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << std::endl);

  reduction[angleEnergyIndex] += energy;
  reduction[virialIndex_XX] += ( force1.x * r12.x + force3.x * r32.x );
  reduction[virialIndex_XY] += ( force1.x * r12.y + force3.x * r32.y );
  reduction[virialIndex_XZ] += ( force1.x * r12.z + force3.x * r32.z );
  reduction[virialIndex_YX] += ( force1.y * r12.x + force3.y * r32.x );
  reduction[virialIndex_YY] += ( force1.y * r12.y + force3.y * r32.y );
  reduction[virialIndex_YZ] += ( force1.y * r12.z + force3.y * r32.z );
  reduction[virialIndex_ZX] += ( force1.z * r12.x + force3.z * r32.x );
  reduction[virialIndex_ZY] += ( force1.z * r12.y + force3.z * r32.y );
  reduction[virialIndex_ZZ] += ( force1.z * r12.z + force3.z * r32.z );

  if (pressureProfileData) {
    BigReal z1 = p[0]->x[localIndex[0]].position.z;
    BigReal z2 = p[1]->x[localIndex[1]].position.z;
    BigReal z3 = p[2]->x[localIndex[2]].position.z;
    int n1 = (int)floor((z1-pressureProfileMin)/pressureProfileThickness);
    int n2 = (int)floor((z2-pressureProfileMin)/pressureProfileThickness);
    int n3 = (int)floor((z3-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(n1, pressureProfileSlabs);
    pp_clamp(n2, pressureProfileSlabs);
    pp_clamp(n3, pressureProfileSlabs);
    int p1 = p[0]->x[localIndex[0]].partition;
    int p2 = p[1]->x[localIndex[1]].partition;
    int p3 = p[2]->x[localIndex[2]].partition;
    int pn = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs, n1, n2, 
                p1, p2, pn,
                force1.x * r12.x, force1.y * r12.y, force1.z * r12.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n3, n2, 
                p3, p2, pn,
                force3.x * r32.x, force3.y * r32.y, force3.z * r32.z,
                pressureProfileData);
  }

 }
}


void AngleElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ANGLE_ENERGY) += data[angleEnergyIndex];
  reduction->item(REDUCTION_BONDED_ENERGY_F) += data[angleEnergyIndex_f];
  reduction->item(REDUCTION_BONDED_ENERGY_TI_1) += data[angleEnergyIndex_ti_1];
  reduction->item(REDUCTION_BONDED_ENERGY_TI_2) += data[angleEnergyIndex_ti_2];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

