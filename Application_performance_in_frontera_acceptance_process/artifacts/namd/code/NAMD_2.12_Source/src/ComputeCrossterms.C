/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeCrossterms.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"

#define FASTER

enum {
  CMAP_TABLE_DIM = 25,
  CMAP_SETUP_DIM = 24,
  CMAP_SPACING = 15,
  CMAP_PHI_0 = -180,
  CMAP_PSI_0 = -180
};

#define INDEX(ncols,i,j)  ((i)*ncols + (j))


// static initialization
int CrosstermElem::pressureProfileSlabs = 0;
int CrosstermElem::pressureProfileAtomTypes = 1;
BigReal CrosstermElem::pressureProfileThickness = 0;
BigReal CrosstermElem::pressureProfileMin = 0;


void CrosstermElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Crossterm** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in CrosstermElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numCrossterms;
  *byatom = mol->crosstermsByAtom;
  *structarray = mol->crossterms;
#endif
}

void CrosstermElem::getParameterPointers(Parameters *p, const CrosstermValue **v) {
  *v = p->crossterm_array;
}

void CrosstermElem::computeForce(CrosstermElem *tuples, int ntuple, BigReal *reduction,
                                BigReal *pressureProfileData)
{
 const Lattice & lattice = tuples[0].p[0]->p->lattice;

 for ( int ituple=0; ituple<ntuple; ++ituple ) {
  const CrosstermElem &tup = tuples[ituple];
  enum { size = 8 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const CrosstermValue * const(&value)(tup.value);

  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << " " <<
	       localIndex[3] << std::endl);

  // Vector r12, r23, r34;	// vector between atoms
  Vector A,B,C,D,E,F;		// cross products
  BigReal rA,rB,rC,rD,rE,rF;	// length of vectors
  BigReal energy=0;	// energy from the angle
  BigReal phi,psi;		// angle between the plans
  double cos_phi,cos_psi;	// cos(phi)
  double sin_phi,sin_psi;	// sin(phi)
  Vector dcosdA,dcosdD;	// Derivative d(cos(phi))/dA
  Vector dcosdB,dcosdE;	// Derivative d(cos(phi))/dB
  Vector dsindC,dsindF;	// Derivative d(sin(phi))/dC
  Vector dsindB,dsindE;	// Derivative d(sin(phi))/dB
  BigReal U,U_phi,U_psi;	// energy constants
  BigReal diff;		// for periodicity
  Force f1(0,0,0),f2(0,0,0),f3(0,0,0);	// force components
  Force f4(0,0,0),f5(0,0,0),f6(0,0,0);	// force components

  //DebugM(3, "::computeForce() -- starting with crossterm type " << crosstermType << std::endl);

  //  Calculate the vectors between atoms
  const Position & pos0 = p[0]->x[localIndex[0]].position;
  const Position & pos1 = p[1]->x[localIndex[1]].position;
  const Position & pos2 = p[2]->x[localIndex[2]].position;
  const Position & pos3 = p[3]->x[localIndex[3]].position;
  const Position & pos4 = p[4]->x[localIndex[4]].position;
  const Position & pos5 = p[5]->x[localIndex[5]].position;
  const Position & pos6 = p[6]->x[localIndex[6]].position;
  const Position & pos7 = p[7]->x[localIndex[7]].position;
  const Vector r12 = lattice.delta(pos0,pos1);
  const Vector r23 = lattice.delta(pos1,pos2);
  const Vector r34 = lattice.delta(pos2,pos3);
  const Vector r56 = lattice.delta(pos4,pos5);
  const Vector r67 = lattice.delta(pos5,pos6);
  const Vector r78 = lattice.delta(pos6,pos7);

  //  Calculate the cross products
  A = cross(r12,r23);
  B = cross(r23,r34);
  C = cross(r23,A);
  D = cross(r56,r67);
  E = cross(r67,r78);
  F = cross(r67,D);

  //  Calculate the distances
  rA = A.length();
  rB = B.length();
  rC = C.length();
  rD = D.length();
  rE = E.length();
  rF = F.length();

  //  Calculate the sin and cos
  cos_phi = A*B/(rA*rB);
  sin_phi = C*B/(rC*rB);
  cos_psi = D*E/(rD*rE);
  sin_psi = F*E/(rF*rE);

  //  Normalize B
  rB = 1.0/rB;
  B *= rB;
  rE = 1.0/rE;
  E *= rE;

  phi= -atan2(sin_phi,cos_phi);
  psi= -atan2(sin_psi,cos_psi);

  if (fabs(sin_phi) > 0.1) {
    //  Normalize A
    rA = 1.0/rA;
    A *= rA;
    dcosdA = rA*(cos_phi*A-B);
    dcosdB = rB*(cos_phi*B-A);
  } else {
    //  Normalize C
    rC = 1.0/rC;
    C *= rC;
    dsindC = rC*(sin_phi*C-B);
    dsindB = rB*(sin_phi*B-C);
  }

  if (fabs(sin_psi) > 0.1) {
    //  Normalize A
    rD = 1.0/rD;
    D *= rD;
    dcosdD = rD*(cos_psi*D-E);
    dcosdE = rE*(cos_psi*E-D);
  } else {
    //  Normalize C
    rF = 1.0/rF;
    F *= rF;
    dsindF = rF*(sin_psi*F-E);
    dsindE = rE*(sin_psi*E-F);
  }

    //  Calculate the energy
{
  const double h = CMAP_SPACING * PI / 180.0;
  const double h_1 = 1.0 / h;
#ifdef FASTER
  const double six_h = 6.0 * h_1;
#endif
  const double phi_0 = CMAP_PHI_0 * PI / 180.0;
  const double psi_0 = CMAP_PSI_0 * PI / 180.0;

  enum { D = CMAP_TABLE_DIM };
  double xa[2], xb[2], dxa[2], dxb[2];
  double ya[2], yb[2], dya[2], dyb[2];
  double t, dx_h, dy_h;
#ifdef FASTER
  double s1, s2, s3, s4, s5;
#endif
  double f = 0, fx = 0, fy = 0;
  const CrosstermData *table = &value->c[0][0];
  int i, j, ilo, jlo, ij;

  /* distance measured in grid points between angle and smallest value */
  dx_h = (phi - phi_0) * h_1;
  dy_h = (psi - psi_0) * h_1;

  /* find smallest numbered grid point in stencil */
  ilo = (int) floor(dx_h);
  if ( ilo < 0 ) ilo = 0; 
  if ( ilo >= CMAP_SETUP_DIM ) ilo = CMAP_SETUP_DIM - 1; 
  jlo = (int) floor(dy_h);
  if ( jlo < 0 ) jlo = 0; 
  if ( jlo >= CMAP_SETUP_DIM ) jlo = CMAP_SETUP_DIM - 1; 

#if !defined(FASTER)

  /* find t for x-dimension and compute xa, xb, dxa, dxb */
  t = dx_h - (double) ilo;
  xa[0] = (1 - t) * (1 - t) * (1 + 2*t);
  xb[0] = h * t * (1 - t) * (1 - t);
  dxa[0] = -6 * t * (1 - t) * h_1;
  dxb[0] = (1 - t) * (1 - 3*t);
  t--;
  xa[1] = (1 + t) * (1 + t) * (1 - 2*t);
  xb[1] = h * t * (1 + t) * (1 + t);
  dxa[1] = -6 * t * (1 + t) * h_1;
  dxb[1] = (1 + t) * (1 + 3*t);

  /* find t for y-dimension and compute ya, yb, dya, dyb */
  t = dy_h - (double) jlo;
  ya[0] = (1 - t) * (1 - t) * (1 + 2*t);
  yb[0] = h * t * (1 - t) * (1 - t);
  dya[0] = -6 * t * (1 - t) * h_1;
  dyb[0] = (1 - t) * (1 - 3*t);
  t--;
  ya[1] = (1 + t) * (1 + t) * (1 - 2*t);
  yb[1] = h * t * (1 + t) * (1 + t);
  dya[1] = -6 * t * (1 + t) * h_1;
  dyb[1] = (1 + t) * (1 + 3*t);

#else

  /* find t for x-dimension and compute xa, xb, dxa, dxb */
  t = dx_h - (double) ilo;
  s1 = 1-t;
  s2 = 2*t;
  s3 = 3*t;
  s4 = t*s1;
  s5 = h*s4;
  xa[0] = s1*s1*(1+s2);
  xa[1] = t*t*(3-s2);
  xb[0] = s5*s1;
  xb[1] = -s5*t;
  dxa[0] = -six_h*s4;
  dxa[1] = -dxa[0];
  dxb[0] = s1*(1-s3);
  dxb[1] = t*(-2+s3);

  /* find t for y-dimension and compute ya, yb, dya, dyb */
  t = dy_h - (double) jlo;
  s1 = 1-t;
  s2 = 2*t;
  s3 = 3*t;
  s4 = t*s1;
  s5 = h*s4;
  ya[0] = s1*s1*(1+s2);
  ya[1] = t*t*(3-s2);
  yb[0] = s5*s1;
  yb[1] = -s5*t;
  dya[0] = -six_h*s4;
  dya[1] = -dya[0];
  dyb[0] = s1*(1-s3);
  dyb[1] = t*(-2+s3);

#endif

  for (i = 0;  i < 2;  i++) {
    for (j = 0;  j < 2;  j++) {
      ij = INDEX(D,i+ilo,j+jlo);

#if !defined(FASTER)

      f += xa[i] * ya[j] * table[ij].d00
        + xb[i] * ya[j] * table[ij].d10
        + xa[i] * yb[j] * table[ij].d01
        + xb[i] * yb[j] * table[ij].d11;

      fx += dxa[i] * ya[j] * table[ij].d00
        + dxb[i] * ya[j] * table[ij].d10
        + dxa[i] * yb[j] * table[ij].d01
        + dxb[i] * yb[j] * table[ij].d11;

      fy += xa[i] * dya[j] * table[ij].d00
        + xb[i] * dya[j] * table[ij].d10
        + xa[i] * dyb[j] * table[ij].d01
        + xb[i] * dyb[j] * table[ij].d11;

#else

      s1=ya[j]*table[ij].d00;
      s2=yb[j]*table[ij].d01;
      s3=ya[j]*table[ij].d10;
      s4=yb[j]*table[ij].d11;

      f+=xa[i]*(s1+s2)+xb[i]*(s3+s4);
      fx+=dxa[i]*(s1+s2)+dxb[i]*(s3+s4);
      fy+=xa[i]*(dya[j]*table[ij].d00+dyb[j]*table[ij].d01)
        +xb[i]*(dya[j]*table[ij].d10+dyb[j]*table[ij].d11);

#endif
    }
  }

  /* return accumulated values */
  U = f * scale;
  U_phi = fx * scale;
  U_psi = fy * scale;

/*
CkPrintf("crossterm %d-%d-%d-%d %d-%d-%d-%d %lf %lf %d %d %lf %lf %lf\n",
    atomID[0], atomID[1], atomID[2], atomID[3],
    atomID[4], atomID[5], atomID[6], atomID[7],
    phi, psi, ilo, jlo, U, U_phi, U_psi);
CkPrintf("%d %d-%d-%d-%d %d-%d-%d-%d\n", CkMyPe(),
   p[0]->patchID, p[1]->patchID, p[2]->patchID, p[3]->patchID,
   p[4]->patchID, p[5]->patchID, p[6]->patchID, p[7]->patchID);
*/

}

    //  Add the energy from this crossterm to the total energy
    energy += U;

    //  Next, we want to calculate the forces.  In order
    //  to do that, we first need to figure out whether the
    //  sin or cos form will be more stable.  For this,
    //  just look at the value of phi
    if (fabs(sin_phi) > 0.1)
    {
      //  use the sin version to avoid 1/cos terms
      U_phi = U_phi/sin_phi;

      f1.x += U_phi*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y += U_phi*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z += U_phi*(r23.x*dcosdA.y - r23.y*dcosdA.x);

      f3.x += U_phi*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y += U_phi*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z += U_phi*(r23.y*dcosdB.x - r23.x*dcosdB.y);

      f2.x += U_phi*(r12.z*dcosdA.y - r12.y*dcosdA.z
               + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y += U_phi*(r12.x*dcosdA.z - r12.z*dcosdA.x
               + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z += U_phi*(r12.y*dcosdA.x - r12.x*dcosdA.y
             + r34.x*dcosdB.y - r34.y*dcosdB.x);
    }
    else
    {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      U_phi = -U_phi/cos_phi;

      f1.x += U_phi*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
      f1.y += U_phi*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
      f1.z += U_phi*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

      f3 += cross(U_phi,dsindB,r23);

      f2.x += U_phi*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
             +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
             +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
             +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y += U_phi*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
             +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
             +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
             +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z += U_phi*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
             +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
             +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
             +dsindB.y*r34.x - dsindB.x*r34.y);
    }

    if (fabs(sin_psi) > 0.1)
    {
      //  use the sin version to avoid 1/cos terms
      U_psi = U_psi/sin_psi;

      f4.x += U_psi*(r67.y*dcosdD.z - r67.z*dcosdD.y);
      f4.y += U_psi*(r67.z*dcosdD.x - r67.x*dcosdD.z);
      f4.z += U_psi*(r67.x*dcosdD.y - r67.y*dcosdD.x);

      f6.x += U_psi*(r67.z*dcosdE.y - r67.y*dcosdE.z);
      f6.y += U_psi*(r67.x*dcosdE.z - r67.z*dcosdE.x);
      f6.z += U_psi*(r67.y*dcosdE.x - r67.x*dcosdE.y);

      f5.x += U_psi*(r56.z*dcosdD.y - r56.y*dcosdD.z
               + r78.y*dcosdE.z - r78.z*dcosdE.y);
      f5.y += U_psi*(r56.x*dcosdD.z - r56.z*dcosdD.x
               + r78.z*dcosdE.x - r78.x*dcosdE.z);
      f5.z += U_psi*(r56.y*dcosdD.x - r56.x*dcosdD.y
             + r78.x*dcosdE.y - r78.y*dcosdE.x);
    }
    else
    {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      U_psi = -U_psi/cos_psi;

      f4.x += U_psi*((r67.y*r67.y + r67.z*r67.z)*dsindF.x
                - r67.x*r67.y*dsindF.y
                - r67.x*r67.z*dsindF.z);
      f4.y += U_psi*((r67.z*r67.z + r67.x*r67.x)*dsindF.y
                - r67.y*r67.z*dsindF.z
                - r67.y*r67.x*dsindF.x);
      f4.z += U_psi*((r67.x*r67.x + r67.y*r67.y)*dsindF.z
                - r67.z*r67.x*dsindF.x
                - r67.z*r67.y*dsindF.y);

      f6 += cross(U_psi,dsindE,r67);

      f5.x += U_psi*(-(r67.y*r56.y + r67.z*r56.z)*dsindF.x
             +(2.0*r67.x*r56.y - r56.x*r67.y)*dsindF.y
             +(2.0*r67.x*r56.z - r56.x*r67.z)*dsindF.z
             +dsindE.z*r78.y - dsindE.y*r78.z);
      f5.y += U_psi*(-(r67.z*r56.z + r67.x*r56.x)*dsindF.y
             +(2.0*r67.y*r56.z - r56.y*r67.z)*dsindF.z
             +(2.0*r67.y*r56.x - r56.y*r67.x)*dsindF.x
             +dsindE.x*r78.z - dsindE.z*r78.x);
      f5.z += U_psi*(-(r67.x*r56.x + r67.y*r56.y)*dsindF.z
             +(2.0*r67.z*r56.x - r56.z*r67.x)*dsindF.x
             +(2.0*r67.z*r56.y - r56.z*r67.y)*dsindF.y
             +dsindE.y*r78.x - dsindE.x*r78.y);
    }

  /* store the forces */
  p[0]->f[localIndex[0]] += f1;
  p[1]->f[localIndex[1]] += f2 - f1;
  p[2]->f[localIndex[2]] += f3 - f2;
  p[3]->f[localIndex[3]] += -f3;
  p[4]->f[localIndex[4]] += f4;
  p[5]->f[localIndex[5]] += f5 - f4;
  p[6]->f[localIndex[6]] += f6 - f5;
  p[7]->f[localIndex[7]] += -f6;

  if ( p[0]->af ) {
    p[0]->af[localIndex[0]] += f1;
    p[1]->af[localIndex[1]] += f2 - f1;
    p[2]->af[localIndex[2]] += f3 - f2;
    p[3]->af[localIndex[3]] += -f3;
    p[4]->af[localIndex[4]] += f4;
    p[5]->af[localIndex[5]] += f5 - f4;
    p[6]->af[localIndex[6]] += f6 - f5;
    p[7]->af[localIndex[7]] += -f6;
  }

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << std::endl);
  reduction[crosstermEnergyIndex] += energy;
  reduction[virialIndex_XX] += ( f1.x * r12.x + f2.x * r23.x + f3.x * r34.x );
  reduction[virialIndex_XY] += ( f1.x * r12.y + f2.x * r23.y + f3.x * r34.y );
  reduction[virialIndex_XZ] += ( f1.x * r12.z + f2.x * r23.z + f3.x * r34.z );
  reduction[virialIndex_YX] += ( f1.y * r12.x + f2.y * r23.x + f3.y * r34.x );
  reduction[virialIndex_YY] += ( f1.y * r12.y + f2.y * r23.y + f3.y * r34.y );
  reduction[virialIndex_YZ] += ( f1.y * r12.z + f2.y * r23.z + f3.y * r34.z );
  reduction[virialIndex_ZX] += ( f1.z * r12.x + f2.z * r23.x + f3.z * r34.x );
  reduction[virialIndex_ZY] += ( f1.z * r12.y + f2.z * r23.y + f3.z * r34.y );
  reduction[virialIndex_ZZ] += ( f1.z * r12.z + f2.z * r23.z + f3.z * r34.z );

  reduction[virialIndex_XX] += ( f4.x * r56.x + f5.x * r67.x + f6.x * r78.x );
  reduction[virialIndex_XY] += ( f4.x * r56.y + f5.x * r67.y + f6.x * r78.y );
  reduction[virialIndex_XZ] += ( f4.x * r56.z + f5.x * r67.z + f6.x * r78.z );
  reduction[virialIndex_YX] += ( f4.y * r56.x + f5.y * r67.x + f6.y * r78.x );
  reduction[virialIndex_YY] += ( f4.y * r56.y + f5.y * r67.y + f6.y * r78.y );
  reduction[virialIndex_YZ] += ( f4.y * r56.z + f5.y * r67.z + f6.y * r78.z );
  reduction[virialIndex_ZX] += ( f4.z * r56.x + f5.z * r67.x + f6.z * r78.x );
  reduction[virialIndex_ZY] += ( f4.z * r56.y + f5.z * r67.y + f6.z * r78.y );
  reduction[virialIndex_ZZ] += ( f4.z * r56.z + f5.z * r67.z + f6.z * r78.z );

  if (pressureProfileData) {
    BigReal z1 = p[0]->x[localIndex[0]].position.z;
    BigReal z2 = p[1]->x[localIndex[1]].position.z;
    BigReal z3 = p[2]->x[localIndex[2]].position.z;
    BigReal z4 = p[3]->x[localIndex[3]].position.z;
    BigReal z5 = p[4]->x[localIndex[4]].position.z;
    BigReal z6 = p[5]->x[localIndex[5]].position.z;
    BigReal z7 = p[6]->x[localIndex[6]].position.z;
    BigReal z8 = p[7]->x[localIndex[7]].position.z;
    int n1 = (int)floor((z1-pressureProfileMin)/pressureProfileThickness);
    int n2 = (int)floor((z2-pressureProfileMin)/pressureProfileThickness);
    int n3 = (int)floor((z3-pressureProfileMin)/pressureProfileThickness);
    int n4 = (int)floor((z4-pressureProfileMin)/pressureProfileThickness);
    int n5 = (int)floor((z5-pressureProfileMin)/pressureProfileThickness);
    int n6 = (int)floor((z6-pressureProfileMin)/pressureProfileThickness);
    int n7 = (int)floor((z7-pressureProfileMin)/pressureProfileThickness);
    int n8 = (int)floor((z8-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(n1, pressureProfileSlabs);
    pp_clamp(n2, pressureProfileSlabs);
    pp_clamp(n3, pressureProfileSlabs);
    pp_clamp(n4, pressureProfileSlabs);
    pp_clamp(n5, pressureProfileSlabs);
    pp_clamp(n6, pressureProfileSlabs);
    pp_clamp(n7, pressureProfileSlabs);
    pp_clamp(n8, pressureProfileSlabs);
    int p1 = p[0]->x[localIndex[0]].partition;
    int p2 = p[1]->x[localIndex[1]].partition;
    int p3 = p[2]->x[localIndex[2]].partition;
    int p4 = p[3]->x[localIndex[3]].partition;
    int p5 = p[4]->x[localIndex[4]].partition;
    int p6 = p[5]->x[localIndex[5]].partition;
    int p7 = p[6]->x[localIndex[6]].partition;
    int p8 = p[7]->x[localIndex[7]].partition;
    int pn = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs, n1, n2, p1, p2, pn,
                f1.x * r12.x, f1.y * r12.y, f1.z * r12.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n2, n3, p2, p3, pn,
                f2.x * r23.x, f2.y * r23.y, f2.z * r23.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n3, n4, p3, p4, pn,
                f3.x * r34.x, f3.y * r34.y, f3.z * r34.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n5, n6, p5, p6, pn,
                f4.x * r56.x, f4.y * r56.y, f4.z * r56.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n6, n7, p6, p7, pn,
                f5.x * r67.x, f5.y * r67.y, f5.z * r67.z,
                pressureProfileData);
    pp_reduction(pressureProfileSlabs, n7, n8, p7, p8, pn,
                f6.x * r78.x, f6.y * r78.y, f6.z * r78.z,
                pressureProfileData);
  }

 }
}

void CrosstermElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_CROSSTERM_ENERGY) += data[crosstermEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_AMD_DIHE,data,virialIndex);
}


/******************************************************************************/

static void lu_decomp_nopivot(double *m, int n);
static void forward_back_sub(double *b, double *m, int n);

void crossterm_setup(CrosstermData *table)
{
  enum { D = CMAP_TABLE_DIM };
  enum { N = CMAP_SETUP_DIM };
  const double h_1 = 1.0 / ( CMAP_SPACING * PI / 180.0) ;
  const double h_2 = h_1 * h_1;
  const double tr_h = 3.0 * h_1;
  int i, j;
  int ij;
  int ijp1p1, ijm1m1, ijp1m1, ijm1p1;
  int ijp2p2, ijm2m2, ijp2m2, ijm2p2;

  /* allocate spline coefficient matrix */
  double* const m = new double[N*N];
  memset(m,0,N*N*sizeof(double));

  /* initialize spline coefficient matrix */
  m[0] = 4;
  for (i = 1;  i < N;  i++) {
    m[INDEX(N,i-1,i)] = 1;
    m[INDEX(N,i,i-1)] = 1;
    m[INDEX(N,i,i)] = 4;
  }
  /* periodic boundary conditions for spline */
  m[INDEX(N,0,N-1)] = 1;
  m[INDEX(N,N-1,0)] = 1;

  /* compute LU-decomposition for this matrix */
  lu_decomp_nopivot(m, N);

  /* allocate vector for solving spline derivatives */
  double* const v = new double[N];
  memset(v,0,N*sizeof(double));

    /* march through rows of table */
    for (i = 0;  i < N;  i++) {

      /* setup RHS vector for solving spline derivatives */
      v[0] = tr_h * (table[INDEX(D,i,1)].d00 - table[INDEX(D,i,N-1)].d00);
      for (j = 1;  j < N;  j++) {
        v[j] = tr_h * (table[INDEX(D,i,j+1)].d00 - table[INDEX(D,i,j-1)].d00);
      }

      /* solve system, returned into vector */
      forward_back_sub(v, m, N);

      /* store values as derivatives wrt differenced table values */
      for (j = 0;  j < N;  j++) {
        table[INDEX(D,i,j)].d01 = v[j];
      }
      table[INDEX(D,i,N)].d01 = v[0];
    }
    for (j = 0;  j <= N;  j++) {
      table[INDEX(D,N,j)].d01 = table[INDEX(D,0,j)].d01;
    }

    /* march through columns of table */
    for (j = 0;  j < N;  j++) {

      /* setup RHS vector for solving spline derivatives */
      v[0] = tr_h * (table[INDEX(D,1,j)].d00 - table[INDEX(D,N-1,j)].d00);
      for (i = 1;  i < N;  i++) {
        v[i] = tr_h * (table[INDEX(D,i+1,j)].d00 - table[INDEX(D,i-1,j)].d00);
      }

      /* solve system, returned into vector */
      forward_back_sub(v, m, N);

      /* store values as derivatives wrt differenced table values */
      for (i = 0;  i < N;  i++) {
        table[INDEX(D,i,j)].d10 = v[i];
      }
      table[INDEX(D,N,j)].d10 = v[0];
    }
    for (i = 0;  i <= N;  i++) {
      table[INDEX(D,i,N)].d10 = table[INDEX(D,i,0)].d10;
    }

    /* march back through rows of table
     *
     * This is CHARMM's approach for calculating mixed partial derivatives,
     * by splining the first derivative values.
     *
     * Here we spline the dx values along y to calculate dxdy derivatives.
     *
     * Test cases show error with CHARMM is within 1e-5 roundoff error.
     */
    for (i = 0;  i < N;  i++) {

      /* setup RHS vector for solving dxy derivatives from dx */
      v[0] = tr_h * (table[INDEX(D,i,1)].d10 - table[INDEX(D,i,N-1)].d10);
      for (j = 1;  j < N;  j++) {
        v[j] = tr_h * (table[INDEX(D,i,j+1)].d10 - table[INDEX(D,i,j-1)].d10);
      }

      /* solve system, returned into vector */
      forward_back_sub(v, m, N);

      /* store values as dxy derivatives wrt differenced table values */
      for (j = 0;  j < N;  j++) {
        table[INDEX(D,i,j)].d11 = v[j];
      }
      table[INDEX(D,i,N)].d11 = v[0];
    }
    for (j = 0;  j <= N;  j++) {
      table[INDEX(D,N,j)].d11 = table[INDEX(D,0,j)].d11;
    }

  /* done with temp storage */
  delete [] m;
  delete [] v;

}


void lu_decomp_nopivot(double *m, int n)
{
  double l_ik;
  int i, j, k;

  for (k = 0;  k < n-1;  k++) {
    for (i = k+1;  i < n;  i++) {
      l_ik = m[INDEX(n,i,k)] / m[INDEX(n,k,k)];
      for (j = k;  j < n;  j++) {
        m[INDEX(n,i,j)] -= l_ik * m[INDEX(n,k,j)];
      }
      m[INDEX(n,i,k)] = l_ik;
    }
  }
}


void forward_back_sub(double *b, double *m, int n)
{
  int i, j;

  /* in place forward elimination (solves Ly=b) using lower triangle of m */
  for (j = 0;  j < n-1;  j++) {
    for (i = j+1;  i < n;  i++) {
      b[i] -= m[INDEX(n,i,j)] * b[j];
    }
  }
  /* in place back substitution (solves Ux=y) using upper triangle of m */
  for (j = n-1;  j >= 0;  j--) {
    b[j] /= m[INDEX(n,j,j)];
    for (i = j-1;  i >= 0;  i--) {
      b[i] -= m[INDEX(n,i,j)] * b[j];
    }
  }
}


