/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <string.h>
//#include <iomanip.h>
#include "common.h"
#include "InfoStream.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyEnums.h"
#include "Vector.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyGroup.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"

#include "NamdTypes.h"
#include "GlobalMaster.h"
#include "GlobalMasterFreeEnergy.h"

// initialize static member variables
// (lambda is the same for all forcing restraints)
double  ARestraint::m_LambdaKf  = 1.0;
double  ARestraint::m_LambdaRef = 0.0;

#if defined(_DEBUG)
void ARestraint::SetCOM(int Index, AVector& Pos) {
//-----------------------------------------------------------------
// only for testing!!!!!!
//-----------------------------------------------------------------
  ASSERT( (Index>=0) && (Index<m_NumGroups) );
  m_pCOMs[Index] = Pos;
}
#endif

ARestraint::ARestraint() {
//-----------------------------------------------------------------
// constructor for base class
//-----------------------------------------------------------------
  m_pGroups = NULL;
  m_pCOMs = NULL;
  m_NumGroups = 0;
}


ARestraint::~ARestraint() {
//-----------------------------------------------------------------
// free space that may have been allocated for Groups and COM's
//-----------------------------------------------------------------
  if (m_pGroups != NULL) {
    ASSERT(m_pCOMs != NULL);
    delete []m_pGroups;
    delete []m_pCOMs;
  }
}


void ARestraint::EarlyExit(char* Str, int AtomID) {
//-----------------------------------------------------------------
// unrecoverable error
//-----------------------------------------------------------------
  char  NumStr[40];

  iout << "FreeEnergy: " << std::endl << endi;
  sprintf(NumStr, "%d", AtomID);
  strcat(Str, " for AtomID: ");
  strcat(Str, NumStr);
  iout << "FreeEnergy: " << Str;
  iout << std::endl << endi;
  NAMD_die("FreeEnergy: Fatal Error with Fixed or Forcing Restraints");
}


void ARestraint::SetGroup(AGroup& Group, int GroupIndex) {
//-----------------------------------------------------------------
// set one group of atoms
//-----------------------------------------------------------------
  ASSERT( (GroupIndex>=0) && (GroupIndex<m_NumGroups) );
  m_pGroups[GroupIndex] = Group;
}


void ARestraint::SetGroups(AGroup& Group1) {
//-----------------------------------------------------------------
// set one group of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 1);
  m_pGroups[0] = Group1;
}


void ARestraint::SetGroups(AGroup& Group1, AGroup& Group2) {
//-----------------------------------------------------------------
// set two groups of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 2);
  m_pGroups[0] = Group1;
  m_pGroups[1] = Group2;
}


void ARestraint::SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3) {
//-----------------------------------------------------------------
// set three groups of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 3);
  m_pGroups[0] = Group1;
  m_pGroups[1] = Group2;
  m_pGroups[2] = Group3;
}


void ARestraint::SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3, AGroup& Group4) {
//-----------------------------------------------------------------
// set four groups of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 4);
  m_pGroups[0] = Group1;
  m_pGroups[1] = Group2;
  m_pGroups[2] = Group3;
  m_pGroups[3] = Group4;
}

  
double ARestraint::GetAngle(AVector& A, AVector& B, AVector& C) {
//-----------------------------------------------------------------
// determine the angle formed by the points A-B-C
//-----------------------------------------------------------------
  double u;
  double a = B.Dist(C);
  double b = A.Dist(C);
  double c = A.Dist(B);

  u = (a*a + c*c - b*b) / (2.0*a*c);
  // protect against acos(<-1.0) and acos(>1.0)
  if (u < -1.0) {u = -1.0;}
  if (u >  1.0) {u =  1.0;}
  return(acos(u));
}


double ARestraint::GetDihe(AVector& A, AVector& B, AVector& C, AVector& D) {
//-----------------------------------------------------------------
// determine the dihedral angle formed by the points A-B-C-D
//-----------------------------------------------------------------
  AVector CD(D - C);
  AVector CB(B - C);
  AVector BC(C - B);
  AVector BA(A - B);
  AVector CDxCB, BCxBA;
  double  top, bot, cos_u, sin_u, Angle;
  AVector topVec;

  CDxCB = CD.cross(CB);
  BCxBA = BC.cross(BA);

  top = CDxCB.dot(BCxBA);
  bot = CDxCB.Dist() * BCxBA.Dist();
  cos_u = top/bot;

  // protect against acos(<-1.0) and acos(>1.0)
  if (cos_u < -1.0) {cos_u = -1.0;}
  if (cos_u >  1.0) {cos_u =  1.0;}

  topVec = CDxCB.cross(BCxBA);
  sin_u = (topVec/bot).dot(CB/CB.Dist());

  // protect against asin(<-1.0) and asin(>1.0)
  if (sin_u < -1.0) {sin_u = -1.0;}
  if (sin_u >  1.0) {sin_u =  1.0;}

  Angle = atan2(sin_u, cos_u);
  return(Angle);
}


void ARestraint::DistributeForce(int WhichGroup, AVector Force,
                                 GlobalMasterFreeEnergy& CFE) {
//----------------------------------------------------------------------
// Distribute Force among the group of atoms specified by WhichGroup
//
// note:  m_pGroups points to an array of Groups
//        m_pGroups[WhichGroup] references one of the Groups
//        m_pGroups[WhichGroup][i] returns an AtomID from the Group
//        (operator[] is defined to return an item from the Group)
//----------------------------------------------------------------------
  int     i, AtomID, NumAtoms, RetVal;
  double  Mass, TotalMass=0;
  AVector SmallForce;
  Vector  NAMD_Vector;

  ASSERT( (WhichGroup>=0) && (WhichGroup<m_NumGroups) );

  // calculate the total mass for the group
  NumAtoms = m_pGroups[WhichGroup].GetNumInGroup();
  for (i=0; i<NumAtoms; i++) {
    AtomID = m_pGroups[WhichGroup][i];
    Mass = CFE.getMass(AtomID);
    if (Mass < 0) {EarlyExit("Negative Mass", AtomID);};
    TotalMass += Mass;
  }

  // distribute Force according to mass of each atom in the group
  for (i=0; i<NumAtoms; i++) {
    AtomID = m_pGroups[WhichGroup][i];
    Mass = CFE.getMass(AtomID);
    if (Mass < 0) {EarlyExit("Negative Mass", AtomID);}
    SmallForce = Force * (Mass/TotalMass);
    // cast SmallForce to a NAMD-type vector (addForce uses Vector class)
    SetEqual(NAMD_Vector, SmallForce);
    RetVal = CFE.addForce(AtomID, NAMD_Vector);
    if (RetVal < 0) {EarlyExit("Can't add Force", AtomID);}
  }
}


void ARestraint::UpdateCOMs(GlobalMasterFreeEnergy& CFE) {
//-----------------------------------------------------------------
// calculate the center-of-mass of each group of atoms
//
// note:  m_pGroups points to an array of Groups
//        m_pGroups[i] references one of the Groups
//        m_pGroups[i][j] returns an AtomID from the Group
//        (operator[] is defined to return an item from the Group)
//-----------------------------------------------------------------
  int      i, j, AtomID, RetVal;
  Vector   NAMD_Vector;
  AVector  COM, Pos;
  double   Mass, TotalMass;

  ASSERT(m_NumGroups > 0);
  // for each group of atoms
  for (i=0; i<m_NumGroups; i++) {
    TotalMass = 0;
    COM.Set(0,0,0);
    // for each atom in the group
    for (j=0; j<m_pGroups[i].GetNumInGroup(); j++) {
      AtomID = m_pGroups[i][j];
      // get its position, weight position with atom's mass
      RetVal = CFE.getPosition(AtomID, NAMD_Vector);
      if (RetVal < 0) {EarlyExit("Can't get Position", AtomID);}
      // cast NAMD_Vector to AVector (getPosition uses Vector class)
      SetEqual(Pos, NAMD_Vector);
      Mass = CFE.getMass(AtomID);
      if (Mass < 0) {EarlyExit("Negative Mass", AtomID);}
      TotalMass += Mass;
      COM += Pos * Mass;
    }
    m_pCOMs[i] = COM / TotalMass;
  }
}


APosRestraint::APosRestraint() {
//-----------------------------------------------------------------
// each APosRestraint restrains 1 group of atoms to a location
//-----------------------------------------------------------------
  m_NumGroups = 1;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new AVector[m_NumGroups];
}


void APosRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the position for this position restraint
//--------------------------------------------------------------------
  double  Distance;
  char    Str[20];

  Distance = GetDistance();
  sprintf(Str, "%7.3f", Distance);

#if defined(_VERBOSE_PMF)
  iout << "Position = ";
  m_pCOMs[0].Out();
  iout << "  Target = ";
  GetPosTarget().Out();
  iout << "  Distance = ";
  iout << Str;
  iout << std::endl << endi;
#else
  m_pCOMs[0].Out();
  iout << "  ";
  GetPosTarget().Out();
  iout << "  ";
  iout << Str;
  iout << " | ";
#endif
}


double APosRestraint::GetE(AVector RefPos, double LambdaKf) {
//--------------------------------------------------------------------
// calculate and return the Energy for this position restraint.
//
//     E = (Kf/2) * (|ri - rref|)**2
//
// where |ri - rref| is the distance between a) the center-of-mass
// of the restrained atoms and b) the reference position.
//
// Note:  COM is calculated before this routine is called.
//--------------------------------------------------------------------
  return( ((m_Kf*LambdaKf)/2.0) * m_pCOMs[0].DistSqr(RefPos) );
}


AVector APosRestraint::GetGrad(int /* WhichGroup */, 
                               AVector RefPos, double LambdaKf) {
//-------------------------------------------------------------------------
// calculate and return the gradient for this position restraint.
//
//     E = (Kf/2) * (|ri - rref|)**2
//
// return:  grad(E)
//
// Notes: COM is calculated before this routine is called.
//        m_pCOMs points to an array of vectors
//        m_pCOMS[0] references the only COM for a position restraint
//        m_pCOMS[0][0] returns the x value from the vector (x,y,z)
//        (operator[] is defined to return an item from the vector)
//-------------------------------------------------------------------------
  // WhichGroup = 0;  // don't care -- there's only 1 atom restrained
  AVector Vec(m_pCOMs[0][0] - RefPos[0],
              m_pCOMs[0][1] - RefPos[1],
              m_pCOMs[0][2] - RefPos[2]);
  return(Vec*m_Kf*LambdaKf);
}


ADistRestraint::ADistRestraint() {
//-----------------------------------------------------------------------
// each ADistRestraint restrains the distance between 2 groups of atoms
//-----------------------------------------------------------------------
  m_NumGroups = 2;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new AVector[m_NumGroups];
}


void ADistRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the distance for this distance restraint
//--------------------------------------------------------------------
  double  Distance;
  char    Str1[20], Str2[20];

  Distance = m_pCOMs[0].Dist(m_pCOMs[1]);
  sprintf(Str1, "%7.3f", Distance);
  Distance = GetDistTarget();
  sprintf(Str2, "%7.3f", Distance);

#if defined(_VERBOSE_PMF)
  iout << "Distance = ";
  iout << Str1;
  iout << "  Target = ";
  iout << Str2;
  iout << std::endl << endi;
#else
  iout << Str1;
  iout << "  ";
  iout << Str2;
  iout << " | ";
#endif
}


double ADistRestraint::GetE(double RefDist, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the Energy for this distance restraint.
//
//     E = (Kf/2) * (di-dref)**2
//
// where di is the distance between 2 centers-of-mass of restrained atoms,
// and dref is the reference distance.
//
// Note:  COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  double Dist, Diff;

  Dist = m_pCOMs[0].Dist(m_pCOMs[1]);
  Diff = Dist - RefDist;
  return( ((m_Kf*LambdaKf)/2.0) * (Diff*Diff) );
}


AVector ADistRestraint::GetGrad(int WhichGroup,
                                double RefDist, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this distance restraint.
//
//     E = (Kf/2) * (di-dref)**2
//
// return:  grad(E)
//
// Notes: COM is calculated before this routine is called.
//        m_pCOMS[0 & 1] reference the COM's of each group of atoms
//---------------------------------------------------------------------------
  double  Dist;
  AVector Vec, A, B;

  ASSERT( (WhichGroup==0) || (WhichGroup==1) );
  if (WhichGroup == 0) {
    A = m_pCOMs[0];
    B = m_pCOMs[1];
  }
  else {
    A = m_pCOMs[1];
    B = m_pCOMs[0];
  }

  Dist = A.Dist(B);
  Vec.Set(A[0]-B[0], A[1]-B[1], A[2]-B[2]);
  Vec *= m_Kf * LambdaKf * (Dist - RefDist) / Dist;
  return(Vec);
}


AnAngleRestraint::AnAngleRestraint() {
//-----------------------------------------------------------------------
// each AnAngleRestraint restrains the angle between 3 groups of atoms
//-----------------------------------------------------------------------
  m_NumGroups = 3;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new AVector[m_NumGroups];
}


void AnAngleRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the angle for this angle restraint
//--------------------------------------------------------------------
  double  Angle;
  char    Str1[20], Str2[20];

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]) * (180/kPi);
  sprintf(Str1, "%8.3f", Angle);
  Angle = GetAngleTarget() * (180/kPi);
  sprintf(Str2, "%8.3f", Angle);

#if defined(_VERBOSE_PMF)
  iout << "Angle = ";
  iout << Str1;
  iout << " degrees";
  iout << "  Target = ";
  iout << Str2;
  iout << " degrees";
  iout << std::endl << endi;
#else
  iout << Str1;
  iout << "  ";
  iout << Str2;
  iout << " | ";
#endif
}


double AnAngleRestraint::GetE(double RefAngle, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the Energy for this angle restraint.
//
//     E = (Kf/2) * (Theta-ThetaRef)**2
//
// where Theta is the angle between 3 centers-of-mass of restrained atoms,
// m_pCOMs[0] -- m_pCOMs[1] -- m_pCOMs[2].
// ThetaRef is the reference angle.
//
// Note:  COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  double Angle, Diff;

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  Diff = Angle - RefAngle;
  return( ((m_Kf*LambdaKf)/2.0) * (Diff*Diff) );
}


AVector AnAngleRestraint::GetGrad(int WhichGroup,
                                  double RefAngle, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this angle restraint.
//
//     E = (Kf/2) * (Theta-ThetaRef)**2
//
// return:  grad(E)
//
// Notes: COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  AVector A, B, C;
  double  Angle;
  double  a, b, c;
  double  u;
  double  Const1, Const2, Const3, Const4, Const5;
  AVector Vec1, Vec2;

  ASSERT( (WhichGroup==0) || (WhichGroup==1) || (WhichGroup==2) );
  A = m_pCOMs[0];
  B = m_pCOMs[1];
  C = m_pCOMs[2];

  a = B.Dist(C);
  b = A.Dist(C);
  c = A.Dist(B);

  u = (a*a + c*c - b*b) / (2.0*a*c);

  // protect against acos(<-1.0), acos(>1.0), sqrt(<0), and divide by 0
  if (u < -1.0)       {u = -1.0;}
  if (u > kAlmostOne) {u = kAlmostOne;}
  Angle = acos(u);

  Const1 = -1.0 / sqrt(1.0 - u*u);
  Const2 = m_Kf * LambdaKf * (Angle-RefAngle) * Const1;

  Const3 = -a/(2.0*c*c) + 1.0/(a+a) + (b*b)/(2.0*a*c*c);
  Const4 = -c/(2.0*a*a) + 1.0/(c+c) + (b*b)/(2.0*c*a*a);
  Const5 = -b/(a*c);

  if (WhichGroup == 0) {
    Vec1 = (A-C) * (Const5/b);
    Vec2 = (A-B) * (Const3/c);
  }
  else if (WhichGroup == 1) {
    Vec1 = (B-A) * (Const3/c);
    Vec2 = (B-C) * (Const4/a);
  }
  else {
    Vec1 = (C-A) * (Const5/b);
    Vec2 = (C-B) * (Const4/a);
  }
  return( (Vec1+Vec2)*Const2);
}


ADiheRestraint::ADiheRestraint() {
//----------------------------------------------------------------------------
// each ADiheRestraint restrains the dihedral angle between 4 groups of atoms
//----------------------------------------------------------------------------
  m_NumGroups = 4;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new AVector[m_NumGroups];
}


void ADiheRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the dihedral angle for this dihedral restraint
//--------------------------------------------------------------------
  double  Dihedral;
  char    Str1[20], Str2[20], Str3[20];

  Dihedral = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]) * (180/kPi);
  sprintf(Str1, "%8.3f", Dihedral);
  Dihedral = GetDiheTarget1() * (180/kPi);
  sprintf(Str2, "%8.3f", Dihedral);
  Dihedral = GetDiheTarget2() * (180/kPi);
  sprintf(Str3, "%8.3f", Dihedral);

#if defined(_VERBOSE_PMF)
  iout << "Dihedral = ";
  iout << Str1;
  iout << " degrees";
  iout << "  Target = ";
  iout << Str2;
  iout << " degrees";
  if (TwoTargets()) {
    iout << " to ";
    iout << Str3;
    iout << " degrees";
  }
  iout << std::endl << endi;
#else
  iout << Str1;
  iout << "  ";
  iout << Str2;
  if (TwoTargets()) {
    iout << ", ";
    iout << Str3;
  }
  iout << " | ";
#endif
}


double ADiheRestraint::GetE(double RefAngle, double Const) {
//---------------------------------------------------------------------------
// calculate and return the Energy for this angle restraint.
//
//    E = (E0/2) * (1 - cos(Chi - ChiRef))
//
// where Chi is the dihedral angle between 4 centers-of-mass of restrained atoms,
// m_pCOMs[0] -- m_pCOMs[1] -- m_pCOMs[2] -- m_pCOMs[3].
// ChiRef is the reference angle.
//
// Note:  COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  double  Angle;

  Angle = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  return( (Const/2.0) * (1.0 - cos(Angle-RefAngle)) );
}


AVector ADiheRestraint::GetGrad(int WhichGroup,
                                double RefAngle, double Const) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this dihedral angle restraint.
//
//    E = (E0/2) * (1 - cos(Chi - ChiRef))
//
// return:  grad(E)
//
// Notes: COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  AVector A, B, C, D;

  ASSERT((WhichGroup==0)||(WhichGroup==1)||(WhichGroup==2)||(WhichGroup==3));

  if ((WhichGroup==0) || (WhichGroup==1)) {
    A = m_pCOMs[0];
    B = m_pCOMs[1];
    C = m_pCOMs[2];
    D = m_pCOMs[3];
  }
  // re-state the problem so the gradient is solved for either atoms 0 or 1
  else {
    A = m_pCOMs[3];
    B = m_pCOMs[2];
    C = m_pCOMs[1];
    D = m_pCOMs[0];
    if (WhichGroup==3) {WhichGroup=0;}
    if (WhichGroup==2) {WhichGroup=1;}
  }

  AVector CD(D - C);
  AVector CB(B - C);
  AVector BC(C - B);
  AVector BA(A - B);
  AVector AC(C - A);
  AVector CDxCB, BCxBA;
  AVector Vec;
  double  phi;
  double  top, bot, u;

  CDxCB = CD.cross(CB);
  BCxBA = BC.cross(BA);

  top = CDxCB.dot(BCxBA);
  bot = CDxCB.Dist() * BCxBA.Dist();

  u = top/bot;
  // protect against acos(<-1.0), acos(>1.0), sqrt(<0), and divide by 0
  if (u < kAlmostMinusOne) {u = kAlmostMinusOne;}
  if (u > kAlmostOne)      {u = kAlmostOne;}

  // get dihedral using atan
  phi = GetDihe(A,B,C,D);

  AVector dP1, dP2, dP3;
  AVector dP4, dP5, dP6;
  ASSERT((WhichGroup==0) || (WhichGroup==1));
  if (WhichGroup==0) {
    dP1.Set( 0,      0,      0    );
    dP2.Set( 0,      0,      0    );
    dP3.Set( 0,      0,      0    );
    dP4.Set( 0,     -BC[2],  BC[1]);
    dP5.Set( BC[2],  0,     -BC[0]);
    dP6.Set(-BC[1],  BC[0],  0    );
  }
  else {
    dP1.Set( 0,     -CD[2],  CD[1]);
    dP2.Set( CD[2],  0,     -CD[0]);
    dP3.Set(-CD[1],  CD[0],  0    );
    dP4.Set( 0,      AC[2], -AC[1]);
    dP5.Set(-AC[2],  0,      AC[0]);
    dP6.Set( AC[1], -AC[0],  0    );
  }

  Vec = gradU(CDxCB, BCxBA, dP1, dP2, dP3, dP4, dP5, dP6);
  Vec *= (Const/2.0) * sin(phi-RefAngle) * (-1.0/sqrt(1.0 - u*u));

  // flip gradient for negative angles
  if (phi < 0) {
    Vec *= -1.0;
  }
  
  return(Vec);
}


AVector ADiheRestraint::gradU(AVector& P1P2P3, AVector& P4P5P6,
                              AVector& dP1,    AVector& dP2,    AVector& dP3,
                              AVector& dP4,    AVector& dP5,    AVector& dP6) {
//----------------------------------------------------------------
// calculate the gradient for ((P1P2P3.dot.P4P5P6)/(mag(P1P2P3)*mag(P4P5P6)))
// P1P2P3 = (P1)i + (P2)j + (P3)k
// P4P5P6 = (P4)i + (P5)j + (P6)k
// dP1 = (d(P1)/dx)i + (d(P1)/dy)j +(d(P1)/dz)k
// dP2 = (d(P2)/dx)i + (d(P2)/dy)j +(d(P2)/dz)k
// dP3 = (d(P3)/dx)i + (d(P3)/dy)j +(d(P3)/dz)k
// dP4 = (d(P4)/dx)i + (d(P4)/dy)j +(d(P4)/dz)k
// dP5 = (d(P5)/dx)i + (d(P5)/dy)j +(d(P5)/dz)k
// dP6 = (d(P6)/dx)i + (d(P6)/dy)j +(d(P6)/dz)k
//----------------------------------------------------------------
  double  Mag123, Mag456, Dot;
  double  Const1, Const2, Const3;
  double  P1, P2, P3, P4, P5, P6;
  AVector RetVec;

  P1 = P1P2P3[0];  P2 = P1P2P3[1];  P3 = P1P2P3[2];
  P4 = P4P5P6[0];  P5 = P4P5P6[1];  P6 = P4P5P6[2];

  Mag123 = P1P2P3.Dist();
  Mag456 = P4P5P6.Dist();
  Dot    = P1P2P3.dot(P4P5P6);

  Const1 =         1.0 / (Mag123*Mag456);
  Const2 = -Dot * (1.0 / (Mag123*Mag456*Mag456*Mag456));
  Const3 = -Dot * (1.0 / (Mag456*Mag123*Mag123*Mag123));

  RetVec = (dP4*P1 + dP1*P4 + dP5*P2 + dP2*P5 + dP6*P3 + dP3*P6) * Const1 +
           (dP4*P4 + dP5*P5 + dP6*P6)                            * Const2 +
           (dP1*P1 + dP2*P2 + dP3*P3)                            * Const3;

  return(RetVec);
}


double AFixedPosRestraint::GetEnergy() {
//--------------------------------------------------------------------
// return the Energy for this fixed position restraint.
//--------------------------------------------------------------------
  return(GetE(m_RefPos));
}


AVector AFixedPosRestraint::GetGradient(int WhichGroup) {
//-------------------------------------------------------------------------
// return the Gradient for this fixed position restraint.
//-------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefPos));
}


double ABoundPosRestraint::GetEnergy() {
//--------------------------------------------------------------------
// calculate and return the Energy for this bound position restraint.
//
// Note:  This is an exception because the form for the E term is
//        different from the other postion restraints.
//--------------------------------------------------------------------
  double  E, Dist, Diff;

  E = 0.0;
  Dist = m_pCOMs[0].Dist(m_RefPos);
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    Diff = Dist - m_RefDist;
    E = (m_Kf/2.0) * (Diff*Diff);
  }
  return(E);
}


AVector ABoundPosRestraint::GetGradient(int /* WhichGroup */) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this bound position restraint.
//
// Note:  This is an exception because the form for the E term is
//        different from the other postion restraints.
//---------------------------------------------------------------------------
  double  Dist;
  AVector Vec;   // Vec is initialized to (0,0,0)

  // WhichGroup = 0;  // don't care -- there's only 1 atom restrained
  Dist = m_pCOMs[0].Dist(m_RefPos);
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    Vec.Set(m_pCOMs[0][0] - m_RefPos[0],
            m_pCOMs[0][1] - m_RefPos[1],
            m_pCOMs[0][2] - m_RefPos[2]);
    Vec *= m_Kf * (Dist - m_RefDist) / Dist;
  }
  return(Vec);
}


double AForcingPosRestraint::GetEnergy() {
//--------------------------------------------------------------------
// return the Energy for this forcing position restraint.
//
// rref = lambda*r1 + (1-lambda)*r0.
// where r0 is the starting position and r1 is the final position
//--------------------------------------------------------------------
  AVector  RefPos;

  RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
  return(GetE(RefPos, m_LambdaKf));
}


AVector AForcingPosRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this forcing position restraint.
//
// rref = lambda*r1 + (1-lambda)*r0.
// where r0 is the starting position and r1 is the final position
//---------------------------------------------------------------------------
  AVector  RefPos;

  RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefPos, m_LambdaKf));
}


double AForcingPosRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this forcing position restraint
//---------------------------------------------------------------------------
  AVector  RefPos;
  double   T1, T2, T3;

  RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
  T1 = (m_pCOMs[0][0] - RefPos[0]) * (m_StartPos[0] - m_StopPos[0]);
  T2 = (m_pCOMs[0][1] - RefPos[1]) * (m_StartPos[1] - m_StopPos[1]);
  T3 = (m_pCOMs[0][2] - RefPos[2]) * (m_StartPos[2] - m_StopPos[2]);
  return( m_Kf * m_LambdaKf * (T1+T2+T3) );
}


double AFixedDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed distance restraint.
//---------------------------------------------------------------------------
  return(GetE(m_RefDist));
}


AVector AFixedDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed distance restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefDist));
}


double ABoundDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound distance restraint.
//---------------------------------------------------------------------------
  double Dist, E;

  E = 0.0;
  Dist = m_pCOMs[0].Dist(m_pCOMs[1]);
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    E = GetE(m_RefDist);
  }
  return(E);
}


AVector ABoundDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound distance restraint.
//---------------------------------------------------------------------------
  double  Dist;
  AVector Vec;

  Dist = m_pCOMs[0].Dist(m_pCOMs[1]);
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    Vec = GetGrad(WhichGroup, m_RefDist);
  }
  return(Vec);
}


double AForcingDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this forcing distance restraint.
//---------------------------------------------------------------------------
  double  RefDist;

  RefDist = m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef);
  return(GetE(RefDist, m_LambdaKf));
}


AVector AForcingDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this forcing distance restraint.
//---------------------------------------------------------------------------
  double  RefDist;
  
  RefDist = m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefDist, m_LambdaKf));
}


double AForcingDistRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this forcing distance restraint
//---------------------------------------------------------------------------
  double  Dist;
  double  RefDist;
  
  Dist = m_pCOMs[0].Dist(m_pCOMs[1]);
  RefDist = m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef);
  return( m_Kf * m_LambdaKf * (Dist-RefDist)*(m_StartDist-m_StopDist) );
}


double AFixedAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed angle restraint.
//---------------------------------------------------------------------------
  return(GetE(m_RefAngle));
}


AVector AFixedAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed angle restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefAngle));
}


double ABoundAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound angle restraint.
//---------------------------------------------------------------------------
  double  E, Angle;

  E = 0.0;
  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  if (((m_Bound==kUpper) && (Angle>m_RefAngle)) ||
      ((m_Bound==kLower) && (Angle<m_RefAngle))) {
    E = GetE(m_RefAngle);
  }
  return(E);
}


AVector ABoundAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound angle restraint
//---------------------------------------------------------------------------
  double  Angle;
  AVector Vec;

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  if (((m_Bound==kUpper) && (Angle>m_RefAngle)) ||
      ((m_Bound==kLower) && (Angle<m_RefAngle))) {
    Vec = GetGrad(WhichGroup, m_RefAngle);
  }
  return(Vec);
}


double AForcingAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this forcing angle restraint.
//---------------------------------------------------------------------------
  double  RefAngle;

  RefAngle = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetE(RefAngle, m_LambdaKf));
}


AVector AForcingAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this forcing angle restraint.
//---------------------------------------------------------------------------
  double  RefAngle;

  RefAngle = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefAngle, m_LambdaKf));
}


double AForcingAngleRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this forcing angle restraint
//---------------------------------------------------------------------------
  double  Angle;
  double  RefAngle;

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  RefAngle = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return( m_Kf * m_LambdaKf * (Angle-RefAngle)*(m_StartAngle-m_StopAngle) );
}


double AFixedDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed dihedral angle restraint.
//---------------------------------------------------------------------------
  return(GetE(m_RefAngle, m_Kf));
}


AVector AFixedDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed dihedral angle restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefAngle, m_Kf));
}


double ABoundDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound dihedral angle restraint.
//---------------------------------------------------------------------------
  double  E, Dihe, Const;

  Const = m_Kf / (1.0 - cos(m_IntervalAngle));
  Dihe = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  // dihedral angle is between LowerAngle and UpperAngle
  if ( (Dihe>m_LowerAngle) && (Dihe<m_UpperAngle) ) {
    E = 0.0;
  }
  // dihedral angle is between LowerAngle and LowerAngle-IntervalAngle
  else if ( (Dihe<m_LowerAngle) && (Dihe>(m_LowerAngle-m_IntervalAngle)) ) {
    E = GetE(m_LowerAngle, Const);
  }
  // dihedral angle is between UpperAngle and UpperAngle+IntervalAngle
  else if ( (Dihe>m_UpperAngle) && (Dihe<(m_UpperAngle+m_IntervalAngle)) ) {
    E = GetE(m_UpperAngle, Const);
  }
  // dihedral angle is more than UpperAngle or less than LowerAngle
  else {
    E = Const;
  }
  return(E);
}


AVector ABoundDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound dihedral angle restraint.
//---------------------------------------------------------------------------
  AVector Vec;
  double  Dihe, Const;

  Const = m_Kf / (1.0 - cos(m_IntervalAngle));
  Dihe = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  // dihedral angle is between LowerAngle and LowerAngle-IntervalAngle
  if ( (Dihe<m_LowerAngle) && (Dihe>(m_LowerAngle-m_IntervalAngle)) ) {
    Vec = GetGrad(WhichGroup, m_LowerAngle, Const);
  }
  // dihedral angle is between UpperAngle and UpperAngle+IntervalAngle
  else if ( (Dihe>m_UpperAngle) && (Dihe<(m_UpperAngle+m_IntervalAngle)) ) {
    Vec = GetGrad(WhichGroup, m_UpperAngle, Const);
  }
  return(Vec);
}


double AForcingDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this forcing dihedral angle restraint.
//---------------------------------------------------------------------------
  double  RefDihe;

  RefDihe = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetE(RefDihe, m_Kf*m_LambdaKf));
}


AVector AForcingDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this forcing dihedral angle restraint.
//---------------------------------------------------------------------------
  double  RefDihe;
  
  RefDihe = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefDihe, m_Kf*m_LambdaKf));
}


double AForcingDiheRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this forcing dihedral angle restraint
//---------------------------------------------------------------------------
  double  Dihe;
  double  RefDihe;
  
  Dihe = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  RefDihe = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return((m_Kf/2)*m_LambdaKf * sin(Dihe-RefDihe) * (m_StartAngle-m_StopAngle));
}

