/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <memory.h>
#include <string.h>
// #include <iomanip.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
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

ARestraintManager::ARestraintManager() {
//------------------------------------------------------------------------
// allocate space for restraint POINTERS
//------------------------------------------------------------------------
  m_ppRestraints = new pRestr[kNumToStart];
  m_NumRestraints = 0;
  m_MaxNum = kNumToStart;
}


ARestraintManager::~ARestraintManager() {
//------------------------------------------------------------------------
// when the RestraintManager goes out of scope,
// free the space that was allocated for each restraint
// (see: ARestraint* GetRestraint(char*, int&)
// then, free the space that was allocated for the pointers
//------------------------------------------------------------------------
  for (int i=0; i<m_NumRestraints; i++) {
    delete m_ppRestraints[i];
  }
  delete []m_ppRestraints;
}


ARestraint* ARestraintManager::operator[] (int Index) {
//------------------------------------------------------------------------
// get a pointer
//------------------------------------------------------------------------
  ASSERT( (Index>=0) && (Index<m_NumRestraints) );
  return(m_ppRestraints[Index]);
}


void ARestraintManager::Add(ARestraint* pRestraint) {
//------------------------------------------------------------------------
// add a pointer to the list.  if there's not enough room, make room.
//------------------------------------------------------------------------
  ARestraint**  ppRestraints;

  // if there's no room for a new pointer
  if (m_NumRestraints == m_MaxNum) {
    // create an array with more space
    m_MaxNum *= kMultiplier;
    ppRestraints = new pRestr[m_MaxNum];
    // fast copy from the full array to the new one (memcpy(dest, src, bytes))
    memcpy(ppRestraints, m_ppRestraints, sizeof(ARestraint*)*m_NumRestraints);
    // return the space used for the full array
    delete []m_ppRestraints;
    // point to the bigger array
    m_ppRestraints = ppRestraints;
  }

  // add the int to the int array
  m_ppRestraints[m_NumRestraints] = pRestraint;
  m_NumRestraints++;
}


void ARestraintManager::UpdateCOMs(GlobalMasterFreeEnergy& CFE) {
//------------------------------------------------------------------------
// update the centers-of-mass in each restraint in the list
// note:  m_ppRestraints points to a list of restraint pointers
//        m_ppRestraints[i] references one of the pointers
//        m_ppRestraints[i]-> accesses a restraint function
//------------------------------------------------------------------------
  for (int i=0; i<m_NumRestraints; i++) {
    m_ppRestraints[i]->UpdateCOMs(CFE);
  }
}


void ARestraintManager::AddForces(GlobalMasterFreeEnergy& CFE) {
//---------------------------------------------------------------------------
// for each restraint, apply restraining force to each COM (center-of-mass).
//---------------------------------------------------------------------------
  int      i, j, NumCOMs;
  AVector  Force;

  // for each restraint
  for (i=0; i<m_NumRestraints; i++) {
    // for each center-of-mass
    NumCOMs = m_ppRestraints[i]->GetNumGroups();
    for (j=0; j<NumCOMs; j++) {
      Force = m_ppRestraints[i]->GetGradient(j);
      // apply restraining force in opposite direction from gradient
      Force *= -1.0;
      m_ppRestraints[i]->DistributeForce(j, Force, CFE);
    }
  }
}


double ARestraintManager::Sum_dU_dLambdas() {
//---------------------------------------------------------------------------
// sum up dU/dLambda from each forcing restraint
//---------------------------------------------------------------------------
  double Sum=0;

  for (int i=0; i<m_NumRestraints; i++) {
    Sum += m_ppRestraints[i]->Get_dU_dLambda();
  }
  return(Sum);
}


Bool_t ARestraintManager::ThereIsAForcingRestraint() {
//---------------------------------------------------------------------------
// return kTrue if there's at least one forcing restraint
//---------------------------------------------------------------------------
  for (int i=0; i<m_NumRestraints; i++) {
    if (m_ppRestraints[i]->IsForcing()) {
      return(kTrue);
    }
  }
  return(kFalse);
}


void ARestraintManager::PrintEnergyInfo() {
//---------------------------------------------------------------------------
// for a restraint, print restraint type and Energy.
//---------------------------------------------------------------------------
#if defined(_VERBOSE_PMF)
  for (int i=0; i<m_NumRestraints; i++) {
    PrintPreInfo(i);
    iout << "Energy = ";
    iout << m_ppRestraints[i]->GetEnergy() << std::endl << endi;
  }
#endif
}


void ARestraintManager::PrintRestraintInfo() {
//---------------------------------------------------------------------------
// for a restraint, print its position, distance, angle, or dihedral angle.
//---------------------------------------------------------------------------
  for (int i=0; i<m_NumRestraints; i++) {
#if defined(_VERBOSE_PMF)
    PrintPreInfo(i);
#endif
    m_ppRestraints[i]->PrintInfo();
  }
#if !defined(_VERBOSE_PMF)
  iout << std::endl << endi;
#endif
}


void ARestraintManager::Print_dU_dLambda_Info() {
//---------------------------------------------------------------------------
// if restraint is a forcing restraint, print dU/dLambda.
//---------------------------------------------------------------------------
#if defined(_VERBOSE_PMF)
  for (int i=0; i<m_NumRestraints; i++) {
    if (m_ppRestraints[i]->IsForcing()) {
      PrintPreInfo(i);
      iout << "dU/dLambda = ";
      iout << m_ppRestraints[i]->Get_dU_dLambda() << std::endl << endi;
    }
  }
#endif
}


void ARestraintManager::PrintPreInfo(int Index) {
//---------------------------------------------------------------------------
// print "Restraint xxx:  Type of Restraint:  "
//---------------------------------------------------------------------------
  char  Str[100];
  char  NumStr[20];

  sprintf(NumStr, "%3d", Index+1);
  iout << "FreeEnergy: " << "Restraint " << NumStr << ":  ";
  m_ppRestraints[Index]->GetStr(Str);
  iout << Str << ":  ";
}

