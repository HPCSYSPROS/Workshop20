/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <memory.h>
//#include <iomanip.h>
#include <stdio.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"
#include "FreeEnergyGroup.h"

ALambdaManager::ALambdaManager() {
//------------------------------------------------------------------------
// make room for some LambdaControl objects.
//------------------------------------------------------------------------
  m_ActiveIndex = 0;
  m_NumObjects = 0;
  m_pPmfBlocks = new ALambdaControl[kLambdaNumToStart];
  m_MaxNum = kLambdaNumToStart;
}


ALambdaManager::~ALambdaManager() {
//------------------------------------------------------------------------
// return borrowed memory to the free store.
//------------------------------------------------------------------------
  delete []m_pPmfBlocks;
}


int  ALambdaManager::GetNum_dU_dLambda() {
//------------------------------------------------------------------------
// get number times dU_dLambda accumulated from active lambda object.
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].GetNum_dU_dLambda());
}


void ALambdaManager::Integrate_MCTI() {
//------------------------------------------------------------------------
// integrate <dU/dLambda> for MCTI
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  (*this)[m_ActiveIndex].Integrate_MCTI();
}


void ALambdaManager::Accumulate(double dU_dLambda) {
//------------------------------------------------------------------------
// add dU_dLambda to the active accumulator (in one of the lambda objects)
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  (*this)[m_ActiveIndex].Accumulate(dU_dLambda);
}


double ALambdaManager::GetAccumulation() {
//------------------------------------------------------------------------
// get the accumulation of dU_dLambda from active lambda object.
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].GetAccumulation());
}


double ALambdaManager::GetIntegration() {
//------------------------------------------------------------------------
// get accumulation of <dU_dLambda> * dLambda from active lambda object.
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].GetIntegration());
}


void ALambdaManager::ZeroAccumulator() {
//------------------------------------------------------------------------
// zero accumulation of dU_dLambda in the active lambda object.
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  (*this)[m_ActiveIndex].ZeroAccumulator();
}


Bool_t ALambdaManager::IsFirstStep() {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if it's time to print restraint information
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsFirstStep());
}


Bool_t ALambdaManager::IsTimeToPrint() {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if it's time to print restraint information
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsTimeToPrint());
}


Bool_t ALambdaManager::IsTimeToPrint_dU_dLambda() {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if it's time to print du/dLambda information
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsTimeToPrint_dU_dLambda());
}


Bool_t ALambdaManager::IsTimeToClearAccumulator() {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if it's time to start accumulating dU/dLambda from zero
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsTimeToClearAccumulator());
}


Bool_t ALambdaManager::IsEndOf_MCTI_Step() {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if this is the last time step of an MCTI step
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsEndOf_MCTI_Step());
}


Bool_t ALambdaManager::IsEndOf_MCTI() {
//------------------------------------------------------------------------
// ASSUMING that the m_ActiveIndex'th LambdaControl is active,
// decide if this is the last time step of an MCTI block
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].IsEndOf_MCTI());
}


void ALambdaManager::PrintLambdaHeader(double dT) {
//------------------------------------------------------------------------
// print header for a new lambda control object
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  (*this)[m_ActiveIndex].PrintLambdaHeader(dT);
}


void ALambdaManager::PrintHeader(double dT) {
//------------------------------------------------------------------------
// print information about current time step
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  (*this)[m_ActiveIndex].PrintHeader(dT);
}


int ALambdaManager::GetNumStepsSoFar() {
//------------------------------------------------------------------------
// return the number of steps taken in the active LambdaControl block
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].GetNumStepsSoFar());
}


int ALambdaManager::GetNumAccumStepsSoFar() {
//------------------------------------------------------------------------
// return the total number of steps dU/dLambda has been accumulated
//------------------------------------------------------------------------
  ASSERT((*this)[m_ActiveIndex].IsActive());
  return((*this)[m_ActiveIndex].GetNumAccumStepsSoFar());
}


void ALambdaManager::PrintSomeSpaces() {
//------------------------------------------------------------------------
// some stuff to make the output look nice
//------------------------------------------------------------------------
#if !defined(_VERBOSE_PMF)
  iout << "                    ";
#endif
}


void ALambdaManager::Print_dU_dLambda_Summary(double Sum_dU_dLambdas) {
//------------------------------------------------------------------------
// print sum of dU/dLambda's for current time-step and
// the accumulation of the above.
//------------------------------------------------------------------------
  char  Str[100];

#if defined(_VERBOSE_PMF)
  iout << "FreeEnergy: ";
  iout << "For all forcing restraints, dU/dLambda  = ";
  iout << Sum_dU_dLambdas << std::endl << endi;
  iout << "FreeEnergy: ";
  iout << "For all forcing restraints, Free Energy = ";
  iout << GetAccumulation();
  iout << " for " << GetNum_dU_dLambda() << " steps" << std::endl << endi;
#else
  sprintf(Str, "%10.2e", GetAccumulation());
  iout << Str << "  ";
  sprintf(Str, "%6d", GetNum_dU_dLambda());
  iout << Str << "  ";
#endif
}


void ALambdaManager::Print_MCTI_Integration() {
//------------------------------------------------------------------------
// print the integral of: <dU/dLambda> * dLambda
//------------------------------------------------------------------------
  iout << "FreeEnergy: ";
  iout << "For MCTI, Free Energy Integral = ";
  iout << GetIntegration();
  iout << " for " << GetNumAccumStepsSoFar() << " steps" << std::endl << endi;
}


Bool_t ALambdaManager::GetLambdas(double& LambdaKf, double& LambdaRef) {
//------------------------------------------------------------------------
// get LambdaKf and LambdaRef from the active LambdaControl
//
// return(kTrue) if an active LambdaControl is found
// return(kFalse) if all LambdaControls have expired
//------------------------------------------------------------------------
  // don't continue if all LamdaControl's have expired
  if (m_ActiveIndex == m_NumObjects) {
    return(kFalse);
  }

  // if the m_ActiveIndex'th LambdaControl is no longer active
  if ( !(*this)[m_ActiveIndex].IsActive()) {
    // move on to the next one
    m_ActiveIndex++;
    // if there is no next object, return KFalse
    if (m_ActiveIndex == m_NumObjects) {
      return(kFalse);
    }
    // otherwise, make sure the next one's active
    else {
      ASSERT( (*this)[m_ActiveIndex].IsActive() );
    }
  }
  // return LambdaKf and LambdaRef from the active LambdaControl
  LambdaKf =  (*this)[m_ActiveIndex].GetLambdaKf();
  LambdaRef = (*this)[m_ActiveIndex].GetLambdaRef();
  return(kTrue);
}


void ALambdaManager::Clear() {
//------------------------------------------------------------------------
// leave memory allocation alone.
//------------------------------------------------------------------------
  m_NumObjects = 0;
}


int ALambdaManager::Add(ALambdaControl& PmfBlock) {
//------------------------------------------------------------------------
// add an object to the list.  if there's not enough room, make room.
// return an index to the added oject.
//------------------------------------------------------------------------
  ALambdaControl*  pPmfBlocks;

  // if there's no room for another object
  if (m_NumObjects == m_MaxNum) {
    // create an array with more space
    m_MaxNum *= kLambdaMultiplier;
    pPmfBlocks = new ALambdaControl[m_MaxNum];
    // copy from the full array to the new one
    for (int i=0; i<m_NumObjects; i++) {
      pPmfBlocks[i] = m_pPmfBlocks[i];
    }
    // return the space used for the full array
    delete []m_pPmfBlocks;
    // point to the bigger array
    m_pPmfBlocks = pPmfBlocks;
  }
  // add the object to the array
  m_pPmfBlocks[m_NumObjects] = PmfBlock;
  m_NumObjects++;
  return(m_NumObjects-1);
}


ALambdaControl& ALambdaManager::operator[] (int Index) {
//------------------------------------------------------------------------
// return an object from this group of objects.
//------------------------------------------------------------------------
  ASSERT((Index>=0) && (Index<m_NumObjects));
  return(m_pPmfBlocks[Index]);
}


int ALambdaManager::GetTotalNumSteps() {
//------------------------------------------------------------------------
// calculate and return the total number of steps needed for all
// pmf and mcti blocks
//------------------------------------------------------------------------
  int  Total, i;
  
  Total = 0;
  for (i=0; i<m_NumObjects; i++) {
    Total += m_pPmfBlocks[i].GetNumSteps();
  }
  return(Total);
}

