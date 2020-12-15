/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <string.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyLambda.h"
#include "Vector.h"
#include "FreeEnergyVector.h"

// initialize static member variables
int  ALambdaControl::m_CurrStep = 0;

ALambdaControl::ALambdaControl() {
//------------------------------------------------------------------------
// initialize member variables
//------------------------------------------------------------------------
  // initialize these for the first block.
  m_StartStep = 0;
  m_LambdaKf  = 1.0;         // applied to Kf
  m_LambdaRef = 0.0;         // applied to reference (pos, dist, angle, dihe)
  m_Sum_dU_dLambda = 0.0;    // for free energy measurement
  m_Num_dU_dLambda = 0;      // "
  m_MCTI_Integration = 0.0;  // "

  // set the rest to illegal values.
  m_Task =           kUnknownTask;
  m_NumSteps =      -1;
  m_NumEquilSteps = -1;
  m_NumAccumSteps = -1;
  m_NumPrintSteps = -1;
  m_NumRepeats =    -1;
  m_StopStep =      -1;
}


void ALambdaControl::Init(ALambdaControl& PriorBlock) {
//------------------------------------------------------------------------
// initialize this object using the settings for the prior object.
//------------------------------------------------------------------------
  m_Task =          PriorBlock.m_Task;
  m_NumSteps =      PriorBlock.m_NumSteps;
  m_NumEquilSteps = PriorBlock.m_NumEquilSteps;
  m_NumAccumSteps = PriorBlock.m_NumAccumSteps;
  m_NumPrintSteps = PriorBlock.m_NumPrintSteps;
  m_NumRepeats =    PriorBlock.m_NumRepeats;
  switch (PriorBlock.m_Task) {

    case kUp:       m_LambdaKf=1.0;  m_LambdaRef=1.0; break;
    case kStepUp:   m_LambdaKf=1.0;  m_LambdaRef=1.0; break;
    case kDown:     m_LambdaKf=1.0;  m_LambdaRef=0.0; break;
    case kStepDown: m_LambdaKf=1.0;  m_LambdaRef=0.0; break;
    case kStop:     m_LambdaKf=1.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;

    case kGrow:     m_LambdaKf=1.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kStepGrow: m_LambdaKf=1.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kFade:     m_LambdaKf=0.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kStepFade: m_LambdaKf=0.0;  m_LambdaRef=PriorBlock.m_LambdaRef; break;
    case kNoGrow:   m_LambdaKf =  PriorBlock.m_LambdaKf;
                    m_LambdaRef = PriorBlock.m_LambdaRef;  break;
    default:        ASSERT(kFalse); break;  //should never get here

  }
  m_StartStep = PriorBlock.GetLastStep();
}


ALambdaControl& ALambdaControl::operator= (ALambdaControl& PmfBlock) {
//------------------------------------------------------------------------
// copy everything from PmfBlock to this block.
//------------------------------------------------------------------------
  m_NumSteps =      PmfBlock.m_NumSteps;
  m_NumEquilSteps = PmfBlock.m_NumEquilSteps;
  m_NumAccumSteps = PmfBlock.m_NumAccumSteps;
  m_NumRepeats =    PmfBlock.m_NumRepeats;
  m_NumPrintSteps = PmfBlock.m_NumPrintSteps;
  m_StartStep =     PmfBlock.m_StartStep;
  m_StopStep =      PmfBlock.m_StopStep;
  m_LambdaKf =      PmfBlock.m_LambdaKf;
  m_LambdaRef =     PmfBlock.m_LambdaRef;
  m_Task =          PmfBlock.m_Task;
  return(*this);
}


void ALambdaControl::Accumulate(double dU_dLambda) {
//------------------------------------------------------------------------
// if lambda is constant, sum dU/dLambda
// if lambda is changing, sum dU/dLambda * dLambda
//------------------------------------------------------------------------

  m_Num_dU_dLambda++;
  switch (m_Task) {
    // lambda is constant
    case kStop:      case kNoGrow:    case kStepUp:
    case kStepDown:  case kStepGrow:  case kStepFade:
      m_Sum_dU_dLambda += dU_dLambda;
      break;
    // lambda is increasing
    case kUp:
    case kGrow:
      m_Sum_dU_dLambda += dU_dLambda / (double)m_NumSteps;
      break;
    // lambda is decreasing
    case kDown:
    case kFade:
      m_Sum_dU_dLambda -= dU_dLambda / (double)m_NumSteps;
      break;
    // should never get here
    default:
      ASSERT(kFalse);
      break;
  }
}


void ALambdaControl::Integrate_MCTI() {
//------------------------------------------------------------------------
// integrate MCTI:  <dU/dLambda> * dLambda
//------------------------------------------------------------------------
  ASSERT(m_Task==kStepUp   || m_Task==kStepDown || 
         m_Task==kStepGrow || m_Task==kStepFade);

  switch (m_Task) {
    // lambda is increasing
    case kStepUp:
    case kStepGrow:
      m_MCTI_Integration += GetAccumulation() / (double)m_NumRepeats;
      break;
    // lambda is decreasing
    case kStepDown:
    case kStepFade:
      m_MCTI_Integration -= GetAccumulation() / (double)m_NumRepeats;
      break;
    // should never get here
    default:
      ASSERT(kFalse);
      break;
  }
}


double ALambdaControl::GetIntegration() {
//------------------------------------------------------------------------
// get MCTI integral.  integral(dU/dLambda> * dLambda
//------------------------------------------------------------------------
  ASSERT(m_Task==kStepUp   || m_Task==kStepDown || 
         m_Task==kStepGrow || m_Task==kStepFade);

  return(m_MCTI_Integration);
}


double ALambdaControl::GetAccumulation() {
//------------------------------------------------------------------------
// if lambda is constant, return (sum/N)
// if lambda is changing, return (sum)
//------------------------------------------------------------------------
  switch (m_Task) {
    // lambda is constant
    case kStop:    case kNoGrow:
    case kStepUp:  case kStepDown:  case kStepGrow:  case kStepFade:
      ASSERT(m_Num_dU_dLambda != 0);
      return(m_Sum_dU_dLambda / (double)m_Num_dU_dLambda);
    // lambda is changing
    case kUp:  case kDown:  case kGrow:  case kFade:
      return(m_Sum_dU_dLambda);
    // should never get here
    default:
      ASSERT(kFalse);
      return(0);
  }
}


int ALambdaControl::GetNumAccumStepsSoFar() {
//------------------------------------------------------------------------
// return the total number of steps dU/dLambda has been accumulated
// this is only called during MCTI
//------------------------------------------------------------------------
  ASSERT(m_Task==kStepUp   || m_Task==kStepDown || 
         m_Task==kStepGrow || m_Task==kStepFade);

  int Count, Maybe;

  Count = (m_CurrStep-m_StartStep) / (m_NumAccumSteps+m_NumEquilSteps);
  Count *= m_NumAccumSteps;
  Maybe = (m_CurrStep-m_StartStep) % (m_NumAccumSteps+m_NumEquilSteps);
  Maybe -= m_NumEquilSteps;
  if (Maybe > 0) {
    Count += Maybe;
  }
  return(Count);
}


int ALambdaControl::GetNumSteps() {
//------------------------------------------------------------------------
// get the number of steps needed for this pmf or mcti block
//------------------------------------------------------------------------
  // make sure m_StopStep is calculated
  GetLastStep();
  
  return( (m_StopStep - m_StartStep) + 1 );
}


Bool_t ALambdaControl::IsLastStep() {
//------------------------------------------------------------------------
// return true if we're on the last step of this pmf or mcti block
//------------------------------------------------------------------------
	if (m_CurrStep == GetLastStep()) {
    return(kTrue);
  }
  else {
    return(kFalse);
  }
}


int ALambdaControl::GetLastStep() {
//------------------------------------------------------------------------
// get the last step of this task
//------------------------------------------------------------------------
  // if it's already calculated, just return it
  if (m_StopStep > 0) {
    return(m_StopStep);
  }
  // otherwise calculate it
  switch (m_Task) {
    case kStepUp:
    case kStepDown:
    case kStepGrow:
    case kStepFade:
      m_StopStep = m_StartStep +
                  (m_NumAccumSteps+m_NumEquilSteps) * m_NumRepeats;
      break;
    default:
      m_StopStep = m_StartStep + m_NumSteps;
      break;
  }
  // and return it
  return(m_StopStep);
}


void ALambdaControl::GetPaddedTaskStr(char* Str) {
//------------------------------------------------------------------------
// get a string that describes this task
//------------------------------------------------------------------------
  switch (m_Task) {
    case kUp:        strcpy(Str, "      Up");      break;
    case kDown:      strcpy(Str, "    Down");      break;
    case kStop:      strcpy(Str, "    Stop");      break;
    case kGrow:      strcpy(Str, "    Grow");      break;
    case kFade:      strcpy(Str, "    Fade");      break;
    case kNoGrow:    strcpy(Str, "  NoGrow");      break;
    case kStepUp:    strcpy(Str, "  StepUp");      break;
    case kStepDown:  strcpy(Str, "StepDown");      break;
    case kStepGrow:  strcpy(Str, "StepGrow");      break;
    case kStepFade:  strcpy(Str, "StepFade");      break;
    default:         strcpy(Str, "Bug Alert!!!");  break;
  }
}


void ALambdaControl::GetTaskStr(char* Str) {
//------------------------------------------------------------------------
// get a string that describes this task
//------------------------------------------------------------------------
  switch (m_Task) {
    case kUp:        strcpy(Str, "Up");            break;
    case kDown:      strcpy(Str, "Down");          break;
    case kStop:      strcpy(Str, "Stop");          break;
    case kGrow:      strcpy(Str, "Grow");          break;
    case kFade:      strcpy(Str, "Fade");          break;
    case kNoGrow:    strcpy(Str, "NoGrow");        break;
    case kStepUp:    strcpy(Str, "StepUp");        break;
    case kStepDown:  strcpy(Str, "StepDown");      break;
    case kStepGrow:  strcpy(Str, "StepGrow");      break;
    case kStepFade:  strcpy(Str, "StepFade");      break;
    default:         strcpy(Str, "Bug Alert!!!");  break;
  }
}


Bool_t ALambdaControl::IsTimeToClearAccumulator() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to start accumulating dU/dLambda from zero.
// (clear the accumulator on the 1st step)
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  switch (m_Task) {
    // for pmf blocks, clear accumulator on first step
    case kUp:      case kDown:    case kStop:
    case kGrow:    case kFade:    case kNoGrow:
      if (m_CurrStep == m_StartStep) {
        RetVal = kTrue;
      }
      break;
    // for mcti blocks, clear accumulator after each equilibration
    case kStepUp:  case kStepDown:  case kStepGrow:  case kStepFade:
      if ( (m_CurrStep-(m_StartStep+m_NumEquilSteps)) % 
             (m_NumEquilSteps+m_NumAccumSteps) == 1) {
        RetVal = kTrue;
      }
      break;
    // should never get here
    default:
      ASSERT(kFalse);
      break;
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsFirstStep() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's the
// first step of the control object.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  if (m_CurrStep == (m_StartStep+1)) {
    RetVal = kTrue;
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsTimeToPrint() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to print out restraint information.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if printing is required
  if (m_NumPrintSteps > 0) {
    // if number-of-steps from StartStep is an even multiple of NumPrintSteps
    // or if it's the last step of this pmf or mcti block,
    // then it's time to print
    if ( IsLastStep() || (((m_CurrStep-m_StartStep)%m_NumPrintSteps)==0) ) {
      RetVal = kTrue;
    }
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsEndOf_MCTI() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if this is
// the last time step of an mcti block.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if an MCTI block is currently active
  if ( (m_Task==kStepUp)   || (m_Task==kStepDown) ||
       (m_Task==kStepGrow) || (m_Task==kStepFade) ) {
    // if this is the last step
    if (m_CurrStep == GetLastStep()) {
      RetVal = kTrue;
    }
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsEndOf_MCTI_Step() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if this is
// the last time step of an mcti step.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if an MCTI block is currently active
  if ( (m_Task==kStepUp)   || (m_Task==kStepDown) ||
       (m_Task==kStepGrow) || (m_Task==kStepFade) ) {
    // if this is the last step of accumulation
    if (((m_CurrStep-m_StartStep)%(m_NumEquilSteps+m_NumAccumSteps))==0) {
      // then this is the last time step of this mcti step
      RetVal = kTrue;
    }
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsTimeToPrint_dU_dLambda() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to print out dU/dLambda information.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if printing is required
  if (m_NumPrintSteps > 0) {
    // if number-of-steps from StartStep is an even multiple of NumPrintSteps
    // or if it's the last step of this pmf or mcti block,
    // then it might be time to print du/dLambda
    if ( IsLastStep() || (((m_CurrStep-m_StartStep)%m_NumPrintSteps)==0) ) {
      // for mcti blocks
      if ((m_Task==kStepUp)   || (m_Task==kStepDown) ||
          (m_Task==kStepGrow) || (m_Task==kStepFade)) {
        // it's only time to print if we're no longer equilibrating
        if ( ((m_CurrStep-m_StartStep-1) % (m_NumEquilSteps+m_NumAccumSteps))
             > m_NumEquilSteps) {
          RetVal = kTrue;
        }
      }
      // for other blocks (up, down, grow, fade) it's time
      else {
        RetVal = kTrue;
      }
    }
  }
  return(RetVal);
}


void ALambdaControl::PrintLambdaHeader(double dT) {
//----------------------------------------------------------------------------
// ASSUMING that this object is currently active, write out the header for
// a new lambda control object
//----------------------------------------------------------------------------
  // double Time;
  
  // calculate current time in femto-seconds
  // Time = (double)m_CurrStep * dT;

iout << "FreeEnergy: ";
#if !defined(_VERBOSE_PMF)
  iout << "nstep  time(ps)  ";
  iout << "    task  lambdaKf  lambdaRef     delta-G  #steps  n*{value  target |}" << std::endl;
  iout << "FreeEnergy: -----  --------  ";
  iout << "--------  --------  ---------  ----------  ------  -------------------" << std::endl;
  iout << endi;
#endif
}


void ALambdaControl::PrintHeader(double dT) {
//----------------------------------------------------------------------------
// ASSUMING that this object is currently active, write out the current time
//----------------------------------------------------------------------------
  double  Time;
  char    Str[100], Str2[100];

  // calculate current time in femto-seconds
  Time = (double)m_CurrStep * dT;

#if defined(_VERBOSE_PMF)
  iout << "FreeEnergy: " << std::endl << endi;
  iout << "FreeEnergy: ";
  iout << "Time Step = "  << m_CurrStep            <<    ",  ";
  iout << "Time = ";
  // write out time in ps
  iout << Time/1000.0     << " ps,  ";
  iout << "Lambda_Kf = "  << m_LambdaKf            <<    ",  ";
  iout << "Lambda_Ref = " << m_LambdaRef           <<     "  ";
  GetTaskStr(Str);
  iout << "(" << Str << ")";
  iout << std::endl << endi;
  iout << "FreeEnergy: ";
  iout << "------------------------------------------------";
  iout << "-------------------";
  iout << std::endl << endi;
#else
  sprintf(Str, "%5d", m_CurrStep);
  // write out time in ps
  sprintf(Str2, "%8.3f", Time/1000.0);
  iout << "FreeEnergy: ";
  iout << Str << "  " << Str2 << "  ";
  GetPaddedTaskStr(Str);
  iout << Str << "  ";
  sprintf(Str, "%8.5f", m_LambdaKf);
  iout << Str << "  ";
  sprintf(Str, "%9.5f", m_LambdaRef);
  iout << Str << "  ";
#endif
}


Bool_t ALambdaControl::IsActive() {
//------------------------------------------------------------------------
// determine if this object is currently active
//------------------------------------------------------------------------
  if ( (m_CurrStep>=m_StartStep) && (m_CurrStep<=GetLastStep()) ) {
    return(kTrue);
  }
  else {
    return(kFalse);
  }
}


double ALambdaControl::GetLambdaKf() {
//------------------------------------------------------------------------
// calculate LambdaKf for grow, fade, stepgrow, and stepfade.
// LambdaKf=1.0, for up, down, stepup, stepdown, and stop.
// for nogrow, LambdaKf is Lambda from the config file.
//------------------------------------------------------------------------
  int     N;
  double  CurrStep  = m_CurrStep;
  double  StartStep = m_StartStep;
  double  NumSteps  = m_NumSteps;
  double  NumEquilSteps = m_NumEquilSteps;
  double  NumAccumSteps = m_NumAccumSteps;
  double  NumRepeats = m_NumRepeats;

  if (IsActive()) {
    switch (m_Task) {
      case kGrow:
        m_LambdaKf = (CurrStep-StartStep)/NumSteps;
        break;
      case kFade:
        m_LambdaKf = 1.0-(CurrStep-StartStep)/NumSteps;
        break;
      case kStepGrow:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaKf = (N+1)/NumRepeats;
        break;
      case kStepFade:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaKf = 1.0 - (N+1)/NumRepeats;
        break;
      case kNoGrow:
        break;              // return prior setting of m_LambdaKf
      default:
        m_LambdaKf=1.0;
    }
  }
  else {
    m_LambdaKf=1.0;
  }
  return(m_LambdaKf);
}


double ALambdaControl::GetLambdaRef() {
//------------------------------------------------------------------------
// calculate LambdaRef for up, down, stepup, and stepdown.
// for stop, LambdaRef is Lambda from the config file.
// for grow, fade, stepgrow, stepfade, and nogrow,
//   LambdaRef is LambdaT from the config file.
//------------------------------------------------------------------------
  int     N;
  double  CurrStep  = m_CurrStep;
  double  StartStep = m_StartStep;
  double  NumSteps  = m_NumSteps;
  double  NumEquilSteps = m_NumEquilSteps;
  double  NumAccumSteps = m_NumAccumSteps;
  double  NumRepeats = m_NumRepeats;

  if (IsActive()) {
    switch (m_Task) {
      case kUp:
        m_LambdaRef = (CurrStep-StartStep)/NumSteps;
        break;
      case kDown:
        m_LambdaRef = 1.0-(CurrStep-StartStep)/NumSteps;
        break;
      case kStepUp:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaRef = (N+1)/NumRepeats;
        break;
      case kStepDown:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        m_LambdaRef = 1.0 - (N+1)/NumRepeats;
      default: 
        break;             // return prior setting of m_LambdaRef
    }
  }
  else {
    m_LambdaRef=0.0;
  }
  return(m_LambdaRef);
}

