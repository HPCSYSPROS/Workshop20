/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//-----------------------------------------------------------------------------
// ALambdaManager contains a (potentially long) list of LambdaControl objects.
// written by David Hurwitz, March to May 1998.
//-----------------------------------------------------------------------------
#if !defined(LAMBDMGR_HPP)
  #define LAMBDMGR_HPP

const int kLambdaNumToStart = 16; // to start, there's room for this num objects.
const int kLambdaMultiplier = 4;  // each time array size is exceeded,
                                  // its size is increased by this many times.

class ALambdaManager {
private:
  ALambdaControl*  m_pPmfBlocks;  // the list of objects
  int  m_NumObjects;              // the number of objects in the list
  int  m_MaxNum;                  // the maximum number of objects allowed in
                                  // the list without allocating more memory
  int  m_ActiveIndex;             // keep track of which LambdaControl is active
  ALambdaControl m_Dummy;

public:
  ALambdaManager();
  ~ALambdaManager();
  ALambdaControl&  operator[] (int index);
  void    Clear();
  int     Add(ALambdaControl& PmfBlock);
  int     GetNumObjects() {return(m_NumObjects);}
  Bool_t  GetLambdas(double& LambdaKf, double& LambdaRef);
  Bool_t  IsTimeToPrint();
  Bool_t  IsFirstStep();
  Bool_t  IsTimeToPrint_dU_dLambda();
  Bool_t  IsTimeToClearAccumulator();
  Bool_t  IsEndOf_MCTI_Step();
  Bool_t  IsEndOf_MCTI();
  int     GetNumStepsSoFar();
  int     GetNumAccumStepsSoFar();
  void    PrintHeader(double dT);
  void    PrintLambdaHeader(double dT);
  void    Print_dU_dLambda_Summary(double Sum_dU_dLambdas);
  void    PrintSomeSpaces();
  void    Print_MCTI_Integration();
  void    IncCurrStep() {m_Dummy.IncCurrStep();}
  int     GetTotalNumSteps();
  void    Integrate_MCTI();
  void    Accumulate(double dU_dLambda);
  double  GetAccumulation();
  double  GetIntegration();
  void    ZeroAccumulator();
  int     GetNum_dU_dLambda();
};

#endif

