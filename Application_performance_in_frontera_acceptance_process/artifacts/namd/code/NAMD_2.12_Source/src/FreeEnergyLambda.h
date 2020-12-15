/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#if !defined(LAMBDA_HPP)
  #define LAMBDA_HPP

class ALambdaControl {
private:

  // don't forget to change operator= if member variables change
  int     m_NumSteps;         // for pmf block
  int     m_NumEquilSteps;    // for mcti block
  int     m_NumAccumSteps;    // "
  int     m_NumRepeats;       // "
  int     m_NumPrintSteps;    // for pmf & mcti blocks
  int     m_StartStep;        // "
  int     m_StopStep;         // "
  double  m_LambdaKf;         // "
  double  m_LambdaRef;        // "
  feptask_t  m_Task;             // "
  double  m_Sum_dU_dLambda;   // for accumulating dU/dLambda
  int     m_Num_dU_dLambda;   // number averaged
  double  m_MCTI_Integration; // for accumulating <dU/dLambda> * dLambda

  static int  m_CurrStep;     // for all pmf & mcti blocks

public:
  ALambdaControl();
  void    Init(ALambdaControl& PriorBlock);
  double  GetLambdaKf();
  double  GetLambdaRef();
  Bool_t  IsActive();
  Bool_t  IsTimeToPrint();
  Bool_t  IsFirstStep();
  Bool_t  IsTimeToPrint_dU_dLambda();
  Bool_t  IsTimeToClearAccumulator();
  Bool_t  IsEndOf_MCTI_Step();
  Bool_t  IsEndOf_MCTI();
  void    PrintHeader(double dT);
  void    PrintLambdaHeader(double dT);
  void    IncCurrStep() {m_CurrStep++;}
  ALambdaControl&  operator= (ALambdaControl& PmfBlock);
  void    GetTaskStr(char* Str);
  void    GetPaddedTaskStr(char* Str);
  void    Integrate_MCTI();
  void    Accumulate(double dU_dLambda);
  double  GetIntegration();
  double  GetAccumulation();
  void    ZeroAccumulator() {
    m_Sum_dU_dLambda = 0.0;
    m_Num_dU_dLambda = 0;
  }

  int    GetNumSteps();
  int    GetNumStepsSoFar()            {return(m_CurrStep-m_StartStep);}
  int    GetNumAccumStepsSoFar();
  int    GetNum_dU_dLambda()           {return(m_Num_dU_dLambda);}
  void   SetNumSteps(int Steps)        {m_NumSteps=Steps;}
  void   SetNumEquilSteps(int Steps)   {m_NumEquilSteps=Steps;}
  void   SetNumAccumSteps(int Steps)   {m_NumAccumSteps=Steps;}
  void   SetNumPrintSteps(int Steps)   {m_NumPrintSteps=Steps;}
  void   SetNumRepeats(int Repeats)    {m_NumRepeats=Repeats;}
  void   SetStartStep(int Step)        {m_StartStep=Step;}
  void   SetStopStep(int Step)         {m_StopStep=Step;}
  void   SetLambdaKf(double LambdaKf)  {m_LambdaKf=LambdaKf;}
  void   SetLambdaRef(double LambdaRef){m_LambdaRef=LambdaRef;}
  void   SetTask(feptask_t Task)          {m_Task=Task;}
  feptask_t GetTask()                     {return(m_Task);}

private:
  Bool_t IsLastStep();
  int    GetLastStep();
};

#endif

