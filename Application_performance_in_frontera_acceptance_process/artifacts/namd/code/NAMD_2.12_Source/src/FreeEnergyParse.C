/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "common.h"
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

#include "FreeEnergyParse.h"


void ProblemParsing(const char* Message, const char* Str, Bool_t Terminate) {
//----------------------------------------------------------------------------
// print this message if there's a problem parsing
//----------------------------------------------------------------------------
  iout << "FreeEnergy: " << std::endl << endi;
  iout << "FreeEnergy: ";
  iout << "Problem parsing input parameters" << std::endl << endi;
  iout << "FreeEnergy: ";
  if (Terminate) {
    iout << "  Error:       " << Message << std::endl << endi;
  }
  else {
    iout << "  Warning:     " << Message << std::endl << endi;
  }
  iout << "FreeEnergy: ";
  iout << "  Read Until:  " << Str << std::endl << endi;
  iout << "FreeEnergy: " << std::endl << endi;
  if (Terminate) {
    NAMD_die("FreeEnergy: Fatal Parsing Error");
  }
}


void CheckParentheses(const char* Str) {
//----------------------------------------------------------------------------
// check for balanced '(' ')' and '{' '}'
//----------------------------------------------------------------------------
  int  ParenthesesCount = 0;
  int  CurlyBracketCount = 0;

  for (unsigned int i=0; i<strlen(Str); i++) {
    if (Str[i] == '(') {ParenthesesCount++;}
    if (Str[i] == ')') {ParenthesesCount--;}
    if (Str[i] == '{') {CurlyBracketCount++;}
    if (Str[i] == '}') {CurlyBracketCount--;}
    if ((ParenthesesCount<0) || (CurlyBracketCount<0)) {
      ProblemParsing("Mismatched Parentheses", Str+i);
    }
  }
}


void ReadInput(char* Str, ARestraintManager& RMgr,
                          ALambdaManager& LMgr,
                          GlobalMasterFreeEnergy& CFE,
                          double dT) {
//----------------------------------------------------------------------------
// parse the input string.  Add restraints to RMgr.  Add PmfBlocks to LMgr.
//----------------------------------------------------------------------------
  int             Count;
  char*           OldStr=NULL;   //make sure it's not equal Str to start
  ALambdaControl  PmfBlock, OldPmfBlock;
  // Bool_t          Terminate;

  // read from Str until can't read anymore
  ToLower(Str);
  CheckParentheses(Str);
  Str += ReadWhite(Str);
  while (OldStr != Str) {
    OldStr = Str;
    // add restraints to restraint manager
    Str += ReadRestraints(Str, RMgr, CFE);
    // read a single PmfBlock
    Count = ReadPmfBlock(Str, PmfBlock, dT);
    if (Count) {
      Str += Count;
      // add it to the Lambda manger
      LMgr.Add(PmfBlock);
      // initialize the default parameters of the next PmfBlock
      OldPmfBlock = PmfBlock;
      PmfBlock.Init(OldPmfBlock);
    }
  }
  Str += ReadWhite(Str);
  if (strlen(Str) > 0) {
    ProblemParsing("Unable to Read Entire Input File", Str, /* Terminate= */ kFalse);
  }
}


int ReadPmfBlock(char* Str, ALambdaControl& PmfBlock, double dT) {
//----------------------------------------------------------------------------
// if str starts with "pmf" or "mcti", read all the specs for this block
// and initialize PmfBlock.
//
// PmfBlock will already have default values for all parameters.
// The parameters that are read will be used to override the defaults.
// LambdaKf and LambdaRef are initialized after all parameters have been
// read, because the meaning of these parameters depends upon the task.
//
// If time-units (fs, ps, ns) are not specified, use ps
//
// return the number of chars to read past this block.
// return 0 if Str does not start with "pmf" or "mcti"
//----------------------------------------------------------------------------
  int     Count, Count1, Count2;
  Bool_t  Finished=kFalse;
  char*   FullString=Str;
  char*   TempStr;
  pmf_t   PmfSpec;
  // if time-units are not specified, then they're in ps.
  TimeUnits_t  TimeUnits;
  TimeUnits_t  DefaultTimeUnits = k_ps;
  // illegal default value.  user will have to specify this.
  feptask_t  Task=kUnknownTask;
  double  Lambda=-1, LambdaT=-1, Time;
  double  Dummy;
  int     NumRepeats=-1;

  const Bool_t  kPrintErrMsg=kTrue; // kNoErrMsg=kFalse;

  // if Str does not begin with "pmf" or "mcti" return 0
  Count1 = ReadWord(Str, "pmf");
  Count2 = ReadWord(Str, "mcti");
  Count = Count1 + Count2;
  if (Count==0) {
    return(0);
  }

  // skip past "pmf" or "mcti"
  Str += Count;
  // skip past "{"
  Str += ReadChar(Str, '{');
  // read spec's until "}" is found or can't read a spec
  do {
    Count = ReadNextPmfSpec(Str, PmfSpec);
    Str += Count;
    if (Count==0) {
      ProblemParsing("Unable to Read PMF Specification", Str);
    }
    // skip past "="
    Str += ReadChar(Str, '=');
    switch (PmfSpec) {
      case kTask:
        TempStr = Str;
        Str += ReadTaskType(Str, Task);
        if (Str == TempStr) {
          ProblemParsing("Can't Read Task", Str);
        }
        PmfBlock.SetTask(Task);
        break;
      case kTime:
        Str += ReadAValue(Str, Time, kPrintErrMsg);
        Str += ReadTimeUnits(Str, TimeUnits, DefaultTimeUnits);
        Time = GetTime(Time, TimeUnits);
        PmfBlock.SetNumSteps((int)(Time/dT));
        break;
      case kLambda:
        Str += ReadAValue(Str, Lambda, kPrintErrMsg);
        break;
      case kLambdaT:
        Str += ReadAValue(Str, LambdaT, kPrintErrMsg);
        break;
      case kPrint:
        Str += ReadAValue(Str, Time, kPrintErrMsg);
        Str += ReadTimeUnits(Str, TimeUnits, DefaultTimeUnits);
        Time = GetTime(Time, TimeUnits);
        PmfBlock.SetNumPrintSteps((int)(Time/dT));
        break;
      case kNoPrint:
        PmfBlock.SetNumPrintSteps(-1);
        break;
      case kEquilTime:
        Str += ReadAValue(Str, Time, kPrintErrMsg);
        Str += ReadTimeUnits(Str, TimeUnits, DefaultTimeUnits);
        Time = GetTime(Time, TimeUnits);
        PmfBlock.SetNumEquilSteps((int)(Time/dT));
        break;
      case kAccumTime:
        Str += ReadAValue(Str, Time, kPrintErrMsg);
        Str += ReadTimeUnits(Str, TimeUnits, DefaultTimeUnits);
        Time = GetTime(Time, TimeUnits);
        PmfBlock.SetNumAccumSteps((int)(Time/dT));
        break;
      case kNumRepeats:
        Str += ReadAValue(Str, Dummy, kPrintErrMsg);
        NumRepeats = (int)Dummy;
        PmfBlock.SetNumRepeats(NumRepeats);
        break;
      default:
        ASSERT(kFalse);
    }
    Count = ReadChar(Str, '}');
    Str += Count;
    // if "}" was read, then we're finished
    if (Count!=0) {Finished=kTrue;}
  }
  while(!Finished);

  // if Task was not specified above, then use PmfBlock's default task
  if (Task==kUnknownTask) {
    Task = PmfBlock.GetTask();
  }
  // if Task wasn't specified earlier, die.
  if (Task==kUnknownTask) {
    ProblemParsing("Must Specify Task", FullString);
  }

  // set LambdaKf and LambdaRef, using Lambda and LambdaT
  // (the former is the internal representation, the latter is the user's)
  // there's no need to set them if they haven't been read in this routine
  switch(Task) {
    case kStop:
      PmfBlock.SetLambdaKf(1.0);
      if (Lambda >= -kALittle) {
        PmfBlock.SetLambdaRef(Lambda);
      }
      break;
    case kNoGrow:
      if (Lambda >= -kALittle) {
        PmfBlock.SetLambdaKf(Lambda);
      }
    case kGrow:
    case kFade:
    case kStepGrow:
    case kStepFade:
      if (LambdaT >= -kALittle) {
        PmfBlock.SetLambdaRef(LambdaT);
      }
      break;
    default:
      ASSERT((Task==kUp)||(Task==kDown)||(Task==kStepUp)||(Task==kStepDown));
      PmfBlock.SetLambdaKf(1.0);
      break;
  }
  
  return(Str-FullString);
}


double GetTime(double Val, TimeUnits_t Units) {
//----------------------------------------------------------------------------
// convert (Val Units) to fs, where Units is either fs, ps, or ns
//----------------------------------------------------------------------------
  switch (Units) {
    case k_fs:  return(Val);
    case k_ps:  return(Val*1000);
    case k_ns:  return(Val*1000000);
    default:
      ASSERT(kFalse);
      return(Val);
  }
}


int ReadTimeUnits(char* Str, TimeUnits_t& Units, TimeUnits_t DefaultUnits) {
//----------------------------------------------------------------------------
// Str should start with one of:
//   "fs", "ps", or "ns"
//
// If Str starts with one of the above, return an identifier for the above.
// if Str does not start with one of the above, set Units to DefaultUnits.
//
// return the number of chars to read past the time units.  return 0 if
// Str does not start with one of the above.
//----------------------------------------------------------------------------
  char*   FullString=Str;
  Bool_t  FoundIt=kFalse;

  if (strncmp(Str,"fs",2)==0) {Units=k_fs; FoundIt=kTrue;}
  if (strncmp(Str,"ps",2)==0) {Units=k_ps; FoundIt=kTrue;}
  if (strncmp(Str,"ns",2)==0) {Units=k_ns; FoundIt=kTrue;}

  if (FoundIt) {
    Str += ReadAlpha(Str);
    Str += ReadWhite(Str);
  }
  else {
    Units = DefaultUnits;
  }
  return(Str-FullString);
}


int ReadTaskType(char* Str, feptask_t& Task) {
//----------------------------------------------------------------------------
// Str should start with a task, one of:
//   "up", "down", "stop", "grow", "fade", "nogrow",
//   "stepup", "stepdown", "stepgrow", "stepfade"
//
// return an identifer for the above.
// return the number of chars to read past the task.  return 0 if Str does
// not start with one of the above.
//----------------------------------------------------------------------------
  char*   FullString=Str;

  Task = kUnknownTask;
  if (strncmp(Str,"up",2)==0)       {Task=kUp;       goto GotIt;}
  if (strncmp(Str,"down",4)==0)     {Task=kDown;     goto GotIt;}
  if (strncmp(Str,"stop",4)==0)     {Task=kStop;     goto GotIt;}
  if (strncmp(Str,"grow",4)==0)     {Task=kGrow;     goto GotIt;}
  if (strncmp(Str,"fade",4)==0)     {Task=kFade;     goto GotIt;}
  if (strncmp(Str,"nogrow",6)==0)   {Task=kNoGrow;   goto GotIt;}
  if (strncmp(Str,"stepup",6)==0)   {Task=kStepUp;   goto GotIt;}
  if (strncmp(Str,"stepdown",8)==0) {Task=kStepDown; goto GotIt;}
  if (strncmp(Str,"stepgrow",8)==0) {Task=kStepGrow; goto GotIt;}
  if (strncmp(Str,"stepfade",8)==0) {Task=kStepFade; goto GotIt;}
  return(0);

GotIt:
  Str += ReadAlpha(Str);
  Str += ReadWhite(Str);
  return(Str-FullString);
}


int ReadNextPmfSpec(char* Str, pmf_t& PmfSpec) {
//----------------------------------------------------------------------------
// Str should start with the next spec for a pmf or mcti block, one of:
//   "task", "time", "lambda", "lambdaT", "print",
//   "equiltime", "accumtime", "numsteps"
//
// Return an identifier of the above.
// Return number of characters, including trailing white space to read past
// this word.  Return NumChars=0 if Str does not begin with one of the above.
//----------------------------------------------------------------------------
  char*  FullString=Str;

  PmfSpec = kUnknownPmf;
  if (strncmp(Str,"task",4)==0)    {PmfSpec=kTask;       goto GotIt;}
  if (strncmp(Str,"time",4)==0)    {PmfSpec=kTime;       goto GotIt;}
  // check for lambdat first, else PmfSpec will be kLambda for "lambdat"
  if (strncmp(Str,"lambdat",7)==0) {PmfSpec=kLambdaT;    goto GotIt;}
  if (strncmp(Str,"lambda",6)==0)  {PmfSpec=kLambda;     goto GotIt;}
  if (strncmp(Str,"print",5)==0)   {PmfSpec=kPrint;      goto GotIt;}
  if (strncmp(Str,"nopr",4)==0)    {PmfSpec=kNoPrint;    goto GotIt;}
  if (strncmp(Str,"equil",5)==0)   {PmfSpec=kEquilTime;  goto GotIt;}
  if (strncmp(Str,"accum",5)==0)   {PmfSpec=kAccumTime;  goto GotIt;}
  if (strncmp(Str,"numstep",7)==0) {PmfSpec=kNumRepeats; goto GotIt;}
  return(0);

GotIt:
  Str += ReadAlpha(Str);
  Str += ReadWhite(Str);
  return(Str-FullString);
}


int ReadRestraints(char* Str, ARestraintManager& AllRestraints,
                              GlobalMasterFreeEnergy& CFE) {
//----------------------------------------------------------------------------
// if Str starts with "urestr", read each restraint spec, create
// a restraint, and put a pointer to the restraint into AllRestraints
//
// return the number of chars to read past all restraints
// return 0 if Str does not start with "urestr"
//----------------------------------------------------------------------------
  int          Count;
  Bool_t       Finished=kFalse;
  char*        FullString=Str;
  ARestraint*  pRestraint;

  Count = ReadWord(Str, "urest");
  // when "urestraint" is found
  if (Count) {
    // skip past "urest..."
    Str += Count;
    // skip past "{"
    Str += ReadChar(Str, '{');
    // read restraints until "}" is found or can't read a restraint
    do {
      pRestraint = GetRestraint(Str, Count, CFE);
      Str += Count;
      // if a restraint could not be read, then we're finished
      if (Count==0) {Finished=kTrue;}
      if (!Finished) {
        AllRestraints.Add(pRestraint);
        Count = ReadChar(Str, '}');
        Str += Count;
        // if "}" was read, then we're finished
        if (Count!=0) {Finished=kTrue;}
      }
    }
    while(!Finished);
  }
  return(Str-FullString);
}


ARestraint* GetRestraint(char* Str, int& NumChars, GlobalMasterFreeEnergy& CFE) {
//----------------------------------------------------------------------------
// read spec's for a restraint, from the input string.
// allocate space for and initialize a restraint object with these specs.
// note:  memory is allocated here, and free'd elsewhere.
//
// return a pointer to this object.
// return the number of characters to read past these specs.
// return NumChars=0 for illegal restraint specs.
//----------------------------------------------------------------------------
  AGroup      Group1, Group2, Group3, Group4;
  ARestraint* pRestraint = NULL;
  restr_t     Restraint;
  Bound_t     Bound;
  int         Count;
  double      Kf;
  double      D,   D0,   D1;
  double      A=0, A0=0, A1=0, A2=0;
  AVector     Pos, Pos0, Pos1;
  char*       FullStr;
  char*       TempStr;

  const Bool_t  kPrintErrMsg=kTrue; // kNoErrMsg=kFalse;

  // save pointer to full string
  FullStr = Str;
  NumChars = 0;

  // get restraint type
  Restraint = ReadNextRestraintType(Str, Count);
  if (Count == 0) {
    ProblemParsing("Can't Read Restraint Type", Str);
    return(pRestraint);
  }

  // skip past restraint type
  ASSERT(Restraint != kUnknownRestr);
  Str += Count;
  
  // read in appropriate number of atoms or groups-of-atoms for
  // this restraint type, put the atoms in Group1 thru Group4
  switch (Restraint) {
    case kDihe:  case kDiheBound:  case kDihePMF:
      Str += AddAtoms(Group4, Str, CFE);
    case kAngle: case kAngleBound: case kAnglePMF:
      Str += AddAtoms(Group3, Str, CFE);
    case kDist:  case kDistBound:  case kDistPMF:
      Str += AddAtoms(Group2, Str, CFE);
    case kPosi:  case kPosiBound:  case kPosiPMF:
      Str += AddAtoms(Group1, Str, CFE);
    default: ;
  }

  // for dihedrals, allow keywords of "barr=", "gap=", OR "kf="
  // for other restraints, just allow "kf="
  TempStr = Str;
  switch(Restraint) {
    case kDihe:
    case kDiheBound:
    case kDihePMF:
      Str += ReadWord(Str, "barr");
      Str += ReadWord(Str, "gap");
    default:
      Str += ReadWord(Str, "kf");
      // make sure the word "barr", "gap", or "kf" was read
      if (Str==TempStr) {
        ProblemParsing("Word Missing: barr, gap, or kf", Str);
      }
      Str += ReadChar(Str, '=');
  }
  // get the Kf value
  Str += ReadAValue(Str, Kf,          kPrintErrMsg);
  
  // read the reference positions, distances or angles
  switch (Restraint) {
    case kPosi:  
      Str += ReadWord(Str, "ref",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadChar(Str, '(');
      Str += ReadAValue(Str, Pos[0],  kPrintErrMsg);
      Str += ReadAValue(Str, Pos[1],  kPrintErrMsg);
      Str += ReadAValue(Str, Pos[2],  kPrintErrMsg);
      Str += ReadChar(Str, ')');
      break;
    case kDist:
      Str += ReadWord(Str, "ref",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, D,       kPrintErrMsg);
      break;
    case kAngle:
      Str += ReadWord(Str, "ref",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A,       kPrintErrMsg);
      break;
    case kDihe:
      Str += ReadWord(Str, "ref",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A,       kPrintErrMsg);
      break;
    case kPosiBound:
      if (ReadBound(Str, Bound) == 0) {ProblemParsing("Missing Bound", Str);}
      Str += ReadWord(Str, "low");
      Str += ReadWord(Str, "hi");
      Str += ReadChar(Str, '=');
      Str += ReadChar(Str, '(');
      Str += ReadAValue(Str, Pos[0],  kPrintErrMsg);
      Str += ReadAValue(Str, Pos[1],  kPrintErrMsg);
      Str += ReadAValue(Str, Pos[2],  kPrintErrMsg);
      Str += ReadAValue(Str, D,       kPrintErrMsg);
      Str += ReadChar(Str, ')');
      break;
    case kDistBound:
      if (ReadBound(Str, Bound) == 0) {ProblemParsing("Missing Bound", Str);}
      Str += ReadWord(Str, "low");
      Str += ReadWord(Str, "hi");
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, D,       kPrintErrMsg);
      break;
    case kAngleBound:
      if (ReadBound(Str, Bound) == 0) {ProblemParsing("Missing Bound", Str);}
      Str += ReadWord(Str, "low");
      Str += ReadWord(Str, "hi");
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A,       kPrintErrMsg);
      break;
    case kDiheBound:
      Str += ReadWord(Str, "low",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A0,      kPrintErrMsg);
      Str += ReadWord(Str, "hi",     kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A1,      kPrintErrMsg);
      Str += ReadWord(Str, "delta",  kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A2,      kPrintErrMsg);
      break;
    case kPosiPMF:
      Str += ReadWord(Str, "low",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadChar(Str, '(');
      Str += ReadAValue(Str, Pos0[0], kPrintErrMsg);
      Str += ReadAValue(Str, Pos0[1], kPrintErrMsg);
      Str += ReadAValue(Str, Pos0[2], kPrintErrMsg);
      Str += ReadChar(Str, ')');
      Str += ReadWord(Str, "hi",     kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadChar(Str, '(');
      Str += ReadAValue(Str, Pos1[0], kPrintErrMsg);
      Str += ReadAValue(Str, Pos1[1], kPrintErrMsg);
      Str += ReadAValue(Str, Pos1[2], kPrintErrMsg);
      Str += ReadChar(Str, ')');
      break;
    case kDistPMF:
      Str += ReadWord(Str, "low",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, D0,      kPrintErrMsg);
      Str += ReadWord(Str, "hi",     kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, D1,      kPrintErrMsg);
      break;
    case kAnglePMF:
      Str += ReadWord(Str, "low",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A0,      kPrintErrMsg);
      Str += ReadWord(Str, "hi",     kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A1,      kPrintErrMsg);
      break;
    case kDihePMF:
      Str += ReadWord(Str, "low",    kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A0,      kPrintErrMsg);
      Str += ReadWord(Str, "hi",     kPrintErrMsg);
      Str += ReadChar(Str, '=');
      Str += ReadAValue(Str, A1,      kPrintErrMsg);
      break;
    default: ;
  }

  // convert degrees to radians
  A  *= (kPi/180);
  A0 *= (kPi/180);
  A1 *= (kPi/180);
  A2 *= (kPi/180);

  // initialize the restraint
  switch (Restraint) {
    case kPosi:
      pRestraint = new AFixedPosRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group1);
      pRestraint->SetRefPos(Pos);
      break;
    case kDist:
      pRestraint = new AFixedDistRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group2, Group1);
      pRestraint->SetRefDist(D);
      break;
    case kAngle:
      pRestraint = new AFixedAngleRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group3, Group2, Group1);
      pRestraint->SetRefAngle(A);
      break;
    case kDihe:
      pRestraint = new AFixedDiheRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group4, Group3, Group2, Group1);
      pRestraint->SetRefAngle(A);
      break;
    case kPosiBound:
      pRestraint = new ABoundPosRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group1);
      pRestraint->SetRefPos(Pos);
      pRestraint->SetRefDist(D);
      pRestraint->SetBound(Bound);
      break;
    case kDistBound:
      pRestraint = new ABoundDistRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group2, Group1);
      pRestraint->SetRefDist(D);
      pRestraint->SetBound(Bound);
      break;
    case kAngleBound:
      pRestraint = new ABoundAngleRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group3, Group2, Group1);
      pRestraint->SetRefAngle(A);
      pRestraint->SetBound(Bound);
      break;
    case kDiheBound:
      pRestraint = new ABoundDiheRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group4, Group3, Group2, Group1);
      pRestraint->SetLowerAngle(A0);
      pRestraint->SetUpperAngle(A1);
      pRestraint->SetIntervalAngle(A2);
      break;
    case kPosiPMF:
      pRestraint = new AForcingPosRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group1);
      pRestraint->SetStartPos(Pos0);
      pRestraint->SetStopPos(Pos1);
      break;
    case kDistPMF:
      pRestraint = new AForcingDistRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group2, Group1);
      pRestraint->SetStartDist(D0);
      pRestraint->SetStopDist(D1);
      break;
    case kAnglePMF:
      pRestraint = new AForcingAngleRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group3, Group2, Group1);
      pRestraint->SetStartAngle(A0);
      pRestraint->SetStopAngle(A1);
      break;
    case kDihePMF:
      pRestraint = new AForcingDiheRestraint;
      pRestraint->SetKf(Kf);
      pRestraint->SetGroups(Group4, Group3, Group2, Group1);
      pRestraint->SetStartAngle(A0);
      pRestraint->SetStopAngle(A1);
      break;
    default: ;
  }
  // calc number of chars to read restraint specs
  NumChars = Str-FullStr;
  return(pRestraint);
}


int ReadBound(char* Str, Bound_t& Bound) {
//----------------------------------------------------------------------------
// Str should start with "low" or "hi".  determine which it is and
// count the number of characters to read past this word + white-space.
// return NumChars=0 if Str does not start with "low" or "hi"
//----------------------------------------------------------------------------
  int  Count;

  Bound = kUnknownBound;
  Count = ReadWord(Str, "low");
  if (Count) {
    Bound=kLower;
    return(Count);
  }

  Count = ReadWord(Str, "hi");
  if (Count) {
    Bound=kUpper;
    return(Count);
  }

  return(Count);  // Count will be 0 if "low" or "hi" wasn't found
}


int ReadAValue(char* Str, double& Value, Bool_t ErrMsg) {
//----------------------------------------------------------------------------
// Str should start with a floating point number.  convert it to a double.
// also, return the number of chars to read past the value + white-space
// return NumChars = 0 if Str does not start with a valid fp number.
// Print an error message if ErrMsg is kTrue, and no value is read.
//----------------------------------------------------------------------------
  int    NumChars;
  char*  NewStr;

  // read f.p. number and trailing white-space
  Value = strtod(Str, &NewStr);
  if (NewStr != Str) {
    NewStr += ReadWhite(NewStr);
  }
  NumChars = NewStr - Str;

  // if no number was read, and ErrMsg is kTrue, print a message
  if ((NumChars==0) && (ErrMsg)) {
    ProblemParsing("Floating Point Number Expected", Str);
  }

  return(NumChars);
}


int ReadChar(char* Str, char Char, Bool_t ErrMsg) {
//----------------------------------------------------------------------------
// Str should start with Char plus trailing white space.
// return the number of chars, including the white space, to read past Char
// return 0 if Str does not start with Char.
// Print an error message if ErrMsg is kTrue, and Char is not read.
//----------------------------------------------------------------------------
  int    NumChars;
  char*  FullString=Str;
  char   Message[64];
  
  // initial part of Message
  strcpy(Message, "Character Missing:  ");

  // read char and trailing white-space
  if (Str[0] == Char) {
    Str += 1;
    Str += ReadWhite(Str);
  }
  NumChars = Str - FullString;

  // if Char was not read, and ErrMsg is kTrue, print a message
  if ((NumChars==0) && (ErrMsg)) {
    // add the character that's missing to Message
    Message[strlen(Message)-1] = Char;
    ProblemParsing(Message, Str);
  }

  return(NumChars);
}


int ReadWord(const char* Str, const char* Word, Bool_t ErrMsg) {
//----------------------------------------------------------------------------
// Str should start with Word plus, perhaps, some extra alphabetic
// characters, plus trailing white space.
//
// (do NOT read past numeric characters which follow Word, so that
//  this routine can be used on the following:  ReadWord(Str, "="),
//  where Str is, for example, "=123 ...")
//
// return the number of chars, including the white space, to read past Word
// return 0 if Str does not start with Word.
//
// Print an error message if ErrMsg is kTrue, and word is not read.
//----------------------------------------------------------------------------
  const char*  FullString=Str;
  int    NumChars, StrLen;
  char   Message[64];
  
  // initial part of message
  strcpy(Message, "Word Missing: ");

  StrLen = strlen(Word);
  if (strncmp(Str, Word, StrLen) == 0) {
    Str += StrLen;
    Str += ReadAlpha(Str);
    Str += ReadWhite(Str);
  }
  NumChars = Str - FullString;

  // if Word was not read, and ErrMsg is kTrue, print a message
  if ((NumChars==0) && (ErrMsg)) {
    strcat(Message, Word);
    ProblemParsing(Message, Str);
  }

  return(NumChars);
}


restr_t ReadNextRestraintType(char* Str, int& NumChars) {
//----------------------------------------------------------------------------
// Str should start with the next restraint type (no leading white space),
// namely one of:
// 
//   "posi",       "dist",       "angle",       "dihe"
//   "posi bound", "dist bound", "angle bound", "dihe bound"
//   "posi pmf",   "dist pmf",   "angle pmf",   "dihe pmf"
//
// figure out which it is.
// also, return the number of characters, including trailing white space.
// return NumChars=0 for an illegal restraint-type
//
// the words "pos*", "dist*", "angle*", "dihe*", "bound*", and "pmf*"
// are all recognized.
//----------------------------------------------------------------------------
  restr_t  RestraintType=kUnknownRestr;
  char*    FullString=Str;

  // check if Str starts with "pos", "dist", "angl", or "dihe"
  if (strncmp(Str,"pos", 3)==0)  {RestraintType=kPosi;  goto GotIt;}
  if (strncmp(Str,"dist",4)==0)  {RestraintType=kDist;  goto GotIt;}
  if (strncmp(Str,"angl",4)==0)  {RestraintType=kAngle; goto GotIt;}
  if (strncmp(Str,"dihe",4)==0)  {RestraintType=kDihe;  goto GotIt;}
  NumChars = 0;
  return(RestraintType);

  // skip to the end of the white space following this word
GotIt:
  Str += 3;
  Str += ReadAlphaNum(Str);
  Str += ReadWhite(Str);

  // check if the next word is "bound", skip to the end of this word
  if (strncmp(Str,"bound",5)==0) {
    switch (RestraintType) {
      case kPosi:   RestraintType=kPosiBound;   break;
      case kDist:   RestraintType=kDistBound;   break;
      case kAngle:  RestraintType=kAngleBound;  break;
      case kDihe:   RestraintType=kDiheBound;   break;
      default: break;
    }
    Str += 5;
    Str += ReadAlphaNum(Str);
  }

  // check if the next word is "pmf", skip to the end of this word
  if (strncmp(Str,"pmf",3)==0) {
    switch (RestraintType) {
      case kPosi:   RestraintType=kPosiPMF;   break;
      case kDist:   RestraintType=kDistPMF;   break;
      case kAngle:  RestraintType=kAnglePMF;  break;
      case kDihe:   RestraintType=kDihePMF;   break;
      default: break;
    }
    Str += 3;
    Str += ReadAlphaNum(Str);
  }

  // skip past trailing white space, calcuate num chars string has been advanced
  Str += ReadWhite(Str);
  NumChars = Str-FullString;
  return(RestraintType);
}


int AddAtoms(AGroup& Group, char* Str, GlobalMasterFreeEnergy& CFE) {
//----------------------------------------------------------------------------
// Str contains specifications for which atoms to add to Group.
// The atoms may be:
//   a) a single atom, b) all atoms of a residue, c) a list of atoms
//   d) all atoms in a list of residues, e) all atoms in a range of residues,
//   e) one or more atomnames in a list of residues, or
//   f) one or more atomnames in a range of residues
// Add the AtomID's for these specified atoms to Group.
// return the number of characters in Str that were read.
//----------------------------------------------------------------------------
  int     NumChars;
  int     RetNumChars = 0;
  Bool_t  GroupMode = kFalse;
  Bool_t  AtomList =  kFalse;
  Bool_t  Finished =  kFalse;
  char*   SavePtr = 0;

  while (!Finished) {
    switch(ReadNextItem(Str, NumChars)) {
      case kStartGroup:
        GroupMode = kTrue;
        break;
      case kEndGroup:
        Finished = kTrue;
        break;
      case kAtom:
        AddAtom(Group, Str, CFE);
        if (!GroupMode) {
          Finished = kTrue;
        }
        break;
      case kAtomName:
      case kAtomNameList:
        AtomList = kTrue;
        SavePtr = Str;
        break;
      case kResidue:
      case kResidueRange:
        if (AtomList) {
          AddAtomsInResidues(Group, SavePtr, Str, CFE);
        }
        else {
          AddResidues(Group, Str, CFE);
        }
        if (!GroupMode) {
          Finished = kTrue;
        }
        break;
      default:
        Finished = kTrue;
        ProblemParsing("Can't Read Atoms", Str);
        break;
    }
    Str += NumChars;
    RetNumChars += NumChars;
  }
  RetNumChars += ReadWhite(Str);
  return(RetNumChars);
}


void AddAtomsInResidues(AGroup& Group, char* AtomNames, char* ResRange,
                        GlobalMasterFreeEnergy& CFE) {
//-------------------------------------------------------------------
// Group contains a list of int's representing AtomID's.
// ResRange should be "(segname, resnum) to (segname, resnum)"
//                 or "(segname, resnum)"
// AtomNames should be "(atomname, atomname, ...):" or "atomname:"
// get the atomID's for each atomname in ResRange, add them to Group.
//-------------------------------------------------------------------
  int   Count, ArrayIndex, i;
  char  AtomNamesArray[21][30];

  // skip to start of first atomname
  if (AtomNames[0] == '(') {
    AtomNames++;
    Count = ReadWhite(AtomNames);
    AtomNames += Count;
  }
  // put each atomname into the array, finish when ':' or ')' is found
  ArrayIndex = 0;
  while ( (AtomNames[0]!=':') && (AtomNames[0]!=')') ) {
    Count = ReadAlphaNum(AtomNames);
    strncpy(AtomNamesArray[ArrayIndex], AtomNames, Count);
    AtomNamesArray[ArrayIndex][Count] = '\0';
    AtomNames += Count;
    Count = ReadWhite(AtomNames);
    AtomNames += Count;
    ArrayIndex++;
  }
  // now add each atomname of Res to Group.
  // if "all" is specified, add all atoms of Res to Group.
  for (i=0; i<ArrayIndex; i++) {
    if (strcmp(AtomNamesArray[i], "all") == 0) {
      AddResidues(Group, ResRange, CFE);
    }
    else {
      AddAtom(Group, ResRange, AtomNamesArray[i], CFE);
    }
  }
}


void AddResidues(AGroup& Group, char* ResRange, GlobalMasterFreeEnergy& CFE) {
//-------------------------------------------------------------------
// Group contains a list of int's representing AtomID's.
// ResRange should be "(segname, resnum) to (segname, resnum)"
//                 or "(segname, resnum)"
// get the atomID's for each atom of ResRange, and add them to Group.
//-------------------------------------------------------------------
  char  SegName[21];
  int   ResNum1, ResNum2, ResNum;
  int   i, NumAtoms, AtomID, RetVal;

  // get start and stop residue numbers
  GetSegName(ResRange, SegName);
  GetResRange(ResRange, ResNum1, ResNum2);

  // for each residue of residue range
  for (ResNum=ResNum1; ResNum<=ResNum2; ResNum++) {
    // for each atom of residue
    NumAtoms = CFE.getNumAtoms(SegName, ResNum);
    if (NumAtoms < 1) { ProblemParsing("No Atoms in Residue", ResRange); }
    for (i=0; i<NumAtoms; i++) {
      // get atomID, register it, add it to Group
      AtomID = CFE.getAtomID(SegName, ResNum, i);
      if (AtomID < 0) { ProblemParsing("Invalid AtomID", ResRange); }
      RetVal = CFE.requestAtom(AtomID);
      if (RetVal < 0) { ProblemParsing("Unable to requestAtom", ResRange); }
      Group.Add(AtomID);
    }
  }
}


void AddAtom(AGroup& Group, char* Atom, GlobalMasterFreeEnergy& CFE) {
//-------------------------------------------------------------------
// Group contains a list of int's representing AtomID's.
// Atom should be "(segname, resnum, atomname)"
// get the atomID for Atom, and add it to Group.
//-------------------------------------------------------------------
  char  AtomName[21];

  GetAtomName(Atom, AtomName);
  AddAtom(Group, Atom, AtomName, CFE);
}


void AddAtom(AGroup& Group, char* ResRange, char* AtomName,
             GlobalMasterFreeEnergy& CFE) {
//-------------------------------------------------------------------
// Group contains a list of int's representing AtomID's.
// ResRange should be "(segname, resnum) to (segname, resnum)"
//                 or "(segname, resnum)"
// AtomName is specified separately.
// get the atomID, and add it to Group.
//-------------------------------------------------------------------
  char  SegName[21];
  int   ResNum, ResNum1, ResNum2, AtomID, RetVal;

  // convert "(segname, resnum1) to (segname, resnum2)"
  //  -> SegName, ResNum1, ResNum2
  GetSegName(ResRange, SegName);
  GetResRange(ResRange, ResNum1, ResNum2);

  // get atomID for each atom in the specified residue range
  // register it, add it to Group
  for (ResNum=ResNum1; ResNum<=ResNum2; ResNum++) {
    AtomID = CFE.getAtomID(SegName, ResNum, AtomName);
    if (AtomID < 0) { ProblemParsing("Invalid AtomID", ResRange); }
    RetVal = CFE.requestAtom(AtomID);
    if (RetVal < 0) { ProblemParsing("Unable to requestAtom", ResRange); }
    Group.Add(AtomID);
  }
}


void GetResRange(char* ResRange, int& ResNum1, int& ResNum2) {
//-------------------------------------------------------------------
// ResRange should be "(segname, resnum1) to (segname, resnum2)"
// return ResNum1 & ResNum2
// if "to" is missing, return resnum1 in both ResNum1 & ResNum2
//-------------------------------------------------------------------
  char SegName1[21], SegName2[21];

  // get start residue number
  GetSegName(ResRange, SegName1);
  GetResNum(ResRange, ResNum1);

  // skip to where "to" should appear
  ResRange += ReadParentheses(ResRange);
  ResRange += ReadWhite(ResRange);

  // if "to" is found
  if (strncmp(ResRange, "to", 2) == 0) {
    //skip to next residue
    ResRange += ReadAlphaNum(ResRange);
    ResRange += ReadWhite(ResRange);
    // get final residue number
    GetSegName(ResRange, SegName2);
    GetResNum(ResRange, ResNum2);
    // do some checks
    if (strcmp(SegName1, SegName2)!=0) {
      ProblemParsing("SegNames Differ", ResRange);
    }
    if (ResNum2 < ResNum1) {
      ProblemParsing("Decreasing Residues", ResRange);
    }
  }

  // otherwise, ResNum2 = ResNum1
  else {
    ResNum2 = ResNum1;
  }
}


int GetSegName(char* Str, char* SegName) {
//-------------------------------------------------------------------
// Str should be (segname, resnum) or (segname, resnum, atomname)
// put segname into SegName
// return the number of characters from start-of-Str thru segname
//-------------------------------------------------------------------
  int    Count;
  char*  FullString=Str;

  if (Str[0] != '(') {ProblemParsing("Missing (", Str);}
  Str += 1;
  Str += ReadWhite(Str);
  Count = ReadAlphaNum(Str);
  if (Count == 0)    {ProblemParsing("Missing Segment Name", Str);}
  strncpy(SegName, Str, Count);
  SegName[Count] = '\0';
  Str += Count;
  return(Str-FullString);
}


int  GetResNum(char* Str, int& ResNum) {
//-------------------------------------------------------------------
// Str should be (segname, resnum) or (segname, resnum, atomname)
// convert resnum to an int and return it
// return the number of characters from start-of-Str thru resnum
//-------------------------------------------------------------------
  int    Count;
  char   SegName[21];
  char*  FullString=Str;

  Str += GetSegName(Str, SegName);
  Str += ReadWhite(Str);
  ResNum = (int) strtol(Str, NULL, 10);
  Count = ReadDigits(Str);
  if (Count == 0) {ProblemParsing("Missing Residue Number", Str);}
  Str += Count;
  return(Str-FullString);
}


int GetAtomName(char* Str, char* AtomName) {
//-------------------------------------------------------------------
// Str should be (segname, resnum, atomname)
// put atomname into AtomName
// return the number of characters from start-of-Str thru ')'
//-------------------------------------------------------------------
  int    Count, ResNum;
  char*  FullString=Str;

  Str += GetResNum(Str, ResNum);
  Str += ReadWhite(Str);
  Count = ReadAlphaNum(Str);
  if (Count == 0)    {ProblemParsing("Missing Atom Name", Str);}
  strncpy(AtomName, Str, Count);
  AtomName[Count] = '\0';
  Str += Count;
  Str += ReadWhite(Str);
  if (Str[0] != ')') {ProblemParsing("Missing )", Str);}
  Str += 1;
  return(Str-FullString);
}


item_t ReadNextItem(char* Str, int& NumChars) {
//-------------------------------------------------------------------
// Figure out what the next item in Str is, and how many characters
// long it is.  The next item should be one of the following:
//   1.  kStartGroup:       Group {
//   2.  kEndGroup:         }
//   3.  kAtomName:         atomname:
//   4.  kAtomNameList:     (atomname, atomname, ... ):
//   5.  kAtom:             (segname, resno, atomname)
//   6.  kResidue:          (segname, resno)
//   7.  kResidueRange:     (segname, resno) to (segname, resno)
// The following assumptions may be made:
//   1. Str is all lowercase
//   2. There are NO leading white-char's
// Return:
//   The length of the next item, plus the white space that follows.
//-------------------------------------------------------------------
  int     Num;
  item_t  RetVal=kUnknownItem;

  Num=IsStartGroup(Str);     if (Num)  {RetVal=kStartGroup;   goto Found;}
  Num=IsEndGroup(Str);       if (Num)  {RetVal=kEndGroup;     goto Found;}
  Num=IsAtomName(Str);       if (Num)  {RetVal=kAtomName;     goto Found;}
  Num=IsAtomNameList(Str);   if (Num)  {RetVal=kAtomNameList; goto Found;}
  Num=IsAtom(Str);           if (Num)  {RetVal=kAtom;         goto Found;}
  Num=IsResidue(Str);        if (Num)  {RetVal=kResidue;      goto Found;}
  Num=IsResidueRange(Str);   if (Num)  {RetVal=kResidueRange; goto Found;}

  // add the white-space after the item to the length of the item.
Found:
  NumChars = Num;
  Str += NumChars;
  NumChars += ReadWhite(Str);
  return(RetVal);
}


int IsStartGroup(char* Str) {
//-------------------------------------------------------------------
// see if Str starts with "group {"
// return:     the number of characters, including white space.
//             0, if Str does not start with "group {"
//-------------------------------------------------------------------
  char*  FullString=Str;

  if (strncmp(Str, "group", 5) == 0) {
    Str += 5;
    Str += ReadWhite(Str);
    if (Str[0] == '{') {
      Str += 1;
      return(Str-FullString);
    }
  }
  return(0);
}


int IsEndGroup(char* Str) {
//-------------------------------------------------------------------
// see if Str starts with "}"
// return:    the number of characters, including white space.
//            0, if Str does not start with "}"
//-------------------------------------------------------------------
  if (Str[0] == '}') {
    return(1);
  }
  return(0);
}


int IsAtomName(char* Str) {
//-------------------------------------------------------------------
// see if Str starts with "atomname:"
// return:    the number of characters, including white space.
//            0, if Str does not start with "atomname:"
//-------------------------------------------------------------------
  int    Count;
  char*  FullString=Str;

  Count = ReadAlphaNum(Str);
  if (Count) {
    Str += Count;
    Str += ReadWhite(Str);
    if (Str[0] == ':') {
      Str += 1;
      return(Str-FullString);
    }
  }
  return(0);
}


int IsAtomNameList(char* Str) {
//------------------------------------------------------------------------
// see if Str starts with "(atomname, atomname, ...):"
// return:    the number of characters, including white space.
//            0, if Str does not start with "(atomname, atomname, ...):"
//------------------------------------------------------------------------
  int    Count;
  char*  FullString=Str;

  // Str will be considered an atom-name-list if it contains the following:
  // '(', anything, ')', ':'
  Count = ReadParentheses(Str);
  if (Count > 0) {
    Str += Count;
    Str += ReadWhite(Str);
    if (Str[0] == ':') {
      Str += 1;
      return(Str-FullString);
    }
  }
  return(0);
}


int IsAtom(char* Str) {
//------------------------------------------------------------------------
// see if Str starts with "(segname, resnum, atomname)"
// return:    the number of characters, including white space.
//            0, if Str does not start with "(segname, resnum, atomname)"
//------------------------------------------------------------------------
  int    Count;
  char*  FullString=Str;

  // if char following the parentheses is ':', this isn't an atom
  if (IsAtomNameList(Str)) {
    return(0);
  }
  // Str must contain the following in sequence to be a legit atom
  // <ws> = optional white-space
  // '(', <ws>, alphanumeric, <ws>, numeric, <ws>, alphanumeric, <ws>, ')'
  if (Str[0] == '(') {
    Str += 1;
    Str += ReadWhite(Str);
    Count = ReadAlphaNum(Str);
    if (Count) {
      Str += Count;
      Str += ReadWhite(Str);
      Count = ReadDigits(Str);
      if (Count) {
        Str += Count;
        Str += ReadWhite(Str);
        Count = ReadAlphaNum(Str);
        if (Count) {
          Str += Count;
          Str += ReadWhite(Str);
          if (Str[0] == ')') {
            Str += 1;
            return(Str-FullString);
          }
        }
      }
    }
  }
  return(0);
}


int IsResidueRange(char* Str) {
//------------------------------------------------------------------------
// see if Str starts with "(segname, resnum) to (segname, resnum)"
// return:    the number of characters, including white space.
//            0, if Str does not start with "(sn, rn) to (sn, rn)"
//------------------------------------------------------------------------
  int    Count;
  char*  FullString=Str;

  // Str must contain the following in sequence to be a legit res range
  // <ws> = optional white-space
  // residue, <ws>, "to", <ws>, residue
  Count = IsAResidue(Str);
  if (Count) {
    Str += Count;
    Str += ReadWhite(Str);
    if (strncmp(Str, "to", 2) == 0) {
      Str += 2;
      Str += ReadWhite(Str);
      Count = IsAResidue(Str);
      if (Count) {
        Str += Count;
        return(Str-FullString);
      }
    }
  }
  return(0);
}


int IsResidue(char* Str) {
//------------------------------------------------------------------------
// see if Str starts with "(segname, resnum)"
// but not "(segname, resnum) to (segname, resnum)"
//------------------------------------------------------------------------
  int Count;

  Count = IsAResidue(Str);
  if (Count) {
    if (IsResidueRange(Str)) {
      return(0);
    }
    else {
      return(Count);
    }
  }
  return(0);
}


int IsAResidue(char* Str) {
//------------------------------------------------------------------------
// see if Str starts with "(segname, resnum)"
// return:    the number of characters, including white space.
//            0, if Str does not start with "(segname, resnum)"
//------------------------------------------------------------------------
  int    Count;
  char*  FullString=Str;

  // if char following the parentheses is ':', this isn't a residue
  if (IsAtomNameList(Str)) {
    return(0);
  }
  // Str must contain the following in sequence to be a legit residue
  // <ws> = optional white-space
  // '(', <ws>, alphanumeric, <ws>, numeric, <ws>, ')'
  if (Str[0] == '(') {
    Str += 1;
    Str += ReadWhite(Str);
    Count = ReadAlphaNum(Str);
    if (Count) {
      Str += Count;
      Str += ReadWhite(Str);
      Count = ReadDigits(Str);
      if (Count) {
        Str += Count;
        Str += ReadWhite(Str);
        if (Str[0] == ')') {
          Str += 1;
          return(Str-FullString);
        }
      }
    }
  }
  return(0);
}


int ReadParentheses(const char* Str) {
//-------------------------------------------------------------------
// count the number of characters from the leading '('
// to the first ')' of Str (inclusive).
//     no leading '('  =>  return  0
//     no closing ')'  =>  return -1
//-------------------------------------------------------------------
  const char* Str2;

  if (Str[0] != '(') {
    return(0);
  }
  Str2 = strchr(Str, ')');
  if (Str2 == NULL) {
    return(-1);
  }
  return((Str2-Str)+1);
}


int ReadAlpha(const char* Str) {
//-------------------------------------------------------------------
// determine the leading number of alphabetic characters in Str.
//-------------------------------------------------------------------
  int  i=0;

  while (1) {
    if (isalpha(Str[i]) || Str[i]=='\'' || Str[i]=='\"' || Str[i] == '*') {
      i++;
    }
    else {
      break;
    }
  }
  return(i);
}


int ReadAlphaNum(const char* Str) {
//-------------------------------------------------------------------
// determine the leading number of alphanumeric characters in Str.
//-------------------------------------------------------------------
  int  i=0;

  while (1) {
    if (isalnum(Str[i]) || Str[i]=='\'' || Str[i]=='\"' || Str[i] == '*') {
      i++;
    }
    else {
      break;
    }
  }
  return(i);
}


int ReadDigits(const char* Str) {
//-------------------------------------------------------------------
// determine the leading number of numeric characters in Str.
//-------------------------------------------------------------------
  int  i=0;

  while (1) {
    if (isdigit(Str[i])) {
      i++;
    }
    else {
      break;
    }
  }
  return(i);
}


int ReadWhite(const char* Str) {
//-------------------------------------------------------------------
// determine the leading number of white characters in Str.
// a white char is:
//   space, tab, linefeed, carriage-return, formfeed,
//   vertical-tab, and newline characters.
// a white char is also (for the sake of this program):
//   comma, semi-colon, and period.
//-------------------------------------------------------------------
  int  i=0;

  // count the number of white-char's.  break at first non-white-char.
  while (1) {
    if (
         (Str[i] == 9)   ||      // tab
         (Str[i] == 10)  ||      // LF
         (Str[i] == 11)  ||      // vertical-tab
         (Str[i] == 12)  ||      // form-feed
         (Str[i] == 13)  ||      // CR
         (Str[i] == 32)  ||      // space
         (Str[i] == ',') ||      // comma
         (Str[i] == ';')         // semi-colon
//         (Str[i] == '.')         // period  (took this out in case of =.4, eg)
       )
    {    
      i++;
    }
    else {
      break;
    }
  }
  return(i);
}


void ToLower(char* Str) {
//-------------------------------------------------------------------
// convert Str to all lower case
//-------------------------------------------------------------------
  for (unsigned int i=0; i<strlen(Str); i++) {
    Str[i] = (char)tolower(Str[i]);
  }
}

