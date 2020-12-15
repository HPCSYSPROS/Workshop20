/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#if !defined(PARSE_HPP)
  #define PARSE_HPP

enum item_t  {kAtom, kResidue, kAtomName, kAtomNameList,
              kResidueRange, kStartGroup, kEndGroup, kUnknownItem};
enum restr_t {kPosi,      kDist,      kAngle,      kDihe,
              kPosiBound, kDistBound, kAngleBound, kDiheBound,
              kPosiPMF,   kDistPMF,   kAnglePMF,   kDihePMF, kUnknownRestr};
enum pmf_t   {kTask, kTime, kLambda, kLambdaT, kPrint, kNoPrint,
              kEquilTime, kAccumTime, kNumRepeats, kUnknownPmf};
enum TimeUnits_t  {k_fs, k_ps, k_ns, kUnknownTime};

// the "main" function for parsing
void   ReadInput(char* Str, ARestraintManager& RMgr,
                            ALambdaManager& LMgr,
                            GlobalMasterFreeEnergy& CFE,
                            double dT);

// get an initialized restraint that's read from the input String
int    ReadRestraints(char* Str, ARestraintManager& RMgr,
                                 GlobalMasterFreeEnergy& CFE);
ARestraint* GetRestraint(char* Str, int& NumChars, GlobalMasterFreeEnergy& CFE);

// for reading pmf/mcti blocks
int    ReadPmfBlock(char* Str, ALambdaControl& PmfBlock, double dT);
int    ReadNextPmfSpec(char* Str, pmf_t& PmfSpec);
int    ReadTaskType(char* Str, feptask_t& Task);
int    ReadTimeUnits(char* Str, TimeUnits_t& Units, TimeUnits_t DefaultUnits);
double GetTime(double Val, TimeUnits_t Units);

// functions for parsing the config file
void    CheckParentheses(const char* Str);
void    ProblemParsing(const char* Message, const char* Str, Bool_t Terminate=kTrue);
void    ToLower(char* Str);
int     ReadWhite(const char* Str);
int     ReadAlphaNum(const char* Str);
int     ReadAlpha(const char* Str);
int     ReadDigits(const char* Str);
int     ReadParentheses(const char* Str);
int     IsStartGroup(char* Str);
int     IsEndGroup(char* Str);
int     IsAtomName(char* Str);
int     IsAtomNameList(char* Str);
int     IsAtom(char* Str);
int     IsAResidue(char* Str);
int     IsResidue(char* Str);
int     IsResidueRange(char* Str);
int     ReadWord(const char* Str, const char* Word, Bool_t ErrMsg=kFalse);
int     ReadChar(char* Str, char Char, Bool_t ErrMsg=kFalse);
int     ReadAValue(char* Str, double& Value, Bool_t ErrMsg=kFalse);
int     ReadBound(char* Str, Bound_t& Bound);
item_t  ReadNextItem(char* Str, int& NumChars);
restr_t ReadNextRestraintType(char* Str, int& NumChars);

// functions for putting the specified atoms into Group's
int  AddAtoms(AGroup& Group, char* Str, GlobalMasterFreeEnergy& CFE);
void AddAtom(AGroup& Group, char* Atom, GlobalMasterFreeEnergy& CFE);
void AddResidues(AGroup& Group, char* ResRange, GlobalMasterFreeEnergy& CFE);
void AddAtomsInResidues(AGroup& Group, char* AtomNames, char* ResRange,
                        GlobalMasterFreeEnergy& CFE);

// functions for helping with the above
void AddAtom(AGroup& Group, char* ResRange, char* AtomName,
             GlobalMasterFreeEnergy& CFE);
void GetResRange(char* ResRange, int& ResNum1, int& ResNum2);
int  GetSegName(char* Str, char* SegName);
int  GetResNum(char* Str, int& ResNum);
int  GetAtomName(char* Str, char* AtomName);

#endif
