/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERTMD_H
#define GLOBALMASTERTMD_H

#include "GlobalMaster.h"
#include "PDBData.h"
#include <map>
//#include <ext/hash_map>
#include <vector>
//namespace std { using namespace __gnu_cxx; }

class GlobalMasterTMD : public GlobalMaster {
public:
  GlobalMasterTMD();
  ~GlobalMasterTMD();

private:
 
  void calculate();
  void parseAtoms(const char *file, int, bool);
  void NewTarget(int);
  int numatoms;
  Vector *atompos;
  Vector *atompos2;
  char *altloc;
  bool qDiffRMSD;
  std::map <int, std::vector<int>  > dmap;
  std::map <int, BigReal > kmap;
  std::map <int, int > altlocmap;
  int numTMDatoms;
  BigReal K;
  BigReal initialRMS, finalRMS;
  int outputFreq;
  int currentStep, firstStep, lastStep;
  BigReal *target;
  BigReal *target2;
  BigReal *weight;
  //int *target_aid;
  // mapping of atom id's to array positions
  //int *aidmap;
};
#endif
