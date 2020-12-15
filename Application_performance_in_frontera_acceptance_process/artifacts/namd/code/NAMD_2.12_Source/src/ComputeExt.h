/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEEXT_H
#define COMPUTEEXT_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;
class ExtForceMsg;
class ComputeExtAtom;

class ComputeExt : public ComputeHomePatches {
public:
  ComputeExt(ComputeID c);
  virtual ~ComputeExt();
  void doWork();
  void saveResults(ExtForceMsg *);

 private:
  SubmitReduction *reduction;

};

#endif

