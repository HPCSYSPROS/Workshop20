/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEFULLDIRECT_H
#define COMPUTEFULLDIRECT_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;

class ComputeFullDirect : public ComputeHomePatches {
private:
  SubmitReduction *reduction;
public:
  ComputeFullDirect(ComputeID c);
  virtual ~ComputeFullDirect();
  void doWork();
};

#endif

