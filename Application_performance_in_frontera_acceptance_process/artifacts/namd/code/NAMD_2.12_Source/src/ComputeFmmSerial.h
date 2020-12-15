/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEFMMSERIAL_H
#define COMPUTEFMMSERIAL_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;
class FmmSerialForceMsg;
class ComputeFmmSerialAtom;

class ComputeFmmSerial : public ComputeHomePatches {
public:
  ComputeFmmSerial(ComputeID c);
  virtual ~ComputeFmmSerial();
  void doWork();
  void saveResults(FmmSerialForceMsg *);

 private:
  SubmitReduction *reduction;

};

#endif

