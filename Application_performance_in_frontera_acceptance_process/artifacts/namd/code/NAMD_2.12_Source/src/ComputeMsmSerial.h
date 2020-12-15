/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMSMSERIAL_H
#define COMPUTEMSMSERIAL_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;
class MsmSerialForceMsg;
class ComputeMsmSerialAtom;

class ComputeMsmSerial : public ComputeHomePatches {
public:
  ComputeMsmSerial(ComputeID c);
  virtual ~ComputeMsmSerial();
  void doWork();
  void saveResults(MsmSerialForceMsg *);

 private:
  SubmitReduction *reduction;

};

#endif

