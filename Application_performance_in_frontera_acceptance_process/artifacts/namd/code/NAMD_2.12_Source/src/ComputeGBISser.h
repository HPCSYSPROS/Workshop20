/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*******************************************************************************
 *******************************************************************************
  This serial version of GBIS is out of date and hasn't been vetted;
  it is to be used as a way of testing new implicit solvent models
  since all atomic coordinates are gathered into a single Compute.
 *******************************************************************************
 ******************************************************************************/

#ifndef COMPUTEGBIS_H
#define COMPUTEGBIS_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;
class GBISForceMsg;
class ComputeGBISAtom;

class ComputeGBISser : public ComputeHomePatches {
public:
  ComputeGBISser(ComputeID c);
  virtual ~ComputeGBISser();
  void doWork();
  void saveResults(GBISForceMsg *);

 private:
  SubmitReduction *reduction;
};

#endif

