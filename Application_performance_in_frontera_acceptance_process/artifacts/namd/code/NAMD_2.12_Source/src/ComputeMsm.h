/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMSM_H
#define COMPUTEMSM_H

#include "Lattice.h"
#include "ComputeMsmMgr.decl.h"
#include "ComputeHomePatches.h"
#include "NamdTypes.h"


class ComputeMsmMgr;
class SubmitReduction;

class MsmInitMsg : public CMessage_MsmInitMsg {
  public:
    ScaledPosition smin, smax;  // needs the extreme positions
};


class ComputeMsm : public ComputeHomePatches {
public:
  ComputeMsm(ComputeID c);
  virtual ~ComputeMsm();
  void doWork();
  void saveResults();

  void setMgr(ComputeMsmMgr *mgr) { myMgr = mgr; }

private:
  SubmitReduction *reduction;
  ComputeMsmMgr *myMgr;  // points to the local MSM manager

  Real qscaling;  // charge scaling constant

  int numLocalPatches;   // total number of local patches to expect
  int cntLocalPatches;   // count local patches into saveResults()
};


#endif // COMPUTEMSM_H
