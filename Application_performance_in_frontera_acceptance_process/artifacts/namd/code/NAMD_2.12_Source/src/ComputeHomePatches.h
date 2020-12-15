/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEHOMEPATCHES_H
#define COMPUTEHOMEPATCHES_H

#include "NamdTypes.h"
#include "common.h"
#include "Compute.h"
#include "HomePatch.h"

#include "Box.h"
#include "OwnerBox.h"

class PatchElem {
  public:
    PatchID patchID;
    HomePatch *p;
    Box<Patch,CompAtom> *positionBox;
    Box<Patch,CompAtom> *avgPositionBox;
    Box<Patch,Results> *forceBox;
    CompAtom *x;
    Results *r;
    Force *f;

  PatchElem() {
    patchID = -1;
    p = NULL;
    positionBox = NULL;
    avgPositionBox = NULL;
    forceBox = NULL;
    x = NULL;
    r = NULL;
    f = NULL;
  }

  PatchElem(PatchID p_param) {
    patchID = p_param;
  }

  PatchElem(HomePatch *p_param, Compute *cid, int useAvgPos) {
    patchID = p_param->getPatchID();
    p = p_param;
    positionBox = p_param->registerPositionPickup(cid);
    if ( useAvgPos ) {
      avgPositionBox = p_param->registerAvgPositionPickup(cid);
    }
    forceBox = p_param->registerForceDeposit(cid);
    x = NULL;
    r = NULL;
    f = NULL;
  }
    
  ~PatchElem() {};

  int operator==(const PatchElem &elem) const {
    return (elem.patchID == patchID);
  }

  int operator<(const PatchElem &elem) const {
    return (patchID < elem.patchID);
  }
};

typedef UniqueSortedArray<PatchElem> ComputeHomePatchList;

class ReductionMgr;

class ComputeHomePatches : public Compute {
protected:
  int useAvgPositions;
  int hasPatchZero;  // used for calculations/data needed only once

  ComputeHomePatchList patchList;

  PatchMap *patchMap;

public:
  ComputeHomePatches(ComputeID c);
  virtual ~ComputeHomePatches();
  virtual void initialize();
  virtual void atomUpdate();
  Flags *getFlags(void) { return &(patchList[0].p->flags); }
};

#endif

