/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Primary class for pairwise force calculations over pairs of patches.
   Takes care of boxes, depositing of forces etc.
*/

#ifndef COMPUTEPPAIR_H
#define COMPUTEPPAIR_H

#include "Compute.h"
#include "PatchTypes.h"

#include "Box.h"
#include "OwnerBox.h"

class Patch;
class Node;
class PatchMap;

class ComputePatchPair : public Compute {

public:
  ComputePatchPair(ComputeID c, PatchID pid[], int t[]);
  virtual ~ComputePatchPair();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();

protected :
  int numAtoms[2];
  CompAtomExt *pExt[2];
  CompAtom* p[2];
  Results* r[2];

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int numWaterAtoms[2];
  #endif

  virtual void doForce(CompAtom* p[2], CompAtomExt* pExt[2], Results* r[2]);
  Patch *patch[2];

// private: // hack for ComputeNonbondedPair::noWork()
  PatchID patchID[2];
  int trans[2];
  Box<Patch,CompAtom> *positionBox[2];
  Box<Patch,Results> *forceBox[2];
};

#endif

