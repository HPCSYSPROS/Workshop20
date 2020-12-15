/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute object which deals with a single patch.
*/

#ifndef COMPUTEPATCH_H
#define COMPUTEPATCH_H

#include "Compute.h"
#include "PatchTypes.h"

#include "Box.h"
#include "OwnerBox.h"

class Patch;
class Node;
class PatchMap;

class ComputePatch : public Compute {

public:
  ComputePatch(ComputeID c, PatchID pid);
  virtual ~ComputePatch();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();

protected :
  int numAtoms;
  CompAtomExt *pExt;
  CompAtom* p;
  Results* r;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int numWaterAtoms;
  #endif

  virtual void doForce(CompAtom* p, CompAtomExt* pExt, Results* r) = 0;
  Patch *patch;

// private:
  PatchID patchID;
  Box<Patch,CompAtom> *positionBox;
  Box<Patch,Results> *forceBox;

};

#endif

