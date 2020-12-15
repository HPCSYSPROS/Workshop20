/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute object which deals with a single patch.
*/

#ifndef COMPUTEHOMEPATCH_H
#define COMPUTEHOMEPATCH_H

#include "Compute.h"
#include "PatchTypes.h"

#include "Box.h"
#include "OwnerBox.h"

class Patch;
class HomePatch;
class Node;
class PatchMap;

class ComputeHomePatch : public Compute {

public:
  ComputeHomePatch(ComputeID c, PatchID pid);
  virtual ~ComputeHomePatch();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();

protected :
  int numAtoms;
  virtual void doForce(FullAtom* p, Results* r) = 0;
  Patch *patch;
  HomePatch *homePatch;

private:
  PatchID patchID;
  Box<Patch,CompAtom> *positionBox;
  Box<Patch,Results> *forceBox;

};

#endif

