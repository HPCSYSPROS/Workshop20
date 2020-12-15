/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeHomePatches.h"
#include "PatchMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

ComputeHomePatches::ComputeHomePatches(ComputeID c) : Compute(c) {
  patchMap = PatchMap::Object();
  useAvgPositions = 0;
  hasPatchZero = 0;
}

ComputeHomePatches::~ComputeHomePatches()
{
  ;
}

void ComputeHomePatches::initialize()
{
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  patchList.resize(0);

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    patchList.add(PatchElem((*ai).patch, this, useAvgPositions));
    if ( (*ai).patch->getPatchID() == 0 ) hasPatchZero = 1;
  }

  setNumPatches(patchList.size());
}

void ComputeHomePatches::atomUpdate()
{
  Compute::atomUpdate();
}

