/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute object which deals with a single patch.
*/

#include "InfoStream.h"
#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputeHomePatch.h"
#include "PatchMap.inl"
#include "HomePatch.h"
#include "Priorities.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

ComputeHomePatch::ComputeHomePatch(ComputeID c, PatchID p) : Compute(c) {
    setNumPatches(1);
    patchID = p;
    patch = NULL;
    homePatch = NULL;
    positionBox = NULL;
    forceBox = NULL;
}

ComputeHomePatch::~ComputeHomePatch() {
  DebugM(4, "~ComputeHomePatch("<<cid<<") numAtoms("<<patchID<<") = " 
    << numAtoms << "\n");
    if (positionBox != NULL) {
      PatchMap::Object()->patch(patchID)->unregisterPositionPickup(this,
	 &positionBox);
    }
    if (forceBox != NULL) {
      PatchMap::Object()->patch(patchID)->unregisterForceDeposit(this,
		&forceBox);
    }
}

void ComputeHomePatch::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

	if (positionBox == NULL) { // We have yet to get boxes
	    if (!(patch = PatchMap::Object()->patch(patchID))) {
	      NAMD_bug("ComputeHomePatch used with unknown patch.");
	    }
            if (!(homePatch = PatchMap::Object()->homePatch(patchID))) {
	      NAMD_bug("ComputeHomePatch used with proxy.");
	    }
	    DebugM(3, "initialize(" << cid <<")  patchid = "<<patch->getPatchID()<<"\n");
	    positionBox = patch->registerPositionPickup(this);
	    forceBox = patch->registerForceDeposit(this);
	}
	numAtoms = patch->getNumAtoms();

  DebugM(3, "initialize("<<cid<<") numAtoms("<<patchID<<") = " 
    << numAtoms  << " patchAddr=" << patch << "\n");
    Compute::initialize();

    int myNode = CkMyPe();
    basePriority = COMPUTE_HOME_PRIORITY + PATCH_PRIORITY(patchID);
}

void ComputeHomePatch::atomUpdate() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?
    numAtoms = patch->getNumAtoms();
}

void ComputeHomePatch::doWork() {
  CompAtom* p;
  Results* r;
  FullAtom* a = homePatch->getAtomList().begin();

  DebugM(3,patchID << ": doWork() called.\n");

  // Open up positionBox, forceBox, and atomBox
  p = positionBox->open();
  r = forceBox->open();

  // Pass pointers to doForce
  doForce(a,r);

  // Close up boxes
  positionBox->close(&p);
  forceBox->close(&r);

  DebugM(2,patchID << ": doWork() completed.\n");
}

