/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute object which deals with a single patch.
*/

#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputePatch.h"
#include "Priorities.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "ComputeMap.h" //needed for checking GBIS type
#include "LdbCoordinator.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

ComputePatch::ComputePatch(ComputeID c, PatchID p) : Compute(c) {
    setNumPatches(1);
    patchID = p;
    patch = NULL;
    positionBox = NULL;
    forceBox = NULL;
    gbisPhase = 3;
}

ComputePatch::~ComputePatch() {
  DebugM(4, "~ComputePatch("<<cid<<") numAtoms("<<patchID<<") = " 
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

void ComputePatch::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

	if (positionBox == NULL) { // We have yet to get boxes
	    if (!(patch = PatchMap::Object()->patch(patchID))) {
	      NAMD_bug("ComputePatch used with unknown patch.");
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
    if ( PatchMap::Object()->node(patchID) != myNode ) {
      basePriority = GB1_COMPUTE_PROXY_PRIORITY + PATCH_PRIORITY(patchID);
      gbisPhasePriority[0] = 0;
      gbisPhasePriority[1] = GB2_COMPUTE_PROXY_PRIORITY-GB1_COMPUTE_PROXY_PRIORITY;//sub GB1_PRIOR
      gbisPhasePriority[2] = COMPUTE_PROXY_PRIORITY-GB1_COMPUTE_PROXY_PRIORITY;
    } else {
      basePriority = GB1_COMPUTE_HOME_PRIORITY + PATCH_PRIORITY(patchID);
      gbisPhasePriority[0] = 0;
      gbisPhasePriority[1] = GB2_COMPUTE_HOME_PRIORITY-GB1_COMPUTE_HOME_PRIORITY;
      gbisPhasePriority[2] = COMPUTE_HOME_PRIORITY-GB1_COMPUTE_HOME_PRIORITY;
    }
}

void ComputePatch::atomUpdate() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?
    numAtoms = patch->getNumAtoms();

    // DMK - Atom Separation (water vs. non-water)
    #if NAMD_SeparateWaters != 0
      numWaterAtoms = patch->getNumWaterAtoms();
    #endif
}

void ComputePatch::doWork() {
  DebugM(3,patchID << ": doWork() called.\n");

#ifndef NAMD_CUDA
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  LdbCoordinator::Object()->startWork(ldObjHandle);
#endif
  if ( (computeType != computeNonbondedSelfType ) ||
       (!patch->flags.doGBIS || gbisPhase == 1) ) {
    // Open up positionBox, forceBox
    p = positionBox->open();
    r = forceBox->open();
    pExt = patch->getCompAtomExtInfo();
  }

  // Pass pointers to doForce
  doForce(p, pExt, r);
// Inform load balancer
#ifndef NAMD_CUDA
  if (patch->flags.doGBIS && (gbisPhase == 1 || gbisPhase == 2)) {
    LdbCoordinator::Object()->pauseWork(ldObjHandle);
  } else {
    LdbCoordinator::Object()->endWork(ldObjHandle);
  }
#endif
  // Close up boxes
  if ( (computeType != computeNonbondedSelfType   ) ||
       (!patch->flags.doGBIS || gbisPhase == 3) ) {
    positionBox->close(&p);
    forceBox->close(&r);
  }
  DebugM(2,patchID << ": doWork() completed.\n");
}

