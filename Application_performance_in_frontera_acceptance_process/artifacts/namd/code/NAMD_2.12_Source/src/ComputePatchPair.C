/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputePatchPair.h"
#include "Priorities.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

ComputePatchPair::ComputePatchPair(ComputeID c, PatchID p[], int t[]) 
    : Compute(c) {

  gbisPhase = 3;
  setNumPatches(2);

  for (int i=0; i<2; i++) {
      patchID[i] = p[i];
      trans[i] = t[i];
      patch[i] = NULL;
      positionBox[i] = NULL;
      forceBox[i] = NULL;
  }
}

ComputePatchPair::~ComputePatchPair() {
  DebugM(4, "~ComputePatchPair("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms("<<patchID[1]<<") = " << numAtoms[1] << "\n" );
  DebugM(4, "~ComputePatchPair("<<cid<<") addr("<<patchID[0]<<") = " 
    << PatchMap::Object()->patch(patchID[0]) << " addr("<<patchID[1]<<") = "
    << PatchMap::Object()->patch(patchID[1]) << "\n");
  for (int i=0; i<2; i++) {
    if (positionBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterPositionPickup(this,
	 &positionBox[i]);
    }
    if (forceBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterForceDeposit(this,
		&forceBox[i]);
    }
  }

}

void ComputePatchPair::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

    for (int i=0; i<2; i++) {
	if (positionBox[i] == NULL) { // We have yet to get boxes
	    if (!(patch[i] = PatchMap::Object()->patch(patchID[i]))) {
	      DebugM(5,"invalid patch(" << patchID[i] 
		   << ")  pointer!\n");
	    }
	    positionBox[i] = patch[i]->registerPositionPickup(this);
	    forceBox[i] = patch[i]->registerForceDeposit(this);
	}
	numAtoms[i] = patch[i]->getNumAtoms();
    }

  DebugM(4, "initialize("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms(" <<patchID[1]<<") = " << numAtoms[1] << "\n" );

    Compute::initialize();

    // proxies are more urgent (lower priority) than patches
    int myNode = CkMyPe();
    int p0 = PATCH_PRIORITY(patchID[0]);
    if ( PatchMap::Object()->node(patchID[0]) == myNode ) {
      p0 += GB1_COMPUTE_HOME_PRIORITY;
    } else {
      p0 += GB1_COMPUTE_PROXY_PRIORITY;
    }
    int p1 = PATCH_PRIORITY(patchID[1]);
    if ( PatchMap::Object()->node(patchID[1]) == myNode ) {
      p1 += GB1_COMPUTE_HOME_PRIORITY;
    } else {
      p1 += GB1_COMPUTE_PROXY_PRIORITY;
    }
    if (p0<p1) { //base phase priorities off of p0
      if ( PatchMap::Object()->node(patchID[0]) == myNode ) {
        gbisPhasePriority[0] = 0;
        gbisPhasePriority[1] = GB2_COMPUTE_HOME_PRIORITY-GB1_COMPUTE_HOME_PRIORITY;
        gbisPhasePriority[2] = COMPUTE_HOME_PRIORITY-GB1_COMPUTE_HOME_PRIORITY;
      } else {
        gbisPhasePriority[0] = 0;
        gbisPhasePriority[1] = GB2_COMPUTE_PROXY_PRIORITY-GB1_COMPUTE_PROXY_PRIORITY;
        gbisPhasePriority[2] = COMPUTE_PROXY_PRIORITY-GB1_COMPUTE_PROXY_PRIORITY;
      }
    } else { //base phase priorities off of p1
      if ( PatchMap::Object()->node(patchID[1]) == myNode ) {
        gbisPhasePriority[0] = 0;
        gbisPhasePriority[1] = GB2_COMPUTE_HOME_PRIORITY-GB1_COMPUTE_HOME_PRIORITY;
        gbisPhasePriority[2] = COMPUTE_HOME_PRIORITY-GB1_COMPUTE_HOME_PRIORITY;
      } else {
        gbisPhasePriority[0] = 0;
        gbisPhasePriority[1] = GB2_COMPUTE_PROXY_PRIORITY-GB1_COMPUTE_PROXY_PRIORITY;
        gbisPhasePriority[2] = COMPUTE_PROXY_PRIORITY-GB1_COMPUTE_PROXY_PRIORITY;
      }
    }
    basePriority = ((p0<p1)?p0:p1);   // most urgent wins

}

void ComputePatchPair::atomUpdate() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

    // DebugM(4,"atomUpdate() - positionBox[0] is " << positionBox[0] << "\n");
    for (int i=0; i<2; i++) {

	numAtoms[i] = patch[i]->getNumAtoms();

        // DMK - Atom Separation (water vs. non-water)
        #if NAMD_SeparateWaters != 0
	  numWaterAtoms[i] = patch[i]->getNumWaterAtoms();
        #endif
    }



    // Compute::atomUpdate();
}


void ComputePatchPair::doForce(CompAtom* p[2], CompAtomExt* pExt[2], Results* r[2])
{
    CkPrintf("ComputePatchPair::doForce() - Dummy eval was sent\n");
    CkPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
}

//---------------------------------------------------------------------
// Where the actual computation is invoked.  doForce is 
// overloaded with specific calculation
//---------------------------------------------------------------------
void ComputePatchPair::doWork() {

#ifndef NAMD_CUDA
  LdbCoordinator::Object()->startWork(ldObjHandle);
#endif
  if ( ( computeType != computeNonbondedPairType ) ||
        (!patch[0]->flags.doGBIS || gbisPhase == 1) ) {
    // Open up positionBox, forceBox, and atomBox
    for (int i=0; i<2; i++) {
      p[i] = positionBox[i]->open();
      r[i] = forceBox[i]->open();
      pExt[i] = patch[i]->getCompAtomExtInfo();
    }
  }

  doForce(p, pExt, r);
 // Inform load balancer
#ifndef NAMD_CUDA
  if (patch[0]->flags.doGBIS && (gbisPhase == 1 || gbisPhase == 2)) {
    LdbCoordinator::Object()->pauseWork(ldObjHandle);
  } else {
    LdbCoordinator::Object()->endWork(ldObjHandle);
  }
#endif

  // Close up boxes
  if ( ( computeType != computeNonbondedPairType ) ||
      (!patch[0]->flags.doGBIS || gbisPhase == 3) ) {
    for (int i=0; i<2; i++) {
      positionBox[i]->close(&p[i]);
      forceBox[i]->close(&r[i]);
    }
  }

}

