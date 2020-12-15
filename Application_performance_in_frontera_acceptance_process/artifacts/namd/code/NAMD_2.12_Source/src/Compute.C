/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Top of Compute hierarchy.  
   enqueueWork() - delivers Compute object itself to queue up for doWork()
   doWork() - called by work queue
*/

#include "main.h"
#include "charm++.h"

#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"

#include "NamdTypes.h"
#include "Box.h"
#include "OwnerBox.h"

#include "Node.h"
#include "Compute.h"
#include "Priorities.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

Compute::Compute(ComputeID c) : gbisPhase(1),basePriority(0), cid(c),
	localWorkMsg(new (PRIORITY_SIZE) LocalWorkMsg) { 
  gbisPhasePriority[0] = 0;
  gbisPhasePriority[1] = 0;
  gbisPhasePriority[2] = 0;
  doAtomUpdate = false;
  computeType = ComputeMap::Object()->type(c);
  ldObjHandle.id.id[0] = -1;
}

Compute::~Compute() {
  delete localWorkMsg;
}

void Compute::enqueueWork() {
  if (!this) { DebugM(4,"This Compute is NULL!!!\n"); }
  if ( ! noWork() ) {
    //gbisPhase = 1; //first phase - this should already be 1
    WorkDistrib::messageEnqueueWork(this);  // should be in ComputeMgr?
  } else {
    //don't enqueue work
  }
}


//---------------------------------------------------------------------
// Signal from patch or proxy that data is ready.
// When all Patches and Proxies needed by this Compute object
// have checked-in, we are ready to enqueueWork()
//---------------------------------------------------------------------
void Compute::patchReady(PatchID patchID, int doneMigration, int seq) { 
  if (doneMigration) { // If any patch has done migration - we must remap
    doAtomUpdate = true; 
  }

  if (numPatches <= 0) {
      DebugM(5,"Compute::patchReady("<<patchID<<")-call not valid!\n");
  } else {
    if (! --patchReadyCounter) {
      patchReadyCounter = numPatches;
      //gbisPhase = 1;
      sequenceNumber = seq;  // breaks CUDA priority if done earlier
      if (doAtomUpdate) {
	atomUpdate();
	doAtomUpdate = false;
      }
      enqueueWork();
    }
  }
}

void Compute::gbisP2PatchReady(PatchID pid, int seq) {

  if (! --patchReadyCounter) {
    patchReadyCounter = numPatches;
    //gbisPhase = 2;
    sequenceNumber = seq;
    enqueueWork();
  }
}

void Compute::gbisP3PatchReady(PatchID pid, int seq) {
  if (! --patchReadyCounter) {
    patchReadyCounter = numPatches;
    //gbisPhase = 3;
    sequenceNumber = seq;
    enqueueWork();
  }
}


int Compute::noWork() {
  return 0;
}

void Compute::doWork() {
    DebugM(5,"Default Compute::doWork() called.\n");
}

void Compute::finishPatch(int) {
    DebugM(5,"Default Compute::finishPatch() called.\n");
}

