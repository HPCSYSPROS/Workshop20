/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Base class for all Compute objects. Almost an abstract
   class except that it does do the basic patchReady()
   countdown.
*/

#ifndef COMPUTE_H
#define COMPUTE_H

#include "main.h"
#include "charm++.h"

#include "NamdTypes.h"

class Node;
class PatchMap;
class LocalWorkMsg;

// Base class for various forms of Compute objects
// for example: <linkto class=ComputeAngles>ComputeAngles</linkto> 
// and <linkto class=ComputeNonbondedExcl>ComputeNonbondedExcl</linkto>
class Compute {
private:
  int patchReadyCounter;
  int numPatches;
  int doAtomUpdate;
  int sequenceNumber;

protected:
  int computeType;
  int basePriority;
  void enqueueWork();
  int gbisPhase;//earlier phases have higher priority
  int gbisPhasePriority[3];//earlier phases have higher priority

public:
  const ComputeID cid;
  LDObjHandle ldObjHandle;

  LocalWorkMsg *const localWorkMsg;
  Compute(ComputeID);
  int type() { return computeType; };

  virtual ~Compute();

  void setNumPatches(int n) { patchReadyCounter = numPatches = n; }
  int getNumPatches() { return (numPatches); };

  // registers for boxes
  virtual void initialize() {};
  // destructor better unregister for boxes!

  virtual void atomUpdate() {};
  virtual void patchReady(PatchID, int doneMigration, int seq);
  virtual int noWork(); // cleans up and returns 1 if no work to do
  virtual void doWork(); // actually does the work if noWork() returns 0
  virtual void finishPatch(int);
  int sequence(void) { return sequenceNumber; }
  int priority(void) { return basePriority+gbisPhasePriority[gbisPhase-1]; }
  int getGBISPhase(void) {return gbisPhase;}

  virtual void gbisP2PatchReady(PatchID, int seq);
  virtual void gbisP3PatchReady(PatchID, int seq);

};

/* For projection's usage: each compute object's work is associated 
 * with a user event in projections. This macro indicates the offset 
 * of the event ID that those compute objects' user events begin with.
 */ 
#define TRACE_COMPOBJ_IDOFFSET 10000

#endif

