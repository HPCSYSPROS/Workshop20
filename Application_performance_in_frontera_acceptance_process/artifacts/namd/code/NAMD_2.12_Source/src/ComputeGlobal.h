/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef COMPUTEGLOBAL_H
#define COMPUTEGLOBAL_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeMgr;
class SubmitReduction;

struct intpair {
  int first, second;
  intpair() {;}
  intpair(int f, int s) : first(f), second(s) {;}
};

inline bool operator<(const intpair &lhs, const intpair &rhs) {
  return lhs.first < rhs.first ? true :
         lhs.first != rhs.first ? false :
         lhs.second < rhs.second;
}

class ComputeGlobal : public ComputeHomePatches {
public:
  ComputeGlobal(ComputeID, ComputeMgr*);
  virtual ~ComputeGlobal();
  void doWork();
  // void recvConfig(ComputeGlobalConfigMsg *);
  void recvResults(ComputeGlobalResultsMsg *);
  // For "loadtotalforces" TCL command
  void saveTotalForces(HomePatch *);

private:
  ComputeMgr *comm;

  void sendData();
  void configure(AtomIDList &newaid, AtomIDList &newgdef);

  AtomIDList aid;
  AtomIDList gdef;  // definitions of groups
  ResizeArray<intpair> gpair;
  
  // (For "loadtotalforces" TCL command)
  // The atom IDs and forces of the requested atoms on the node
  // after force evaluation. "fid" could be slightly different
  // from "aid", since the latter is after atom migration.
  AtomIDList fid;
  ForceList totalForce;
  ForceList groupTotalForce;
  int numGroupsRequested;

  Force **forcePtrs;
  FullAtom **atomPtrs;
  
  int forceSendEnabled; // are total forces received?
  int forceSendActive; // are total forces received this step?
  int gfcount;  // count of atoms contributing to group forces
  char *isRequested;  // whether this atom is requested by the TCL script
  int isRequestedAllocSize;  // size of array
  int endRequested;  // starting at this point assume not requested
  int dofull;  // whether "Results::slow" force will exist

  int firsttime;
  SubmitReduction *reduction;
};

#endif

