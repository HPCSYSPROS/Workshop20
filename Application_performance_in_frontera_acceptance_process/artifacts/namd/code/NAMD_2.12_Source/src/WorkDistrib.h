/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef _WORKDISTRIB_H
#define _WORKDISTRIB_H

#include "charm++.h"

#include "main.h"

#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ComputeMap.h"
#include "WorkDistrib.decl.h"

class Node;
class Compute;
class Molecule;

// For Compute objects to enqueue themselves when ready to compute
class LocalWorkMsg : public CMessage_LocalWorkMsg
{
public:
  Compute *compute;
};

class FinishWorkMsg : public CMessage_FinishWorkMsg
{
public:
  Compute *compute;
  int data;
};

enum { maxPatchDepends = 126 };

class PatchMapMsg;
class ComputeMapMsg;
class ComputeMapChangeMsg;

class WorkDistrib : public CBase_WorkDistrib
{
public:
  WorkDistrib();
  ~WorkDistrib(void);

  // static void messageMovePatchDone();
  // void movePatchDone();

  static void messageEnqueueWork(Compute *);
  static void messageFinishCUDA(Compute *);
  static void messageFinishMIC(Compute *);
  void enqueueWork(LocalWorkMsg *msg);
  void enqueueExcls(LocalWorkMsg *msg);
  void enqueueBonds(LocalWorkMsg *msg);
  void enqueueAngles(LocalWorkMsg *msg);
  void enqueueDihedrals(LocalWorkMsg *msg);
  void enqueueImpropers(LocalWorkMsg *msg);
  void enqueueThole(LocalWorkMsg *msg);  // Drude model
  void enqueueAniso(LocalWorkMsg *msg);  // Drude model
  void enqueueCrossterms(LocalWorkMsg *msg);
  // JLai
  void enqueueGromacsPair(LocalWorkMsg *msg);
  // End of JLai
  void enqueuePme(LocalWorkMsg *msg);
  void enqueueSelfA1(LocalWorkMsg *msg);
  void enqueueSelfA2(LocalWorkMsg *msg);
  void enqueueSelfA3(LocalWorkMsg *msg);
  void enqueueSelfB1(LocalWorkMsg *msg);
  void enqueueSelfB2(LocalWorkMsg *msg);
  void enqueueSelfB3(LocalWorkMsg *msg);
  void enqueueWorkA1(LocalWorkMsg *msg);
  void enqueueWorkA2(LocalWorkMsg *msg);
  void enqueueWorkA3(LocalWorkMsg *msg);
  void enqueueWorkB1(LocalWorkMsg *msg);
  void enqueueWorkB2(LocalWorkMsg *msg);
  void enqueueWorkB3(LocalWorkMsg *msg);
  void enqueueWorkC(LocalWorkMsg *msg);
  void enqueueCUDA(LocalWorkMsg *msg);
  void enqueueCUDAP2(LocalWorkMsg *msg);
  void enqueueCUDAP3(LocalWorkMsg *msg);
  void finishCUDAPatch(FinishWorkMsg *msg);
  void finishCUDA(LocalWorkMsg *msg);
  void finishCUDAP2(LocalWorkMsg *msg);
  void finishCUDAP3(LocalWorkMsg *msg);
  void enqueueMIC(LocalWorkMsg *msg);
  void finishMIC(LocalWorkMsg *msg);
  void enqueueLCPO(LocalWorkMsg *msg);

  void mapComputes(void);
  void sendPatchMap(void);
  void sendComputeMap(void);
  void saveComputeMapChanges(int,CkGroupID);
  void recvComputeMapChanges(ComputeMapChangeMsg *);
  void doneSaveComputeMap(CkReductionMsg *);

  FullAtomList *createAtomLists(const char *basename=0);
  void createHomePatches(void);
  void distributeHomePatches(void);

  void reinitAtoms(const char *basename=0);
  void patchMapInit(void);
  void assignNodeToPatch(void);

  void savePatchMap(PatchMapMsg *msg);
  void saveComputeMap(ComputeMapMsg *msg);
  inline void setPatchMapArrived(bool s) {patchMapArrived=s;}

#ifdef MEM_OPT_VERSION
  void fillAtomListForOnePatch(int pid, FullAtomList &alist);
  void random_velocities_parallel(BigReal Temp,InputAtomList &inAtoms);
#endif

  static int peOrderingInit;           // used during startup
  static int *peDiffuseOrdering;       // pes in diffuse order
  static int *peDiffuseOrderingIndex;  // index of pe in diffuse order
  static int *peCompactOrdering;       // pes in compact order
  static int *peCompactOrderingIndex;  // index of pe in compact order

  struct pe_sortop_diffuse {
    inline bool operator() (int a, int b) const {
      const int *index = WorkDistrib::peDiffuseOrderingIndex;
      return ( index[a] < index[b] );
    }
  };
  struct pe_sortop_compact {
    inline bool operator() (int a, int b) const {
      const int *index = WorkDistrib::peCompactOrderingIndex;
      return ( index[a] < index[b] );
    }
  };

  static void sortPmePes(int *pmepes, int xdim, int ydim);

  static void peOrderingReady();

  // MIC-Specific
  static void send_initHostDeviceLDB();
  /* entry */ void initHostDeviceLDB();
  static void send_contributeHostDeviceLDB(int peSetLen, int * peSet);
  /* entry */ void contributeHostDeviceLDB(int peSetLen, int * peSet);
  static void send_setDeviceLDBParams(int dt, int hs, int sp1, int pp1, int pp2);
  /* entry */ void setDeviceLDBParams(int dt, int hs, int sp1, int pp1, int pp2);

  static void buildNodeAwarePeOrdering(void);

private:
  void mapComputeNonbonded(void);
  void mapComputeLCPO(void);
  void mapComputeNode(ComputeType);
  void mapComputeHomePatches(ComputeType);
  void mapComputeHomeTuples(ComputeType);
  void mapComputePatch(ComputeType);
  void assignPatchesToLowestLoadNode(void);
  void assignPatchesRecursiveBisection(void);
  void assignPatchesRoundRobin(void);
  void assignPatchesSpaceFillingCurve(void);
  void assignPatchesBitReversal(void);
  int  assignPatchesTopoGridRecBisection();

  void sortNodesAndAssign(int *assignedNode, int baseNodes = 0);
  void velocities_from_PDB(const char *filename, 
			   Vector *v, int totalAtoms);
  void velocities_from_binfile(const char *fname, Vector *vels, int n);
  void random_velocities(BigReal Temp, Molecule *structure,
			 Vector *v, int totalAtoms);
  void remove_com_motion(Vector *vel, Molecule *structure, int n);

  bool patchMapArrived;
  bool computeMapArrived;

  int saveComputeMapReturnEP;
  CkGroupID saveComputeMapReturnChareID;
};

#endif /* WORKDISTRIB_H */

