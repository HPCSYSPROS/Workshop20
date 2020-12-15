/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/PatchMgr.h,v $
 * $Author: jim $
 * $Date: 2015/03/03 17:54:14 $
 * $Revision: 1.1032 $
 *****************************************************************************/

#ifndef PATCHMGR_H
#define PATCHMGR_H

#include "charm++.h"

#include "NamdTypes.h"
#include "SortedArray.h"
#include "HomePatch.h"
#include "HomePatchList.h"
#include "BOCgroup.h"
#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "PatchMgr.decl.h"

#if USE_TOPOMAP 
#include "TopoManager.h"
#endif

class HomePatch;

class MovePatchesMsg : public CMessage_MovePatchesMsg {
public:
    NodeID  fromNodeID;
    PatchID pid;
    FullAtomList atom;

    MovePatchesMsg(void) { ; }

    MovePatchesMsg(PatchID n, FullAtomList &a) : pid(n)
    {
      atom.swap(a);
      fromNodeID = CkMyPe();
    }

  // pack and unpack functions
  static void* pack(MovePatchesMsg *msg);
  static MovePatchesMsg* unpack(void *ptr);
};

class MoveAtomMsg : public CMessage_MoveAtomMsg {
public:
  int atomid;
  int moveto;
  Vector coord;
};

class MoveAllByMsg : public CMessage_MoveAllByMsg {
public:
  Vector offset;
};

class SetLatticeMsg : public CMessage_SetLatticeMsg {
public:
  Lattice lattice;
};


class CheckpointAtomsReqMsg : public CMessage_CheckpointAtomsReqMsg {
public:
  int task;
  int pid;
  int replica;
  int pe;
  char *key;
};

class CheckpointAtomsMsg : public CMessage_CheckpointAtomsMsg {
public:
  Lattice lattice;
  int task;
  int pid;
  int replica;
  int pe;
  int berendsenPressure_count;
  int numAtoms;
  FullAtom *atoms;
  char *key;
};

extern "C" {
  void recvCheckpointReq_handler(envelope*);
  void recvCheckpointLoad_handler(envelope*);
  void recvCheckpointStore_handler(envelope*);
  void recvCheckpointAck_handler(envelope*);
}


class ExchangeAtomsReqMsg : public CMessage_ExchangeAtomsReqMsg {
public:
  int pid;
  int dstpe;
};

class ExchangeAtomsMsg : public CMessage_ExchangeAtomsMsg {
public:
  Lattice lattice;
  int pid;
  int numAtoms;
  FullAtom *atoms;
};

extern "C" {
  void recvExchangeReq_handler(envelope*);
  void recvExchangeMsg_handler(envelope*);
}

// PatchMgr creates and manages homepatches. There exist one instance of 
// PatchMgr on each node (derived from Charm++ Group).  // That is, when a new operator causes creation of one instance on each node. 
// In addition to creation of homepatches, it handles the atom redistribution
// at the end of each cycle (i.e., atoms can move from patch to patch at the
// cycle boundaries).
struct MovePatch 
{
    MovePatch(PatchID p=-1, NodeID n=-1) : nodeID(n), pid(p) {};
    ~MovePatch() {};

    NodeID nodeID;
    PatchID pid;

    int operator<(MovePatch m) {
      return ( nodeID < m.nodeID );
    }

    int operator==(MovePatch m) {
      return ( nodeID == m.nodeID );
    }
};

typedef SortedArray<MovePatch> MovePatchList;
typedef ResizeArrayIter<MovePatch> MovePatchListIter;

class PatchMgr : public CBase_PatchMgr
{

public:
  PatchMgr();
  ~PatchMgr();

  static PatchMgr* Object() { return CkpvAccess(PatchMgr_instance); }
  
  void createHomePatch(PatchID pid, FullAtomList &a);

  //atomCnt is the number of atoms patch pid has
  void preCreateHomePatch(PatchID pid, int atomCnt);

  void movePatch(PatchID, NodeID);
  void sendMovePatches();
  void recvMovePatches(MovePatchesMsg *msg);

  void sendAtoms(PatchID pid, FullAtomList &a);
  void recvAtoms(MovePatchesMsg *msg);

  // void ackMovePatches(AckMovePatchesMsg *msg);

  HomePatch *homePatch(PatchID pid) {
     return homePatches.find(HomePatchElem(pid))->patch;
  } 

  // void sendMigrationMsg(PatchID, MigrationInfo);
  void sendMigrationMsgs(PatchID, MigrationInfo*, int);
  // void recvMigrateAtoms(MigrateAtomsMsg *);
  void recvMigrateAtomsCombined(MigrateAtomsCombinedMsg *);

  void moveAtom(MoveAtomMsg *msg);
  void moveAllBy(MoveAllByMsg *msg);
  void setLattice(SetLatticeMsg *msg);

  void sendCheckpointReq(int pid, int remote, const char *key, int task);
  void recvCheckpointReq(CheckpointAtomsReqMsg *msg);
  void sendCheckpointLoad(CheckpointAtomsMsg *msg, int dst, int dstpe);
  void recvCheckpointLoad(CheckpointAtomsMsg *msg);
  void sendCheckpointStore(CheckpointAtomsMsg *msg, int dst, int dstpe);
  void recvCheckpointStore(CheckpointAtomsMsg *msg);
  void sendCheckpointAck(int pid, int dst, int dstpe);
  void recvCheckpointAck(CheckpointAtomsReqMsg *msg);

  void sendExchangeReq(int pid, int src);
  void recvExchangeReq(ExchangeAtomsReqMsg *msg);
  void sendExchangeMsg(ExchangeAtomsMsg *msg, int dst, int dstpe);
  void recvExchangeMsg(ExchangeAtomsMsg *msg);

private:
  friend class PatchMap;
  PatchMap *patchMap;

  int numAllPatches;
  int numHomePatches;

  int recvCheckpointReq_index;
  int recvCheckpointLoad_index;
  int recvCheckpointStore_index;
  int recvCheckpointAck_index;

  int recvExchangeReq_index;
  int recvExchangeMsg_index;

  // an array of patch pointers residing on this node
  HomePatchList homePatches;

  // an array of patches to move off this node
  MovePatchList move;
  int ackMovePending;

  // data for combining migration messages
  MigrateAtomsCombinedMsg ** combineMigrationMsgs;
  ResizeArray<int> combineMigrationDestPes;
  int migrationCountdown;

public:
  void setHomePatchFixedAtomNum(int patchId, int numFixed){
      HomePatch *thisHomePatch = patchMap->homePatch(patchId);
      thisHomePatch->setNumFixedAtoms(numFixed);
  }

  void sendOneHomePatch(int patchId, int nodeId);
};



#endif /* PATCHMGR_H */

