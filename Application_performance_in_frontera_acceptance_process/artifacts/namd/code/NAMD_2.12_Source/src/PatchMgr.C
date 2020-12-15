/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "PatchMgr.decl.h"
#include "PatchMgr.h"

#include "NamdTypes.h"
//#include "Compute.h"
#include "HomePatch.h"
#include "PatchMap.h"
#include "AtomMap.h"

#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "main.decl.h"
#include "main.h"

#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"
#include "Node.h"
#include "SimParameters.h"

#include "packmsg.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"


// BOC constructor
PatchMgr::PatchMgr()
{
    // CkPrintf("[%d] PatchMgr Created\n", CkMyPe());

    // Singleton pattern
    if (CkpvAccess(PatchMgr_instance) == NULL) {
	CkpvAccess(PatchMgr_instance) = this;
    } else {
	NAMD_bug("PatchMgr instanced twice on same processor!");
    }

    // Get PatchMap singleton started
    patchMap = PatchMap::Instance();
    patchMap->registerPatchMgr(this);

    recvCheckpointReq_index = CmiRegisterHandler((CmiHandler)recvCheckpointReq_handler);
    recvCheckpointLoad_index = CmiRegisterHandler((CmiHandler)recvCheckpointLoad_handler);
    recvCheckpointStore_index = CmiRegisterHandler((CmiHandler)recvCheckpointStore_handler);
    recvCheckpointAck_index = CmiRegisterHandler((CmiHandler)recvCheckpointAck_handler);

    recvExchangeReq_index = CmiRegisterHandler((CmiHandler)recvExchangeReq_handler);
    recvExchangeMsg_index = CmiRegisterHandler((CmiHandler)recvExchangeMsg_handler);

    // Message combining initialization
    migrationCountdown = 0;
    combineMigrationMsgs = new MigrateAtomsCombinedMsg*[CkNumPes()];
    int numPes = CkNumPes();
    for ( int i = 0; i < numPes; ++i ) combineMigrationMsgs[i] = 0;
}

PatchMgr::~PatchMgr()
{
    HomePatchListIter hi(homePatches);
    for ( hi = hi.begin(); hi != hi.end(); hi++) {
      HomePatchElem* elem = homePatches.find(HomePatchElem(hi->pid));
      delete elem->patch;
    }
    delete [] combineMigrationMsgs;
}

void PatchMgr::preCreateHomePatch(PatchID pid, int atomCnt){
    HomePatch *patch = new HomePatch(pid, atomCnt);
    homePatches.load(HomePatchElem(pid, patch));
    patchMap->registerPatch(pid, patch);
}

void PatchMgr::createHomePatch(PatchID pid, FullAtomList &a) 
{
    HomePatch *patch = new HomePatch(pid, a);
    homePatches.load(HomePatchElem(pid, patch));
    patchMap->registerPatch(pid, patch);
}


// Add a HomePatch to a list of patches to be moved 
// HomePatches are actually moved by invoking sendMovePatches() below
void PatchMgr::movePatch(PatchID pid, NodeID nodeID) 
{
    move.load(MovePatch(pid,nodeID));
}

void PatchMgr::sendOneHomePatch(int patchId, int nodeId){
    HomePatch *p = homePatch(patchId);
    patchMap->unregisterPatch(patchId, p);

    MovePatchesMsg *msg = new MovePatchesMsg(patchId, p->atom);

    // Deleting the HomePatchElem will call a destructor for clean up
    // but the msg elements are safe since they use a container template
    // that uses ref counting.
    delete p;
    homePatches.del(HomePatchElem(patchId)); 

    if ( msg->atom.shared() ) NAMD_bug("shared message array in PatchMgr::sendOneHomePatch");

    // Sending to PatchMgr::recvMovePatches on remote node
    CProxy_PatchMgr cp(thisgroup);
    cp[nodeId].recvMovePatches(msg);
}

// Uses list constructed by movePatch() and dispatches
// HomePatch(es) to new nodes
void PatchMgr::sendMovePatches() 
{
    if (! move.size())
	return;

    MovePatchListIter m(move);
    for ( m = m.begin(); m != m.end(); m++) {
      HomePatch *p = homePatch(m->pid);
      patchMap->unregisterPatch(m->pid, p);

      MovePatchesMsg *msg = new MovePatchesMsg(m->pid, p->atom);

      // Deleting the HomePatchElem will call a destructor for clean up
      // but the msg elements are safe since they use a container template
      // that uses ref counting.
      delete p;
      homePatches.del(HomePatchElem(m->pid)); 

      if ( msg->atom.shared() ) NAMD_bug("shared message array in PatchMgr::sendMovePatches");

      // Sending to PatchMgr::recvMovePatches on remote node
      CProxy_PatchMgr cp(thisgroup);
      cp[m->nodeID].recvMovePatches(msg);
    }
    move.resize(0);
}

void PatchMgr::recvMovePatches(MovePatchesMsg *msg) {
    // Make a new HomePatch
    createHomePatch(msg->pid, msg->atom);
    delete msg;

    // Tell sending PatchMgr we received MovePatchMsg
//    AckMovePatchesMsg *ackmsg = 
//      new AckMovePatchesMsg;
//    CSendMsgBranch(PatchMgr,ackMovePatches, ackmsg, thisgroup, msg->fromNodeID);
}
    

//void PatchMgr::ackMovePatches(AckMovePatchesMsg *msg)
//{
//    delete msg;
//    if (! --ackMovePending) 
//	WorkDistrib::messageMovePatchDone();
//}


void PatchMgr::sendAtoms(PatchID pid, FullAtomList &a) {

      MovePatchesMsg *msg = new MovePatchesMsg(pid, a);

      if ( msg->atom.shared() ) NAMD_bug("shared message array in PatchMgr::sendAtoms");

      CProxy_PatchMgr cp(thisgroup);
      cp[patchMap->node(pid)].recvAtoms(msg);

}

void PatchMgr::recvAtoms(MovePatchesMsg *msg) {
    patchMap->homePatch(msg->pid)->reinitAtoms(msg->atom);
    delete msg;
}

// Called by HomePatch to migrate atoms off to new patches
// Message combining occurs here
void PatchMgr::sendMigrationMsgs(PatchID src, MigrationInfo *m, int numMsgs) {
/*
  for (int i=0; i < numMsgs; i++) {
    PatchMgr::Object()->sendMigrationMsg(src, m[i]);
  }
*/
  if ( ! migrationCountdown )  // (re)initialize
  {
    // DebugM(3,"migrationCountdown (re)initialize\n");
    numHomePatches = patchMap->numHomePatches();
    migrationCountdown = numHomePatches;
    combineMigrationDestPes.resize(0);
  }
  for (int i=0; i < numMsgs; i++) {  // buffer messages
    int destNodeID = m[i].destNodeID;
    if ( 1 ) // destNodeID != CkMyPe() )
    {
      if ( ! combineMigrationMsgs[destNodeID] )
      {
        combineMigrationMsgs[destNodeID] = new MigrateAtomsCombinedMsg();
        combineMigrationDestPes.add(destNodeID);
      }
      combineMigrationMsgs[destNodeID]->add(src,m[i].destPatchID,m[i].mList);
    }
    else
    {
	// for now buffer local messages too
    }
  }
  migrationCountdown -= 1;
  // DebugM(3,"migrationCountdown = " << migrationCountdown << "\n");
  if ( ! migrationCountdown )  // send out combined messages
  {
    int n = combineMigrationDestPes.size();
    for ( int i = 0; i < n; ++i ) {
        int destNodeID = combineMigrationDestPes[i];
	DebugM(3,"Sending MigrateAtomsCombinedMsg to node " << destNodeID << "\n");
        CProxy_PatchMgr cp(thisgroup);
        cp[destNodeID].recvMigrateAtomsCombined(combineMigrationMsgs[destNodeID]);
        combineMigrationMsgs[destNodeID] = 0;
    }
  }
}

void PatchMgr::recvMigrateAtomsCombined (MigrateAtomsCombinedMsg *msg)
{
  DebugM(3,"Received MigrateAtomsCombinedMsg with " << msg->srcPatchID.size() << " messages.\n");
  msg->distribute();
  delete msg;
}

void PatchMgr::moveAtom(MoveAtomMsg *msg) {
  LocalID lid = AtomMap::Object()->localID(msg->atomid);
  if ( lid.pid != notUsed ) {
    HomePatch *hp = patchMap->homePatch(lid.pid);
    if ( hp ) {
      FullAtom &a = hp->atom[lid.index];
      if ( msg->moveto ) {
        a.fixedPosition = msg->coord;
      } else {
        a.fixedPosition = hp->lattice.reverse_transform(a.position,a.transform);
        a.fixedPosition += msg->coord;
      }
      a.position = hp->lattice.apply_transform(a.fixedPosition,a.transform);
    }
  }
  delete msg;
}

void PatchMgr::moveAllBy(MoveAllByMsg *msg) {
  // loop over homePatches, moving every atom
  for (HomePatchElem *elem = homePatches.begin(); elem != homePatches.end(); elem++) {
    HomePatch *hp = elem->patch;
    for (int i=0; i<hp->getNumAtoms(); i++) {
      FullAtom &a = hp->atom[i];
      a.fixedPosition = hp->lattice.reverse_transform(a.position,a.transform);
      a.fixedPosition += msg->offset;
      a.position = hp->lattice.apply_transform(a.fixedPosition,a.transform);
    }
  }
  delete msg;
}

void PatchMgr::setLattice(SetLatticeMsg *msg) {
  // loop over homePatches, setting the lattice to the new value.
  for (HomePatchElem *elem = homePatches.begin(); elem != homePatches.end(); elem++) {
    HomePatch *hp = elem->patch;
    hp->lattice = msg->lattice;
  }
  // Must also do this for SimParameters in order for pressure profile to work!
  Node::Object()->simParameters->lattice = msg->lattice;
}


// initiating replica
void PatchMgr::sendCheckpointReq(int pid, int remote, const char *key, int task) {
  CheckpointAtomsReqMsg *msg = new (1+strlen(key),0) CheckpointAtomsReqMsg;
  msg->pid = pid;
  msg->pe = CkMyPe();
  msg->replica = CmiMyPartition();
  msg->task = task;
  strcpy(msg->key,key);
  envelope *env = UsrToEnv(CheckpointAtomsReqMsg::pack(msg));
  CmiSetHandler(env,recvCheckpointReq_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(CkMyPe(),remote,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(CkMyPe(),env->getTotalsize(),(char*)env);
#endif
}

// responding replica
extern "C" {
  void recvCheckpointReq_handler(envelope *env) {
    PatchMgr::Object()->recvCheckpointReq(CheckpointAtomsReqMsg::unpack(EnvToUsr(env)));
  }
}

// responding replica
void PatchMgr::recvCheckpointReq(CheckpointAtomsReqMsg *msg) {
  int patchnode = patchMap->node(msg->pid);
  if ( CkMyPe() != patchnode ) {
    thisProxy[patchnode].recvCheckpointReq(msg);
  } else {
    HomePatch *hp = patchMap->homePatch(msg->pid);
    if ( ! hp ) NAMD_bug("null HomePatch pointer in PatchMgr::recvCheckpointReq");
    hp->recvCheckpointReq(msg->task, msg->key, msg->replica, msg->pe);
    delete msg;
  }
}


// responding replica
void PatchMgr::sendCheckpointLoad(CheckpointAtomsMsg *msg, int dst, int dstpe) {
  envelope *env = UsrToEnv(CheckpointAtomsMsg::pack(msg));
  CmiSetHandler(env,recvCheckpointLoad_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(dstpe,dst,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(dstpe,env->getTotalsize(),(char*)env);
#endif
}

// initiating replica
extern "C" {
  void recvCheckpointLoad_handler(envelope *env) {
    PatchMgr::Object()->recvCheckpointLoad(CheckpointAtomsMsg::unpack(EnvToUsr(env)));
  }
}

// initiating replica
void PatchMgr::recvCheckpointLoad(CheckpointAtomsMsg *msg) {
  HomePatch *hp = patchMap->homePatch(msg->pid);
  hp->recvCheckpointLoad(msg);
}


// initiating replica
void PatchMgr::sendCheckpointStore(CheckpointAtomsMsg *msg, int dst, int dstpe) {
  envelope *env = UsrToEnv(CheckpointAtomsMsg::pack(msg));
  CmiSetHandler(env,recvCheckpointStore_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(dstpe,dst,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(dstpe,env->getTotalsize(),(char*)env);
#endif
}

// responding replica
extern "C" {
  void recvCheckpointStore_handler(envelope *env) {
    PatchMgr::Object()->recvCheckpointStore(CheckpointAtomsMsg::unpack(EnvToUsr(env)));
  }
}

// responding replica
void PatchMgr::recvCheckpointStore(CheckpointAtomsMsg *msg) {
  HomePatch *hp = patchMap->homePatch(msg->pid);
  hp->recvCheckpointStore(msg);
}


// responding replica
void PatchMgr::sendCheckpointAck(int pid, int dst, int dstpe) {
  CheckpointAtomsReqMsg *msg = new CheckpointAtomsReqMsg;
  msg->pid = pid;
  envelope *env = UsrToEnv(CheckpointAtomsReqMsg::pack(msg));
  CmiSetHandler(env,recvCheckpointAck_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(dstpe,dst,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(dstpe,env->getTotalsize(),(char*)env);
#endif
}

// initiating replica
extern "C" {
  void recvCheckpointAck_handler(envelope *env) {
    PatchMgr::Object()->recvCheckpointAck(CheckpointAtomsReqMsg::unpack(EnvToUsr(env)));
  }
}

// initiating replica
void PatchMgr::recvCheckpointAck(CheckpointAtomsReqMsg *msg) {
  HomePatch *hp = patchMap->homePatch(msg->pid);
  if ( ! hp ) NAMD_bug("null HomePatch pointer in PatchMgr::recvCheckpointAck");
  hp->recvCheckpointAck();
  delete msg;
}


void PatchMgr::sendExchangeReq(int pid, int src) {
  ExchangeAtomsReqMsg *msg = new ExchangeAtomsReqMsg;
  msg->pid = pid;
  msg->dstpe = CkMyPe();
  envelope *env = UsrToEnv(ExchangeAtomsReqMsg::pack(msg));
  CmiSetHandler(env,recvExchangeReq_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(CkMyPe(),src,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(CkMyPe(),env->getTotalsize(),(char*)env);
#endif
}

extern "C" {
  void recvExchangeReq_handler(envelope *env) {
    PatchMgr::Object()->recvExchangeReq(ExchangeAtomsReqMsg::unpack(EnvToUsr(env)));
  }
}

void PatchMgr::recvExchangeReq(ExchangeAtomsReqMsg *msg) {
  int patchnode = patchMap->node(msg->pid);
  if ( CkMyPe() != patchnode ) {
    thisProxy[patchnode].recvExchangeReq(msg);
  } else {
    HomePatch *hp = patchMap->homePatch(msg->pid);
    if ( ! hp ) NAMD_bug("null HomePatch pointer in PatchMgr::recvExchangeReq");
    hp->recvExchangeReq(msg->dstpe);
    delete msg;
  }
}

void PatchMgr::sendExchangeMsg(ExchangeAtomsMsg *msg, int dst, int dstpe) {
  envelope *env = UsrToEnv(ExchangeAtomsMsg::pack(msg));
  CmiSetHandler(env,recvExchangeMsg_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(dstpe,dst,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(dstpe,env->getTotalsize(),(char*)env);
#endif
}

extern "C" {
  void recvExchangeMsg_handler(envelope *env) {
    PatchMgr::Object()->recvExchangeMsg(ExchangeAtomsMsg::unpack(EnvToUsr(env)));
  }
}

void PatchMgr::recvExchangeMsg(ExchangeAtomsMsg *msg) {
  HomePatch *hp = patchMap->homePatch(msg->pid);
  hp->recvExchangeMsg(msg);
}

PACK_MSG(MovePatchesMsg,
  PACK(fromNodeID);
  PACK(pid);
  PACK_RESIZE(atom);
)


#include "PatchMgr.def.h"

