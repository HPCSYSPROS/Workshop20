/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef BOCGROUP_H
#define BOCGROUP_H

#include "charm++.h"
#include "ckmulticast.h"
#if     CMK_SMP && USE_CKLOOP
#include "CkLoopAPI.h"
#endif

class BOCgroup {
public:
  CkGroupID workDistrib;
  CkGroupID patchMgr;
  CkGroupID proxyMgr;
  CkGroupID computeMgr;
  CkGroupID computePmeMgr;
  CkGroupID nodePmeMgr;
#ifdef NAMD_CUDA
  CkGroupID computePmeCUDAMgr;
  CkGroupID computeCUDAMgr;
#endif
  //  CkGroupID delegateMgr;
#ifdef OPENATOM_VERSION
  CkGroupID computeMoaMgr;
#endif // OPENATOM_VERSION
  CkGroupID computeExtMgr;
  CkGroupID computeQMMgr;
  CkGroupID computeGBISserMgr;
  CkGroupID computeFmmSerialMgr;
  CkGroupID computeMsmSerialMgr;
  CkGroupID computeMsmMsaMgr;
  CkGroupID computeMsmMgr;
  CkGroupID multicastMgr;  // Charm CkMulticast library module
  CkGroupID reductionMgr;
  CkGroupID collectionMgr;
  CkGroupID broadcastMgr;
  CkGroupID ldbCoordinator;
  CkGroupID sync;
  CkGroupID node;
  CkGroupID ioMgr;
  #ifdef USE_NODEPATCHMGR
  CkGroupID nodeProxyMgr;
  #endif
  
#if     CMK_SMP && USE_CKLOOP
  CProxy_FuncCkLoop ckLoop;
#endif

  CkGroupID dataExchanger;
};

#endif /* BOCGROUP_H */


