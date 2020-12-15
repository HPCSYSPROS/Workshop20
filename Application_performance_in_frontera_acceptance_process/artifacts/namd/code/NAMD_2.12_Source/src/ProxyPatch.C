/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Lattice.h"
#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "main.decl.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "AtomMap.h"
#include "PatchMap.h"
#include "Priorities.h"

#define MIN_DEBUG_LEVEL 2
//#define  DEBUGM
#include "Debug.h"

ProxyPatch::ProxyPatch(PatchID pd) : 
  Patch(pd), proxyMsgBufferStatus(PROXYMSGNOTBUFFERED), 
  curProxyMsg(NULL), prevProxyMsg(NULL)
{
  DebugM(4, "ProxyPatch(" << pd << ") at " << this << "\n");
  ProxyMgr::Object()->registerProxy(patchID);
  numAtoms = -1;
  parent = -1;

#ifndef NODEAWARE_PROXY_SPANNINGTREE
  nChild = 0;
  child = new int[proxySpanDim];
#endif

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  localphs = 0;
#ifdef REMOVE_PROXYRESULTMSG_EXTRACOPY
  int msgstart = sizeof(envelope)+sizeof(ProxyResultVarsizeMsg);
#else
  int msgstart = sizeof(envelope)+sizeof(ProxyResultMsg);
#endif
  localphs = CmiCreatePersistent(PatchMap::Object()->node(patchID), 30000, msgstart);
  ntreephs = 0;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
  treephs = NULL;
#else
  treephs = new PersistentHandle[proxySpanDim];
#endif
#endif

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = -1;
  #endif
  
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    depositLock = CmiCreateLock();
  #endif
}

ProxyPatch::~ProxyPatch()
{
  DebugM(4, "ProxyPatch(" << patchID << ") deleted at " << this << "\n");
  ProxyMgr::Object()->unregisterProxy(patchID);

  // ProxyPatch may be freed because of load balancing if the compute object
  // it corresponds to no longer exist on this specific processor.
  CmiAssert(prevProxyMsg!=NULL);
  if(prevProxyMsg!=NULL) {
// #ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
//       AtomMap::Object()->unregisterIDs(patchID,positionPtrBegin, positionPtrEnd);
// #else
      atomMapper->unregisterIDsCompAtomExt(pExt.begin(),pExt.end());
// #endif      
#if ! CMK_PERSISTENT_COMM || ! USE_PERSISTENT_TREE
      delete prevProxyMsg;
#endif
      prevProxyMsg = NULL;
  }


#ifdef NODEAWARE_PROXY_SPANNINGTREE
  #ifdef USE_NODEPATCHMGR
  delete [] nodeChildren;  
  #endif
#endif
  delete [] child;

  p.resize(0);
  pExt.resize(0);

  lcpoType.resize(0);

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  CmiDestoryPersistent(localphs);
  localphs = 0;
  for (int i=0; i<ntreephs; i++)  CmiDestoryPersistent(treephs[i]);
  delete [] treephs;
#endif
}

void ProxyPatch::boxClosed(int box) {
  ProxyGBISP1ResultMsg *msg1;
  ProxyGBISP2ResultMsg *msg2;
 
  if (box == 1) { // force Box
    // Note: delay the deletion of proxyDataMsg (of the 
    // current step) until the next step. This is done 
    // for the sake of atom migration (ProxyDataMsg) 
    // as the ProxyPatch has to  unregister the atoms 
    // of the previous step in the AtomMap data structure 
    // also denotes end of gbis phase 3
    sendResults();
  } else if ( box == 5) {//end phase 1
    //this msg should only have nonzero atoms if flags.doNonbonded
    int msgAtoms = (flags.doNonbonded) ? numAtoms : 0;
    msg1 = new (msgAtoms,PRIORITY_SIZE) ProxyGBISP1ResultMsg;
    for (int i = 0; i < msgAtoms; i++) {
      msg1->psiSum[i] = psiSum[i];
    }
    msg1->patch = patchID;
    msg1->psiSumLen = msgAtoms;
    msg1->origPe = CkMyPe();
    SET_PRIORITY(msg1,flags.sequence,GB1_PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
    ProxyMgr::Object()->sendResult(msg1);
  } else if ( box == 8) {//end phase 2
    //this msg should only have nonzero atoms if flags.doFullElectrostatics
    int msgAtoms = (flags.doFullElectrostatics) ? numAtoms : 0;
    msg2 = new (msgAtoms,PRIORITY_SIZE) ProxyGBISP2ResultMsg;
    for (int i = 0; i < msgAtoms; i++) {
      msg2->dEdaSum[i] = dEdaSum[i];
    }
    msg2->patch = patchID;
    msg2->dEdaSumLen = msgAtoms;
    msg2->origPe = CkMyPe();
    SET_PRIORITY(msg2,flags.sequence,GB2_PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
    ProxyMgr::Object()->sendResult(msg2);
  } else if (box == 9) {
    //nothing
  } else if (box == 10) {
    // LCPO do nothing
  }


  if ( ! --boxesOpen ) {
    DebugM(2,patchID << ": " << "Checking message buffer.\n");    
    
    if(proxyMsgBufferStatus == PROXYALLMSGBUFFERED) {
          CmiAssert(curProxyMsg != NULL);
          DebugM(3,"Patch " << patchID << " processing buffered proxy ALL data.\n");
          receiveAll(curProxyMsg);          
    }else if(proxyMsgBufferStatus == PROXYDATAMSGBUFFERED) {
          CmiAssert(curProxyMsg != NULL);
          DebugM(3,"Patch " << patchID << " processing buffered proxy data.\n");
          receiveData(curProxyMsg);
    }
  } else {
       DebugM(3,"ProxyPatch " << patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

//each timestep
void ProxyPatch::receiveData(ProxyDataMsg *msg)
{
  DebugM(3, "receiveData(" << patchID << ")\n");

  //delete the ProxyDataMsg of the previous step
  delete prevProxyMsg;
  prevProxyMsg = NULL;

  if ( boxesOpen )
  {
      proxyMsgBufferStatus = PROXYDATAMSGBUFFERED;
    // store message in queue (only need one element, though)
    curProxyMsg = msg;
    return;
  }

  //Reuse position arrays inside proxyDataMsg --Chao Mei
  curProxyMsg = msg;
  prevProxyMsg = curProxyMsg;
  flags = msg->flags;

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  if ( ((int64)msg->positionList) % 32 ) { // not aligned
    p.resize(msg->plLen);
    positionPtrBegin = p.begin();
    memcpy(positionPtrBegin, msg->positionList, sizeof(CompAtom)*(msg->plLen));
  } else { // aligned
    positionPtrBegin = msg->positionList;
  }
  positionPtrEnd = positionPtrBegin + msg->plLen;
  if ( ((int64)positionPtrBegin) % 32 ) NAMD_bug("ProxyPatch::receiveData positionPtrBegin not 32-byte aligned");
#else
  p.resize(msg->plLen);
  memcpy(p.begin(), msg->positionList, sizeof(CompAtom)*(msg->plLen));
#endif

// DMK
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  cudaAtomPtr = msg->cudaAtomList;
#endif
  
  avgPositionPtrBegin = msg->avgPositionList;
  avgPositionPtrEnd = msg->avgPositionList + msg->avgPlLen;
  
  // BEGIN LA
  velocityPtrBegin = msg->velocityList;
  velocityPtrEnd = msg->velocityList + msg->vlLen;
  // END LA

  if ( numAtoms == -1 ) { // for new proxies since receiveAtoms is not called
      //numAtoms = p.size();
      numAtoms = msg->plLen;

      //Retrieve the CompAtomExt list
      CmiAssert(msg->plExtLen!=0);
      pExt.resize(msg->plExtLen);
      memcpy(pExt.begin(), msg->positionExtList, sizeof(CompAtomExt)*(msg->plExtLen));


    // DMK - Atom Separation (water vs. non-water)
    #if NAMD_SeparateWaters != 0
      numWaterAtoms = msg->numWaterAtoms;
    #endif

    positionsReady(1);
  } else {
    positionsReady(0);
  }
}

//every doMigration
void ProxyPatch::receiveAll(ProxyDataMsg *msg)
{
  DebugM(3, "receiveAll(" << patchID << ")\n");

  if ( boxesOpen )
  {
    proxyMsgBufferStatus = PROXYALLMSGBUFFERED;    
    curProxyMsg = msg;
    return;
  }  

  //The prevProxyMsg has to be deleted after this if-statement because
  // positionPtrBegin points to the space inside the prevProxyMsg
  if(prevProxyMsg!=NULL) {
// #ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
//       AtomMap::Object()->unregisterIDs(patchID,positionPtrBegin,positionPtrEnd);
// #else
      atomMapper->unregisterIDsCompAtomExt(pExt.begin(), pExt.end());
// #endif
  }
  //Now delete the ProxyDataMsg of the previous step
#if ! CMK_PERSISTENT_COMM || ! USE_PERSISTENT_TREE
  delete prevProxyMsg;
#endif
  curProxyMsg = msg;
  prevProxyMsg = curProxyMsg;

  flags = msg->flags;

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  if ( ((int64)msg->positionList) % 32 ) { // not aligned
    p.resize(msg->plLen);
    positionPtrBegin = p.begin();
    memcpy(positionPtrBegin, msg->positionList, sizeof(CompAtom)*(msg->plLen));
  } else { // aligned
    positionPtrBegin = msg->positionList;
  }
  positionPtrEnd = positionPtrBegin + msg->plLen;
  if ( ((int64)positionPtrBegin) % 32 ) NAMD_bug("ProxyPatch::receiveAll positionPtrBegin not 32-byte aligned");
#else
  p.resize(msg->plLen);
  memcpy(p.begin(), msg->positionList, sizeof(CompAtom)*(msg->plLen));
#endif

// DMK
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  cudaAtomPtr = msg->cudaAtomList;
#endif

  numAtoms = msg->plLen;
  //numAtoms = p.size();
  
  avgPositionPtrBegin = msg->avgPositionList;
  avgPositionPtrEnd = msg->avgPositionList + msg->avgPlLen;
  
  // BEGIN LA
  velocityPtrBegin = msg->velocityList;
  velocityPtrEnd = msg->velocityList + msg->vlLen;
  // END LA

  if (flags.doGBIS) {
    intRad.resize(numAtoms*2);
    for (int i = 0; i < numAtoms*2;i++) {
      intRad[i] = msg->intRadList[i];
    }
  }

  if (flags.doLCPO) {
    lcpoType.resize(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
      lcpoType[i] = msg->lcpoTypeList[i];
    }
  }

  //We cannot reuse the CompAtomExt list inside the msg because
  //the information is needed at every step. In the current implementation
  //scheme, the ProxyDataMsg msg will be deleted for every step.
  //In order to keep this information, we have to do the extra copy. But
  //this overhead is amortized among the steps that atoms don't migrate
  // --Chao Mei
  pExt.resize(msg->plExtLen);
  memcpy(pExt.begin(), msg->positionExtList, sizeof(CompAtomExt)*(msg->plExtLen));

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = msg->numWaterAtoms;
  #endif

  positionsReady(1);
}

void ProxyPatch::sendResults(void)
{
  DebugM(3, "sendResults(" << patchID << ")\n");
  register int i = 0;
  register ForceList::iterator f_i, f_e, f2_i;
  for ( i = Results::normal + 1 ; i <= flags.maxForceMerged; ++i ) {
    f_i = f[Results::normal].begin(); f_e = f[Results::normal].end();
    f2_i = f[i].begin();
    for ( ; f_i != f_e; ++f_i, ++f2_i ) *f_i += *f2_i;
    f[i].resize(0);
  }
  for ( i = flags.maxForceUsed + 1; i < Results::maxNumForces; ++i )
    f[i].resize(0);

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  CmiUsePersistentHandle(&localphs, 1);
#endif

  if (proxyRecvSpanning == 0) {
#ifdef REMOVE_PROXYRESULTMSG_EXTRACOPY
    ProxyResultVarsizeMsg *msg = ProxyResultVarsizeMsg::getANewMsg(CkMyPe(), patchID, PRIORITY_SIZE, f); 
#else
    ProxyResultMsg *msg = new (PRIORITY_SIZE) ProxyResultMsg;    
    msg->node = CkMyPe();
    msg->patch = patchID;
    for ( i = 0; i < Results::maxNumForces; ++i ) 
      msg->forceList[i] = &(f[i]);
#endif
    SET_PRIORITY(msg,flags.sequence,PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
    //sending results to HomePatch
    ProxyMgr::Object()->sendResults(msg);
  }
  else {
    ProxyCombinedResultMsg *msg = new (PRIORITY_SIZE) ProxyCombinedResultMsg;
    SET_PRIORITY(msg,flags.sequence,
		PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
    msg->nodes.add(CkMyPe());
    msg->patch = patchID;
    for ( i = 0; i < Results::maxNumForces; ++i ) 
      msg->forceList[i] = &(f[i]);
    //sending results to HomePatch
    ProxyMgr::Object()->sendResults(msg);
  }
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  CmiUsePersistentHandle(NULL, 0);
#endif
}

#ifdef NODEAWARE_PROXY_SPANNINGTREE
void ProxyPatch::setSpanningTree(int p, int *c, int n) { 
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE && ! defined(USE_NODEPATCHMGR)
  if (ntreephs!=0) {
      for (int i=0; i<ntreephs; i++)  CmiDestoryPersistent(treephs[i]);
      delete [] treephs;
  }
  treephs = NULL;
  if (n) {
      treephs = new PersistentHandle[n];
      for (int i=0; i<n; i++) {
           treephs[i] = CmiCreatePersistent(c[i], 27000, sizeof(envelope)+sizeof(ProxyDataMsg));
      }
  }
  ntreephs = n;
#endif
  parent=p; nChild = n; nWait = 0;
  delete [] child;
  if(n==0) {
      child = NULL;
      return;
  }
  child = new int[n];
  for (int i=0; i<n; i++) child[i] = c[i];

  #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("ProxyPatch[%d] has %d children: ", patchID, nChild);
    for(int i=0; i<nChild; i++)
        dft->writeTrace("%d ", child[i]);
    dft->writeTrace("\n");
    dft->closeTrace();
  #endif
//CkPrintf("setSpanningTree: [%d:%d] %d %d:%d %d\n", CkMyPe(), patchID, parent, nChild, child[0], child[1]);
}

#ifdef USE_NODEPATCHMGR
void ProxyPatch::setSTNodeChildren(int numNids, int *nids){
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  if (numNodeChild!=0) {
      for (int i=0; i<numNodeChild; i++)  CmiDestoryPersistent(treephs[i]);
      delete [] treephs;
  }
  treephs = NULL;
  if (numNids) {
      treephs = new PersistentHandle[numNids];
      for (int i=0; i<numNids; i++) {
           treephs[i] = CmiCreateNodePersistent(nids[i], 27000, sizeof(envelope)+sizeof(ProxyDataMsg));
      }
  }
  ntreephs = numNids;
#endif
    numNodeChild = numNids;
    delete [] nodeChildren;
    if(numNids==0) {
        nodeChildren = NULL;
        return;
    }
    nodeChildren = new int[numNids];
    for(int i=0; i<numNids; i++) nodeChildren[i] = nids[i]; 
}
#endif

#else //branch for NODEAWARE_PROXY_SPANNINGTREE not defined
void ProxyPatch::setSpanningTree(int p, int *c, int n) { 
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  if (ntreephs!=0) {
      for (int i=0; i<ntreephs; i++)  CmiDestoryPersistent(treephs[i]);
  }
  for (int i=0; i<n; i++) {
       treephs[i] = CmiCreatePersistent(c[i], 27000, sizeof(envelope)+sizeof(ProxyDataMsg));
  }
  ntreephs = n;
#endif
  parent=p; nChild = n; nWait = 0;
  for (int i=0; i<n; i++) child[i] = c[i];
//CkPrintf("setSpanningTree: [%d:%d] %d %d:%d %d\n", CkMyPe(), patchID, parent, nChild, child[0], child[1]);
}
#endif

int ProxyPatch::getSpanningTreeChild(int *c) { 
  for (int i=0; i<nChild; i++) c[i] = child[i];
  return nChild;
}

ProxyCombinedResultMsg *ProxyPatch::depositCombinedResultMsg(ProxyCombinedResultMsg *msg) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CmiLock(depositLock);
#endif
  nWait++;
  if (nWait == 1) msgCBuffer = msg;
  else {
    NodeIDList::iterator n_i, n_e;
    n_i = msg->nodes.begin();
    n_e = msg->nodes.end();
    for (; n_i!=n_e; ++n_i) msgCBuffer->nodes.add(*n_i);
    for ( int k = 0; k < Results::maxNumForces; ++k )
    {
    register ForceList::iterator r_i;
    r_i = msgCBuffer->forceList[k]->begin();
    register ForceList::iterator f_i, f_e;
    f_i = msg->forceList[k]->begin();
    f_e = msg->forceList[k]->end();
    //    for ( ; f_i != f_e; ++f_i, ++r_i ) *r_i += *f_i;

    int nf = f_e - f_i;
#ifdef ARCH_POWERPC
#pragma disjoint (*f_i, *r_i)
#pragma unroll(4)
#endif
    for (int count = 0; count < nf; count++) {
      r_i[count].x += f_i[count].x;      
      r_i[count].y += f_i[count].y;      
      r_i[count].z += f_i[count].z;
    }

    }
    delete msg;
  }
//CkPrintf("[%d:%d] wait: %d of %d (%d %d %d)\n", CkMyPe(), patchID, nWait, nChild+1, parent, child[0],child[1]);
  if (nWait == nChild + 1) {
    nWait = 0;
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    CmiUnlock(depositLock);
#endif
    
    return msgCBuffer;
  }

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CmiUnlock(depositLock);
#endif

  return NULL;
}

//receive data after phase 1 to begin phase 2
void ProxyPatch::receiveData(ProxyGBISP2DataMsg *msg) {
  memcpy(bornRad.begin(), msg->bornRad, sizeof(Real)*numAtoms);
  delete msg;
  Patch::gbisP2Ready();
}

void ProxyPatch::receiveData(ProxyGBISP3DataMsg *msg) {
  memcpy(dHdrPrefix.begin(), msg->dHdrPrefix, sizeof(Real)*msg->dHdrPrefixLen);
  delete msg;
  Patch::gbisP3Ready();
}

ProxyCombinedResultMsg *ProxyPatch::depositCombinedResultRawMsg(ProxyCombinedResultRawMsg *msg) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CmiLock(depositLock);
#endif
  nWait++;
  if (nWait == 1) msgCBuffer = ProxyCombinedResultMsg::fromRaw(msg);
  else {
    for (int i=0; i<msg->nodeSize; i++) msgCBuffer->nodes.add(msg->nodes[i]);

    register char* isNonZero = msg->isForceNonZero;
	register Force* f_i = msg->forceArr;
	for ( int k = 0; k < Results::maxNumForces; ++k )
    {
		register ForceList::iterator r_i;
		r_i = msgCBuffer->forceList[k]->begin();
        int nf = msg->flLen[k];

#ifdef ARCH_POWERPC
#pragma disjoint (*f_i, *r_i)
#endif
		for (int count = 0; count < nf; count++) {
			if(*isNonZero){
				r_i[count].x += f_i->x;
				r_i[count].y += f_i->y;
				r_i[count].z += f_i->z;
				f_i++;
			}
			isNonZero++;
		}
    }
    delete msg;
  }
//CkPrintf("[%d:%d] wait: %d of %d (%d %d %d)\n", CkMyPe(), patchID, nWait, nChild+1, parent, child[0],child[1]);
  if (nWait == nChild + 1) {
    nWait = 0;
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    CmiUnlock(depositLock);
#endif

    return msgCBuffer;
  }

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CmiUnlock(depositLock);
#endif

  return NULL;
}



