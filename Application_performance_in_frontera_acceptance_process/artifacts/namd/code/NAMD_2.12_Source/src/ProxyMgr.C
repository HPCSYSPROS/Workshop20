/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "main.h"
#include "BOCgroup.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "PatchMap.inl"
#include "ProxyPatch.h"
#include "ComputeMap.h"
#include "HomePatch.h"
#include <string.h>
#include "ProcessorPrivate.h"
#include "packmsg.h"
#include "Priorities.h"
#ifndef _NO_ALLOCA_H
#include <alloca.h>
#endif
#ifndef _NO_MALLOC_H
#include <malloc.h>
#endif

#include <map>
#include <vector>
#include <algorithm>

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
#include "qd.h"
#endif

#include "ComputeNonbondedMICKernel.h"
#include "SimParameters.h"
#include "Node.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 2
#include "Debug.h"

#define ALLOCA(TYPE,NAME,SIZE) TYPE *NAME = (TYPE *) alloca((SIZE)*sizeof(TYPE))

int proxySendSpanning	= 0;
int proxyRecvSpanning	= 0;
//"proxySpanDim" is a configuration parameter as "proxyTreeBranchFactor" in configuration file
int proxySpanDim	= 4;
int inNodeProxySpanDim = 16;

PACK_MSG(ProxySpanningTreeMsg,
  PACK(patch);
  PACK(node);
  PACK_RESIZE(tree);
)

void* ProxyResultMsg::pack(ProxyResultMsg *msg) {

  int msg_size = 0;
  msg_size += sizeof(msg->node);
  msg_size += sizeof(msg->patch);

  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j]->size();
    msg_size += sizeof(array_size);
    msg_size += array_size * sizeof(char);    
    msg_size = ALIGN_8 (msg_size);
    Force* f = msg->forceList[j]->begin();
    int nonzero_count = 0;
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) { ++nonzero_count; }
    }
    msg_size += nonzero_count * sizeof(Vector);
  }

  void *msg_buf = CkAllocBuffer(msg,msg_size);
  char *msg_cur = (char *)msg_buf;

  CmiMemcpy((void*)msg_cur,(void*)(&(msg->node)),sizeof(msg->node));
  msg_cur += sizeof(msg->node);
  CmiMemcpy((void*)msg_cur,(void*)(&(msg->patch)),sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j]->size();
    *(int *) msg_cur = array_size;
    msg_cur += sizeof(int);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    msg_cur = (char *)ALIGN_8 (msg_cur);
    Vector *farr = (Vector *)msg_cur;
    Force* f = msg->forceList[j]->begin();

    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) {
        nonzero[i] = 1;
	farr->x = f[i].x;
	farr->y = f[i].y;
	farr->z = f[i].z;
	farr ++;
      } else {
        nonzero[i] = 0;
      }
    }
    msg_cur = (char *) farr;	  
  }

  delete msg;
  return msg_buf;
}

ProxyResultMsg* ProxyResultMsg::unpack(void *ptr) {

  void *vmsg = CkAllocBuffer(ptr,sizeof(ProxyResultMsg));
  ProxyResultMsg *msg = new (vmsg) ProxyResultMsg;
  char *msg_cur = (char*)ptr;

  CmiMemcpy((void*)(&(msg->node)),(void*)msg_cur,sizeof(msg->node));
  msg_cur += sizeof(msg->node);
  CmiMemcpy((void*)(&(msg->patch)),(void*)msg_cur,sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = *(int *) msg_cur;
    msg_cur += sizeof(array_size);
    msg->forceList[j] = &(msg->forceListInternal[j]);
    msg->forceList[j]->resize(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);    
    msg_cur = (char *)ALIGN_8 (msg_cur);
    Vector* farr = (Vector *) msg_cur;
    Force* f = msg->forceList[j]->begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( nonzero[i] ) {
	f[i].x = farr->x;
	f[i].y = farr->y;
	f[i].z = farr->z;
	farr++;
      } else {
        f[i].x = 0.;  f[i].y = 0.;  f[i].z = 0.;
      }
    }    
    msg_cur = (char *) farr;
  }

  CkFreeMsg(ptr);
  return msg;
}

ProxyResultVarsizeMsg *ProxyResultVarsizeMsg::getANewMsg(NodeID nid, PatchID pid, int prioSize, ForceList *fls){

    //1. decide the length of forceArr and iszero field.
    int tmpLen[Results::maxNumForces];
    int iszeroLen = 0;
    for (int i=0; i<Results::maxNumForces; i++){
        tmpLen[i] = fls[i].size();
        iszeroLen += tmpLen[i];
    }
    char *tmpIszero = new char[iszeroLen];
    char *iszeroPtr = tmpIszero;
    int fArrLen = 0;
    for(int i=0; i<Results::maxNumForces; i++) {        
        Force *fiPtr = fls[i].begin();
        for(int j=0; j<tmpLen[i]; j++, fiPtr++, iszeroPtr++) {         
            if(fiPtr->x!=0.0 || fiPtr->y!=0.0 || fiPtr->z!=0) {
                *iszeroPtr=0;
                fArrLen++;
            }else{
                *iszeroPtr=1;
            }            
        }
    }

    //2. Ready to create the msg, and set all fields
    ProxyResultVarsizeMsg *retmsg = new(fArrLen, iszeroLen, prioSize)ProxyResultVarsizeMsg;
    retmsg->node = nid;
    retmsg->patch = pid;
    memcpy(retmsg->flLen, tmpLen, sizeof(int)*Results::maxNumForces);
    iszeroPtr = tmpIszero;
    Force *forcePtr = retmsg->forceArr;
    for(int i=0; i<Results::maxNumForces; i++) {        
        Force *fiPtr = fls[i].begin();
        for(int j=0; j<tmpLen[i]; j++, fiPtr++, iszeroPtr++) {
            if((*iszeroPtr)!=1) {
                forcePtr->x = fiPtr->x;
                forcePtr->y = fiPtr->y;
                forcePtr->z = fiPtr->z;
                forcePtr++;
            }            
        }
    }
    memcpy(retmsg->isZero, tmpIszero, sizeof(char)*iszeroLen);
    delete [] tmpIszero;
    return retmsg;
}

ProxyNodeAwareSpanningTreeMsg *ProxyNodeAwareSpanningTreeMsg::getANewMsg(PatchID pid, NodeID nid, proxyTreeNode *tree, int size){
    int numAllPes = 0;
    for(int i=0; i<size; i++) {
        numAllPes += tree[i].numPes;
    }
    ProxyNodeAwareSpanningTreeMsg *retmsg = new(size, numAllPes, 0) ProxyNodeAwareSpanningTreeMsg;
    retmsg->patch = pid;
    retmsg->procID = nid;
    retmsg->numNodesWithProxies = size;
    int *pAllPes = retmsg->allPes;
    for(int i=0; i<size; i++) {
        retmsg->numPesOfNode[i] = tree[i].numPes;
        for(int j=0; j<tree[i].numPes; j++) {
            *pAllPes = tree[i].peIDs[j];
            pAllPes++;
        }
    }
    return retmsg;
}

//Only available when macro PROCTRACE_DEBUG is defined
void ProxyNodeAwareSpanningTreeMsg::printOut(char *tag){
#ifdef PROCTRACE_DEBUG
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    const char *patchname = "ProxyPatch";
    if(procID == CkMyPe()) patchname = "HomePatch";
    dft->writeTrace("%s: %s[%d] on proc %d node %d has ST (src %d) with %d nodes\n", 
                    tag, patchname, patch, CkMyPe(), CkMyNode(), procID, numNodesWithProxies);
    if(numNodesWithProxies==0) {
        dft->closeTrace();
        return;
    }
    dft->writeTrace("%s: ===%d===pes/node: ", tag, patch);
    for(int i=0; i<numNodesWithProxies; i++) {
        dft->writeTrace("%d ", numPesOfNode[i]);
    }
    dft->writeTrace("\n%s: ===%d===pe list: ", tag, patch);
    int *p = allPes;
    for(int i=0; i<numNodesWithProxies; i++) {
        for(int j=0; j<numPesOfNode[i]; j++) {
            dft->writeTrace("%d ", *p);
            p++;
        }
    }
    dft->writeTrace("\n");    
    dft->closeTrace();
#endif
}

// for spanning tree
ProxyCombinedResultRawMsg* ProxyCombinedResultMsg::toRaw(ProxyCombinedResultMsg *msg) {
  int totalFLLen=0;
  int nonzero_count = 0;
  int nodeSize = msg->nodes.size();
  for (int j = 0; j < Results::maxNumForces; ++j ) {
        int array_size = msg->forceList[j]->size();
    totalFLLen +=  array_size;
    Force* f = msg->forceList[j]->begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) { ++nonzero_count; }
    }
  }

  ProxyCombinedResultRawMsg *msg_buf = new(nodeSize, totalFLLen, nonzero_count, PRIORITY_SIZE)ProxyCombinedResultRawMsg;
  //Copy envelope stuff
  {
         envelope *oenv = UsrToEnv(msg);
         envelope *nenv = UsrToEnv(msg_buf);
         CmiMemcpy(nenv->getPrioPtr(), oenv->getPrioPtr(), nenv->getPrioBytes());
  }

  msg_buf->nodeSize = nodeSize;
  for (int i=0; i<nodeSize; i++) {
    msg_buf->nodes[i] = msg->nodes[i];
  }
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  msg_buf->destPe = msg->destPe;
  #if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
  msg_buf->isFromImmMsgCall = msg->isFromImmMsgCall;
  #endif
  #endif
  msg_buf->patch = msg->patch;

  Force *farr = msg_buf->forceArr;
  char *isNonZeroPtr = msg_buf->isForceNonZero;
  for ( int j = 0; j < Results::maxNumForces; ++j ) {
        int array_size = msg->forceList[j]->size();
    msg_buf->flLen[j] = array_size;
    Force* f = msg->forceList[j]->begin();
    for ( int i = 0; i < array_size; ++i , isNonZeroPtr++) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) {
        *isNonZeroPtr = 1;
                farr->x  =  f[i].x;
        farr->y  =  f[i].y;
        farr->z  =  f[i].z;
        farr ++;
      } else {
        *isNonZeroPtr = 0;
      }
    }
  }
  delete msg;
  return msg_buf;
}

ProxyCombinedResultMsg* ProxyCombinedResultMsg::fromRaw(ProxyCombinedResultRawMsg *ptr) {

  //CkPrintf("[%d]: unpacking: plainData=%p\n", CkMyPe(), ptr->plainData);      

  void *vmsg = CkAllocBuffer(ptr,sizeof(ProxyCombinedResultMsg));
  ProxyCombinedResultMsg *msg = new (vmsg) ProxyCombinedResultMsg;

  for (int i=0; i<ptr->nodeSize; i++) {
    msg->nodes.add(ptr->nodes[i]);
  }
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  msg->destPe = ptr->destPe;
  #if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
  msg->isFromImmMsgCall = ptr->isFromImmMsgCall;
  #endif
  #endif
  msg->patch = ptr->patch;

  char *nonzero = ptr->isForceNonZero;
  Force* farr = ptr->forceArr;

  for ( int j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = ptr->flLen[j];
    msg->forceList[j] = &(msg->forceListInternal[j]);
    msg->forceList[j]->resize(array_size);
    Force* f = msg->forceList[j]->begin();

    for ( int i = 0; i < array_size; ++i, nonzero++ ) {
      if ( *nonzero ) {
                f[i].x = farr->x;
                f[i].y = farr->y;
                f[i].z = farr->z;
                farr++;
      } else {
        f[i].x = 0.;  f[i].y = 0.;  f[i].z = 0.;
      }
    }
  }

  delete ptr;
  return msg;
}

// class static
int ProxyMgr::nodecount = 0;

ProxyMgr::ProxyMgr() { 
  if (CkpvAccess(ProxyMgr_instance)) {
    NAMD_bug("Tried to create ProxyMgr twice.");
  }
  CkpvAccess(ProxyMgr_instance) = this;
}

ProxyMgr::~ProxyMgr() { 
  removeProxies();
  CkpvAccess(ProxyMgr_instance) = NULL;
}


void ProxyMgr::setSendSpanning() {
  if(CkMyRank()!=0) return; 
  proxySendSpanning = 1;
}

int ProxyMgr::getSendSpanning() {
  return proxySendSpanning;
}

void ProxyMgr::setRecvSpanning() {
  if(CkMyRank()!=0) return;
  proxyRecvSpanning = 1;
}

int ProxyMgr::getRecvSpanning() {
  return proxyRecvSpanning;
}

void ProxyMgr::setProxyTreeBranchFactor(int dim){
    if(CkMyRank()!=0) return;
    proxySpanDim = dim;
}

ProxyTree &ProxyMgr::getPtree() {
  return ptree;
}

void ProxyMgr::removeProxies(void)
{
  ProxySetIter pi(proxySet);
  for ( pi = pi.begin(); pi != pi.end(); pi++)
  {
    delete pi->proxyPatch;
  }
  proxySet.clear();
}

void ProxyMgr::removeUnusedProxies(void)
{
  ResizeArray<PatchID> toDelete;
  ProxySetIter pi(proxySet);
  for ( pi = pi.begin(); pi != pi.end(); pi++)
  {
    if ( pi->proxyPatch->getNumComputes() == 0 ) {
      toDelete.add(pi->patchID);
      //fprintf(stderr, "Proxy Deleted Patch %d Proc %d", pi->patchID, CkMyPe());
    }
  }
  PatchID *pidi = toDelete.begin();
  for ( ; pidi != toDelete.end(); ++pidi ) {
    removeProxy(*pidi);
  }
}

// Figure out which proxies we need and create them
void ProxyMgr::createProxies(void)
{
  // Delete the old proxies.
  removeProxies();

  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  int myNode = CkMyPe();
  enum PatchFlag { Unknown, Home, NeedProxy };
  int *patchFlag = new int[numPatches]; 
  int i, j;

  // Note all home patches.
  for ( i = 0; i < numPatches; ++i )
  {
    patchFlag[i] = ( patchMap->node(i) == myNode ) ? Home : Unknown;
  }

  // Add all upstream neighbors.
  PatchID neighbors[PatchMap::MaxOneAway];
  PatchIDList basepids;
  patchMap->basePatchIDList(myNode,basepids);
  for ( i = 0; i < basepids.size(); ++i )
  {
    if ( patchMap->node(basepids[i]) != myNode ) {
	patchFlag[basepids[i]] = NeedProxy;
    }
    int numNeighbors = patchMap->upstreamNeighbors(basepids[i],neighbors);
    for ( j = 0; j < numNeighbors; ++j )
    {
      if ( ! patchFlag[neighbors[j]] ) {
	patchFlag[neighbors[j]] = NeedProxy;
      }
    }
  }

  ComputeMap *computeMap = ComputeMap::Object();

  // Check all patch-based compute objects.
  int nc = computeMap->numComputes();
  for ( i = 0; i < nc; ++i )
  {
#if defined(NAMD_CUDA)
    ComputeType t = computeMap->type(i);
    if ( t == computeNonbondedSelfType || t == computeNonbondedPairType )
      continue;
#elif defined(NAMD_MIC)
    ComputeType t = computeMap->type(i);
    if ( computeMap->directToDevice(i) != 0 ) { continue; } // NOTE: Compute for device will take care of requiring the patch
#endif
    if ( computeMap->node(i) != myNode ) 
      continue;
    int numPid = computeMap->numPids(i);
    for ( j = 0; j < numPid; ++j )
    {
      int pid = computeMap->pid(i,j);
      if ( ! patchFlag[pid] ) {
	patchFlag[pid] = NeedProxy;
      }
    }
  }
  // Create proxy list
  for ( i = 0; i < numPatches; ++i ) {
    if ( patchFlag[i] == NeedProxy )
    { // create proxy patch
      ProxyPatch *proxy = new ProxyPatch(i);
      proxySet.add(ProxyElem(i, proxy));
      patchMap->registerPatch(i, proxy);
    }
  }
  delete[] patchFlag;
}

void
ProxyMgr::createProxy(PatchID pid) {
  Patch *p = PatchMap::Object()->patch(pid);
  if (!p) {
     DebugM(4,"createProxy("<<pid<<")\n");
     ProxyPatch *proxy = new ProxyPatch(pid);
     proxySet.add(ProxyElem(pid,proxy));
     PatchMap::Object()->registerPatch(pid,proxy);
  }
  else {
     DebugM(4,"createProxy("<<pid<<") found " << p->getPatchID() << "\n");
  }
    
}

void
ProxyMgr::removeProxy(PatchID pid) {
  ProxyElem *p = proxySet.find(ProxyElem(pid));
  if (p) { 
    PatchMap::Object()->unregisterPatch(pid,p->proxyPatch);
    delete p->proxyPatch;
    proxySet.del(ProxyElem(pid));
    // iout << iINFO << "Removing unused proxy " << pid << " on " << iPE << ".\n" << endi;
  }
}
  
void
ProxyMgr::registerProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  RegisterProxyMsg *msg = new RegisterProxyMsg;

  msg->node=CkMyPe();
  msg->patch = pid;

  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp[node].recvRegisterProxy(msg);
}

void
ProxyMgr::recvRegisterProxy(RegisterProxyMsg *msg) {
  HomePatch *homePatch = PatchMap::Object()->homePatch(msg->patch);
  homePatch->registerProxy(msg); // message deleted in registerProxy()
}

void
ProxyMgr::unregisterProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  UnregisterProxyMsg *msg = new UnregisterProxyMsg;

  msg->node=CkMyPe();
  msg->patch = pid;

  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp[node].recvUnregisterProxy(msg);
}

void
ProxyMgr::recvUnregisterProxy(UnregisterProxyMsg *msg) {
  HomePatch *homePatch = PatchMap::Object()->homePatch(msg->patch);
  homePatch->unregisterProxy(msg); // message deleted in registerProxy()
}

void 
ProxyMgr::buildProxySpanningTree()
{
  PatchIDList pids;
  if (!CkMyPe()) iout << iINFO << "Building spanning tree ... send: " << proxySendSpanning << " recv: " << proxyRecvSpanning 
      << " with branch factor " << proxySpanDim <<"\n" << endi;
  PatchMap::Object()->homePatchIDList(pids);
  for (int i=0; i<pids.size(); i++) {
    HomePatch *home = PatchMap::Object()->homePatch(pids[i]);
    if (home == NULL) CkPrintf("ERROR: homepatch NULL\n");
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    home->buildNodeAwareSpanningTree();
#else
    home->buildSpanningTree();
#endif
  }
}

void 
ProxyMgr::buildProxySpanningTree2()
{
#if 0
  //The homePatchIDList is an expensive as it goes through
  //every patch ids.
  PatchIDList pids;
  PatchMap::Object()->homePatchIDList(pids);
  for (int i=0; i<pids.size(); i++) {
    HomePatch *home = PatchMap::Object()->homePatch(pids[i]);
    if (home == NULL) CkPrintf("ERROR: homepatch NULL\n");
    home->sendProxies();
  }
#else
  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  HomePatchListIter iter(*hpl);
  for(iter=iter.begin(); iter!=iter.end(); iter++) {
	  HomePatch *home = iter->patch;
	  home->sendProxies();
  }
#endif
}

void 
ProxyMgr::sendProxies(int pid, int *list, int n)
{
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp[0].recvProxies(pid, list, n);
}

//The value defines the max number of intermediate proxies (acting
//as the node to relay proxy msgs to children) allowed to reside 
//on a physical node for proxy spanning tree
#define MAX_INTERNODE 1

//Only for debug
static void outputProxyTree(ProxyTree &ptree, int np){
	FILE *ofp = fopen("patch_proxylist.txt", "w");
	std::vector<int> plist;
	for(int i=0; i<np; i++) {
		fprintf(ofp, "%d: ", i);
		int listlen = ptree.proxylist[i].size();
		fprintf(ofp, "#%d ", listlen);
		plist.clear();
		for(int j=0; j<listlen; j++) {
			plist.push_back(ptree.proxylist[i][j]);
		}
		std::sort(plist.begin(), plist.end());
		for(int j=0; j<listlen; j++) {
			fprintf(ofp, "%d ", plist[j]);
		}
		fprintf(ofp, "\n");
	}
	fclose(ofp);
}

// only on PE 0
void 
ProxyMgr::recvProxies(int pid, int *list, int n)
{
  int nPatches = PatchMap::Object()->numPatches();
  if (ptree.proxylist == NULL)
    ptree.proxylist = new NodeIDList[nPatches];
  ptree.proxylist[pid].resize(n);
  for (int i=0; i<n; i++)
    ptree.proxylist[pid][i] = list[i];
  ptree.proxyMsgCount ++;
  if (ptree.proxyMsgCount == nPatches) {
	//debug
	//outputProxyTree(ptree, nPatches);

    ptree.proxyMsgCount = 0;
    // building and sending of trees is done in two steps now
    // so that the building step can be shifted to the load balancer
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    buildNodeAwareSpanningTree0();
#else
    buildSpanningTree0();    
#endif
    sendSpanningTrees();
  }
}

void ProxyMgr::recvPatchProxyInfo(PatchProxyListMsg *msg){
	int nPatches = PatchMap::Object()->numPatches();
	if(ptree.proxylist == NULL) ptree.proxylist = new NodeIDList[nPatches];
	CmiAssert(msg->numPatches == nPatches);
	int peIdx = 0;
	for(int i=0; i<nPatches; i++) {
		int pid = msg->patchIDs[i];
		int plen = msg->proxyListLen[i];
		ptree.proxylist[pid].resize(plen);
		for(int j=0; j<plen; j++) {
			ptree.proxylist[pid][j] = msg->proxyPEs[peIdx++];
		}		
	}
	delete msg;
	
	//debug
	//outputProxyTree(ptree, nPatches);

	ptree.proxyMsgCount = 0;
    // building and sending of trees is done in two steps now
    // so that the building step can be shifted to the load balancer
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    buildNodeAwareSpanningTree0();
#else
    buildSpanningTree0();
#endif
    sendSpanningTrees();
}

//
// XXX static and global variables are unsafe for shared memory builds.
// The global and static vars should be eliminated.  
// Unfortunately, the routines that use these below are actually 
// in use in NAMD.
//
extern double *cpuloads;
static int *procidx = NULL;
static double averageLoad = 0.0;

static int compLoad(const void *a, const void *b)
{
  int i1 = *(int *)a;
  int i2 = *(int *)b;
  double d1 = cpuloads[i1];
  double d2 = cpuloads[i2];
  if (d1 < d2) 
    return 1;
  else if (d1 == d2) 
    return 0;
  else 
    return -1;
  // sort from high to low
}

static void processCpuLoad()
{
  int i;
  if (!procidx) {
    procidx = new int[CkNumPes()];
  }
  for (i=0; i<CkNumPes(); i++) procidx[i] = i;
  qsort(procidx, CkNumPes(), sizeof(int), compLoad);

  double averageLoad = 0.0;
  for (i=0; i<CkNumPes(); i++) averageLoad += cpuloads[i];
  averageLoad /= CkNumPes();
//  iout << "buildSpanningTree1: no intermediate node on " << procidx[0] << " " << procidx[1] << endi;

}

static int noInterNode(int p)
{
  int exclude = 0;
  if(CkNumPes()<1025)
    exclude = 5;
  else if(CkNumPes()<4097)
    exclude = 10;
  else if(CkNumPes()<8193)
    exclude = 40;
  else if(CkNumPes()<16385)
    exclude = 40;
  else
    exclude = 80;
  for (int i=0; i<exclude; i++) if (procidx[i] == p) return 1;
//  if (cpuloads[p] > averageLoad) return 1;
  return 0;
}

#ifdef NODEAWARE_PROXY_SPANNINGTREE
//only on PE 0
void ProxyMgr::buildNodeAwareSpanningTree0(){
	CkPrintf("Info: build node-aware spanning tree with send: %d, recv: %d with branch factor %d\n", 
			 proxySendSpanning, proxyRecvSpanning, proxySpanDim);
    int numPatches = PatchMap::Object()->numPatches();
    if (ptree.naTrees == NULL) ptree.naTrees = new proxyTreeNodeList[numPatches];
    for (int pid=0; pid<numPatches; pid++)     
        buildSinglePatchNodeAwareSpanningTree(pid, ptree.proxylist[pid], ptree.naTrees[pid]);
       

    //Debug
    //printf("#######################Naive ST#######################\n");
    //printProxySpanningTree();

    //Now the naive spanning tree has been constructed and stored in oneNATree;
    //Afterwards, some optimizations on this naive spanning tree could be done.
    //except the first element as the tree root always contains the processor
    //that has home patch

    //1st Optimization: reduce intermediate nodes as much as possible. In details,
    //the optimal case is that on a single physical smp node, there should be no
    //two proxies who act as the intermediate nodes to pass information to childrens
    //in the spanning tree. E.g, for patch A's proxy spanning tree, it has a node X as
    //its intermediate node. However, for patch B's, it also has a node X as its intermediate
    //node. We should avoid this situation as node X becomes the bottleneck as it has twice
    //amount of work to process now.
    //Step1: foward to the first patch that has proxies
    //Now proxyNodeMap records the info that how many intermediate nodes on a node
        //each element indiates the number of proxies residing on this node    
    int pid=0;
    for(;pid<numPatches; pid++) {
        if(ptree.proxylist[pid].size()>0) break;
    }
    if(pid==numPatches) {
        return;
    }
    int *proxyNodeMap = new int[CkNumNodes()];
    memset(proxyNodeMap, 0, sizeof(int)*CkNumNodes());
   {
    proxyTreeNodeList &onePatchT = ptree.naTrees[pid];
    //If a node is an intermediate node, then its idx should satisfy
    //idx*proxySpanDim + 1 < onePatchT.size()
    int lastInterNodeIdx = (onePatchT.size()-2)/proxySpanDim;
    for(int i=1; i<lastInterNodeIdx; i++) { //excluding the root node
        int nid = onePatchT.item(i).nodeID;
        proxyNodeMap[nid]++;
    }
   }
    //Step2: iterate over each patch's proxy spanning tree to adjust
    //the tree node positions. The bad thing here is that it may involve
    //many memory allocations and deallocation for small-size (~100bytes)
    //chunks.
    pid++; //advance to the next patch
    for(; pid<numPatches; pid++) {
        if(ptree.proxylist[pid].size()==0) continue;
        proxyTreeNodeList &onePatchT = ptree.naTrees[pid];
        int lastInterNodeIdx = (onePatchT.size()-2)/proxySpanDim;
        for(int i=1; i<=lastInterNodeIdx; i++) {
            int nid = onePatchT.item(i).nodeID;
            if(proxyNodeMap[nid]<MAX_INTERNODE) {
                proxyNodeMap[nid]++;
                continue;
            }
            //the position is occupied, so search the children
            //nodes to see whether there's one to swap this node
            //if not found, find the first position that has smallest
            //amount of nodes.
            int leastIdx = -1;
            int leastAmount = ~(1<<31);
            //iterate children nodes
            int swapPos;
            for(swapPos=lastInterNodeIdx+1; swapPos<onePatchT.size(); swapPos++) {
                int chiNId = onePatchT.item(swapPos).nodeID;
                if(proxyNodeMap[chiNId]<MAX_INTERNODE) {
                    break;
                }
                if(proxyNodeMap[chiNId]<leastAmount) {
                    leastAmount = proxyNodeMap[chiNId];
                    leastIdx = swapPos;
                }
            }
            if(swapPos==onePatchT.size()) {
                CmiAssert(leastIdx!=-1); //because the above loop at least executes once
                //indicate we cannot find a physical node which
                //still allows the intermediate proxy.
                swapPos = leastIdx;
            }
            //swap the current proxy tree node "i" with node "swapPos"
            proxyTreeNode *curNode = &onePatchT.item(i);
            proxyTreeNode *swapNode = &onePatchT.item(swapPos);
            proxyNodeMap[swapNode->nodeID]++; //update the proxyNodeMap record
            int tmp = curNode->nodeID;
            curNode->nodeID = swapNode->nodeID;
            swapNode->nodeID = tmp;
            tmp = curNode->numPes;
            curNode->numPes = swapNode->numPes;
            swapNode->numPes = tmp;
            int *tmpPes = curNode->peIDs;
            curNode->peIDs = swapNode->peIDs;
            swapNode->peIDs = tmpPes;
        }
    }
    delete [] proxyNodeMap;    

    //Debug
    //printf("#######################After 1st optimization#######################\n");
    //printProxySpanningTree();

    //2nd optimization: similar to the 1st optimization but now thinking in
    //the core level. If we cannot avoid place two intermediate proxy
    //on the same node, we'd better to place them in different cores inside
    //the node
    if(CmiMyNodeSize()==1) {
        //No need to perform the second optimization as every node has only 1 core
        return;
    }
    //Step1: forward to the first patch that has proxies
    pid=0;
    for(;pid<numPatches; pid++) {
        if(ptree.proxylist[pid].size()>0) break;
    }
    if(pid==numPatches) {
        return;
    }
    int *proxyCoreMap = new int[CkNumPes()];
    memset(proxyCoreMap, 0, sizeof(int)*CkNumPes());
   {
    proxyTreeNodeList &onePatchT = ptree.naTrees[pid];
    //If a node is an intermediate node, then its idx should satisfy
    //idx*proxySpanDim + 1 < onePatchT.size()
    int lastInterNodeIdx = (onePatchT.size()-2)/proxySpanDim;
    for(int i=1; i<lastInterNodeIdx; i++) { //excluding the root node
        int rootProcID = onePatchT.item(i).peIDs[0];
        proxyCoreMap[rootProcID]++;
    }
   }
    //Step2: iterate over each patch's proxy spanning tree to adjust
    //the root's position of intermediate proxies.
    pid++; //advance to the next patch
    for(; pid<numPatches; pid++) {
        if(ptree.proxylist[pid].size()==0) continue;
        proxyTreeNodeList &onePatchT = ptree.naTrees[pid];
        int lastInterNodeIdx = (onePatchT.size()-2)/proxySpanDim;
        for(int i=1; i<=lastInterNodeIdx; i++) {
            proxyTreeNode *curNode = &onePatchT.item(i);
            int rootProcID = curNode->peIDs[0];
            if(curNode->numPes==1 || proxyCoreMap[rootProcID]<MAX_INTERNODE){
                //if this node contains only 1 core, then we have to leave it as it is
                //because there are no other cores in the same node that could be used to
                //adjust
                proxyCoreMap[rootProcID]++;
                continue;
            }
            
            //foound more than MAX_INTERNODE intermediate proxies on the same core,
            //adjust the root id of the core of this proxy tree node
            int leastIdx = -1;
            int leastAmount = ~(1<<31);
            //iterate children nodes
            int swapPos;
            
            for(swapPos=1; swapPos<curNode->numPes; swapPos++) {
                int otherCoreID = curNode->peIDs[swapPos];
                if(proxyCoreMap[otherCoreID]<MAX_INTERNODE) {
                    break;
                }
                if(proxyCoreMap[otherCoreID]<leastAmount) {
                    leastAmount = proxyCoreMap[otherCoreID];
                    leastIdx = swapPos;
                }
            }
            if(swapPos==curNode->numPes) {
	        CmiAssert(leastIdx!=-1); //because the above loop body must execute at least once
                //indicate we cannot find a physical node which
                //still allows the intermediate proxy.
                swapPos = leastIdx;
            }
            int tmp = curNode->peIDs[swapPos];
            curNode->peIDs[swapPos] = curNode->peIDs[0];
            curNode->peIDs[0] = tmp;
            proxyCoreMap[tmp]++;
        }      
    }

    delete proxyCoreMap;

    //Debug
    //printf("#######################After 2nd optimization#######################\n");
    //printProxySpanningTree();
}

void ProxyMgr::buildSinglePatchNodeAwareSpanningTree(PatchID pid, NodeIDList &proxyList, 
                                                     proxyTreeNodeList &ptnTree){       
    int numProxies = proxyList.size();
    if (numProxies == 0) {
        //CkPrintf ("This is sheer evil in building node-aware spanning tree!\n\n");            
        return;
    }        
 
    //usually the #proxies is at most 62 (considering 2-away in all dimensions)
    //so the access in proxyNodeMap and proxyTreeIdx is at most log2(62)=8 if
    //all proxies are in different nodes
    //could be better than a CkNumNodes() array in that cache perf. is better
    //because of the reduced memory footprint -Chao Mei
    std::map<int, int> proxyNodeMap; //<node id, numProxies>    
    std::vector<int> proxyNodeIDs;
    std::map<int, int> proxyTreeIdx; //<node id, idx in proxyNodeIDs>
    
    //the processor id of home patch
    int hpProcID = PatchMap::Object()->node(pid);
    int hpNodeID = CkNodeOf(hpProcID);
    proxyNodeMap[hpNodeID]=1;
    proxyTreeIdx[hpNodeID]=0;
    proxyNodeIDs.push_back(hpNodeID);
    //proxyNodeList[0] = hpNodeID;
    //int numNodesWithProxies = 1;
    
    for(int i=0; i<numProxies; i++) {
        int procId = proxyList[i];
        int nodeId = CkNodeOf(procId);
        std::map<int, int>::iterator it=proxyNodeMap.find(nodeId);
        if(it==proxyNodeMap.end()) {
            proxyNodeMap[nodeId] = 1;
            proxyTreeIdx[nodeId] = proxyNodeIDs.size();
            proxyNodeIDs.push_back(nodeId);
        }else{
            proxyNodeMap[nodeId]++;
        }        
    }
    proxyTreeNodeList &oneNATree = ptnTree;   // spanning tree
    int numNodesWithProxies = proxyNodeIDs.size();
    oneNATree.resize(numNodesWithProxies);
    //initialize oneNATree
    for(int i=0; i<numNodesWithProxies; i++) {
        proxyTreeNode *oneNode = &oneNATree.item(i);
        delete [] oneNode->peIDs;
        oneNode->nodeID = proxyNodeIDs[i];
        oneNode->peIDs = new int[proxyNodeMap[oneNode->nodeID]];                        
        oneNode->numPes = 0; //initially set to zero as used for incrementing later
    }
    
    //set up the tree root which contains the home patch processor
    proxyTreeNode *rootnode = &oneNATree.item(0);
    rootnode->peIDs[0] = hpProcID;
    rootnode->numPes++;
    
    for(int i=0; i<numProxies; i++) {
        int procId = proxyList[i];
        int nodeId = CkNodeOf(procId);
        int idxInTree = proxyTreeIdx[nodeId];
        CmiAssert(idxInTree>=0 && idxInTree<numNodesWithProxies);
        proxyTreeNode *oneNode = &oneNATree.item(idxInTree);
        oneNode->peIDs[oneNode->numPes] = procId;
        oneNode->numPes++;
    }
}
#else //branch of NODEAWARE_PROXY_SPANNINGTREE
// only on PE 0
void 
ProxyMgr::buildSpanningTree0()
{
	CkPrintf("Info: build spanning tree with send: %d, recv: %d with branch factor %d\n", 
			 proxySendSpanning, proxyRecvSpanning, proxySpanDim);

  int i;

  processCpuLoad();

  int *numPatchesOnNode = new int[CkNumPes()];
  int numNodesWithPatches = 0;
  for (i=0; i<CkNumPes(); i++) numPatchesOnNode[i] = 0;
  int numPatches = PatchMap::Object()->numPatches();
  for (i=0; i<numPatches; i++) {
    int node = PatchMap::Object()->node(i);
    numPatchesOnNode[node]++;
    if (numPatchesOnNode[node] == 1)
      numNodesWithPatches ++;
  }
  int patchNodesLast =
    ( numNodesWithPatches < ( 0.7 * CkNumPes() ) );
  int *ntrees = new int[CkNumPes()];
  for (i=0; i<CkNumPes(); i++) ntrees[i] = 0;
  if (ptree.trees == NULL) ptree.trees = new NodeIDList[numPatches];
  for (int pid=0; pid<numPatches; pid++) 
  {
    int numProxies = ptree.proxylist[pid].size();
    if (numProxies == 0) {
      //CkPrintf ("This is sheer evil!\n\n");
      //ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, NULL, 0);
      delete [] ntrees;
      delete [] numPatchesOnNode;
      return;
    }
    NodeIDList &tree = ptree.trees[pid];   // spanning tree
    NodeIDList oldtree;  oldtree.swap(tree);
    tree.resize(numProxies+1);
    tree.setall(-1);
    tree[0] = PatchMap::Object()->node(pid);
    int s=1, e=numProxies;
    int nNonPatch = 0;
    int treesize = 1;
    int pp;

    // keep tree persistent for non-intermediate nodes
    for (pp=0; pp<numProxies; pp++) {
      int p = ptree.proxylist[pid][pp];
      int oldindex = oldtree.find(p);
      if (oldindex != -1 && oldindex <= numProxies) {
        int isIntermediate = (oldindex*proxySpanDim+1 <= numProxies);
        if (!isIntermediate) {
          tree[oldindex] = p;
        }
        else if (ntrees[p] < MAX_INTERNODE) {
          tree[oldindex] = p;
          ntrees[p] ++;
        }
      }
    }

    for (pp=0; pp<numProxies; pp++) {
      int p = ptree.proxylist[pid][pp];              // processor number
      if (tree.find(p) != -1) continue;        // already used
      treesize++;
      if (patchNodesLast && numPatchesOnNode[p] ) {
        while (tree[e] != -1) { e--; if (e==-1) e = numProxies; }
        tree[e] = p;
        int isIntermediate = (e*proxySpanDim+1 <= numProxies);
        if (isIntermediate) ntrees[p]++;
      }
      else {
        while (tree[s] != -1) { s++; if (s==numProxies+1) s = 1; }
        int isIntermediate = (s*proxySpanDim+1 <= numProxies);
        if (isIntermediate && (ntrees[p] >= MAX_INTERNODE || noInterNode(p))) {   // TOO MANY INTERMEDIATE TREES
        //if (isIntermediate && ntrees[p] >= MAX_INTERNODE)    // TOO MANY INTERMEDIATE TREES
          while (tree[e] != -1) { e--; if (e==-1) e = numProxies; }
          tree[e] = p;
          isIntermediate = (e*proxySpanDim+1 <= numProxies);
          if (isIntermediate) ntrees[p]++;
        }
        else {
          tree[s] = p;
          nNonPatch++;
          if (isIntermediate) ntrees[p]++;
        }
      }
    }
    // send homepatch's proxy tree
    if(ptree.sizes)
      ptree.sizes[pid] = treesize;
    //ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, &tree[0], treesize);
  }
  /*for (i=0; i<CkNumPes(); i++) {
    if (ntrees[i] > MAX_INTERNODE) iout << "Processor " << i << "has (guess) " << ntrees[i] << " intermediate nodes." << endi;
  }*/
  delete [] ntrees;
  delete [] numPatchesOnNode;
}
#endif

void ProxyMgr::sendSpanningTrees()
{
  int numPatches = PatchMap::Object()->numPatches();
  for (int pid=0; pid<numPatches; pid++) {
    int numProxies = ptree.proxylist[pid].size();
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    if (numProxies == 0)
      ProxyMgr::Object()->sendNodeAwareSpanningTreeToHomePatch(pid, NULL, 0);
    else {
      ProxyMgr::Object()->sendNodeAwareSpanningTreeToHomePatch(pid, ptree.naTrees[pid].begin(), ptree.naTrees[pid].size());
    }
#else
    if (numProxies == 0)
      ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, NULL, 0);
    else {
      ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, ptree.trees[pid].begin(), ptree.trees[pid].size());
    }
#endif
  }
}

void ProxyMgr::sendSpanningTreeToHomePatch(int pid, int *tree, int n)
{
  CProxy_ProxyMgr cp(thisgroup);
  cp[PatchMap::Object()->node(pid)].recvSpanningTreeOnHomePatch(pid, tree, n);
}

void ProxyMgr::recvSpanningTreeOnHomePatch(int pid, int *tree, int n)
{
  HomePatch *p = PatchMap::Object()->homePatch(pid);
  p->recvSpanningTree(tree, n);
}

void ProxyMgr::sendNodeAwareSpanningTreeToHomePatch(int pid, proxyTreeNode *tree, int n)
{
  CProxy_ProxyMgr cp(thisgroup);
  ProxyNodeAwareSpanningTreeMsg *msg = ProxyNodeAwareSpanningTreeMsg::getANewMsg(pid, -1, tree, n);
  cp[PatchMap::Object()->node(pid)].recvNodeAwareSpanningTreeOnHomePatch(msg);
}

void ProxyMgr::recvNodeAwareSpanningTreeOnHomePatch(ProxyNodeAwareSpanningTreeMsg *msg)
{
  HomePatch *p = PatchMap::Object()->homePatch(msg->patch);
  p->recvNodeAwareSpanningTree(msg);
  delete msg;
}

void 
ProxyMgr::sendSpanningTree(ProxySpanningTreeMsg *msg) {
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp[msg->tree[0]].recvSpanningTree(msg);
}

void ProxyMgr::sendNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg){
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  int pe = msg->allPes[0]; //the root procID

#if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
  DebugFileTrace *dft = DebugFileTrace::Object();
  dft->openTrace();
  dft->writeTrace("PMgr::sndST: from proc %d for patch[%d]\n", pe, msg->patch);
  dft->closeTrace();
#endif

  cp[pe].recvNodeAwareSpanningTree(msg);
}

void 
ProxyMgr::recvSpanningTree(ProxySpanningTreeMsg *msg) {
  int size = msg->tree.size();
  int *child = new int[proxySpanDim];
  int nChild = 0;
  int i;
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  for (i=0; i<proxySpanDim; i++) {
    if (size > i+1) { child[i] = msg->tree[i+1]; nChild++; }
  }
  if (!PatchMap::Object()->homePatch(msg->patch)) {
    proxy->setSpanningTree(msg->node, child, nChild);
  }

  // build subtree and pass down
 if (nChild > 0) {

  nodecount ++;
  //if (nodecount > MAX_INTERNODE) 
  //  iout << "Processor " << CkMyPe() << "has (actual) " << nodecount << " intermediate nodes." << endi;

//CkPrintf("[%d] %d:(%d) %d %d %d %d %d\n", CkMyPe(), msg->patch, size, msg->tree[0], msg->tree[1], msg->tree[2], msg->tree[3], msg->tree[4]);
  NodeIDList *tree = new NodeIDList[proxySpanDim];
  int level = 1, index=1;
  int done = 0;
  while (!done) {
    for (int n=0; n<nChild; n++) {
      if (done) break;
      for (int j=0; j<level; j++) {
       if (index >= size) { done = 1; break; }
       tree[n].add(msg->tree[index]);
       index++;
      }
    }
    level *=proxySpanDim;
  }

  ProxyMgr *proxyMgr = ProxyMgr::Object();
  for (i=0; i<proxySpanDim; i++) {
    if (tree[i].size()) {
      ProxySpanningTreeMsg *cmsg = new ProxySpanningTreeMsg;
      cmsg->patch = msg->patch;
      cmsg->node = CkMyPe();
      cmsg->tree.copy(tree[i]);  // copy data for thread safety
      proxyMgr->sendSpanningTree(cmsg);
    }
  }

  delete [] tree;
 }

  delete [] child;
  delete msg;
}

//The "msg" represents the subtree rooted at this proc
void ProxyMgr::recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg){
#if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("PMgr::recvST0 for patch[%d] with #nodes=%d\n", msg->patch, msg->numNodesWithProxies);
    dft->closeTrace();
    msg->printOut("PMgr::recvST");
#endif

    //This function is divided into three parts. The tree root is msg->allPes[0]
    //1. set up its own immediate childrens
    int treesize = msg->numNodesWithProxies; //at least include one as its internal procs    
    int iNChild = msg->numPesOfNode[0]-1; //number of internal children
    int eNChild = treesize-1; //number of external children

    CmiAssert(treesize>0);
    //use the same way of computing children in HomePatch::setupChildrenFromProxySpanningTree    
    eNChild = (proxySpanDim>eNChild)?eNChild:proxySpanDim;
    int iSlots = proxySpanDim-eNChild;
    if(iNChild>0) {
        if(iSlots==0){
            //at least having one internal child
            iNChild = 1;    
        }else{
            iNChild = (iSlots>iNChild)?iNChild:iSlots;
        }
    }    
    int numChild = iNChild + eNChild;
    if(numChild==0){
        //Indicating this proxy is a leaf in the spanning tree
        ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
        proxy->setSpanningTree(msg->procID, NULL, 0);
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
		//When using NODEPATCHMGR, the proc-level is a flat list attached to the node
		//while the node-level spanning tree obeys the branch factor.
		//As a result, when passing down spanning trees, if this proc is on the same node
		//of its parent, then the NodeProxyMgr has already been set by its parent. There's
		//no need resetting here. However, the nodeChildren attached to this proxy has
		//to be set to NULL. -Chao Mei
		int onSameNode = (CkMyNode() == CkNodeOf(msg->procID));
        //set up proxyInfo inside NodeProxyMgr
        if(!onSameNode && !PatchMap::Object()->homePatch(msg->patch)){
            //only when this processor contains a proxy patch of "msg->patch"
            //is the patch registeration in NodeProxyMgr needed,
            //and itself needs to be registered
            CProxy_NodeProxyMgr pm(CkpvAccess(BOCclass_group).nodeProxyMgr);
            NodeProxyMgr *npm = pm[CkMyNode()].ckLocalBranch();        
            npm->registerPatch(msg->patch, msg->numPesOfNode[0], msg->allPes);            
        }
        //set children in terms of node ids
        proxy->setSTNodeChildren(0, NULL);       
#endif
        delete msg;
        return;
    }

    nodecount++;
    //if (nodecount > MAX_INTERNODE) 
    //  iout << "Processor " << CkMyPe() << "has (actual) " << nodecount << " intermediate nodes." << endi;

    if(!PatchMap::Object()->homePatch(msg->patch)){
        //the home patch of this spanning tree has been already set
        //in HomePatch::setupChildrenFromProxySpanningTree, so we
        //only care about the children setup for proxy patches here
        ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
        ALLOCA(int,children,numChild);
        //add external children
        int *p = msg->allPes + msg->numPesOfNode[0];
        for(int i=0; i<eNChild; i++) {
            children[i] = *p;
            p += msg->numPesOfNode[i+1];
        }
        //add internal children
        for(int i=eNChild, j=1; i<numChild; i++, j++) {
            children[i] = msg->allPes[j]; 
        }
        proxy->setSpanningTree(msg->procID, children, numChild);

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
		int onSameNode = (CkMyNode() == CkNodeOf(msg->procID));
		if(!onSameNode) {
			//set up proxyInfo inside NodeProxyMgr
			CProxy_NodeProxyMgr pm(CkpvAccess(BOCclass_group).nodeProxyMgr);
			NodeProxyMgr *npm = pm[CkMyNode()].ckLocalBranch();        
			npm->registerPatch(msg->patch, msg->numPesOfNode[0], msg->allPes);        
	
			//set children in terms of node ids
			ALLOCA(int,nodeChildren,eNChild+1);
			p = msg->allPes + msg->numPesOfNode[0];
			for(int i=0; i<eNChild; i++) {
				nodeChildren[i] = CkNodeOf(*p);
				p += msg->numPesOfNode[i+1];
			}
			//the last entry always stores the node id that contains this proxy
			nodeChildren[eNChild] = CkNodeOf(msg->allPes[0]);
			proxy->setSTNodeChildren(eNChild+1, nodeChildren);
		} else {
			proxy->setSTNodeChildren(0, NULL);
		}
#endif
    }

    //2. send msgs for the tree to children proxies
    ResizeArray<int> *exTreeChildSize = new ResizeArray<int>[numChild];
    ResizeArray<int *> *exTreeChildPtr = new ResizeArray<int *>[numChild];    

    //2a. first processing children of external nodes
    if(eNChild > 0) {    
        int nodesToCnt = 1; //the number of children each root (current root's 
                            //immedidate external nodes) has in each level
        int pos = 1; //track the iteration over msg->numPesOfNode and skip the current root
        int *pePtr = msg->allPes + msg->numPesOfNode[0];
        int done = 0;
        while(!done) {
            for(int childID=0; childID<eNChild; childID++) {
                //iterate nodes on each level
                for(int i=0; i<nodesToCnt; i++) {
                    int cursize = msg->numPesOfNode[pos];
                    exTreeChildSize[childID].add(cursize);
                    exTreeChildPtr[childID].add(pePtr);
                    pos++;
                    pePtr += cursize; 
                    if(pos==msg->numNodesWithProxies) {
                        done = 1;
                        break;
                    }
                }
                if(done) break;                         
            }
            nodesToCnt *= proxySpanDim;
        }
    }

    //2b. secondly processing children on the same node
    if(iNChild>0) {
        int nodesToCnt = 1; //the number of children each root (current root's 
                            //immedidate internal child proxies) has in each level
        int pos = 1; //track the iteration over proxies on the same node and skip the current root
        int *pePtr = msg->allPes+1; //skip the root
        int done = 0;
        while(!done) {
            for(int childID=eNChild; childID<numChild; childID++) {
                //iterate nodes on each level
                for(int i=0; i<nodesToCnt; i++) {                    
                    exTreeChildSize[childID].add(1);
                    exTreeChildPtr[childID].add(pePtr);
                    pos++;
                    pePtr++; 
                    if(pos==msg->numPesOfNode[0]) {
                        done = 1;
                        break;
                    }
                }
                if(done) break;                         
            }
            nodesToCnt *= proxySpanDim;
        }
    }
          
    for(int i=0; i<numChild; i++) {                
        ResizeArray<int> *allSizes = &exTreeChildSize[i];
        ResizeArray<int *> *allPtrs = &exTreeChildPtr[i];
        int totalNodes = allSizes->size();
        int totalPes = 0;
        for(int j=0; j<totalNodes; j++) totalPes += allSizes->item(j);
        ProxyNodeAwareSpanningTreeMsg *cmsg = new(totalNodes, totalPes, 0) ProxyNodeAwareSpanningTreeMsg;
        cmsg->patch = msg->patch;
        cmsg->procID = CkMyPe();
        cmsg->numNodesWithProxies = totalNodes;
        int *pAllPes = cmsg->allPes;
        for(int j=0; j<totalNodes; j++) {
            int numPes = allSizes->item(j);
            cmsg->numPesOfNode[j] = numPes;
            memcpy(pAllPes, allPtrs->item(j), sizeof(int)*numPes);
            pAllPes += numPes;
        }
        #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
        cmsg->printOut("sndChi:");
        #endif
        ProxyMgr::Object()->sendNodeAwareSpanningTree(cmsg);
    }
    delete [] exTreeChildSize;
    delete [] exTreeChildPtr;  
    delete msg;
}

void ProxyMgr::recvNodeAwareSTParent(int patch, int parent){
#if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("PMgr::recvSTParent: for ProxyPatch[%d], parent is %d\n", patch, parent);
    dft->closeTrace();
#endif
    ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(patch);
    CmiAssert(proxy!=NULL);
    proxy->setSpanningTree(parent, NULL, 0);
}

void ProxyMgr::sendResults(ProxyResultVarsizeMsg *msg) {
    CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
    NodeID node = PatchMap::Object()->node(msg->patch);
    CmiEnableUrgentSend(1);
    cp[node].recvResults(msg);
    CmiEnableUrgentSend(0);
}

void ProxyMgr::recvResults(ProxyResultVarsizeMsg *msg) {
    HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
    home->receiveResults(msg); // delete done in HomePatch::receiveResults()
}

void ProxyMgr::sendResults(ProxyResultMsg *msg) {
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  NodeID node = PatchMap::Object()->node(msg->patch);
  CmiEnableUrgentSend(1);
  cp[node].recvResults(msg);
  CmiEnableUrgentSend(0);
}

void ProxyMgr::recvResults(ProxyResultMsg *msg) {
  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  home->receiveResults(msg); // delete done in HomePatch::receiveResults()
}

//sendResults is a direct function call, not an entry method
void ProxyMgr::sendResults(ProxyCombinedResultMsg *msg) {
  ProxyPatch *patch = (ProxyPatch *)PatchMap::Object()->patch(msg->patch);
  ProxyCombinedResultMsg *ocMsg = patch->depositCombinedResultMsg(msg);
  if (ocMsg) {
	ProxyCombinedResultRawMsg *cMsg = ProxyCombinedResultMsg::toRaw(ocMsg);        
    int destPe = patch->getSpanningTreeParent();
    CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
    CmiAssert(destPe!=CkMyPe());
    //if(destPe != CkMyPe()) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
      /*CkPrintf("ready to call node::recvImmRes on pe[%d] to dest[%d]\n", CkMyPe(), destPe);
      fflush(stdout);*/

      cMsg->destPe = destPe;
      CProxy_NodeProxyMgr cnp(CkpvAccess(BOCclass_group).nodeProxyMgr);
      cnp[CkNodeOf(destPe)].recvImmediateResults(cMsg);
#else
      cp[destPe].recvImmediateResults(cMsg);
#endif
    //}
    //else{
    ////IT SHOULD NEVER BE ENTERED
    //  cp[destPe].recvResults(cMsg);
    //}
  }
}

void ProxyMgr::recvResults(ProxyCombinedResultRawMsg *omsg) {
	ProxyCombinedResultRawMsg *msg = omsg;

//Chao Mei: hack for QD in case of SMP with immediate msg
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
    if(proxyRecvSpanning && msg->isFromImmMsgCall){
//    CkPrintf("qdcreate called on pe[%d]\n", CkMyPe());
//    fflush(stdout);
        //To compensate for the counter loss for message creation
        //inside the process of immediate message on comm thread
        CkpvAccess(_qd)->create();
    }
#endif

  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  if (home) {
    //printf("Home got a message\n");
    home->receiveResults(msg); // delete done in HomePatch::receiveResults()
  }
  else {
    NAMD_bug("ProxyMgr should receive result message on home processor");
  }
}

void ProxyMgr::recvImmediateResults(ProxyCombinedResultRawMsg *omsg) {
  HomePatch *home = PatchMap::Object()->homePatch(omsg->patch);
  if (home) {
    CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);        
    CmiEnableUrgentSend(1);
    cp[CkMyPe()].recvResults(omsg);
    CmiEnableUrgentSend(0);
  }
  else {
    ProxyPatch *patch = (ProxyPatch *)PatchMap::Object()->patch(omsg->patch);
	ProxyCombinedResultMsg *ocMsg = patch->depositCombinedResultRawMsg(omsg);
    if (ocMsg) {
		CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
		ProxyCombinedResultRawMsg *cMsg = ProxyCombinedResultMsg::toRaw(ocMsg);		
		cp[patch->getSpanningTreeParent()].recvImmediateResults(cMsg);
    }
  }
}

void NodeProxyMgr::recvImmediateResults(ProxyCombinedResultRawMsg *omsg){
    ProxyCombinedResultRawMsg *msg = omsg;
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    //CkPrintf("recvImmRes called on comm thread%d pe[%d]\n", CkMyRank()==CmiMyNodeSize(), CkMyPe());
    //fflush(stdout);

    int destRank = CkRankOf(msg->destPe);
    PatchMap *pmap = localPatchMaps[destRank];
    HomePatch *home = pmap->homePatch(msg->patch);
    if (home) {
#if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
        msg->isFromImmMsgCall = (CkMyRank()==CkMyNodeSize());
#endif
        CProxy_ProxyMgr cp(localProxyMgr);
        CmiEnableUrgentSend(1);
        cp[msg->destPe].recvResults(msg);
        CmiEnableUrgentSend(0);
/*
        char *srcfrom = "Isfrom";
        if(CkMyRank()!=CmiMyNodeSize()) srcfrom="Notfrom";
      CkPrintf("%s comm thread from pe[%d]\n", srcfrom, CkMyPe());
      fflush(stdout);
*/
    }
    else {
        ProxyPatch *patch = (ProxyPatch *)pmap->patch(msg->patch);
        ProxyCombinedResultMsg *ocMsg = patch->depositCombinedResultRawMsg(msg); 
        if (ocMsg) {
            CProxy_NodeProxyMgr cnp(thisgroup);
            ocMsg->destPe = patch->getSpanningTreeParent();
			ProxyCombinedResultRawMsg *cMsg = ProxyCombinedResultMsg::toRaw(ocMsg);
            cnp[CkNodeOf(cMsg->destPe)].recvImmediateResults(cMsg);
        }
    }
#endif
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg, int pcnt, int *pids) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    if(proxySendSpanning == 1) {
        CProxy_NodeProxyMgr cnp(CkpvAccess(BOCclass_group).nodeProxyMgr);
        for(int i=0; i<pcnt-1; i++) {
            ProxyDataMsg *copymsg = (ProxyDataMsg *)CkCopyMsg((void **)&msg);
            cnp[pids[i]].recvImmediateProxyData(copymsg);
        }
        cnp[pids[pcnt-1]].recvImmediateProxyData(msg);
        return;
    }
#endif
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp.recvImmediateProxyData(msg,pcnt,pids);
}

void 
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
//Chao Mei: hack for QD in case of SMP with immediate msg
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
    if(proxySendSpanning && msg->isFromImmMsgCall){
//    CkPrintf("qdcreate called on pe[%d]\n", CkMyPe());
//    fflush(stdout);
	//To compensate for the counter loss for message creation
	//inside the process of immediate message on comm thread
	CkpvAccess(_qd)->create();
    }
#endif
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::recvImmediateProxyData(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);  
  if (proxySendSpanning == 1) {
    // copy the message and send to spanning children
    //int *pids = (int*)alloca(proxy->getSpanningTreeNChild()*sizeof(int));
    //int npid = proxy->getSpanningTreeChild(pids);
    int npid = proxy->getSpanningTreeNChild();
    int *pids = (int *)proxy->getSpanningTreeChildPtr();
    if (npid) {        
        ProxyDataMsg *newmsg = (ProxyDataMsg *)CkCopyMsg((void **)&msg);     
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
        int ntreephs;
        PersistentHandle *treephs = proxy->getSpanningTreePhs(ntreephs);
        CmiAssert(treephs && ntreephs == npid);
        CmiUsePersistentHandle(treephs, ntreephs);
#endif
        ProxyMgr::Object()->sendProxyData(newmsg,npid,pids);
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
        CmiUsePersistentHandle(NULL, 0);
#endif
      #if 0
      //ChaoMei: buggy code??? the spanning tree doesn't always have 2 levels
      //At the second level of the tree immediate messages are not needed
      CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
      cp.recvProxyData(newmsg,npid,pids);
      #endif
    }
  }
  /* send to self via EP method to preserve priority */
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp[CkMyPe()].recvProxyData(msg);
}

void NodeProxyMgr::recvImmediateProxyData(ProxyDataMsg *msg) {    
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    CProxy_ProxyMgr cp(localProxyMgr);
    proxyTreeNode *ptn = proxyInfo[msg->patch];
    CmiAssert(ptn->numPes!=0);

    //re-send msg to this nodes's children nodes.
    //only the first pe of a node of node-aware ST should contain children nodes
    int rank = CkRankOf(ptn->peIDs[0]);
    PatchMap *pmap = localPatchMaps[rank];
    ProxyPatch *ppatch = (ProxyPatch *)pmap->patch(msg->patch);

    int npid = ppatch->getSTNNodeChild();
    int *pids = ppatch->getSTNodeChildPtr();
    if(npid>0) {        
        //only needs to send to other nodes, so check the last entry of pids.
        //This is because the data for proxies on the same node have been sent
        //later in this function by NodeProxyMgr.
        if(pids[npid-1]==CkMyNode()) npid--;
    }    
    CProxy_NodeProxyMgr cnp(thisgroup);
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    if (npid) {
        int ntreephs;
        PersistentHandle *treephs = ppatch->getSpanningTreePhs(ntreephs);
        CmiAssert(treephs && ntreephs >= npid);
        CmiUsePersistentHandle(treephs, ntreephs);
    }
#endif
    for(int i=0; i<npid; i++) {
        ProxyDataMsg *copymsg = (ProxyDataMsg *)CkCopyMsg((void **)&msg);
        cnp[pids[i]].recvImmediateProxyData(copymsg);
    }    
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    CmiUsePersistentHandle(NULL, 0);
#endif

    //re-send msg to it's internal cores
#if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
    msg->isFromImmMsgCall = (CkMyRank()==CkMyNodeSize());
#endif
    cp.recvProxyData(msg, ptn->numPes, ptn->peIDs);
#else
    CkAbort("Bad execution path to NodeProxyMgr::recvImmediateProxyData\n");
#endif
}

void
ProxyMgr::sendProxyAll(ProxyDataMsg *msg, int pcnt, int *pids) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    if(proxySendSpanning == 1) {
        CProxy_NodeProxyMgr cnp(CkpvAccess(BOCclass_group).nodeProxyMgr);
        for(int i=0; i<pcnt-1; i++) {
            ProxyDataMsg *copymsg = (ProxyDataMsg *)CkCopyMsg((void **)&msg);
            cnp[pids[i]].recvImmediateProxyAll(copymsg);
        }
        cnp[pids[pcnt-1]].recvImmediateProxyAll(msg);
        return;
    }
#endif
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp.recvImmediateProxyAll(msg,pcnt,pids);
}

void 
ProxyMgr::recvProxyAll(ProxyDataMsg *msg) {
//Chao Mei: hack for QD in case of SMP with immediate msg
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
    if(proxySendSpanning && msg->isFromImmMsgCall){
//    CkPrintf("qdcreate called on pe[%d]\n", CkMyPe());
//    fflush(stdout);
	//To compensate for the counter loss for message creation
	//inside the process of immediate message on comm thread
	CkpvAccess(_qd)->create();
    }
#endif

  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAll(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::recvImmediateProxyAll(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
  DebugFileTrace *dft = DebugFileTrace::Object();
  dft->openTrace();
  dft->writeTrace("PMgr::recvImmPAll for patch[%d]\n", msg->patch);
  CmiAssert(proxy!=NULL);
  dft->writeTrace("PMgr::recvImmPAll assertion OK for patch[%d]\n", msg->patch);
  dft->closeTrace();
  #endif
  if (proxySendSpanning == 1) {
    // copy the message and send to spanning children
    //int *pids = (int*)alloca(proxy->getSpanningTreeNChild()*sizeof(int));
    //int npid = proxy->getSpanningTreeChild(pids);
    int npid = proxy->getSpanningTreeNChild();
    int *pids = (int *)proxy->getSpanningTreeChildPtr();
    if (npid) {
        ProxyDataMsg *newmsg = (ProxyDataMsg *)CkCopyMsg((void **)&msg);      
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
        int ntreephs;
        PersistentHandle *treephs = proxy->getSpanningTreePhs(ntreephs);
        CmiAssert(treephs && ntreephs == npid);
        CmiUsePersistentHandle(treephs, ntreephs);
#endif
        ProxyMgr::Object()->sendProxyAll(newmsg,npid,pids);
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
        CmiUsePersistentHandle(NULL, 0);
#endif
    }
  }
  /* send to self via EP method to preserve priority */
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  cp[CkMyPe()].recvProxyAll(msg);
}

void NodeProxyMgr::recvImmediateProxyAll(ProxyDataMsg *msg) {    
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    CProxy_ProxyMgr cp(localProxyMgr);
    proxyTreeNode *ptn = proxyInfo[msg->patch];
    CmiAssert(ptn->numPes!=0);
    #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    //This could be executed on comm thd.
    printf("NodePMgr::recvImmPAll for patch[%d] on node %d rank %d, prepare to send proc ", msg->patch, CkMyNode(), CkMyRank());
    for(int i=0; i<ptn->numPes; i++) {
        printf("%d, ", ptn->peIDs[i]);
    }
    printf("\n");
    fflush(stdout);
    #endif

    //re-send msg to this nodes's children nodes.
    //only the first pe of a node of node-aware ST should contain children nodes
    int rank = CkRankOf(ptn->peIDs[0]);
    PatchMap *pmap = localPatchMaps[rank];
    ProxyPatch *ppatch = (ProxyPatch *)pmap->patch(msg->patch);

    int npid = ppatch->getSTNNodeChild();
    int *pids = ppatch->getSTNodeChildPtr();
    if(npid>0) {        
        //only needs to send to other nodes, so check the last entry of pids.
        //This is because the data for proxies on the same node have been sent
        //later in this function by NodeProxyMgr.
        if(pids[npid-1]==CkMyNode()) npid--;
    }
    
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    if (npid) {
        int ntreephs;
        PersistentHandle *treephs = ppatch->getSpanningTreePhs(ntreephs);
        CmiAssert(treephs && ntreephs >= npid);
        CmiUsePersistentHandle(treephs, ntreephs);
    }
#endif
    CProxy_NodeProxyMgr cnp(thisgroup);
    for(int i=0; i<npid; i++) {
        ProxyDataMsg *copymsg = (ProxyDataMsg *)CkCopyMsg((void **)&msg);
        cnp[pids[i]].recvImmediateProxyAll(copymsg);
    }    
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    CmiUsePersistentHandle(NULL, 0);
#endif

    //re-send msg to it's internal cores
#if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
    msg->isFromImmMsgCall = (CkMyRank()==CkMyNodeSize());
#endif
    cp.recvProxyAll(msg, ptn->numPes, ptn->peIDs);
#else
    CkAbort("Bad execution path to NodeProxyMgr::recvImmediateProxyData\n");
#endif
}

void ProxyMgr::printProxySpanningTree(){
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    int numPatches = PatchMap::Object()->numPatches();
    for(int i=0; i<numPatches; i++) {
        proxyTreeNodeList &oneList = ptree.naTrees[i];
        printf("ST tree for HomePatch[%d]: #nodes = %d\n", i, oneList.size()); 
        if(ptree.proxylist[i].size()==0) continue;
        printf("===%d=== pes/node: ", i);
        for(int j=0; j<oneList.size(); j++) {
            printf("%d ", oneList.item(j).numPes);
        }
        printf("\n");
        printf("===%d=== pe ids: ", i);
        for(int j=0; j<oneList.size(); j++) {
            for(int k=0; k<oneList.item(j).numPes; k++) {
                printf("%d ", oneList.item(j).peIDs[k]);
            }            
        }
        printf("\n");
    }    
    fflush(stdout);  
#else
    int numPatches = PatchMap::Object()->numPatches();
    for(int i=0; i<numPatches; i++) {
        NodeIDList &oneList = ptree.trees[i];
        printf("ST tree for HomePatch[%d]: #nodes = %d\n", i, oneList.size()); 
        if(ptree.proxylist[i].size()==0) continue;        
        printf("===%d=== pe ids: ", i);
        for(int j=0; j<oneList.size(); j++) {            
            printf("%d ", oneList.item(j));            
        }
        printf("\n");
    }    
    fflush(stdout);  
#endif
}

void NodeProxyMgr::registerPatch(int patchID, int numPes, int *pes){
    if(proxyInfo[patchID]) {
        delete proxyInfo[patchID];
    }
    if(numPes == 0) {
        proxyInfo[patchID] = NULL;
    }else{
        proxyInfo[patchID] = new proxyTreeNode(CkNodeOf(pes[0]),numPes,pes);
    }
}

void ProxyMgr::sendResult(ProxyGBISP1ResultMsg *msg) { //pp -r> hp
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  NodeID node = PatchMap::Object()->node(msg->patch);
  CmiEnableUrgentSend(1);
  cp[node].recvResult(msg);
  CmiEnableUrgentSend(0);
}
void ProxyMgr::recvResult(ProxyGBISP1ResultMsg *msg) { //pp -r> hp
  HomePatch *homePatch = PatchMap::Object()->homePatch(msg->patch);
  homePatch->receiveResult(msg); // message deleted in registerProxy()
}
void ProxyMgr::recvData(  ProxyGBISP2DataMsg *msg) {  //hp -d> pp
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms() ?
}
void ProxyMgr::sendResult(ProxyGBISP2ResultMsg *msg) { //pp -r> hp
  CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
  NodeID node = PatchMap::Object()->node(msg->patch);
  CmiEnableUrgentSend(1);
  cp[node].recvResult(msg);
  CmiEnableUrgentSend(0);
}
void ProxyMgr::recvResult(ProxyGBISP2ResultMsg *msg) { //pp -r> hp
  HomePatch *homePatch = PatchMap::Object()->homePatch(msg->patch);
  homePatch->receiveResult(msg); // message deleted in registerProxy()
}
void ProxyMgr::recvData(  ProxyGBISP3DataMsg *msg) {   //hp -d> pp
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms() ?
}

PatchProxyListMsg *PatchProxyListMsg::createPatchProxyListMsg(PatchProxyListMsg **bufs, int bufSize, ProxyListInfo *info, int size){
	//1. compute the total patches this node manages, and the total length of all proxy lists
	int totalPatches = 0;
	int totalProxies = 0;
	for(int i=0; i<bufSize; i++) {
		PatchProxyListMsg *one = bufs[i];
		totalPatches += one->numPatches;
		for(int j=0; j<one->numPatches; j++) totalProxies += one->proxyListLen[j];
	}
	totalPatches += size;
	for(int i=0; i<size; i++) {
		totalProxies += info[i].numProxies;
	}

	PatchProxyListMsg *msg = new(totalPatches, totalPatches, totalProxies, 0)PatchProxyListMsg(totalPatches);
	int msgPatchIdx = 0;
	int msgProxyPeIdx = 0;
	for(int i=0; i<bufSize; i++) {
		PatchProxyListMsg *one = bufs[i];
		int curPeIdx = 0;
		for(int j=0; j<one->numPatches; j++) {
			msg->patchIDs[msgPatchIdx] = one->patchIDs[j];
			int curListLen = one->proxyListLen[j];
			msg->proxyListLen[msgPatchIdx++] = curListLen;
			memcpy(msg->proxyPEs+msgProxyPeIdx, one->proxyPEs+curPeIdx, sizeof(int)*curListLen);
			curPeIdx += curListLen;
			msgProxyPeIdx += curListLen;
		}
	}
	for(int i=0; i<size; i++) {
		msg->patchIDs[msgPatchIdx] = info[i].patchID;
		int curListLen = info[i].numProxies;
		msg->proxyListLen[msgPatchIdx++] = curListLen;
		memcpy(msg->proxyPEs+msgProxyPeIdx, info[i].proxyList, sizeof(int)*curListLen);
		msgProxyPeIdx += curListLen;
	}
	return msg;
}

#define HOMEPATCH_TREE_BRFACTOR 2
void NodeProxyMgr::createSTForHomePatches(PatchMap *pmap){
	//We use implicit tree construction for all home patches
	std::vector<int> nodesWithPatches; //record the id of node that has home patches
	int myNodeIdx = -1; //the index into the above vector of this node
	for(int nodeId=0; nodeId<CkNumNodes(); ++nodeId) {
		int hpCnt = 0;
                int firstPe = CkNodeFirst(nodeId);
                int endPe = firstPe + CkNodeSize(nodeId);
		for(int pe=firstPe; pe < endPe; ++pe) {
			hpCnt += pmap->numPatchesOnNode(pe);
		}
		if(hpCnt==0) continue;

		nodesWithPatches.push_back(nodeId);
		if(CkMyNode() == nodeId) {
			//on my node
			myNodeIdx = nodesWithPatches.size()-1;
			numHomePatches = hpCnt;
			homepatchRecved = 0;
			localProxyLists = new ProxyListInfo[hpCnt];
			memset(localProxyLists, 0, sizeof(ProxyListInfo)*hpCnt);			
		}
	}

	if(myNodeIdx==-1){
		//there's no home patches on this node
		//just set to a value that doesn't make sense in spanning tree.
		parentNode = -2; 
		numKidNodes = 0;
		kidRecved = 0;
		return;
	}
	
	//calculate parent
	if(myNodeIdx == 0) {
		parentNode = -1;
	}else{
		int parentIdx = (myNodeIdx-1)/HOMEPATCH_TREE_BRFACTOR;
		parentNode = nodesWithPatches[parentIdx];
	}

	//calculate kids
	numKidNodes = 0;
	int totalNodes = nodesWithPatches.size();
	for(int i=1; i<=HOMEPATCH_TREE_BRFACTOR; i++) {
		int kidId = myNodeIdx*HOMEPATCH_TREE_BRFACTOR+i;
		if(kidId >= totalNodes) break;
		numKidNodes++;
	}
	if(numKidNodes!=0) {
		remoteProxyLists = new PatchProxyListMsg *[numKidNodes];
	}
	kidRecved = 0;

	//CkPrintf("Node[%d] has %d homepatches with parent=%d and %d kids \n", CkMyNode(), numHomePatches, parentNode, numKidNodes);
}

void NodeProxyMgr::sendProxyList(int pid, int *plist, int size){
	int insertIdx; //indexed from 0
	CmiLock(localDepositLock);
	insertIdx = homepatchRecved++; //ensure the atomic increment

	localProxyLists[insertIdx].patchID = pid;
	localProxyLists[insertIdx].numProxies = size;
	localProxyLists[insertIdx].proxyList = plist;

	if(insertIdx == (numHomePatches-1)) {
		//all local home patches have contributed
		contributeToParent();
	}
	CmiUnlock(localDepositLock);
}

void NodeProxyMgr::sendProxyListInfo(PatchProxyListMsg *msg){
	int insertIdx; //indexed from 0
	CmiLock(localDepositLock);
	insertIdx = kidRecved++;
	
	remoteProxyLists[insertIdx] = msg;
	if(insertIdx == (numKidNodes-1)) {
		//all kids have contributed;
		contributeToParent();
	}
	CmiUnlock(localDepositLock);
}

void NodeProxyMgr::contributeToParent(){
	if(homepatchRecved!=numHomePatches || kidRecved != numKidNodes) return;

	homepatchRecved = 0;
	kidRecved = 0;
	//construct the msg
	PatchProxyListMsg *msg = PatchProxyListMsg::createPatchProxyListMsg(remoteProxyLists, numKidNodes, localProxyLists, numHomePatches);
	if(parentNode == -1) {
		//send to proxy mgr on PE[0] as this is the root node
		CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
		cp[0].recvPatchProxyInfo(msg);
	}else{
		CProxy_NodeProxyMgr cnp(thisgroup);
		cnp[parentNode].sendProxyListInfo(msg);
	}
	for(int i=0; i<numKidNodes; i++) {
		delete remoteProxyLists[i];
	}
}

#include "ProxyMgr.def.h"

