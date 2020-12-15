
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
    Sync will ensure that all homepatches finished updating before Computes starts and all proxies finished updating themselves.
*/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <stdio.h>

#include "InfoStream.h"
#include "Patch.h"
#include "PatchMap.h"
#include "ProxyMgr.h"
#include "Compute.h"
#include "ComputeMap.h"

#include "Sync.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

// make sure all HomePatches get their positions data and sendProxyData to 
// their proxies before computes get positionsReady.

// useProxySync will make sure all proxies get updated before computes' 
// positionsReady triggered so that real computation starts.
// when these two combined, it will make sure that homepatch get all its force 
// and positions data, and proxies receive its updated data before all 
// computes start.

static int eventHoldComputes;
static int eventReleaseComputes;

Sync::Sync(): INCREASE(600), step(0), counter(0), homeReady(0)
{
    if (CkpvAccess(Sync_instance) == NULL) {
        CkpvAccess(Sync_instance) = this;
    } else {
	NAMD_bug("Sync instanced twice on same processor!");
    }
    capacity = INCREASE;
    clist = new _clist[capacity];
    cnum = 0;
    nPatcheReady = 0;
    numPatches = -1;
    eventHoldComputes = traceRegisterUserEvent("Sync::holdComputes", 133);
    eventReleaseComputes = traceRegisterUserEvent("Sync::releaseComputes", 134);
}

Sync::~Sync()
{
  delete [] clist;
}

void Sync::openSync(void)
{
  int reportPe = 1;
  while ( 2 * reportPe < CkNumPes() ) reportPe *= 2;
  step = -1;
  useSync = 1;
  if ( PatchMap::Object()->numPatches() >= 4 * CkNumPes() ) useSync = 0;
  if ( CmiNumNodes() < 2 ) useSync = 0;
  if ( CmiNumPhysicalNodes() < 2 ) useSync = 0;
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  useSync = 0;
#endif
  useProxySync = 0;
  if (useSync) {
    // if use proxy spanning tree, proxy sync is forced
    if (!useProxySync && (proxySendSpanning || proxyRecvSpanning)
        && PatchMap::Object()->numPatches() < 4 * CkNumPes() ) {
      // If on BG/P, useProxySync should not be turned on for better performance
      #if ! (CMK_BLUEGENEQ || CMK_BLUEGENEP)
      // CmiPrintf("[%d] useProxySync is turned on. \n", CkMyPe());
      useProxySync = 1;
      #endif
    }
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    // immediate messages can be processed by any PE
    if (CkMyNodeSize() > 2) useProxySync = 0;
#endif
    // no proxies on this node, no need to use proxy sync.
    if (useProxySync && ProxyMgr::Object()->numProxies() == 0) {
      // CmiPrintf("[%d] useProxySync is turned off because no proxy. \n", CkMyPe());
      useProxySync = 0;
    }
    // if no proxy sync and no home patch, then disable home patch sync as well
    if (!useProxySync && PatchMap::Object()->numHomePatches() == 0) useSync = 0;
  }
  if(CkMyPe() == reportPe)
    iout << iINFO << "useSync: " << useSync << " useProxySync: " << useProxySync << "\n" << endi;
}    

// called from Patch::positionsReady()
int Sync::holdComputes(PatchID pid, Compute **cbegin, Compute **cend, int doneMigration, int seq)
{
  if (!useSync) return 0;
  if (step < 0) step = seq;
  if (!useProxySync) {
    // only hold when homepatches are not ready
    PatchMap *patchMap = PatchMap::Object();
    if (homeReady && seq == step) {
      nPatcheReady++;
      triggerCompute();
      return 0;
    }
  }
  traceUserEvent(eventHoldComputes);

  int slot = 0;
  for (; slot < cnum; slot++)
     if (clist[slot].pid == -1) break;
  if (slot == cnum) {
    cnum++;
    // table is full, expand the list
    if (cnum == capacity) {
      capacity += INCREASE;
      struct _clist *tmp = new _clist[capacity];
      memcpy(tmp, clist, cnum*sizeof(_clist));
      delete [] clist;
      clist = tmp;
      //CmiPrintf("[%d] Info:: Sync buffer overflow and expanded!\n", CkMyPe());
    }
  }

  clist[slot].cbegin = cbegin;
  clist[slot].cend = cend;
  clist[slot].pid = pid;
  clist[slot].doneMigration  = doneMigration;
  clist[slot].step = seq;

//  CkPrintf("REG[%d]: patch:%d step:%d-%d slot:%d\n", CkMyPe(), pid, patchMap->patch(pid)->flags.sequence, step, slot);

  if (clist[slot].step == step) {
      nPatcheReady++;
      triggerCompute();
  }
  return 1;
}

// called from HomePatch::positionsReady()
void Sync::PatchReady(void)
{
 if ( useSync ) {
  counter ++;
  triggerCompute();
 }
}

void Sync::releaseComputes()
{
  PatchMap *patchMap = PatchMap::Object();

  traceUserEvent(eventReleaseComputes);

  for (int i= 0; i<cnum; i++) {
    int &pid = clist[i].pid;
    if (pid == -1) continue;
    if (clist[i].step != step) {
      continue;
    }
    //         CkPrintf(" %d-%d-%d ",
    //	 clist[i].pid, clist[i].step,
    //      patchMap->patch(pid)->flags.sequence);

    Compute **cend = clist[i].cend;
    for(Compute **cid = clist[i].cbegin; cid != cend; cid++) {
      (*cid)->patchReady(pid,clist[i].doneMigration,step);
    }
    pid = -1;
  }
//  CkPrintf("\n");
}

void Sync::triggerCompute()
{
  PatchMap *patchMap = PatchMap::Object();
  const int numHomePatches = patchMap->numHomePatches();

  if (numPatches == -1) 
    numPatches = ProxyMgr::Object()->numProxies() + numHomePatches;

// if (CkMyPe()<=8) CkPrintf("SYNC[%d]: PATCHREADY:%d %d patches:%d %d\n", CkMyPe(), counter, numHomePatches, nPatcheReady, numPatches);
//  CkPrintf("SYNC[%d]: PATCHREADY:%d %d patches:%d %d\n", CkMyPe(), counter, PatchMap::Object()->numHomePatches(), nPatcheReady, numPatches);

  if (homeReady == 0 && counter >= numHomePatches) {
    homeReady = 1;
 // if (CkMyPe()<=8) CkPrintf("HOMEREADY[%d]\n", CkMyPe());
    if (!useProxySync)  releaseComputes();
  }

  if (homeReady && nPatcheReady == numPatches)
  {
// if (CkMyPe()<=8) CkPrintf("TRIGGERED[%d]\n", CkMyPe());
//     CkPrintf("TRIGGERED[%d]\n", CkMyPe());
    if (useProxySync) releaseComputes();

    // reset counter
    numPatches = -1;
    step++;
    nPatcheReady = 0;
    for (int i= 0; i<cnum; i++) {
      if (clist[i].pid != -1 && clist[i].step == step) ++nPatcheReady;
    }
    homeReady = 0;
    if ( numHomePatches ) {
      counter -= numHomePatches;
      if (counter >= numHomePatches) triggerCompute();
    }
  }
}


#if 0
// hack for SUN MPCC
// force compiler to instantiate ArrayElementT
void *frightenCompilerIntoInstantiatingHack(void) {
  return new ArrayElementT<int>;
}
#endif


#include "Sync.def.h"
