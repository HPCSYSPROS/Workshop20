/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/TorusLB.C,v $
 * $Author: jim $
 * $Date: 2013/08/30 18:18:20 $
 * $Revision: 1.35 $
 *****************************************************************************/
 
/** \file TorusLB.C
 *  Author: Abhinav S Bhatele
 *  Date Created: June 05th, 2007 
 *
 *  Replacement for AlgSeven.C
 */

#include "TorusLB.h"
#include "ProxyMgr.h"
#define SHRINK_INNER_BRICK 1

TorusLB::TorusLB(computeInfo *cs, patchInfo *pas, processorInfo *pes, int ncs, 
int npas, int npes) : RefineTorusLB(cs, pas, pes, ncs, npas, npes, 0)
{
  strategyName = "TorusLB";
  strategy();
  //if ( computeMax() <= origMaxLoad ) {
  //  binaryRefine();
  //  printLoads();
  //}
  // CREATE THE SPANNING TREE IN THE LOAD BALANCER
  //if(proxySendSpanning || proxyRecvSpanning)
  //  createSpanningTree();
}

TorusLB::~TorusLB() { }

void TorusLB::strategy() {
  int index;
  // compute the average load by (compute load + background load) / numPesAvailable
  computeAverage();
  // two heaps of self and pair computes
  makeTwoHeaps();

  const int beginGroup = processors[0].Id;
  const int endGroup = beginGroup + P;
#define INGROUP(PROC) ((PROC) >= beginGroup && (PROC) < endGroup)

  computeInfo *c;
  processorInfo *p, *minp;
  Iterator nextP;
  overLoad = 1.2;

  for(int I=0; I<numComputes; I++) {

  c = (computeInfo *) computePairHeap->deleteMax();
  if ( ! c ) c = (computeInfo *) computeSelfHeap->deleteMax(); 

  if(c->processor != -1) continue; // go to the next compute
  if(!c) CkAbort("TorusLB: Compute Heap empty!\n");

  for(int j=0; j<6; j++) {
    bestPe[j] = 0;
    goodPe[j] = 0;
    badPe[j] = 0;
  }

  // Look at pes which have the compute's patches

  // HYBRID check if processor is in local group
#define SELECT_REALPE(X) if INGROUP((X)) { \
  selectPes(&processors[(X) - beginGroup], c); \
  }

  const int realPe1 = patches[c->patch1].processor;
  SELECT_REALPE(realPe1)

  const int realPe2 = patches[c->patch2].processor;
  if ( realPe2 != realPe1 ) {
    SELECT_REALPE(realPe2)
  }

  // Try the processors which have the patches' proxies
  p = (processorInfo *)(patches[c->patch1].proxiesOn.iterator((Iterator *)&nextP));
  while(p) {						// patch 1
    if INGROUP(p->Id) selectPes(p, c);
    p = (processorInfo *)(patches[c->patch1].proxiesOn.next((Iterator *)&nextP));
  } 

  p = (processorInfo *)(patches[c->patch2].proxiesOn.iterator((Iterator *)&nextP));
  while(p) {						// patch 2
    if INGROUP(p->Id) selectPes(p, c);
    p = (processorInfo *)(patches[c->patch2].proxiesOn.next((Iterator *)&nextP));
  }

  // see if we have found a processor to place the compute on
  p = 0;
  if((p = bestPe[5])
#if USE_TOPOMAP
  || (p = goodPe[5])
#endif
  || (p = bestPe[4])
#if USE_TOPOMAP
  || (p = goodPe[4])
#endif
  || (p = bestPe[3])
#if USE_TOPOMAP
  || (p = goodPe[3])
#endif
  || (p = bestPe[1])
#if USE_TOPOMAP
  || (p = goodPe[1])
#endif
  || (p = bestPe[2])
#if USE_TOPOMAP
  || (p = goodPe[2])
#endif
  || (p = bestPe[0])
#if USE_TOPOMAP
  || (p = goodPe[0])
#endif
  ) {
    assign(c, p);
    continue;
  }

    // Try all pes on the nodes of the home patches
    if ( CmiNumNodes() > 1 ) {  // else not useful
      double minLoad = overLoad * averageLoad;
      minp = 0;
      int realNode1 = CmiNodeOf(realPe1);
      int nodeSize = CmiNodeSize(realNode1);
      if ( nodeSize > 1 ) {  // else did it already
        int firstpe = CmiNodeFirst(realNode1);
        for ( int rpe = firstpe; rpe < firstpe+nodeSize; ++rpe ) {
          if INGROUP(rpe) {
            p = &processors[rpe - beginGroup];
            if ( p->available && ( p->load + c->load < minLoad ) ) {
              minLoad = p->load + c->load;
              minp = p;
            }
          }
        }
      }
      if ( realPe2 != realPe1 ) {
        int realNode2 = CmiNodeOf(realPe2);
        if ( realNode2 != realNode1 ) {  // else did it already
          nodeSize = CmiNodeSize(realNode2);
          if ( nodeSize > 1 ) {
            int firstpe = CmiNodeFirst(realNode2);
            for ( int rpe = firstpe; rpe < firstpe+nodeSize; ++rpe ) {
              if INGROUP(rpe) {
                p = &processors[rpe - beginGroup];
                if ( p->available && ( p->load + c->load < minLoad ) ) {
                  minLoad = p->load + c->load;
                  minp = p;
                }
              }
            }
          }
        }
      }
      if(minp) {
        assign(c, minp);
        continue;
      }
    }

    // Try all pes on the physical nodes of the home patches
    if ( ( CmiNumPhysicalNodes() > 1 ) &&
         ( CmiNumPhysicalNodes() < CmiNumNodes() ) ) {  // else not useful
      double minLoad = overLoad * averageLoad;
      minp = 0;
      int realNode1 = CmiPhysicalNodeID(realPe1);
      int *rpelist;
      int nodeSize;
      CmiGetPesOnPhysicalNode(realNode1, &rpelist, &nodeSize);
      if ( nodeSize > 1 ) {  // else did it already
        for ( int ipe = 0; ipe < nodeSize; ++ipe ) {
          int rpe = rpelist[ipe];
          if INGROUP(rpe) {
            p = &processors[rpe - beginGroup];
            if ( p->available && ( p->load + c->load < minLoad ) ) {
              minLoad = p->load + c->load;
              minp = p;
            }
          }
        }
      }
      if ( realPe2 != realPe1 ) {
        int realNode2 = CmiPhysicalNodeID(realPe2);
        if ( realNode2 != realNode1 ) {  // else did it already
          CmiGetPesOnPhysicalNode(realNode2, &rpelist, &nodeSize);
          if ( nodeSize > 1 ) {  // else did it already
            for ( int ipe = 0; ipe < nodeSize; ++ipe ) {
              int rpe = rpelist[ipe];
              if INGROUP(rpe) {
                p = &processors[rpe - beginGroup];
                if ( p->available && ( p->load + c->load < minLoad ) ) {
                  minLoad = p->load + c->load;
                  minp = p;
                }
              }
            }
          }
        }
      }
      if(minp) {
        assign(c, minp);
        continue;
      }
    }

 
  int found = 0;
#if USE_TOPOMAP
  // If no processor found, go through the whole list in a topological fashion
  // first try the inner brick
  int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM, t1, t2;
  int dimNX, dimNY, dimNZ, dimNT;
  double minLoad;
  p1 = patches[c->patch1].processor;
  p2 = patches[c->patch2].processor;

  tmgr.rankToCoordinates(p1, x1, y1, z1, t1);
  tmgr.rankToCoordinates(p2, x2, y2, z2, t2);
  dimNX = tmgr.getDimNX();
  dimNY = tmgr.getDimNY();
  dimNZ = tmgr.getDimNZ();
  dimNT = tmgr.getDimNT();

  brickDim(x1, x2, dimNX, xm, xM);
  brickDim(y1, y2, dimNY, ym, yM);
  brickDim(z1, z2, dimNZ, zm, zM);

  // to shrink the inner brick by some hops
#if 0
  xm=xm+SHRINK_INNER_BRICK;
  ym=ym+SHRINK_INNER_BRICK;
  zm=zm+SHRINK_INNER_BRICK;

  xM=xM-SHRINK_INNER_BRICK;
  yM=yM-SHRINK_INNER_BRICK;
  zM=zM-SHRINK_INNER_BRICK;
#endif

  // first go over the processors inside the brick and choose the least 
  // overloaded one
  p = 0; minp = 0;
  minLoad = overLoad * averageLoad;
  for(int i=xm; i<=xM; i++)
    for(int j=ym; j<=yM; j++)
      for(int k=zm; k<=zM; k++)
        for(int l=0; l<dimNT; l++)
        {
          pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
          if ( ! INGROUP(pe) ) continue;
          p = &processors[pe - beginGroup];
          if(c->load + p->load < minLoad) { 
            minLoad = c->load + p->load;
            minp = p;
	    found = 1;
          }
        }

  // if no success, then go through the remaining torus (outer brick)
  // and pick the first underloaded one
  minLoad = overLoad * averageLoad;
  if(found == 0) {
    p = 0; minp = 0;
    for(int i=xM+1; i<xm+dimNX; i++)
      for(int j=0; j<dimNY; j++)
	for(int k=0; k<dimNZ; k++)
          for(int l=0; l<dimNT; l++)
          {
            pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
            if ( ! INGROUP(pe) ) continue;
            p = &processors[pe - beginGroup];
            if(c->load + p->load < minLoad) { 
              minp = p;
	      found = 1; break;
            }
          }
  }

  if(found == 0) {
    for(int j=yM+1; j<ym+dimNY; j++)
      for(int i=xm; i<=xM; i++)
	for(int k=0; k<dimNZ; k++)
          for(int l=0; l<dimNT; l++)
          {
            pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
            if ( ! INGROUP(pe) ) continue;
            p = &processors[pe - beginGroup];
            if(c->load + p->load < minLoad) { 
              minp = p;
	      found = 1; break;
            }
          }
  }

  if(found == 0) {
    for(int k=zM+1; k<zm+dimNZ; k++)
      for(int i=xm; i<=xM; i++)
        for(int j=ym; j<=yM; j++)
          for(int l=0; l<dimNT; l++)
          {
            pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
            if ( ! INGROUP(pe) ) continue;
            p = &processors[pe - beginGroup];
            if(c->load + p->load < minLoad) { 
              minp = p;
	      found = 1; break;
            }
          }
  }

  
  if(found == 1) {
    assign(c, minp);
    continue;
  }

#endif /* USE_TOPOMAP */

  if(found == 0) {
    heapIterator nextp;
    processorInfo *p = (processorInfo *)(pes->iterator((heapIterator *) &nextp));
    while (p) {
      selectPes(p, c);
      p = (processorInfo *)(pes->next(&nextp));
    }
    p = 0;
    if((p = bestPe[5])
#if USE_TOPOMAP
    || (p = goodPe[5])
#endif
    || (p = bestPe[4])
#if USE_TOPOMAP
    || (p = goodPe[4])
#endif
    || (p = bestPe[3])
#if USE_TOPOMAP
    || (p = goodPe[3])
#endif
    || (p = bestPe[1])
#if USE_TOPOMAP
    || (p = goodPe[1])
#endif
    || (p = bestPe[2])
#if USE_TOPOMAP
    || (p = goodPe[2])
#endif
    || (p = bestPe[0])
#if USE_TOPOMAP
    || (p = goodPe[0])
#endif
    ) {
      assign(c, p);
      found = 1;
      continue;
    }
  }

  if(found == 0) {
    p = 0;
    if((p = badPe[5])
    || (p = badPe[4])
    || (p = badPe[3])
    || (p = badPe[1])
    || (p = badPe[2])
    || (p = badPe[0])) {
      assign(c, p);
      found = 1;
      continue;
    }
  }

  if(found == 0) {
     CkPrintf("TorusLB: No receiver found average %f overload %f\n", averageLoad, overLoad);
     CkAbort("TorusLB: No receiver found\n");
  }
 
  } // end of computes for-loop

  printLoads(3);
}

void TorusLB::selectPes(processorInfo *p, computeInfo *c) {
  if (p->available == false)
    return;

  // find the position in bestPe/goodPe to place this pair
  // HP HP HP HP HP HP
  // 02 11 20 01 10 00
  //  5  4  3  2  1  0
  int numPatches, numProxies, /* badForComm, */ index;
  numAvailable(c, p, &numPatches, &numProxies, 0 /* &badForComm */); 
  int numEither = numPatches + numProxies;
  index = (numEither*(numEither+1))/2 + numProxies;

#if USE_TOPOMAP
  int x, y, z, t;
  int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM, t1, t2;
  int dimNX, dimNY, dimNZ, dimNT;
  double minLoad;
  p1 = patches[c->patch1].processor;
  p2 = patches[c->patch2].processor;

  tmgr.rankToCoordinates(p1, x1, y1, z1, t1);
  tmgr.rankToCoordinates(p2, x2, y2, z2, t2);
  dimNX = tmgr.getDimNX();
  dimNY = tmgr.getDimNY();
  dimNZ = tmgr.getDimNZ();
  dimNT = tmgr.getDimNT();

  brickDim(x1, x2, dimNX, xm, xM);
  brickDim(y1, y2, dimNY, ym, yM);
  brickDim(z1, z2, dimNZ, zm, zM);
#endif

  if (p->load + c->load < overLoad * averageLoad) {
  // replace only if the new processor is underloaded
#if USE_TOPOMAP
    tmgr.rankToCoordinates(p->Id, x, y, z, t);
    int wB = withinBrick(x, y, z, xm, xM, dimNX, ym, yM, dimNY, zm, zM, dimNZ);
    if (wB) {
      // if within the brick formed by the patch processors
#endif
      // or the non-topology case
      processorInfo* &oldp = bestPe[index];
      if (!(oldp) || p->load < oldp->load )
        oldp = p;
#if USE_TOPOMAP
    } else {
      // if outside the brick formed by the patch processors
      processorInfo* &oldp = goodPe[index];
      double loadDiff = 0.0;
      
      if (!(oldp)) {
	// replace if there is no processor at that position
        oldp = p;
      }
      else {
	// get the load difference if the processor exixts
	loadDiff = oldp->load - p->load;
	if ((loadDiff > 0.4) || (loadDiff > 0.0 && (tmgr.getHopsBetweenRanks(p->Id, p1) + 
	    tmgr.getHopsBetweenRanks(p->Id, p2) < tmgr.getHopsBetweenRanks(oldp->Id, p1) + 
	    tmgr.getHopsBetweenRanks(oldp->Id, p2)))) {
	  // if weights are similar, look at the hops
          oldp = p;
	}
      }
    }
#endif
  } else {
    // for the first placement, we must find a processor to place
    // the compute on, so choose a bad processor if necessary
    processorInfo* &oldp = badPe[index];
    if (!(oldp) || p->load < oldp->load )
      oldp = p;
  }
}

