/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/AlgSeven.C,v $
 * $Author: jim $
 * $Date: 2013/06/07 22:34:36 $
 * $Revision: 1.60 $
 *****************************************************************************/

#include "common.h"
#include "InfoStream.h"
#include "Node.h"
#include "Alg7.h"

#define TINYLOAD 0.0005

Alg7::Alg7(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
  Rebalancer(computeArray, patchArray, 
	     processorArray, nComps, 
	     nPatches, nPes)
{
strategyName = "Alg7";
strategy();
}

extern int isPmeProcessor(int);

void Alg7::togrid(processorInfo* goodP[3][3][2], processorInfo* poorP[3][3][2],
			processorInfo *p, computeInfo *c) {
      if(p->available == false) return;

      int nPatches, nProxies, badForComm;
      numAvailable(c,p,&nPatches,&nProxies,&badForComm);

      if (c->load + p->load < overLoad*averageLoad) {
        processorInfo* &altp = goodP[nPatches][nProxies][badForComm];	

#if USE_TOPOMAP 
	if(!altp)
	  altp = p;
	else {
	  //Find processors that are patch neighbors on the BGL torus
	  int neighbor = 0, neighbor_alt = 0;
	  
	  /*
	    if((tmgr->isNeighbor(altp->Id, patches[c->patch1].processor) ||
	    tmgr->isNeighbor(altp->Id, patches[c->patch2].processor)))
	    neighbor_alt = 1;
	    
	    if(tmgr->isNeighbor(p->Id, patches[c->patch1].processor) ||
	    tmgr->isNeighbor(p->Id, patches[c->patch2].processor))
	    neighbor = 1;
	  */
	  
	  if(tmgr.areNeighbors(altp->Id, patches[c->patch1].processor,
				    patches[c->patch2].processor, 4))
	    neighbor_alt = 1;
	  
	  if(tmgr.areNeighbors(p->Id, patches[c->patch1].processor, 
				    patches[c->patch2].processor, 4))
	    neighbor = 1;
	  
	  if(neighbor_alt == 1 && neighbor == 1) {
	    //Both are neighbors, only replace if lighter
	    if (p->load < altp->load ) {
	      altp = p;
	    }
	  }
	  else if(neighbor_alt == 0 && neighbor == 1)
	    //Previous was not a neighbor, kick him out
	    altp = p;
	  else if(neighbor_alt == 1 && neighbor == 0)
	    ;      //Give preference to good neighbors
	  else {
	    //Both not neighbors, choose nearby node to minimize hop bytes
	    /*
	      if (!altp || p->load < altp->load ) {
	      altp = p;
	      }
	    */

	    int alt_dist = 0, dist = 0;	    
	    int ax,ay,az, x,y,z, p1x,p1y,p1z, p2x,p2y,p2z;
	    
	    tmgr.rankToCoordinates(altp->Id, ax,ay,az);
	    tmgr.rankToCoordinates(p->Id, x,y,z);
	    tmgr.rankToCoordinates(patches[c->patch1].processor, p1x, p1y, p1z);
	    tmgr.rankToCoordinates(patches[c->patch2].processor, p2x, p2y, p2z);
 
	    alt_dist = abs(p1x - ax) + abs(p2x - ax) +
	      abs(p1y - ay) + abs(p1z - az) +
	      abs(p2y - ay) + abs(p2z - az);
	    
	    dist = abs(p1x - x) + abs(p2x - x) +
	      abs(p1y - y) + abs(p1z - z) +
	      abs(p2y - y) + abs(p2z - z);
	    
	    if(alt_dist > dist)
	      altp = p;	  
	  }
	}
#else 
        if (!altp || p->load < altp->load ) {	
	  altp = p;
        }
#endif	  
      }

      {
        processorInfo* &altp = poorP[nPatches][nProxies][badForComm];
        if (!altp || p->load < altp->load ) {
	  altp = p;
        }
      }
}

void Alg7::strategy()
{
  // double bestSize0, bestSize1, bestSize2;
  computeInfo *c;
  int numAssigned;
  processorInfo* goodP[3][3][2];  // goodP[# of real patches][# of proxies]
  processorInfo* poorP[3][3][2];  // fallback option

  double startTime = CmiWallTimer();

  // iout << iINFO << "calling makeHeaps. \n";
  adjustBackgroundLoadAndComputeAverage();
  makeHeaps();
  // iout << iINFO << "Before assignment\n" << endi;
  // printLoads();

  /*
  int numOverloaded = 0;
  for (int ip=0; ip<P; ip++) {
    if ( processors[ip].backgroundLoad > averageLoad ) ++numOverloaded;
  }
  if ( numOverloaded ) {
    iout << iWARN << numOverloaded
      << " processors are overloaded due to background load.\n" << endi;
  }
  */
	      
  numAssigned = 0;

  //   for (int i=0; i<numPatches; i++)
  //     { std::cout << "(" << patches[i].Id << "," << patches[i].processor ;}
  overLoad = 1.2;
  for (int ic=0; ic<numComputes; ic++) {

    // place computes w/ patches on heavily background loaded nodes first
    // place pair before self, because self is more flexible
    c = (computeInfo *) computeBgPairHeap->deleteMax();
    if ( ! c ) c = (computeInfo *) computeBgSelfHeap->deleteMax();
    if ( ! c ) c = (computeInfo *) computePairHeap->deleteMax();
    if ( ! c ) c = (computeInfo *) computeSelfHeap->deleteMax();

    if (c->processor != -1) continue; // skip to the next compute;

    if ( ! c ) NAMD_bug("Alg7: computesHeap empty!");
    int i,j,k;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++) {
        for(k=0;k<2;k++) {
	  goodP[i][j][k]=0;
	  poorP[i][j][k]=0;
        }
      }

    // first try for at least one proxy
    {
      Iterator nextProc;
      processorInfo *p;

      p = &processors[patches[c->patch1].processor];
      togrid(goodP, poorP, p, c);

      p = &processors[patches[c->patch2].processor];
      togrid(goodP, poorP, p, c);

      p = (processorInfo *)patches[c->patch1].
                            proxiesOn.iterator((Iterator *)&nextProc);
      while (p) {
        togrid(goodP, poorP, p, c);
        p = (processorInfo *)patches[c->patch1].
                            proxiesOn.next((Iterator*)&nextProc);
      }

      p = (processorInfo *)patches[c->patch2].
                            proxiesOn.iterator((Iterator *)&nextProc);
      while (p) {
        togrid(goodP, poorP, p, c);
        p = (processorInfo *)patches[c->patch2].
                            proxiesOn.next((Iterator*)&nextProc);
      }
      p = 0;
      // prefer to place compute with existing proxies over home patches
      if ((p = goodP[0][2][0])    // No home, two proxies
       || (p = goodP[1][1][0])    // One home, one proxy
       || (p = goodP[2][0][0])    // Two home, no proxies
       || (p = goodP[0][1][0])    // No home, one proxy
       || (p = goodP[1][0][0])    // One home, no proxies
       || (p = goodP[0][0][0])    // No home, no proxies
       || (p = goodP[0][1][1])    // No home, one proxy
       || (p = goodP[1][0][1])    // One home, no proxies
       || (p = goodP[0][0][1])    // No home, no proxies
         ) {
        assign(c,p); numAssigned++;
        continue;
      }
    }

    // no luck, do it the long way

    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    while (p) {
      togrid(goodP, poorP, p, c);
      p = (processorInfo *) pes->next(&nextProcessor);
    }

    //    if (numAssigned >= 0) {  Else is commented out below

    p = 0;
      // prefer to place compute with existing proxies over home patches
      if ((p = goodP[0][2][0])    // No home, two proxies
       || (p = goodP[1][1][0])    // One home, one proxy
       || (p = goodP[2][0][0])    // Two home, no proxies
       || (p = goodP[0][1][0])    // No home, one proxy
       || (p = goodP[1][0][0])    // One home, no proxies
       || (p = goodP[0][0][0])    // No home, no proxies
       || (p = goodP[0][1][1])    // No home, one proxy
       || (p = goodP[1][0][1])    // One home, no proxies
       || (p = goodP[0][0][1])    // No home, no proxies
       ) {
      assign(c,p); numAssigned++;
   } else if (   // overloaded processors
          (p = poorP[0][2][0])    // No home, two proxies
       || (p = poorP[1][1][0])    // One home, one proxy
       || (p = poorP[2][0][0])    // Two home, no proxies
       || (p = poorP[0][1][0])    // No home, one proxy
       || (p = poorP[1][0][0])    // One home, no proxies
       || (p = poorP[0][0][0])    // No home, no proxies
       || (p = poorP[0][1][1])    // No home, one proxy
       || (p = poorP[1][0][1])    // One home, no proxies
       || (p = poorP[0][0][1])    // No home, no proxies
       ) {
      //iout << iWARN << "overload assign to " << p->Id << "\n" << endi;
      assign(c,p); numAssigned++;
    } else {
      NAMD_bug("*** Alg 7 No receiver found 1 ***");
      break;
    }
  }

  printLoads();

  if ( computeMax() <= origMaxLoad ) {
    // binary-search refinement procedure
    multirefine(1.05);
    printLoads();
  }

}

