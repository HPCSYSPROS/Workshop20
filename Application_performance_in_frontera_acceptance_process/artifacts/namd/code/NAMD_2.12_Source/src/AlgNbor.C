/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "common.h"
#include "InfoStream.h"
#include "AlgNbor.h"

#define TINYLOAD 0.0005

AlgNbor::AlgNbor(int pe, computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes, int nNbs) :
  Rebalancer(computeArray, patchArray, 
	     processorArray, nComps, 
	     nPatches, nPes)
{
strategyName = "AlgNbor";
mype = pe;
nNbors = nNbs;
strategy();
}

void AlgNbor::strategy()
{
  int i, j;
  double startTime = CmiWallTimer();
  // CmiPrintf("strategy started on %d\n", CmiMyPe());

  // make initial assignment as it is now
  for (i=0; i<numComputes; i++) {
        assign((computeInfo *) &(computes[i]),
		 (processorInfo *) &(processors[computes[i].oldProcessor]));
  }

  // black-red mode balancing
  if (processors[mype].available == false) return;

  // calculate avarage load for neighbors
  double myload = processors[mype].load;
  double avgload = 0.0;
  for (i=0; i<P; i++) {
    if (processors[i].Id >= 0) {
//CmiPrintf("[%d] %d:%f\n", CmiMyPe(), i, processors[i].load);
      avgload += processors[i].load;
    }
  }
  avgload /= (nNbors+1);

  if (myload <= avgload) return;

  CmiPrintf("[%d]:Myload: %f, avrage load: %f. \n", mype, myload, avgload);

  IRSet *lightProcessors = new IRSet();
  for (i=0; i<P; i++) {
    if (processors[i].Id >= 0 &&  i != mype)
      if (processors[i].load < processors[mype].load)
	lightProcessors->insert((InfoRecord *) &(processors[i]));
  }

  int done = 0;
  while (processors[mype].load > avgload)
  {
    processorInfo* goodP[3][3];
    computeInfo* goodCompute[3][3];
    double goodSize[3][3];
    processorInfo* bestP;

    Iterator nextProcessor;
    if (lightProcessors->hasElements() == 0) break;
    processorInfo *p = (processorInfo *)lightProcessors->iterator((Iterator *) &nextProcessor);
    for(i=0; i < 3; i++)
        for(j=0; j<3; j++) {
          goodP[i][j] = 0;
          goodCompute[i][j] = 0;
          goodSize[i][j] = 0.;
        }
    while (p) 
    {
      Iterator nextCompute;
      nextCompute.id = 0;
      if (processors[mype].computeSet.hasElements() == 0) break;
      computeInfo *c = (computeInfo *)
		processors[mype].computeSet.iterator((Iterator *)&nextCompute);
      while (c) {
      
	if (c->load + p->load < processors[mype].load - c->load) {
          int nPatches, nProxies, badForComm;
          numAvailable(c,p,&nPatches,&nProxies,&badForComm);

	  if ( c->load > goodSize[nPatches][nProxies] ) {
	    goodSize[nPatches][nProxies] = c->load;
	    goodCompute[nPatches][nProxies] = c;
	    goodP[nPatches][nProxies] = p;
          }
	}

        nextCompute.id++;
        c = (computeInfo *) processors[mype].computeSet.
	  next((Iterator *)&nextCompute);
      }
      p = (processorInfo *)
	 lightProcessors->next((Iterator *) &nextProcessor);
    }
    if (goodCompute[2][0]) {
         deAssign(goodCompute[2][0], &processors[mype]);
         assign(goodCompute[2][0], goodP[2][0]);
         bestP = goodP[2][0];
    } else if (goodCompute[1][1]) {
         deAssign(goodCompute[1][1], &processors[mype]);
         assign(goodCompute[1][1], goodP[1][1]);
         bestP = goodP[1][1];
    } else if (goodCompute[0][2]) {
         deAssign(goodCompute[0][2], &processors[mype]);
         assign(goodCompute[0][2], goodP[0][2]);
         bestP = goodP[0][2];
    } else if (goodCompute[1][0]) {
         deAssign(goodCompute[1][0], &processors[mype]);
         assign(goodCompute[1][0], goodP[1][0]);
         bestP = goodP[1][0];
    } else if (goodCompute[0][1]) {
         deAssign(goodCompute[0][1], &processors[mype]);
         assign(goodCompute[0][1], goodP[0][1]);
         bestP = goodP[0][1];
    } else if (goodCompute[0][0]) {
         deAssign(goodCompute[0][0], &processors[mype]);
         assign(goodCompute[0][0], goodP[0][0]);
         bestP = goodP[0][0];
    } else {
         iout << iINFO << "AlgNbor: No receiver found" << "\n" << endi;
         break;
    }
    if (bestP->load > processors[mype].load) 
      lightProcessors->remove(bestP);
  }

/*
  do {
    processorInfo *p;
    computeInfo *c;

    p = (processorInfo *)procs.deleteMin();
    if (p == NULL) break;
    bool objfound = false;
    do {
      c = (computeInfo *)objs.deleteMax();
      if (c == NULL) break;

      double new_p_load = p->load + c->load;
      double my_new_load = myload - c->load;
//CmiPrintf("new load: new_p_load:%e my_new_load:%e c:%e\n", new_p_load, my_new_load, c->load);
      if (new_p_load < my_new_load) {
	objfound = true;
      } 
    } while (!objfound);
    if (!objfound) break;
    deAssign(c, &processors[mype]);
    assign(c, p);
    myload -= c->load;
    procs.insert(p);
  } while (myload > avgload);
  */
/*
  // double bestSize0, bestSize1, bestSize2;
  computeInfo *c;
  int numAssigned;
  processorInfo* goodP[3][3];  // goodP[# of real patches][# of proxies]
  processorInfo* poorP[3][3];  // fallback option


  //   iout << iINFO  << "calling makeHeaps. \n";
  computeAverage();
  makeHeaps();
  //   iout << iINFO
  //	<< "Before assignment\n" << endi;
  //   printLoads();
	      
  numAssigned = 0;

  //   for (int i=0; i<numPatches; i++)
  //     { std::cout << "(" << patches[i].Id << "," << patches[i].processor ;}
  overLoad = 1.2;
  for (int ic=0; ic<numComputes; ic++) {
    c = (computeInfo *) computesHeap->deleteMax();
    if ( ! c ) NAMD_bug("AlgNbor: computesHeap empty!");
    if (c->processor != -1) continue; // skip to the next compute;
    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    int i,j;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++) {
	goodP[i][j]=0;
	poorP[i][j]=0;
      }
    while (p) {
      int nPatches = numPatchesAvail(c,p);
      int nProxies = numProxiesAvail(c,p);
      if (nPatches < 0 || nPatches > 2)
	iout << iERROR << "Too many patches: " << nPatches << "\n" << endi;
      if (nProxies < 0 || nProxies > 2)
	iout << iERROR << "Too many proxies: " << nProxies << "\n" << endi;

      if (!goodP[nPatches][nProxies] ||
	    (p->load < goodP[nPatches][nProxies]->load)) {
        if (c->load + p->load < overLoad*averageLoad) {
	  goodP[nPatches][nProxies] = p;
        }
      }
      if (!poorP[nPatches][nProxies] ||
	    (p->load < poorP[nPatches][nProxies]->load)) {
	poorP[nPatches][nProxies] = p;   // fallback
      }
      p = (processorInfo *) pes->next(&nextProcessor);
    }

    //    if (numAssigned >= 0) {  Else is commented out below

    p = 0;
    if ((p = goodP[2][0])    // Two home, no proxies
     || (p = goodP[1][1])    // One home, one proxy
     || (p = goodP[0][2])    // No home, two proxies
     || (p = goodP[1][0])    // One home, no proxies
     || (p = goodP[0][1])    // No home, one proxy
     || (p = goodP[0][0])    // No home, no proxies
     || (p = poorP[2][0])    // Two home, no proxies, overload
     || (p = poorP[1][1])    // One home, one proxy, overload
     || (p = poorP[0][2])    // No home, two proxies, overload
     || (p = poorP[1][0])    // One home, no proxies, overload
     || (p = poorP[0][1])    // No home, one proxy, overload
     || (p = poorP[0][0])    // No home, no proxies, overload
       ) {
      assign(c,p); numAssigned++;
    } else {
      NAMD_bug("*** Alg 7 No receiver found 1 ***");
      break;
    }

  }


  // binary-search refinement procedure
  multirefine();
*/

  CmiPrintf("AlgNbor finish time: %f.\n", CmiWallTimer()-startTime);
}















