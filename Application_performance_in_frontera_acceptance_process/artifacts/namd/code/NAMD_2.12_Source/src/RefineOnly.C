/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "RefineOnly.h"

RefineOnly::RefineOnly(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
Rebalancer(computeArray, patchArray, 
	   processorArray, nComps, 
	   nPatches, nPes)
{
  strategyName = "Refine";
  strategy();
#if 0
    if(proxySendSpanning || proxyRecvSpanning) {
      decrSTLoad();
      computeAverage();
      createSpanningTree();
      incrSTLoad();
      for(int i=0; i<P; i++)
        delete [] processors[i].proxyUsage; 
      InitProxyUsage();
      multirefine(); 
      printLoads();
      //createSpanningTree();
    }
#endif

}

void RefineOnly::strategy()
{ 
  // iout << iINFO << "numComputes: " << numComputes << "\n" << endi;
  for (int i=0; i<numComputes; i++)
    assign((computeInfo *) &(computes[i]),
	   (processorInfo *) &(processors[computes[i].oldProcessor]));
	 
  // merge refinement code and move to Rebalancer
  multirefine();

  //  while (!refine())
  //    overLoad += .01;

  computeAverage();
#if 0
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO 
       << "After load balancing (predicted stats):\n" << endi;
#endif
  printLoads();
#if 0
       << "------------------------------------------------------------\n"
       << endi;
#endif
  double max = computeMax();

}
