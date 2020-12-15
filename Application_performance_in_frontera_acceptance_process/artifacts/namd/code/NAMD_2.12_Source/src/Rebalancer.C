/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/Rebalancer.C,v $
 * $Author: jim $
 * $Date: 2013/08/30 21:43:01 $
 * $Revision: 1.100 $
 *****************************************************************************/

#include "InfoStream.h"
#include "Node.h"
#include "Rebalancer.h"
#include "ProxyMgr.h"
#include "PatchMap.h"
#include "LdbCoordinator.h"
#include "memusage.h"
#include <iomanip>

#define ST_NODE_LOAD 		0.005
#define PROXY_LOAD              0.001
#define COMPUTE_LOAD            0.00005

Rebalancer::Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
      processorInfo *processorArray, int nComps, int nPatches, int nPes)
{
   bytesPerAtom = 32;
   strategyName = "None";
   computes = computeArray;
   patches =  patchArray;
   processors =  processorArray;
   numComputes = nComps;
   numPatches = nPatches;
   P = nPes;
   pes = NULL;
   computePairHeap = NULL;
   computeSelfHeap = NULL;
   computeBgPairHeap = NULL;
   computeBgSelfHeap = NULL;
   overLoad = 0.;
   numPesAvailable = 0;
   firstAssignInRefine = 0;
   collMsg = 0;

  const int beginGroup = processors[0].Id;
  const int endGroup = beginGroup + P;
#define INGROUP(PROC) ((PROC) >= beginGroup && (PROC) < endGroup)

   int i;
   int index;
   for (i=0; i<P; i++)
   {
      // For testing only...
      // processors[i].backgroundLoad = 0;
      // End of test section
      processors[i].load = processors[i].backgroundLoad;
      processors[i].computeLoad = 0;
      if (processors[i].available) {
        numPesAvailable += 1;
      }
   }

   for (i=0; i<nPatches; i++) {
     // Only for those patches which are in my group (hierarchical case)
     if INGROUP(patches[i].processor) {
       index = patches[i].processor - beginGroup;
       if (!patches[i].proxiesOn.find(&(processors[index]))) {
       	 patches[i].proxiesOn.unchecked_insert(&(processors[index]));
       	 processors[index].proxies.unchecked_insert(&(patches[i]));
       }
       processors[index].patchSet.unchecked_insert(&patches[i]);
     }
   }		          

   InitProxyUsage();

   for (i=0; i<numComputes; i++)
      computeArray[i].processor = -1;

   for (i=0; i < numComputes; i++) {
     // Only for those computes which are in my group (hierarchical case)
     if INGROUP(computes[i].oldProcessor) {
       index = computes[i].oldProcessor - beginGroup;
       processors[index].computeLoad += computes[i].load;
     }
   }

   // Added 4-29-98: Temporarily adds the compute load to the background
   // load so that the correct value for the total load can be displayed.
   float *temploads = new float[P];
   for(i=0; i<P; i++)
   {
      temploads[i] = processors[i].load;
      processors[i].load += processors[i].computeLoad;
   }

   origMaxLoad = computeMax();

   // iout << iINFO << "Initial load" << "\n";
   printLoads(1);

   for(i=0;i<P; i++)
   {
      processors[i].load = temploads[i];
      processors[i].computeLoad = 0;
   }
   
   delete [] temploads;

   // int count1=0, count2=0;
   // for (i=0; i<nPatches; i++)
   // {
   //    if (patches[i].proxiesOn.numElements() <= 1)
   //    count1++;
   //    else count2++;
   // }		          
   // iout << iINFO << "Count1 = " << count1
   //      << "Count2 = " << count2
   //      << "\n" << std::endl;
   // 
   // for (i=0; i <P; i++) 
   // {
   //    iout << iINFO << "\n proxies on proc. " << i << " are for patches:";
   //    processorArray[i].proxies->print();
   // }
   // 
   // iout << iINFO <<"\n" << endi;
   // strategy();

   // for (i=0; i<nPatches; i++)
   // {
   //    iout << "patch " << i << " on processor " << patches[i].processor << "\n" << endi;
   // }
}

Rebalancer::~Rebalancer()
{
  if ( computeMax() > origMaxLoad ) {
   if ( P == CkNumPes() ) {
   iout << "LDB:";
   if ( P != CkNumPes() ) {
     int w = 1;   
     int maxinw = 10;
     while ( maxinw < CkNumPes() ) {
       ++w;
       maxinw = 10*maxinw;
     }
     iout << " PES " <<
             std::setw(w) << std::right << processors[0].Id << "-" <<
             std::setw(w) << std::left  << processors[P-1].Id <<
             std::right;
   }
   iout << " Reverting to original mapping\n" << endi;
   fflush(stdout);
   } else {  // P != CkNumPes()
     if ( ! collMsg ) NAMD_bug("Rebalancer::~Rebalancer() collMsg null.");
     collMsg->finalAvgPeLoad = collMsg->initAvgPeLoad;
     collMsg->finalMaxPeLoad = collMsg->initMaxPeLoad;
     collMsg->finalTotalProxies = collMsg->initTotalProxies;
     collMsg->finalMaxPeProxies = collMsg->initMaxPeProxies;
     collMsg->finalMaxPatchProxies = collMsg->initMaxPatchProxies;
     collMsg->reverted = 1;
   }
   const int beginGroup = processors[0].Id;
   const int endGroup = beginGroup + P;
   for (int i=0; i < numComputes; i++) {
     // Only for those computes which are in my group (hierarchical case)
     if INGROUP(computes[i].oldProcessor) {
       computes[i].processor = computes[i].oldProcessor;
     }
   }
  }

  if ( P != CkNumPes() ) {
    if ( ! collMsg ) NAMD_bug("Rebalancer::~Rebalancer() collMsg null.");
    LdbCoordinator::Object()->sendCollectLoads(collMsg);
    collMsg = 0;
  }

  //for(int i=0; i<P; i++)
  //  delete [] processors[i].proxyUsage;
   delete pes;
   delete computePairHeap;
   delete computeSelfHeap;
   delete computeBgPairHeap;
   delete computeBgSelfHeap;
}

// Added 4-29-98: array proxyUsage on each processor keeps track of 
// how many computes are accessing each proxy on the processor.  If
// no computes are accessing it, the proxy can be removed in DeAssign
void Rebalancer::InitProxyUsage()
{
   int i;
   numProxies = 0;

   for(i=0; i<P; i++) {
     //processors[i].proxyUsage = new unsigned char [numPatches];
     //for(int j=0; j<numPatches; j++)
     //{
     //  processors[i].proxyUsage[j] = 0;
     //}

      Iterator nextCompute;
      nextCompute.id = 0;

      computeInfo *c = (computeInfo *)
         processors[i].computeSet.iterator((Iterator *)&nextCompute);

      while(c)
      {
	/* int n1 = */ //processors[i].proxyUsage[c->patch1]++;
	proxyUsage.increment (i, c->patch1); 
	/* int n2 = */ //processors[i].proxyUsage[c->patch2]++;
	proxyUsage.increment (i, c->patch2); 

         // iout << iINFO  
         // << "Assigning compute " << c->Id << " with work = " << c->load 
         // << " to processor " << processors[i].Id << "\n"
         // << "\tproxyUsage[" << c->patch1 << "]: " << n1 << " --> " << n1+1 << "\n";
         // << "\tproxyUsage[" << c->patch2 << "]: " << n2 << " --> " << n2+1 << "\n";
         // << std::endl;

         nextCompute.id++;
         c = (computeInfo *) processors[i].computeSet.next((Iterator *)&nextCompute);
      }
   }

  for (i=0; i<numPatches; i++)
  {
      numProxies += ( patches[i].proxiesOn.numElements() - 1 );
      Iterator nextProc;
      processorInfo *p = (processorInfo *)patches[i].proxiesOn.iterator((Iterator *)&nextProc);
      while (p) {
	  //p->proxyUsage[i] += 1;
	  proxyUsage.increment (p->Id, i);
          p = (processorInfo *)patches[i].proxiesOn.next((Iterator*)&nextProc);
      }
  }

}


void Rebalancer::strategy()
{
   iout << iINFO << "Strategy not implemented for the base class.\n" << "\n";
}

void Rebalancer::makeHeaps()
{
   int i, j;

   delete pes;
   pes = new minHeap(P+2);
   for (i=0; i<P; i++)
      pes->insert((InfoRecord *) &(processors[i]));

   delete computePairHeap;
   delete computeSelfHeap;
   delete computeBgPairHeap;
   delete computeBgSelfHeap;

   double bgLoadLimit = 0.5 * averageLoad;
   /*
   iout << iINFO << "Background load limit = " << bgLoadLimit << "\n";
   for (i=0; i<P; i++)
     if ( processors[i].backgroundLoad > bgLoadLimit )
       iout << iINFO << "Processor " << i << " background load = "
            << processors[i].backgroundLoad << "\n";
   iout << endi;
   */

   int numSelfComputes, numPairComputes, numBgSelfComputes, numBgPairComputes;

   while ( 1 ) {
    numSelfComputes = 0;
    numPairComputes = 0;
    numBgSelfComputes = 0;
    numBgPairComputes = 0;
    for (i=0; i<numComputes; i++) {
     int pa1 = computes[i].patch1;
     int pa2 = computes[i].patch2;
     if ( pa1 == pa2 ) {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit) {
          ++numBgSelfComputes;
        } else {
          ++numSelfComputes;
        }
     } else {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit
          || processors[patches[pa2].processor].backgroundLoad > bgLoadLimit) {
          ++numBgPairComputes;
        } else {
          ++numPairComputes;
        }
     }
    }

    int numBgComputes = numBgPairComputes + numBgSelfComputes;

    /*if ( numBgComputes ) {
        iout << iINFO << numBgComputes << " of " << numComputes
        << " computes have background load > " << bgLoadLimit << "\n" << endi;
    }*/

    if ( numBgComputes < 0.3 * numComputes ) break;
    else bgLoadLimit += 0.1 * averageLoad;
   }

   computePairHeap = new maxHeap(numPairComputes+2);
   computeSelfHeap = new maxHeap(numSelfComputes+2);
   computeBgPairHeap = new maxHeap(numBgPairComputes+2);
   computeBgSelfHeap = new maxHeap(numBgSelfComputes+2);

   for (i=0; i<numComputes; i++) {
     int pa1 = computes[i].patch1;
     int pa2 = computes[i].patch2;
     if ( pa1 == pa2 ) {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit) {
          computeBgSelfHeap->insert( (InfoRecord *) &(computes[i]));
        } else {
          computeSelfHeap->insert( (InfoRecord *) &(computes[i]));
        }
     } else {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit
          || processors[patches[pa2].processor].backgroundLoad > bgLoadLimit) {
          computeBgPairHeap->insert( (InfoRecord *) &(computes[i]));
        } else {
          computePairHeap->insert( (InfoRecord *) &(computes[i]));
        }
     }
   }

/*
   delete computePairHeap;
   delete computeSelfHeap;

   int numSelfComputes = 0;
   for (i=0; i<numComputes; i++)
      if ( computes[i].patch1 == computes[i].patch2 ) ++numSelfComputes;

   computeSelfHeap = new maxHeap(numSelfComputes+2);
   computePairHeap = new maxHeap(numComputes-numSelfComputes+2);

   for (i=0; i<numComputes; i++)
      if ( computes[i].patch1 == computes[i].patch2 )
         computeSelfHeap->insert( (InfoRecord *) &(computes[i]));
      else
         computePairHeap->insert( (InfoRecord *) &(computes[i]));
*/
}

void Rebalancer::makeTwoHeaps()
{
   int i, j;

   delete pes;
   pes = new minHeap(P+2);
   for (i=0; i<P; i++)
      pes->insert((InfoRecord *) &(processors[i]));

   delete computePairHeap;
   delete computeSelfHeap;
   delete computeBgPairHeap;
   delete computeBgSelfHeap;

   int numSelfComputes, numPairComputes;

   numSelfComputes = 0;
   numPairComputes = 0;
   for (i=0; i<numComputes; i++) {
     int pa1 = computes[i].patch1;
     int pa2 = computes[i].patch2;
     if (pa1 == pa2)
       ++numSelfComputes;
     else
       ++numPairComputes;
   }

   computePairHeap = new maxHeap(numPairComputes+2);
   computeSelfHeap = new maxHeap(numSelfComputes+2);

   for (i=0; i<numComputes; i++) {
     int pa1 = computes[i].patch1;
     int pa2 = computes[i].patch2;
     if ( pa1 == pa2 )
       computeSelfHeap->insert( (InfoRecord *) &(computes[i]));
     else
       computePairHeap->insert( (InfoRecord *) &(computes[i]));
   }
}

// not safe with hybrid balancer
//void Rebalancer::assign(computeInfo *c, int processor)
//{
//   assign(c, &(processors[processor]));
//}

void Rebalancer::assign(computeInfo *c, processorInfo *p)
{
   c->processor = p->Id;
   p->computeSet.unchecked_insert((InfoRecord *) c);
#if COMPUTE_CORRECTION
   if(firstAssignInRefine)
     p->computeLoad += c->load + COMPUTE_LOAD;
   else
#endif
     p->computeLoad += c->load;
     
   p->load = p->computeLoad + p->backgroundLoad;
   patchInfo* patch1 = (patchInfo *) &(patches[c->patch1]);
   patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);

   if (!patch1->proxiesOn.find(p)) {
     p->proxies.unchecked_insert(patch1); 
     patch1->proxiesOn.unchecked_insert(p); 
     numProxies++;
#if PROXY_CORRECTION
     if(firstAssignInRefine) {
       processors[p->Id].load += PROXY_LOAD;
       processors[p->Id].backgroundLoad += PROXY_LOAD;
     }
#endif
   }

   if (!patch2->proxiesOn.find(p)) {
     p->proxies.unchecked_insert(patch2); 
     patch2->proxiesOn.unchecked_insert(p);
     numProxies++;
#if PROXY_CORRECTION
     if(firstAssignInRefine) {
       processors[p->Id].load += PROXY_LOAD;
       processors[p->Id].backgroundLoad += PROXY_LOAD;
     }
#endif
   }
   
   // 4-29-98: Added the following code to keep track of how many proxies
   // on each processor are being used by a compute on that processor
   /* int n1 = */ //p->proxyUsage[c->patch1]++;
   proxyUsage.increment (p->Id, c->patch1);
   /* int n2 = */ //p->proxyUsage[c->patch2]++;
   proxyUsage.increment (p->Id, c->patch2);

   // iout << iINFO  
   // << "Assigning compute " << c->Id << " with work = " << c->load 
   // << " to processor " << p->Id << "\n"
   // << "\tproxyUsage[" << c->patch1 << "]: " << n1 << " --> " << n1+1 << "\n"
   // << "\tproxyUsage[" << c->patch2 << "]: " << n2 << " --> " << n2+1 << "\n"
   // << std::endl;

#if 0
   iout << "Assign " << c->Id << " patches " << c->patch1 << " " << c->patch2
        << " load " << c->load << " to " << p->Id << " new load "
        << p->load << " background " << p->backgroundLoad
        << " nPatches " << nPatches << " nProxies " << nProxies;
   if ( nPatches + nProxies < 2 ) iout << " addProxy";
   if ( badForComm ) iout << " badForComm";
   iout << "\n" << endi;
#endif
}

void  Rebalancer::deAssign(computeInfo *c, processorInfo *p)
{
   if (!p->computeSet.remove(c))  {
      iout << iINFO << "ERROR: Rebalancer tried to deAssign an object that is not on the processor.\n" << endi;
      return;
   }

   double temp_load = 0.0;

   c->processor = -1;
   p->computeLoad -= c->load;
   CmiAssert(p->computeLoad >= 0.0);
   temp_load = p->load - c->load;
   p->load = p->computeLoad + p->backgroundLoad;
   CmiAssert( fabs(temp_load - p->load) < 0.001 );

   // 4-29-98: Added the following code to keep track of how many proxies 
   // on each processor are being used by a compute on that processor.
   // If no computes are using the proxy, it should be removed if it is not
   // on the processor that its patch is on.
   /* int n1 = */ //p->proxyUsage[c->patch1]--;
   proxyUsage.decrement (p->Id, c->patch1);
   /* int n2 = */ //p->proxyUsage[c->patch2]--;
   proxyUsage.decrement (p->Id, c->patch2);

   // iout << iINFO
   // << "De-assigning compute " << c->Id << " from processor " << p->Id << "\n"
   // << "\tproxyUsage[" << c->patch1 << "]: " << n1 << " --> " << n1-1 << "\n"
   // << "\tproxyUsage[" << c->patch2 << "]: " << n2 << " --> " << n2-1 << "\n"
   // << std::endl;

   //if(p->proxyUsage[c->patch1] <= 0 && p->Id != patches[c->patch1].processor)
   if(proxyUsage.getVal(p->Id, c->patch1) <= 0 && p->Id != patches[c->patch1].processor)
   {
      // iout << iINFO 
      // << "REMOVING PROXY " << c->patch1 << " FROM PROCESSOR " << p->Id 
      // << std::endl << endl;

      patchInfo* patch1 = (patchInfo *) &(patches[c->patch1]);
      p->proxies.remove(patch1);
      patch1->proxiesOn.remove(p);
      numProxies--;
#if PROXY_CORRECTION
      if(firstAssignInRefine) {
	processors[p->Id].load -= PROXY_LOAD;
	processors[p->Id].backgroundLoad -= PROXY_LOAD;
	if(processors[p->Id].backgroundLoad < 0.0) {
	  processors[p->Id].backgroundLoad = 0.0;
	  processors[p->Id].load += PROXY_LOAD;
	}
      }
#endif
   }
   
   //if(p->proxyUsage[c->patch2] <= 0 && p->Id != patches[c->patch2].processor)
   if(proxyUsage.getVal(p->Id, c->patch2) <= 0 && p->Id != patches[c->patch2].processor)
   {
      // iout << iINFO
      // << "REMOVING PROXY " << c->patch1 << " FROM PROCESSOR " << p->Id 
      // << std::endl << endl;

      patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);
      p->proxies.remove(patch2);
      patch2->proxiesOn.remove(p);
      numProxies--;
#if PROXY_CORRECTION
      if(firstAssignInRefine) {
	processors[p->Id].load -= PROXY_LOAD;
	processors[p->Id].backgroundLoad -= PROXY_LOAD;
	if(processors[p->Id].backgroundLoad < 0.0) {
	  processors[p->Id].backgroundLoad = 0.0;
	  processors[p->Id].load += PROXY_LOAD;
	}
      }
#endif
   }
}

void Rebalancer::refine_togrid(pcgrid &grid, double thresholdLoad,
			processorInfo *p, computeInfo *c) {

  if(p->available == false) return;

  if ( c->load + p->load < thresholdLoad) {
    int nPatches, nProxies, badForComm;
    numAvailable(c,p,&nPatches,&nProxies,&badForComm);

    // if ( badForComm ) return;

    pcpair *pair = &grid[nPatches][nProxies][badForComm];

    if (! pair->c) {
      pair->c = c;
      pair->p = p;
    } else {
      double newval = p->load - c->load;
      if ( c->load + p->load < averageLoad ) {
         newval -= averageLoad;
      }
      double oldval = pair->p->load - pair->c->load;
      if ( pair->c->load + pair->p->load < averageLoad ) {
         oldval -= averageLoad;
      }
      if (newval < oldval) {
	pair->c = c;
	pair->p = p;
      }
    }
  }
}

int Rebalancer::refine()
{
   int finish = 1;
   int no_new_proxies = 0;  // set to true if new proxies are futile
   maxHeap *heavyProcessors = new maxHeap(P);

   IRSet *lightProcessors = new IRSet();
   int i;
   double thresholdLoad = overLoad * averageLoad;
   for (i=0; i<P; i++)
   {
      // iout << iINFO << "\n Computes on processor " << i << " ";
      // processors[i].computeSet->print();
      // iout << iINFO << "\n" << endi;
      if (processors[i].load > thresholdLoad)
         heavyProcessors->insert((InfoRecord *) &(processors[i]));
      else lightProcessors->insert((InfoRecord *) &(processors[i]));
   }

#if LDB_DEBUG
   iout << "\nBefore Refinement Summary" << "\n";
   printSummary();
#endif

   int done = 0;
   while (!done)
   {
      // processorInfo *donor = (processorInfo *) heavyProcessors->deleteMax();
      /* Keep selecting new donors, until we find one with some compute to
       * migrate
       */
/*
      computeInfo* c=0;
      while (donor && !c) {
        Iterator nextCompute;
        nextCompute.id = 0;
        c = (computeInfo *) donor->
            computeSet.iterator((Iterator *)&nextCompute);
        if (!c) {
          iout << iINFO << "Ignoring donor " << donor->Id
               << " because no computes\n" << endi;
	  donor = (processorInfo*)heavyProcessors->deleteMax();
        }
      };
*/

      processorInfo *donor;
      while (donor = (processorInfo*)heavyProcessors->deleteMax()) {
	if (donor->computeSet.hasElements()) break;
        if ( ! no_new_proxies ) {
          /*
          iout << iINFO << "Most-loaded processor " << donor->Id
               << " (" << donor->patchSet.numElements() << " patches, "
               << donor->proxies.numElements() << " proxies)"
               << " has no migratable work.\n" << endi;
          */
          no_new_proxies = 1;  // New proxies would not improve load balance.
        }
      }
  
      if (!donor) break;  // No donors found at all! Give up 

      pcgrid grid;
#define REASSIGN(GRID) if (GRID.c) { \
	   deAssign(GRID.c, donor); \
           assign(GRID.c, GRID.p); \
           bestP = GRID.p; \
        }

      // try for at least one proxy
      {
        Iterator nextCompute;
        nextCompute.id = 0;
        computeInfo *c = (computeInfo *) 
           donor->computeSet.iterator((Iterator *)&nextCompute);
        while (c)
        {
          Iterator nextProc;
          processorInfo *p;

	  p = &processors[patches[c->patch1].processor];
	  refine_togrid(grid, thresholdLoad, p, c);

	  if (c->patch1 != c->patch2)
	  {
	  p = &processors[patches[c->patch2].processor];
	  refine_togrid(grid, thresholdLoad, p, c);
	  }

	  p = (processorInfo *)patches[c->patch1].
				proxiesOn.iterator((Iterator *)&nextProc);
          while (p) {
	    refine_togrid(grid, thresholdLoad, p, c);
            p = (processorInfo *)patches[c->patch1].
				proxiesOn.next((Iterator*)&nextProc);
          }

	  if (c->patch1 != c->patch2) 
	  {
	  p = (processorInfo *)patches[c->patch2].
				proxiesOn.iterator((Iterator *)&nextProc);
          while (p) {
	    refine_togrid(grid, thresholdLoad, p, c);
            p = (processorInfo *)patches[c->patch2].
				proxiesOn.next((Iterator*)&nextProc);
          }
	  }

          nextCompute.id++;
          c = (computeInfo *) donor->computeSet.
	    next((Iterator *)&nextCompute);
        }
        processorInfo* bestP = 0;
        // prefer proxies to home patches
	REASSIGN(grid[0][2][0])
	else REASSIGN(grid[1][1][0])
	else REASSIGN(grid[2][0][0])
        else if ( no_new_proxies ) { finish = 0; break; }
	else REASSIGN(grid[0][1][0])
	else REASSIGN(grid[1][0][0])
	else REASSIGN(grid[0][0][0])
	// else REASSIGN(grid[0][1][1])
	// else REASSIGN(grid[1][0][1])
	// else REASSIGN(grid[0][0][1])
        if (bestP) {
	  if (bestP->load > averageLoad) lightProcessors->remove(bestP);
	  if (donor->load > thresholdLoad)
		heavyProcessors->insert((InfoRecord *) donor);
	  else lightProcessors->insert((InfoRecord *) donor);
	  continue;
        }
      }

      if ( no_new_proxies ) iout << iINFO
         << "ERROR: Rebalancer::refine() algorithm is broken.\n" << endi;

      // no luck, do it the long way

      //find the best pair (c,receiver)
      Iterator nextProcessor;
      processorInfo *p = (processorInfo *) 
      lightProcessors->iterator((Iterator *) &nextProcessor);

      while (p)
      {
         Iterator nextCompute;
         nextCompute.id = 0;
         computeInfo *c = (computeInfo *) 
            donor->computeSet.iterator((Iterator *)&nextCompute);
         while (c)
         {
#if USE_TOPOMAP
           int flag = tmgr.areNeighbors(p->Id, patches[c->patch1].processor, 
				     patches[c->patch2].processor, 8);
	   if(flag)
#endif
	     {	     
	       refine_togrid(grid, thresholdLoad, p, c);
	     }
	   nextCompute.id++;
	   c = (computeInfo *) donor->computeSet.
	     next((Iterator *)&nextCompute);
         }
         p = (processorInfo *) 
	   lightProcessors->next((Iterator *) &nextProcessor);
      }

      //we have narrowed the choice to 6 candidates.
      // prefer proxies to home patches
      {
        processorInfo* bestP = 0;
	REASSIGN(grid[0][2][0])
	else REASSIGN(grid[1][1][0])
	else REASSIGN(grid[2][0][0])
	else REASSIGN(grid[0][1][0])
	else REASSIGN(grid[1][0][0])
	else REASSIGN(grid[0][0][0])
	// else REASSIGN(grid[0][1][1])
	// else REASSIGN(grid[1][0][1])
	// else REASSIGN(grid[0][0][1])
	else { finish = 0; break; }
	if (bestP->load > averageLoad) lightProcessors->remove(bestP);
	if (donor->load > thresholdLoad)
		heavyProcessors->insert((InfoRecord *) donor);
	else lightProcessors->insert((InfoRecord *) donor);
      }

   }  

#if LDB_DEBUG
   iout << "After Refinement Summary" << "\n";
   printSummary();

   if (!finish) {
     iout << iINFO << "Refine: No solution found for overLoad = " 
	  << overLoad << "\n" << endi;
   }
#endif

   delete heavyProcessors;
   delete lightProcessors;

   return finish;
}

// this binary search refinement procedure assume you already assigned computes
// to their processors before calling this!!
void Rebalancer::multirefine(double overload_start)
{
  // The New refinement procedure.  This is identical to the code in
  // RefineOnly.C, and probably should be merged with that code to form
  // a binary-search function

  double avg = computeAverage();
  double max = computeMax();

#if LDB_DEBUG
  iout << "******** Processors with background load > average load ********" << "\n";
#endif

  int numOverloaded = 0;
  for (int ip=0; ip<P; ip++) {
    if ( processors[ip].backgroundLoad > averageLoad ) {
      ++numOverloaded;
#if LDB_DEBUG
      iout << iINFO << "Info about proc " << ip << ": Load: " << processors[ip].load << " Bg Load: " << processors[ip].backgroundLoad << " Compute Load: " << processors[ip].computeLoad << " No of computes: " << processors[ip].computeSet.numElements() << "\n";
#endif
    }
  }
  if ( numOverloaded ) {
    iout << iWARN << numOverloaded
      << " processors are overloaded due to high background load.\n" << endi;
  }
#if LDB_DEBUG
  iout << "******** Processor List Ends ********" << "\n\n";
#endif

  const double overloadStep = 0.01;
  const double overloadStart = overload_start;       //1.05;
  double dCurOverload = max / avg;
  
  int minOverload = 0;   //Min overload should be 1.05 ?
  int maxOverload = (int)((dCurOverload - overloadStart)/overloadStep + 1);
  double dMinOverload = minOverload * overloadStep + overloadStart;
  double dMaxOverload = maxOverload * overloadStep + overloadStart;

#if LDB_DEBUG 
  iout << iINFO
       << "Balancing from " << minOverload << " = " << dMinOverload 
       << " to " << maxOverload << "=" << dMaxOverload 
       << " dCurOverload=" << dCurOverload << " max=" << max << " avg=" << avg
       << "\n" << endi;
#endif

  int curOverload;
  int refineDone = 0;

  overLoad = dMinOverload;
  if (refine())
    refineDone = 1;
  else {
    overLoad = dMaxOverload;
    if (!refine()) {
      iout << iINFO << "ERROR: Could not refine at max overload\n" << endi;
      refineDone = 1;
    }
  }

  // Scan up, until we find a refine that works
  while (!refineDone) {
    if (maxOverload - minOverload <= 1)
      refineDone = 1;
    else {
      curOverload = (maxOverload + minOverload ) / 2;

      overLoad = curOverload * overloadStep + overloadStart;
#if LDB_DEBUG 
      iout << iINFO << "Testing curOverload " << curOverload 
	   << "=" << overLoad << " [min,max]=" 
	   << minOverload << ", " << maxOverload
	   << "\n" << endi;
#endif     
      if (refine())
	maxOverload = curOverload;
      else
	minOverload = curOverload;
    }
  }

}

void Rebalancer::printResults()
{
  iout << iINFO << "ready to print result \n" << "\n";
}


void Rebalancer::printLoads(int phase)  // 0=nocollective, 1=initial, 2=proxies, 3=final
{

   int i, total = 0, numBytes = 0;
   double max;
   int maxproxies = 0;
   int maxpatchproxies = 0;
   double avgBgLoad =0.0;

   for (i=0; i<P; i++) {
      int nproxies = processors[i].proxies.numElements() - 
			processors[i].patchSet.numElements();
      total += nproxies;
      if ( nproxies > maxproxies ) maxproxies = nproxies;
      avgBgLoad += processors[i].backgroundLoad;
      Iterator p;
      int count = 0;
    
      patchInfo *patch = (patchInfo *) processors[i].patchSet.iterator(&p);
      while (patch)
      {
         int myProxies;
         myProxies = patch->proxiesOn.numElements()-1;
         if ( myProxies > maxpatchproxies ) maxpatchproxies = myProxies;
         numBytes += myProxies *patch->numAtoms*bytesPerAtom;
         count += myProxies;
         patch = (patchInfo *)processors[i].patchSet.next(&p);
      }
   }

   avgBgLoad /= P;
   computeAverage();
   max = computeMax();

  if ( P == CkNumPes() ) {
   iout << "LDB:";
   if ( P != CkNumPes() ) {
     int w = 1;   
     int maxinw = 10;
     while ( maxinw < CkNumPes() ) {
       ++w;
       maxinw = 10*maxinw;
     }
     iout << " PES " <<
             std::setw(w) << std::right << processors[0].Id << "-" <<
             std::setw(w) << std::left  << processors[P-1].Id <<
             std::right;
   }
   iout << " TIME " << CmiWallTimer() << " LOAD: AVG " << averageLoad 
     << " MAX " << max << "  PROXIES: TOTAL " << total << " MAXPE " << 
     maxproxies << " MAXPATCH " << maxpatchproxies << " " << strategyName 
     << " MEM: " << memusage_MB() << " MB\n" << endi;
   fflush(stdout);
  }

   if ( P != CkNumPes() ) {  // collect stats on pe 0
     switch ( phase ) {
     case 0:  // no collective
       NAMD_bug("Rebalancer::printLoads(0) called with hybrid balancer.");
     break;
     case 1:  // initial
       if ( collMsg ) NAMD_bug("Rebalancer::printLoads(1) collMsg not null.");
       collMsg = new CollectLoadsMsg;
       collMsg->reverted = 0;
       collMsg->firstPe = processors[0].Id;
       collMsg->lastPe = processors[P-1].Id;
       collMsg->initTime = CmiWallTimer();
       collMsg->initMemory = memusage_MB();
       collMsg->initAvgPeLoad = averageLoad;
       collMsg->initMaxPeLoad = max;
       collMsg->initTotalProxies = total;
       collMsg->initMaxPeProxies = maxproxies;
       collMsg->initMaxPatchProxies = maxpatchproxies;
     break;
     case 2:  // proxies (optional)
       if ( ! collMsg ) NAMD_bug("Rebalancer::printLoads(2) collMsg null.");
       collMsg->initAvgPeLoad = averageLoad;
       collMsg->initMaxPeLoad = max;
       collMsg->initTotalProxies = total;
       collMsg->initMaxPeProxies = maxproxies;
       collMsg->initMaxPatchProxies = maxpatchproxies;
     break;
     case 3:  // final
       if ( ! collMsg ) NAMD_bug("Rebalancer::printLoads(3) collMsg null.");
       collMsg->finalTime = CmiWallTimer();
       collMsg->finalMemory = memusage_MB();
       collMsg->finalAvgPeLoad = averageLoad;
       collMsg->finalMaxPeLoad = max;
       collMsg->finalTotalProxies = total;
       collMsg->finalMaxPeProxies = maxproxies;
       collMsg->finalMaxPatchProxies = maxpatchproxies;
       strncpy(collMsg->strategyName,strategyName,15);
       collMsg->strategyName[15] = 0;
     break;
     default:
       NAMD_bug("Rebalancer::printLoads() called with unknown phase.");
     }
   }

}

void Rebalancer::printSummary()
{
   int i;
   // After refining, compute min, max and avg processor load
   double total = processors[0].load;
   double min = processors[0].load;
   int min_proc = 0;
   double max = processors[0].load;
   int max_proc = 0;
   for (i=1; i<P; i++) {
     total += processors[i].load;
     if (processors[i].load < min) {
       min = processors[i].load;
       min_proc = i;
     }
     if (processors[i].load > max) {
       max = processors[i].load;
       max_proc = i;
     }
   }
   iout << iINFO << "  min = " << min << " processor " << min_proc << "\n";
   iout << iINFO << "  max = " << max << " processor " << max_proc << "\n";
   iout << iINFO << "  total = " << total << " average = " << total/P << "\n";
   iout << iINFO << "Info about most overloaded processor " << max_proc << ": Load: " << processors[max_proc].load << " Bg Load: " << processors[max_proc].backgroundLoad << " Compute Load: " << processors[max_proc].computeLoad << " No of computes: " << processors[max_proc].computeSet.numElements() << " No. of proxies: " << processors[max_proc].proxies.numElements() << "\n" << endi; 
}

double Rebalancer::computeAverage()
{
   int i;
   double total = 0.;
   for (i=0; i<numComputes; i++)
      total += computes[i].load;

   for (i=0; i<P; i++) {
      if (processors[i].available) {
        total += processors[i].backgroundLoad;
      }
   }
  
   if (numPesAvailable == 0) {
     CmiPrintf("Warning: no processors available for load balancing!\n");
     averageLoad = 0.0;
   }
   else 
     averageLoad = total/numPesAvailable;
   return averageLoad;
}

void Rebalancer::adjustBackgroundLoadAndComputeAverage()
{
  // useful for AlgSeven when some loads start out as zero

   if (numPesAvailable == 0) {
     computeAverage();  // because otherwise someone will forget
     return;
   }
  
   int i;
   double bgtotal = 0.;
   for (i=0; i<P; i++) {
      if (processors[i].available) {
        bgtotal += processors[i].backgroundLoad;
      }
   }
   double bgavg = bgtotal / numPesAvailable;

   int nadjusted = 0;
   for (i=0; i<P; i++) {
      if (processors[i].available) {
        double bgload = processors[i].backgroundLoad;
        if ( bgload < bgavg ) {
          processors[i].backgroundLoad = bgavg;
          ++nadjusted;
        }
      }
   }
   // iout << iINFO << "Adjusted background load on " << nadjusted
   //     << " nodes.\n" << endi;

   computeAverage();  // because otherwise someone will forget
}

double Rebalancer::computeMax()
{
   int i;
   double max = processors[0].load;
   for (i=1; i<P; i++)
   {
      if (processors[i].load > max) 
         max = processors[i].load;
   }
   return max;
}

int Rebalancer::isAvailableOn(patchInfo *patch, processorInfo *p)
{
   return  patch->proxiesOn.find(p);
}

void Rebalancer::numAvailable(computeInfo *c, processorInfo *p,
           int *nPatches, int *nProxies, int *isBadForCommunication)
{
   // return the number of proxy/home patches available on p for c (0, 1, 2)
   int realPe, index;
   int patch_count = 0;
   int proxy_count = 0;

  const int beginGroup = processors[0].Id;
  const int endGroup = beginGroup + P;

   patchInfo &pa1 = patches[c->patch1];
   patchInfo &pa2 = patches[c->patch2];
   int pa1_avail = 1;
   int pa2_avail = 1;

   if (pa1.processor == p->Id) {
     patch_count++;
   } else if ( pa1.proxiesOn.find(p) ) {
     proxy_count++;
   } else {
     pa1_avail = 0;
   }

   // self computes get one patch for free here
   if (c->patch1 == c->patch2 || pa2.processor == p->Id) {
     patch_count++;
   } else if ( pa2.proxiesOn.find(p) ) {
     proxy_count++;
   } else {
     pa2_avail = 0;
   }

   *nPatches = patch_count;
   *nProxies = proxy_count;

  if ( isBadForCommunication ) {  // skip work if pointer is null
   int bad = 0;

   if ( patch_count + proxy_count < 2 ) {
     double bgLoadLimit = 1.2 * averageLoad;
     if ( p->backgroundLoad > bgLoadLimit ) bad = 1;
     else {
       int proxiesPerPeLimit = numProxies / numPesAvailable + 3;
       if ( proxiesPerPeLimit < 6 ) proxiesPerPeLimit = 6;

       if ( p->proxies.numElements() > proxiesPerPeLimit ) bad = 1;

       int proxiesPerPatchLimit = numProxies / numPatches + 3;
       if ( proxiesPerPatchLimit < 6 ) proxiesPerPatchLimit = 6;

       if ( ! bad && ! pa1_avail ) {
	 // HYBRID check for range in local group
	 realPe = pa1.processor;
	 if INGROUP(realPe) {
	   index = realPe - beginGroup;
           //BACKUP if ( processors[pa1.processor].backgroundLoad > bgLoadLimit) bad = 1;
           if (processors[index].backgroundLoad > bgLoadLimit) bad = 1;
           else if ( pa1.proxiesOn.numElements() > proxiesPerPatchLimit ) bad = 1;
	 } else bad = 1;  // patch has proxies we don't know about
       }

       if ( ! bad && ! pa2_avail ) {
	 // HYBRID check for range in local group
	 realPe = pa2.processor;
	 if INGROUP(realPe) {
	   index = realPe - beginGroup;
	   // BACKUP if ( processors[pa2.processor].backgroundLoad > bgLoadLimit) bad = 1;
	   if ( processors[index].backgroundLoad > bgLoadLimit) bad = 1;
	   else if ( pa2.proxiesOn.numElements() > proxiesPerPatchLimit ) bad = 1;
	 } else bad = 1;  // patch has proxies we don't know about
       }

     }
   }

   *isBadForCommunication = bad;
  }
}

void Rebalancer::createSpanningTree() {
  ProxyTree &pt = ProxyMgr::Object()->getPtree();
  Iterator nextP;
  processorInfo *p;
#ifndef NODEAWARE_PROXY_SPANNINGTREE
  if(pt.sizes==NULL)
    pt.sizes = new int[numPatches];
#endif
 
  if (pt.proxylist == NULL)
    pt.proxylist = new NodeIDList[numPatches];
  for(int i=0; i<numPatches; i++)
  {
    pt.proxylist[i].resize(patches[i].proxiesOn.numElements());
    nextP.id = 0;
    p = (processorInfo *)(patches[i].proxiesOn.iterator((Iterator *)&nextP));
    int j = 0;
    while(p) {
      //if (p->Id < 0)
      //  printf ("Inserting proxy on -ve processor %d for patch %d\n", p->Id, i);

      if (p->Id == (PatchMap::Object()->node(i))) {
        p = (processorInfo *)(patches[i].proxiesOn.next((Iterator *)&nextP));
        continue;
      }

      pt.proxylist[i][j] = p->Id;
      nextP.id++;
      p = (processorInfo *)(patches[i].proxiesOn.next((Iterator *)&nextP));
      j++;
    }
    pt.proxylist[i].resize(j);
  }
  CkPrintf("Done intialising\n");
#ifdef NODEAWARE_PROXY_SPANNINGTREE
  ProxyMgr::Object()->buildNodeAwareSpanningTree0();
#else
  ProxyMgr::Object()->buildSpanningTree0();
#endif
}

void Rebalancer::decrSTLoad() {
  int pe;
  ProxyTree &pt = ProxyMgr::Object()->getPtree();
  for(int i=0; i<numPatches; i++)
    for(int j=1; j<pt.proxylist[i].size() && j<proxySpanDim; j++) {
      pe = pt.proxylist[i][j];
      processors[pe].load -= ST_NODE_LOAD;
      processors[pe].backgroundLoad -= ST_NODE_LOAD;
      if(processors[pe].load < 0.0)
        processors[pe].load = 0.0;
      if(processors[pe].backgroundLoad < 0.0)
        processors[pe].backgroundLoad = 0.0;
    }
}

void Rebalancer::incrSTLoad() {
  int pe;
  ProxyTree &pt = ProxyMgr::Object()->getPtree();
  for(int i=0; i<numPatches; i++)
    for(int j=1; j<pt.proxylist[i].size() && j<proxySpanDim; j++) {
      pe = pt.proxylist[i][j];
      processors[pe].load += ST_NODE_LOAD;
      processors[pe].backgroundLoad += ST_NODE_LOAD;
    }
}

