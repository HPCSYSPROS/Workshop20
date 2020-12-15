/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/Rebalancer.h,v $
 * $Author: jim $
 * $Date: 2013/08/30 18:18:20 $
 * $Revision: 1.35 $
 *****************************************************************************/

/** \file Rebalancer.h
 *  This class is a super class for all the load balancing classes 
 *  which define new strategies. It has common functions like makeHeaps,
 *  makeTwoHeaps, refine_togrid and others which are used by different
 *  LDBs.
 *
 *  NOTE: Confusingly, this class also houses the strategy for RefineOnly
 *  in the function "refine".
 */

#ifndef REFINEONLY_DEFS_H
#define REFINEONLY_DEFS_H

#include "elements.h"
#include "heap.h"
#if USE_TOPOMAP 
#include "TopoManager.h"
#endif
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"

#define LDB_DEBUG		0	// for verbose LDB output
#define PROXY_CORRECTION        0
#define COMPUTE_CORRECTION	0

#include "ckhashtable.h"

class ProxyUsageKey {
 protected:
  int      processor;
  int      patch;
  
 public:
  ProxyUsageKey (int pe, int patch) {
    this->processor = pe;
    this->patch     = patch;
  }

  CkHashCode hash (void) const {
    return (patch << 16) + processor;
  }

  static CkHashCode  staticHash (const void *x, size_t size) {
    return ((ProxyUsageKey *)x)->hash();
  }

  int compare (const ProxyUsageKey &in) const {
    if ((in.patch == patch) && (in.processor == processor))
      return 1;
    
    return 0;
  }
   
  static int staticCompare (const void *a, const void *b, size_t size) {
    return ((ProxyUsageKey *)a)->compare(* (ProxyUsageKey *)b);
  }
};

class ProxyUsage {
 protected:
  CkHashtableT <ProxyUsageKey, int>  htable;
  
 public:
  
  ProxyUsage () : htable (1217, 0.5) {}   //pass in a large prime close to 
                                          //1k processors

  void increment (int pe, int patch) {
    ProxyUsageKey  pkey (pe, patch);

    int val = htable.get (pkey);
    htable.put (pkey) =  val + 1;      
  }

  void decrement (int pe, int patch) {
    ProxyUsageKey  pkey (pe, patch);
    
    int val = htable.get (pkey);
    CkAssert (val > 0);
    val --;

    if (val == 0)
      htable.remove (pkey);
    else 
      htable.put (pkey) = val; 
  }
     
  int getVal (int pe, int patch) {
    ProxyUsageKey  pkey (pe, patch);  
    return htable.get (pkey);
  }
};
  
class CollectLoadsMsg;

class Rebalancer {
public:
  struct pcpair {
    processorInfo *p;
    computeInfo *c;
    pcpair() : p(0),c(0) {;}
    void reset () { p = 0; c = 0; }
  };

  Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
             processorInfo *processorArray,
             int nComps, int nPatches, int nPes);
  ~Rebalancer();

protected: 
  int bytesPerAtom;
  typedef pcpair pcgrid[3][3][2];
  ProxyUsage  proxyUsage;
  const char *strategyName;
  computeInfo *computes;
  patchInfo *patches;
  processorInfo *processors;
  minHeap *pes;
  maxHeap *computePairHeap;
  maxHeap *computeSelfHeap;
  maxHeap *computeBgPairHeap;
  maxHeap *computeBgSelfHeap;
  int P;
  int numPatches;
  int numComputes;
  int numProxies;
  int numPesAvailable;
  double averageLoad;
  double origMaxLoad;
  int firstAssignInRefine;
#if USE_TOPOMAP
  TopoManager tmgr;
#endif

  int isAvailableOn(patchInfo *patch, processorInfo *p);
  void numAvailable(computeInfo *c, processorInfo *p,
	      int *nPatches, int *nProxies, int *isBadForCommunication);
  void refine_togrid(pcgrid &grid, double thresholdLoad,
	      processorInfo *p, computeInfo *c);
  void strategy();
  void makeHeaps();
  void makeTwoHeaps();
  void assign(computeInfo *c, processorInfo *pRec);
  void assign(computeInfo *c, int p);
  void deAssign(computeInfo *c, processorInfo *pRec);
  int refine();
  void multirefine(double overload_start=1.02);
  void printSummary();
  void printResults();
  void printLoads(int phase=0);
  CollectLoadsMsg *collMsg;
  double computeAverage();
  void adjustBackgroundLoadAndComputeAverage();
  double computeMax();
  double overLoad;
  void createSpanningTree();
  void decrSTLoad();
  void incrSTLoad();
  void InitProxyUsage();

  /** \function brickDim
   *  This function returns the coordinates of the inner brick
   *  between any two points on the torus. The coordinates need
   *  to be seen modulo the dimension in that direction
  */
  inline void brickDim(int a, int b, int dim, int &min, int &max) {
    int x1, x2, x3, x4, temp, i;
    if(a < b)
      { x1 = a; x2 = b; } 
    else
      { x1 = b; x2 = a; }

    x3 = x2 - x1;
    x4 = dim - x3;
    if(x3 < x4) {
      min = x1; max = x2;
    } else {
      min = x2; max = x1 + dim;
    }
  }

  /** \function withinBrick
   *  This function returns whether a particular coordinate is
   *  within the region defined by xm, xM, ym, yM, zm, zM
   */
  inline int withinBrick(int x, int y, int z, int xm, int xM, int dimX,
			 int ym, int yM, int dimY, int zm, int zM, int dimZ) {
    int wbX, wbY, wbZ;
    if( ((x >= xm) && (x <= xM)) || ((x < xm) && (x+dimX <= xM)) ) wbX = 1; else return 0;
    if( ((y >= ym) && (y <= yM)) || ((y < ym) && (y+dimY <= yM)) ) wbY = 1; else return 0;
    if( ((z >= zm) && (z <= zM)) || ((z < zm) && (z+dimZ <= zM)) ) wbZ = 1; else return 0;

    //if( wbX && wbY && wbZ)
      return 1;
  }

};

#endif
