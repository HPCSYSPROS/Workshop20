/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "common.h"
#include "InfoStream.h"
#include "AlgRecBisection.h"

//#define DEBUG

//
//   load balancing computes using recursion of bisection algorithm
//
//
AlgRecBisection::AlgRecBisection(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
  Rebalancer(computeArray, patchArray, 
	     processorArray, nComps, 
	     nPatches, nPes)
{
strategyName = "Rob";
strategy();
}


void AlgRecBisection::rec_divide(int n, Partition &p)
{
  int i;
  int midpos;
  int n1, n2;
  double load1, currentload;
  int maxdir, count;
  Partition p1, p2;

#ifdef DEBUG
  CmiPrintf("rec_divide: partition n:%d count:%d load:%f (%d %d %d, %d %d %d)\n", n, p.count, p.load, p.origin[0], p.origin[1], p.origin[2], p.corner[0], p.corner[1], p.corner[2]);
#endif

  if (n==1) {
    partitions[currentp++] = p;
    return;
  }
/*
  if (p.origin.x==p.corner.x && p.origin.y==p.corner.y && p.origin.z==p.corner.z) 
     NAMD_die("AlgRecBisection failed in recursion.\n"); 
*/

  n2 = n/2;
  n1 = n-n2;

  load1 = (1.0*n1/n) * p.load;

  p1 = p;
  p1.refno = ++refno;
  p2 = p;
  p2.refno = ++refno;

  // determine the best division direction
  int maxSpan=-1;
  maxdir = XDIR;
  for (i=XDIR; i<=ZDIR; i++) {
    int myspan = p.corner[i] - p.origin[i];
    if (myspan > maxSpan) {
      maxdir = i;
      maxSpan = myspan;
    }
  }

  // other two dimensions
  int dir2 = (maxdir+1)%3;
  int dir3 = (maxdir+2)%3;

  currentload = 0.0;
  count = 0;
  midpos = p.origin[maxdir];
  for (i=0; i<numComputes; i++) {
    // not belong to this partition
    if (computeLoad[vArray[maxdir][i].id].refno != p.refno) continue;
    if (vArray[maxdir][i].v<p.origin[maxdir]) continue;
    if (vArray[maxdir][i].v>p.corner[maxdir]) break;

    int cid = vArray[maxdir][i].id;	// this compute ID
    // check if this compute is within the partition
    if ( computeLoad[cid].v[dir2] >= p.origin[dir2] &&
	 computeLoad[cid].v[dir2] <= p.corner[dir2] &&
	 computeLoad[cid].v[dir3] >= p.origin[dir3] &&
	 computeLoad[cid].v[dir3] <= p.corner[dir3]  ) {
      // this compute is set to the first partition
      if (currentload <= load1) {
	computeLoad[cid].refno = p1.refno;
        currentload += computeLoad[cid].load;
        count ++;
	midpos = computeLoad[cid].v[maxdir];
      }
      else {	// or the next partition
/*
	double c = currentload + computeLoad[cid].load;
	if (c - load1 < load1 - currentload) {
	  computeLoad[cid].refno = p1.refno;
	  currentload = c;
	  count ++;
	  midpos = computeLoad[cid].v[maxdir];
	}
	else
*/
	computeLoad[cid].refno = p2.refno;
      }
    }
  }
//  CmiPrintf("X:cur:%d, prev:%d load:%f %f\n", cur, prev, currentload, prevload);
#ifdef DEBUG
  CmiPrintf("DIR:%d %d load:%f\n", maxdir, midpos, currentload);
#endif


  p1.corner[maxdir] = midpos;
  p2.origin[maxdir] = midpos;

  p1.load = currentload;
  p1.count = count;
  p2.load = p.load - p1.load;
  p2.count = p.count - p1.count;
#ifdef DEBUG
  CmiPrintf("p1: n:%d count:%d load:%f\n", n1, p1.count, p1.load);
  CmiPrintf("p2: n:%d count:%d load:%f\n", n2, p2.count, p2.load);
#endif
  CmiAssert(! (p1 == p2));
  rec_divide(n1, p1);
  rec_divide(n2, p2);
}

void AlgRecBisection::setVal(int x, int y, int z)
{
  int i;
  for (i=0; i<numComputes; i++) {
    computeLoad[i].tv = 1000000*computeLoad[i].v[x]+
			1000*computeLoad[i].v[y]+
			computeLoad[i].v[z];
  }
#if 0
  CmiPrintf("original:%d\n", x);
  for (i=0; i<numComputes; i++) 
    CmiPrintf("%d ", computeLoad[i].tv);
  CmiPrintf("\n");
#endif
}

int AlgRecBisection::sort_partition(int x, int p, int r)
{
  int mid = computeLoad[vArray[x][p].id].tv;
  int i= p;
  int j= r;
  while (1) {
    while (computeLoad[vArray[x][j].id].tv > mid && j>i) j--;
    while (computeLoad[vArray[x][i].id].tv < mid && i<j) i++;
    if (i<j) {
      if (computeLoad[vArray[x][i].id].tv == computeLoad[vArray[x][j].id].tv)
      {
	if (computeLoad[vArray[x][i].id].tv != mid) NAMD_die("my god!\n");
	if (i-p < r-j) i++;
	else j--;
	continue;
      }
      VecArray tmp = vArray[x][i];
      vArray[x][i] = vArray[x][j];
      vArray[x][j] = tmp;
    }
    else
      return j;
  }
}

void AlgRecBisection::qsort(int x, int p, int r)
{
  if (p<r) {
    int q = sort_partition(x, p, r);
//CmiPrintf("midpoint: %d %d %d\n", p,q,r);
    qsort(x, p, q-1);
    qsort(x, q+1, r);
  }
}

void AlgRecBisection::quicksort(int x)
{
  int y = (x+1)%3;
  int z = (x+2)%3;
  setVal(x, y, z);
  qsort(x, 0, numComputes-1);

#if 0
  CmiPrintf("result:%d\n", x);
  for (int i=0; i<numComputes; i++) 
    CmiPrintf("%d ", computeLoad[vArray[x][i].id].tv);
  CmiPrintf("\n");
#endif
}

void AlgRecBisection::mapPartitionsToNodes()
{
  int i,j,k;
#if 0
  for (i=0; i<P; i++) partitions[i].node = i;
#else
  PatchMap *patchMap = PatchMap::Object();

  int **pool = new int *[P];
  for (i=0; i<P; i++) pool[i] = new int[P];
  for (i=0; i<P; i++) for (j=0; j<P; j++) pool[i][j] = 0;

  // sum up the number of nodes that patches of computes are on
  for (i=0; i<numComputes; i++)
  {
    for (j=0; j<P; j++)
      if (computeLoad[i].refno == partitions[j].refno) 
      {
	//int node1 = patchMap->node(computes[i].patch1);
	//int node2 = patchMap->node(computes[i].patch2);
        int node1 = patches[computes[i].patch1].processor;
        int node2 = patches[computes[i].patch2].processor;
	pool[j][node1]++;
	pool[j][node2]++;
      }
  }
#ifdef DEBUG
  for (i=0; i<P; i++) {
    for (j=0; j<P; j++) CmiPrintf("%d ", pool[i][j]);
    CmiPrintf("\n");
  }
#endif
  while (1)
  {
    int index=-1, node=0, eager=-1;
    for (j=0; j<P; j++) {
      if (partitions[j].node != -1) continue;
      int wantmost=-1, maxnodes=-1;
      for (k=0; k<P; k++) if (pool[j][k] > maxnodes && !partitions[k].mapped) {wantmost=k; maxnodes = pool[j][k];}
      if (maxnodes > eager) {
	index = j; node = wantmost; eager = maxnodes;
      }
    }
    if (eager == -1) break;
    partitions[index].node = node;
    partitions[node].mapped = 1;
  }

  for (i=0; i<P; i++) delete [] pool[i];
  delete [] pool;
#endif

  CmiPrintf("partitions to nodes mapping: ");
  for (i=0; i<P; i++) CmiPrintf("%d ", partitions[i].node);
  CmiPrintf("\n");
}

void AlgRecBisection::strategy()
{
  int i,j;

  PatchMap *patchMap = PatchMap::Object();

  // create computeLoad and calculate tentative computes coordinates
  computeLoad = new ComputeLoad[numComputes];
  for (i=XDIR; i<=ZDIR; i++) vArray[i] = new VecArray[numComputes];

  CmiPrintf("AlgRecBisection: numComputes:%d\n", numComputes);

  int size_x = patchMap->gridsize_a();
  int size_y = patchMap->gridsize_b();
  int size_z = patchMap->gridsize_c();

  // v[0] = XDIR  v[1] = YDIR v[2] = ZDIR
  // vArray[XDIR] is an array holding the x vector for all computes
  for (i=0; i<numComputes; i++) {
    int pid1 = computes[i].patch1;
    int pid2 = computes[i].patch2;
    computeLoad[i].id = i;
    int a1 = patchMap->index_a(pid1);
    int b1 = patchMap->index_b(pid1);
    int c1 = patchMap->index_c(pid1);
    int a2 = patchMap->index_a(pid2);
    int b2 = patchMap->index_b(pid2);
    int c2 = patchMap->index_c(pid1);
    computeLoad[i].v[0] = a1 + a2;
    computeLoad[i].v[1] = b1 + b2;
    computeLoad[i].v[2] = c1 + c2;
//    CmiPrintf("(%d %d %d)", computeLoad[i].x, computeLoad[i].y, computeLoad[i].z);
    computeLoad[i].load = computes[i].load;
    computeLoad[i].refno = 0;

    for (j=XDIR; j<=ZDIR; j++) {
      vArray[j][i].id = i;
      vArray[j][i].v = computeLoad[i].v[j];
    }
  }
//  CmiPrintf("\n");

  double t = CmiWallTimer();

  quicksort(XDIR);
  quicksort(YDIR);
  quicksort(ZDIR);
  CmiPrintf("qsort time: %f\n", CmiWallTimer() - t);

  npartition = P;
  partitions = new Partition[npartition];

  top_partition.origin[XDIR] = 0;
  top_partition.origin[YDIR] = 0;
  top_partition.origin[ZDIR] = 0;
  top_partition.corner[XDIR] = 2*(size_x-1);
  top_partition.corner[YDIR] = 2*(size_y-1);
  top_partition.corner[ZDIR] = 2*(size_z-1);

  top_partition.refno = 0;
  top_partition.load = 0.0;
  top_partition.count = numComputes;
  for (i=0; i<numComputes; i++) top_partition.load += computes[i].load;

  currentp = 0;
  refno = 0;

  // recursively divide
  rec_divide(npartition, top_partition);

  CmiPrintf("After partitioning: \n");
  for (i=0; i<P; i++) {
    CmiPrintf("[%d] (%d,%d,%d) (%d,%d,%d) load:%f count:%d\n", i, partitions[i].origin[0], partitions[i].origin[1], partitions[i].origin[2], partitions[i].corner[0], partitions[i].corner[1], partitions[i].corner[2], partitions[i].load, partitions[i].count);
  }

  for (i=0; i<numComputes; i++) computes[i].processor=-1;

  // mapping partitions to nodes
  mapPartitionsToNodes();

  // this is for debugging
  int *num = new int[P];
  for (i=0; i<P; i++) num[i] = 0;

  for (i=0; i<numComputes; i++)
  {
    for (j=0; j<P; j++)
      if (computeLoad[i].refno == partitions[j].refno) 
        { computes[computeLoad[i].id].processor = partitions[j].node; num[j]++; break; }
  }

  for (i=0; i<P; i++)
    if (num[i] != partitions[i].count) 
      NAMD_die("AlgRecBisection: Compute counts don't agree!\n");

  delete [] num;

  for (i=0; i<numComputes; i++) {
    if (computes[i].processor == -1) NAMD_bug("AlgRecBisection failure!\n");
  }


  delete [] computeLoad;
  for (i=0; i<3; i++) delete [] vArray[i];
  delete [] partitions;

  // use refinement
  for (i=0; i<numComputes; i++)
    assign((computeInfo *) &(computes[i]),
	   (processorInfo *) &(processors[computes[i].processor]));
	 
  printLoads();

#if 1
  multirefine();
#else
  printSummary();
#endif

  printLoads();

  CmiPrintf("AlgRecBisection finished time: %f\n", CmiWallTimer() - t);

}


