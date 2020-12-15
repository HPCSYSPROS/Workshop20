/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h>
#include <stdio.h>

#include "InfoStream.h"
#include "ComputeMap.h"
#include "Compute.h"

#include "charm++.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

#include "ComputeNonbondedMICKernel.h"

ComputeMap* ComputeMap::instance;

// Singleton method
ComputeMap *ComputeMap::Instance() {
  if (instance == 0) {
    instance = new ComputeMap;	// this is never deleted
  }
  return instance;
}


//----------------------------------------------------------------------
ComputeMap::ComputeMap(void)
{
  nComputes=0;
  computePtrs=0;
}

//----------------------------------------------------------------------
ComputeMap::~ComputeMap(void)
{
  delete [] computePtrs;
}

void
ComputeMap::checkMap(void)
{
  int computeCount = nComputes;
  for (int i=0; i<nComputes; i++) {
    if (computePtrs[i]) {
      computeCount++;
      if (! (computePtrs[i]->cid == i)) {
	DebugM(4, "ComputeID("<<computePtrs[i]->cid<<") != ComputeID("
	  << i <<")\n");
      }
    }
  }
  DebugM(4, "Compute Count = " << computeCount << "\n");
}

void ComputeMap::pack (ComputeData *buffer)
{
  DebugM(4,"Packing ComputeMap\n");
  memcpy(buffer, computeData.begin(), nComputes * sizeof(ComputeData));
}

void ComputeMap::unpack (int n, ComputeData *ptr)
{
  DebugM(4,"Unpacking ComputeMap\n");

  if ( nComputes && n != nComputes ) {
    NAMD_bug("number of computes in new ComputeMap has changed!\n");
  }

  nComputes = n;
  computeData.resize(nComputes);
  memcpy(computeData.begin(), ptr, nComputes * sizeof(ComputeData));
}

void ComputeMap::initPtrs() {
  if ( ! computePtrs ) {
    computePtrs = new Compute*[nComputes];
    memset(computePtrs, 0, nComputes*sizeof(Compute*));
  }
}

void ComputeMap::extendPtrs() {
  if ( ! computePtrs ) NAMD_bug("ComputeMap::extendPtrs() 1");
  int oldN = nComputes;
  nComputes = computeData.size();
  if ( nComputes > oldN ) {
    Compute **oldPtrs = computePtrs;
    computePtrs = new Compute*[nComputes];
    memcpy(computePtrs, oldPtrs, oldN*sizeof(Compute*));
    memset(computePtrs+oldN, 0, (nComputes-oldN)*sizeof(Compute*));
    delete [] oldPtrs;
  }
}

//----------------------------------------------------------------------
int ComputeMap::numPids(ComputeID cid)
{
    return computeData[cid].numPids;
}

//----------------------------------------------------------------------
int ComputeMap::pid(ComputeID cid,int i)
{
    return computeData[cid].pids[i].pid;
}

int ComputeMap::trans(ComputeID cid,int i)
{
    return computeData[cid].pids[i].trans;
}

//----------------------------------------------------------------------
ComputeType ComputeMap::type(ComputeID cid)
{
  if (nComputes)
    return computeData[cid].type;
  else return computeErrorType;
}

//----------------------------------------------------------------------
int ComputeMap::partition(ComputeID cid)
{
  if (nComputes)
    return computeData[cid].partition;
  else return computeErrorType;
}
//----------------------------------------------------------------------
int ComputeMap::numPartitions(ComputeID cid)
{
  if (nComputes)
    return computeData[cid].numPartitions;
  else return computeErrorType;
}

//----------------------------------------------------------------------
int ComputeMap::allocateCids()
{
  nComputes = 0;
  computeData.resize(500);
  computeData.resize(0);

  return 0;
}

//----------------------------------------------------------------------
ComputeID ComputeMap::storeCompute(int inode, int maxPids, 
				   ComputeType type, 
				   int partition,int numPartitions)
{
  if (maxPids > numPidsAllocated) {
    NAMD_bug("ComputeMap::storeCompute called with maxPids > numPidsAllocated");
  }

  int cid;

  cid = nComputes;
  nComputes++;
  computeData.resize(nComputes);

  computeData[cid].node=inode;

  computeData[cid].type = type;
  computeData[cid].partition = partition;
  computeData[cid].numPartitions = numPartitions;

  computeData[cid].numPids = 0;

  #if defined(NAMD_MIC)
    // Initially in MIC runs, all computes are mapped to the host.  The host vs
    //   device LDB scheme will change this mapping later.
    computeData[cid].directToDevice = 0;
  #endif

  return cid;
}

//----------------------------------------------------------------------
ComputeID ComputeMap::cloneCompute(ComputeID src, int partition)
{
  const int cid = computeData.size();
  computeData.resize(cid+1);

  computeData[cid] = computeData[src];
  computeData[cid].partition = partition;
  computeData[cid].node = -1;

  return cid;
}

//----------------------------------------------------------------------
void ComputeMap::newPid(ComputeID cid, PatchID pid, int trans)
{
  computeData[cid].pids[computeData[cid].numPids].pid=pid;
  computeData[cid].pids[computeData[cid].numPids].trans=trans;
  computeData[cid].numPids++;
}

//----------------------------------------------------------------------
void ComputeMap::printComputeMap(void)
{
  DebugM(2,"---------------------------------------");
  DebugM(2,"---------------------------------------\n");

  DebugM(2,"nComputes = " << nComputes << '\n');
  DebugM(2,"nAllocated = " << nComputes << '\n');
  for(int i=0; i < nComputes; i++)
  {
    DebugM(2,"Compute " << i << '\n');
    DebugM(2,"  node = " << computeData[i].node << '\n');
    DebugM(2,"  numPids = " << computeData[i].numPids << '\n');
    for(int j=0; j < computeData[i].numPids; j++)
    {
      DebugM(2,computeData[i].pids[j].pid);
      if (!((j+1) % 6))
	DebugM(2,'\n');
    }
    DebugM(2,"\n---------------------------------------");
    DebugM(2,"---------------------------------------\n");

  }

char *fname;
#ifdef MEM_OPT_VERSION
fname = "computeMap.opt";
#else
fname = "computeMap.orig";
#endif

  FILE *ofp = fopen(fname, "w");
  fprintf(ofp,"---------------------------------------");
  fprintf(ofp,"---------------------------------------\n");

  fprintf(ofp,"nComputes = %d\n", nComputes);
  fprintf(ofp,"nAllocated = %d\n", nComputes);
  for(int i=0; i < nComputes; i++)
  {
    fprintf(ofp,"Compute %d\n", i);
    fprintf(ofp,"  node = %d\n",computeData[i].node);
    fprintf(ofp,"  numPids = %d\n",computeData[i].numPids);
	fprintf(ofp,"  type = %d\n",computeData[i].type);
    for(int j=0; j < computeData[i].numPids; j++)
    {
      fprintf(ofp,"%d ",computeData[i].pids[j].pid);
      if (!((j+1) % 6))
	fprintf(ofp,"\n");
    }
    fprintf(ofp,"\n---------------------------------------");
    fprintf(ofp,"---------------------------------------\n");

  }

fclose(ofp);

}

void ComputeMap::saveComputeMap(const char *fname)
{
  static int count = 0;
  char f[128];
  sprintf(f, "%s.%d", fname, count++);
  FILE *fp = fopen(f, "w");
  CmiAssert(fp != NULL);
  fprintf(fp, "%d\n", nComputes);
  for(int i=0; i < nComputes; i++)
  {
    fprintf(fp, "%d\n", computeData[i].node);
  }
  fclose(fp);
  CkPrintf("ComputeMap has been stored in %s.\n", f);
}

void ComputeMap::loadComputeMap(const char *fname)
{
  FILE *fp = fopen(fname, "r");
  CmiAssert(fp != NULL);
  int n;
  fscanf(fp, "%d\n", &n);
  CmiAssert(n == nComputes);
  for(int i=0; i < nComputes; i++)
  {
    fscanf(fp, "%d\n", &computeData[i].node);
  }
  fclose(fp);
}

//----------------------------------------------------------------------
#if defined(NAMD_MIC)

void ComputeMap::setDirectToDevice(const ComputeID cid, const int d) {
  if (cid < 0 || cid >= nComputes) {
    NAMD_bug("ComputeMap::setDirectToDevice() called with an invalid cid value");
  }
  computeData[cid].directToDevice = ((d == 0) ? (0) : (1));
}

int ComputeMap::directToDevice(const ComputeID cid) const {
  if (cid < 0 || cid >= nComputes) {
    NAMD_bug("ComputeMap::directToDevice() called with an invalid cid value");
  }
  return computeData[cid].directToDevice;
}

#endif // defined(NAMD_MIC)
