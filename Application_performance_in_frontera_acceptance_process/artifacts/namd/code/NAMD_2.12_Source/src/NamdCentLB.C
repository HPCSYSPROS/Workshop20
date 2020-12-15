/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdCentLB.C,v $
 * $Author: jim $
 * $Date: 2015/01/16 21:36:18 $
 * $Revision: 1.124 $
 *****************************************************************************/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include "InfoStream.h"
#include "NamdCentLB.h"
#include "NamdCentLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

// #define DUMP_LDBDATA 1
// #define LOAD_LDBDATA 1

double *cpuloads = NULL;

void CreateNamdCentLB() {
  // CkPrintf("[%d] creating NamdCentLB %d\n",CkMyPe(),loadbalancer);
  loadbalancer = CProxy_NamdCentLB::ckNew();
  // CkPrintf("[%d] created NamdCentLB %d\n",CkMyPe(),loadbalancer);
  if (CkMyRank() == 0 && cpuloads == NULL) {    
    cpuloads = new double[CkNumPes()];
    CmiAssert(cpuloads != NULL);
    for (int i=0; i<CkNumPes(); i++) cpuloads[i] = 0.0;
  }
}

NamdCentLB *AllocateNamdCentLB() {
  return new NamdCentLB((CkMigrateMessage*)NULL);
}

/**
 * Migratable Object Constructor.
 */
NamdCentLB::NamdCentLB(CkMigrateMessage *msg): CentralLB(msg) {
  processorArray = 0;
  patchArray = 0;
  computeArray = 0;
} 

NamdCentLB::NamdCentLB(): CentralLB(CkLBOptions(-1))
{
  //  if (CkMyPe()==0)
  //   CkPrintf("[%d] NamdCentLB created\n",CkMyPe());
  processorArray = 0;
  patchArray = 0;
  computeArray = 0;
}

/*
NamdCentLB::~NamdCentLB()
{
  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;
}
*/

bool NamdCentLB::QueryBalanceNow(int _step)
{
  //  CkPrintf("[%d] Balancing on step %d\n",CkMyPe(),_step);
  if ( LdbCoordinator::Object()->takingLdbData ) {
    return true;
  } else {
    return false;
  }
}

bool NamdCentLB::QueryDumpData()
{
#if 0
  if (LdbCoordinator::Object()->ldbCycleNum == 1)  return true;
  if (LdbCoordinator::Object()->ldbCycleNum == 2)  return true;
#endif
  return false;
}

CLBMigrateMsg* NamdCentLB::Strategy(LDStats* stats)
{
  //  CkPrintf("LDB: All statistics received at %f, %f\n",
  //  CmiTimer(),CmiWallTimer());

  int numProcessors = stats->nprocs();
  int numPatches = PatchMap::Object()->numPatches();
  ComputeMap *computeMap = ComputeMap::Object();
  const int numComputes = computeMap->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  // these sizes should never change
  if ( ! processorArray ) processorArray = new processorInfo[numProcessors];
  if ( ! patchArray ) patchArray = new patchInfo[numPatches];
  if ( ! computeArray ) computeArray = new computeInfo[numComputes];

  int nMoveableComputes = buildData(stats);

#if LDB_DEBUG
#define DUMP_LDBDATA 1
#define LOAD_LDBDATA 1
#endif

#if DUMP_LDBDATA 
  dumpDataASCII("ldbd_before", numProcessors, numPatches, nMoveableComputes);
#elif LOAD_LDBDATA
  loadDataASCII("ldbd_before.5", numProcessors, numPatches, nMoveableComputes);
  // CkExit();
#endif

  double averageLoad = 0.;
  double avgCompute;
  {
   int i;
   double total = 0.;
   double maxCompute = 0.;
   int maxi = 0;
   for (i=0; i<nMoveableComputes; i++) {
      double load = computeArray[i].load;
      total += load;
      if ( load > maxCompute ) { maxCompute = load;  maxi = i; }
   }
   avgCompute = total / nMoveableComputes;

    int P = stats->nprocs();
   int numPesAvailable = 0;
   for (i=0; i<P; i++) {
      if (processorArray[i].available) {
        ++numPesAvailable;
        total += processorArray[i].backgroundLoad;
      }
   }
   if (numPesAvailable == 0)
     NAMD_die("No processors available for load balancing!\n");

   averageLoad = total/numPesAvailable;
   CkPrintf("LDB: Largest compute %d load %f is %.1f%% of average load %f\n",
            computeArray[maxi].handle.id.id[0],
            maxCompute, 100. * maxCompute / averageLoad, averageLoad);
   CkPrintf("LDB: Average compute %f is %.1f%% of average load %f\n",
            avgCompute, 100. * avgCompute / averageLoad, averageLoad);
  }

  if ( step() == 1 ) {
    // compute splitting only
    // partitions are stored as char but mostly limited by
    // high load noise at low outer-loop iteration counts
    int maxParts = 10;
#ifdef NAMD_CUDA
//split LCPO compute very small, else CUDA compute is delayed
    if (simParams->LCPOOn) {
      maxParts = 20;
    }
#endif
    int totalAddedParts = 0;
    double maxCompute = averageLoad / 10.;
    if ( maxCompute < 2. * avgCompute ) maxCompute = 2. * avgCompute;
    if ( simParams->ldbRelativeGrainsize > 0. ) {
      maxCompute = averageLoad * simParams->ldbRelativeGrainsize;
    }
    CkPrintf("LDB: Partitioning computes with target load %f\n", maxCompute);
    double maxUnsplit = 0.;
    for (int i=0; i<nMoveableComputes; i++) {
      computeArray[i].processor = computeArray[i].oldProcessor;
      const int cid = computeArray[i].handle.id.id[0];
      const double load = computeArray[i].load;
      if ( computeMap->numPartitions(cid) == 0 ) {
        if ( load > maxUnsplit ) maxUnsplit = load;
        continue;
      }
      int nparts = (int) ceil(load / maxCompute);
      if ( nparts > maxParts ) nparts = maxParts;
      if ( nparts < 1 ) nparts = 1;
      if ( 0 && nparts > 1 ) {
        CkPrintf("LDB: Partitioning compute %d with load %f by %d\n",
                  cid, load, nparts);
      }
      computeMap->setNewNumPartitions(cid,nparts);
      totalAddedParts += nparts - 1;
    }
    CkPrintf("LDB: Increased migratable compute count from %d to %d\n",
              nMoveableComputes,nMoveableComputes+totalAddedParts);
    CkPrintf("LDB: Largest unpartitionable compute is %f\n", maxUnsplit);
  } else if (simParams->ldbStrategy == LDBSTRAT_DEFAULT) { // default
    if (step() < 4)
      TorusLB(computeArray, patchArray, processorArray,
	          nMoveableComputes, numPatches, numProcessors);
    else
      RefineTorusLB(computeArray, patchArray, processorArray,
                  nMoveableComputes, numPatches, numProcessors, 1);
  } else if (simParams->ldbStrategy == LDBSTRAT_COMPREHENSIVE) {
    TorusLB(computeArray, patchArray, processorArray,
	          nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_REFINEONLY) {
    RefineTorusLB(computeArray, patchArray, processorArray,
                  nMoveableComputes, numPatches, numProcessors, 1);
  } else if (simParams->ldbStrategy == LDBSTRAT_OLD) {
    if (step() < 4)
      Alg7(computeArray, patchArray, processorArray,
	          nMoveableComputes, numPatches, numProcessors);
    else
      RefineOnly(computeArray, patchArray, processorArray, 
                  nMoveableComputes, numPatches, numProcessors);
  }

#if LDB_DEBUG && USE_TOPOMAP
  TopoManager tmgr;
  int pe1, pe2, pe3, hops=0;
  /* This is double counting the hops
  for(int i=0; i<nMoveableComputes; i++)
  {
    pe1 = computeArray[i].processor;
    pe2 = patchArray[computeArray[i].patch1].processor;
    pe3 = patchArray[computeArray[i].patch2].processor;
    hops += tmgr.getHopsBetweenRanks(pe1, pe2);
    if(computeArray[i].patch1 != computeArray[i].patch2)
      hops += tmgr.getHopsBetweenRanks(pe1, pe3);  
  }*/
  for (int i=0; i<numPatches; i++)  {
    //int num = patchArray[i].proxiesOn.numElements();
    pe1 = patchArray[i].processor;
    Iterator nextProc;
    processorInfo *p = (processorInfo *)patchArray[i].proxiesOn.iterator((Iterator *)&nextProc);
    while (p) {
      pe2 = p->Id;
      hops += tmgr.getHopsBetweenRanks(pe1, pe2);
      p = (processorInfo *)patchArray[i].proxiesOn.next((Iterator*)&nextProc);
    }
  }
  CkPrintf("Load Balancing: Number of Hops: %d\n", hops);
#endif

#if DUMP_LDBDATA
  dumpDataASCII("ldbd_after", numProcessors, numPatches, nMoveableComputes);
#elif LOAD_LDBDATA
  dumpDataASCII("ldbd_after.5", numProcessors, numPatches, nMoveableComputes);
  // loadDataASCII("ldbd_after", numProcessors, numPatches, nMoveableComputes);
  // CkExit();
#endif

  // For error checking:
  // Count up computes, to see if somebody doesn't have any computes
  int i;
#if 0
  int* computeCount = new int[numProcessors];
  for(i=0; i<numProcessors; i++)
    computeCount[i]=0;
  for(i=0; i<nMoveableComputes; i++)
    computeCount[computeArray[i].processor]++;
  for(i=0; i<numProcessors; i++) {
    if (computeCount[i]==0)
      iout << iINFO <<"Warning: Processor " << i 
	   << " has NO moveable computes.\n" << endi;
  }
  delete [] computeCount;
#endif
  
  CkVec<MigrateInfo *> migrateInfo;
  for(i=0;i<nMoveableComputes;i++) {
    if (computeArray[i].processor != computeArray[i].oldProcessor) {
      //      CkPrintf("[%d] Obj %d migrating from %d to %d\n",
      //               CkMyPe(),computeArray[i].handle.id.id[0],
      //	       computeArray[i].processor,computeArray[i].oldProcessor);
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      migrateMe->from_pe = computeArray[i].oldProcessor;
      migrateMe->to_pe = computeArray[i].processor;
      migrateInfo.insertAtEnd(migrateMe);

      // sneak in updates to ComputeMap
      computeMap->setNewNode(computeArray[i].handle.id.id[0],
	 			computeArray[i].processor);
    }
  }
  
  int migrate_count=migrateInfo.length();
  // CkPrintf("NamdCentLB migrating %d elements\n",migrate_count);
  CLBMigrateMsg* msg = new(migrate_count,CkNumPes(),CkNumPes(),0) CLBMigrateMsg;

  msg->n_moves = migrate_count;
  for(i=0; i < migrate_count; i++) {
    MigrateInfo* item = migrateInfo[i];
    msg->moves[i] = *item;
    delete item;
    migrateInfo[i] = 0;
  }

  for (i=0; i<numProcessors; i++) {
    cpuloads[i] = processorArray[i].load;
  }

  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;

  processorArray = NULL;
  patchArray = NULL;
  computeArray = NULL;
  
  return msg;
};

#ifndef WIN32

void NamdCentLB::dumpDataASCII(char *file, int numProcessors,
			       int numPatches, int numComputes)
{
  char filename[128];
  sprintf(filename, "%s.%d", file, step());
  FILE* fp = fopen(filename,"w");
  if (fp == NULL){
     perror("dumpLDStatsASCII");
     return;
  }
  CkPrintf("***** DUMP data to file: %s ***** \n", filename);
  fprintf(fp,"%d %d %d\n",numProcessors,numPatches,numComputes);

  int i;
  for(i=0;i<numProcessors;i++) {
    processorInfo* p = processorArray + i;
    fprintf(fp,"%d %e %e %e %e\n",p->Id,p->load,p->backgroundLoad,p->computeLoad,p->idleTime);
  }

  for(i=0;i < numPatches; i++) {
    patchInfo* p = patchArray + i;
    fprintf(fp,"%d %e %d %d\n",p->Id,p->load,p->processor,p->numAtoms);
  }
    
  for(i=0; i < numComputes; i++) {
    computeInfo* c = computeArray + i;
    fprintf(fp,"%d %e %d %d %d %d",c->Id,c->load,c->patch1,c->patch2,
	    c->processor,c->oldProcessor);
    fprintf(fp, "\n");
  }

  // dump patchSet
  for (i=0; i< numProcessors; i++) {
      int num = processorArray[i].proxies.numElements();
      fprintf(fp, "%d %d: ", i, num);
      Iterator nextProxy;
      patchInfo *p = (patchInfo *)processorArray[i].proxies.
	iterator((Iterator *)&nextProxy);
      while (p) {
          fprintf(fp, "%d ", p->Id);
          p = (patchInfo *)processorArray[i].proxies.
	    next((Iterator*)&nextProxy);
      }
      fprintf(fp, "\n");
  }
  // dump proxiesOn
  for (i=0; i<numPatches; i++)  {
    int num = patchArray[i].proxiesOn.numElements();
    fprintf(fp, "%d %d: ", i, num);
      Iterator nextProc;
      processorInfo *p = (processorInfo *)patchArray[i].proxiesOn.
	iterator((Iterator *)&nextProc);
      while (p) {
	fprintf(fp, "%d ", p->Id);
	p = (processorInfo *)patchArray[i].proxiesOn.
	  next((Iterator*)&nextProc);
      }
      fprintf(fp, "\n");
  }

  fclose(fp);
  //CkExit();
}

void NamdCentLB::loadDataASCII(char *file, int &numProcessors,
			       int &numPatches, int &numComputes)
{
  char filename[128];
  //sprintf(filename, "%s.%d", file, step());
  sprintf(filename, "%s", file);

  CkPrintf("***** Load ascii data from file: %s ***** \n", filename);

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
     perror("loadDataASCII");
     return;
  }

  fscanf(fp,"%d %d %d",&numProcessors,&numPatches,&numComputes);

  printf("numProcs: %d numPatches: %d numComputes: %d\n", numProcessors,numPatches, numComputes);

  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;
  processorArray = new processorInfo[numProcessors];
  patchArray = new patchInfo[numPatches];
  computeArray = new computeInfo[numComputes];

  int i;
  for(i=0;i<numProcessors;i++) {
    processorInfo* p = processorArray + i;
    fscanf(fp,"%d %le %le %le", &p->Id, &p->load, &p->backgroundLoad, &p->computeLoad);
    fscanf(fp,"%le\n", &p->idleTime);
    if (p->Id != i) CmiAbort("Reading processorArray error!");
//    p->backgroundLoad = 0.0;
  }

  for(i=0;i < numPatches; i++) {
    patchInfo* p = patchArray + i;
    fscanf(fp,"%d %le %d %d\n",&p->Id,&p->load,&p->processor,&p->numAtoms);
    if (p->Id != i || p->processor > numProcessors || p->processor < 0) 
      CmiAbort("Reading patchArray error!");
  }
    
  for(i=0; i < numComputes; i++) {
    computeInfo* c = computeArray + i;
    fscanf(fp,"%d %le %d %d %d %d",&c->Id,&c->load,&c->patch1,&c->patch2,
	    &c->processor,&c->oldProcessor);

    if (c->patch1 < 0 || c->patch1 > numPatches || c->patch2 < 0 || c->patch2 > numPatches)
      CmiAbort("Reading computeArray error!");
  // printf("%d %e %d %d %d %d\n", c->Id,c->load,c->patch1,c->patch2,c->processor,c->oldProcessor);
  }

  // dump patchSet
  for (i=0; i< numProcessors; i++) {
      int num, curp;
      fscanf(fp,"%d %d: ",&curp, &num);
      if(curp != i)
	CmiAbort("Reading patchsSet error!");
      for (int j=0; j<num; j++) {
          int id;
          fscanf(fp,"%d",&id);
          processorArray[i].proxies.unchecked_insert(&patchArray[id]);
      }
  }
  // dump proxiesOn
  for (i=0; i<numPatches; i++)  {
      int num, curp;
      fscanf(fp,"%d %d: ",&curp, &num);
      if(curp != i)
	CmiAbort("Reading proxiesOn error!");
      for (int j=0; j<num; j++) {
          int id;
	  fscanf(fp,"%d",&id);
          patchArray[i].proxiesOn.insert(&processorArray[id]);
      }
  }

  fclose(fp);
}
#endif

extern int isPmeProcessor(int); 
#ifdef MEM_OPT_VERSION
extern int isOutputProcessor(int); 
#endif
#if defined(NAMD_MIC)
extern int isMICProcessor(int);
#endif

int NamdCentLB::buildData(LDStats* stats)
{
  int n_pes = stats->nprocs();

  PatchMap* patchMap = PatchMap::Object();
  ComputeMap* computeMap = ComputeMap::Object();
  const SimParameters* simParams = Node::Object()->simParameters;

  BigReal bgfactor = simParams->ldbBackgroundScaling;
  BigReal pmebgfactor = simParams->ldbPMEBackgroundScaling;
  BigReal homebgfactor = simParams->ldbHomeBackgroundScaling;
  int pmeOn = simParams->PMEOn;
  int unLoadPme = simParams->ldbUnloadPME;
  int pmeBarrier = simParams->PMEBarrier;
  int unLoadZero = simParams->ldbUnloadZero;
  int unLoadOne = simParams->ldbUnloadOne;
  int unLoadIO= simParams->ldbUnloadOutputPEs;
  int i;
  for (i=0; i<n_pes; ++i) {
    processorArray[i].Id = i;
    processorArray[i].available = true;
    if ( pmeOn && isPmeProcessor(i) ) {
      processorArray[i].backgroundLoad = pmebgfactor * stats->procs[i].bg_walltime;
    } else if (patchMap->numPatchesOnNode(i) > 0) {
      processorArray[i].backgroundLoad = homebgfactor * stats->procs[i].bg_walltime;
    } else {
      processorArray[i].backgroundLoad = bgfactor * stats->procs[i].bg_walltime;
    }
    processorArray[i].idleTime = stats->procs[i].idletime;
    processorArray[i].load = processorArray[i].computeLoad = 0.0;
  }

/* *********** this code is defunct *****************
#if 0
  double bgfactor = 1.0 + 1.0 * CkNumPes()/1000.0;
  if ( bgfactor > 2.0 ) bgfactor = 2.0;
  iout << iINFO << "Scaling background load by " << bgfactor << ".\n" << endi;
  int i;
  for (i=0; i<n_pes; i++) {
    processorArray[i].Id = i;
    processorArray[i].backgroundLoad = bgfactor * stats[i].bg_walltime;
  }

  double bg_weight = 0.7;

  int i;
  for (i=0; i<n_pes; i++) {
    processorArray[i].Id = i;
    if (patchMap->numPatchesOnNode(i) > 0)
      processorArray[i].backgroundLoad = bg_weight * stats->procs[i].bg_walltime;
    else 
      processorArray[i].backgroundLoad = stats[i].bg_walltime;
  }
  
  //Modification to reduce the coputeload on PME processors
  const SimParameters* simParams = Node::Object()->simParameters;  
  
  // CkPrintf("BACKGROUND LOAD\n");
  if(simParams->PMEOn) {
    double bgfactor = 1.0 + 1.0 * CkNumPes()/1000.0;
    if ( bgfactor > 2.0 ) bgfactor = 2.0;
    for (i=0; i<n_pes; i++) {
      // CkPrintf("BG[%d] =  %5.5lf,", i, processorArray[i].backgroundLoad);
      if(isPmeProcessor(i)) {
	processorArray[i].backgroundLoad *= bgfactor;
      }
      // CkPrintf("%5.5lf;  ", processorArray[i].backgroundLoad);
    }
  }
  // CkPrintf("\n");
#endif  
*********** end of defunct code *********** */

  if (unLoadZero) processorArray[0].available = false;
  if (unLoadOne) processorArray[1].available = false;

  // if all pes are Pme, disable this flag
  if (pmeOn && unLoadPme) {
    for (i=0; i<n_pes; i++) {
      if (!isPmeProcessor(i))  break;
    }
    if (i == n_pes) {
      iout << iINFO << "Turned off unLoadPme flag!\n"  << endi;
      unLoadPme = 0;
    }
  }
  
  if (pmeOn && unLoadPme) {
    for (i=0; i<n_pes; i++) {
      if ((pmeBarrier && i==0) || isPmeProcessor(i)) 
	processorArray[i].available = false;
    }
  }
  // if all pes are output, disable this flag
#ifdef MEM_OPT_VERSION

  if (unLoadIO) {
      if (simParams->numoutputprocs == n_pes) {
	  iout << iINFO << "Turned off unLoadIO flag!\n"  << endi;
	  unLoadIO = 0;
      }
  }
  if (unLoadIO){
    iout << iINFO << "Testing for output processors!\n"  << endi;
      for (i=0; i<n_pes; i++) {
	  if (isOutputProcessor(stats->procs[i].pe)) 
	    {
	      //	      iout << iINFO << "Removed output PE "<< stats->procs[i].pe <<" from available list!\n"  << endi;
	      processorArray[i].available = false;
	    }
	  else
	    {
	      //	      iout << iINFO << "Nonoutput PE "<< stats->procs[i].pe <<" is in available list!\n"  << endi;
	    }
      }
  }
#endif

  // Unload PEs driving MIC devices, if need be
  #if defined(NAMD_MIC)
    if (simParams->mic_unloadMICPEs != 0) {
      for (i = 0; i < n_pes; i++) {
        if (isMICProcessor(i) != 0) { processorArray[i].available = false; }
      }
    }
  #endif

  int nMoveableComputes=0;
  int nProxies = 0;		// total number of estimated proxies
  int nIdleComputes = 0;

  int j;
  for (j=0; j < stats->n_objs; j++) {
      const LDObjData &this_obj = stats->objData[j];
      int frompe = stats->from_proc[j];

      // filter out non-NAMD managed objects (like PME array)
      if (this_obj.omID().id.idx != 1) {
        // CkPrintf("non-NAMD object %d on pe %d with walltime %lf\n",
        // this_obj.id().id[0], stats->from_proc[j], this_obj.wallTime);
        processorArray[stats->from_proc[j]].backgroundLoad += this_obj.wallTime;
        continue;
      }

      if (this_obj.id().id[1] == -2) { // Its a patch
	const int pid = this_obj.id().id[0];
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

	patchArray[pid].Id = pid;
	patchArray[pid].numAtoms = 0;
	patchArray[pid].processor = stats->from_proc[j];
	const int numProxies = 
#if USE_TOPOMAP
	requiredProxiesOnProcGrid(pid,neighborNodes);
#else
	requiredProxies(pid, neighborNodes);
#endif

        nProxies += numProxies;

	for (int k=0; k<numProxies; k++) {
	  processorArray[neighborNodes[k]].proxies.unchecked_insert(&patchArray[pid]);
	  patchArray[pid].proxiesOn.unchecked_insert(&processorArray[neighborNodes[k]]);
	}
	processorArray[stats->from_proc[j]].backgroundLoad += this_obj.wallTime;
      } else if (this_obj.id().id[1] == -3) { // Its a bonded compute
	processorArray[stats->from_proc[j]].backgroundLoad += this_obj.wallTime;
      } else if (this_obj.migratable) { // Its a compute
       if ( this_obj.wallTime == 0. ) { // don't migrate idle computes
         ++nIdleComputes;
       } else {
	const int cid = this_obj.id().id[0];
	const int p0 = computeMap->pid(cid,0);

	// For self-interactions, just return the same pid twice
	int p1;
	if (computeMap->numPids(cid) > 1)
	  p1 = computeMap->pid(cid,1);
	else p1 = p0;
	computeArray[nMoveableComputes].Id = cid;
	computeArray[nMoveableComputes].oldProcessor = stats->from_proc[j];
	processorArray[stats->from_proc[j]].computeLoad += this_obj.wallTime;
	computeArray[nMoveableComputes].processor = -1;
	computeArray[nMoveableComputes].patch1 = p0;
	computeArray[nMoveableComputes].patch2 = p1;
	computeArray[nMoveableComputes].handle = this_obj.handle;
	computeArray[nMoveableComputes].load = this_obj.wallTime;
	nMoveableComputes++;
       }
      } else {
	processorArray[stats->from_proc[j]].backgroundLoad += this_obj.wallTime;
      }
    }

   if ( nIdleComputes )
     CkPrintf("LDB: %d computes have load of zero\n", nIdleComputes);

/* *********** this code is defunct *****************
#if 0
  int averageProxy = nProxies / n_pes;
  CkPrintf("total proxies: %d, avervage: %d\n", nProxies, averageProxy);
  for (i=0; i<n_pes; i++) {
    // too many proxies on this node, weight the background load
    int proxies = processorArray[i].proxies.numElements();
    if (proxies > averageProxy) {
      double factor = 1.0*(proxies-averageProxy)/nProxies;
      processorArray[i].backgroundLoad *= (1.0 + factor);
      CkPrintf("On [%d]: too many proxies: %d, increased bg load by %f\n", i, nProxies, factor);
    }
  }
#endif
*********** end of defunct code *********** */

  for (i=0; i<n_pes; i++) {
    processorArray[i].load = processorArray[i].backgroundLoad + processorArray[i].computeLoad;
  }
  stats->clear();
  return nMoveableComputes;
}

// Figure out which proxies we will definitely create on other
// nodes, without regard for non-bonded computes.  This code is swiped
// from ProxyMgr, and changes there probable need to be propagated here.

int NamdCentLB::requiredProxies(PatchID id, int neighborNodes[])
{
  PatchMap* patchMap = PatchMap::Object();
  int myNode = patchMap->node(id);
  int nProxyNodes = 0;

#define IF_NEW_NODE \
    int j; \
    for ( j=0; j<nProxyNodes && neighborNodes[j] != proxyNode; ++j ); \
    if ( j == nProxyNodes )

  PatchID neighbors[1 + PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
  neighbors[0] = id;
  int numNeighbors = 1 + patchMap->downstreamNeighbors(id,neighbors+1);
  for ( int i = 0; i < numNeighbors; ++i ) {
    const int proxyNode = patchMap->basenode(neighbors[i]);
    if ( proxyNode != myNode ) {
      IF_NEW_NODE {
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
    }
  }

  // Distribute initial default proxies across empty processors.
  // This shouldn't be necessary, but may constrain the load balancer
  // and avoid placing too many proxies on a single processor.  -JCP

  // This code needs to be turned off when the creation of ST is
  // shifted to the load balancers -ASB

#if 1
  int numPes = CkNumPes();
  int numPatches = patchMap->numPatches();
  int emptyNodes = numPes - numPatches;
  if ( emptyNodes > numPatches ) {
    int nodesPerPatch = nProxyNodes + 1 + (emptyNodes-1) / numPatches;
    int maxNodesPerPatch = PatchMap::MaxOneAway + PatchMap::MaxTwoAway;
    if ( nodesPerPatch > maxNodesPerPatch ) nodesPerPatch = maxNodesPerPatch;
    int proxyNode = (myNode + 1) % numPes;
    while ( nProxyNodes < nodesPerPatch &&
                        ! patchMap->numPatchesOnNode(proxyNode) ) {
      if ( proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
      proxyNode = (proxyNode + 1) % numPes;
    }
    proxyNode = (myNode - 1 + numPes) % numPes;
    while ( nProxyNodes < nodesPerPatch &&
                        ! patchMap->numPatchesOnNode(proxyNode) ) {
      if ( proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
      proxyNode = (proxyNode - 1 + numPes) % numPes;
    }
    proxyNode = (myNode + 1) % numPes;
    int count = 0;
    while ( nProxyNodes < nodesPerPatch ) {
      if ( ! patchMap->numPatchesOnNode(proxyNode) && proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
      proxyNode = (proxyNode + 1) % numPes;
      count ++; if (count == numPes) break;   // we looped all
    }
  } else {
    int proxyNode = myNode - 1;
    if ( proxyNode >= 0 && ! patchMap->numPatchesOnNode(proxyNode) ) {
      if ( proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
    }
    proxyNode = myNode + 1;
    if ( proxyNode < numPes && ! patchMap->numPatchesOnNode(proxyNode) ) {
      if ( proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
    }
  }
#endif

  return nProxyNodes;
}

#if USE_TOPOMAP 
// Figure out which proxies we will definitely create on other nodes,
// without regard for non-bonded computes.  This code is swiped from
// ProxyMgr, and changes there probable need to be propagated here.
// The proxies are placed on nearby processors on the 3d-grid along
// the X, Y, Z and T dimensions

int NamdCentLB::requiredProxiesOnProcGrid(PatchID id, int neighborNodes[])
{
  enum proxyHere { No, Yes };
  int numPes = CkNumPes();
  proxyHere *proxyNodes = new proxyHere[numPes];
  int nProxyNodes;
  int i, j, k, l;

  int xsize = 0, ysize = 0, zsize = 0, tsize = 0;
  int my_x = 0, my_y = 0, my_z = 0, my_t = 0;

  PatchMap* patchMap = PatchMap::Object();
  int myNode = patchMap->node(id);
    
  TopoManager tmgr;
  xsize = tmgr.getDimNX();
  ysize = tmgr.getDimNY();
  zsize = tmgr.getDimNZ();
  tsize = tmgr.getDimNT();
  
  tmgr.rankToCoordinates(myNode, my_x, my_y, my_z, my_t);
  
  if(xsize * ysize * zsize * tsize != CkNumPes()) {
    delete [] proxyNodes;
    return requiredProxies(id, neighborNodes);
  }  

  // Note all home patches.
  for ( i = 0; i < numPes; ++i )
  {
    proxyNodes[i] = No;
  }
  nProxyNodes = 0;

  // Check all two-away neighbors.
  // This is really just one-away neighbors, since 
  // two-away always returns zero: RKB
  PatchID neighbors[1 + PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

  // Assign a proxy to all your neighbors. But dont increment counter
  // because these have to be there anyway.
  neighbors[0] = id;  
  int numNeighbors = 1 + patchMap->downstreamNeighbors(id,neighbors+1);
  
  // Small Flag chooses between different loadbalancing schemes.
  // Small Flag == true, patches are close to each other
  // false, patches are far from each other
  bool smallFlag = false;
  double pnodes = CkNumPes();
  pnodes *= 0.25;    
  smallFlag = (patchMap->numPatches() > pnodes )?1:0;

  //If there are lot of patches its likely they will all be neighbors, 
  //so all we need to do is to place proxies on downstream patches.
  //if (smallFlag) {
  for ( i = 1; i < numNeighbors; ++i )
    {
      int proxyNode = patchMap->basenode(neighbors[i]);
      
      if (proxyNode != myNode)
	if (proxyNodes[proxyNode] == No)
	  {
	    proxyNodes[proxyNode] = Yes;
	    neighborNodes[nProxyNodes] = proxyNode;
	    nProxyNodes++;
	  }
    }
  //}
 
  if (step() > 2) {
    delete [] proxyNodes;
    return nProxyNodes;
  }
 
  // Place numPesPerPatch proxies on the 3d torus neighbors of a processor

  int numPatches = patchMap->numPatches();
  int emptyNodes = numPes - numPatches;
  //if ( emptyNodes > numPatches ) {
  
  int nodesPerPatch = nProxyNodes + 4 * (emptyNodes-1) / numPatches + 1;
  int proxyNode = 0 ;
  int proxy_x=0, proxy_y=0, proxy_z=0;
  
  //Choose from the 26 neighbors of mynode.
  //CkAssert(nodesPerPatch - nProxyNodes <= 26);  
  //Too few patches otherwise, try twoaway?
  
  for(k=-1; k<= 1; k++) {
    proxy_z = (my_z + k + zsize) % zsize;
    for(j=-1; j <= 1; j++) {
      proxy_y = (my_y + j + ysize) % ysize;
      for(i = -1; i <= 1; i++) {
	proxy_x = (my_x + i + xsize) % xsize;
	for(l = 0; l < tsize; l++) {
	  if(i == 0 && j == 0 && k == 0 && l == 0)
	    continue;

	  proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z, l);

	  if((! patchMap->numPatchesOnNode(proxyNode) || !smallFlag) &&
	     proxyNodes[proxyNode] == No) {
	    proxyNodes[proxyNode] = Yes;
	    neighborNodes[nProxyNodes] = proxyNode;
	    nProxyNodes++;
	  }
	  
	  if(nProxyNodes >= nodesPerPatch || 
	     nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	    break;
	} // end for

	if(nProxyNodes >= nodesPerPatch || 
	   nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	  break;
      } // end for
      
      if(nProxyNodes >= nodesPerPatch || 
	 nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	break;	  
    } // end for

    if(nProxyNodes >= nodesPerPatch || 
       nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
      break;	  
  } // end for

#if 1
  if(!smallFlag) {
    for(k=-2; k<= 2; k+=2) {
      proxy_z = (my_z + k + zsize) % zsize;
      for(j=-2; j <= 2; j+=2) {
	proxy_y = (my_y + j + ysize) % ysize;
	for(i = -2; i <= 2; i+=2) {
	  proxy_x = (my_x + i + xsize) % xsize;
	  for(l = 0; l < tsize; l++) {
	    if(i == 0 && j == 0 && k == 0 && l == 0)
	      continue;
	  
	    proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z, l);
	  
	    if((! patchMap->numPatchesOnNode(proxyNode) || !smallFlag) &&
	       proxyNodes[proxyNode] == No) {
	      proxyNodes[proxyNode] = Yes;
	      neighborNodes[nProxyNodes] = proxyNode;
	      nProxyNodes++;
	    }
	    
	    if(nProxyNodes >= nodesPerPatch || 
	       nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	      break;
	  } // end for

	  if(nProxyNodes >= nodesPerPatch || 
	     nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	    break;
	} // end for
	
	if(nProxyNodes >= nodesPerPatch || 
	   nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	  break;	  
      } // end for

      if(nProxyNodes >= nodesPerPatch || 
	 nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	break;	  
    } // end for
  }

#else
  #if 0
  const SimParameters* params = Node::Object()->simParameters;

  if(!smallFlag) {
    //Add two-away proxies
    if(patchMap->numaway_a() == 2) {
      proxy_y = (my_y + 2) % ysize;
      proxy_x = my_x  % xsize;
      proxy_z = my_z  % zsize;
      
      proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z);
      if(proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
      nProxyNodes++;
      }
      
      proxy_y = (my_y - 2 + ysize) % ysize;
      proxy_x = my_x  % xsize;
      proxy_z = my_z % zsize;
      
      proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z);
      if(proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
    }
    
    //Add two away proxies
    if(patchMap->numaway_b() == 2) {
      proxy_y = my_y  % ysize;
      proxy_x = my_x  % xsize;
      proxy_z = (my_z + 2) % zsize;
      
      proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z);
      if(proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
      
      proxy_y = my_y  % ysize;
      proxy_x = my_x  % xsize;
      proxy_z = (my_z - 2 + zsize) % zsize;
      
      proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z);
      if(proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
    }
    
    //Add two away proxies
    if(patchMap->numaway_c() == 2) {
      proxy_y = my_y  % ysize;
      proxy_x = (my_x + 2) % xsize;
      proxy_z = my_z  % zsize;
      
      proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z);
      if(proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
      nProxyNodes++;
      }
      
      proxy_y = my_y  % ysize;
      proxy_x = (my_x  - 2 + xsize) % xsize;
      proxy_z = my_z % zsize;
      
      proxyNode = tmgr.coordinatesToRank(proxy_x, proxy_y, proxy_z);
      if(proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
    }
  }
  #endif
#endif
  
  // CkPrintf("Returning %d proxies\n", nProxyNodes);

  delete [] proxyNodes;
  return nProxyNodes;
}

#endif
