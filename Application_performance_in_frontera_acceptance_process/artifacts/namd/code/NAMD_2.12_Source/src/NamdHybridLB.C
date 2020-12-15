/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdHybridLB.C,v $
 * $Author: jim $
 * $Date: 2013/08/29 02:20:57 $
 * $Revision: 1.39 $
 *****************************************************************************/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include "InfoStream.h"
#include "NamdHybridLB.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

class SplitComputesMsg : public CMessage_SplitComputesMsg {
public:
  double maxUnsplit;
  double averageLoad;
  double avgCompute;
  double maxCompute;
  int maxComputeId;
  int nMoveableComputes;
  int numPesAvailable;
  int n;
  int *cid;
  float *load;
};

#include "NamdHybridLB.def.h"

// #define DUMP_LDBDATA 1
// #define LOAD_LDBDATA 1

extern int isPmeProcessor(int); 
#ifdef MEM_OPT_VERSION
extern int isOutputProcessor(int); 
#endif
// Load array defined in NamdCentLB.C
extern double *cpuloads;

/**
 * Creates the chare array for the hybrid load balancer.
 */ 
void CreateNamdHybridLB() {
	CProxy_NamdHybridLB::ckNew();

	// creating an array to store the loads of all processors
	// to be used with proxy spanning tree
	if (CkMyPe() == 0 && cpuloads == NULL) {
		cpuloads = new double[CkNumPes()];
		CmiAssert(cpuloads != NULL);
		for (int i=0; i<CkNumPes(); i++) cpuloads[i] = 0.0;
	}
}

/**
 * @brief Default constructor.
 */
NamdHybridLB::NamdHybridLB(): HybridBaseLB(CkLBOptions(-1))
{
  // setting the name
  lbname = (char *)"NamdHybridLB";

  delete tree;        // delete the tree built from the base class
  const SimParameters* simParams = Node::Object()->simParameters;
  if (CkNumPes() <= simParams->hybridGroupSize)  {
    tree = new TwoLevelTree;   // similar to centralized load balancing
  }
  else {
    tree = new ThreeLevelTree(simParams->hybridGroupSize);
    initTree();
    // can only do shrink strategy on levels > 1
    statsStrategy = SHRINK_NULL;
  }

  // initializing thisProxy
  thisProxy = CProxy_NamdHybridLB(thisgroup);
  
  // initializing the central LB
  centralLB = AllocateNamdCentLB();

  // initializing the dummy LB
  dummyLB = AllocateNamdDummyLB();

  // assigning initial values to variables
  from_procs = NULL;
  computeArray = NULL;
  patchArray = NULL;
  processorArray = NULL;
  updateCount = 0;
  splitCount = 0;
  splitComputesMsgs = 0;
  updateFlag = false;
  collectFlag = false;

}

/**
 * @brief Function used to discover if load balancer can be run at that point.
 *
 * It is called from HybridBase every time AtSync method is called.
 */
bool NamdHybridLB::QueryBalanceNow(int _step){ 
  if ( LdbCoordinator::Object()->takingLdbData ) {
	  return true;
  } else {
	  return false;
  } 
}

bool NamdHybridLB::QueryDumpData() {
#if 0                                                                                             
  if (LdbCoordinator::Object()->ldbCycleNum == 1)  return true;                                
  if (LdbCoordinator::Object()->ldbCycleNum == 2)  return true;                                
#endif                                                                                            
  return false;                                                                                
}

#if 0
/**
 *  Runs the load balancing strategy with shrinking the load information.
 *  Note: now, it is just calling the Strategy, but eventually will have
 *  its own code.
 */
LBVectorMigrateMsg* NamdHybridLB::VectorStrategy(LDStats* stats){
  CkPrintf("[%d] Using Vector Strategy to balance the load\n",CkMyPe());
  LBVectorMigrateMsg* msg = new(0,0) LBVectorMigrateMsg;
  msg->n_moves = 0;
  msg->level = currentLevel;
  return msg;
}
#endif

/*
 * Runs the load balancing strategy
 */
CLBMigrateMsg* NamdHybridLB::Strategy(LDStats* stats)
{
	int i;
  	// CkPrintf("[%d] NamdHybridLB at Strategy\n",CkMyPe());
	
	// calling the centralLB for level 1		
	if(currentLevel == 1){
		LevelData *lData = levelData[currentLevel];
		CLBMigrateMsg *msg;
		LocalLBInfoMsg *newMsg;
		msg = GrpLevelStrategy(stats);

		// creating a new message to send to its parent
		newMsg = new(msg->n_moves,endPE-startPE+1) LocalLBInfoMsg;
		newMsg->n_moves = msg->n_moves;
		newMsg->startPE = startPE;
		newMsg->endPE = endPE;
		for(i=0; i<msg->n_moves; i++){
			newMsg->moves[i] = msg->moves[i];
		}
		for(i=0; i<endPE-startPE+1; i++){
			newMsg->cpuloads[i] = peLoads[i];
		}
		delete [] peLoads;
		thisProxy[0].UpdateLocalLBInfo(newMsg);
		return msg;
	}else{
		dummyLB->work(stats);
		return createMigrateMsg(stats);
	}
}

/**
 * Updates the compute map with the migration information from its children.
 */
void NamdHybridLB::UpdateLocalLBInfo(LocalLBInfoMsg *msg){
	int children;
	int i;

	// getting the number of children
	children = tree->numNodes(currentLevel);
	// CkPrintf("[%d] Updating compute map, total %d\n",CkMyPe(),siblings);

	// getting the compute map to insert the changes coming from the children
	ComputeMap *computeMap = ComputeMap::Object();

	// traversing the set of moves in msg
	for(i=0; i<msg->n_moves; i++){
	    if (msg->moves[i].to_pe != -1)
		computeMap->setNewNode(msg->moves[i].obj.id.id[0],msg->moves[i].to_pe);	
	}

	// CODING
	// updating cpuloads array
	for(i=msg->startPE; i<=msg->endPE; i++){
		cpuloads[i] = msg->cpuloads[i-msg->startPE];
	}

	// checking if all children have sent the update
	updateCount++;
	if(updateCount == children){
		updateCount = 0;
		updateFlag = true;
		 // CkPrintf("[%d] UPDATE READY\n",CkMyPe());		
	}

	// checking if the collect info is ready
	if(updateFlag && collectFlag){
		updateFlag = false;
		collectFlag = false;	
		thisProxy[parent_backup].CollectInfo(loc_backup, n_backup, fromlevel_backup);
	}

        delete msg;
}

void NamdHybridLB::splitComputes(SplitComputesMsg *msg) {
  const int children = tree->numNodes(1);

  if ( ! splitComputesMsgs ) {
    splitComputesMsgs = new SplitComputesMsg*[children];
  }
  
  splitComputesMsgs[splitCount] = msg;

  if ( ++splitCount == children ) {
    splitCount = 0;

    const SimParameters* simParams = Node::Object()->simParameters;
    ComputeMap *computeMap = ComputeMap::Object();

    double maxUnsplit = 0.;
    double averageLoad = 0.;
    double avgCompute = 0.;
    double maxCompute = 0.;
    int maxComputeId = -1;
    int nMoveableComputes = 0;
    int numPesAvailable = 0;
    int nToSplit = 0;

    for ( int j=0; j < children; ++j ) {
      SplitComputesMsg *msg = splitComputesMsgs[j];
      if ( msg->maxUnsplit > maxUnsplit ) { maxUnsplit = msg->maxUnsplit; }
      if ( msg->maxCompute > maxCompute ) { maxCompute = msg->maxCompute; maxComputeId = msg->maxComputeId; }
      averageLoad += msg->averageLoad * msg->numPesAvailable;
      numPesAvailable += msg->numPesAvailable;
      avgCompute += msg->avgCompute * msg->nMoveableComputes;
      nMoveableComputes += msg->nMoveableComputes;
      nToSplit += msg->n;
    }
    
    averageLoad /= numPesAvailable;
    avgCompute /= nMoveableComputes;

    CkPrintf("LDB: Largest compute %d load %f is %.1f%% of average load %f\n",
            maxComputeId, maxCompute, 100. * maxCompute / averageLoad, averageLoad);
    CkPrintf("LDB: Average compute %f is %.1f%% of average load %f\n",
            avgCompute, 100. * avgCompute / averageLoad, averageLoad);

    if ( ! nToSplit ) {
      for ( int j=0; j < children; ++j ) {
        delete splitComputesMsgs[j];
      }
    } else {
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
      maxCompute = averageLoad / 10.;
      if ( maxCompute < 2. * avgCompute ) maxCompute = 2. * avgCompute;
      if ( simParams->ldbRelativeGrainsize > 0. ) {
        maxCompute = averageLoad * simParams->ldbRelativeGrainsize;
      }
      CkPrintf("LDB: Partitioning computes with target load %f\n", maxCompute);

      for ( int j=0; j < children; ++j ) {
        SplitComputesMsg *msg = splitComputesMsgs[j];
        for (int i=0; i < msg->n; ++i) {
          int nparts = (int) ceil(msg->load[i] / maxCompute);
          if ( nparts > maxParts ) nparts = maxParts;
          if ( nparts < 1 ) nparts = 1;
          computeMap->setNewNumPartitions(msg->cid[i],nparts);
          totalAddedParts += nparts - 1;
        }
        delete msg;
      }

      CkPrintf("LDB: Increased migratable compute count from %d to %d\n",
              nMoveableComputes,nMoveableComputes+totalAddedParts);
      CkPrintf("LDB: Largest unpartitionable compute is %f\n", maxUnsplit);
    }
  }
}


/**
 * This function implements a strategy similar to the one used in the 
 * centralized case in NamdCentLB.
 */
CLBMigrateMsg* NamdHybridLB::GrpLevelStrategy(LDStats* stats) {
  int numProcessors = stats->nprocs();	// number of processors at group level
  int numPatches = PatchMap::Object()->numPatches();
  ComputeMap *computeMap = ComputeMap::Object();
  const int numComputes = computeMap->numComputes();
  const int numGroupComputes = stats->n_migrateobjs;
  const SimParameters* simParams = Node::Object()->simParameters;

  if ( ! processorArray ) processorArray = new processorInfo[numProcessors];
  // these data structures are global and need to be distributed
  if ( ! patchArray ) patchArray = new patchInfo[numPatches];
  if ( ! computeArray ) computeArray = new computeInfo[numGroupComputes];
  if ( ! from_procs ) from_procs = new int[numGroupComputes];

  int nMoveableComputes = buildData(stats);
  CmiAssert(nMoveableComputes <= numGroupComputes);


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
  double maxCompute;
  int maxComputeId;
  int numPesAvailable;
  {
   int i;
   double total = 0.;
   maxCompute = 0.;
   int maxi = 0;
   for (i=0; i<nMoveableComputes; i++) {
      double load = computeArray[i].load;
      total += load;
      if ( load > maxCompute ) { maxCompute = load;  maxi = i; }
   }
   avgCompute = total / nMoveableComputes;
   maxComputeId = computeArray[maxi].handle.id.id[0];

    int P = stats->nprocs();
   numPesAvailable = 0;
   for (i=0; i<P; i++) {
      if (processorArray[i].available) {
        ++numPesAvailable;
        total += processorArray[i].backgroundLoad;
      }
   }
   if (numPesAvailable == 0)
     NAMD_die("No processors available for load balancing!\n");

   averageLoad = total/numPesAvailable;
  }

  int i_split = 0;
  double maxUnsplit = 0.;

  if ( step() == 1 ) {
    for (int i=0; i<nMoveableComputes; i++) {
      const int cid = computeArray[i].handle.id.id[0];
      if ( computeMap->numPartitions(cid) == 0 ) {
        const double load = computeArray[i].load;
        if ( load > maxUnsplit ) maxUnsplit = load;
        continue;
      }
      ++i_split;
    }
  }

  {
    SplitComputesMsg *msg = new(i_split,i_split) SplitComputesMsg;
    msg->maxUnsplit = maxUnsplit;
    msg->averageLoad = averageLoad;
    msg->avgCompute = avgCompute;
    msg->maxCompute = maxCompute;
    msg->maxComputeId = maxComputeId;
    msg->nMoveableComputes = nMoveableComputes;
    msg->numPesAvailable = numPesAvailable;
    msg->n = i_split;

    if ( step() == 1 ) {
      i_split = 0;
      for (int i=0; i<nMoveableComputes; i++) {
        computeArray[i].processor = computeArray[i].oldProcessor;
        const int cid = computeArray[i].handle.id.id[0];
        if ( computeMap->numPartitions(cid) == 0 ) {
          continue;
        }
        msg->cid[i_split] = cid;
        msg->load[i_split] = computeArray[i].load;
        ++i_split;
      }
    }

    thisProxy[0].splitComputes(msg);
  }

  if ( step() == 1 ) {
    // compute splitting only
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
    NAMD_die("Old load balancer strategy is not compatible with hybrid balancer.");
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
    if (computeArray[i].processor != from_procs[i]+stats->procs[0].pe) {
      /* CkPrintf("[%d] Obj %d migrating from %d (%d) to %d\n",
                     CkMyPe(),computeArray[i].handle.id.id[0],
			 from_procs[i], computeArray[i].oldProcessor, computeArray[i].processor); */
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      //migrateMe->from_pe = computeArray[i].oldProcessor;
      int frompe = from_procs[i];
      if (frompe == numProcessors)
        frompe = -1;
      else
        frompe = frompe + stats->procs[0].pe;
      migrateMe->from_pe = frompe;
      migrateMe->to_pe = computeArray[i].processor;
      if (frompe == -1) {
          // don't know yet which processor this compute belongs to, but
	  // inform receiver
        LDObjData obj;
        obj.handle = computeArray[i].handle;
        thisProxy[computeArray[i].processor].ObjMigrated(obj, NULL, 0, currentLevel-1);
      } 
      migrateInfo.insertAtEnd(migrateMe);

      // sneak in updates to ComputeMap
      //ERASE CkPrintf("%d setting %d to processor %d\n",CkMyPe(),computeArray[i].handle.id.id[0],computeArray[i].processor);
      computeMap->setNewNode(computeArray[i].handle.id.id[0],
				computeArray[i].processor);
    }
  }
  // CkPrintf("LOAD BALANCING READY %d\n",CkMyPe()); 

  LBMigrateMsg* msg;
  msg = createMigrateMsg(migrateInfo, numProcessors);

  peLoads = new double [numProcessors]; 
  startPE = processorArray[0].Id;
  endPE = processorArray[numProcessors-1].Id;
  // CkPrintf("[%d] numProcessors=%d, %d to %d\n",CkMyPe(),numProcessors,processorArray[0].Id,processorArray[numProcessors-1].Id);
  for (i=0; i<numProcessors; i++) {
	peLoads[i] = processorArray[i].load;
  }


  delete [] from_procs;
  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;

  from_procs = NULL;
  processorArray = NULL;
  patchArray = NULL;
  computeArray = NULL;
  
  return msg;

}

void NamdHybridLB::dumpDataASCII(char *file, int numProcessors,
                               int numPatches, int numComputes)
{
  char filename[128];
  sprintf(filename, "%s_%d.%d", file, CkMyPe(), step());
  FILE* fp = fopen(filename,"w");
  if (fp == NULL){
     perror("dumpLDStatsASCII");
     return;
  }
  // CkPrintf("***** DUMP data to file: %s ***** \n", filename);
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


/**
 * @brief Builds the data structures required for the load balancing strategies in NAMD.
 */ 
int NamdHybridLB::buildData(LDStats* stats) {
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
  // traversing the list of processors and getting their load information
  int i, pe_no;
  for (i=0; i<n_pes; ++i) {
    pe_no = stats->procs[i].pe;

    // BACKUP processorArray[i].Id = i; 
    processorArray[i].Id = pe_no;               // absolute pe number
    processorArray[i].available = true;
    // BACKUP if ( pmeOn && isPmeProcessor(i) )
    if ( pmeOn && isPmeProcessor(pe_no) ) {
      processorArray[i].backgroundLoad = pmebgfactor * stats->procs[i].bg_walltime;
    // BACKUP } else if (patchMap->numPatchesOnNode(i) > 0) {
    } else if (patchMap->numPatchesOnNode(pe_no) > 0) {
      processorArray[i].backgroundLoad = homebgfactor * stats->procs[i].bg_walltime;
    } else {
      processorArray[i].backgroundLoad = bgfactor * stats->procs[i].bg_walltime;
    }
    processorArray[i].idleTime = stats->procs[i].idletime;
    processorArray[i].load = processorArray[i].computeLoad = 0.0;
  }

  // If I am group zero, then offload processor 0 and 1 in my group
  if(stats->procs[0].pe == 0) {
    if(unLoadZero) processorArray[0].available = false;
    if(unLoadOne) processorArray[1].available = false;
  }

  // if all pes are Pme, disable this flag
  if (pmeOn && unLoadPme) {
    for (i=0; i<n_pes; i++) {
      if(!isPmeProcessor(stats->procs[i].pe))  break;
    }
    if (i == n_pes) {
      iout << iINFO << "Turned off unLoadPme flag!\n"  << endi;
      unLoadPme = 0;
    }
  }

  if (pmeOn && unLoadPme) {
    for (i=0; i<n_pes; i++) {
      if ((pmeBarrier && i==0) || isPmeProcessor(stats->procs[i].pe)) 
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
      for (i=0; i<n_pes; i++) {
	  if (isOutputProcessor(stats->procs[i].pe)) 
	      processorArray[i].available = false;
      }
  }
#endif

  // need to go over all patches to get all required proxies
  int numPatches = patchMap->numPatches();
  int totalLocalProxies = 0;
  int totalProxies = 0;
  for ( int pid=0; pid<numPatches; ++pid ) {
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

	patchArray[pid].Id = pid;
	patchArray[pid].numAtoms = 0;
	patchArray[pid].processor = patchMap->node(pid);

	const int numProxies = 
#if 0 // USE_TOPOMAP - this function needs to be there for the hybrid case
	requiredProxiesOnProcGrid(pid,neighborNodes);
#else
	requiredProxies(pid, neighborNodes);
#endif

        int numLocalProxies = 0;
	for (int k=0; k<numProxies; k++) {
		if( (neighborNodes[k] >= stats->procs[0].pe) && (neighborNodes[k] <= stats->procs[n_pes-1].pe) ){
			++numLocalProxies;
			int index = neighborNodes[k] - stats->procs[0].pe;
  			processorArray[index].proxies.unchecked_insert(&patchArray[pid]);
  			patchArray[pid].proxiesOn.unchecked_insert(&processorArray[index]);
		}
	}
#if 0
	if ( numLocalProxies ) {
	    CkPrintf("LDB Pe %d patch %d has %d local of %d total proxies\n",
		CkMyPe(), pid, numLocalProxies, numProxies);
	}
#endif
	totalLocalProxies += numLocalProxies;
	totalProxies += numProxies;
  }
#if 0
  CkPrintf("LDB Pe %d has %d local of %d total proxies\n",
		CkMyPe(), totalLocalProxies, totalProxies);
#endif
  
  int nMoveableComputes=0;
  int index;

  int j;

  // this loop goes over only the objects in this group
  for(j=0; j < stats->n_objs; j++) {
	const LDObjData &this_obj = stats->objData[j];
      	int frompe = stats->from_proc[j];

	// filter out non-NAMD managed objects (like PME array)
      	if (this_obj.omID().id.idx != 1) {
                // CmiAssert(frompe>=0 && frompe<n_pes);
                // CkPrintf("non-NAMD object %d on pe %d with walltime %lf\n",
                // this_obj.id().id[0], frompe + stats->procs[0].pe, this_obj.wallTime);
		processorArray[frompe].backgroundLoad += this_obj.wallTime;
        	continue;
	}

      	if (this_obj.id().id[1] == -2) { // Its a patch
		// handled above to get required proxies from all patches
		processorArray[frompe].backgroundLoad += this_obj.wallTime;
	} else if (this_obj.id().id[1] == -3) { // Its a bonded compute
		processorArray[frompe].backgroundLoad += this_obj.wallTime;
	} else if (this_obj.migratable && this_obj.wallTime != 0.) { // Its a compute

		const int cid = this_obj.id().id[0];
		const int p0 = computeMap->pid(cid,0);

		// For self-interactions, just return the same pid twice
		int p1;
		if (computeMap->numPids(cid) > 1)
	  		p1 = computeMap->pid(cid,1);
			else p1 = p0;
			computeArray[nMoveableComputes].Id = cid;
			//BACKUP computeArray[nMoveableComputes].oldProcessor = stats->from_proc[j];
			if (frompe >= n_pes) {  // from outside
CkPrintf("assigning random old processor...this looks broken\n");
			  computeArray[nMoveableComputes].oldProcessor = CrnRand()%n_pes + stats->procs[0].pe;     // random
			}
			else {
			  computeArray[nMoveableComputes].oldProcessor = frompe + stats->procs[0].pe;
			}
			from_procs[nMoveableComputes] = frompe;

			//BACKUP2 index = stats->from_proc[j] - stats->procs[0].pe;
			//BACKUP processorArray[stats->from_proc[j]].computeLoad += this_obj.wallTime;
			int index = computeArray[nMoveableComputes].oldProcessor - stats->procs[0].pe; 
			processorArray[index].computeLoad += this_obj.wallTime;
			computeArray[nMoveableComputes].processor = -1;
			computeArray[nMoveableComputes].patch1 = p0;
			computeArray[nMoveableComputes].patch2 = p1;
			computeArray[nMoveableComputes].handle = this_obj.handle;
			computeArray[nMoveableComputes].load = this_obj.wallTime;
			nMoveableComputes++;
      	}
  }

  	for (i=0; i<n_pes; i++) {
	  processorArray[i].load = processorArray[i].backgroundLoad + processorArray[i].computeLoad;
  	}
  	stats->clear();
  	return nMoveableComputes;
}


int NamdHybridLB::requiredProxies(PatchID id, int neighborNodes[])
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
  int numNodes = CkNumPes();
  int numPatches = patchMap->numPatches();
  int emptyNodes = numNodes - numPatches;
  if ( emptyNodes > numPatches ) {
    int nodesPerPatch = nProxyNodes + 1 + (emptyNodes-1) / numPatches;
    int maxNodesPerPatch = PatchMap::MaxOneAway + PatchMap::MaxTwoAway;
    if ( nodesPerPatch > maxNodesPerPatch ) nodesPerPatch = maxNodesPerPatch;
    int proxyNode = (myNode + 1) % numNodes;
    while ( nProxyNodes < nodesPerPatch &&
			! patchMap->numPatchesOnNode(proxyNode) ) {
      if ( proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
      proxyNode = (proxyNode + 1) % numNodes;
    }
    proxyNode = (myNode - 1 + numNodes) % numNodes;
    while ( nProxyNodes < nodesPerPatch &&
			! patchMap->numPatchesOnNode(proxyNode) ) {
      if ( proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
      proxyNode = (proxyNode - 1 + numNodes) % numNodes;
    }
    proxyNode = (myNode + 1) % numNodes;
    int count = 0;
    while ( nProxyNodes < nodesPerPatch ) {
      if ( ! patchMap->numPatchesOnNode(proxyNode) && proxyNode != myNode ) {
        IF_NEW_NODE {
          neighborNodes[nProxyNodes] = proxyNode;
          nProxyNodes++;
        }
      }
      proxyNode = (proxyNode + 1) % numNodes;
      count ++; if (count == numNodes) break;   // we looped all
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
    if ( proxyNode < numNodes && ! patchMap->numPatchesOnNode(proxyNode) ) {
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

#if 0
void NamdHybridLB::CollectInfo(Location *loc, int n, int fromlevel)
{
   int atlevel = fromlevel + 1;
   LevelData *lData = levelData[atlevel];
   lData->info_recved++;

   CkVec<Location> &matchedObjs = lData->matchedObjs;
CmiAssert(0);

   // sort into mactched and unmatched list
   std::map<LDObjKey, int> &unmatchedObjs = lData->unmatchedObjs;
   for (int i=0; i<n; i++) {
     std::map<LDObjKey, int>::iterator iter = unmatchedObjs.find(loc[i].key);
     if (iter != unmatchedObjs.end()) {
       CmiAssert(iter->second != -1 || loc[i].loc != -1);
       if (loc[i].loc == -1) loc[i].loc = iter->second;
       matchedObjs.push_back(loc[i]);
       unmatchedObjs.erase(iter);
     }
     else
       unmatchedObjs[loc[i].key] = loc[i].loc;
   }

//  DEBUGF(("[%d] level %d has %d unmatched and %d matched. \n", CkMyPe(), atlevel, unmatchedObjs.size(), matchedObjs.size()));

   if (lData->info_recved == lData->nChildren) {
     lData->info_recved = 0;
     if (_lb_args.debug() > 1)
         CkPrintf("[%d] CollectInfo at level %d started at %f\n",
	        CkMyPe(), atlevel, CkWallTimer());
     if (lData->parent != -1) {

		// NAMD specific
		CkVec<Location> unmatchedbuf;
   		for(std::map<LDObjKey, int>::const_iterator it = unmatchedObjs.begin(); it != unmatchedObjs.end(); ++it){
    		unmatchedbuf.push_back(Location(it->first, it->second));
   		}
		// checking if update of ComputeMap is ready before calling parent
		if(CkMyPe() == 0){
			if(updateFlag){
				updateFlag = false;
				collectFlag = false;
				thisProxy[lData->parent].CollectInfo(unmatchedbuf.getVec(), unmatchedbuf.size(), atlevel);
			}else{
				CkPrintf("[%d] COMPUTEMAP UPDATE NOT READY\n",CkMyPe());
				collectFlag = true;
				parent_backup = lData->parent;
				loc_backup = unmatchedbuf.getVec();
				n_backup = unmatchedbuf.size();
				fromlevel_backup = atlevel;
			}
		}else{
			// send only unmatched ones up the tree
			thisProxy[lData->parent].CollectInfo(unmatchedbuf.getVec(), unmatchedbuf.size(), atlevel);
		}

     }
     else { // root
       // we should have all answers now
       CmiAssert(unmatchedObjs.size() == 0);
       // start send match list down
       thisProxy.PropagateInfo(matchedObjs.getVec(), matchedObjs.size(), atlevel, lData->nChildren, lData->children);
       lData->statsData->clear();
     }
   }
}
#endif
