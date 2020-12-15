
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include "InfoStream.h"
#include "NamdNborLB.h"
#include "NamdNborLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

void CreateNamdNborLB()
{
  // CkPrintf("[%d] creating NamdNborLB %d\n",CkMyPe(),loadbalancer);
  CProxy_NamdNborLB::ckNew();
  // CkPrintf("[%d] created NamdNborLB %d\n",CkMyPe(),loadbalancer);
}

NamdNborLB::NamdNborLB(): NeighborLB(CkLBOptions(-1))
{
  //  if (CkMyPe()==0)
  //   CkPrintf("[%d] NamdNborLB created\n",CkMyPe());
  processorArray = 0;
  patchArray = 0;
  computeArray = 0;
  act = 0;
  ldbNum = 0;
}

/*
NamdNborLB::~NamdNborLB()
{
  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;
}
*/

int NamdNborLB::max_neighbors() {
  return 4;
}

int NamdNborLB::num_neighbors() {
  return numNbors;
}

void NamdNborLB::neighbors(int* _n) {
#if 0
    const int me = CkMyPe();
    const int npe = CkNumPes();
    if (npe > 1)
      _n[0] = (me + npe - 1) % npe;
    if (npe > 2)
      _n[1] = (me + 1) % npe;

    int bigstep = (npe - 1) / 3 + 1;
    if (bigstep == 1) bigstep++;

    if (npe > 3)
      _n[2] = (me + bigstep) % npe;
    if (npe > 4)
      _n[3] = (me + npe - bigstep) % npe;
#else	// 2D mesh or Torus mesh
#define SEQ(x, y) ((x)*yDim + (y))
#define WRAP   0
    numNbors = 0;
    int yDim = (int)sqrt((double)CkNumPes());
    int xDim = CkNumPes() / yDim;
    if (CkNumPes() % yDim) xDim++;
    int x = CkMyPe()/yDim;
    int y = CkMyPe()%yDim;
    int x1, y1, s;
    // CmiPrintf("[%d]info: %d %d %d %d\n", CkMyPe(), xDim, yDim, x,y);

    x1=x; y1 = y-1;
#if WRAP
    if (y1==-1) y1=yDim-1;
    if (SEQ(x1, y1) >= CkNumPes()) s = CkNumPes()-1;
    else s = SEQ(x1, y1);
    if (s != CkMyPe()) _n[numNbors++] = s;
#else
    if (y1 != -1)  _n[numNbors++] = SEQ(x1, y1);
#endif

    x1=x; y1=y+1;
#if WRAP
    if (y1 == yDim || SEQ(x1,y1) >= CkNumPes()) y1=0;
    s = SEQ(x1, y1);
    if (s != _n[numNbors-1] && s != CkMyPe()) _n[numNbors++] = s;
#else
    if (y1 == yDim || SEQ(x1,y1) >= CkNumPes()) ;
    else _n[numNbors++] = SEQ(x1, y1);
#endif

    y1=y; x1=x-1;
#if WRAP
    if (x1==-1) x1=xDim-1;
    if (SEQ(x1, y1) >= CkNumPes()) x1--;
    s = SEQ(x1, y1);
    if (s != CkMyPe()) _n[numNbors++] = s;
#else
    if (x1!=-1) _n[numNbors++] = SEQ(x1, y1);
#endif

    y1=y; x1=x+1;
#if WRAP
    if (x1==xDim || SEQ(x1,y1) >= CkNumPes()) x1=0;
    s = SEQ(x1, y1);
    if (s != _n[numNbors-1] && s != CkMyPe()) _n[numNbors++] = s;
#else
    if (x1==xDim || SEQ(x1,y1) >= CkNumPes()) ;
    else _n[numNbors++] = SEQ(x1,y1);
#endif
    // CmiPrintf("[%d] %d neighbors: %d %d %d %d\n", CkMyPe(), numNbors, _n[0], _n[1], _n[2], _n[3]);
    act = (x+y)%2;
#endif

};

bool NamdNborLB::QueryBalanceNow(int _step)
{
  // CkPrintf("[%d] QueryBalanceNow on step %d: %d\n",CkMyPe(),_step, LdbCoordinator::Object()->takingLdbData);
  if ( LdbCoordinator::Object()->takingLdbData ) {
    return true;
  } else {
    return false;
  }
}

bool NamdNborLB::QueryMigrateStep(int _step)
{
  // CmiPrintf("[%d] QueryMigrateStep %d %d.\n", CkMyPe(), _step, act);
  return (act+ldbNum)%2 == 0;
}

NLBMigrateMsg* NamdNborLB::Strategy(NborBaseLB::LDStats* stats, int count)
{
  //  CkPrintf("LDB:[%d] All statistics received at %f, %f\n",
  //  CkMyPe(), CmiTimer(),CmiWallTimer());
  int i,j;
  int mype = CkMyPe();

  ldbNum ++;

  const int numProcessors = CkNumPes();
  const int numPatches = PatchMap::Object()->numPatches();
  const int numComputes = ComputeMap::Object()->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  int nMoveableComputes = 0;
  for (i=0; i < count+1; i++) {
    LDStats &thisLDStats = ((i==count)?myStats:stats[i]);
    for (j=0; j < thisLDStats.n_objs; j++) {
      const LDObjData &this_obj = thisLDStats.objData[j];
      if (this_obj.omID().id.idx != 1) continue;
      if (this_obj.id().id[1] == -2) continue;
      if (this_obj.migratable)  nMoveableComputes++;
    }
  }
  // CmiPrintf("%d nMoveableComputes: %d\n", CkMyPe(), nMoveableComputes);

  // these sizes should never change
  processorArray = new processorInfo[numProcessors];
  patchArray = new patchInfo[numPatches];
//  if ( ! computeArray ) computeArray = new computeInfo[nMoveableComputes];
  computeArray = new computeInfo[nMoveableComputes];

  nMoveableComputes = buildData(stats,count);

  //CmiPrintf("AlgNbor begin on %d\n", CkMyPe());
  AlgNbor(mype, computeArray,patchArray,processorArray,
			nMoveableComputes, numPatches, numProcessors, count);
  //CmiPrintf("AlgNbor end on %d\n", CkMyPe());

  CkVec<MigrateInfo *> migrateInfo;
  for(i=0;i<nMoveableComputes;i++) {
    if (computeArray[i].oldProcessor == mype)
    if (computeArray[i].processor != computeArray[i].oldProcessor) {
      // CkPrintf("[%d] Obj %d migrating from %d to %d\n",
      //          CkMyPe(),computeArray[i].handle.id.id[0],
      //       computeArray[i].processor,computeArray[i].oldProcessor);
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      migrateMe->from_pe = computeArray[i].oldProcessor;
      migrateMe->to_pe = computeArray[i].processor;
      migrateInfo.insertAtEnd(migrateMe);
    }
  }
  
  int migrate_count=migrateInfo.length();
  CkPrintf("NamdNborLB [%d] migrating %d elements\n", CkMyPe(), migrate_count);
  NLBMigrateMsg* msg = new(migrate_count,CkNumPes(),CkNumPes(),0) NLBMigrateMsg;
  msg->n_moves = migrate_count;
  for(i=0; i < migrate_count; i++) {
    MigrateInfo* item = migrateInfo[i];
    msg->moves[i] = *item;
    delete item;
    migrateInfo[i] = 0;
  }

  delete [] patchArray;
  delete [] computeArray;
  /*
  for(i=0; i<numProcessors; i++)
      delete [] processorArray[i].proxyUsage;
  */
  delete [] processorArray;

  return msg;
};


int NamdNborLB::buildData(NborBaseLB::LDStats* stats, int count)
{
  PatchMap* patchMap = PatchMap::Object();
  ComputeMap* computeMap = ComputeMap::Object();
  double bg_weight = 0.7;

  int i;
  for (i=0; i<CkNumPes(); i++) {
    processorArray[i].load = 0.0;
    processorArray[i].backgroundLoad = 0.0;
    processorArray[i].available = false;
    if (i == CkMyPe()) {
      processorArray[i].Id = i;
      processorArray[i].available = myStats.move;
    if (patchMap->numPatches() > 0)
      processorArray[i].backgroundLoad = myStats.bg_walltime*bg_weight;
    else 
      processorArray[i].backgroundLoad = myStats.bg_walltime;
      continue;
    }
    int peslot = NeighborIndex(i);
    if (peslot != -1) {
    processorArray[i].Id = i;
      processorArray[i].available = stats[peslot].move;
    if (patchMap->numPatches() > 0)
      processorArray[i].backgroundLoad = bg_weight * stats[peslot].bg_walltime;
    else 
      processorArray[i].backgroundLoad = stats[peslot].bg_walltime;
    }
    else 
      processorArray[i].Id = -2;
  }

  int nMoveableComputes=0;
  int nProxies = 0;		// total number of estimated proxies
  for (i=0; i<patchMap->numPatches(); i++) {
	patchArray[i].Id = i;
	patchArray[i].numAtoms = 0;
	patchArray[i].processor = patchMap->node(i);
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
	const int numProxies = requiredProxies(i,neighborNodes);
        nProxies += numProxies;

	for (int k=0; k<numProxies; k++) {
	  if (NeighborIndex(neighborNodes[k]) != -1) {
	    processorArray[neighborNodes[k]].proxies.unchecked_insert(&patchArray[i]);
	    patchArray[i].proxiesOn.unchecked_insert(&processorArray[neighborNodes[k]]);
	  }
	}
  }
  for (i=0; i < count+1; i++) {
    int j;
    LDStats &thisLDStats = ((i==count)?myStats:stats[i]);
    for (j=0; j < thisLDStats.n_objs; j++) {
      const LDObjData &this_obj = thisLDStats.objData[j];
      // filter out non-NAMD managed objects (like PME array)
      if (this_obj.omID().id.idx != 1) continue;
      if (this_obj.id().id[1] == -2) { // Its a patch
/*
	const int pid = this_obj.id.id[0];
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

	patchArray[pid].Id = pid;
	patchArray[pid].numAtoms = 0;
	patchArray[pid].processor = i;
	const int numProxies = requiredProxies(pid,neighborNodes);
        nProxies += numProxies;

	for (int k=0; k<numProxies; k++) {
	  processorArray[neighborNodes[k]].proxies.unchecked_insert(&patchArray[pid]);
	  patchArray[pid].proxiesOn.unchecked_insert(&processorArray[neighborNodes[k]]);
	}
*/
      } else if (this_obj.migratable) { // Its a compute
	const int cid = this_obj.id().id[0];
	const int p0 = computeMap->pid(cid,0);

	// For self-interactions, just return the same pid twice
	int p1;
	if (computeMap->numPids(cid) > 1)
	  p1 = computeMap->pid(cid,1);
	else p1 = p0;
	computeArray[nMoveableComputes].Id = cid;
	computeArray[nMoveableComputes].oldProcessor = thisLDStats.from_pe;
	computeArray[nMoveableComputes].processor = -1;
	computeArray[nMoveableComputes].patch1 = p0;
	computeArray[nMoveableComputes].patch2 = p1;
	computeArray[nMoveableComputes].handle = this_obj.handle;
	computeArray[nMoveableComputes].load = this_obj.wallTime;
	nMoveableComputes++;
      }
    }
  }
  return nMoveableComputes;
}

// Figure out which proxies we will definitely create on other
// nodes, without regard for non-bonded computes.  This code is swiped
// from ProxyMgr, and changes there probable need to be propagated here.

int NamdNborLB::requiredProxies(PatchID id, int neighborNodes[])
{
  enum proxyHere { No, Yes };
  int numNodes = CkNumPes();
  proxyHere *proxyNodes = new proxyHere[numNodes];
  int nProxyNodes;
  int i;

  // Note all home patches.
  for ( i = 0; i < numNodes; ++i )
  {
    proxyNodes[i] = No;
  }
  nProxyNodes=0;

  // Check all two-away neighbors.
  // This is really just one-away neighbors, since 
  // two-away always returns zero: RKB
  PatchID neighbors[1 + PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

  PatchMap* patchMap = PatchMap::Object();

  int myNode = patchMap->node(id);
  neighbors[0] = id;
  int numNeighbors = 1 + patchMap->downstreamNeighbors(id,neighbors+1);
  for ( i = 0; i < numNeighbors; ++i )
  {
    const int proxyNode = patchMap->basenode(neighbors[i]);
    if (proxyNode != myNode)
      if (proxyNodes[proxyNode] == No)
      {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
  }

  delete [] proxyNodes;
  return nProxyNodes;
}

