
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HomePatch owns the actual atoms of a Patch of space
   Proxy(s) get messages via ProxyMgr from HomePatch(es)
   to update lists of atoms and their coordinates
   HomePatch(es) also have a Sequencer bound to them

   superclass: 	Patch		
*/

#include "time.h"
#include <math.h>
#include "charm++.h"
#include "qd.h"

#include "SimParameters.h"
#include "HomePatch.h"
#include "AtomMap.h"
#include "Node.h"
#include "PatchMap.inl"
#include "main.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "Migration.h"
#include "Molecule.h"
#include "PatchMgr.h"
#include "Sequencer.h"
#include "Broadcasts.h"
#include "LdbCoordinator.h"
#include "ReductionMgr.h"
#include "Sync.h"
#include "Random.h"
#include "Priorities.h"
#include "ComputeNonbondedUtil.h"
#include "ComputeGBIS.inl"
#include "Priorities.h"
#include "SortAtoms.h"

#include "ComputeQM.h"
#include "ComputeQMMgr.decl.h"

//#define PRINT_COMP
#define TINY 1.0e-20;
#define MAXHGS 10
#define MIN_DEBUG_LEVEL 2
//#define DEBUGM
#include "Debug.h"

#include <vector>
#include <algorithm>
using namespace std;

typedef int HGArrayInt[MAXHGS];
typedef BigReal HGArrayBigReal[MAXHGS];
typedef zVector HGArrayVector[MAXHGS];
typedef BigReal HGMatrixBigReal[MAXHGS][MAXHGS];
typedef zVector HGMatrixVector[MAXHGS][MAXHGS];

#include "ComputeNonbondedMICKernel.h"

int average(CompAtom *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial);

void mollify(CompAtom *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab);


// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0

// Macro to test if a hydrogen group represents a water molecule.
// NOTE: This test is the same test in Molecule.C for setting the
//   OxygenAtom flag in status.
// hgtype should be the number of atoms in a water hydrogen group
// It must now be set based on simulation parameters because we might
// be using tip4p

// DJH: This will give false positive for full Drude model,
//      e.g. O D H is not water but has hgs==3
#define IS_HYDROGEN_GROUP_WATER(hgs, mass)                 \
  ((hgs >= 3) && ((mass >= 14.0) && (mass <= 18.0)))

#endif


HomePatch::HomePatch(PatchID pd, int atomCnt) : Patch(pd)
// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0
  ,tempAtom()
#endif
{
  settle_initialized = 0;

  doAtomUpdate = true;
  rattleListValid = false;

  exchange_msg = 0;
  exchange_req = -1;

  //tracking the end of gbis phases
  numGBISP1Arrived = 0;
  numGBISP2Arrived = 0;
  numGBISP3Arrived = 0;
  phase1BoxClosedCalled = false;
  phase2BoxClosedCalled = false;
  phase3BoxClosedCalled = false;

  min.x = PatchMap::Object()->min_a(patchID);
  min.y = PatchMap::Object()->min_b(patchID);
  min.z = PatchMap::Object()->min_c(patchID);
  max.x = PatchMap::Object()->max_a(patchID);
  max.y = PatchMap::Object()->max_b(patchID);
  max.z = PatchMap::Object()->max_c(patchID);
  center = 0.5*(min+max);

  int aAway = PatchMap::Object()->numaway_a();
  if ( PatchMap::Object()->periodic_a() ||
       PatchMap::Object()->gridsize_a() > aAway + 1 ) {
    aAwayDist = (max.x - min.x) * aAway;
  } else {
    aAwayDist = Node::Object()->simParameters->patchDimension;
  }
  int bAway = PatchMap::Object()->numaway_b();
  if ( PatchMap::Object()->periodic_b() ||
       PatchMap::Object()->gridsize_b() > bAway + 1 ) {
    bAwayDist = (max.y - min.y) * bAway;
  } else {
    bAwayDist = Node::Object()->simParameters->patchDimension;
  }
  int cAway = PatchMap::Object()->numaway_c();
  if ( PatchMap::Object()->periodic_c() ||
       PatchMap::Object()->gridsize_c() > cAway + 1 ) {
    cAwayDist = (max.z - min.z) * cAway;
  } else {
    cAwayDist = Node::Object()->simParameters->patchDimension;
  }

  migrationSuspended = false;
  allMigrationIn = false;
  marginViolations = 0;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
  inMigration = false;
  numMlBuf = 0;
  flags.sequence = -1;
  flags.maxForceUsed = -1;

  numAtoms = atomCnt;
  replacementForces = 0;

  SimParameters *simParams = Node::Object()->simParameters;
  doPairlistCheck_newTolerance = 
	0.5 * ( simParams->pairlistDist - simParams->cutoff );


  numFixedAtoms = 0;
  //if ( simParams->fixedAtomsOn ) {
  //  for ( int i = 0; i < numAtoms; ++i ) {
  //    numFixedAtoms += ( atom[i].atomFixed ? 1 : 0 );
  //  }
  //}

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  ptnTree.resize(0);
  /*children = NULL;
  numChild = 0;*/
#else
  child =  new int[proxySpanDim];
  nChild = 0;	// number of proxy spanning tree children
#endif

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  nphs = 0;
  localphs = NULL;
  isProxyChanged = 0;
#endif


  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0

    // Create the scratch memory for separating atoms
    tempAtom.resize(numAtoms);
    numWaterAtoms = -1;

  #endif
  // Handle unusual water models here
  if (simParams->watmodel == WAT_TIP4) init_tip4();
  else if (simParams->watmodel == WAT_SWM4) init_swm4();

  isNewProxyAdded = 0;

}

HomePatch::HomePatch(PatchID pd, FullAtomList &al) : Patch(pd)
// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0
  ,tempAtom()
#endif
{ 
  atom.swap(al);
  settle_initialized = 0;

  doAtomUpdate = true;
  rattleListValid = false;

  exchange_msg = 0;
  exchange_req = -1;

  numGBISP1Arrived = 0;
  numGBISP2Arrived = 0;
  numGBISP3Arrived = 0;
  phase1BoxClosedCalled = false;
  phase2BoxClosedCalled = false;
  phase3BoxClosedCalled = false;

  min.x = PatchMap::Object()->min_a(patchID);
  min.y = PatchMap::Object()->min_b(patchID);
  min.z = PatchMap::Object()->min_c(patchID);
  max.x = PatchMap::Object()->max_a(patchID);
  max.y = PatchMap::Object()->max_b(patchID);
  max.z = PatchMap::Object()->max_c(patchID);
  center = 0.5*(min+max);

  int aAway = PatchMap::Object()->numaway_a();
  if ( PatchMap::Object()->periodic_a() ||
       PatchMap::Object()->gridsize_a() > aAway + 1 ) {
    aAwayDist = (max.x - min.x) * aAway;
  } else {
    aAwayDist = Node::Object()->simParameters->patchDimension;
  }
  int bAway = PatchMap::Object()->numaway_b();
  if ( PatchMap::Object()->periodic_b() ||
       PatchMap::Object()->gridsize_b() > bAway + 1 ) {
    bAwayDist = (max.y - min.y) * bAway;
  } else {
    bAwayDist = Node::Object()->simParameters->patchDimension;
  }
  int cAway = PatchMap::Object()->numaway_c();
  if ( PatchMap::Object()->periodic_c() ||
       PatchMap::Object()->gridsize_c() > cAway + 1 ) {
    cAwayDist = (max.z - min.z) * cAway;
  } else {
    cAwayDist = Node::Object()->simParameters->patchDimension;
  }

  migrationSuspended = false;
  allMigrationIn = false;
  marginViolations = 0;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
  inMigration = false;
  numMlBuf = 0;
  flags.sequence = -1;
  flags.maxForceUsed = -1;

  numAtoms = atom.size();
  replacementForces = 0;

  SimParameters *simParams = Node::Object()->simParameters;
  doPairlistCheck_newTolerance = 
	0.5 * ( simParams->pairlistDist - simParams->cutoff );


  numFixedAtoms = 0;
  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      numFixedAtoms += ( atom[i].atomFixed ? 1 : 0 );
    }
  }

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  ptnTree.resize(0);
  /*children = NULL;
  numChild = 0;*/
#else
  child =  new int[proxySpanDim];
  nChild = 0;	// number of proxy spanning tree children
#endif

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  nphs = 0;
  localphs = NULL;
  isProxyChanged = 0;
#endif


  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0

    // Create the scratch memory for separating atoms
    tempAtom.resize(numAtoms);
    numWaterAtoms = -1;

    // Separate the current list of atoms
    separateAtoms();

  #endif
    
  // Handle unusual water models here
  if (simParams->watmodel == WAT_TIP4) init_tip4();
  else if (simParams->watmodel == WAT_SWM4) init_swm4();

  isNewProxyAdded = 0;
}

void HomePatch::write_tip4_props() {
  printf("Writing r_om and r_ohc: %f | %f\n", r_om, r_ohc);
}

void HomePatch::init_tip4() {
  // initialize the distances needed for the tip4p water model
  Molecule *mol = Node::Object()->molecule;
  r_om = mol->r_om;
  r_ohc = mol->r_ohc;
}


void ::HomePatch::init_swm4() {
  // initialize the distances needed for the SWM4 water model
  Molecule *mol = Node::Object()->molecule;
  r_om = mol->r_om;
  r_ohc = mol->r_ohc;
}


void HomePatch::reinitAtoms(FullAtomList &al) {
  atomMapper->unregisterIDsFullAtom(atom.begin(),atom.end());

  atom.swap(al);
  numAtoms = atom.size();

  doAtomUpdate = true;
  rattleListValid = false;

  if ( ! numNeighbors ) atomMapper->registerIDsFullAtom(atom.begin(),atom.end());

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0

    // Reset the numWaterAtoms value
    numWaterAtoms = -1;

    // Separate the atoms
    separateAtoms();

  #endif
}

// Bind a Sequencer to this HomePatch
void HomePatch::useSequencer(Sequencer *sequencerPtr)
{ sequencer=sequencerPtr; }

// start simulation over this Patch of atoms
void HomePatch::runSequencer(void)
{ sequencer->run(); }

void HomePatch::readPatchMap() {
  // iout << "Patch " << patchID << " has " << proxy.size() << " proxies.\n" << endi;
  PatchMap *p = PatchMap::Object();
  PatchID nnPatchID[PatchMap::MaxOneAway];

  patchMigrationCounter = numNeighbors 
    = PatchMap::Object()->oneAwayNeighbors(patchID, nnPatchID);
  DebugM( 1, "NumNeighbors for pid " <<patchID<<" is "<< numNeighbors << "\n");
  int n;
  for (n=0; n<numNeighbors; n++) {
    realInfo[n].destNodeID = p->node(realInfo[n].destPatchID = nnPatchID[n]);
     DebugM( 1, " nnPatchID=" <<nnPatchID[n]<<" nnNodeID="<< realInfo[n].destNodeID<< "\n");
    realInfo[n].mList.resize(0);
  }

  // Make mapping from the 3x3x3 cube of pointers to real migration info
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
      {
	int pid =  p->pid(p->index_a(patchID)+i-1, 
	    p->index_b(patchID)+j-1, p->index_c(patchID)+k-1);
	if (pid < 0) {
	   DebugM(5, "ERROR, for patchID " << patchID <<" I got neigh pid = " << pid << "\n");
	}
	if (pid == patchID && ! (
		( (i-1) && p->periodic_a() ) ||
		( (j-1) && p->periodic_b() ) ||
		( (k-1) && p->periodic_c() ) )) {
	  mInfo[i][j][k] = NULL;
	}
	else {
	  // Does not work as expected for periodic with only two patches.
	  // Also need to check which image we want, but OK for now.  -JCP
	  for (n = 0; n<numNeighbors; n++) {
	    if (pid == realInfo[n].destPatchID) {
	      mInfo[i][j][k] = &realInfo[n];
	      break;
	    }
	  }
	  if (n == numNeighbors) { // disaster! 
	    DebugM(4,"BAD News, I could not find PID " << pid << "\n");
	  }
	}
      }

  DebugM(1,"Patch("<<patchID<<") # of neighbors = " << numNeighbors << "\n");
}

HomePatch::~HomePatch()
{
    atomMapper->unregisterIDsFullAtom(atom.begin(),atom.end());
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    ptnTree.resize(0);
    #ifdef USE_NODEPATCHMGR
    delete [] nodeChildren;    
    #endif
#endif
  delete [] child;
}


void HomePatch::boxClosed(int box) {
  // begin gbis
  if (box == 5) {// end of phase 1
    phase1BoxClosedCalled = true;
    if (!psiSumBox.isOpen() && numGBISP1Arrived == proxy.size()) {
      if (flags.doGBIS && flags.doNonbonded) {
        sequencer->awaken();
      }
    } else {
      //need to wait until proxies arrive before awakening
    }
  } else if (box == 6) {// intRad
    //do nothing
  } else if (box == 7) {// bornRad
    //do nothing
  } else if (box == 8) {// end of phase 2
    phase2BoxClosedCalled = true;
    //if no proxies, AfterP1 can't be called from receive
    //so it will be called from here
    if (!dEdaSumBox.isOpen() && numGBISP2Arrived == proxy.size()) {
      if (flags.doGBIS && flags.doNonbonded) {
        sequencer->awaken();
      }
    } else {
      //need to wait until proxies arrive before awakening
    }
  } else if (box == 9) {
    //do nothing
  } else if (box == 10) {
    //lcpoType Box closed: do nothing
  } else {
    //do nothing
  }
  // end gbis

  if ( ! --boxesOpen ) {
    if ( replacementForces ) {
      for ( int i = 0; i < numAtoms; ++i ) {
        if ( replacementForces[i].replace ) {
          for ( int j = 0; j < Results::maxNumForces; ++j ) { f[j][i] = 0; }
          f[Results::normal][i] = replacementForces[i].force;
        }
      }
      replacementForces = 0;
    }
    DebugM(1,patchID << ": " << CthSelf() << " awakening sequencer "
	<< sequencer->thread << "(" << patchID << ") @" << CmiTimer() << "\n");
    // only awaken suspended threads.  Then say it is suspended.

    phase3BoxClosedCalled = true;
    if (flags.doGBIS) {
      if (flags.doNonbonded) {
        sequencer->awaken();
      } else {
        if (numGBISP1Arrived == proxy.size() &&
          numGBISP2Arrived == proxy.size() &&
          numGBISP3Arrived == proxy.size()) {
          sequencer->awaken();//all boxes closed and all proxies arrived
        }
      }
    } else {//non-gbis awaken
      sequencer->awaken();
    }
  } else {
    DebugM(1,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void HomePatch::registerProxy(RegisterProxyMsg *msg) {
  DebugM(4, "registerProxy("<<patchID<<") - adding node " <<msg->node<<"\n");
  proxy.add(msg->node);
  forceBox.clientAdd();

  isNewProxyAdded = 1;
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  isProxyChanged = 1;
#endif

  Random((patchID + 37) * 137).reorder(proxy.begin(),proxy.size());
  delete msg;
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  isProxyChanged = 1;
#endif
  int n = msg->node;
  NodeID *pe = proxy.begin();
  for ( ; *pe != n; ++pe );
  forceBox.clientRemove();
  proxy.del(pe - proxy.begin());
  delete msg;
}

#if USE_TOPOMAP && USE_SPANNING_TREE

int HomePatch::findSubroots(int dim, int* subroots, int psize, int* pidscopy){
  int nChild = 0;
  int cones[6][proxySpanDim*proxySpanDim+proxySpanDim];
  int conesizes[6] = {0,0,0,0,0,0};
  int conecounters[6] = {0,0,0,0,0,0};
  int childcounter = 0;
  nChild = (psize>proxySpanDim)?proxySpanDim:psize;
  TopoManager tmgr;
  for(int i=0;i<psize;i++){
    int cone = tmgr.getConeNumberForRank(pidscopy[i]);
    cones[cone][conesizes[cone]++] = pidscopy[i];
  }

  while(childcounter<nChild){
    for(int i=0;i<6;i++){
      if(conecounters[i]<conesizes[i]){
        subroots[childcounter++] = cones[i][conecounters[i]++];
      }
    }
  }
  for(int i=nChild;i<proxySpanDim;i++)
    subroots[i] = -1;
  return nChild;
}
#endif // USE_TOPOMAP 

static int compDistance(const void *a, const void *b)
{
  int d1 = abs(*(int *)a - CkMyPe());
  int d2 = abs(*(int *)b - CkMyPe());
  if (d1 < d2) 
    return -1;
  else if (d1 == d2) 
    return 0;
  else 
    return 1;
}

void HomePatch::sendProxies()
{
#if USE_NODEPATCHMGR
	CProxy_NodeProxyMgr pm(CkpvAccess(BOCclass_group).nodeProxyMgr);
	NodeProxyMgr *npm = pm[CkMyNode()].ckLocalBranch();
	npm->sendProxyList(patchID, proxy.begin(), proxy.size());
#else
	ProxyMgr::Object()->sendProxies(patchID, proxy.begin(), proxy.size());
#endif
}

#ifdef NODEAWARE_PROXY_SPANNINGTREE
void HomePatch::buildNodeAwareSpanningTree(void){
#if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
	DebugFileTrace *dft = DebugFileTrace::Object();
	dft->openTrace();
	dft->writeTrace("HomePatch[%d] has %d proxy on proc[%d] node[%d]\n", patchID, proxy.size(), CkMyPe(), CkMyNode());
	dft->writeTrace("Proxies are: ");
	for(int i=0; i<proxy.size(); i++) dft->writeTrace("%d(%d), ", proxy[i], CkNodeOf(proxy[i]));
	dft->writeTrace("\n");
	dft->closeTrace();
#endif
 
    //build the naive spanning tree for this home patch
    if(! proxy.size()) {
        //this case will not happen in practice.
        //In debugging state where spanning tree is enforced, then this could happen
        //Chao Mei        
       return;
    }
    ProxyMgr::buildSinglePatchNodeAwareSpanningTree(patchID, proxy, ptnTree);
    //optimize on the naive spanning tree

    //setup the children
    setupChildrenFromProxySpanningTree();
    //send down to children
    sendNodeAwareSpanningTree();
}

void HomePatch::setupChildrenFromProxySpanningTree(){
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    isProxyChanged = 1;
#endif
    if(ptnTree.size()==0) {
        nChild = 0;
        delete [] child;
        child = NULL;
        #ifdef USE_NODEPATCHMGR
        numNodeChild = 0;
        delete [] nodeChildren;
        nodeChildren = NULL;        
        #endif
        return;
    }
    proxyTreeNode *rootnode = &ptnTree.item(0);
    CmiAssert(rootnode->peIDs[0] == CkMyPe());
    //set up children
    //1. add external children (the first proc inside the proxy tree node)    
    //2. add internal children (with threshold that would enable spanning    
    int internalChild = rootnode->numPes-1;
    int externalChild = ptnTree.size()-1;
    externalChild = (proxySpanDim>externalChild)?externalChild:proxySpanDim;
    int internalSlots = proxySpanDim-externalChild;
    if(internalChild>0){
      if(internalSlots==0) {
         //at least having one internal child
        internalChild = 1;
      }else{
        internalChild = (internalSlots>internalChild)?internalChild:internalSlots;
      }
    }
    
    nChild = externalChild+internalChild;
    CmiAssert(nChild>0);

    //exclude the root node        
    delete [] child;
    child = new int[nChild];    

    for(int i=0; i<externalChild; i++) {
        child[i] = ptnTree.item(i+1).peIDs[0];
    }
    for(int i=externalChild, j=1; i<nChild; i++, j++) {
        child[i] = rootnode->peIDs[j];
    }

#ifdef USE_NODEPATCHMGR
    //only register the cores that have proxy patches. The HomePach's core
    //doesn't need to be registered.
    CProxy_NodeProxyMgr pm(CkpvAccess(BOCclass_group).nodeProxyMgr);
    NodeProxyMgr *npm = pm[CkMyNode()].ckLocalBranch();
    if(rootnode->numPes==1){
        npm->registerPatch(patchID, 0, NULL);        
    }
    else{
        npm->registerPatch(patchID, rootnode->numPes-1, &rootnode->peIDs[1]);        
    }

    //set up childrens in terms of node ids
    numNodeChild = externalChild;
    if(internalChild) numNodeChild++;
    delete [] nodeChildren;
    nodeChildren = new int[numNodeChild];    
    for(int i=0; i<externalChild; i++) {
        nodeChildren[i] = ptnTree.item(i+1).nodeID;        
    }
    //the last entry always stores this node id if there are proxies
    //on other cores of the same node for this patch.
    if(internalChild)
        nodeChildren[numNodeChild-1] = rootnode->nodeID;
#endif
    
#if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("HomePatch[%d] has %d children: ", patchID, nChild);
    for(int i=0; i<nChild; i++)
        dft->writeTrace("%d ", child[i]);
    dft->writeTrace("\n");
    dft->closeTrace();
#endif   
}
#endif

#ifdef NODEAWARE_PROXY_SPANNINGTREE
//This is not an entry method, but takes an argument of message type
void HomePatch::recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg){
    //set up the whole tree ptnTree
    int treesize = msg->numNodesWithProxies;    
    ptnTree.resize(treesize);    
    int *pAllPes = msg->allPes;
    for(int i=0; i<treesize; i++) {
        proxyTreeNode *oneNode = &ptnTree.item(i);
        delete [] oneNode->peIDs;
        oneNode->numPes = msg->numPesOfNode[i];
        oneNode->nodeID = CkNodeOf(*pAllPes);
        oneNode->peIDs = new int[oneNode->numPes];
        for(int j=0; j<oneNode->numPes; j++) {
            oneNode->peIDs[j] = *pAllPes;
            pAllPes++;
        }
    }
    //setup children
    setupChildrenFromProxySpanningTree();
    //send down to children
    sendNodeAwareSpanningTree();
}

void HomePatch::sendNodeAwareSpanningTree(){
    if(ptnTree.size()==0) return;    
    ProxyNodeAwareSpanningTreeMsg *msg = 
        ProxyNodeAwareSpanningTreeMsg::getANewMsg(patchID, CkMyPe(), ptnTree.begin(), ptnTree.size());
   
    #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    msg->printOut("HP::sendST");
    #endif

    CmiAssert(CkMyPe() == msg->allPes[0]);
    ProxyMgr::Object()->sendNodeAwareSpanningTree(msg);

}
#else
void HomePatch::recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg){}
void HomePatch::sendNodeAwareSpanningTree(){}
#endif

#ifndef NODEAWARE_PROXY_SPANNINGTREE
// recv a spanning tree from load balancer
void HomePatch::recvSpanningTree(int *t, int n)
{
  int i;
  nChild=0;
  tree.resize(n);
  for (i=0; i<n; i++) {
    tree[i] = t[i];
  }

  for (i=1; i<=proxySpanDim; i++) {
    if (tree.size() <= i) break;
    child[i-1] = tree[i];
    nChild++;
  }

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
  isProxyChanged = 1;
#endif

  // send down to children
  sendSpanningTree();
}

void HomePatch::sendSpanningTree()
{
  if (tree.size() == 0) return;
  ProxySpanningTreeMsg *msg = new ProxySpanningTreeMsg;
  msg->patch = patchID;
  msg->node = CkMyPe();
  msg->tree.copy(tree);  // copy data for thread safety
  ProxyMgr::Object()->sendSpanningTree(msg);  
}
#else
void HomePatch::recvSpanningTree(int *t, int n){}
void HomePatch::sendSpanningTree(){}
#endif

#ifndef NODEAWARE_PROXY_SPANNINGTREE
void HomePatch::buildSpanningTree(void)
{
  nChild = 0;
  int psize = proxy.size();
  if (psize == 0) return;
  NodeIDList oldtree;  oldtree.swap(tree);
  int oldsize = oldtree.size();
  tree.resize(psize + 1);
  tree.setall(-1);
  tree[0] = CkMyPe();
  int s=1, e=psize+1;
  NodeIDList::iterator pli;
  int patchNodesLast =
    ( PatchMap::Object()->numNodesWithPatches() < ( 0.7 * CkNumPes() ) );
  int nNonPatch = 0;
#if 1
    // try to put it to the same old tree
  for ( pli = proxy.begin(); pli != proxy.end(); ++pli )
  {
    int oldindex = oldtree.find(*pli);
    if (oldindex != -1 && oldindex < psize) {
      tree[oldindex] = *pli;
    }
  }
  s=1; e=psize;
  for ( pli = proxy.begin(); pli != proxy.end(); ++pli )
  {
    if (tree.find(*pli) != -1) continue;    // already assigned
    if ( patchNodesLast && PatchMap::Object()->numPatchesOnNode(*pli) ) {
      while (tree[e] != -1) { e--; if (e==-1) e = psize; }
      tree[e] = *pli;
    } else {
      while (tree[s] != -1) { s++; if (s==psize+1) s = 1; }
      tree[s] = *pli;
      nNonPatch++;
    }
  }
#if 1
  if (oldsize==0 && nNonPatch) {
    // first time, sort by distance
    qsort(&tree[1], nNonPatch, sizeof(int), compDistance);
  }
#endif

  //CkPrintf("home: %d:(%d) %d %d %d %d %d\n", patchID, tree.size(),tree[0],tree[1],tree[2],tree[3],tree[4]);

#if USE_TOPOMAP && USE_SPANNING_TREE
  
  //Right now only works for spanning trees with two levels
  int *treecopy = new int [psize];
  int subroots[proxySpanDim];
  int subsizes[proxySpanDim];
  int subtrees[proxySpanDim][proxySpanDim];
  int idxes[proxySpanDim];
  int i = 0;

  for(i=0;i<proxySpanDim;i++){
    subsizes[i] = 0;
    idxes[i] = i;
  }
  
  for(i=0;i<psize;i++){
    treecopy[i] = tree[i+1];
  }
  
  TopoManager tmgr;
  tmgr.sortRanksByHops(treecopy,nNonPatch);
  tmgr.sortRanksByHops(treecopy+nNonPatch,
						psize-nNonPatch);  
  
  /* build tree and subtrees */
  nChild = findSubroots(proxySpanDim,subroots,psize,treecopy);
  delete [] treecopy;
  
  for(int i=1;i<psize+1;i++){
    int isSubroot=0;
    for(int j=0;j<nChild;j++)
      if(tree[i]==subroots[j]){
        isSubroot=1;
	break;
      }
    if(isSubroot) continue;
    
    int bAdded = 0;
    tmgr.sortIndexByHops(tree[i], subroots,
						  idxes, proxySpanDim);
    for(int j=0;j<proxySpanDim;j++){
      if(subsizes[idxes[j]]<proxySpanDim){
        subtrees[idxes[j]][(subsizes[idxes[j]])++] = tree[i];
	bAdded = 1; 
        break;
      }
    }
    if( psize > proxySpanDim && ! bAdded ) {
      NAMD_bug("HomePatch BGL Spanning Tree error: Couldn't find subtree for leaf\n");
    }
  }

#else /* USE_TOPOMAP && USE_SPANNING_TREE */
  
  for (int i=1; i<=proxySpanDim; i++) {
    if (tree.size() <= i) break;
    child[i-1] = tree[i];
    nChild++;
  }
#endif
#endif
  
#if 0
  // for debugging
  CkPrintf("[%d] Spanning tree for %d with %d children %d nNonPatch %d\n", CkMyPe(), patchID, psize, nNonPatch);
  for (int i=0; i<psize+1; i++) {
    CkPrintf("%d ", tree[i]);
  }
  CkPrintf("\n");
#endif
  // send to children nodes
  sendSpanningTree();
}
#endif


void HomePatch::receiveResults(ProxyResultVarsizeMsg *msg) {

    numGBISP3Arrived++;
    DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
    int n = msg->node;
    Results *r = forceBox.clientOpen();

    char *iszeroPtr = msg->isZero;
    Force *msgFPtr = msg->forceArr;

    for ( int k = 0; k < Results::maxNumForces; ++k )
    {
      Force *rfPtr = r->f[k];
      for(int i=0; i<msg->flLen[k]; i++, rfPtr++, iszeroPtr++) {
          if((*iszeroPtr)!=1) {
              *rfPtr += *msgFPtr;
              msgFPtr++;
          }
      }      
    }
    forceBox.clientClose();
    delete msg;
}

void HomePatch::receiveResults(ProxyResultMsg *msg) {
  numGBISP3Arrived++;
  DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
  int n = msg->node;
  Results *r = forceBox.clientOpen();
  for ( int k = 0; k < Results::maxNumForces; ++k )
  {
    Force *f = r->f[k];
    register ForceList::iterator f_i, f_e;
    f_i = msg->forceList[k]->begin();
    f_e = msg->forceList[k]->end();
    for ( ; f_i != f_e; ++f_i, ++f ) *f += *f_i;
  }
  forceBox.clientClose();
  delete msg;
}

void HomePatch::receiveResults(ProxyCombinedResultRawMsg* msg)
{
    numGBISP3Arrived++;
  DebugM(4, "patchID("<<patchID<<") receiveRes() #nodes("<<msg->nodeSize<<")\n");
    Results *r = forceBox.clientOpen(msg->nodeSize);
      register char* isNonZero = msg->isForceNonZero;
      register Force* f_i = msg->forceArr;
      for ( int k = 0; k < Results::maxNumForces; ++k )
      {
        Force *f = r->f[k];
		int nf = msg->flLen[k];
#ifdef ARCH_POWERPC
#pragma disjoint (*f_i, *f)
#endif
        for (int count = 0; count < nf; count++) {
          if(*isNonZero){
			f[count].x += f_i->x;
			f[count].y += f_i->y;
			f[count].z += f_i->z;
			f_i++;
          }
          isNonZero++;
        }
      }
    forceBox.clientClose(msg->nodeSize);

    delete msg;
}

void HomePatch::qmSwapAtoms()
{
    // This is used for LSS in QM/MM simulations.
    // Changes atom labels so that we effectively exchange solvent
    // molecules between classical and quantum modes.
    
    SimParameters *simParams = Node::Object()->simParameters;
    int numQMAtms = Node::Object()->molecule->get_numQMAtoms();
    const Real * const qmAtomGroup = Node::Object()->molecule->get_qmAtomGroup() ;
    const int *qmAtmIndx = Node::Object()->molecule->get_qmAtmIndx() ;
    Real *qmAtmChrg = Node::Object()->molecule->get_qmAtmChrg() ;
    
    ComputeQMMgr *mgrP = CProxy_ComputeQMMgr::ckLocalBranch(
            CkpvAccess(BOCclass_group).computeQMMgr) ;
    
    FullAtom *a_i = atom.begin();
    
    for (int i=0; i<numAtoms; ++i ) { 
        
        LSSSubsDat *subP = lssSubs(mgrP).find( LSSSubsDat(a_i[i].id) ) ;
        
        if ( subP != NULL ) {
            a_i[i].id = subP->newID ;
            a_i[i].vdwType = subP->newVdWType ;
            
            // If we are swappign a classical atom with a QM one, the charge
            // will need extra handling.
            if (qmAtomGroup[subP->newID] > 0 && simParams->PMEOn) {
                // We make sure that the last atom charge calculated for the 
                // QM atom being transfered here is available for PME
                // in the next step.
                
                // Loops over all QM atoms (in all QM groups) comparing their 
                // global indices (the sequential atom ID from NAMD).
                for (int qmIter=0; qmIter<numQMAtms; qmIter++) {
                    
                    if (qmAtmIndx[qmIter] == subP->newID) {
                        qmAtmChrg[qmIter] = subP->newCharge;
                        break;
                    }
                    
                }
                
                // For QM atoms, the charge in the full atom structure is zero.
                // Electrostatic interactions between QM atoms and their 
                // environment is handled in ComputeQM.
                a_i[i].charge = 0;
            }
            else {
                // If we are swappign a QM atom with a Classical one, only the charge
                // in the full atomstructure needs updating, since it used to be zero.
                a_i[i].charge = subP->newCharge ;
            }
        }
    }
    
    return;
}

void HomePatch::positionsReady(int doMigration)
{    
  SimParameters *simParams = Node::Object()->simParameters;

  flags.sequence++;

  if (!patchMapRead) { readPatchMap(); }
      
  if (numNeighbors && ! simParams->staticAtomAssignment) {
    if (doMigration) {
      rattleListValid = false;
      doAtomMigration();
    } else {
      doMarginCheck();
    }
  }
  
  if (doMigration && simParams->qmLSSOn)
      qmSwapAtoms();

#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  if ( doMigration ) {
    int n = numAtoms;
    FullAtom *a_i = atom.begin();
    #if defined(NAMD_CUDA) || (defined(NAMD_MIC) && MIC_SORT_ATOMS != 0)
    int *ao = new int[n];
    int nfree;
    if ( simParams->fixedAtomsOn && ! simParams->fixedAtomsForces ) {
      int k = 0;
      int k2 = n;
      for ( int j=0; j<n; ++j ) {
        // put fixed atoms at end
        if ( a_i[j].atomFixed ) ao[--k2] = j;
        else ao[k++] = j;
      }
      nfree = k;
    } else {
      nfree = n;
      for ( int j=0; j<n; ++j ) {
        ao[j] = j;
      }
    }

    sortAtomsForCUDA(ao,a_i,nfree,n);
  
    for ( int i=0; i<n; ++i ) { 
      a_i[i].sortOrder = ao[i];
    }
    delete [] ao;
    #else
      for (int i = 0; i < n; ++i) {
        a_i[i].sortOrder = i;
      }
    #endif
  }

  { 
    const double charge_scaling = sqrt(COULOMB * ComputeNonbondedUtil::scaling * ComputeNonbondedUtil::dielectric_1);
    const Vector ucenter = lattice.unscale(center);
    const BigReal ucenter_x = ucenter.x;
    const BigReal ucenter_y = ucenter.y;
    const BigReal ucenter_z = ucenter.z;
    const int n = numAtoms;
    #if defined(NAMD_MIC) && (MIC_HANDCODE_FORCE_SOA_VS_AOS == 0)
      int n_16 = n;
      n_16 = (n + 15) & (~15);
      cudaAtomList.resize(n_16);
      CudaAtom *ac = cudaAtomPtr = cudaAtomList.begin();
      mic_position_t *atom_x = ((mic_position_t*)ac) + (0 * n_16);
      mic_position_t *atom_y = ((mic_position_t*)ac) + (1 * n_16);
      mic_position_t *atom_z = ((mic_position_t*)ac) + (2 * n_16);
      mic_position_t *atom_q = ((mic_position_t*)ac) + (3 * n_16);
    #else
    cudaAtomList.resize(n);
    CudaAtom *ac = cudaAtomPtr = cudaAtomList.begin();
    #endif
    const FullAtom *a = atom.begin();
    for ( int k=0; k<n; ++k ) {
      #if defined(NAMD_MIC) && (MIC_HANDCODE_FORCE_SOA_VS_AOS == 0)
        int j = a[k].sortOrder;
        atom_x[k] = a[j].position.x - ucenter_x;
        atom_y[k] = a[j].position.y - ucenter_y;
        atom_z[k] = a[j].position.z - ucenter_z;
        atom_q[k] = charge_scaling * a[j].charge;
      #else
      int j = a[k].sortOrder;
      ac[k].x = a[j].position.x - ucenter_x;
      ac[k].y = a[j].position.y - ucenter_y;
      ac[k].z = a[j].position.z - ucenter_z;
      ac[k].q = charge_scaling * a[j].charge;
      #endif
    }
  }
#else
  doMigration = doMigration && numNeighbors;
#endif
  doMigration = doMigration || ! patchMapRead;

  doMigration = doMigration || doAtomUpdate;
  doAtomUpdate = false;

  // Workaround for oversize groups
  doGroupSizeCheck();

  // Copy information needed by computes and proxys to Patch::p.
  p.resize(numAtoms);
  CompAtom *p_i = p.begin();
  pExt.resize(numAtoms);
  CompAtomExt *pExt_i = pExt.begin();
  FullAtom *a_i = atom.begin();
  int i; int n = numAtoms;
  for ( i=0; i<n; ++i ) { 
    p_i[i] = a_i[i]; 
    pExt_i[i] = a_i[i];
  }

  // Measure atom movement to test pairlist validity
  doPairlistCheck();

  if (flags.doMolly) mollyAverage();
  // BEGIN LA
  if (flags.doLoweAndersen) loweAndersenVelocities();
  // END LA

    if (flags.doGBIS) {
      //reset for next time step
      numGBISP1Arrived = 0;
      phase1BoxClosedCalled = false;
      numGBISP2Arrived = 0;
      phase2BoxClosedCalled = false;
      numGBISP3Arrived = 0;
      phase3BoxClosedCalled = false;
      if (doMigration || isNewProxyAdded)
        setGBISIntrinsicRadii();
    }

  if (flags.doLCPO) {
    if (doMigration || isNewProxyAdded) {
      setLcpoType();
    }
  }

  // Must Add Proxy Changes when migration completed!
  NodeIDList::iterator pli;
  int *pids = NULL;
  int pidsPreAllocated = 1;
  int npid;
  if (proxySendSpanning == 0) {
    npid = proxy.size();
    pids = new int[npid];
    pidsPreAllocated = 0;
    int *pidi = pids;
    int *pide = pids + proxy.size();
    int patchNodesLast =
      ( PatchMap::Object()->numNodesWithPatches() < ( 0.7 * CkNumPes() ) );
    for ( pli = proxy.begin(); pli != proxy.end(); ++pli )
    {
      if ( patchNodesLast && PatchMap::Object()->numPatchesOnNode(*pli) ) {
        *(--pide) = *pli;
      } else {
        *(pidi++) = *pli;
      }
    }
  }
  else {
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    #ifdef USE_NODEPATCHMGR
    npid = numNodeChild;
    pids = nodeChildren;
    #else
    npid = nChild;
    pids = child;
    #endif
#else
    npid = nChild;
    pidsPreAllocated = 0;
    pids = new int[proxySpanDim];
    for (int i=0; i<nChild; i++) pids[i] = child[i];
#endif
  }
  if (npid) { //have proxies
    int seq = flags.sequence;
    int priority = PROXY_DATA_PRIORITY + PATCH_PRIORITY(patchID);
    //begin to prepare proxy msg and send it
    int pdMsgPLLen = p.size();
    int pdMsgAvgPLLen = 0;
    if(flags.doMolly) {
        pdMsgAvgPLLen = p_avg.size();
    }
    // BEGIN LA
    int pdMsgVLLen = 0;
    if (flags.doLoweAndersen) {
	pdMsgVLLen = v.size();
    }
    // END LA

    int intRadLen = 0;
    if (flags.doGBIS && (doMigration || isNewProxyAdded)) {
	    intRadLen = numAtoms * 2;
    }

    //LCPO
    int lcpoTypeLen = 0;
    if (flags.doLCPO && (doMigration || isNewProxyAdded)) {
	    lcpoTypeLen = numAtoms;
    }

    int pdMsgPLExtLen = 0;
    if(doMigration || isNewProxyAdded) {
        pdMsgPLExtLen = pExt.size();
    }

    int cudaAtomLen = 0;
#ifdef NAMD_CUDA
    cudaAtomLen = numAtoms;
#endif

    #ifdef NAMD_MIC
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        cudaAtomLen = numAtoms;
      #else
        cudaAtomLen = (numAtoms + 15) & (~15);
      #endif
    #endif

    ProxyDataMsg *nmsg = new (pdMsgPLLen, pdMsgAvgPLLen, pdMsgVLLen, intRadLen,
      lcpoTypeLen, pdMsgPLExtLen, cudaAtomLen, PRIORITY_SIZE) ProxyDataMsg; // BEGIN LA, END LA

    SET_PRIORITY(nmsg,seq,priority);
    nmsg->patch = patchID;
    nmsg->flags = flags;
    nmsg->plLen = pdMsgPLLen;                
    //copying data to the newly created msg
    memcpy(nmsg->positionList, p.begin(), sizeof(CompAtom)*pdMsgPLLen);
    nmsg->avgPlLen = pdMsgAvgPLLen;        
    if(flags.doMolly) {
        memcpy(nmsg->avgPositionList, p_avg.begin(), sizeof(CompAtom)*pdMsgAvgPLLen);
    }
    // BEGIN LA
    nmsg->vlLen = pdMsgVLLen;
    if (flags.doLoweAndersen) {
	memcpy(nmsg->velocityList, v.begin(), sizeof(CompAtom)*pdMsgVLLen);
    }
    // END LA

    if (flags.doGBIS && (doMigration || isNewProxyAdded)) {
      for (int i = 0; i < numAtoms * 2; i++) {
        nmsg->intRadList[i] = intRad[i];
      }
    }

    if (flags.doLCPO && (doMigration || isNewProxyAdded)) {
      for (int i = 0; i < numAtoms; i++) {
        nmsg->lcpoTypeList[i] = lcpoType[i];
      }
    }

    nmsg->plExtLen = pdMsgPLExtLen;
    if(doMigration || isNewProxyAdded){     
        memcpy(nmsg->positionExtList, pExt.begin(), sizeof(CompAtomExt)*pdMsgPLExtLen);
    }

// DMK
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
    memcpy(nmsg->cudaAtomList, cudaAtomPtr, sizeof(CudaAtom)*cudaAtomLen);
#endif
    
#if NAMD_SeparateWaters != 0
    //DMK - Atom Separation (water vs. non-water)
    nmsg->numWaterAtoms = numWaterAtoms;
#endif

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
    nmsg->isFromImmMsgCall = 0;
#endif
    
    #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("HP::posReady: for HomePatch[%d], sending proxy msg to: ", patchID);
    for(int i=0; i<npid; i++) {
        dft->writeTrace("%d ", pids[i]);
    }
    dft->writeTrace("\n");
    dft->closeTrace();
    #endif

#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    if (isProxyChanged || localphs == NULL)
    {
//CmiPrintf("[%d] Build persistent: isProxyChanged: %d %p\n", CkMyPe(), isProxyChanged, localphs);
     //CmiAssert(isProxyChanged);
     if (nphs) {
       for (int i=0; i<nphs; i++) {
         CmiDestoryPersistent(localphs[i]);
       }
       delete [] localphs;
     }
     localphs = new PersistentHandle[npid];
     int persist_size = sizeof(envelope) + sizeof(ProxyDataMsg) + sizeof(CompAtom)*(pdMsgPLLen+pdMsgAvgPLLen+pdMsgVLLen) + intRadLen*sizeof(Real) + lcpoTypeLen*sizeof(int) + sizeof(CompAtomExt)*pdMsgPLExtLen + sizeof(CudaAtom)*cudaAtomLen + PRIORITY_SIZE/8 + 2048;
     for (int i=0; i<npid; i++) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
       if (proxySendSpanning)
           localphs[i] = CmiCreateNodePersistent(pids[i], persist_size, sizeof(envelope)+sizeof(ProxyDataMsg));
       else
#endif
       localphs[i] = CmiCreatePersistent(pids[i], persist_size, sizeof(envelope)+sizeof(ProxyDataMsg));
     }
     nphs = npid;
    }
    CmiAssert(nphs == npid && localphs != NULL);
    CmiUsePersistentHandle(localphs, nphs);
#endif
    if(doMigration || isNewProxyAdded) {
        ProxyMgr::Object()->sendProxyAll(nmsg,npid,pids);
    }else{
        ProxyMgr::Object()->sendProxyData(nmsg,npid,pids);
    }
#if CMK_PERSISTENT_COMM && USE_PERSISTENT_TREE
    CmiUsePersistentHandle(NULL, 0);
#endif
    isNewProxyAdded = 0;
  }
  isProxyChanged = 0;
  if(!pidsPreAllocated) delete [] pids;
  DebugM(4, "patchID("<<patchID<<") doing positions Ready\n");

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  positionPtrBegin = p.begin();
  positionPtrEnd = p.end();
#endif

  if(flags.doMolly) {
      avgPositionPtrBegin = p_avg.begin();
      avgPositionPtrEnd = p_avg.end();
  }
  // BEGIN LA
  if (flags.doLoweAndersen) {
      velocityPtrBegin = v.begin();
      velocityPtrEnd = v.end();
  }
  // END LA

  Patch::positionsReady(doMigration);

  patchMapRead = 1;

  // gzheng
  Sync::Object()->PatchReady();
}

void HomePatch::replaceForces(ExtForce *f)
{
  replacementForces = f;
}

void HomePatch::saveForce(const int ftag)
{
  f_saved[ftag].resize(numAtoms);
  for ( int i = 0; i < numAtoms; ++i )
  {
    f_saved[ftag][i] = f[ftag][i];
  }
}


#undef DEBUG_REDISTRIB_FORCE 
#undef DEBUG_REDISTRIB_FORCE_VERBOSE
//#define DEBUG_REDISTRIB_FORCE
/*
 * Redistribute forces from lone pair site onto its host atoms.
 *
 * Atoms are "labeled" i, j, k, l, where atom i is the lone pair.
 * Positions of atoms are ri, rj, rk, rl.
 * Forces of atoms are fi, fj, fk, fl; updated forces are returned.
 * Accumulate updates to the virial.
 *
 * The forces on the atoms are updated so that:
 *   - the force fi on the lone pair site is 0
 *   - the net force fi+fj+fk+fl is conserved
 *   - the net torque cross(ri,fi)+cross(rj,fj)+cross(rk,fk)+cross(rl,fl)
 *     is conserved
 *
 * If "midpt" is true (nonzero), then use the midpoint of rk and rl
 * (e.g. rk and rl are the hydrogen atoms for water).  Otherwise use
 * planes defined by ri,rj,rk and rj,rk,rl.
 *
 * Having "midpt" set true corresponds in CHARMM to having a negative
 * distance parameter in Lphost.
 *
 * Use FIX_FOR_WATER below to fix the case that occurs when the lone pair
 * site for water lies on the perpendicular bisector of rk and rl, making
 * the cross(r,s) zero.
 */
#define FIX_FOR_WATER
void HomePatch::redistrib_lp_force(
    Vector& fi, Vector& fj, Vector& fk, Vector& fl,
    const Vector& ri, const Vector& rj, const Vector& rk, const Vector& rl,
    Tensor *virial, int midpt) {
#ifdef DEBUG_REDISTRIB_FORCE
  Vector foldnet, toldnet;  // old net force, old net torque
  foldnet = fi + fj + fk + fl;
  toldnet = cross(ri,fi) + cross(rj,fj) + cross(rk,fk) + cross(rl,fl);
#endif
  Vector fja(0), fka(0), fla(0);

  Vector r = ri - rj;
  BigReal invr2 = r.rlength();
  invr2 *= invr2;
  BigReal fdot = (fi*r) * invr2;
  Vector fr = r * fdot;

  fja += fr;

  Vector s, t;
  if (midpt) {
    s = rj - 0.5*(rk + rl);
    t = 0.5*(rk - rl);
  }
  else {
    s = rj - rk;
    t = rk - rl;
  }
  BigReal invs2 = s.rlength();
  invs2 *= invs2;

  Vector p = cross(r,s);
#if !defined(FIX_FOR_WATER)
  BigReal invp = p.rlength();
#else
  BigReal p2 = p.length2();  // fix division by zero above
#endif

  Vector q = cross(s,t);
  BigReal invq = q.rlength();

#if !defined(FIX_FOR_WATER)
  BigReal fpdot = (fi*p) * invp;
  Vector fp = p * fpdot;
  Vector ft = fi - fr - fp;
#else
  BigReal fpdot;
  Vector fp, ft;
  if (p2 < 1e-6) {  // vector is near zero, assume no fp contribution to force
    fpdot = 0;
    fp = 0;
    ft = fi - fr;
  }
  else {
    fpdot = (fi*p) / sqrt(p2);
    fp = p * fpdot;
    ft = fi - fr - fp;
  }
#endif

  fja += ft;
  Vector v = cross(r,ft);  // torque
  ft = cross(s,v) * invs2;
  fja -= ft;

  if (midpt) {
    fka += 0.5 * ft;
    fla += 0.5 * ft;
  }
  else {
    fka += ft;
  }

  BigReal srdot = (s*r) * invs2;
  Vector rr = r - s*srdot;
  BigReal rrdot = rr.length();
  BigReal stdot = (s*t) * invs2;
  Vector tt = t - s*stdot;
  BigReal invtt = tt.rlength();
  BigReal fact = rrdot*fpdot*invtt*invq;
  Vector fq = q * fact;

  fla += fq;
  fja += fp*(1+srdot) + fq*stdot;

  ft = fq*(1+stdot) + fp*srdot;

  if (midpt) {
    fka += -0.5*ft;
    fla += -0.5*ft;
  }
  else {
    fka -= ft;
  }

  if (virial) {
    Tensor va = outer(fja,rj);
    va += outer(fka,rk);
    va += outer(fla,rl);
    va -= outer(fi,ri);
    *virial += va;
  }

  fi = 0;  // lone pair has zero force
  fj += fja;
  fk += fka;
  fl += fla;

#ifdef DEBUG_REDISTRIB_FORCE
#define TOL_REDISTRIB  1e-4
  Vector fnewnet, tnewnet;  // new net force, new net torque
  fnewnet = fi + fj + fk + fl;
  tnewnet = cross(ri,fi) + cross(rj,fj) + cross(rk,fk) + cross(rl,fl);
  Vector fdiff = fnewnet - foldnet;
  Vector tdiff = tnewnet - toldnet;
  if (fdiff.length2() > TOL_REDISTRIB*TOL_REDISTRIB) {
    printf("Error:  force redistribution for water exceeded tolerance:  "
        "fdiff=(%f, %f, %f)\n", fdiff.x, fdiff.y, fdiff.z);
  }
  if (tdiff.length2() > TOL_REDISTRIB*TOL_REDISTRIB) {
    printf("Error:  torque redistribution for water exceeded tolerance:  "
        "tdiff=(%f, %f, %f)\n", tdiff.x, tdiff.y, tdiff.z);
  }
#endif
}


/* Redistribute forces from the massless lonepair charge particle onto
 * the other atoms of the water.
 *
 * This is done using the same algorithm as charmm uses for TIP4P lonepairs.
 *
 * Pass by reference the forces (O H1 H2 LP) to be modified,
 * pass by constant reference the corresponding positions,
 * and a pointer to virial.
 */
void HomePatch::redistrib_lp_water_force(
    Vector& f_ox, Vector& f_h1, Vector& f_h2, Vector& f_lp,
    const Vector& p_ox, const Vector& p_h1, const Vector& p_h2,
    const Vector& p_lp, Tensor *virial) {

#ifdef DEBUG_REDISTRIB_FORCE 
  // Debug information to check against results at end

  // total force and torque relative to origin
  Vector totforce, tottorque;
  totforce = f_ox + f_h1 + f_h2 + f_lp;
  tottorque = cross(f_ox, p_ox) + cross(f_h1, p_h1) + cross(f_h2, p_h2);
  //printf("Torque without LP is %f/%f/%f\n",
  //    tottorque.x, tottorque.y, tottorque.z);
  tottorque += cross(f_lp, p_lp);
  //printf("Torque with LP is %f/%f/%f\n",
  //    tottorque.x, tottorque.y, tottorque.z);
#endif

  // accumulate force adjustments
  Vector fad_ox(0), fad_h(0);

  // Calculate the radial component of the force and add it to the oxygen
  Vector r_ox_lp = p_lp - p_ox;
  BigReal invlen2_r_ox_lp = r_ox_lp.rlength();
  invlen2_r_ox_lp *= invlen2_r_ox_lp;
  BigReal rad_factor = (f_lp * r_ox_lp) * invlen2_r_ox_lp;
  Vector f_rad = r_ox_lp * rad_factor;

  fad_ox += f_rad;

  // Calculate the angular component
  Vector r_hcom_ox = p_ox - ( (p_h1 + p_h2) * 0.5 );
  Vector r_h2_h1_2 = (p_h1 - p_h2) * 0.5; // half of r_h2_h1

  // deviation from collinearity of charge site
  //Vector r_oop = cross(r_ox_lp, r_hcom_ox);
  //
  // vector out of o-h-h plane
  //Vector r_perp = cross(r_hcom_ox, r_h2_h1_2);

  // Here we assume that Ox/Lp/Hcom are linear
  // If you want to correct for deviations, this is the place

//  printf("Deviation from linearity for ox %i: %f/%f/%f\n", oxind, r_oop.x, r_oop.y, r_oop.z);

  Vector f_ang = f_lp - f_rad; // leave the angular component

  // now split this component onto the other atoms
  BigReal len_r_ox_lp = r_ox_lp.length();
  BigReal invlen_r_hcom_ox = r_hcom_ox.rlength();
  BigReal oxcomp = (r_hcom_ox.length() - len_r_ox_lp) * invlen_r_hcom_ox;
  BigReal hydcomp = 0.5 * len_r_ox_lp * invlen_r_hcom_ox;

  fad_ox += (f_ang * oxcomp);
  fad_h += (f_ang * hydcomp);  // adjustment for both hydrogens

  // Add virial contributions
  if (virial) {
    Tensor vir = outer(fad_ox, p_ox);
    vir += outer(fad_h, p_h1);
    vir += outer(fad_h, p_h2);
    vir -= outer(f_lp, p_lp);
    *virial += vir;
  }

  //Vector zerovec(0.0, 0.0, 0.0);
  f_lp = 0;
  f_ox += fad_ox;
  f_h1 += fad_h;
  f_h2 += fad_h;

#ifdef DEBUG_REDISTRIB_FORCE 
  // Check that the total force and torque come out right
  Vector newforce, newtorque;
  newforce = f_ox + f_h1 + f_h2;
  newtorque = cross(f_ox, p_ox) + cross(f_h1, p_h1) + cross(f_h2, p_h2);
  Vector fdiff = newforce - totforce;
  Vector tdiff = newtorque - tottorque;
  BigReal error = fdiff.length();
  if (error > 0.0001) {
     printf("Error:  Force redistribution for water "
         "exceeded force tolerance:  error=%f\n", error);
  }
#ifdef DEBUG_REDISTRIB_FORCE_VERBOSE
  printf("Error in net force:  %f\n", error);
#endif

  error = tdiff.length();
  if (error > 0.0001) {
     printf("Error:  Force redistribution for water "
         "exceeded torque tolerance:  error=%f\n", error);
  }
#ifdef DEBUG_REDISTRIB_FORCE_VERBOSE
  printf("Error in net torque:  %f\n", error);
#endif
#endif /* DEBUG */
}

void HomePatch::reposition_lonepair(
    Vector& ri, const Vector& rj, const Vector& rk, const Vector& rl,
    Real distance, Real angle, Real dihedral)
{
  if ( (rj-rk).length2() > 100. || (rj-rl).length2() > 100. ) {
    iout << iWARN << "Large distance between lonepair reference atoms: ("
      << rj << ") (" << rk << ") (" << rl << ")\n" << endi;
  }
  BigReal r, t, p, cst, snt, csp, snp, invlen;
  Vector v, w, a, b, c;

  if (distance >= 0) {
    v = rk;
    r = distance;
  }
  else {
    v = 0.5*(rk + rl);
    r = -distance;
  }

  t = angle;
  p = dihedral;
  cst = cos(t);
  snt = sin(t);
  csp = cos(p);
  snp = sin(p);
  a = v - rj;
  b = rl - v;
  invlen = a.rlength();
  a *= invlen;
  c = cross(b, a);
  invlen = c.rlength();
  c *= invlen;
  b = cross(a, c);
  w.x = r*cst;
  w.y = r*snt*csp;
  w.z = r*snt*snp;
  ri.x = rj.x + w.x*a.x + w.y*b.x + w.z*c.x;
  ri.y = rj.y + w.x*a.y + w.y*b.y + w.z*c.y;
  ri.z = rj.z + w.x*a.z + w.y*b.z + w.z*c.z;
}

void HomePatch::reposition_all_lonepairs(void) {
  // ASSERT: simParams->lonepairs == TRUE
  for (int i=0;  i < numAtoms;  i++) {
    if (atom[i].mass < 0.01) {
      // found a lone pair
      AtomID aid = atom[i].id;  // global atom ID of lp
      Lphost *lph = Node::Object()->molecule->get_lphost(aid);  // its lphost
      if (lph == NULL) {
        char errmsg[512];
        sprintf(errmsg, "reposition lone pairs: "
            "no Lphost exists for LP %d\n", aid);
        NAMD_die(errmsg);
      }
      LocalID j = AtomMap::Object()->localID(lph->atom2);
      LocalID k = AtomMap::Object()->localID(lph->atom3);
      LocalID l = AtomMap::Object()->localID(lph->atom4);
      if (j.pid != patchID || k.pid != patchID || l.pid != patchID) {
        char errmsg[512];
        sprintf(errmsg, "reposition lone pairs: "
            "LP %d has some Lphost atom off patch\n", aid);
        NAMD_die(errmsg);
      }
      // reposition this lone pair
      reposition_lonepair(atom[i].position, atom[j.index].position,
          atom[k.index].position, atom[l.index].position,
          lph->distance, lph->angle, lph->dihedral);
    }
  }
}

void HomePatch::swm4_omrepos(Vector *ref, Vector *pos, Vector *vel,
    BigReal invdt) {
  // Reposition lonepair (Om) particle of Drude SWM4 water.
  // Same comments apply as to tip4_omrepos(), but the ordering of atoms
  // is different: O, D, LP, H1, H2.
  pos[2] = pos[0] + (0.5 * (pos[3] + pos[4]) - pos[0]) * (r_om / r_ohc);
  // Now, adjust velocity of particle to get it to appropriate place
  // during next integration "drift-step"
  if (invdt != 0) {
    vel[2] = (pos[2] - ref[2]) * invdt;
  }
  // No virial correction needed since lonepair is massless
}

void HomePatch::tip4_omrepos(Vector* ref, Vector* pos, Vector* vel, BigReal invdt) {
  /* Reposition the om particle of a tip4p water
   * A little geometry shows that the appropriate position is given by
   * R_O + (1 / 2 r_ohc) * ( 0.5 (R_H1 + R_H2) - R_O ) 
   * Here r_om is the distance from the oxygen to Om site, and r_ohc
   * is the altitude from the oxygen to the hydrogen center of mass
   * Those quantities are precalculated upon initialization of HomePatch
   *
   * Ordering of TIP4P atoms: O, H1, H2, LP.
   */

  //printf("rom/rohc are %f %f and invdt is %f\n", r_om, r_ohc, invdt);
  //printf("Other positions are: \n  0: %f %f %f\n  1: %f %f %f\n  2: %f %f %f\n", pos[0].x, pos[0].y, pos[0].z, pos[1].x, pos[1].y, pos[1].z, pos[2].x, pos[2].y, pos[2].z);
  pos[3] = pos[0] + (0.5 * (pos[1] + pos[2]) - pos[0]) * (r_om / r_ohc); 
  //printf("New position for lp is %f %f %f\n", pos[3].x, pos[3].y, pos[3].z);

  // Now, adjust the velocity of the particle to get it to the appropriate place
  if (invdt != 0) {
    vel[3] = (pos[3] - ref[3]) * invdt;
  }

  // No virial correction needed, since this is a massless particle
  return;
}

void HomePatch::redistrib_lonepair_forces(const int ftag, Tensor *virial) {
  // ASSERT: simParams->lonepairs == TRUE
  ForceList *f_mod = f;
  for (int i = 0;  i < numAtoms;  i++) {
    if (atom[i].mass < 0.01) {
      // found a lone pair
      AtomID aid = atom[i].id;  // global atom ID of lp
      Lphost *lph = Node::Object()->molecule->get_lphost(aid);  // its lphost
      if (lph == NULL) {
        char errmsg[512];
        sprintf(errmsg, "redistrib lone pair forces: "
            "no Lphost exists for LP %d\n", aid);
        NAMD_die(errmsg);
      }
      LocalID j = AtomMap::Object()->localID(lph->atom2);
      LocalID k = AtomMap::Object()->localID(lph->atom3);
      LocalID l = AtomMap::Object()->localID(lph->atom4);
      if (j.pid != patchID || k.pid != patchID || l.pid != patchID) {
        char errmsg[512];
        sprintf(errmsg, "redistrib lone pair forces: "
            "LP %d has some Lphost atom off patch\n", aid);
        NAMD_die(errmsg);
      }
      // redistribute forces from this lone pair
      int midpt = (lph->distance < 0);
      redistrib_lp_force(f_mod[ftag][i], f_mod[ftag][j.index],
          f_mod[ftag][k.index], f_mod[ftag][l.index],
          atom[i].position, atom[j.index].position,
          atom[k.index].position, atom[l.index].position, virial, midpt);
    }
  }
}

void HomePatch::redistrib_swm4_forces(const int ftag, Tensor *virial) {
  // Loop over the patch's atoms and apply the appropriate corrections
  // to get all forces off of lone pairs
  ForceList *f_mod = f;
  for (int i = 0;  i < numAtoms;  i++) {
    if (atom[i].mass < 0.01) {
      // found lonepair
      redistrib_lp_water_force(f_mod[ftag][i-2], f_mod[ftag][i+1],
          f_mod[ftag][i+2], f_mod[ftag][i],
          atom[i-2].position, atom[i+1].position,
          atom[i+2].position, atom[i].position, virial);
    }
  }
}

void HomePatch::redistrib_tip4p_forces(const int ftag, Tensor* virial) {
  // Loop over the patch's atoms and apply the appropriate corrections
  // to get all forces off of lone pairs
  // Atom ordering:  O H1 H2 LP

  ForceList *f_mod =f;
  for (int i=0; i<numAtoms; i++) {
    if (atom[i].mass < 0.01) {
      // found lonepair
      redistrib_lp_water_force(f_mod[ftag][i-3], f_mod[ftag][i-2],
          f_mod[ftag][i-1], f_mod[ftag][i],
          atom[i-3].position, atom[i-2].position,
          atom[i-1].position, atom[i].position, virial);
    }
  }
}


void HomePatch::addForceToMomentum(const BigReal timestep, const int ftag,
							const int useSaved)
{
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  ForceList *f_use = (useSaved ? f_saved : f);

  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      if ( atom[i].atomFixed ) {
        atom[i].velocity = 0;
      } else {
        BigReal recip_val = ( atom[i].mass > 0. ? dt * namd_reciprocal( atom[i].mass ) : 0.); 
        atom[i].velocity += f_use[ftag][i] * recip_val;
      }
    }
  } else {
    FullAtom *atom_arr  = atom.begin();
    const Force    *force_arr = f_use[ftag].const_begin();
#ifdef ARCH_POWERPC
#pragma disjoint (*force_arr, *atom_arr)
#endif
    for ( int i = 0; i < numAtoms; ++i ) {
      if (atom[i].mass == 0.) continue;
      BigReal recip_val = ( atom[i].mass > 0. ? dt * namd_reciprocal( atom[i].mass ) : 0.); 
      //printf("Taking reciprocal of mass %f\n", atom[i].mass);
      atom_arr[i].velocity.x += force_arr[i].x * recip_val;
      atom_arr[i].velocity.y += force_arr[i].y * recip_val;
      atom_arr[i].velocity.z += force_arr[i].z * recip_val;
    }
  }
}

void HomePatch::addForceToMomentum3(const BigReal timestep1, const int ftag1, const int useSaved1,
    const BigReal timestep2, const int ftag2, const int useSaved2,
    const BigReal timestep3, const int ftag3, const int useSaved3)
{
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt1 = timestep1 / TIMEFACTOR;
  const BigReal dt2 = timestep2 / TIMEFACTOR;
  const BigReal dt3 = timestep3 / TIMEFACTOR;
  ForceList *f_use1 = (useSaved1 ? f_saved : f);
  ForceList *f_use2 = (useSaved2 ? f_saved : f);
  ForceList *f_use3 = (useSaved3 ? f_saved : f);

  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      if ( atom[i].atomFixed ) {
        atom[i].velocity = 0;
      } else {
        BigReal recip_val = ( atom[i].mass > 0. ? namd_reciprocal( atom[i].mass ) : 0.); 
        atom[i].velocity += (f_use1[ftag1][i]*dt1 + f_use2[ftag2][i]*dt2 + f_use3[ftag3][i]*dt3)*recip_val;
      }
    }
  } else {
    FullAtom *atom_arr  = atom.begin();
    const Force *force_arr1 = f_use1[ftag1].const_begin();
    const Force *force_arr2 = f_use2[ftag2].const_begin();
    const Force *force_arr3 = f_use3[ftag3].const_begin();
#ifdef ARCH_POWERPC
#pragma disjoint (*force_arr1, *atom_arr)
#pragma disjoint (*force_arr2, *atom_arr)
#pragma disjoint (*force_arr3, *atom_arr)
#endif
    for ( int i = 0; i < numAtoms; ++i ) {
      if (atom[i].mass == 0.) continue;
      BigReal recip_val = ( atom[i].mass > 0. ? namd_reciprocal( atom[i].mass ) : 0.); 
      //printf("Taking reciprocal of mass %f\n", atom[i].mass);
      atom_arr[i].velocity.x += (force_arr1[i].x*dt1 + force_arr2[i].x*dt2 + force_arr3[i].x*dt3) * recip_val;
      atom_arr[i].velocity.y += (force_arr1[i].y*dt1 + force_arr2[i].y*dt2 + force_arr3[i].y*dt3) * recip_val;
      atom_arr[i].velocity.z += (force_arr1[i].z*dt1 + force_arr2[i].z*dt2 + force_arr3[i].z*dt3) * recip_val;
    }
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      if ( ! atom[i].atomFixed ) atom[i].position += atom[i].velocity * dt;
    }
  } else {
    FullAtom *atom_arr  = atom.begin();
    for ( int i = 0; i < numAtoms; ++i ) {
      atom_arr[i].position.x  +=  atom_arr[i].velocity.x * dt;
      atom_arr[i].position.y  +=  atom_arr[i].velocity.y * dt;
      atom_arr[i].position.z  +=  atom_arr[i].velocity.z * dt;
    }
  }
}

int HomePatch::hardWallDrude(const BigReal timestep, Tensor *virial,
    SubmitReduction *ppreduction)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal kbt=BOLTZMANN*simParams->drudeTemp;
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const BigReal dt = timestep / TIMEFACTOR;
  const BigReal invdt = (dt == 0.) ? 0. : 1.0 / dt; // precalc 1/dt
  int i, ia, ib, j;
  int dieOnError = simParams->rigidDie;
  Tensor wc;  // constraint virial
  BigReal idz, zmin, delta_T, maxtime=timestep,v_Bond;
  int nslabs;

  // start data for hard wall boundary between drude and its host atom
  // static int Count=0;
  int Idx;
  double r_wall, r_wall_SQ, rab, rab_SQ, dr, mass_a, mass_b, mass_sum;
  Vector v_ab, vb_1, vp_1, vb_2, vp_2, new_vel_a, new_vel_b, new_pos_a, new_pos_b, *new_pos, *new_vel;
  double dot_v_r_1, dot_v_r_2;
  double vb_cm, dr_a, dr_b;
  // end data for hard wall boundary between drude and its host atom

  // start calculation of hard wall boundary between drude and its host atom
  if (simParams->drudeHardWallOn) {
    if (ppreduction) {
      nslabs = simParams->pressureProfileSlabs;
      idz = nslabs/lattice.c().z;
      zmin = lattice.origin().z - 0.5*lattice.c().z;
    }

    r_wall = simParams->drudeBondLen;
    r_wall_SQ = r_wall*r_wall;
    // Count++;
    for (i=1; i<numAtoms; i++)	{
      if ( (atom[i].mass > 0.01) && ((atom[i].mass < 1.0)) ) { // drude particle
        ia = i-1;
        ib = i;

        v_ab = atom[ib].position - atom[ia].position;
        rab_SQ = v_ab.x*v_ab.x + v_ab.y*v_ab.y + v_ab.z*v_ab.z;

        if (rab_SQ > r_wall_SQ)	{  // to impose the hard wall constraint
          rab = sqrt(rab_SQ);
          if ( (rab > (2.0*r_wall)) && dieOnError ) {  // unexpected situation
            iout << iERROR << "HardWallDrude> "
              << "The drude is too far away from atom "
              << (atom[ia].id + 1) << " d = " << rab << "!\n" << endi;
            return -1;  // triggers early exit
          }

          v_ab.x /= rab;
          v_ab.y /= rab;
          v_ab.z /= rab;

          if ( fixedAtomsOn && atom[ia].atomFixed ) {  // the heavy atom is fixed
            if (atom[ib].atomFixed) {  // the drude is fixed too
              continue;
            }
            else {  // only the heavy atom is fixed
              dot_v_r_2 = atom[ib].velocity.x*v_ab.x
                + atom[ib].velocity.y*v_ab.y + atom[ib].velocity.z*v_ab.z;
              vb_2 = dot_v_r_2 * v_ab;
              vp_2 = atom[ib].velocity - vb_2;

              dr = rab - r_wall;
              if(dot_v_r_2 == 0.0) {
              	delta_T = maxtime;
              }
              else {
                delta_T = dr/fabs(dot_v_r_2); // the time since the collision occurs
                if(delta_T > maxtime ) delta_T = maxtime; // make sure it is not crazy
              }

              dot_v_r_2 = -dot_v_r_2*sqrt(kbt/atom[ib].mass)/fabs(dot_v_r_2);

              vb_2 = dot_v_r_2 * v_ab;

              new_vel_a = atom[ia].velocity;
              new_vel_b = vp_2 + vb_2;

              dr_b = -dr + delta_T*dot_v_r_2;  // L = L_0 + dT *v_new, v was flipped

              new_pos_a = atom[ia].position;
              new_pos_b = atom[ib].position + dr_b*v_ab; // correct the position
            }
          }
          else {
            mass_a = atom[ia].mass;
            mass_b = atom[ib].mass;
            mass_sum = mass_a+mass_b;

            dot_v_r_1 = atom[ia].velocity.x*v_ab.x
              + atom[ia].velocity.y*v_ab.y + atom[ia].velocity.z*v_ab.z;
            vb_1 = dot_v_r_1 * v_ab;
            vp_1 = atom[ia].velocity - vb_1;

            dot_v_r_2 = atom[ib].velocity.x*v_ab.x
              + atom[ib].velocity.y*v_ab.y + atom[ib].velocity.z*v_ab.z;
            vb_2 = dot_v_r_2 * v_ab;
            vp_2 = atom[ib].velocity - vb_2;

            vb_cm = (mass_a*dot_v_r_1 + mass_b*dot_v_r_2)/mass_sum;

            dot_v_r_1 -= vb_cm;
            dot_v_r_2 -= vb_cm;

            dr = rab - r_wall;

            if(dot_v_r_2 == dot_v_r_1) {
            	delta_T = maxtime;
            }
            else {
              delta_T = dr/fabs(dot_v_r_2 - dot_v_r_1);  // the time since the collision occurs
              if(delta_T > maxtime ) delta_T = maxtime; // make sure it is not crazy
            }
            
            // the relative velocity between ia and ib. Drawn according to T_Drude
            v_Bond = sqrt(kbt/mass_b);

            // reflect the velocity along bond vector and scale down
            dot_v_r_1 = -dot_v_r_1*v_Bond*mass_b/(fabs(dot_v_r_1)*mass_sum);
            dot_v_r_2 = -dot_v_r_2*v_Bond*mass_a/(fabs(dot_v_r_2)*mass_sum);

            dr_a = dr*mass_b/mass_sum + delta_T*dot_v_r_1;
            dr_b = -dr*mass_a/mass_sum + delta_T*dot_v_r_2;

            new_pos_a = atom[ia].position + dr_a*v_ab;	// correct the position
            new_pos_b = atom[ib].position + dr_b*v_ab;
            // atom[ia].position += (dr_a*v_ab);  // correct the position
            // atom[ib].position += (dr_b*v_ab);
            
            dot_v_r_1 += vb_cm;
            dot_v_r_2 += vb_cm;

            vb_1 = dot_v_r_1 * v_ab;
            vb_2 = dot_v_r_2 * v_ab;

            new_vel_a = vp_1 + vb_1;
            new_vel_b = vp_2 + vb_2;
          }

          int ppoffset, partition;
          if ( invdt == 0 ) {
            atom[ia].position = new_pos_a;
            atom[ib].position = new_pos_b;
          }
          else if ( virial == 0 ) {
            atom[ia].velocity = new_vel_a;
            atom[ib].velocity = new_vel_b;
          }
          else {
            for ( j = 0; j < 2; j++ ) {
              if (j==0) {  // atom ia, heavy atom
                Idx = ia;
                new_pos = &new_pos_a;
                new_vel = &new_vel_a;
              }
              else if (j==1) {  // atom ib, drude
                Idx = ib;
                new_pos = &new_pos_b;
                new_vel = &new_vel_b;
              }
              Force df = (*new_vel - atom[Idx].velocity) *
                ( atom[Idx].mass * invdt );
              Tensor vir = outer(df, atom[Idx].position);
              wc += vir;
              atom[Idx].velocity = *new_vel;
              atom[Idx].position = *new_pos;

              if (ppreduction) {
                if (!j) {
                  BigReal z = new_pos->z;
                  int partition = atom[Idx].partition;
                  int slab = (int)floor((z-zmin)*idz);
                  if (slab < 0) slab += nslabs;
                  else if (slab >= nslabs) slab -= nslabs;
                  ppoffset = 3*(slab + nslabs*partition);
                }
                ppreduction->item(ppoffset  ) += vir.xx;
                ppreduction->item(ppoffset+1) += vir.yy;
                ppreduction->item(ppoffset+2) += vir.zz;
              }

            }
          }
        }				
      }
    }

    // if ( (Count>10000) && (Count%10==0) ) {
    //   v_ab = atom[1].position - atom[0].position;
    //   rab_SQ = v_ab.x*v_ab.x + v_ab.y*v_ab.y + v_ab.z*v_ab.z;
    //   iout << "DBG_R: " << Count << "  " << sqrt(rab_SQ) << "\n" << endi;
    // }

  }

  // end calculation of hard wall boundary between drude and its host atom

  if ( dt && virial ) *virial += wc;

  return 0;
}

void HomePatch::buildRattleList() {

  SimParameters *simParams = Node::Object()->simParameters;
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;

  // Re-size to containt numAtoms elements
  velNew.resize(numAtoms);
  posNew.resize(numAtoms);

  // Size of a hydrogen group for water
  int wathgsize = 3;
  int watmodel = simParams->watmodel;
  if (watmodel == WAT_TIP4) wathgsize = 4;
  else if (watmodel == WAT_SWM4) wathgsize = 5;

  // Initialize the settle algorithm with water parameters
  // settle1() assumes all waters are identical,
  // and will generate bad results if they are not.
  // XXX this will move to Molecule::build_atom_status when that 
  // version is debugged
  if ( ! settle_initialized ) {
    for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
      // find a water
      if (atom[ig].rigidBondLength > 0) {
        int oatm;
        if (simParams->watmodel == WAT_SWM4) {
          oatm = ig+3;  // skip over Drude and Lonepair
          //printf("ig=%d  mass_ig=%g  oatm=%d  mass_oatm=%g\n",
          //    ig, atom[ig].mass, oatm, atom[oatm].mass);
        }
        else {
          oatm = ig+1;
          // Avoid using the Om site to set this by mistake
          if (atom[ig].mass < 0.5 || atom[ig+1].mass < 0.5) {
            oatm += 1;
          }
        }

        // initialize settle water parameters
        settle1init(atom[ig].mass, atom[oatm].mass, 
                    atom[ig].rigidBondLength,
                    atom[oatm].rigidBondLength,
                    settle_mOrmT, settle_mHrmT, settle_ra,
                    settle_rb, settle_rc, settle_rra);
        settle_initialized = 1;
        break; // done with init
      }
    }
  }
  
  Vector ref[10];
  BigReal rmass[10];
  BigReal dsq[10];
  int fixed[10];
  int ial[10];
  int ibl[10];

  int numSettle = 0;
  int numRattle = 0;
  int posRattleParam = 0;

  settleList.clear();
  rattleList.clear();
  noconstList.clear();
  rattleParam.clear();

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) {
      // only one atom in group
      noconstList.push_back(ig);
      continue;
    }
    int anyfixed = 0;
    for (int i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      rmass[i] = (atom[ig+i].mass > 0. ? 1. / atom[ig+i].mass : 0.);
      fixed[i] = ( fixedAtomsOn && atom[ig+i].atomFixed );
      if ( fixed[i] ) {
        anyfixed = 1;
        rmass[i] = 0.;
      }
    }
    int icnt = 0;
    BigReal tmp = atom[ig].rigidBondLength;
    if (tmp > 0.0) {  // for water
      if (hgs != wathgsize) {
        char errmsg[256];
        sprintf(errmsg, "Water molecule starting with atom %d contains %d atoms "
                         "but the specified water model requires %d atoms.\n",
                         atom[ig].id+1, hgs, wathgsize);
        NAMD_die(errmsg);
      }
      // Use SETTLE for water unless some of the water atoms are fixed,
      if (useSettle && !anyfixed) {
        // Store to Settle -list
        settleList.push_back(ig);
        continue;
      }
      if ( !(fixed[1] && fixed[2]) ) {
        dsq[icnt] = tmp * tmp;
        ial[icnt] = 1;
        ibl[icnt] = 2;
        ++icnt;
      }
    } // if (tmp > 0.0)
    for (int i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = atom[ig+i].rigidBondLength ) > 0 ) {
        if ( !(fixed[0] && fixed[i]) ) {
          dsq[icnt] = tmp * tmp;
          ial[icnt] = 0;
          ibl[icnt] = i;
          ++icnt;
        }
      }
    }
    if ( icnt == 0 ) {
      // no constraints
      noconstList.push_back(ig);
      continue;  
    }
    // Store to Rattle -list
    RattleList rattleListElem;
    rattleListElem.ig  = ig;
    rattleListElem.icnt = icnt;
    rattleList.push_back(rattleListElem);
    for (int i = 0; i < icnt; ++i ) {
      int a = ial[i];
      int b = ibl[i];
      RattleParam rattleParamElem;
      rattleParamElem.ia = a;
      rattleParamElem.ib = b;
      rattleParamElem.dsq = dsq[i];
      rattleParamElem.rma = rmass[a];
      rattleParamElem.rmb = rmass[b];
      rattleParam.push_back(rattleParamElem);
    }
  }

}

void HomePatch::addRattleForce(const BigReal invdt, Tensor& wc) {
  for (int ig = 0; ig < numAtoms; ++ig ) {
    Force df = (velNew[ig] - atom[ig].velocity) * ( atom[ig].mass * invdt );
    Tensor vir = outer(df, atom[ig].position);
    wc += vir;
    f[Results::normal][ig] += df;
    atom[ig].velocity = velNew[ig];
  }
}

int HomePatch::rattle1(const BigReal timestep, Tensor *virial, 
  SubmitReduction *ppreduction) {

  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->watmodel != WAT_TIP3 || ppreduction) {
    // Call old rattle1 -method instead
    return rattle1old(timestep, virial, ppreduction);
  }

  if (!rattleListValid) {
    buildRattleList();
    rattleListValid = true;
  }

  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;
  const BigReal dt = timestep / TIMEFACTOR;
  const BigReal invdt = (dt == 0.) ? 0. : 1.0 / dt; // precalc 1/dt
  const BigReal tol2 = 2.0 * simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;

  Vector ref[10];  // reference position
  Vector pos[10];  // new position
  Vector vel[10];  // new velocity

  // Manual un-roll
  int n = (settleList.size()/2)*2;
  for (int j=0;j < n;j+=2) {
    int ig;
    ig = settleList[j];
    for (int i = 0; i < 3; ++i ) {
      ref[i] = atom[ig+i].position;
      pos[i] = atom[ig+i].position + atom[ig+i].velocity * dt;
    }
    ig = settleList[j+1];
    for (int i = 0; i < 3; ++i ) {
      ref[i+3] = atom[ig+i].position;
      pos[i+3] = atom[ig+i].position + atom[ig+i].velocity * dt;
    }
    settle1_SIMD<2>(ref, pos,
      settle_mOrmT, settle_mHrmT, settle_ra,
      settle_rb, settle_rc, settle_rra);

    ig = settleList[j];
    for (int i = 0; i < 3; ++i ) {
      velNew[ig+i] = (pos[i] - ref[i])*invdt;
      posNew[ig+i] = pos[i];
    }
    ig = settleList[j+1];
    for (int i = 0; i < 3; ++i ) {
      velNew[ig+i] = (pos[i+3] - ref[i+3])*invdt;
      posNew[ig+i] = pos[i+3];
    }

  }

  if (settleList.size() % 2) {
    int ig = settleList[settleList.size()-1];
    for (int i = 0; i < 3; ++i ) {
      ref[i] = atom[ig+i].position;
      pos[i] = atom[ig+i].position + atom[ig+i].velocity * dt;
    }
    settle1_SIMD<1>(ref, pos,
            settle_mOrmT, settle_mHrmT, settle_ra,
            settle_rb, settle_rc, settle_rra);        
    for (int i = 0; i < 3; ++i ) {
      velNew[ig+i] = (pos[i] - ref[i])*invdt;
      posNew[ig+i] = pos[i];
    }
  }

  int posParam = 0;
  for (int j=0;j < rattleList.size();++j) {

    BigReal refx[10];
    BigReal refy[10];
    BigReal refz[10];

    BigReal posx[10];
    BigReal posy[10];
    BigReal posz[10];

    int ig = rattleList[j].ig;
    int icnt = rattleList[j].icnt;
    int hgs = atom[ig].hydrogenGroupSize;
    for (int i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      pos[i] = atom[ig+i].position;
      if (!(fixedAtomsOn && atom[ig+i].atomFixed))
        pos[i] += atom[ig+i].velocity * dt;
      refx[i] = ref[i].x;
      refy[i] = ref[i].y;
      refz[i] = ref[i].z;
      posx[i] = pos[i].x;
      posy[i] = pos[i].y;
      posz[i] = pos[i].z;
    }


    bool done;
    bool consFailure;
    if (icnt == 1) {
      rattlePair<1>(&rattleParam[posParam],
        refx, refy, refz,
        posx, posy, posz);
      done = true;
      consFailure = false;
    } else {
      rattleN(icnt, &rattleParam[posParam],
        refx, refy, refz,
        posx, posy, posz,
        tol2, maxiter,
        done, consFailure);
    }

    // Advance position in rattleParam
    posParam += icnt;

    for (int i = 0; i < hgs; ++i ) {
      pos[i].x = posx[i];
      pos[i].y = posy[i];
      pos[i].z = posz[i];
    }

    for (int i = 0; i < hgs; ++i ) {
      velNew[ig+i] = (pos[i] - ref[i])*invdt;
      posNew[ig+i] = pos[i];
    }

    if ( consFailure ) {
      if ( dieOnError ) {
        iout << iERROR << "Constraint failure in RATTLE algorithm for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
        return -1;  // triggers early exit
      } else {
        iout << iWARN << "Constraint failure in RATTLE algorithm for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
      }
    } else if ( ! done ) {
      if ( dieOnError ) {
        iout << iERROR << "Exceeded RATTLE iteration limit for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
        return -1;  // triggers early exit
      } else {
        iout << iWARN << "Exceeded RATTLE iteration limit for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
      }
    }

  }

  // Finally, we have to go through atoms that are not involved in rattle just so that we have
  // their positions and velocities up-to-date in posNew and velNew
  for (int j=0;j < noconstList.size();++j) {
    int ig = noconstList[j];
    int hgs = atom[ig].hydrogenGroupSize;
    for (int i = 0; i < hgs; ++i ) {
      velNew[ig+i] = atom[ig+i].velocity;
      posNew[ig+i] = atom[ig+i].position;
    }
  }

  if ( invdt == 0 ) {
    for (int ig = 0; ig < numAtoms; ++ig )
      atom[ig].position = posNew[ig];
  } else if ( virial == 0 ) {
    for (int ig = 0; ig < numAtoms; ++ig )
      atom[ig].velocity = velNew[ig];
  } else {
    Tensor wc;  // constraint virial
    addRattleForce(invdt, wc);
    *virial += wc;
  }

  return 0;
}

//  RATTLE algorithm from Allen & Tildesley
int HomePatch::rattle1old(const BigReal timestep, Tensor *virial, 
    SubmitReduction *ppreduction)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;
  const BigReal dt = timestep / TIMEFACTOR;
  const BigReal invdt = (dt == 0.) ? 0. : 1.0 / dt; // precalc 1/dt
  BigReal tol2 = 2.0 * simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsq[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector pos[10];  // new position
  Vector vel[10];  // new velocity
  Vector netdp[10];  // total momentum change from constraint
  BigReal rmass[10];  // 1 / mass
  int fixed[10];  // is atom fixed?
  Tensor wc;  // constraint virial
  BigReal idz, zmin;
  int nslabs;

  // Initialize the settle algorithm with water parameters
  // settle1() assumes all waters are identical,
  // and will generate bad results if they are not.
  // XXX this will move to Molecule::build_atom_status when that 
  // version is debugged
  if ( ! settle_initialized ) {
    for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
      // find a water
      if (atom[ig].rigidBondLength > 0) {
        int oatm;
        if (simParams->watmodel == WAT_SWM4) {
          oatm = ig+3;  // skip over Drude and Lonepair
          //printf("ig=%d  mass_ig=%g  oatm=%d  mass_oatm=%g\n",
          //    ig, atom[ig].mass, oatm, atom[oatm].mass);
        }
        else {
          oatm = ig+1;
          // Avoid using the Om site to set this by mistake
          if (atom[ig].mass < 0.5 || atom[ig+1].mass < 0.5) {
            oatm += 1;
          }
        }

        // initialize settle water parameters
        settle1init(atom[ig].mass, atom[oatm].mass, 
                    atom[ig].rigidBondLength,
                    atom[oatm].rigidBondLength,
                    settle_mOrmT, settle_mHrmT, settle_ra,
                    settle_rb, settle_rc, settle_rra);
        settle_initialized = 1;
        break; // done with init
      }
    }
  }

  if (ppreduction) {
    nslabs = simParams->pressureProfileSlabs;
    idz = nslabs/lattice.c().z;
    zmin = lattice.origin().z - 0.5*lattice.c().z;
  }

  // Size of a hydrogen group for water
  int wathgsize = 3;
  int watmodel = simParams->watmodel;
  if (watmodel == WAT_TIP4) wathgsize = 4;
  else if (watmodel == WAT_SWM4) wathgsize = 5;
  
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    int anyfixed = 0;
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      pos[i] = atom[ig+i].position;
      vel[i] = atom[ig+i].velocity;
      rmass[i] = (atom[ig+i].mass > 0. ? 1. / atom[ig+i].mass : 0.);
      //printf("rmass of %i is %f\n", ig+i, rmass[i]);
      fixed[i] = ( fixedAtomsOn && atom[ig+i].atomFixed );
      //printf("fixed status of %i is %i\n", i, fixed[i]);
      // undo addVelocityToPosition to get proper reference coordinates
      if ( fixed[i] ) { anyfixed = 1; rmass[i] = 0.; } else pos[i] += vel[i] * dt;
    }
    int icnt = 0;
    if ( ( tmp = atom[ig].rigidBondLength ) > 0 ) {  // for water
      if (hgs != wathgsize) {
        char errmsg[256];
        sprintf(errmsg, "Water molecule starting with atom %d contains %d atoms "
                         "but the specified water model requires %d atoms.\n",
                         atom[ig].id+1, hgs, wathgsize);
        NAMD_die(errmsg);
      }
      // Use SETTLE for water unless some of the water atoms are fixed,
      if (useSettle && !anyfixed) {
        if (simParams->watmodel == WAT_SWM4) {
          // SWM4 ordering:  O D LP H1 H2
          // do swap(O,LP) and call settle with subarray O H1 H2
          // swap back after we return
          Vector lp_ref = ref[2];
          Vector lp_pos = pos[2];
          Vector lp_vel = vel[2];
          ref[2] = ref[0];
          pos[2] = pos[0];
          vel[2] = vel[0];
          settle1(ref+2, pos+2, vel+2, invdt,
                  settle_mOrmT, settle_mHrmT, settle_ra,
                  settle_rb, settle_rc, settle_rra);
          ref[0] = ref[2];
          pos[0] = pos[2];
          vel[0] = vel[2];
          ref[2] = lp_ref;
          pos[2] = lp_pos;
          vel[2] = lp_vel;
          // determine for LP updated pos and vel
          swm4_omrepos(ref, pos, vel, invdt);
        }
        else {
            settle1(ref, pos, vel, invdt,
                    settle_mOrmT, settle_mHrmT, settle_ra,
                    settle_rb, settle_rc, settle_rra);            
          if (simParams->watmodel == WAT_TIP4) {
            tip4_omrepos(ref, pos, vel, invdt);
          }
        }

        // which slab the hydrogen group will belong to
        // for pprofile calculations.
        int ppoffset, partition;
        if ( invdt == 0 ) for ( i = 0; i < wathgsize; ++i ) {
          atom[ig+i].position = pos[i];
        } else if ( virial == 0 ) for ( i = 0; i < wathgsize; ++i ) {
          atom[ig+i].velocity = vel[i];
        } else for ( i = 0; i < wathgsize; ++i ) {
          Force df = (vel[i] - atom[ig+i].velocity) * ( atom[ig+i].mass * invdt );
          Tensor vir = outer(df, ref[i]);
          wc += vir;
          f[Results::normal][ig+i] += df;
          atom[ig+i].velocity = vel[i];
          if (ppreduction) {
            // put all the atoms from a water in the same slab.  Atom 0
            // should be the parent atom.
            if (!i) {
              BigReal z = pos[i].z;
              partition = atom[ig].partition;
              int slab = (int)floor((z-zmin)*idz);
              if (slab < 0) slab += nslabs;
              else if (slab >= nslabs) slab -= nslabs;
              ppoffset = 3*(slab + nslabs*partition);
            }
            ppreduction->item(ppoffset  ) += vir.xx;
            ppreduction->item(ppoffset+1) += vir.yy;
            ppreduction->item(ppoffset+2) += vir.zz;
          }
        }
        continue;
      }
      if ( !(fixed[1] && fixed[2]) ) {
        dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = atom[ig+i].rigidBondLength ) > 0 ) {
        if ( !(fixed[0] && fixed[i]) ) {
          dsq[icnt] = tmp * tmp;  ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
        }
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    for ( i = 0; i < hgs; ++i ) {
      netdp[i] = 0.;
    }
    int done;
    int consFailure;
    for ( iter = 0; iter < maxiter; ++iter ) {
//if (iter > 0) CkPrintf("iteration %d\n", iter);
      done = 1;
      consFailure = 0;
      for ( i = 0; i < icnt; ++i ) {
        int a = ial[i];  int b = ibl[i];
        Vector pab = pos[a] - pos[b];
	      BigReal pabsq = pab.x*pab.x + pab.y*pab.y + pab.z*pab.z;
        BigReal rabsq = dsq[i];
        BigReal diffsq = rabsq - pabsq;
        if ( fabs(diffsq) > (rabsq * tol2) ) {
          Vector &rab = refab[i];
          BigReal rpab = rab.x*pab.x + rab.y*pab.y + rab.z*pab.z;
          if ( rpab < ( rabsq * 1.0e-6 ) ) {
            done = 0;
            consFailure = 1;
            continue;
          }
          BigReal rma = rmass[a];
          BigReal rmb = rmass[b];
          BigReal gab = diffsq / ( 2.0 * ( rma + rmb ) * rpab );
          Vector dp = rab * gab;
          pos[a] += rma * dp;
          pos[b] -= rmb * dp;
          if ( invdt != 0. ) {
            dp *= invdt;
            netdp[a] += dp;
            netdp[b] -= dp;
          }
          done = 0;
        }
      }
      if ( done ) break;
    }

    if ( consFailure ) {
      if ( dieOnError ) {
        iout << iERROR << "Constraint failure in RATTLE algorithm for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
        return -1;  // triggers early exit
      } else {
        iout << iWARN << "Constraint failure in RATTLE algorithm for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
      }
    } else if ( ! done ) {
      if ( dieOnError ) {
        iout << iERROR << "Exceeded RATTLE iteration limit for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
        return -1;  // triggers early exit
      } else {
        iout << iWARN << "Exceeded RATTLE iteration limit for atom "
        << (atom[ig].id + 1) << "!\n" << endi;
      }
    }

    // store data back to patch
    int ppoffset, partition;
    if ( invdt == 0 ) for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].position = pos[i];
    } else if ( virial == 0 ) for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].velocity = vel[i] + rmass[i] * netdp[i];
    } else for ( i = 0; i < hgs; ++i ) {
      Force df = netdp[i] * invdt;
      Tensor vir = outer(df, ref[i]);
      wc += vir;
      f[Results::normal][ig+i] += df;
      atom[ig+i].velocity = vel[i] + rmass[i] * netdp[i];
      if (ppreduction) {
        if (!i) {
          BigReal z = pos[i].z;
          int partition = atom[ig].partition;
          int slab = (int)floor((z-zmin)*idz);
          if (slab < 0) slab += nslabs;
          else if (slab >= nslabs) slab -= nslabs;
          ppoffset = 3*(slab + nslabs*partition);
        }
        ppreduction->item(ppoffset  ) += vir.xx;
        ppreduction->item(ppoffset+1) += vir.yy;
        ppreduction->item(ppoffset+2) += vir.zz;
      }
    }
  }
  if ( dt && virial ) *virial += wc;

  return 0;
}

//  RATTLE algorithm from Allen & Tildesley
void HomePatch::rattle2(const BigReal timestep, Tensor *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;
  const BigReal dt = timestep / TIMEFACTOR;
  Tensor wc;  // constraint virial
  BigReal tol = simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsqi[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector vel[10];  // new velocity
  BigReal rmass[10];  // 1 / mass
  BigReal redmass[10];  // reduced mass
  int fixed[10];  // is atom fixed?
  
  // Size of a hydrogen group for water
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  //  CkPrintf("In rattle2!\n");
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    //    CkPrintf("ig=%d\n",ig);
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    int anyfixed = 0;
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      vel[i] = atom[ig+i].velocity;
      rmass[i] = atom[ig+i].mass > 0. ? 1. / atom[ig+i].mass : 0.;
      fixed[i] = ( fixedAtomsOn && atom[ig+i].atomFixed );
      if ( fixed[i] ) { anyfixed = 1; rmass[i] = 0.; }
    }
    int icnt = 0;
    if ( ( tmp = atom[ig].rigidBondLength ) > 0 ) {  // for water
      if (hgs != wathgsize) {
        NAMD_bug("Hydrogen group error caught in rattle2().");
      }
      // Use SETTLE for water unless some of the water atoms are fixed,
      if (useSettle && !anyfixed) {
        if (simParams->watmodel == WAT_SWM4) {
          // SWM4 ordering:  O D LP H1 H2
          // do swap(O,LP) and call settle with subarray O H1 H2
          // swap back after we return
          Vector lp_ref = ref[2];
          // Vector lp_vel = vel[2];
          ref[2] = ref[0];
          vel[2] = vel[0];
          settle2(atom[ig].mass, atom[ig+3].mass, ref+2, vel+2, dt, virial);
          ref[0] = ref[2];
          vel[0] = vel[2];
          ref[2] = lp_ref;
          // vel[2] = vel[0];  // set LP vel to O vel
        }
        else {
          settle2(atom[ig].mass, atom[ig+1].mass, ref, vel, dt, virial);
          if (simParams->watmodel == WAT_TIP4) {
            vel[3] = vel[0];
          }
        }
        for (i=0; i<hgs; i++) {
          atom[ig+i].velocity = vel[i];
        }
        continue;
      }
      if ( !(fixed[1] && fixed[2]) ) {
	redmass[icnt] = 1. / (rmass[1] + rmass[2]);
	dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    //    CkPrintf("Loop 2\n");
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = atom[ig+i].rigidBondLength ) > 0 ) {
        if ( !(fixed[0] && fixed[i]) ) {
	  redmass[icnt] = 1. / (rmass[0] + rmass[i]);
	  dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 0;
	  ibl[icnt] = i;  ++icnt;
	}
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    //    CkPrintf("Loop 3\n");
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    //    CkPrintf("Loop 4\n");
    int done;
    for ( iter = 0; iter < maxiter; ++iter ) {
      done = 1;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector vab = vel[a] - vel[b];
	Vector &rab = refab[i];
	BigReal rabsqi = dsqi[i];
	BigReal rvab = rab.x*vab.x + rab.y*vab.y + rab.z*vab.z;
	if ( (fabs(rvab) * dt * rabsqi) > tol ) {
	  Vector dp = rab * (-rvab * redmass[i] * rabsqi);
	  wc += outer(dp,rab);
	  vel[a] += rmass[a] * dp;
	  vel[b] -= rmass[b] * dp;
	  done = 0;
	}
      }
      if ( done ) break;
      //if (done) { if (iter > 0) CkPrintf("iter=%d\n", iter); break; }
    }
    if ( ! done ) {
      if ( dieOnError ) {
	NAMD_die("Exceeded maximum number of iterations in rattle2().");
      } else {
	iout << iWARN <<
	  "Exceeded maximum number of iterations in rattle2().\n" << endi;
      }
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].velocity = vel[i];
    }
  }
  //  CkPrintf("Leaving rattle2!\n");
  // check that there isn't a constant needed here!
  *virial += wc / ( 0.5 * dt );

}


//  Adjust gradients for minimizer
void HomePatch::minimize_rattle2(const BigReal timestep, Tensor *virial, bool forces)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Force *f1 = f[Results::normal].begin();
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;
  const BigReal dt = timestep / TIMEFACTOR;
  Tensor wc;  // constraint virial
  BigReal tol = simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsqi[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector vel[10];  // new velocity
  BigReal rmass[10];  // 1 / mass
  BigReal redmass[10];  // reduced mass
  int fixed[10];  // is atom fixed?
  
  // Size of a hydrogen group for water
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  //  CkPrintf("In rattle2!\n");
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    //    CkPrintf("ig=%d\n",ig);
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    int anyfixed = 0;
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      vel[i] = ( forces ? f1[ig+i] : atom[ig+i].velocity );
      rmass[i] = 1.0;
      fixed[i] = ( fixedAtomsOn && atom[ig+i].atomFixed );
      if ( fixed[i] ) { anyfixed = 1; rmass[i] = 0.; }
    }
    int icnt = 0;
    if ( ( tmp = atom[ig].rigidBondLength ) > 0 ) {  // for water
      if (hgs != wathgsize) {
        NAMD_bug("Hydrogen group error caught in rattle2().");
      }
      // Use SETTLE for water unless some of the water atoms are fixed,
      if (useSettle && !anyfixed) {
        if (simParams->watmodel == WAT_SWM4) {
          // SWM4 ordering:  O D LP H1 H2
          // do swap(O,LP) and call settle with subarray O H1 H2
          // swap back after we return
          Vector lp_ref = ref[2];
          // Vector lp_vel = vel[2];
          ref[2] = ref[0];
          vel[2] = vel[0];
          settle2(1.0, 1.0, ref+2, vel+2, dt, virial);
          ref[0] = ref[2];
          vel[0] = vel[2];
          ref[2] = lp_ref;
          // vel[2] = vel[0];  // set LP vel to O vel
        }
        else {
          settle2(1.0, 1.0, ref, vel, dt, virial);
          if (simParams->watmodel == WAT_TIP4) {
            vel[3] = vel[0];
          }
        }
        for (i=0; i<hgs; i++) {
          ( forces ? f1[ig+i] : atom[ig+i].velocity ) = vel[i];
        }
        continue;
      }
      if ( !(fixed[1] && fixed[2]) ) {
	redmass[icnt] = 1. / (rmass[1] + rmass[2]);
	dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    //    CkPrintf("Loop 2\n");
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = atom[ig+i].rigidBondLength ) > 0 ) {
        if ( !(fixed[0] && fixed[i]) ) {
	  redmass[icnt] = 1. / (rmass[0] + rmass[i]);
	  dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 0;
	  ibl[icnt] = i;  ++icnt;
	}
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    //    CkPrintf("Loop 3\n");
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    //    CkPrintf("Loop 4\n");
    int done;
    for ( iter = 0; iter < maxiter; ++iter ) {
      done = 1;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector vab = vel[a] - vel[b];
	Vector &rab = refab[i];
	BigReal rabsqi = dsqi[i];
	BigReal rvab = rab.x*vab.x + rab.y*vab.y + rab.z*vab.z;
	if ( (fabs(rvab) * dt * rabsqi) > tol ) {
	  Vector dp = rab * (-rvab * redmass[i] * rabsqi);
	  wc += outer(dp,rab);
	  vel[a] += rmass[a] * dp;
	  vel[b] -= rmass[b] * dp;
	  done = 0;
	}
      }
      if ( done ) break;
      //if (done) { if (iter > 0) CkPrintf("iter=%d\n", iter); break; }
    }
    if ( ! done ) {
      if ( dieOnError ) {
	NAMD_die("Exceeded maximum number of iterations in rattle2().");
      } else {
	iout << iWARN <<
	  "Exceeded maximum number of iterations in rattle2().\n" << endi;
      }
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      ( forces ? f1[ig+i] : atom[ig+i].velocity ) = vel[i];
    }
  }
  //  CkPrintf("Leaving rattle2!\n");
  // check that there isn't a constant needed here!
  *virial += wc / ( 0.5 * dt );

}


// BEGIN LA
void HomePatch::loweAndersenVelocities()
{
    DebugM(2, "loweAndersenVelocities\n");
    Molecule *mol = Node::Object()->molecule;
    SimParameters *simParams = Node::Object()->simParameters;
    v.resize(numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
	//v[i] = p[i];
	// co-opt CompAtom structure to pass velocity and mass information
	v[i].position = atom[i].velocity;
	v[i].charge = atom[i].mass;
    }
    DebugM(2, "loweAndersenVelocities\n");
}

void HomePatch::loweAndersenFinish()
{
    DebugM(2, "loweAndersenFinish\n");
    v.resize(0);
}
// END LA

//LCPO
void HomePatch::setLcpoType() {
  Molecule *mol = Node::Object()->molecule;
  const int *lcpoParamType = mol->getLcpoParamType();

  lcpoType.resize(numAtoms);
  for (int i = 0; i < numAtoms; i++) {
    lcpoType[i] = lcpoParamType[pExt[i].id];
  }
}

//set intrinsic radii of atom when doMigration
void HomePatch::setGBISIntrinsicRadii() {
  intRad.resize(numAtoms*2);
  intRad.setall(0);
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Real offset = simParams->coulomb_radius_offset;
  for (int i = 0; i < numAtoms; i++) {
    Real rad = MassToRadius(atom[i].mass);//in ComputeGBIS.inl
    Real screen = MassToScreen(atom[i].mass);//same
    intRad[2*i+0] = rad - offset;//r0
    intRad[2*i+1] = screen*(rad - offset);//s0
  }
}

//compute born radius after phase 1, before phase 2
void HomePatch::gbisComputeAfterP1() {

  SimParameters *simParams = Node::Object()->simParameters;
  BigReal alphaMax = simParams->alpha_max;
  BigReal delta = simParams->gbis_delta;
  BigReal beta = simParams->gbis_beta;
  BigReal gamma = simParams->gbis_gamma;
  BigReal coulomb_radius_offset = simParams->coulomb_radius_offset;

  BigReal rhoi;
  BigReal rhoi0;
  //calculate bornRad from psiSum
  for (int i = 0; i < numAtoms; i++) {
    rhoi0 = intRad[2*i];
    rhoi = rhoi0+coulomb_radius_offset;
    psiFin[i] += psiSum[i];
    psiFin[i] *= rhoi0;
    bornRad[i]=1/(1/rhoi0-1/rhoi*tanh(psiFin[i]*(delta+psiFin[i]*(-beta+gamma*psiFin[i]))));
    bornRad[i] = (bornRad[i] > alphaMax) ? alphaMax : bornRad[i];
#ifdef PRINT_COMP
    CkPrintf("BORNRAD(%04i)[%04i] = % .4e\n",flags.sequence,pExt[i].id,bornRad[i]);
#endif
  }

  gbisP2Ready();
}

//compute dHdrPrefix after phase 2, before phase 3
void HomePatch::gbisComputeAfterP2() {

  SimParameters *simParams = Node::Object()->simParameters;
  BigReal delta = simParams->gbis_delta;
  BigReal beta = simParams->gbis_beta;
  BigReal gamma = simParams->gbis_gamma;
  BigReal epsilon_s = simParams->solvent_dielectric;
  BigReal epsilon_p = simParams->dielectric;
  BigReal epsilon_s_i = 1/simParams->solvent_dielectric;
  BigReal epsilon_p_i = 1/simParams->dielectric;
  BigReal coulomb_radius_offset = simParams->coulomb_radius_offset;
  BigReal kappa = simParams->kappa;
  BigReal fij, expkappa, Dij, dEdai, dedasum;
  BigReal rhoi, rhoi0, psii, nbetapsi;
  BigReal gammapsi2, tanhi, daidr;
  for (int i = 0; i < numAtoms; i++) {
    //add diagonal dEda term
    dHdrPrefix[i] += dEdaSum[i];//accumulated from proxies
    fij = bornRad[i];//inf
    expkappa = exp(-kappa*fij);//0
    Dij = epsilon_p_i - expkappa*epsilon_s_i;//dielectric term
    //calculate dHij prefix
    dEdai = -0.5*COULOMB*atom[i].charge*atom[i].charge
                  *(kappa*epsilon_s_i*expkappa-Dij/fij)/bornRad[i];
    dHdrPrefix[i] += dEdai;
    dedasum = dHdrPrefix[i];

    rhoi0 = intRad[2*i];
    rhoi = rhoi0+coulomb_radius_offset;
    psii = psiFin[i];
    nbetapsi = -beta*psii;
    gammapsi2 = gamma*psii*psii;
    tanhi = tanh(psii*(delta+nbetapsi+gammapsi2));
    daidr = bornRad[i]*bornRad[i]*rhoi0/rhoi*(1-tanhi*tanhi)
           * (delta+nbetapsi+nbetapsi+gammapsi2+gammapsi2+gammapsi2);
    dHdrPrefix[i] *= daidr;//dHdrPrefix previously equaled dEda
#ifdef PRINT_COMP
    CkPrintf("DHDR(%04i)[%04i] = % .4e\n",flags.sequence,pExt[i].id,dHdrPrefix[i]);
#endif
  }
  gbisP3Ready();
}

//send born radius to proxies to begin phase 2
void HomePatch::gbisP2Ready() {
  if (proxy.size() > 0) {
    CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
    for (int i = 0; i < proxy.size(); i++) {
      int node = proxy[i];
      ProxyGBISP2DataMsg *msg=new(numAtoms,PRIORITY_SIZE) ProxyGBISP2DataMsg;
      msg->patch = patchID;
      msg->origPe = CkMyPe();
      memcpy(msg->bornRad,bornRad.begin(),numAtoms*sizeof(Real));
      msg->destPe = node;
      int seq = flags.sequence;
      int priority = GB2_PROXY_DATA_PRIORITY + PATCH_PRIORITY(patchID);
      SET_PRIORITY(msg,seq,priority);
      cp[node].recvData(msg);
    }
  }
  Patch::gbisP2Ready();
}

//send dHdrPrefix to proxies to begin phase 3
void HomePatch::gbisP3Ready() {
  if (proxy.size() > 0) {
    CProxy_ProxyMgr cp(CkpvAccess(BOCclass_group).proxyMgr);
    //only nonzero message should be sent for doFullElec
    int msgAtoms = (flags.doFullElectrostatics) ? numAtoms : 0;
    for (int i = 0; i < proxy.size(); i++) {
      int node = proxy[i];
      ProxyGBISP3DataMsg *msg = new(msgAtoms,PRIORITY_SIZE) ProxyGBISP3DataMsg;
      msg->patch = patchID;
      msg->dHdrPrefixLen = msgAtoms;
      msg->origPe = CkMyPe();
      memcpy(msg->dHdrPrefix, dHdrPrefix.begin(), msgAtoms*sizeof(Real));
      msg->destPe = node;
      int seq = flags.sequence;
      int priority = GB3_PROXY_DATA_PRIORITY + PATCH_PRIORITY(patchID);
      SET_PRIORITY(msg,seq,priority);
      cp[node].recvData(msg);
    }
  }
  Patch::gbisP3Ready();
}

//receive proxy results from phase 1
void HomePatch::receiveResult(ProxyGBISP1ResultMsg *msg) {
  ++numGBISP1Arrived;
    for ( int i = 0; i < msg->psiSumLen; ++i ) {
      psiFin[i] += msg->psiSum[i];
    }
  delete msg;

  if (flags.doNonbonded) {
    //awaken if phase 1 done
    if (phase1BoxClosedCalled == true &&
        numGBISP1Arrived==proxy.size() ) {
        sequencer->awaken();
    }
  } else {
    //awaken if all phases done on noWork step
    if (boxesOpen == 0 &&
        numGBISP1Arrived == proxy.size() &&
        numGBISP2Arrived == proxy.size() &&
        numGBISP3Arrived == proxy.size()) {
      sequencer->awaken();
    }
  }
}

//receive proxy results from phase 2
void HomePatch::receiveResult(ProxyGBISP2ResultMsg *msg) {
  ++numGBISP2Arrived;
  //accumulate dEda
  for ( int i = 0; i < msg->dEdaSumLen; ++i ) {
    dHdrPrefix[i] += msg->dEdaSum[i];
  }
  delete msg;

  if (flags.doNonbonded) {
    //awaken if phase 2 done
    if (phase2BoxClosedCalled == true &&
        numGBISP2Arrived==proxy.size() ) {
      sequencer->awaken();
    }
  } else {
    //awaken if all phases done on noWork step
    if (boxesOpen == 0 &&
        numGBISP1Arrived == proxy.size() &&
        numGBISP2Arrived == proxy.size() &&
        numGBISP3Arrived == proxy.size()) {
      sequencer->awaken();
    }
  }
}

//  MOLLY algorithm part 1
void HomePatch::mollyAverage()
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  BigReal tol = simParams->mollyTol;
  int maxiter = simParams->mollyIter;
  int i, iter;
  HGArrayBigReal dsq;
  BigReal tmp;
  HGArrayInt ial, ibl;
  HGArrayVector ref;  // reference position
  HGArrayVector refab;  // reference vector
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  CompAtom *avg;  // averaged position
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  //  iout<<iINFO << "mollyAverage: "<<std::endl<<endi;
  p_avg.resize(numAtoms);
  for ( i=0; i<numAtoms; ++i ) p_avg[i] = p[i];

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = atom[ig+i].position;
	  rmass[i] = 1. / atom[ig+i].mass;
	  fixed[i] = ( simParams->fixedAtomsOn && atom[ig+i].atomFixed );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	avg = &(p_avg[ig]);
	int icnt = 0;

  if ( ( tmp = atom[ig].rigidBondLength ) > 0 ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyAverage().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
    if ( ( tmp = atom[ig+i].rigidBondLength ) > 0 ) {
	    if ( !(fixed[0] && fixed[i]) ) {
	      dsq[icnt] =  tmp * tmp;  ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	    }
	  }
	}
	if ( icnt == 0 ) continue;  // no constraints
	numLambdas += icnt;
	molly_lambda.resize(numLambdas);
	lambda = &(molly_lambda[numLambdas - icnt]);
	for ( i = 0; i < icnt; ++i ) {
	  refab[i] = ref[ial[i]] - ref[ibl[i]];
	}
	//	iout<<iINFO<<"hgs="<<hgs<<" m="<<icnt<<std::endl<<endi;
	iter=average(avg,ref,lambda,hgs,icnt,rmass,dsq,ial,ibl,refab,tol,maxiter);
	if ( iter == maxiter ) {
	  iout << iWARN << "Exceeded maximum number of iterations in mollyAverage().\n"<<endi;
	}
  }

  // for ( i=0; i<numAtoms; ++i ) {
  //    if ( ( p_avg[i].position - p[i].position ).length2() > 1.0 ) {
  //      iout << iERROR << "MOLLY moved atom " << (p[i].id + 1) << " from "
  //        << p[i].position << " to " << p_avg[i].position << "\n" << endi;
  //    }
  // }

}


//  MOLLY algorithm part 2
void HomePatch::mollyMollify(Tensor *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Tensor wc;  // constraint virial
  int i;
  HGArrayInt ial, ibl;
  HGArrayVector ref;  // reference position
  CompAtom *avg;  // averaged position
  HGArrayVector refab;  // reference vector
  HGArrayVector force;  // new force
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if (hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = atom[ig+i].position;
	  force[i] = f[Results::slow][ig+i];
	  rmass[i] = 1. / atom[ig+i].mass;
	  fixed[i] = ( simParams->fixedAtomsOn && atom[ig+i].atomFixed );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	int icnt = 0;
	// c-ji I'm only going to mollify water for now
  if ( atom[ig].rigidBondLength > 0 ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyMollify().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
    if ( atom[ig+i].rigidBondLength > 0 ) {
	    if ( !(fixed[0] && fixed[i]) ) {
	      ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	    }
	  }
	}

	if ( icnt == 0 ) continue;  // no constraints
	lambda = &(molly_lambda[numLambdas]);
	numLambdas += icnt;
	for ( i = 0; i < icnt; ++i ) {
	  refab[i] = ref[ial[i]] - ref[ibl[i]];
	}
	avg = &(p_avg[ig]);
	mollify(avg,ref,lambda,force,hgs,icnt,rmass,ial,ibl,refab);
	// store data back to patch
	for ( i = 0; i < hgs; ++i ) {
	  wc += outer(force[i]-f[Results::slow][ig+i],ref[i]);
	  f[Results::slow][ig+i] = force[i];
	}
  }
  // check that there isn't a constant needed here!
  *virial += wc;
  p_avg.resize(0);
}

void HomePatch::checkpoint(void) {
  checkpoint_atom.copy(atom);
  checkpoint_lattice = lattice;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    checkpoint_numWaterAtoms = numWaterAtoms;
  #endif
}

void HomePatch::revert(void) {
  atomMapper->unregisterIDsFullAtom(atom.begin(),atom.end());

  atom.copy(checkpoint_atom);
  numAtoms = atom.size();
  lattice = checkpoint_lattice;

  doAtomUpdate = true;
  rattleListValid = false;

  if ( ! numNeighbors ) atomMapper->registerIDsFullAtom(atom.begin(),atom.end());

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = checkpoint_numWaterAtoms;
  #endif
}

void HomePatch::exchangeCheckpoint(int scriptTask, int &bpc) {  // initiating replica
  SimParameters *simParams = Node::Object()->simParameters;
  checkpoint_task = scriptTask;
  const int remote = simParams->scriptIntArg1;
  const char *key = simParams->scriptStringArg1;
  PatchMgr::Object()->sendCheckpointReq(patchID, remote, key, scriptTask);
}

void HomePatch::recvCheckpointReq(int task, const char *key, int replica, int pe) {  // responding replica
  if ( task == SCRIPT_CHECKPOINT_FREE ) {
    if ( ! checkpoints.count(key) ) {
      NAMD_die("Unable to free checkpoint, requested key was never stored.");
    }
    delete checkpoints[key];
    checkpoints.erase(key);
    PatchMgr::Object()->sendCheckpointAck(patchID, replica, pe);
    return;
  }
  CheckpointAtomsMsg *msg;
  if ( task == SCRIPT_CHECKPOINT_LOAD || task == SCRIPT_CHECKPOINT_SWAP ) {
    if ( ! checkpoints.count(key) ) {
      NAMD_die("Unable to load checkpoint, requested key was never stored.");
    }
    checkpoint_t &cp = *checkpoints[key];
    msg = new (cp.numAtoms,1,0) CheckpointAtomsMsg;
    msg->lattice = cp.lattice;
    msg->berendsenPressure_count = cp.berendsenPressure_count;
    msg->numAtoms = cp.numAtoms;
    memcpy(msg->atoms,cp.atoms.begin(),cp.numAtoms*sizeof(FullAtom));
  } else {
    msg = new (0,1,0) CheckpointAtomsMsg;
  }
  msg->pid = patchID;
  msg->replica = CmiMyPartition();
  msg->pe = CkMyPe();
  PatchMgr::Object()->sendCheckpointLoad(msg, replica, pe);
}

void HomePatch::recvCheckpointLoad(CheckpointAtomsMsg *msg) {  // initiating replica
  SimParameters *simParams = Node::Object()->simParameters;
  const int remote = simParams->scriptIntArg1;
  const char *key = simParams->scriptStringArg1;
  if ( checkpoint_task == SCRIPT_CHECKPOINT_FREE ) {
    NAMD_bug("HomePatch::recvCheckpointLoad called during checkpointFree.");
  }
  if ( msg->replica != remote ) {
    NAMD_bug("HomePatch::recvCheckpointLoad message from wrong replica.");
  }
  if ( checkpoint_task == SCRIPT_CHECKPOINT_STORE || checkpoint_task == SCRIPT_CHECKPOINT_SWAP ) {
    CheckpointAtomsMsg *newmsg = new (numAtoms,1+strlen(key),0) CheckpointAtomsMsg;
    strcpy(newmsg->key,key);
    newmsg->lattice = lattice;
    newmsg->berendsenPressure_count = sequencer->berendsenPressure_count;
    newmsg->pid = patchID;
    newmsg->pe = CkMyPe();
    newmsg->replica = CmiMyPartition();
    newmsg->numAtoms = numAtoms;
    memcpy(newmsg->atoms,atom.begin(),numAtoms*sizeof(FullAtom));
    PatchMgr::Object()->sendCheckpointStore(newmsg, remote, msg->pe);
  }
  if ( checkpoint_task == SCRIPT_CHECKPOINT_LOAD || checkpoint_task == SCRIPT_CHECKPOINT_SWAP ) {
    atomMapper->unregisterIDsFullAtom(atom.begin(),atom.end());
    lattice = msg->lattice;
    sequencer->berendsenPressure_count = msg->berendsenPressure_count;
    numAtoms = msg->numAtoms;
    atom.resize(numAtoms);
    memcpy(atom.begin(),msg->atoms,numAtoms*sizeof(FullAtom));
    doAtomUpdate = true;
    rattleListValid = false;
    if ( ! numNeighbors ) atomMapper->registerIDsFullAtom(atom.begin(),atom.end());
  }
  if ( checkpoint_task == SCRIPT_CHECKPOINT_LOAD ) {
    recvCheckpointAck();
  }
  delete msg;
}

void HomePatch::recvCheckpointStore(CheckpointAtomsMsg *msg) {  // responding replica
  if ( ! checkpoints.count(msg->key) ) {
    checkpoints[msg->key] = new checkpoint_t;
  }
  checkpoint_t &cp = *checkpoints[msg->key];
  cp.lattice = msg->lattice;
  cp.berendsenPressure_count = msg->berendsenPressure_count;
  cp.numAtoms = msg->numAtoms;
  cp.atoms.resize(cp.numAtoms);
  memcpy(cp.atoms.begin(),msg->atoms,cp.numAtoms*sizeof(FullAtom));
  PatchMgr::Object()->sendCheckpointAck(patchID, msg->replica, msg->pe);
  delete msg;
}

void HomePatch::recvCheckpointAck() {  // initiating replica
  CkpvAccess(_qd)->process();
}


void HomePatch::exchangeAtoms(int scriptTask) {
  SimParameters *simParams = Node::Object()->simParameters;
  // CkPrintf("exchangeAtoms %d %d %d %d\n", CmiMyPartition(), scriptTask, (int)(simParams->scriptArg1), (int)(simParams->scriptArg2));
  if ( scriptTask == SCRIPT_ATOMSEND || scriptTask == SCRIPT_ATOMSENDRECV ) {
    exchange_dst = (int) simParams->scriptArg1;
    // create and save outgoing message
    exchange_msg = new (numAtoms,0) ExchangeAtomsMsg;
    exchange_msg->lattice = lattice;
    exchange_msg->pid = patchID;
    exchange_msg->numAtoms = numAtoms;
    memcpy(exchange_msg->atoms,atom.begin(),numAtoms*sizeof(FullAtom));
    if ( exchange_req >= 0 ) {
      recvExchangeReq(exchange_req);
    }
  }
  if ( scriptTask == SCRIPT_ATOMRECV || scriptTask == SCRIPT_ATOMSENDRECV ) {
    exchange_src = (int) simParams->scriptArg2;
    PatchMgr::Object()->sendExchangeReq(patchID, exchange_src);
  }
}

void HomePatch::recvExchangeReq(int req) {
  exchange_req = req;
  if ( exchange_msg ) {
    // CkPrintf("recvExchangeReq %d %d\n", CmiMyPartition(), exchange_dst);
    PatchMgr::Object()->sendExchangeMsg(exchange_msg, exchange_dst, exchange_req);
    exchange_msg = 0;
    exchange_req = -1;
    CkpvAccess(_qd)->process();
  }
}

void HomePatch::recvExchangeMsg(ExchangeAtomsMsg *msg) {
  // CkPrintf("recvExchangeMsg %d %d\n", CmiMyPartition(), exchange_src);
  atomMapper->unregisterIDsFullAtom(atom.begin(),atom.end());
  lattice = msg->lattice;
  numAtoms = msg->numAtoms;
  atom.resize(numAtoms);
  memcpy(atom.begin(),msg->atoms,numAtoms*sizeof(FullAtom));
  delete msg;
  CkpvAccess(_qd)->process();
  doAtomUpdate = true;
  rattleListValid = false;
  if ( ! numNeighbors ) atomMapper->registerIDsFullAtom(atom.begin(),atom.end());
}

void HomePatch::submitLoadStats(int timestep)
{
  LdbCoordinator::Object()->patchLoad(patchID,numAtoms,timestep);
}


void HomePatch::doPairlistCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;

  if ( numAtoms == 0 || ! flags.usePairlists ) {
    flags.pairlistTolerance = 0.;
    flags.maxAtomMovement = 99999.;
    return;
  }

  int i; int n = numAtoms;
  CompAtom *p_i = p.begin();

  if ( flags.savePairlists ) {
    flags.pairlistTolerance = doPairlistCheck_newTolerance;
    flags.maxAtomMovement = 0.;
    doPairlistCheck_newTolerance *= (1. - simParams->pairlistShrink);
    doPairlistCheck_lattice = lattice;
    doPairlistCheck_positions.resize(numAtoms);
    CompAtom *psave_i = doPairlistCheck_positions.begin();
    for ( i=0; i<n; ++i ) { psave_i[i] = p_i[i]; }
    return;
  }

  Lattice &lattice_old = doPairlistCheck_lattice;
  Position center_cur = lattice.unscale(center);
  Position center_old = lattice_old.unscale(center);
  Vector center_delta = center_cur - center_old;
  
  // find max deviation to corner (any neighbor shares a corner)
  BigReal max_cd = 0.;
  for ( i=0; i<2; ++i ) {
    for ( int j=0; j<2; ++j ) {
      for ( int k=0; k<2; ++k ) {
	ScaledPosition corner(	i ? min.x : max.x ,
				j ? min.y : max.y ,
				k ? min.z : max.z );
	Vector corner_delta =
		lattice.unscale(corner) - lattice_old.unscale(corner);
        corner_delta -= center_delta;
	BigReal cd = corner_delta.length2();
        if ( cd > max_cd ) max_cd = cd;
      }
    }
  }
  max_cd = sqrt(max_cd);

  // find max deviation of atoms relative to center
  BigReal max_pd = 0.;
  CompAtom *p_old_i = doPairlistCheck_positions.begin();
  for ( i=0; i<n; ++i ) {
    Vector p_delta = p_i[i].position - p_old_i[i].position;
    p_delta -= center_delta;
    BigReal pd = p_delta.length2();
    if ( pd > max_pd ) max_pd = pd;
  }
  max_pd = sqrt(max_pd);

  BigReal max_tol = max_pd + max_cd;

  flags.maxAtomMovement = max_tol;

  // if ( max_tol > flags.pairlistTolerance ) iout << "tolerance " << max_tol << " > " << flags.pairlistTolerance << "\n" << endi;

  if ( max_tol > ( (1. - simParams->pairlistTrigger) *
				doPairlistCheck_newTolerance ) ) {
    doPairlistCheck_newTolerance *= (1. + simParams->pairlistGrow);
  }

  if ( max_tol > doPairlistCheck_newTolerance ) {
    doPairlistCheck_newTolerance = max_tol / (1. - simParams->pairlistTrigger);
  }

}

void HomePatch::doGroupSizeCheck()
{
  if ( ! flags.doNonbonded ) return;

  SimParameters *simParams = Node::Object()->simParameters;
  BigReal hgcut = 0.5 * simParams->hgroupCutoff;  hgcut *= hgcut;
  BigReal maxrad2 = 0.;

  FullAtomList::iterator p_i = atom.begin();
  FullAtomList::iterator p_e = atom.end();

  while ( p_i != p_e ) {
    const int hgs = p_i->hydrogenGroupSize;
    if ( ! hgs ) break;  // avoid infinite loop on bug
    int ngs = hgs;
    if ( ngs > 5 ) ngs = 5;  // limit to at most 5 atoms per group
    BigReal x = p_i->position.x;
    BigReal y = p_i->position.y;
    BigReal z = p_i->position.z;
    int i;
    for ( i = 1; i < ngs; ++i ) {  // limit spatial extent
      p_i[i].nonbondedGroupSize = 0;
      BigReal dx = p_i[i].position.x - x;
      BigReal dy = p_i[i].position.y - y;
      BigReal dz = p_i[i].position.z - z;
      BigReal r2 = dx * dx + dy * dy + dz * dz;
      if ( r2 > hgcut ) break;
      else if ( r2 > maxrad2 ) maxrad2 = r2;
    }
    p_i->nonbondedGroupSize = i;
    for ( ; i < hgs; ++i ) {
      p_i[i].nonbondedGroupSize = 1;
    }
    p_i += hgs;
  }

  if ( p_i != p_e ) {
    NAMD_bug("hydrogenGroupSize is zero in HomePatch::doGroupSizeCheck");
  }

  flags.maxGroupRadius = sqrt(maxrad2);

}

void HomePatch::doMarginCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;

  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal sysdimc = lattice.c_r().unit() * lattice.c();

  BigReal minSize = simParams->patchDimension - simParams->margin;

  if ( ( aAwayDist*sysdima < minSize*0.9999 ) ||
       ( bAwayDist*sysdimb < minSize*0.9999 ) ||
       ( cAwayDist*sysdimc < minSize*0.9999 ) ) {

    NAMD_die("Periodic cell has become too small for original patch grid!\n"
      "Possible solutions are to restart from a recent checkpoint,\n"
      "increase margin, or disable useFlexibleCell for liquid simulation.");
  }

  BigReal cutoff = simParams->cutoff;

  BigReal margina = 0.5 * ( aAwayDist - cutoff / sysdima );
  BigReal marginb = 0.5 * ( bAwayDist - cutoff / sysdimb );
  BigReal marginc = 0.5 * ( cAwayDist - cutoff / sysdimc );

  if ( (margina < -0.0001) || (marginb < -0.0001) || (marginc < -0.0001) ) {
    NAMD_die("Periodic cell has become too small for original patch grid!\n"
      "There are probably many margin violations already reported.\n"
      "Possible solutions are to restart from a recent checkpoint,\n"
      "increase margin, or disable useFlexibleCell for liquid simulation.");
  }

  BigReal minx = min.x - margina;
  BigReal miny = min.y - marginb;
  BigReal minz = min.z - marginc;
  BigReal maxx = max.x + margina;
  BigReal maxy = max.y + marginb;
  BigReal maxz = max.z + marginc;

  int xdev, ydev, zdev;
  int problemCount = 0;

  FullAtomList::iterator p_i = atom.begin();
  FullAtomList::iterator p_e = atom.end();
  for ( ; p_i != p_e; ++p_i ) {

    ScaledPosition s = lattice.scale(p_i->position);

    // check if atom is within bounds
    if (s.x < minx) xdev = 0;
    else if (maxx <= s.x) xdev = 2; 
    else xdev = 1;

    if (s.y < miny) ydev = 0;
    else if (maxy <= s.y) ydev = 2; 
    else ydev = 1;

    if (s.z < minz) zdev = 0;
    else if (maxz <= s.z) zdev = 2; 
    else zdev = 1;

    if (mInfo[xdev][ydev][zdev]) { // somewhere else to be
	++problemCount;
    }

  }

  marginViolations = problemCount;
  // if ( problemCount ) {
  //     iout << iERROR <<
  //       "Found " << problemCount << " margin violations!\n" << endi;
  // } 

}


void
HomePatch::doAtomMigration()
{
  int i;

  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList.resize(0);
  }

  // Purge the AtomMap
  atomMapper->unregisterIDsFullAtom(atom.begin(),atom.end());

  // Determine atoms that need to migrate

  BigReal minx = min.x;
  BigReal miny = min.y;
  BigReal minz = min.z;
  BigReal maxx = max.x;
  BigReal maxy = max.y;
  BigReal maxz = max.z;

  int xdev, ydev, zdev;
  int delnum = 0;

  FullAtomList::iterator atom_i = atom.begin();
  FullAtomList::iterator atom_e = atom.end();

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    FullAtomList::iterator atom_first = atom_i;
    int numLostWaterAtoms = 0;
  #endif

  while ( atom_i != atom_e ) {
    if ( atom_i->migrationGroupSize ) {
      Position pos = atom_i->position;
      if ( atom_i->migrationGroupSize != atom_i->hydrogenGroupSize ) {
        int mgs = atom_i->migrationGroupSize;
        int c = 1;
        for ( int j=atom_i->hydrogenGroupSize; j<mgs;
				j+=(atom_i+j)->hydrogenGroupSize ) {
          pos += (atom_i+j)->position;
          ++c;
        }
        pos *= 1./c;
        // iout << "mgroup " << atom_i->id << " at " << pos << "\n" << endi;
      }

      ScaledPosition s = lattice.scale(pos);

      // check if atom is within bounds
      if (s.x < minx) xdev = 0;
      else if (maxx <= s.x) xdev = 2;
      else xdev = 1;

      if (s.y < miny) ydev = 0;
      else if (maxy <= s.y) ydev = 2;
      else ydev = 1;

      if (s.z < minz) zdev = 0;
      else if (maxz <= s.z) zdev = 2;
      else zdev = 1;

    }

    if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                    // Don't migrate if destination is myself

      // See if we have a migration list already
      MigrationList &mCur = mInfo[xdev][ydev][zdev]->mList;
      DebugM(3,"Migrating atom " << atom_i->id << " from patch "
		<< patchID << " with position " << atom_i->position << "\n");
      mCur.add(*atom_i);

      ++delnum;


      // DMK - Atom Separation (water vs. non-water)
      #if NAMD_SeparateWaters != 0
        // Check to see if this atom is part of a water molecule.  If
        //   so, numWaterAtoms needs to be adjusted to refect the lost of
        //   this atom.
        // NOTE: The atom separation code assumes that if the oxygen
        //   atom of the hydrogen group making up the water molecule is
        //   migrated to another HomePatch, the hydrogens will also
        //   move!!!
        int atomIndex = atom_i - atom_first;
        if (atomIndex < numWaterAtoms)
          numLostWaterAtoms++;
      #endif


    } else {

      if ( delnum ) { *(atom_i-delnum) = *atom_i; }

    }

    ++atom_i;
  }

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms -= numLostWaterAtoms;
  #endif


  int delpos = numAtoms - delnum;
  DebugM(4,"numAtoms " << numAtoms << " deleted " << delnum << "\n");
  atom.del(delpos,delnum);

  numAtoms = atom.size();

  PatchMgr::Object()->sendMigrationMsgs(patchID, realInfo, numNeighbors);

  // signal depositMigration() that we are inMigration mode
  inMigration = true;

  // Drain the migration message buffer
  for (i=0; i<numMlBuf; i++) {
     DebugM(1, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
     depositMigration(msgbuf[i]);
  }
  numMlBuf = 0;
     
  if (!allMigrationIn) {
    DebugM(3,"All Migrations NOT in, we are suspending patch "<<patchID<<"\n");
    migrationSuspended = true;
    sequencer->suspend();
    migrationSuspended = false;
  }
  allMigrationIn = false;

  inMigration = false;
  marginViolations = 0;
}

void 
HomePatch::depositMigration(MigrateAtomsMsg *msg)
{

  if (!inMigration) { // We have to buffer changes due to migration
		      // until our patch is in migration mode
    msgbuf[numMlBuf++] = msg;
    return;
  } 


  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0


    // Merge the incoming list of atoms with the current list of
    //   atoms.  Note that mergeSeparatedAtomList() will apply any
    //   required transformations to the incoming atoms as it is
    //   separating them.
    mergeAtomList(msg->migrationList);


  #else

    // Merge the incoming list of atoms with the current list of
    // atoms.  Apply transformations to the incoming atoms as they are
    // added to this patch's list.
    {
      MigrationList &migrationList = msg->migrationList;
      MigrationList::iterator mi;
      Transform mother_transform;
      for (mi = migrationList.begin(); mi != migrationList.end(); mi++) {
        DebugM(1,"Migrating atom " << mi->id << " to patch "
		  << patchID << " with position " << mi->position << "\n"); 
        if ( mi->migrationGroupSize ) {
          if ( mi->migrationGroupSize != mi->hydrogenGroupSize ) {
            Position pos = mi->position;
            int mgs = mi->migrationGroupSize;
            int c = 1;
            for ( int j=mi->hydrogenGroupSize; j<mgs;
                                j+=(mi+j)->hydrogenGroupSize ) {
              pos += (mi+j)->position;
              ++c;
            }
            pos *= 1./c;
            // iout << "mgroup " << mi->id << " at " << pos << "\n" << endi;
            mother_transform = mi->transform;
            pos = lattice.nearest(pos,center,&mother_transform);
            mi->position = lattice.reverse_transform(mi->position,mi->transform);
            mi->position = lattice.apply_transform(mi->position,mother_transform);
            mi->transform = mother_transform;
          } else {
            mi->position = lattice.nearest(mi->position,center,&(mi->transform));
            mother_transform = mi->transform;
          }
        } else {
          mi->position = lattice.reverse_transform(mi->position,mi->transform);
          mi->position = lattice.apply_transform(mi->position,mother_transform);
          mi->transform = mother_transform;
        }
        atom.add(*mi);
      }
    }


  #endif // if (NAMD_SeparateWaters != 0)


  numAtoms = atom.size();
  delete msg;

  DebugM(3,"Counter on " << patchID << " = " << patchMigrationCounter << "\n");
  if (!--patchMigrationCounter) {
    DebugM(3,"All Migrations are in for patch "<<patchID<<"\n");
    allMigrationIn = true;
    patchMigrationCounter = numNeighbors;
    if (migrationSuspended) {
      DebugM(3,"patch "<<patchID<<" is being awakened\n");
      migrationSuspended = false;
      sequencer->awaken();
      return;
    }
    else {
       DebugM(3,"patch "<<patchID<<" was not suspended\n");
    }
  }
}



// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0

// This function will separate waters from non-waters in the current
//   atom list (regardless of whether or not the atom list is has been
//   sorted yet or not).
void HomePatch::separateAtoms() {

  // Basic Idea:  Iterate through all the atoms in the current list
  //   of atoms.  Pack the waters in the current atoms list and move
  //   the non-waters to the scratch list.  Once the atoms have all
  //   been separated, copy the non-waters to the end of the waters.
  // NOTE:  This code does not assume that the atoms list has been
  //   separated in any manner.

  // NOTE: Sanity check - Doesn't look like the default constructor actually
  //   adds any atoms but does set numAtoms. ???
  if (atom.size() < 0) return;  // Nothing to do.

  // Resize the scratch FullAtomList (tempAtom)
  tempAtom.resize(numAtoms);  // NOTE: Worst case: all non-water

  // Define size of a water hydrogen group
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  // Iterate through all the atoms
  int i = 0;
  int waterIndex = 0;
  int nonWaterIndex = 0;
  while (i < numAtoms) {

    FullAtom &atom_i = atom[i];
    Mass mass = atom_i.mass;
    int hgs = atom_i.hydrogenGroupSize; 
    // Check to see if this hydrogen group is a water molecule
    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {

      // Move this hydrogen group up in the current atom list
      if (waterIndex != i) {
        atom[waterIndex    ] = atom[i    ];  // Oxygen
        atom[waterIndex + 1] = atom[i + 1];  // Hydrogen
        atom[waterIndex + 2] = atom[i + 2];  // Hydrogen
        if (wathgsize > 3) atom[waterIndex + 3] = atom[i + 3];  // lonepair
        if (wathgsize > 4) atom[waterIndex + 4] = atom[i + 4];  // drude
          // actual Drude water is arranged:  O D LP H H
      }

      waterIndex += wathgsize;
      i += wathgsize;

    } else {

      // Move this hydrogen group into non-water (scratch) atom list
      for (int j = 0; j < hgs; j++)
        tempAtom[nonWaterIndex + j] = atom[i + j];

      nonWaterIndex += hgs;
      i += hgs;
    }

  } // end iterating through atoms

  // Iterate through the non-water (scratch) atom list, adding the
  //   atoms to the end of the atom list.
  // NOTE: This could be done with a straight memcpy if the internal
  //   data structures of ResizeArray could be accessed directly.
  //   Or, perhaps add a member to ResizeArray that can add a consecutive
  //   list of elements starting at a particular index (would be cleaner).
  for (i = 0; i < nonWaterIndex; i++)
    atom[waterIndex + i] = tempAtom[i];

  // Set numWaterAtoms
  numWaterAtoms = waterIndex;
}


// This function will merge the given list of atoms (not assumed to
//   be separated) with the current list of atoms (already assumed
//   to be separated).
// NOTE: This function applies the transformations to the incoming
//   atoms as it is separating them.
void HomePatch::mergeAtomList(FullAtomList &al) {

  // Sanity check
  if (al.size() <= 0) return;  // Nothing to do

  const int orig_atomSize = atom.size();
  const int orig_alSize = al.size();

  // Resize the atom list (will eventually hold contents of both lists)
  atom.resize(orig_atomSize + orig_alSize); // NOTE: Will have contents of both


  #if 0  // version where non-waters are moved to scratch first

  
  // Basic Idea:  The current list is separated already so copy the
  //   non-water atoms out of it into the scratch atom array.  Then
  //   separate the incoming/given list (al), adding the waters to the
  //   end of the waters in atom list and non-waters to the end of the
  //   scratch list.  At this point, all waters are in atom list and all
  //   non-waters are in the scratch list so just copy the scratch list
  //   to the end of the atom list.
  // NOTE: If al is already separated and the number of waters in it
  //   is know, could simply move the non-waters in atoms back by that
  //   amount and directly copy the waters in al into the created gap
  //   and the non-waters in al to the end.  Leave this as an
  //   optimization for later since I'm not sure if this would actually
  //   do better as the combining code (for combining migration
  //   messages) would also have to merge the contents of the atom lists
  //   they carry.  Generally speaking, there is probably a faster way
  //   to do this, but this will get it working.

  // Copy all the non-waters in the current atom list into the
  //   scratch atom list.
  const int orig_atom_numNonWaters = orig_atomSize - numWaterAtoms;
  tempAtom.resize(orig_atom_numNonWaters + al.size()); // NOTE: Worst case
  for (int i = 0; i < orig_atom_numNonWaters; i++)
    tempAtom[i] = atom[numWaterAtoms + i];

  // Separate the contents of the given atom list (applying the
  // transforms as needed)
  int atom_waterIndex = numWaterAtoms;
  int atom_nonWaterIndex = orig_atom_numNonWaters;
  int i = 0;
  while (i < orig_alSize) {

    FullAtom &atom_i = al[i];
    int hgs = atom_i.hydrogenGroupSize;
    if ( hgs != atom_i.migrationGroupSize ) {
      NAMD_bug("HomePatch::mergeAtomList() not updated for migration groups!");
    }
    Mass mass = atom_i.mass;

    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {

      // Apply the transforms

      // Oxygen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogen (@ +1)
      al[i+1].position = lattice.reverse_transform(al[i+1].position, al[i+1].transform);
      al[i+1].position = lattice.apply_transform(al[i+1].position, mother_transform);
      al[i+1].transform = mother_transform;

      // Hydrogen (@ +2)
      al[i+2].position = lattice.reverse_transform(al[i+2].position, al[i+2].transform);
      al[i+2].position = lattice.apply_transform(al[i+2].position, mother_transform);
      al[i+2].transform = mother_transform;

      // Add to the end of the waters in the current list of atoms
      atom[atom_waterIndex    ] = al[i    ];
      atom[atom_waterIndex + 1] = al[i + 1];
      atom[atom_waterIndex + 2] = al[i + 2];

      atom_waterIndex += 3;
      i += 3;

    } else {

      // Apply the transforms

      // Non-Hydrogen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogens (@ +1 -> +(hgs-1))
      for (int j = 1; j < hgs; j++) {
        al[i+j].position = lattice.reverse_transform(al[i+j].position, al[i+j].transform);
        al[i+j].position = lattice.apply_transform(al[i+j].position, mother_transform);
        al[i+j].transform = mother_transform;
      }

      // Add to the end of the non-waters (scratch) atom list
      for (int j = 0; j < hgs; j++)
        tempAtom[atom_nonWaterIndex + j] = al[i + j];

      atom_nonWaterIndex += hgs;
      i += hgs;
    }

  } // end while iterating through given atom list

  // Copy all the non-waters to the end of the current atom list
  for (int i = 0; i < atom_nonWaterIndex; i++)
    atom[atom_waterIndex + i] = tempAtom[i];

  // Set numWaterAtoms and numAtoms
  numWaterAtoms = atom_waterIndex;
  numAtoms = atom.size();


  #else


  // Basic Idea:  Count the number of water atoms in the incoming atom
  //   list then move the non-waters back in the current atom list to
  //   make room for the incoming waters.  Once there is room in the
  //   current list, separate the incoming list as the atoms are being
  //   added to the current list.
  // NOTE:  Since the incoming atom list is likely to be small,
  //   iterating over its hydrogen groups twice should not be too bad.
  // NOTE:  This code assumes the current list is already separated,
  //   the incoming list may not be separated, and the transforms are
  //   applied to the incoming atoms as the separation occurs.

  // size of a water hydrogen group
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  // Count the number of waters in the given atom list
  int al_numWaterAtoms = 0;
  int i = 0;
  while (i < orig_alSize) {

    FullAtom &atom_i = al[i];
    int hgs = atom_i.hydrogenGroupSize;
    Mass mass = atom_i.mass;

    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {
      al_numWaterAtoms += wathgsize;
    }

    i += hgs;
  }

  // Move all of the non-waters in the current atom list back (to a
  //   higher index) by the number of waters in the given list.
  if (al_numWaterAtoms > 0) {
    for (i = orig_atomSize - 1; i >= numWaterAtoms; i--) {
      atom[i + al_numWaterAtoms] = atom[i];
    }
  }

  // Separte the atoms in the given atom list.  Perform the
  //   transformations on them and then add them to the appropriate
  //   location in the current atom list.
  int atom_waterIndex = numWaterAtoms;
  int atom_nonWaterIndex = orig_atomSize + al_numWaterAtoms;
  i = 0;
  while (i < orig_alSize) {

    FullAtom &atom_i = al[i];
    int hgs = atom_i.hydrogenGroupSize;
    if ( hgs != atom_i.migrationGroupSize ) {
      NAMD_bug("HomePatch::mergeAtomList() not updated for migration groups!");
    }
    Mass mass = atom_i.mass;

    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {

      // Apply the transforms

      // Oxygen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogen (@ +1)
      al[i+1].position = lattice.reverse_transform(al[i+1].position, al[i+1].transform);
      al[i+1].position = lattice.apply_transform(al[i+1].position, mother_transform);
      al[i+1].transform = mother_transform;

      // Hydrogen (@ +2)
      al[i+2].position = lattice.reverse_transform(al[i+2].position, al[i+2].transform);
      al[i+2].position = lattice.apply_transform(al[i+2].position, mother_transform);
      al[i+2].transform = mother_transform;

      // Add to the end of the waters in the current list of atoms
      atom[atom_waterIndex    ] = al[i    ];
      atom[atom_waterIndex + 1] = al[i + 1];
      atom[atom_waterIndex + 2] = al[i + 2];

      if (wathgsize > 3) atom[atom_waterIndex + 3] = al[i + 3];

      atom_waterIndex += wathgsize;
      i += wathgsize;

    } else {

      // Apply the transforms

      // Non-Hydrogen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogens (@ +1 -> +(hgs-1))
      for (int j = 1; j < hgs; j++) {
        al[i+j].position = lattice.reverse_transform(al[i+j].position, al[i+j].transform);
        al[i+j].position = lattice.apply_transform(al[i+j].position, mother_transform);
        al[i+j].transform = mother_transform;
      }

      // Add to the end of the non-waters (scratch) atom list
      for (int j = 0; j < hgs; j++)
        atom[atom_nonWaterIndex + j] = al[i + j];

      atom_nonWaterIndex += hgs;
      i += hgs;
    }

  } // end while iterating through given atom list

  // Set numWaterAtoms and numAtoms
  numWaterAtoms = atom_waterIndex;
  numAtoms = atom_nonWaterIndex;

  #endif
}

#endif



inline void lubksb(HGMatrixBigReal &a, int n, HGArrayInt &indx,
                                              HGArrayBigReal &b)
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii >= 0)
			for (j=ii;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


inline void ludcmp(HGMatrixBigReal &a, int n, HGArrayInt &indx, BigReal *d)
{

	int i,imax,j,k;
	double big,dum,sum,temp;
	HGArrayBigReal vv;
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) NAMD_die("Singular matrix in routine ludcmp\n");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}


inline void G_q(const HGArrayVector &refab, HGMatrixVector &gqij,
     const int n, const int m, const HGArrayInt &ial, const HGArrayInt &ibl) {
  int i; 
  // step through the rows of the matrix
  for(i=0;i<m;i++) {
    gqij[i][ial[i]]=2.0*refab[i];
    gqij[i][ibl[i]]=-gqij[i][ial[i]];
  }
};


// c-ji code for MOLLY 7-31-99
int average(CompAtom *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial) {
  //  input:  n = length of hydrogen group to be averaged (shaked)
  //          q[n] = vector array of original positions
  //          m = number of constraints
  //          imass[n] = inverse mass for each atom
  //          length2[m] = square of reference bond length for each constraint
  //          ial[m] = atom a in each constraint 
  //          ibl[m] = atom b in each constraint 
  //          refab[m] = vector of q_ial(i) - q_ibl(i) for each constraint
  //          tolf = function error tolerance for Newton's iteration
  //          ntrial = max number of Newton's iterations
  //  output: lambda[m] = double array of lagrange multipliers (used by mollify)
  //          qtilde[n] = vector array of averaged (shaked) positions

  int k,k1,i,j;
  BigReal errx,errf,d,tolx;

  HGArrayInt indx;
  HGArrayBigReal p;
  HGArrayBigReal fvec;
  HGMatrixBigReal fjac;
  HGArrayVector avgab;
  HGMatrixVector grhs;
  HGMatrixVector auxrhs;
  HGMatrixVector glhs;

  //  iout <<iINFO << "average: n="<<n<<" m="<<m<<std::endl<<endi;
  tolx=tolf; 
  
  // initialize lambda, globalGrhs

  for (i=0;i<m;i++) {
    lambda[i]=0.0;
  }

  // define grhs, auxrhs for all iterations
  // grhs= g_x(q)
  //
  G_q(refab,grhs,n,m,ial,ibl);
  for (k=1;k<=ntrial;k++) {
    //    usrfun(qtilde,q0,lambda,fvec,fjac,n,water); 
    HGArrayBigReal gij;
    // this used to be main loop of usrfun
    // compute qtilde given q0, lambda, IMASSes
    {
      BigReal multiplier;
      HGArrayVector tmp;
      for (i=0;i<m;i++) {
	multiplier = lambda[i];
	// auxrhs = M^{-1}grhs^{T}
	for (j=0;j<n;j++) {
	  auxrhs[i][j]=multiplier*imass[j]*grhs[i][j];
	}
      }
      for (j=0;j<n;j++) {
	//      tmp[j]=0.0;      
	for (i=0;i<m;i++) {
	  tmp[j]+=auxrhs[i][j];
	}
      }
 
      for (j=0;j<n;j++) {
	qtilde[j].position=q[j]+tmp[j];
      }
      //      delete [] tmp;
    }
  
    for ( i = 0; i < m; i++ ) {
      avgab[i] = qtilde[ial[i]].position - qtilde[ibl[i]].position;
    }

    //  iout<<iINFO << "Calling Jac" << std::endl<<endi;
    //  Jac(qtilde, q0, fjac,n,water);
    {
      //  Vector glhs[3*n+3];

      HGMatrixVector grhs2;

      G_q(avgab,glhs,n,m,ial,ibl);
#ifdef DEBUG0
      iout<<iINFO << "G_q:" << std::endl<<endi;
      for (i=0;i<m;i++) {
	iout<<iINFO << glhs[i*n+0] << " " << glhs[i*n+1] << " " << glhs[i*n+2] << std::endl<<endi;
      }
#endif
      //      G_q(refab,grhs2,m,ial,ibl);
      // update with the masses
      for (j=0; j<n; j++) { // number of atoms
	for (i=0; i<m;i++) { // number of constraints
	  grhs2[i][j] = grhs[i][j]*imass[j];
	}
      }

      // G_q(qtilde) * M^-1 G_q'(q0) =
      // G_q(qtilde) * grhs'
      for (i=0;i<m;i++) { // number of constraints
	for (j=0;j<m;j++) { // number of constraints
	  fjac[i][j] = 0; 
	  for (k1=0;k1<n;k1++) {
	    fjac[i][j] += glhs[i][k1]*grhs2[j][k1]; 
	  }
	}
      }
#ifdef DEBUG0  
      iout<<iINFO << "glhs" <<endi;
      for(i=0;i<9;i++) {
	iout<<iINFO << glhs[i] << ","<<endi;
      }
      iout<<iINFO << std::endl<<endi;
      for(i=0;i<9;i++) {
	iout<<iINFO << grhs2[i] << ","<<endi;
      }
      iout<<iINFO << std::endl<<endi;
#endif
      //      delete[] grhs2;
    }
    // end of Jac calculation
#ifdef DEBUG0
    iout<<iINFO << "Jac" << std::endl<<endi;
    for (i=0;i<m;i++) 
      for (j=0;j<m;j++)
	iout<<iINFO << fjac[i][j] << " "<<endi;
    iout<< std::endl<<endi;
#endif
    // calculate constraints in gij for n constraints this being a water
    //  G(qtilde, gij, n, water);
    for (i=0;i<m;i++) {
      gij[i]=avgab[i]*avgab[i]-length2[i];
    }
#ifdef DEBUG0
    iout<<iINFO << "G" << std::endl<<endi;
    iout<<iINFO << "( "<<endi;
    for(i=0;i<m-1;i++) {
      iout<<iINFO << gij[i] << ", "<<endi;
    }
    iout<<iINFO << gij[m-1] << ")" << std::endl<<endi;
#endif
    // fill the return vector
    for(i=0;i<m;i++) {
      fvec[i] = gij[i];
    }
    // free up the constraints
    //    delete[] gij;
    // continue Newton's iteration    
    errf=0.0;
    for (i=0;i<m;i++) errf += fabs(fvec[i]);
#ifdef DEBUG0
    iout<<iINFO << "errf: " << errf << std::endl<<endi;
#endif
    if (errf <= tolf) {
      break;
    }
    for (i=0;i<m;i++) p[i] = -fvec[i];
    //    iout<<iINFO << "Doing dcmp in average " << std::endl<<endi;
    ludcmp(fjac,m,indx,&d);
    lubksb(fjac,m,indx,p);

    errx=0.0;
    for (i=0;i<m;i++) {
      errx += fabs(p[i]);
    }
    for (i=0;i<m;i++)  
      lambda[i] += p[i];

#ifdef DEBUG0
    iout<<iINFO << "lambda:" << lambda[0] 
	 << " " << lambda[1] << " " << lambda[2] << std::endl<<endi;
    iout<<iINFO << "errx: " << errx << std::endl<<endi;
#endif
    if (errx <= tolx) break;
#ifdef DEBUG0
    iout<<iINFO << "Qtilde:" << std::endl<<endi;
    iout<<iINFO << qtilde[0].position << " " << qtilde[1].position << " " << qtilde[2].position << std::endl<<endi; 
#endif
  }
#ifdef DEBUG
  iout<<iINFO << "LAMBDA:" << lambda[0] << " " << lambda[1] << " " << lambda[2] << std::endl<<endi;
#endif

  return k; // 
}

void mollify(CompAtom *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab) {
  int i,j,k;
  BigReal d;
  HGMatrixBigReal fjac;
  Vector zero(0.0,0.0,0.0);
  
  HGArrayVector tmpforce;
  HGArrayVector tmpforce2;
  HGArrayVector y;
  HGMatrixVector grhs;
  HGMatrixVector glhs;
  HGArrayBigReal aux;
  HGArrayInt indx;

  for(i=0;i<n;i++) {
    tmpforce[i]=imass[i]*force[i];
  }

  HGMatrixVector grhs2;
  HGArrayVector avgab;

  for ( i = 0; i < m; i++ ) {
	avgab[i] = qtilde[ial[i]].position - qtilde[ibl[i]].position;
  }

  G_q(avgab,glhs,n,m,ial,ibl);
  G_q(refab,grhs,n,m,ial,ibl);
  // update with the masses
  for (j=0; j<n; j++) { // number of atoms
	for (i=0; i<m;i++) { // number of constraints
	  grhs2[i][j] = grhs[i][j]*imass[j];
	}
  }

  // G_q(qtilde) * M^-1 G_q'(q0) =
  // G_q(qtilde) * grhs'
  for (i=0;i<m;i++) { // number of constraints
	for (j=0;j<m;j++) { // number of constraints
	  fjac[j][i] = 0; 
	  for (k=0;k<n;k++) {
	    fjac[j][i] += glhs[i][k]*grhs2[j][k]; 
	  }
	}
  }

  // aux=gqij*tmpforce
  //  globalGrhs::computeGlobalGrhs(q0,n,water);
  //  G_q(refab,grhs,m,ial,ibl);
  for(i=0;i<m;i++) {
    aux[i]=0.0;
    for(j=0;j<n;j++) {
      aux[i]+=grhs[i][j]*tmpforce[j];
    }
  }

  ludcmp(fjac,m,indx,&d);
  lubksb(fjac,m,indx,aux);

  for(j=0;j<n;j++) {
    y[j] = zero;
    for(i=0;i<m;i++) {
      y[j] += aux[i]*glhs[i][j];
    }
  }
  for(i=0;i<n;i++) {
    y[i]=force[i]-y[i];
  }
    
  // gqq12*y
  for(i=0;i<n;i++) {
    tmpforce2[i]=imass[i]*y[i];
  }

  // here we assume that tmpforce is initialized to zero.
  for (i=0;i<n;i++) {
    tmpforce[i]=zero;
  }
  
  for (j=0;j<m;j++) {
    Vector tmpf = 2.0*lambda[j]*(tmpforce2[ial[j]]-tmpforce2[ibl[j]]);
    tmpforce[ial[j]] += tmpf;
    tmpforce[ibl[j]] -= tmpf;
  }
  // c-ji the other bug for 2 constraint water was this line (2-4-99)
  //  for(i=0;i<m;i++) {
  for(i=0;i<n;i++) {
    force[i]=tmpforce[i]+y[i];
  }

}
