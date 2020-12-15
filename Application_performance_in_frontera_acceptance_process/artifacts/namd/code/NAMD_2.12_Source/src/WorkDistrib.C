/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/WorkDistrib.C,v $
 * $Author: jim $
 * $Date: 2016/09/29 20:31:47 $
 * $Revision: 1.1290 $
 *****************************************************************************/

/** \file WorkDistrib.C
 *  Currently, WorkDistrib generates the layout of the Patches,
 *  directs the construction and distribution of Computes and
 *  associates Computes with Patches.
 */

#include <stdio.h>
 
#include "InfoStream.h"
#include "Communicate.h"
#include "ProcessorPrivate.h"
#include "BOCgroup.h"
#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"
#include "Lattice.h"
#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "main.decl.h"
#include "main.h"
#include "Node.h"
#include "PatchMgr.h"
#include "PatchMap.inl"
#include "NamdTypes.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "Compute.h"
#include "ComputeMap.h"
#include "RecBisection.h"
#include "Random.h"
#include "varsizemsg.h"
#include "ProxyMgr.h"
#include "Priorities.h"
#include "SortAtoms.h"
#include <algorithm>
#include "TopoManager.h"
#include "ComputePmeCUDAMgr.h"

#include "DeviceCUDA.h"
#ifdef NAMD_CUDA
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;
#endif

//#define DEBUGM
#define MIN_DEBUG_LEVEL 2
#include "Debug.h"
#ifdef MEM_OPT_VERSION
extern int isOutputProcessor(int); 
#endif
class ComputeMapChangeMsg : public CMessage_ComputeMapChangeMsg
{
public:

  int numNewNodes;
  int numNewNumPartitions;
  int *newNodes;
  char *newNumPartitions;

//  VARSIZE_DECL(ComputeMapChangeMsg);
};

/*
VARSIZE_MSG(ComputeMapChangeMsg,
  VARSIZE_ARRAY(newNodes);
)
*/

static int randtopo;

static void build_ordering(void *) {
  WorkDistrib::buildNodeAwarePeOrdering();
}

void topo_getargs(char **argv) {
  randtopo = CmiGetArgFlag(argv, "+randtopo");
  if ( CkMyPe() >= CkNumPes() ) return;
  CcdCallOnCondition(CcdTOPOLOGY_AVAIL, (CcdVoidFn)build_ordering, (void*)0);
}

static int eventMachineProgress;

//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib()
{
  CkpvAccess(BOCclass_group).workDistrib = thisgroup;
  patchMapArrived = false;
  computeMapArrived = false;

#if CMK_SMP
#define MACHINE_PROGRESS
#else
#define MACHINE_PROGRESS { traceUserEvent(eventMachineProgress);  CmiMachineProgressImpl(); }
  if ( CkMyNodeSize() > 1 ) NAMD_bug("CkMyNodeSize() > 1 for non-smp build");
  eventMachineProgress = traceRegisterUserEvent("CmiMachineProgressImpl",233);
#endif
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{ }

static int compare_bit_reversed(int a, int b) {
  int d = a ^ b;
  int c = 1;
  if ( d ) while ( ! (d & c) ) {
    c = c << 1;
  }
  return (a & c) - (b & c);
}

static bool less_than_bit_reversed(int a, int b) {
  int d = a ^ b;
  int c = 1;
  if ( d ) while ( ! (d & c) ) {
    c = c << 1;
  }
  return d && (b & c);
}

struct pe_sortop_bit_reversed {
  int *rankInPhysOfNode;
  pe_sortop_bit_reversed(int *r) : rankInPhysOfNode(r) {}
  inline bool operator() (int a, int b) const {
    int c = compare_bit_reversed(CmiRankOf(a),CmiRankOf(b));
    if ( c < 0 ) return true;
    if ( c > 0 ) return false;
    c = compare_bit_reversed(
        rankInPhysOfNode[CmiNodeOf(a)],rankInPhysOfNode[CmiNodeOf(b)]);
    if ( c < 0 ) return true;
    if ( c > 0 ) return false;
    c = compare_bit_reversed(CmiPhysicalNodeID(a),CmiPhysicalNodeID(b));
    return ( c < 0 );
  }
};

int WorkDistrib::peOrderingInit;
int* WorkDistrib::peDiffuseOrdering;
int* WorkDistrib::peDiffuseOrderingIndex;
int* WorkDistrib::peCompactOrdering;
int* WorkDistrib::peCompactOrderingIndex;

#ifdef NAMD_CUDA
extern void cuda_initialize();
#endif

void mic_initialize();

void WorkDistrib::peOrderingReady() {
  //CkPrintf("WorkDistrib::peOrderingReady on %d\n", CkMyPe());
#ifdef NAMD_CUDA
  cuda_initialize();
#endif
#ifdef NAMD_MIC
  mic_initialize();
#endif
}

void WorkDistrib::buildNodeAwarePeOrdering() {

 CmiMemLock();
 if ( ! peOrderingInit ) {
  //CkPrintf("WorkDistrib::buildNodeAwarePeOrdering on pe %d\n", CkMyPe());

  const int numPhys = CmiNumPhysicalNodes();
  const int numNode = CmiNumNodes();
  const int numPe = CmiNumPes();
  ResizeArray<int> numNodeInPhys(numPhys);
  ResizeArray<int> rankInPhysOfNode(numNode);

  peDiffuseOrdering = new int[numPe];
  peDiffuseOrderingIndex = new int[numPe];
  peCompactOrdering = new int[numPe];
  peCompactOrderingIndex = new int[numPe];

  int k = 0;
  for ( int ph=0; ph<numPhys; ++ph ) {
    int *pes, npes;
    CmiGetPesOnPhysicalNode(ph, &pes, &npes);
    for ( int i=0; i<npes; ++i, ++k ) {
      peCompactOrdering[k] = pes[i];
    }
    numNodeInPhys[ph] = 0;
    for ( int i=0, j=0; i<npes; i += CmiNodeSize(CmiNodeOf(pes[i])), ++j ) {
      rankInPhysOfNode[CmiNodeOf(pes[i])] = j;
      numNodeInPhys[ph] += 1;
    }
  }

  if ( randtopo && numPhys > 2 ) {
    if ( ! CkMyNode() ) {
      iout << iWARN << "RANDOMIZING PHYSICAL NODE ORDERING\n" << endi;
    }
    ResizeArray<int> randPhysOrder(numPhys);
    for ( int j=0; j<numPhys; ++j ) {
      randPhysOrder[j] = j;
    }
    Random(314159265).reorder(randPhysOrder.begin()+2, numPhys-2);
    for ( int j=0, k=0; j<numPhys; ++j ) {
      const int ph = randPhysOrder[j];
      int *pes, npes;
      CmiGetPesOnPhysicalNode(ph, &pes, &npes);
      for ( int i=0; i<npes; ++i, ++k ) {
        peCompactOrdering[k] = pes[i];
      }
    }
  }

  for ( int i=0; i<numPe; ++i ) {
    peDiffuseOrdering[i] = i;
  }
  std::sort(peDiffuseOrdering, peDiffuseOrdering+numPe,
    pe_sortop_bit_reversed(rankInPhysOfNode.begin()));

  for ( int i=0; i<numPe; ++i ) {
    peDiffuseOrderingIndex[peDiffuseOrdering[i]] = i;
    peCompactOrderingIndex[peCompactOrdering[i]] = i;
  }

  if ( 0 && CmiMyNode() == 0 ) for ( int i=0; i<numPe; ++i ) {
    CkPrintf("order %5d %5d %5d %5d %5d\n", i,
      peDiffuseOrdering[i],
      peDiffuseOrderingIndex[i],
      peCompactOrdering[i],
      peCompactOrderingIndex[i]);
  }

  peOrderingInit = 1;
 }
 CmiMemUnlock();
 peOrderingReady();

}

struct pe_sortop_coord_x {
  ScaledPosition *spos;
  pe_sortop_coord_x(ScaledPosition *s) : spos(s) {}
  inline bool operator() (int a, int b) const {
    return ( spos[a].x < spos[b].x );
  }
};

struct pe_sortop_coord_y {
  ScaledPosition *spos;
  pe_sortop_coord_y(ScaledPosition *s) : spos(s) {}
  inline bool operator() (int a, int b) const {
    return ( spos[a].y < spos[b].y );
  }
};

static void recursive_bisect_coord(
    int x_begin, int x_end, int y_begin, int y_end,
    int *pe_begin, ScaledPosition *coord,
    int *result, int ydim
  ) {
  int x_len = x_end - x_begin;
  int y_len = y_end - y_begin;
  if ( x_len == 1 && y_len == 1 ) {
    // done, now put this pe in the right place
    if ( 0 ) CkPrintf("pme %5d %5d on pe %5d at %f %f\n", x_begin, y_begin, *pe_begin,
      coord[*pe_begin].x, coord[*pe_begin].y);
    result[x_begin*ydim + y_begin] = *pe_begin;
    return;
  }
  int *pe_end = pe_begin + x_len * y_len;
  if ( x_len >= y_len ) {
    std::sort(pe_begin, pe_end, pe_sortop_coord_x(coord));
    int x_split = x_begin + x_len / 2;
    int* pe_split = pe_begin + (x_split - x_begin) * y_len;
    //CkPrintf("x_split %5d %5d %5d\n", x_begin, x_split, x_end);
    recursive_bisect_coord(x_begin, x_split, y_begin, y_end, pe_begin, coord, result, ydim);
    recursive_bisect_coord(x_split, x_end, y_begin, y_end, pe_split, coord, result, ydim);
  } else {
    std::sort(pe_begin, pe_end, pe_sortop_coord_y(coord));
    int y_split = y_begin + y_len / 2;
    int* pe_split = pe_begin + (y_split - y_begin) * x_len;
    //CkPrintf("y_split %5d %5d %5d\n", y_begin, y_split, y_end);
    recursive_bisect_coord(x_begin, x_end, y_begin, y_split, pe_begin, coord, result, ydim);
    recursive_bisect_coord(x_begin, x_end, y_split, y_end, pe_split, coord, result, ydim);
  }
}

void WorkDistrib::sortPmePes(int *pmepes, int xdim, int ydim) {
  int numpes = CkNumPes();
  ResizeArray<int> count(numpes);
  ResizeArray<ScaledPosition> sumPos(numpes);
  ResizeArray<ScaledPosition> avgPos(numpes);
  for ( int i=0; i<numpes; ++i ) {
    count[i] = 0;
    sumPos[i] = 0;
    avgPos[i] = 0;
  }
  PatchMap *patchMap = PatchMap::Object();
  for ( int i=0, npatches=patchMap->numPatches(); i<npatches; ++i ) {
    int pe = patchMap->node(i);
    count[pe] += 1;
    sumPos[pe] += patchMap->center(i);
  }
  const int npmepes = xdim*ydim;
  ResizeArray<int> sortpes(npmepes);
  for ( int i=0; i<npmepes; ++i ) {
    int pe = sortpes[i] = pmepes[i];
    int cnt = count[pe];
    ScaledPosition sum = sumPos[pe];
    if ( cnt == 0 ) {
      // average over node
      int node = CkNodeOf(pe);
      int nsize = CkNodeSize(node);
      int pe2 = CkNodeFirst(node);
      for ( int j=0; j<nsize; ++j, ++pe2 )  {
        cnt += count[pe2];
        sum += sumPos[pe2];
      }
    }
    if ( cnt == 0 ) {
      // average over physical node
      int node = CmiPhysicalNodeID(pe);
      int nsize, *nlist;
      CmiGetPesOnPhysicalNode(node, &nlist, &nsize);
      for ( int j=0; j<nsize; ++j )  {
        int pe2 = nlist[j];
        cnt += count[pe2];
        sum += sumPos[pe2];
      }
    }
    if ( cnt ) {
      avgPos[pe] = sum / cnt;
    }
  }
  recursive_bisect_coord(0, xdim, 0, ydim, sortpes.begin(), avgPos.begin(), pmepes, ydim);
}


//----------------------------------------------------------------------
void WorkDistrib::saveComputeMapChanges(int ep, CkGroupID chareID)
{
  saveComputeMapReturnEP = ep;
  saveComputeMapReturnChareID = chareID;

  ComputeMapChangeMsg *mapMsg = new (0, 0, 0) ComputeMapChangeMsg;
  CProxy_WorkDistrib(thisgroup).recvComputeMapChanges(mapMsg);

/*
    // store the latest compute map
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->storeComputeMap) {
    computeMap->saveComputeMap(simParams->computeMapFilename);
    CkPrintf("ComputeMap has been stored in %s.\n", simParams->computeMapFilename);
  }
*/
}

void WorkDistrib::recvComputeMapChanges(ComputeMapChangeMsg *msg) {

  delete msg;

  ComputeMap *computeMap = ComputeMap::Object();

  int i;
  int nc = computeMap->numComputes();
  
  if ( ! CkMyPe() ) { // send
    // CkPrintf("At %f on %d WorkDistrib::recvComputeMapChanges %d\n", CmiWallTimer(), CkMyPe(), nc);
    MOStream *msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, COMPUTEMAPTAG, BUFSIZE);
    msg->put(nc);
    for (i=0; i<nc; i++) {
      int data = computeMap->newNode(i);
      msg->put(data);
    }
    msg->put(nc);
    for (i=0; i<nc; i++) {
      char data = computeMap->newNumPartitions(i);
      msg->put(data);
    }
    msg->put(nc);
    msg->end();
    delete msg;
    // CkPrintf("At %f on %d done WorkDistrib::recvComputeMapChanges %d\n", CmiWallTimer(), CkMyPe(), nc);
  } else if ( ! CkMyRank() ) { // receive
    // if ( CkMyNode() == 1 ) CkPrintf("At %f on %d WorkDistrib::recvComputeMapChanges %d\n", CmiWallTimer(), CkMyPe(), nc);
    MIStream *msg = CkpvAccess(comm)->newInputStream(0, COMPUTEMAPTAG);
    msg->get(i);
    if ( i != nc ) NAMD_bug("WorkDistrib::recvComputeMapChanges check 1 failed\n");
    for (i=0; i<nc; i++) {
      int data;
      msg->get(data);
      computeMap->setNewNode(i,data);
    }
    msg->get(i);
    if ( i != nc ) NAMD_bug("WorkDistrib::recvComputeMapChanges check 2 failed\n");
    for (i=0; i<nc; i++) {
      char data;
      msg->get(data);
      computeMap->setNewNumPartitions(i,data);
    }
    msg->get(i);
    if ( i != nc ) NAMD_bug("WorkDistrib::recvComputeMapChanges check 3 failed\n");
    delete msg;
    // if ( CkMyNode() == 1 ) CkPrintf("At %f on %d done WorkDistrib::recvComputeMapChanges %d\n", CmiWallTimer(), CkMyPe(), nc);
  }

  CkCallback cb(CkIndex_WorkDistrib::doneSaveComputeMap(NULL), 0, thisgroup);
  contribute(0, NULL, CkReduction::random, cb);
}

void WorkDistrib::doneSaveComputeMap(CkReductionMsg *msg) {
  delete msg;

  CkSendMsgBranch(saveComputeMapReturnEP, CkAllocMsg(0,0,0), 0, saveComputeMapReturnChareID);
}

#ifdef MEM_OPT_VERSION
//All basic info already exists for each atom inside the FullAtomList because
//it is loaded when reading the binary per-atom file. This function will fill
//the info regarding transform, nonbondedGroupSize etc. Refer to
//WorkDistrib::createAtomLists
void WorkDistrib::fillAtomListForOnePatch(int pid, FullAtomList &alist){
  PatchMap *patchMap = PatchMap::Object();

  ScaledPosition center(0.5*(patchMap->min_a(pid)+patchMap->max_a(pid)),
                          0.5*(patchMap->min_b(pid)+patchMap->max_b(pid)),
                          0.5*(patchMap->min_c(pid)+patchMap->max_c(pid)));

    int n = alist.size();
    FullAtom *a = alist.begin();
/* 
    //Those options are not supported in MEM_OPT_VERSIOn -Chao Mei 
//Modifications for alchemical fep
    Bool alchFepOn = params->alchFepOn;
    Bool alchThermIntOn = params->alchThermIntOn;
//fepe
    Bool lesOn = params->lesOn;
    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);
*/
    SimParameters *params = Node::Object()->simParameters;
    const Lattice lattice = params->lattice;
    Transform mother_transform;
    for(int j=0; j < n; j++)
    {
      int aid = a[j].id;
      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates
      
      // a[j].fixedPosition = a[j].position;  ParallelIOMgr stores ref coord here.

      if ( a[j].migrationGroupSize ) {
       if ( a[j].migrationGroupSize != a[j].hydrogenGroupSize ) {
            Position pos = a[j].position;
            int mgs = a[j].migrationGroupSize;
            int c = 1;

            for ( int k=a[j].hydrogenGroupSize; k<mgs;
                                k+=a[j+k].hydrogenGroupSize ) {
              pos += a[j+k].position;
              ++c;
            }

            pos *= 1./c;
            mother_transform = a[j].transform;  // should be 0,0,0
            pos = lattice.nearest(pos,center,&mother_transform);
            a[j].position = lattice.apply_transform(a[j].position,mother_transform);
            a[j].transform = mother_transform;        
       } else {
        a[j].position = lattice.nearest(a[j].position, center, &(a[j].transform));
        mother_transform = a[j].transform;
       }
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }

/* 
    //Those options are not supported in MEM_OPT_VERSIOn -Chao Mei 
//Modifications for alchemical fep
      if ( alchOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      } 
      else {
        a[j].partition = 0;
      }
//fepe
*/
      a[j].partition = 0;

      //set langevinParams based on atom status
      if(params->langevinOn) {
        BigReal bval = params->langevinDamping;
        if(!params->langevinHydrogen && 
           ((a[j].status & HydrogenAtom)!=0)) {
          bval = 0;
        }else if ((a[j].status & LonepairAtom)!=0) {
          bval = 0;
        }else if ((a[j].status & DrudeAtom)!=0) {
          bval = params->langevinDamping;
        }
        a[j].langevinParam = bval;
      }

    }

    int size, allfixed;
    for(int j=0; j < n; j+=size) {
      size = a[j].hydrogenGroupSize;
      if ( ! size ) {
        NAMD_bug("Mother atom with hydrogenGroupSize of 0!");
      }
      allfixed = 1;
      for (int k = 0; k < size; ++k ) {
        allfixed = ( allfixed && (a[j+k].atomFixed) );
      }
      for (int k = 0; k < size; ++k ) {
        a[j+k].groupFixed = allfixed ? 1 : 0;
      }
    }

    if ( params->outputPatchDetails ) {
      int patchId = pid;
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(int j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " on proc " << CkMyPe()
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }

}

void WorkDistrib::random_velocities_parallel(BigReal Temp,InputAtomList &inAtoms)
{
  int i, j;             //  Loop counter
  BigReal kbT;          //  Boltzman constant * Temp
  BigReal randnum;      //  Random number from -6.0 to 6.0
  BigReal kbToverM;     //  sqrt(Kb*Temp/Mass)
  SimParameters *simParams = Node::Object()->simParameters;
  Bool lesOn = simParams->lesOn;
  Random vel_random(simParams->randomSeed);
  int lesReduceTemp = lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  kbT = Temp*BOLTZMANN;
  int count=0;
  int totalAtoms = inAtoms.size();
  for(i=0;i<totalAtoms;i++)
  {
    Real atomMs=inAtoms[i].mass;

    if (atomMs <= 0.) {
      kbToverM = 0.;
    } else {
      /*
       * lesOn is not supported in MEM_OPT_VERSION, so the original assignment
       * is simplified. --Chao Mei
       */
      //kbToverM = sqrt(kbT *
        //( lesOn && structure->get_fep_type(aid) ? tempFactor : 1.0 ) /
        //                  atomMs );
      kbToverM = sqrt(kbT * 1.0 / atomMs);
    }
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    inAtoms[i].velocity.x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    inAtoms[i].velocity.y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;
    
    inAtoms[i].velocity.z = randnum*kbToverM;
  }
}
#endif

//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
FullAtomList *WorkDistrib::createAtomLists(const char *basename)
{
  int i;
  StringList *current;	//  Pointer used to retrieve configuration items
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
  SimParameters *params = node->simParameters;
  Molecule *molecule = node->molecule;
  PDB *pdb = node->pdb;

  int numPatches = patchMap->numPatches();
  int numAtoms = pdb->num_atoms();

  Vector *positions = new Position[numAtoms];
  Vector *velocities = new Velocity[numAtoms];

 if ( basename ) {
  if ( params->binaryOutput ) {
    read_binary_file((std::string(basename)+".coor").c_str(), positions, numAtoms);
    read_binary_file((std::string(basename)+".vel").c_str(), velocities, numAtoms);
  } else {
    PDB coorpdb((std::string(basename)+".coor").c_str());
    if ( coorpdb.num_atoms() != numAtoms ) {
      NAMD_die("Incorrect atom count in coordinate pdb file");
    }
    coorpdb.get_all_positions(positions);
    velocities_from_PDB((std::string(basename)+".vel").c_str(), velocities, numAtoms);
  }
 } else {
  pdb->get_all_positions(positions);

  if ( params->initialTemp < 0.0 ) {
    Bool binvels=FALSE;

    //  Reading the veolcities from a PDB
    current = node->configList->find("velocities");

    if (current == NULL) {
      current = node->configList->find("binvelocities");
      binvels = TRUE;
    }

    if (!binvels) {
      velocities_from_PDB(current->data, velocities, numAtoms);
    }
    else {
      velocities_from_binfile(current->data, velocities, numAtoms);
    }
  }
  else {
    // Random velocities for a given temperature
    random_velocities(params->initialTemp, molecule, velocities, numAtoms);
  }
 }

  //  If COMMotion == no, remove center of mass motion
  if (!(params->comMove)) {
    remove_com_motion(velocities, molecule, numAtoms);
  }

  FullAtomList *atoms = new FullAtomList[numPatches];

  const Lattice lattice = params->lattice;

    if ( params->staticAtomAssignment ) {
      FullAtomList sortAtoms;
      for ( i=0; i < numAtoms; i++ ) {
        HydrogenGroupID &h = molecule->hydrogenGroup[i];
        if ( ! h.isMP ) continue;
        FullAtom a;
        a.id = i;
        a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
        a.position = positions[h.atomID];
        sortAtoms.add(a);
      } 
      int *order = new int[sortAtoms.size()];
      for ( i=0; i < sortAtoms.size(); i++ ) {
        order[i] = i;
      } 
      int *breaks = new int[numPatches];
      sortAtomsForPatches(order,breaks,sortAtoms.begin(),
                        sortAtoms.size(),numAtoms,
                        patchMap->gridsize_c(),
                        patchMap->gridsize_b(),
                        patchMap->gridsize_a());

      i = 0;
      for ( int pid = 0; pid < numPatches; ++pid ) {
        int iend = breaks[pid];
        for ( ; i<iend; ++i ) {
          FullAtom &sa = sortAtoms[order[i]];
          int mgs = sa.migrationGroupSize;
/*
CkPrintf("patch %d (%d %d %d) has group %d atom %d size %d at %.2f %.2f %.2f\n",
          pid, patchMap->index_a(pid), patchMap->index_b(pid), 
          patchMap->index_c(pid), order[i], sa.id, mgs,
          sa.position.x, sa.position.y, sa.position.z);
*/
          for ( int k=0; k<mgs; ++k ) {
            HydrogenGroupID &h = molecule->hydrogenGroup[sa.id + k];
            int aid = h.atomID;
            FullAtom a;
            a.id = aid;
            a.position = positions[aid];
            a.velocity = velocities[aid];
            a.vdwType = molecule->atomvdwtype(aid);
            a.status = molecule->getAtoms()[aid].status;
            a.langevinParam = molecule->langevin_param(aid);
            a.hydrogenGroupSize = h.isGP ? h.atomsInGroup : 0;
            a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
            if(params->rigidBonds != RIGID_NONE) {
              a.rigidBondLength = molecule->rigid_bond_length(aid);
            }else{
              a.rigidBondLength = 0.0;
            }
            atoms[pid].add(a);
          }
        }
CkPrintf("patch %d (%d %d %d) has %d atoms\n",
          pid, patchMap->index_a(pid), patchMap->index_b(pid), 
          patchMap->index_c(pid), atoms[pid].size());
      }
      delete [] order;
      delete [] breaks;
    } else
    {
    // split atoms into patches based on migration group and position
    int aid, pid=0;
    for(i=0; i < numAtoms; i++)
      {
      // Assign atoms to patches without splitting hydrogen groups.
      // We know that the hydrogenGroup array is sorted with group parents
      // listed first.  Thus, only change the pid if an atom is a group parent.
      HydrogenGroupID &h = molecule->hydrogenGroup[i];
      aid = h.atomID;
      FullAtom a;
      a.id = aid;
      a.position = positions[aid];
      a.velocity = velocities[aid];
      a.vdwType = molecule->atomvdwtype(aid);
      a.status = molecule->getAtoms()[aid].status;
      a.langevinParam = molecule->langevin_param(aid);
      a.hydrogenGroupSize = h.isGP ? h.atomsInGroup : 0;
      a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
      if(params->rigidBonds != RIGID_NONE) {
        a.rigidBondLength = molecule->rigid_bond_length(aid);
      }else{
        a.rigidBondLength = 0.0;
      }
      if (h.isMP) {
	pid = patchMap->assignToPatch(positions[aid],lattice);
      } // else: don't change pid
      atoms[pid].add(a);
      }
    }

  delete [] positions;
  delete [] velocities;

  for(i=0; i < numPatches; i++)
  {
    ScaledPosition center(0.5*(patchMap->min_a(i)+patchMap->max_a(i)),
			  0.5*(patchMap->min_b(i)+patchMap->max_b(i)),
			  0.5*(patchMap->min_c(i)+patchMap->max_c(i)));

    int n = atoms[i].size();
    FullAtom *a = atoms[i].begin();
    int j;
//Modifications for alchemical fep
    Bool alchOn = params->alchOn;
//fepe
    Bool lesOn = params->lesOn;
  
    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);

    Transform mother_transform;
    for(j=0; j < n; j++)
    {
      int aid = a[j].id;

      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates

      a[j].atomFixed = molecule->is_atom_fixed(aid) ? 1 : 0;
      a[j].fixedPosition = a[j].position;

      if ( a[j].migrationGroupSize ) {
       if ( a[j].migrationGroupSize != a[j].hydrogenGroupSize ) {
            Position pos = a[j].position;
            int mgs = a[j].migrationGroupSize;
            int c = 1;
            for ( int k=a[j].hydrogenGroupSize; k<mgs;
                                k+=a[j+k].hydrogenGroupSize ) {
              pos += a[j+k].position;
              ++c;
            }
            pos *= 1./c;
            mother_transform = a[j].transform;  // should be 0,0,0
            pos = lattice.nearest(pos,center,&mother_transform);
            a[j].position = lattice.apply_transform(a[j].position,mother_transform);
            a[j].transform = mother_transform;
       } else {
        a[j].position = lattice.nearest(
		a[j].position, center, &(a[j].transform));
        mother_transform = a[j].transform;
       }
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }

      a[j].mass = molecule->atommass(aid);
      a[j].charge = molecule->atomcharge(aid);

//Modifications for alchemical fep
      if ( alchOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      } 
      else {
        a[j].partition = 0;
      }
//fepe

    }

    int size, allfixed, k;
    for(j=0; j < n; j+=size) {
      size = a[j].hydrogenGroupSize;
      if ( ! size ) {
        NAMD_bug("Mother atom with hydrogenGroupSize of 0!");
      }
      allfixed = 1;
      for ( k = 0; k < size; ++k ) {
        allfixed = ( allfixed && (a[j+k].atomFixed) );
      }
      for ( k = 0; k < size; ++k ) {
        a[j+k].groupFixed = allfixed ? 1 : 0;
      }
    }

    if ( params->outputPatchDetails ) {
      int patchId = i;
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }
  }

  return atoms;

}

//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
void WorkDistrib::createHomePatches(void)
{
  int i;
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  FullAtomList *atoms = createAtomLists();
    
#ifdef MEM_OPT_VERSION
/*  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  node->molecule->delEachAtomSigs();
  node->molecule->delMassChargeSpace();
*/
#endif

  int maxAtoms = -1;
  int maxPatch = -1;
  for(i=0; i < numPatches; i++) {
    int numAtoms = atoms[i].size();
    if ( numAtoms > maxAtoms ) { maxAtoms = numAtoms; maxPatch = i; }
  }
  iout << iINFO << "LARGEST PATCH (" << maxPatch <<
	") HAS " << maxAtoms << " ATOMS\n" << endi;

  for(i=0; i < numPatches; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(3,"Created " << i << " patches so far.\n");
    }

    patchMgr->createHomePatch(i,atoms[i]);
  }

  delete [] atoms;
}

void WorkDistrib::distributeHomePatches() {
  // ref BOC
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
  // ref singleton
  PatchMap *patchMap = PatchMap::Object();

  // Move patches to the proper node
  for(int i=0;i < patchMap->numPatches(); i++)
  {
    if (patchMap->node(i) != node->myid() )
    {
      DebugM(3,"patchMgr->movePatch("
	<< i << "," << patchMap->node(i) << ")\n");
      patchMgr->movePatch(i,patchMap->node(i));
    }
  }
  patchMgr->sendMovePatches();
}

void WorkDistrib::reinitAtoms(const char *basename) {

  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  FullAtomList *atoms = createAtomLists(basename);

  for(int i=0; i < numPatches; i++) {
    patchMgr->sendAtoms(i,atoms[i]);
  }

  delete [] atoms;

}


//----------------------------------------------------------------------

class PatchMapMsg : public CMessage_PatchMapMsg {
  public:
    char *patchMapData;
};

void WorkDistrib::sendPatchMap(void)
{
  if ( CkNumPes() == 1 ) {
    patchMapArrived = true;
    return;
  }

  //Automatically enable spanning tree
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  if( ( PatchMap::Object()->numPatches() <= CkNumPes()/4
#ifdef NODEAWARE_PROXY_SPANNINGTREE 
      || CkNumPes() > CkNumNodes()
      ) && ( CkNumNodes() > 1
#endif
    ) && params->isSendSpanningTreeUnset() )
    ProxyMgr::Object()->setSendSpanning();

#ifdef NODEAWARE_PROXY_SPANNINGTREE 
  if ( CkNumPes() > CkNumNodes() && CkNumNodes() > 1
        && params->isRecvSpanningTreeUnset() )
    ProxyMgr::Object()->setRecvSpanning();
#endif

  int size = PatchMap::Object()->packSize();

  PatchMapMsg *mapMsg = new (size, 0) PatchMapMsg;

  PatchMap::Object()->pack(mapMsg->patchMapData, size);

  CProxy_WorkDistrib workProxy(thisgroup);
  workProxy[0].savePatchMap(mapMsg);
}

// saveMaps() is called when the map message is received
void WorkDistrib::savePatchMap(PatchMapMsg *msg)
{
  // Use a resend to forward messages before processing.  Otherwise the
  // map distribution is slow on many CPUs.  We need to use a tree
  // rather than a broadcast because some implementations of broadcast
  // generate a copy of the message on the sender for each recipient.
  // This is because MPI doesn't allow re-use of an outstanding buffer.

  if ( CkMyRank() ) patchMapArrived = true;

  if ( patchMapArrived && CkMyPe() ) {
    PatchMap::Object()->unpack(msg->patchMapData);

    //Automatically enable spanning tree
    CProxy_Node nd(CkpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    SimParameters *params = node->simParameters;
    if( ( PatchMap::Object()->numPatches() <= CkNumPes()/4
#ifdef NODEAWARE_PROXY_SPANNINGTREE 
        || CkNumPes() > CkNumNodes()
        ) && ( CkNumNodes() > 1
#endif
      ) && params->isSendSpanningTreeUnset() )
      ProxyMgr::Object()->setSendSpanning();

#ifdef NODEAWARE_PROXY_SPANNINGTREE 
    if ( CkNumPes() > CkNumNodes() && CkNumNodes() > 1
          && params->isRecvSpanningTreeUnset() )
      ProxyMgr::Object()->setRecvSpanning();
#endif
  }

  if ( patchMapArrived ) {
    if ( CkMyRank() + 1 < CkNodeSize(CkMyNode()) ) {
      ((CProxy_WorkDistrib(thisgroup))[CkMyPe()+1]).savePatchMap(msg);
    } else {
      delete msg;
    }
    return;
  }

  patchMapArrived = true;

  int self = CkMyNode();
  int range_begin = 0;
  int range_end = CkNumNodes();
  while ( self != range_begin ) {
    ++range_begin;
    int split = range_begin + ( range_end - range_begin ) / 2;
    if ( self < split ) { range_end = split; }
    else { range_begin = split; }
  }
  int send_near = self + 1;
  int send_far = send_near + ( range_end - send_near ) / 2;

  int pids[3];
  int npid = 0;
  if ( send_far < range_end ) pids[npid++] = CkNodeFirst(send_far);
  if ( send_near < send_far ) pids[npid++] = CkNodeFirst(send_near);
  pids[npid++] = CkMyPe();  // always send the message to ourselves
  CProxy_WorkDistrib(thisgroup).savePatchMap(msg,npid,pids);
}


class ComputeMapMsg : public CMessage_ComputeMapMsg {
  public:
    int nComputes;
    ComputeMap::ComputeData *computeMapData;
};

void WorkDistrib::sendComputeMap(void)
{
  if ( CkNumNodes() == 1 ) {
    computeMapArrived = true;
    ComputeMap::Object()->initPtrs();
    return;
  }

  int size = ComputeMap::Object()->numComputes();

  ComputeMapMsg *mapMsg = new (size, 0) ComputeMapMsg;

  mapMsg->nComputes = size;
  ComputeMap::Object()->pack(mapMsg->computeMapData);

  CProxy_WorkDistrib workProxy(thisgroup);
  workProxy[0].saveComputeMap(mapMsg);
}

// saveMaps() is called when the map message is received
void WorkDistrib::saveComputeMap(ComputeMapMsg *msg)
{
  // Use a resend to forward messages before processing.  Otherwise the
  // map distribution is slow on many CPUs.  We need to use a tree
  // rather than a broadcast because some implementations of broadcast
  // generate a copy of the message on the sender for each recipient.
  // This is because MPI doesn't allow re-use of an outstanding buffer.

  if ( CkMyRank() ) {
    NAMD_bug("WorkDistrib::saveComputeMap called on non-rank-zero pe");
  }

  if ( computeMapArrived && CkMyPe() ) {
    ComputeMap::Object()->unpack(msg->nComputes, msg->computeMapData);
  }

  if ( computeMapArrived ) {
    delete msg;
    ComputeMap::Object()->initPtrs();
    return;
  }

  computeMapArrived = true;

  int self = CkMyNode();
  int range_begin = 0;
  int range_end = CkNumNodes();
  while ( self != range_begin ) {
    ++range_begin;
    int split = range_begin + ( range_end - range_begin ) / 2;
    if ( self < split ) { range_end = split; }
    else { range_begin = split; }
  }
  int send_near = self + 1;
  int send_far = send_near + ( range_end - send_near ) / 2;

  int pids[3];
  int npid = 0;
  if ( send_far < range_end ) pids[npid++] = CkNodeFirst(send_far);
  if ( send_near < send_far ) pids[npid++] = CkNodeFirst(send_near);
  pids[npid++] = CkMyPe();  // always send the message to ourselves
  CProxy_WorkDistrib(thisgroup).saveComputeMap(msg,npid,pids);
}


//----------------------------------------------------------------------
void WorkDistrib::patchMapInit(void)
{
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  Lattice lattice = params->lattice;

  BigReal patchSize = params->patchDimension;

#ifndef MEM_OPT_VERSION
  const int totalAtoms = node->pdb->num_atoms();
#else
  const int totalAtoms = node->molecule->numAtoms;
#endif

  ScaledPosition xmin, xmax;

  double maxNumPatches = 1.e9;  // need to adjust fractional values
  if ( params->minAtomsPerPatch > 0 )
    maxNumPatches = totalAtoms / params->minAtomsPerPatch;

  DebugM(3,"Mapping patches\n");
  if ( lattice.a_p() && lattice.b_p() && lattice.c_p() ) {
    xmin = 0.;  xmax = 0.;
  }
  else if ( params->FMAOn || params->MSMOn || params->FMMOn ) {
  // Need to use full box for FMA to match NAMD 1.X results.
#if 0
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r());
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r());
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r());
#endif
    node->pdb->find_extremes(lattice);
    node->pdb->get_extremes(xmin, xmax);
#if 0
    printf("+++ center=%.4f %.4f %.4f\n",
        lattice.origin().x, lattice.origin().y, lattice.origin().z);
    printf("+++ xmin=%.4f  xmax=%.4f\n", xmin.x, xmax.x);
    printf("+++ ymin=%.4f  ymax=%.4f\n", xmin.y, xmax.y);
    printf("+++ zmin=%.4f  zmax=%.4f\n", xmin.z, xmax.z);
#endif
  // Otherwise, this allows a small number of stray atoms.
  }
  else {
#if 0
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r(),0.9);
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r(),0.9);
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r(),0.9);
#endif
    node->pdb->find_extremes(lattice, 1.0);
    node->pdb->get_extremes(xmin, xmax);
    iout << iINFO << "ORIGINAL ATOMS MINMAX IS " << xmin << "  " << xmax << "\n" << endi;
    double frac = ( (double)totalAtoms - 10000. ) / (double)totalAtoms;
    if ( frac < 0.9 ) { frac = 0.9; }
    node->pdb->find_extremes(lattice, frac);
    node->pdb->get_extremes(xmin, xmax);
    iout << iINFO << "ADJUSTED ATOMS MINMAX IS " << xmin << "  " << xmax << "\n" << endi;
  }

#if 0
  BigReal origin_shift;
  origin_shift = lattice.a_r() * lattice.origin();
  xmin.x -= origin_shift;
  xmax.x -= origin_shift;
  origin_shift = lattice.b_r() * lattice.origin();
  xmin.y -= origin_shift;
  xmax.y -= origin_shift;
  origin_shift = lattice.c_r() * lattice.origin();
  xmin.z -= origin_shift;
  xmax.z -= origin_shift;
#endif

  // SimParameters default is -1 for unset
  int twoAwayX = params->twoAwayX;
  int twoAwayY = params->twoAwayY;
  int twoAwayZ = params->twoAwayZ;

  // SASA implementation is not compatible with twoAway patches
  if (params->LCPOOn && patchSize < 32.4) {
    if ( twoAwayX > 0 || twoAwayY > 0 || twoAwayZ > 0 ) {
      iout << iWARN << "Ignoring twoAway[XYZ] due to LCPO SASA implementation.\n" << endi;
    }
    twoAwayX = twoAwayY = twoAwayZ = 0;
  }

  // if you think you know what you're doing go right ahead
  if ( twoAwayX > 0 ) maxNumPatches = 1.e9;
  if ( twoAwayY > 0 ) maxNumPatches = 1.e9;
  if ( twoAwayZ > 0 ) maxNumPatches = 1.e9;
  if ( params->maxPatches > 0 ) {
      maxNumPatches = params->maxPatches;
      iout << iINFO << "LIMITING NUMBER OF PATCHES TO " <<
                                maxNumPatches << "\n" << endi;
  }

  int numpes = CkNumPes();
  SimParameters *simparam = Node::Object()->simParameters;
  if(simparam->simulateInitialMapping) {
    numpes = simparam->simulatedPEs;
    delete [] patchMap->nPatchesOnNode;
    patchMap->nPatchesOnNode = new int[numpes];
    memset(patchMap->nPatchesOnNode, 0, numpes*sizeof(int));	
  }

#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  // for CUDA be sure there are more patches than pes

  int numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  if ( numPatches < numpes && twoAwayX < 0 ) {
    twoAwayX = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < numpes && twoAwayY < 0 ) {
    twoAwayY = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < numpes && twoAwayZ < 0 ) {
    twoAwayZ = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < numpes ) {
    #if defined(NAMD_MIC)
    NAMD_die("MIC-enabled NAMD requires at least one patch per thread.");
    #endif
  }
  if ( numPatches % numpes && numPatches <= 1.4 * numpes ) {
    int exactFit = numPatches - numPatches % numpes;
    int newNumPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,exactFit,params->staticAtomAssignment,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
    if ( newNumPatches == exactFit ) {
      iout << iINFO << "REDUCING NUMBER OF PATCHES TO IMPROVE LOAD BALANCE\n" << endi;
      maxNumPatches = exactFit;
    }
  }

  patchMap->makePatches(xmin,xmax,lattice,patchSize,maxNumPatches,
	params->staticAtomAssignment, params->replicaUniformPatchGrids, params->LCPOOn,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);

#else

  int availPes = numpes;
  if ( params->noPatchesOnZero && numpes > 1 ) {
      availPes -= 1;
      if(params->noPatchesOnOne && numpes > 2)
        availPes -= 1;
  }
#ifdef MEM_OPT_VERSION
  if(params->noPatchesOnOutputPEs && numpes - params->numoutputprocs >2)
    {
      availPes -= params->numoutputprocs;
      if ( params->noPatchesOnZero && numpes > 1 && isOutputProcessor(0)){
	availPes++;
      }
      if ( params->noPatchesOnOne && numpes > 2 && isOutputProcessor(1)){
	availPes++;
      }
    }
#endif

  int numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  if ( ( numPatches > (0.3*availPes) || numPatches > maxNumPatches
       ) && twoAwayZ < 0 ) {
    twoAwayZ = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( ( numPatches > (0.6*availPes) || numPatches > maxNumPatches
       ) && twoAwayY < 0 ) {
    twoAwayY = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( ( numPatches > availPes || numPatches > maxNumPatches
       ) && twoAwayX < 0 ) {
    twoAwayX = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,1.e9,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( numPatches > availPes && numPatches <= (1.4*availPes) && availPes <= maxNumPatches ) {
    int newNumPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,availPes,params->staticAtomAssignment,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
    if ( newNumPatches <= availPes && numPatches <= (1.4*newNumPatches) ) {
      iout << iINFO << "REDUCING NUMBER OF PATCHES TO IMPROVE LOAD BALANCE\n" << endi;
      maxNumPatches = availPes;
    }
  }

  patchMap->makePatches(xmin,xmax,lattice,patchSize,maxNumPatches,
	params->staticAtomAssignment, params->replicaUniformPatchGrids, params->LCPOOn,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);

#endif

}


//----------------------------------------------------------------------
void WorkDistrib::assignNodeToPatch()
{
  PatchMap *patchMap = PatchMap::Object();
  int nNodes = Node::Object()->numNodes();
  SimParameters *simparam = Node::Object()->simParameters;
  if(simparam->simulateInitialMapping) {
	  nNodes = simparam->simulatedPEs;
  }

#if (CMK_BLUEGENEP | CMK_BLUEGENEL) && USE_TOPOMAP 
  TopoManager tmgr;
  int numPes = tmgr.getDimNX() * tmgr.getDimNY() * tmgr.getDimNZ();
  if (numPes > patchMap->numPatches() && (assignPatchesTopoGridRecBisection() > 0)) {
    CkPrintf ("Blue Gene/L topology partitioner finished successfully \n");
  }
  else
#endif
  assignPatchesSpaceFillingCurve();	  
  
  int *nAtoms = new int[nNodes];
  int numAtoms=0;
  int i;
  for(i=0; i < nNodes; i++)
    nAtoms[i] = 0;

  for(i=0; i < patchMap->numPatches(); i++)
  {
    //    iout << iINFO << "Patch " << i << " has " 
    //	 << patchMap->patch(i)->getNumAtoms() << " atoms and "
    //	 << patchMap->patch(i)->getNumAtoms() * 
    //            patchMap->patch(i)->getNumAtoms() 
    //	 << " pairs.\n" << endi;	 
#ifdef MEM_OPT_VERSION
      numAtoms += patchMap->numAtoms(i);
      nAtoms[patchMap->node(i)] += patchMap->numAtoms(i);	  
#else
    if (patchMap->patch(i)) {
      numAtoms += patchMap->patch(i)->getNumAtoms();
      nAtoms[patchMap->node(i)] += patchMap->patch(i)->getNumAtoms();	  
    }
#endif
  }

  if ( numAtoms != Node::Object()->molecule->numAtoms ) {
    for(i=0; i < nNodes; i++)
      iout << iINFO << nAtoms[i] << " atoms assigned to node " << i << "\n" << endi;
    iout << iINFO << "Assigned " << numAtoms << " atoms but expected " << Node::Object()->molecule->numAtoms << "\n" << endi;
    NAMD_die("Incorrect atom count in WorkDistrib::assignNodeToPatch\n");
  }

  delete [] nAtoms;
 
  //  PatchMap::Object()->printPatchMap();
}

//----------------------------------------------------------------------
// void WorkDistrib::assignPatchesSlices() 
// {
//   int pid; 
//   int assignedNode = 0;
//   PatchMap *patchMap = PatchMap::Object();
//   Node *node = CLocalBranch(Node, CkpvAccess(BOCclass_group).node);

//   int *numAtoms = new int[node->numNodes()];
//   for (int i=0; i<node->numNodes(); i++) {
//     numAtoms[i] = 0;
//   }

//   // Assign patch to node with least atoms assigned.
//   for(pid=0; pid < patchMap->numPatches(); pid++) {
//     assignedNode = 0;
//     for (i=1; i < node->numNodes(); i++) {
//       if (numAtoms[i] < numAtoms[assignedNode]) assignedNode = i;
//     }
//     patchMap->assignNode(pid, assignedNode);
//     numAtoms[assignedNode] += patchMap->patch(pid)->getNumAtoms();

//     /*
//     iout << iINFO << "Patch (" << pid << ") has " 
//       << patchMap->patch(pid)->getNumAtoms() 
//       << " atoms:  Assigned to Node(" << assignedNode << ")\n" 
//       << endi;
//     */
//   }

//   delete[] numAtoms;
// }

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesToLowestLoadNode() 
{
  int pid; 
  int assignedNode = 0;
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;
  int ncpus = node->numNodes();
  if(simParams->simulateInitialMapping) {
	  ncpus = simParams->simulatedPEs;
  }

  int *load = new int[ncpus];
  int *assignedNodes = new int[patchMap->numPatches()];
  for (int i=0; i<ncpus; i++) {
    load[i] = 0;
  }
  CkPrintf("assignPatchesToLowestLoadNode\n");
  int defaultNode = 0;
  if ( simParams->noPatchesOnZero && ncpus > 1 ){
    defaultNode = 1;
    if( simParams->noPatchesOnOne && ncpus > 2)
      defaultNode = 2;
  }
  // Assign patch to node with least atoms assigned.
  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode = defaultNode;
    for (int i=assignedNode + 1; i < ncpus; i++) {
      if (load[i] < load[assignedNode]) assignedNode = i;
    }
    assignedNodes[pid] = assignedNode;
#ifdef MEM_OPT_VERSION
    load[assignedNode] += patchMap->numAtoms(pid) + 1;
#else
    load[assignedNode] += patchMap->patch(pid)->getNumAtoms() + 1;
#endif
  }

  delete[] load;
  sortNodesAndAssign(assignedNodes);
  delete[] assignedNodes;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesBitReversal() 
{
  int pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simparam = node->simParameters;

  int ncpus = node->numNodes();
  if(simparam->simulateInitialMapping) {
	  ncpus = simparam->simulatedPEs;
  }
  int npatches = patchMap->numPatches();
  if ( ncpus <= npatches )
    NAMD_bug("WorkDistrib::assignPatchesBitReversal called improperly");

  SortableResizeArray<int> seq(ncpus);
  // avoid using node 0 (reverse of 0 is 0 so start at 1)
  for ( int i = 1; i < ncpus; ++i ) {
    seq[i-1] = peDiffuseOrdering[i];
  }

  // extract and sort patch locations
  sortNodesAndAssign(seq.begin());
  if ( ncpus > 2*npatches ) sortNodesAndAssign(seq.begin()+npatches, 1);
}

//----------------------------------------------------------------------
struct nodesort {
  int node;
  int a_total;
  int b_total;
  int c_total;
  int npatches;
  nodesort() : node(-1),a_total(0),b_total(0),c_total(0),npatches(0) { ; }
  int operator==(const nodesort &o) const {
    float a1 = ((float)a_total)/((float)npatches);
    float a2 = ((float)o.a_total)/((float)o.npatches);
    float b1 = ((float)b_total)/((float)npatches);
    float b2 = ((float)o.b_total)/((float)o.npatches);
    float c1 = ((float)c_total)/((float)npatches);
    float c2 = ((float)o.c_total)/((float)o.npatches);
    return ((a1 == a2) && (b1 == b2) && (c1 == c2));
  }
  int operator<(const nodesort &o) const {
    float a1 = ((float)a_total)/((float)npatches);
    float a2 = ((float)o.a_total)/((float)o.npatches);
    float b1 = ((float)b_total)/((float)npatches);
    float b2 = ((float)o.b_total)/((float)o.npatches);
    float c1 = ((float)c_total)/((float)npatches);
    float c2 = ((float)o.c_total)/((float)o.npatches);
    return ( (a1 < a2) || ((a1 == a2) && (b1 < b2)) ||
		((a1 == a2) && (b1 == b2) && (c1 < c2)) );
  }
};

void WorkDistrib::sortNodesAndAssign(int *assignedNode, int baseNodes) {
  // if baseNodes is zero (default) then set both nodes and basenodes
  // if baseNodes is nonzero then this is a second call to set basenodes only
  int i, pid; 
  PatchMap *patchMap = PatchMap::Object();
  int npatches = patchMap->numPatches();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  int nnodes = node->numNodes();
  SimParameters *simparam = node->simParameters;
  if(simparam->simulateInitialMapping) {
	  nnodes = simparam->simulatedPEs;
  }

  ResizeArray<nodesort> allnodes(nnodes);
  for ( i=0; i < nnodes; ++i ) {
    allnodes[i].node = i;
  }
  for ( pid=0; pid<npatches; ++pid ) {
    // iout << pid << " " << assignedNode[pid] << "\n" << endi;
    allnodes[assignedNode[pid]].npatches++;
    allnodes[assignedNode[pid]].a_total += patchMap->index_a(pid);
    allnodes[assignedNode[pid]].b_total += patchMap->index_b(pid);
    allnodes[assignedNode[pid]].c_total += patchMap->index_c(pid);
  }
  SortableResizeArray<nodesort> usednodes(nnodes);
  usednodes.resize(0);
  for ( i=0; i < nnodes; ++i ) {
    if ( allnodes[i].npatches ) usednodes.add(allnodes[i]);
  }
  usednodes.sort();
  int i2 = 0;
  for ( i=0; i < nnodes; ++i ) {
    int pe = peCompactOrdering[i];
    if ( allnodes[pe].npatches ) allnodes[usednodes[i2++].node].node = pe;
  }

  for ( pid=0; pid<npatches; ++pid ) {
    // iout << pid << " " <<  allnodes[assignedNode[pid]].node << "\n" << endi;
    if ( ! baseNodes ) {
      patchMap->assignNode(pid, allnodes[assignedNode[pid]].node);	
    }
    patchMap->assignBaseNode(pid, allnodes[assignedNode[pid]].node);	
  }
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRoundRobin() 
{
  int pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simparam = node->simParameters;
  int ncpus = node->numNodes();
  if(simparam->simulateInitialMapping) {
	  ncpus = simparam->simulatedPEs;
  }
  int *assignedNode = new int[patchMap->numPatches()];

  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode[pid] = pid % ncpus;
  }

  sortNodesAndAssign(assignedNode);
  delete [] assignedNode;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRecursiveBisection() 
{
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  SimParameters *simParams = Node::Object()->simParameters;
  int numNodes = Node::Object()->numNodes();
  if(simParams->simulateInitialMapping) {
	  numNodes = simParams->simulatedPEs;
  }
  
  int usedNodes = numNodes;
  int unusedNodes = 0;
  CkPrintf("assignPatchesRecursiveBisection\n");
  if ( simParams->noPatchesOnZero && numNodes > 1 ){
    usedNodes -= 1;
    if(simParams->noPatchesOnOne && numNodes > 2)
      usedNodes -= 1;
  }
  unusedNodes = numNodes - usedNodes;
  RecBisection recBisec(usedNodes,PatchMap::Object());
  if ( recBisec.partition(assignedNode) ) {
    if ( unusedNodes !=0 ) {
      for ( int i=0; i<patchMap->numPatches(); ++i ) {
        assignedNode[i] += unusedNodes;
      }
    }
    sortNodesAndAssign(assignedNode);
    delete [] assignedNode;
  } else {
    //free the array here since a same array will be allocated
    //in assignPatchesToLowestLoadNode function, thus reducting
    //temporary memory usage
    delete [] assignedNode; 
    
    iout << iWARN 
	 << "WorkDistrib: Recursive bisection fails, "
	 << "invoking space-filling curve algorithm\n";
    assignPatchesSpaceFillingCurve();
  }
}

// class to re-order dimensions in decreasing size
struct TopoManagerWrapper {
  TopoManager tmgr;
  int a_dim, b_dim, c_dim, d_dim, e_dim;
  int a_rot, b_rot, c_rot, d_rot, e_rot;
  int a_mod, b_mod, c_mod, d_mod, e_mod;
  int fixpe(int pe) {  // compensate for lame fallback topology information
    return CmiGetFirstPeOnPhysicalNode(CmiPhysicalNodeID(pe));
  }
  TopoManagerWrapper() {
#if CMK_BLUEGENEQ
    int na=tmgr.getDimNA();
    int nb=tmgr.getDimNB();
    int nc=tmgr.getDimNC();
    int nd=tmgr.getDimND();
    int ne=tmgr.getDimNE();
#else
    int na=tmgr.getDimNX();
    int nb=tmgr.getDimNY();
    int nc=tmgr.getDimNZ();
    int nd=1;
    int ne=1;
#endif
    ResizeArray<int> a_flags(na);
    ResizeArray<int> b_flags(nb);
    ResizeArray<int> c_flags(nc);
    ResizeArray<int> d_flags(nd);
    ResizeArray<int> e_flags(ne);
    for ( int i=0; i<na; ++i ) { a_flags[i] = 0; }
    for ( int i=0; i<nb; ++i ) { b_flags[i] = 0; }
    for ( int i=0; i<nc; ++i ) { c_flags[i] = 0; }
    for ( int i=0; i<nd; ++i ) { d_flags[i] = 0; }
    for ( int i=0; i<ne; ++i ) { e_flags[i] = 0; }
    int npes = CkNumPes();
    for ( int pe=0; pe<npes; ++pe ) {
      int a,b,c,d,e,t;
#if CMK_BLUEGENEQ
      tmgr.rankToCoordinates(fixpe(pe),a,b,c,d,e,t);
#else
      tmgr.rankToCoordinates(fixpe(pe),a,b,c,t);
      d=0; e=0;
#endif
      if ( a < 0 || a >= na ) NAMD_bug("inconsistent torus topology!");
      if ( b < 0 || b >= nb ) NAMD_bug("inconsistent torus topology!");
      if ( c < 0 || c >= nc ) NAMD_bug("inconsistent torus topology!");
      if ( d < 0 || d >= nd ) NAMD_bug("inconsistent torus topology!");
      if ( e < 0 || e >= ne ) NAMD_bug("inconsistent torus topology!");
      a_flags[a] = 1;
      b_flags[b] = 1;
      c_flags[c] = 1;
      d_flags[d] = 1;
      e_flags[e] = 1;
    }
    iout << iINFO << "TORUS A SIZE " << na << " USING";
    for ( int i=0; i<na; ++i ) { if ( a_flags[i] ) iout << " " << i; }
    iout << "\n" << endi;
    iout << iINFO << "TORUS B SIZE " << nb << " USING";
    for ( int i=0; i<nb; ++i ) { if ( b_flags[i] ) iout << " " << i; }
    iout << "\n" << endi;
    iout << iINFO << "TORUS C SIZE " << nc << " USING";
    for ( int i=0; i<nc; ++i ) { if ( c_flags[i] ) iout << " " << i; }
    iout << "\n" << endi;
#if CMK_BLUEGENEQ
    iout << iINFO << "TORUS D SIZE " << nd << " USING";
    for ( int i=0; i<nd; ++i ) { if ( d_flags[i] ) iout << " " << i; }
    iout << "\n" << endi;
    iout << iINFO << "TORUS E SIZE " << ne << " USING";
    for ( int i=0; i<ne; ++i ) { if ( e_flags[i] ) iout << " " << i; }
    iout << "\n" << endi;
#endif
    // find most compact representation of our subset
    a_rot = b_rot = c_rot = d_rot = e_rot = 0;
    a_mod = na; b_mod = nb; c_mod = nc; d_mod = nd; e_mod = ne;
#if CMK_BLUEGENEQ
    if ( tmgr.absA(na) == 0 ) // torus
#else
    if ( tmgr.absX(na) == 0 ) // torus
#endif
      for ( int i=0, gaplen=0, gapstart=0; i<2*na; ++i ) {
        if ( a_flags[i%na] ) gapstart = i+1;
        else if ( i - gapstart >= gaplen ) {
          a_rot = 2*na-i-1; gaplen = i - gapstart;
        }
      }
#if CMK_BLUEGENEQ
    if ( tmgr.absB(nb) == 0 ) // torus
#else
    if ( tmgr.absY(nb) == 0 ) // torus
#endif
      for ( int i=0, gaplen=0, gapstart=0; i<2*nb; ++i ) {
        if ( b_flags[i%nb] ) gapstart = i+1;
        else if ( i - gapstart >= gaplen ) {
          b_rot = 2*nb-i-1; gaplen = i - gapstart;
        }
      }
#if CMK_BLUEGENEQ
    if ( tmgr.absC(nc) == 0 ) // torus
#else
    if ( tmgr.absZ(nc) == 0 ) // torus
#endif
      for ( int i=0, gaplen=0, gapstart=0; i<2*nc; ++i ) {
        if ( c_flags[i%nc] ) gapstart = i+1;
        else if ( i - gapstart >= gaplen ) {
          c_rot = 2*nc-i-1; gaplen = i - gapstart;
        }
      }
#if CMK_BLUEGENEQ
    if ( tmgr.absD(nd) == 0 ) // torus
      for ( int i=0, gaplen=0, gapstart=0; i<2*nd; ++i ) {
        if ( d_flags[i%nd] ) gapstart = i+1;
        else if ( i - gapstart >= gaplen ) {
          d_rot = 2*nd-i-1; gaplen = i - gapstart;
        }
      }
    if ( tmgr.absE(ne) == 0 ) // torus
      for ( int i=0, gaplen=0, gapstart=0; i<2*ne; ++i ) {
        if ( e_flags[i%ne] ) gapstart = i+1;
        else if ( i - gapstart >= gaplen ) {
          e_rot = 2*ne-i-1; gaplen = i - gapstart;
        }
      }
#endif
    // order dimensions by length
    int a_min=na, a_max=-1;
    int b_min=nb, b_max=-1;
    int c_min=nc, c_max=-1;
    int d_min=nd, d_max=-1;
    int e_min=ne, e_max=-1;
    for ( int pe=0; pe<npes; ++pe ) {
      int a,b,c,d,e,t;
#if CMK_BLUEGENEQ
      tmgr.rankToCoordinates(fixpe(pe),a,b,c,d,e,t);
#else
      tmgr.rankToCoordinates(fixpe(pe),a,b,c,t);
      d=0; e=0;
#endif
      a = (a+a_rot)%a_mod;
      b = (b+b_rot)%b_mod;
      c = (c+c_rot)%c_mod;
      d = (d+d_rot)%d_mod;
      e = (e+e_rot)%e_mod;
      if ( a < a_min ) a_min = a;
      if ( b < b_min ) b_min = b;
      if ( c < c_min ) c_min = c;
      if ( d < d_min ) d_min = d;
      if ( e < e_min ) e_min = e;
      if ( a > a_max ) a_max = a;
      if ( b > b_max ) b_max = b;
      if ( c > c_max ) c_max = c;
      if ( d > d_max ) d_max = d;
      if ( e > e_max ) e_max = e;
    }
    int a_len = a_max - a_min + 1;
    int b_len = b_max - b_min + 1;
    int c_len = c_max - c_min + 1;
    int d_len = d_max - d_min + 1;
    int e_len = e_max - e_min + 1;
    int lensort[5];
    lensort[0] = (a_len << 3) + 0;
    lensort[1] = (b_len << 3) + 1;
    lensort[2] = (c_len << 3) + 2;
    lensort[3] = (d_len << 3) + 3;
    lensort[4] = (e_len << 3) + 4;
    // CkPrintf("TopoManagerWrapper lensort before %d %d %d %d %d\n", lensort[0] & 7, lensort[1] & 7, lensort[2] & 7, lensort[3] & 7, lensort[4] & 7);
    std::sort(lensort, lensort+5);
    // CkPrintf("TopoManagerWrapper lensort after %d %d %d %d %d\n", lensort[0] & 7, lensort[1] & 7, lensort[2] & 7, lensort[3] & 7, lensort[4] & 7);
    for ( int i=0; i<5; ++i ) { if ( (lensort[i] & 7) == 0 ) a_dim = 4-i; }
    for ( int i=0; i<5; ++i ) { if ( (lensort[i] & 7) == 1 ) b_dim = 4-i; }
    for ( int i=0; i<5; ++i ) { if ( (lensort[i] & 7) == 2 ) c_dim = 4-i; }
    for ( int i=0; i<5; ++i ) { if ( (lensort[i] & 7) == 3 ) d_dim = 4-i; }
    for ( int i=0; i<5; ++i ) { if ( (lensort[i] & 7) == 4 ) e_dim = 4-i; }
#if 0
    if ( a_len >= b_len && a_len >= c_len ) {
      a_dim = 0;
      if ( b_len >= c_len ) {
        b_dim = 1; c_dim = 2;
      } else {
        b_dim = 2; c_dim = 1;
      }
    } else if ( b_len >= a_len && b_len >= c_len ) {
      b_dim = 0;
      if ( a_len >= c_len ) {
        a_dim = 1; c_dim = 2;
      } else {
        a_dim = 2; c_dim = 1;
      }
    } else { // c is longest
      c_dim = 0;
      if ( a_len >= b_len ) {
        a_dim = 1; b_dim = 2;
      } else {
        a_dim = 2; b_dim = 1;
      }
    }
#endif
    iout << iINFO << "TORUS MINIMAL MESH SIZE IS " << a_len << " BY " << b_len << " BY " << c_len
#if CMK_BLUEGENEQ
    << " BY " << d_len << " BY " << e_len
#endif
    << "\n" << endi;
    // CkPrintf("TopoManagerWrapper dims %d %d %d %d %d\n", a_dim, b_dim, c_dim, d_dim, e_dim);
  }
  void coords(int pe, int *crds) {
    int a,b,c,d,e,t;
#if CMK_BLUEGENEQ
    tmgr.rankToCoordinates(fixpe(pe),a,b,c,d,e,t);
#else
    tmgr.rankToCoordinates(fixpe(pe),a,b,c,t);
    d=0; e=0;
#endif
    if ( a_dim < 3 ) crds[a_dim] = (a+a_rot)%a_mod;
    if ( b_dim < 3 ) crds[b_dim] = (b+b_rot)%b_mod;
    if ( c_dim < 3 ) crds[c_dim] = (c+c_rot)%c_mod;
    if ( d_dim < 3 ) crds[d_dim] = (d+d_rot)%d_mod;
    if ( e_dim < 3 ) crds[e_dim] = (e+e_rot)%e_mod;
  }
  int coord(int pe, int dim) {
    int crds[3];
    coords(pe,crds);
    return crds[dim];
  }
  struct pe_sortop_topo {
    TopoManagerWrapper &tmgr;
    const int *sortdims;
    pe_sortop_topo(TopoManagerWrapper &t, int *d) : tmgr(t), sortdims(d) {}
    bool operator() (int pe1, int pe2) const {
      int crds1[3], crds2[3];
      tmgr.coords(pe1,crds1);
      tmgr.coords(pe2,crds2);
      for ( int i=0; i<3; ++i ) {
        int d = sortdims[i];
        if ( crds1[d] != crds2[d] ) return ( crds1[d] < crds2[d] );
      }
      const int *index = WorkDistrib::peCompactOrderingIndex;
      return ( index[pe1] < index[pe2] );
    }
  };
  int* sortAndSplit(int *node_begin, int *node_end, int splitdim) {
    if ( node_begin == node_end ) return node_begin;
    int tmins[3], tmaxs[3], tlens[3], sortdims[3];
    coords(*node_begin, tmins);
    coords(*node_begin, tmaxs);
    for ( int *peitr = node_begin; peitr != node_end; ++peitr ) {
      int tvals[3];
      coords(*peitr, tvals);
      for ( int i=0; i<3; ++i ) {
        if ( tvals[i] < tmins[i] ) tmins[i] = tvals[i];
        if ( tvals[i] > tmaxs[i] ) tmaxs[i] = tvals[i];
      }
    }
    for ( int i=0; i<3; ++i ) {
      tlens[i] = tmaxs[i] - tmins[i];
    }
    sortdims[0] = splitdim;
    for ( int i=0, j=0; i<3; ++i ) {
      if ( i != splitdim ) sortdims[++j] = i;
    }
    if ( tlens[sortdims[1]] < tlens[sortdims[2]] ) {
      int tmp = sortdims[1];
      sortdims[1] = sortdims[2];
      sortdims[2] = tmp;
    }
    std::sort(node_begin,node_end,pe_sortop_topo(*this,sortdims));
    int *nodes = node_begin;
    int nnodes = node_end - node_begin;
    int i_split = 0;
#if 0
    int c_split = coord(nodes[0],splitdim);
    for ( int i=0; i<nnodes; ++i ) {
      if ( coord(nodes[i],splitdim) != c_split ) {
        int mid = (nnodes+1)/2;
        if ( abs(i-mid) < abs(i_split-mid) ) {
          i_split = i;
          c_split = coord(i,splitdim);
        }
        else break;
      }
    }
#endif
    for ( int i=0; i<nnodes; ++i ) {
      if ( ! CmiPeOnSamePhysicalNode(nodes[i_split],nodes[i]) ) {
        int mid = (nnodes+1)/2;
        if ( abs(i-mid) < abs(i_split-mid) ) i_split = i;
        else break;
      }
    }
    return ( node_begin + i_split );
  }
};

struct patch_sortop_curve_a {
  PatchMap *pmap;
  patch_sortop_curve_a(PatchMap *m) : pmap(m) {}
  inline bool operator() (int p1, int p2) const {
    int a1 = pmap->index_a(p1);
    int a2 = pmap->index_a(p2);
    if ( a1 < a2 ) return true;
    if ( a1 > a2 ) return false;
    int dir = ( (a1 & 1) ? -1 : 1 );
    int b1 = pmap->index_b(p1);
    int b2 = pmap->index_b(p2);
    if ( b1 * dir < b2 * dir ) return true;
    if ( b1 * dir > b2 * dir ) return false;
    dir *= ( (b1 & 1) ? -1 : 1 );
    int c1 = pmap->index_c(p1);
    int c2 = pmap->index_c(p2);
    if ( c1 * dir < c2 * dir ) return true;
    return false;
  }
};

struct patch_sortop_curve_b {
  PatchMap *pmap;
  patch_sortop_curve_b(PatchMap *m) : pmap(m) {}
  inline bool operator() (int p1, int p2) const {
    int a1 = pmap->index_b(p1);
    int a2 = pmap->index_b(p2);
    if ( a1 < a2 ) return true;
    if ( a1 > a2 ) return false;
    int dir = ( (a1 & 1) ? -1 : 1 );
    int b1 = pmap->index_a(p1);
    int b2 = pmap->index_a(p2);
    if ( b1 * dir < b2 * dir ) return true;
    if ( b1 * dir > b2 * dir ) return false;
    dir *= ( (b1 & 1) ? -1 : 1 );
    int c1 = pmap->index_c(p1);
    int c2 = pmap->index_c(p2);
    if ( c1 * dir < c2 * dir ) return true;
    return false;
  }
};

struct patch_sortop_curve_c {
  PatchMap *pmap;
  patch_sortop_curve_c(PatchMap *m) : pmap(m) {}
  inline bool operator() (int p1, int p2) const {
    int a1 = pmap->index_c(p1);
    int a2 = pmap->index_c(p2);
    if ( a1 < a2 ) return true;
    if ( a1 > a2 ) return false;
    int dir = ( (a1 & 1) ? -1 : 1 );
    int b1 = pmap->index_a(p1);
    int b2 = pmap->index_a(p2);
    if ( b1 * dir < b2 * dir ) return true;
    if ( b1 * dir > b2 * dir ) return false;
    dir *= ( (b1 & 1) ? -1 : 1 );
    int c1 = pmap->index_b(p1);
    int c2 = pmap->index_b(p2);
    if ( c1 * dir < c2 * dir ) return true;
    return false;
  }
};

static void recursive_bisect_with_curve(
  int *patch_begin, int *patch_end,
  int *node_begin, int *node_end,
  double *patchLoads,
  double *sortedLoads,
  int *assignedNode,
  TopoManagerWrapper &tmgr
  ) {

  SimParameters *simParams = Node::Object()->simParameters;
  PatchMap *patchMap = PatchMap::Object();
  int *patches = patch_begin;
  int npatches = patch_end - patch_begin;
  int *nodes = node_begin;
  int nnodes = node_end - node_begin;

  // assign patch loads
  double totalRawLoad = 0;
  for ( int i=0; i<npatches; ++i ) {
    int pid=patches[i];
#ifdef MEM_OPT_VERSION
    double load = patchMap->numAtoms(pid) + 10;      
#else
    double load = patchMap->patch(pid)->getNumAtoms() + 10;
#endif
    patchLoads[pid] = load;
    sortedLoads[i] = load;
    totalRawLoad += load;
  }
  std::sort(sortedLoads,sortedLoads+npatches);

  // limit maxPatchLoad to adjusted average load per node
  double sumLoad = 0;
  double maxPatchLoad = 1;
  for ( int i=0; i<npatches; ++i ) {
    double load = sortedLoads[i];
    double total = sumLoad + (npatches-i) * load;
    if ( nnodes * load > total ) break;
    sumLoad += load;
    maxPatchLoad = load;
  }
  double totalLoad = 0;
  for ( int i=0; i<npatches; ++i ) {
    int pid=patches[i];
    if ( patchLoads[pid] > maxPatchLoad ) patchLoads[pid] = maxPatchLoad;
    totalLoad += patchLoads[pid];
  }
  if ( nnodes * maxPatchLoad > totalLoad )
    NAMD_bug("algorithm failure in WorkDistrib recursive_bisect_with_curve()");

  int a_len, b_len, c_len;
  int a_min, b_min, c_min;
  { // find dimensions
    a_min = patchMap->index_a(patches[0]);
    b_min = patchMap->index_b(patches[0]);
    c_min = patchMap->index_c(patches[0]);
    int a_max = a_min;
    int b_max = b_min;
    int c_max = c_min;
    for ( int i=1; i<npatches; ++i ) {
      int a = patchMap->index_a(patches[i]);
      int b = patchMap->index_b(patches[i]);
      int c = patchMap->index_c(patches[i]);
      if ( a < a_min ) a_min = a;
      if ( b < b_min ) b_min = b;
      if ( c < c_min ) c_min = c;
      if ( a > a_max ) a_max = a;
      if ( b > b_max ) b_max = b;
      if ( c > c_max ) c_max = c;
    }
    a_len = a_max - a_min;
    b_len = b_max - b_min;
    c_len = c_max - c_min;
  }

  int *node_split = node_begin;

  if ( simParams->disableTopology ) ; else
  if ( a_len >= b_len && a_len >= c_len ) {
    node_split = tmgr.sortAndSplit(node_begin,node_end,0);
  } else if ( b_len >= a_len && b_len >= c_len ) {
    node_split = tmgr.sortAndSplit(node_begin,node_end,1);
  } else if ( c_len >= a_len && c_len >= b_len ) {
    node_split = tmgr.sortAndSplit(node_begin,node_end,2);
  }

  if ( node_split == node_begin ) {  // unable to split torus
    // make sure physical nodes are together
    std::sort(node_begin, node_end, WorkDistrib::pe_sortop_compact());
    // find physical node boundary to split on
    int i_split = 0;
    for ( int i=0; i<nnodes; ++i ) {
      if ( ! CmiPeOnSamePhysicalNode(nodes[i_split],nodes[i]) ) {
        int mid = (nnodes+1)/2;
        if ( abs(i-mid) < abs(i_split-mid) ) i_split = i;
        else break;
      }
    }
    node_split = node_begin + i_split;
  }

  if ( node_split == node_begin ) {
    if ( simParams->verboseTopology ) {
      int crds[3];
      tmgr.coords(*node_begin, crds);
      CkPrintf("WorkDistrib: physnode %5d pe %5d node %5d at %5d %5d %5d from %5d %5d %5d has %5d patches %5d x %5d x %5d load %7f pes %5d\n",
               CmiPhysicalNodeID(*node_begin), *node_begin,
               CkNodeOf(*node_begin), crds[0], crds[1], crds[2],
               a_min, b_min, c_min, npatches,
               a_len+1, b_len+1, c_len+1, totalRawLoad, nnodes);
    }

    // final sort along a to minimize pme message count
    std::sort(patch_begin,patch_end,patch_sortop_curve_a(patchMap));

    // walk through patches in sorted order
    int *node = node_begin;
    sumLoad = 0;
    for ( int i=0; i < npatches; ++i ) {
      int pid = patches[i];
      assignedNode[pid] = *node;
      sumLoad += patchLoads[pid];
      double targetLoad = totalLoad *
        ((double)(node-node_begin+1) / (double)nnodes);
      if ( 0 ) CkPrintf("assign %5d node %5d patch %5d %5d %5d load %7f target %7f\n",
                i, *node,
                patchMap->index_a(pid),
                patchMap->index_b(pid),
                patchMap->index_c(pid),
                sumLoad, targetLoad);
      double extra = ( i+1 < npatches ? 0.5 * patchLoads[patches[i+1]] : 0 );
      if ( node+1 < node_end && sumLoad + extra >= targetLoad ) { ++node; }
    }

    return;
  }

  if ( a_len >= b_len && a_len >= c_len ) {
    if ( 0 ) CkPrintf("sort a\n");
    std::sort(patch_begin,patch_end,patch_sortop_curve_a(patchMap));
  } else if ( b_len >= a_len && b_len >= c_len ) {
    if ( 0 ) CkPrintf("sort b\n");
    std::sort(patch_begin,patch_end,patch_sortop_curve_b(patchMap));
  } else if ( c_len >= a_len && c_len >= b_len ) {
    if ( 0 ) CkPrintf("sort c\n");
    std::sort(patch_begin,patch_end,patch_sortop_curve_c(patchMap));
  }

  int *patch_split;
  { // walk through patches in sorted order
    int *node = node_begin;
    sumLoad = 0;
    for ( patch_split = patch_begin;
          patch_split != patch_end && node != node_split;
          ++patch_split ) {
      sumLoad += patchLoads[*patch_split];
      double targetLoad = totalLoad *
        ((double)(node-node_begin+1) / (double)nnodes);
      if ( 0 ) CkPrintf("test %5d node %5d patch %5d %5d %5d load %7f target %7f\n",
                patch_split - patch_begin, *node,
                patchMap->index_a(*patch_split),
                patchMap->index_b(*patch_split),
                patchMap->index_c(*patch_split),
                sumLoad, targetLoad);
      double extra = ( patch_split+1 != patch_end ? 0.5 * patchLoads[*(patch_split+1)] : 0 );
      if ( node+1 < node_end && sumLoad + extra >= targetLoad ) { ++node; }
    }
    double targetLoad = totalLoad *
      ((double)(node_split-node_begin) / (double)nnodes);
    if ( 0 ) CkPrintf("split node %5d/%5d patch %5d/%5d load %7f target %7f\n",
              node_split-node_begin, nnodes,
              patch_split-patch_begin, npatches,
              sumLoad, targetLoad);
  }

  // recurse
  recursive_bisect_with_curve(
    patch_begin, patch_split, node_begin, node_split,
    patchLoads, sortedLoads, assignedNode, tmgr);
  recursive_bisect_with_curve(
    patch_split, patch_end, node_split, node_end,
    patchLoads, sortedLoads, assignedNode, tmgr);
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesSpaceFillingCurve() 
{
  TopoManagerWrapper tmgr;
  PatchMap *patchMap = PatchMap::Object();
  const int numPatches = patchMap->numPatches();
  int *assignedNode = new int[numPatches];
  ResizeArray<double> patchLoads(numPatches);
  SortableResizeArray<double> sortedLoads(numPatches);
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  if(simParams->simulateInitialMapping) {
          NAMD_die("simulateInitialMapping not supported by assignPatchesSpaceFillingCurve()");
	  numNodes = simParams->simulatedPEs;
  }

  ResizeArray<int> patchOrdering(numPatches);
  for ( int i=0; i<numPatches; ++i ) {
    patchOrdering[i] = i;
  }

  ResizeArray<int> nodeOrdering(numNodes);
  nodeOrdering.resize(0);
  for ( int i=0; i<numNodes; ++i ) {
    int pe = peDiffuseOrdering[(i+1)%numNodes];  // avoid 0 if possible
    if ( simParams->noPatchesOnZero && numNodes > 1 ) {
      if ( pe == 0 ) continue;
      if(simParams->noPatchesOnOne && numNodes > 2) {
        if ( pe == 1 ) continue;
      }
    }  
#ifdef MEM_OPT_VERSION
    if(simParams->noPatchesOnOutputPEs && numNodes-simParams->numoutputprocs >2) {
      if ( isOutputProcessor(pe) ) continue;
    }
#endif
    nodeOrdering.add(pe);
    if ( 0 ) CkPrintf("using pe %5d\n", pe);
  }

  int *node_begin = nodeOrdering.begin();
  int *node_end = nodeOrdering.end();
  if ( nodeOrdering.size() > numPatches ) {
    node_end = node_begin + numPatches;
  }
  std::sort(node_begin, node_end, pe_sortop_compact());

  int *basenode_begin = node_begin;
  int *basenode_end = node_end;
  if ( nodeOrdering.size() > 2*numPatches ) {
    basenode_begin = node_end;
    basenode_end = basenode_begin + numPatches;
    std::sort(basenode_begin, basenode_end, pe_sortop_compact());
  }

  if ( simParams->disableTopology ) {
    iout << iWARN << "IGNORING TORUS TOPOLOGY DURING PATCH PLACEMENT\n" << endi;
  }

  recursive_bisect_with_curve(
    patchOrdering.begin(), patchOrdering.end(),
    node_begin, node_end,
    patchLoads.begin(), sortedLoads.begin(), assignedNode, tmgr);

  std::sort(node_begin, node_end, pe_sortop_compact());

  int samenodecount = 0;

  for ( int pid=0; pid<numPatches; ++pid ) {
    int node = assignedNode[pid];
    patchMap->assignNode(pid, node);
    int nodeidx = std::lower_bound(node_begin, node_end, node,
                                   pe_sortop_compact()) - node_begin;
    int basenode = basenode_begin[nodeidx];
    patchMap->assignBaseNode(pid, basenode);
    if ( CmiPeOnSamePhysicalNode(node,basenode) ) ++samenodecount;
  }

  iout << iINFO << "Placed " << (samenodecount*100./numPatches) << "% of base nodes on same physical node as patch\n" << endi;

  delete [] assignedNode; 
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputes(void)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  DebugM(3,"Mapping computes\n");

  computeMap->allocateCids();

  // Handle full electrostatics
  if ( node->simParameters->fullDirectOn )
    mapComputeHomePatches(computeFullDirectType);
  if ( node->simParameters->FMAOn )
#ifdef DPMTA
    mapComputeHomePatches(computeDPMTAType);
#else
    NAMD_die("This binary does not include DPMTA (FMA).");
#endif
  if ( node->simParameters->PMEOn ) {
#ifdef DPME
    if ( node->simParameters->useDPME )
      mapComputeHomePatches(computeDPMEType);
    else {
      if (node->simParameters->useOptPME) {
	mapComputeHomePatches(optPmeType);
	if ( node->simParameters->pressureProfileEwaldOn )
	  mapComputeHomePatches(computeEwaldType);
      }
      else {
	mapComputeHomePatches(computePmeType);
	if ( node->simParameters->pressureProfileEwaldOn )
	  mapComputeHomePatches(computeEwaldType);
      }
    }
#else
    if (node->simParameters->useOptPME) {
      mapComputeHomePatches(optPmeType);
      if ( node->simParameters->pressureProfileEwaldOn )
	mapComputeHomePatches(computeEwaldType);
    }
    else {      
#ifdef NAMD_CUDA
      if (node->simParameters->usePMECUDA) {
        mapComputePatch(computePmeCUDAType);
      } else 
#endif
      {
        mapComputePatch(computePmeType);
      }
      if ( node->simParameters->pressureProfileEwaldOn )
	mapComputeHomePatches(computeEwaldType);
    }
#endif
  }

  if ( node->simParameters->globalForcesOn ) {
    DebugM(2,"adding ComputeGlobal\n");
    mapComputeHomePatches(computeGlobalType);
  }

  if ( node->simParameters->extForcesOn )
    mapComputeHomePatches(computeExtType);

  if ( node->simParameters->qmForcesOn )
    mapComputeHomePatches(computeQMType);

  if ( node->simParameters->GBISserOn )
    mapComputeHomePatches(computeGBISserType);

  if ( node->simParameters->MsmSerialOn )
    mapComputeHomePatches(computeMsmSerialType);
#ifdef CHARM_HAS_MSA
  else if ( node->simParameters->MSMOn )
    mapComputeHomePatches(computeMsmMsaType);
#else
  else if ( node->simParameters->MSMOn )
    mapComputeHomePatches(computeMsmType);
#endif

  if ( node->simParameters->FMMOn )
    mapComputeHomePatches(computeFmmType);

#ifdef NAMD_CUDA
  if (node->simParameters->useCUDA2) {
    mapComputeNode(computeNonbondedCUDA2Type);
  } else {
    mapComputeNode(computeNonbondedCUDAType);
  }
  mapComputeHomeTuples(computeExclsType);
  mapComputePatch(computeSelfExclsType);
#endif

#ifdef NAMD_MIC
  mapComputeNode(computeNonbondedMICType);
#endif

  mapComputeNonbonded();

  if ( node->simParameters->LCPOOn ) {
    mapComputeLCPO();
  }

  // If we're doing true pair interactions, no need for bonded terms.
  // But if we're doing within-group interactions, we do need them.
  if ( !node->simParameters->pairInteractionOn || 
      node->simParameters->pairInteractionSelf) { 
    mapComputeHomeTuples(computeBondsType);
    mapComputeHomeTuples(computeAnglesType);
    mapComputeHomeTuples(computeDihedralsType);
    mapComputeHomeTuples(computeImpropersType);
    mapComputeHomeTuples(computeCrosstermsType);
    mapComputePatch(computeSelfBondsType);
    mapComputePatch(computeSelfAnglesType);
    mapComputePatch(computeSelfDihedralsType);
    mapComputePatch(computeSelfImpropersType);
    mapComputePatch(computeSelfCrosstermsType);
  }

  if ( node->simParameters->goGroPair ) {
      // JLai
      mapComputeHomeTuples(computeGromacsPairType);
      mapComputePatch(computeSelfGromacsPairType);
    // End of JLai
  }

  if ( node->simParameters->drudeOn ) {
    mapComputeHomeTuples(computeTholeType);
    mapComputePatch(computeSelfTholeType);
    mapComputeHomeTuples(computeAnisoType);
    mapComputePatch(computeSelfAnisoType);
  }

  if ( node->simParameters->eFieldOn )
    mapComputePatch(computeEFieldType);
  /* BEGIN gf */
  if ( node->simParameters->mgridforceOn )
    mapComputePatch(computeGridForceType);
  /* END gf */
  if ( node->simParameters->stirOn )
    mapComputePatch(computeStirType);
  if ( node->simParameters->sphericalBCOn )
    mapComputePatch(computeSphericalBCType);
  if ( node->simParameters->cylindricalBCOn )
    mapComputePatch(computeCylindricalBCType);
  if ( node->simParameters->tclBCOn ) {
    mapComputeHomePatches(computeTclBCType);
  }
  if ( node->simParameters->constraintsOn )
    mapComputePatch(computeRestraintsType);
  if ( node->simParameters->consForceOn )
    mapComputePatch(computeConsForceType);
  if ( node->simParameters->consTorqueOn )
    mapComputePatch(computeConsTorqueType);

    // store the latest compute map
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->storeComputeMap) {
    computeMap->saveComputeMap(simParams->computeMapFilename);
  }
    // override mapping decision
  if (simParams->loadComputeMap) {
    computeMap->loadComputeMap(simParams->computeMapFilename);
    CkPrintf("ComputeMap has been loaded from %s.\n", simParams->computeMapFilename);
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomeTuples(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int numNodes = node->numNodes();
  SimParameters *simparam = node->simParameters;
  if(simparam->simulateInitialMapping) {
	  numNodes = simparam->simulatedPEs;
  }

  char *isBaseNode = new char[numNodes];
  memset(isBaseNode,0,numNodes*sizeof(char));

  int numPatches = patchMap->numPatches();
  for(int j=0; j<numPatches; j++) {
    isBaseNode[patchMap->basenode(j)] = 1;
  }

  for(int i=0; i<numNodes; i++) {
    if ( isBaseNode[i] ) {
      computeMap->storeCompute(i,0,type);
    }
  }

  delete [] isBaseNode;
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomePatches(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int numNodes = node->numNodes();
  SimParameters *simparam = node->simParameters;
  if(simparam->simulateInitialMapping) {
	  numNodes = simparam->simulatedPEs;
  }

  for(int i=0; i<numNodes; i++) {
    if ( patchMap->numPatchesOnNode(i) ) {
      computeMap->storeCompute(i,0,type);
    }
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputePatch(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID i;
  ComputeID cid;

  for(i=0; i<patchMap->numPatches(); i++)
  {
    cid=computeMap->storeCompute(patchMap->node(i),1,type);
    computeMap->newPid(cid,i);
    patchMap->newCid(i,cid);
  }

}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeNode(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID i;
  ComputeID cid;

  int ncpus = CkNumPes();
  SimParameters *simparam = Node::Object()->simParameters;
  if(simparam->simulateInitialMapping) {
	  ncpus = simparam->simulatedPEs;
  }

  for(int i=0; i<ncpus; i++) {
    computeMap->storeCompute(i,0,type);
  }

}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeNonbonded(void)
{
  // For each patch, create 1 electrostatic object for self-interaction.
  // Then create 1 for each 1-away and 2-away neighbor which has a larger
  // pid.

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = Node::Object()->simParameters;
  int ncpus = CkNumPes();
  int nodesize = CkMyNodeSize();
  if(simParams->simulateInitialMapping) {
	  ncpus = simParams->simulatedPEs;
	  nodesize = simParams->simulatedNodeSize;
  }

  PatchID oneAway[PatchMap::MaxOneOrTwoAway];
  PatchID oneAwayDownstream[PatchMap::MaxOneOrTwoAway];
  int oneAwayTrans[PatchMap::MaxOneOrTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;
  double partScaling = 1.0;
  if ( ncpus < patchMap->numPatches() ) {
    partScaling = ((double)ncpus) / ((double)patchMap->numPatches());
  }

  for(i=0; i<patchMap->numPatches(); i++) // do the self 
  {

   int numPartitions = 1;
#if 0
   if ( simParams->ldBalancer == LDBAL_HYBRID ) {
#ifdef  MEM_OPT_VERSION    
    int64 numFixed = patchMap->numFixedAtoms(i);
    int64 numAtoms = patchMap->numAtoms(i);
#else
    int64 numFixed = patchMap->patch(i)->getNumFixedAtoms();  // avoid overflow
    int64 numAtoms = patchMap->patch(i)->getNumAtoms();
#endif

    int divide = node->simParameters->numAtomsSelf;
    if (divide > 0) {
      numPartitions = (int) ( partScaling * ( 0.5 +
        (numAtoms*numAtoms-numFixed*numFixed) / (double)(2*divide*divide) ) );
    }
    if (numPartitions < 1) numPartitions = 1;
    if ( numPartitions > node->simParameters->maxSelfPart )
			numPartitions = node->simParameters->maxSelfPart;
    // self-interaction
    DebugM(4,"Mapping " << numPartitions << " ComputeNonbondedSelf objects for patch " << i << "\n");
//    iout <<"Self numPartitions = " <<numPartitions <<" numAtoms " <<numAtoms <<std::endl;
   }
#endif

    // DMK - NOTE - For MIC builds (i.e. NAMD_MIC is defined), it is assumed that self computes are
    //   mapped to the PE their associated patch is on. If the code below should change, making that
    //   untrue, MIC builds should be special cased so that assumption still holds (or the host vs
    //   device load balancing scheme should be modified).  (See the comment in the function
    //   mic_assignComputes() in ComputeNonbondedMIC.C for more details.)
    for(int partition=0; partition < numPartitions; partition++)
    {
      cid=computeMap->storeCompute(patchMap->node(i),1,
				   computeNonbondedSelfType,
				   partition,numPartitions);
      computeMap->newPid(cid,i);
      patchMap->newCid(i,cid);
    }
  }

  for(int p1=0; p1 <patchMap->numPatches(); p1++) // do the pairs
  {
    // this only returns half of neighbors, which is what we want
    numNeighbors=patchMap->oneOrTwoAwayNeighbors(p1,oneAway,oneAwayDownstream,oneAwayTrans);
    for(j=0;j<numNeighbors;j++)
    {
	int p2 = oneAway[j];
	int dsp = oneAwayDownstream[j];

      int numPartitions = 1;
#if 0
      if ( simParams->ldBalancer == LDBAL_HYBRID ) {
#ifdef  MEM_OPT_VERSION        
        int64 numAtoms1 = patchMap->numAtoms(p1);
        int64 numAtoms2 = patchMap->numAtoms(p2);
        int64 numFixed1 = patchMap->numFixedAtoms(p1);
        int64 numFixed2 = patchMap->numFixedAtoms(p2);
#else
        int64 numAtoms1 = patchMap->patch(p1)->getNumAtoms();
        int64 numAtoms2 = patchMap->patch(p2)->getNumAtoms();
        int64 numFixed1 = patchMap->patch(p1)->getNumFixedAtoms();
        int64 numFixed2 = patchMap->patch(p2)->getNumFixedAtoms();
#endif


        const int t2 = oneAwayTrans[j];
        const int adim = patchMap->gridsize_a();
        const int bdim = patchMap->gridsize_b();
        const int cdim = patchMap->gridsize_c();
        const int nax = patchMap->numaway_a();  // 1 or 2
        const int nay = patchMap->numaway_b();  // 1 or 2
        const int naz = patchMap->numaway_c();  // 1 or 2
        const int ia1 = patchMap->index_a(p1);
        const int ia2 = patchMap->index_a(p2) + adim * Lattice::offset_a(t2);
        const int ib1 = patchMap->index_b(p1);
        const int ib2 = patchMap->index_b(p2) + bdim * Lattice::offset_b(t2);
        const int ic1 = patchMap->index_c(p1);
        const int ic2 = patchMap->index_c(p2) + cdim * Lattice::offset_c(t2);

        if ( abs(ia2-ia1) > nax ||
             abs(ib2-ib1) > nay ||
             abs(ic2-ic1) > naz )
          NAMD_bug("Bad patch distance in WorkDistrib::mapComputeNonbonded");

	int distance = 3;
 	if ( ia1 == ia2 ) --distance;
 	else if ( ia1 == ia2 + nax - 1 ) --distance;
 	else if ( ia1 + nax - 1 == ia2 ) --distance;
 	if ( ib1 == ib2 ) --distance;
 	else if ( ib1 == ib2 + nay - 1 ) --distance;
 	else if ( ib1 + nay - 1 == ib2 ) --distance;
 	if ( ic1 == ic2 ) --distance;
 	else if ( ic1 == ic2 + naz - 1 ) --distance;
 	else if ( ic1 + naz - 1 == ic2 ) --distance;
	int divide = 0;
	if ( distance == 0 ) {
	  divide = node->simParameters->numAtomsSelf2;
        } else if (distance == 1) {
	  divide = node->simParameters->numAtomsPair;
	} else {
	  divide = node->simParameters->numAtomsPair2;
	}
	if (divide > 0) {
          numPartitions = (int) ( partScaling * ( 0.5 +
	    (numAtoms1*numAtoms2-numFixed1*numFixed2)/(double)(divide*divide) ) );
	}
        if ( numPartitions < 1 ) numPartitions = 1;
        if ( numPartitions > node->simParameters->maxPairPart )
			numPartitions = node->simParameters->maxPairPart;
//	if ( numPartitions > 1 ) iout << "Mapping " << numPartitions << " ComputeNonbondedPair objects for patches " << p1 << "(" << numAtoms1 << ") and " << p2 << "(" << numAtoms2 << ")\n" << endi;
      }
#endif
		for(int partition=0; partition < numPartitions; partition++)
		{
		  cid=computeMap->storeCompute( patchMap->basenode(dsp),
			2,computeNonbondedPairType,partition,numPartitions);
		  computeMap->newPid(cid,p1);
		  computeMap->newPid(cid,p2,oneAwayTrans[j]);
		  patchMap->newCid(p1,cid);
		  patchMap->newCid(p2,cid);
		}
    }
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeLCPO(void) {
  //iterate over all needed objects

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = Node::Object()->simParameters;
  int ncpus = CkNumPes();
  int nodesize = CkMyNodeSize();
  const int maxPatches = 8;

  int numPatchesInOctet;
  PatchID patchesInOctet[maxPatches];
  int oneAwayTrans[maxPatches];

  //partitioned after 1st timestep
  int numPartitions = 1;

  PatchID i;
  ComputeID cid;

  // one octet per patch
  for(i=0; i<patchMap->numPatches(); i++) {
    numPatchesInOctet =
        patchMap->getPatchesInOctet(i, patchesInOctet, oneAwayTrans);

		for(int partition=0; partition < numPartitions; partition++) {
      cid=computeMap->storeCompute(patchMap->node(i),
          numPatchesInOctet,
				  computeLCPOType,
				  partition,
          numPartitions);
      for (int p = 0; p < numPatchesInOctet; p++) {
        computeMap->newPid(cid, patchesInOctet[p], oneAwayTrans[p]);
      }
      for (int p = 0; p < numPatchesInOctet; p++) {
        patchMap->newCid(patchesInOctet[p],cid);
      }
    } // for partitions 
  } // for patches
} // mapComputeLCPO

//----------------------------------------------------------------------
void WorkDistrib::messageEnqueueWork(Compute *compute) {
  LocalWorkMsg *msg = compute->localWorkMsg;
  int seq = compute->sequence();
  int gbisPhase = compute->getGBISPhase();

  if ( seq < 0 ) {
    NAMD_bug("compute->sequence() < 0 in WorkDistrib::messageEnqueueWork");
  } else {
    SET_PRIORITY(msg,seq,compute->priority());
  }

  msg->compute = compute; // pointer is valid since send is to local Pe
  int type = compute->type();
  int cid = compute->cid;

  CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
  switch ( type ) {
  case computeExclsType:
  case computeSelfExclsType:
    wdProxy[CkMyPe()].enqueueExcls(msg);
    break;
  case computeBondsType:
  case computeSelfBondsType:
    wdProxy[CkMyPe()].enqueueBonds(msg);
    break;
  case computeAnglesType:
  case computeSelfAnglesType:
    wdProxy[CkMyPe()].enqueueAngles(msg);
    break;
  case computeDihedralsType:
  case computeSelfDihedralsType:
    wdProxy[CkMyPe()].enqueueDihedrals(msg);
    break;
  case computeImpropersType:
  case computeSelfImpropersType:
    wdProxy[CkMyPe()].enqueueImpropers(msg);
    break;
  case computeTholeType:
  case computeSelfTholeType:
    wdProxy[CkMyPe()].enqueueThole(msg);
    break;
  case computeAnisoType:
  case computeSelfAnisoType:
    wdProxy[CkMyPe()].enqueueAniso(msg);
    break;
  case computeCrosstermsType:
  case computeSelfCrosstermsType:
    wdProxy[CkMyPe()].enqueueCrossterms(msg);
    break;
  // JLai
  case computeGromacsPairType:
  case computeSelfGromacsPairType:
    wdProxy[CkMyPe()].enqueueGromacsPair(msg);
    break;    
  // End of JLai
  case computeLCPOType:
    wdProxy[CkMyPe()].enqueueLCPO(msg);
    break;
  case computeNonbondedSelfType:
    switch ( seq % 2 ) {
    case 0:
      //wdProxy[CkMyPe()].enqueueSelfA(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueSelfA1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueSelfA2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueSelfA3(msg);
           break;
      }
      break;
    case 1:
      //wdProxy[CkMyPe()].enqueueSelfB(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueSelfB1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueSelfB2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueSelfB3(msg);
           break;
      }
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueSelf case statement error!");
    }
    break;
  case computeNonbondedPairType:
    switch ( seq % 2 ) {
    case 0:
      //wdProxy[CkMyPe()].enqueueWorkA(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueWorkA1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueWorkA2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueWorkA3(msg);
           break;
      }
      break;
    case 1:
      //wdProxy[CkMyPe()].enqueueWorkB(msg);
      switch ( gbisPhase ) {
         case 1:
           wdProxy[CkMyPe()].enqueueWorkB1(msg);
           break;
         case 2:
           wdProxy[CkMyPe()].enqueueWorkB2(msg);
           break;
         case 3:
           wdProxy[CkMyPe()].enqueueWorkB3(msg);
           break;
      }
      break;
    case 2:
      wdProxy[CkMyPe()].enqueueWorkC(msg);
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueWork case statement error!");
    }
    break;
  case computeNonbondedCUDAType:
#ifdef NAMD_CUDA
  case computeNonbondedCUDA2Type:
//     CkPrintf("WorkDistrib[%d]::CUDA seq=%d phase=%d\n", CkMyPe(), seq, gbisPhase);
    //wdProxy[CkMyPe()].enqueueCUDA(msg);
    switch ( gbisPhase ) {
       case 1:
         wdProxy[CkMyPe()].enqueueCUDA(msg);
         break;
       case 2:
         wdProxy[CkMyPe()].enqueueCUDAP2(msg);
         break;
       case 3:
         wdProxy[CkMyPe()].enqueueCUDAP3(msg);
         break;
    }
#else
    msg->compute->doWork();  MACHINE_PROGRESS
#endif
    break;
  case computeNonbondedMICType:
#ifdef NAMD_MIC
    wdProxy[CkMyPe()].enqueueMIC(msg);
#endif
    break;
  case computePmeType:
    // CkPrintf("PME %d %d %x\n", CkMyPe(), seq, compute->priority());
    wdProxy[CkMyPe()].enqueuePme(msg);
    break;
#ifdef NAMD_CUDA
  case computePmeCUDAType:
    wdProxy[CkMyPe()].enqueuePme(msg);
    break;
#endif
  case optPmeType:
    // CkPrintf("PME %d %d %x\n", CkMyPe(), seq, compute->priority());
#ifdef NAMD_CUDA
    wdProxy[CkMyPe()].enqueuePme(msg);
#else
    msg->compute->doWork();  MACHINE_PROGRESS
#endif
    break;
  default:
    wdProxy[CkMyPe()].enqueueWork(msg);
  }
}

//----------------------------------------------------------------------
void WorkDistrib::messageFinishCUDA(Compute *compute) {
  LocalWorkMsg *msg = compute->localWorkMsg;
  int seq = compute->sequence();
  int gbisPhase = compute->getGBISPhase();

  if ( seq < 0 ) {
    NAMD_bug("compute->sequence() < 0 in WorkDistrib::messageEnqueueWork");
  } else {
    SET_PRIORITY(msg,seq,compute->priority());
  }

  msg->compute = compute; // pointer is valid since send is to local Pe
  CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);

#ifdef NAMD_CUDA
    //wdProxy[CkMyPe()].finishCUDA(msg);
    switch ( gbisPhase ) {
       case 1:
         wdProxy[CkMyPe()].finishCUDA(msg);
         break;
       case 2:
         wdProxy[CkMyPe()].finishCUDAP2(msg);
         break;
       case 3:
         wdProxy[CkMyPe()].finishCUDAP3(msg);
         break;
    }
#else
    msg->compute->doWork();  MACHINE_PROGRESS
#endif
}

//----------------------------------------------------------------------
void WorkDistrib::messageFinishMIC(Compute *compute) {
  LocalWorkMsg *msg = compute->localWorkMsg;
  int seq = compute->sequence();
  int gbisPhase = compute->getGBISPhase();

  if ( seq < 0 ) {
    NAMD_bug("compute->sequence() < 0 in WorkDistrib::messageFinishMIC");
  } else {
    SET_PRIORITY(msg,seq,compute->priority());
  }

  msg->compute = compute; // pointer is valid since send is to local Pe
  CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);

#ifdef NAMD_MIC
    wdProxy[CkMyPe()].finishMIC(msg);
#else
    msg->compute->doWork();  MACHINE_PROGRESS
#endif
}

void WorkDistrib::enqueueWork(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueExcls(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueBonds(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueAngles(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueDihedrals(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueImpropers(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueThole(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueAniso(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueCrossterms(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

// JLai
void WorkDistrib::enqueueGromacsPair(LocalWorkMsg *msg) {
  msg->compute->doWork();
  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("\nWorkDistrib LocalWorkMsg recycling failed! Check enqueueGromacsPair from WorkDistrib.C\n");
}
// End of JLai

void WorkDistrib::enqueuePme(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueLCPO(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfA1(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfA2(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfA3(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueSelfB1(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfB2(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueSelfB3(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkA1(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkA2(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkA3(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkB1(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkB2(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueWorkB3(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}



void WorkDistrib::enqueueWorkC(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueCUDA(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  // ComputeNonbondedCUDA *c = msg->compute;
  // if ( c->localWorkMsg != msg && c->localWorkMsg2 != msg )
  //   NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::enqueueCUDAP2(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
}
void WorkDistrib::enqueueCUDAP3(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
}

void WorkDistrib::finishCUDAPatch(FinishWorkMsg *msg) {
  msg->compute->finishPatch(msg->data);
}

void WorkDistrib::finishCUDA(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
  // ComputeNonbondedCUDA *c = msg->compute;
  // if ( c->localWorkMsg != msg && c->localWorkMsg2 != msg )
  //   NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}
void WorkDistrib::finishCUDAP2(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
}
void WorkDistrib::finishCUDAP3(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
}

void WorkDistrib::enqueueMIC(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
}
void WorkDistrib::finishMIC(LocalWorkMsg *msg) {
  msg->compute->doWork();  MACHINE_PROGRESS
}


//**********************************************************************
//
//			FUNCTION velocities_from_PDB
//
//   INPUTS:
//      v - Array of vectors to populate
//	filename - name of the PDB filename to read in
//
//	This function reads in a set of initial velocities from a
//      PDB file.  It places the velocities into the array of Vectors
//      passed to it.
//
//***********************************************************************/

void WorkDistrib::velocities_from_PDB(const char *filename, 
				      Vector *v, int totalAtoms)
{
  PDB *v_pdb;		//  PDB info from velocity PDB
  int i;

  //  Read the PDB
  v_pdb = new PDB(filename);
  if ( v_pdb == NULL )
  {
    NAMD_die("memory allocation failed in Node::velocities_from_PDB");
  }

  //  Make sure the number of velocities read in matches
  //  the number of atoms we have
  if (v_pdb->num_atoms() != totalAtoms)
  {
    char err_msg[129];

    sprintf(err_msg, "FOUND %d COORDINATES IN VELOCITY PDB!!",
	    v_pdb->num_atoms());

    NAMD_die(err_msg);
  }

  //  Get the entire list of atom info and loop through
  //  them assigning the velocity vector for each one
  v_pdb->get_all_positions(v);

  for (i=0; i<totalAtoms; i++)
  {
    v[i].x *= PDBVELINVFACTOR;
    v[i].y *= PDBVELINVFACTOR;
    v[i].z *= PDBVELINVFACTOR;
  }

  delete v_pdb;
}
//		END OF FUNCTION velocities_from_PDB

//**********************************************************************
//
// 			FUNCTION velocities_from_binfile
//
//    INPUTS:
// 	fname - File name to write velocities to
//	n - Number of atoms in system
//	vels - Array of velocity vectors
//					
//	This function writes out the velocities in binary format.  This is
//     done to preserve accuracy between restarts of namd.
//
//**********************************************************************

void WorkDistrib::velocities_from_binfile(const char *fname, Vector *vels, int n)
{
  read_binary_file(fname,vels,n);
}
//               END OF FUNCTION velocities_from_binfile

//**********************************************************************
//
//			FUNCTION random_velocities
//
//   INPUTS:
//	v - array of vectors to populate
//	Temp - Temperature to acheive
//
//	This function assigns a random velocity distribution to a
//   simulation to achieve a desired initial temperature.  The method
//   used here was stolen from the program X-PLOR.
//
//**********************************************************************

void WorkDistrib::random_velocities(BigReal Temp,Molecule *structure,
				    Vector *v, int totalAtoms)
{
  int i, j;		//  Loop counter
  BigReal kbT;		//  Boltzman constant * Temp
  BigReal randnum;	//  Random number from -6.0 to 6.0
  BigReal kbToverM;	//  sqrt(Kb*Temp/Mass)
  SimParameters *simParams = Node::Object()->simParameters;
  Bool lesOn = simParams->lesOn;
  Random vel_random(simParams->randomSeed);

  int lesReduceTemp = lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  kbT = Temp*BOLTZMANN;

  //  Loop through all the atoms and assign velocities in
  //  the x, y and z directions for each one
  for (i=0; i<totalAtoms; i++)
  {
    if (structure->atommass(i) <= 0.) {
      kbToverM = 0.;
    } else {
      kbToverM = sqrt(kbT *
        ( lesOn && structure->get_fep_type(i) ? tempFactor : 1.0 ) /
			  structure->atommass(i) );
    }

    //  The following comment was stolen from X-PLOR where
    //  the following section of code was adapted from.
    
    //  This section generates a Gaussian random
    //  deviate of 0.0 mean and standard deviation RFD for
    //  each of the three spatial dimensions.
    //  The algorithm is a "sum of uniform deviates algorithm"
    //  which may be found in Abramowitz and Stegun,
    //  "Handbook of Mathematical Functions", pg 952.
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[i].x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[i].y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;
    
    v[i].z = randnum*kbToverM;
  }

  if ( simParams->drudeOn ) for (i=0; i<totalAtoms; i++) {
    if ( structure->is_drude(i) ) {
      v[i] = v[structure->get_mother_atom(i)];  // zero is good enough
    }
  }
}
/*			END OF FUNCTION random_velocities		*/

//**********************************************************************
//
//			FUNCTION remove_com_motion
//
//   INPUTS:
//	vel - Array of initial velocity vectors
//
//	This function removes the center of mass motion from a molecule.
//
//**********************************************************************

void WorkDistrib::remove_com_motion(Vector *vel, Molecule *structure, int n)
{
  Vector mv(0,0,0);		//  Sum of (mv)_i
  BigReal totalMass=0; 	//  Total mass of system
  int i;			//  Loop counter

  //  Loop through and compute the net momentum
  for (i=0; i<n; i++)
  {
    BigReal mass = structure->atommass(i);
    mv += mass * vel[i];
    totalMass += mass;
  }

  mv /= totalMass;

  iout << iINFO << "REMOVING COM VELOCITY "
	<< ( PDBVELFACTOR * mv ) << "\n" << endi;

  for (i=0; i<n; i++) { vel[i] -= mv; }

}
/*			END OF FUNCTION remove_com_motion		*/

#if USE_TOPOMAP 

//Specifically designed for BGL and other 3d Tori architectures
//Partition Torus and Patch grid together using recursive bisection.
int WorkDistrib::assignPatchesTopoGridRecBisection() {
  
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  if(simParams->simulateInitialMapping) {
	  numNodes = simParams->simulatedPEs;
  }

  int usedNodes = numNodes;
  CkPrintf("assignPatchesTopoGridRecBisection\n");
  if ( simParams->noPatchesOnZero && numNodes > 1 ) {
    usedNodes -= 1;
    if ( simParams->noPatchesOnOne && numNodes > 2 )
      usedNodes -= 1;
  }
  RecBisection recBisec(patchMap->numPatches(), PatchMap::Object());
  
  int xsize = 0, ysize = 0, zsize = 0;
  
  // Right now assumes a T*** (e.g. TXYZ) mapping
  TopoManager tmgr;
  xsize = tmgr.getDimNX();
  ysize = tmgr.getDimNY();
  zsize = tmgr.getDimNZ();
  
  //Fix to not assign patches to processor 0
  int rc = recBisec.partitionProcGrid(xsize, ysize, zsize, assignedNode);
 
  delete [] assignedNode;

  return rc;
}
#endif


#if defined(NAMD_MIC)
  extern void mic_hostDeviceLDB();
  extern void mic_contributeHostDeviceLDB(int idLen, int * id);
  extern void mic_setDeviceLDBParams(int dt, int hs, int sp1, int pp1, int pp2);
#endif

void WorkDistrib::send_initHostDeviceLDB() {
  #if defined(NAMD_MIC)
    CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
    wdProxy.initHostDeviceLDB();
  #endif
}

void WorkDistrib::initHostDeviceLDB() {
  #if defined(NAMD_MIC)
    mic_hostDeviceLDB();
  #endif
}

void WorkDistrib::send_contributeHostDeviceLDB(int peSetLen, int * peSet) {
  #if defined(NAMD_MIC)
    CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
    wdProxy[0].contributeHostDeviceLDB(peSetLen, peSet);
  #endif
}

void WorkDistrib::contributeHostDeviceLDB(int peSetLen, int * peSet) {
  #if defined(NAMD_MIC)
    mic_contributeHostDeviceLDB(peSetLen, peSet);
  #endif
}

void WorkDistrib::send_setDeviceLDBParams(int dt, int hs, int sp1, int pp1, int pp2) {
  #if defined(NAMD_MIC)
    CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
    wdProxy.setDeviceLDBParams(dt, hs, sp1, pp1, pp2);
  #endif
}

void WorkDistrib::setDeviceLDBParams(int dt, int hs, int sp1, int pp1, int pp2) {
  #if defined(NAMD_MIC)
    mic_setDeviceLDBParams(dt, hs, sp1, pp1, pp2);
  #endif
}


#include "WorkDistrib.def.h"

