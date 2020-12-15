/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PATCHMAP_H
#define PATCHMAP_H

#include "NamdTypes.h"
#include "HomePatchList.h"
#include "Lattice.h"
#include "ProcessorPrivate.h"

#include <vector>


class Patch;
class PatchMgr;
class HomePatch;
template<class Type> class ObjectArena;

class PatchMap
{
public:
  static PatchMap *Instance();
  inline static PatchMap *Object() { return CkpvAccess(PatchMap_instance); }
  inline static PatchMap *ObjectOnPe(int pe) {
    return CkpvAccessOther(PatchMap_instance, CmiRankOf(pe));
  }

  int sizeGrid(ScaledPosition xmin, ScaledPosition xmax,
			const Lattice &lattice, BigReal patchSize,
			double maxNumPatches, int staticAtomAssignment,
			int asplit, int bsplit, int csplit);
  void makePatches(ScaledPosition xmin, ScaledPosition xmax,
			const Lattice &lattice, BigReal patchSize,
			double maxNumPatches, int staticAtomAssignment,
			int replicaUniformPatchGrids, int lcpo,
			int asplit, int bsplit, int csplit);
  void checkMap();

  ~PatchMap(void);

  enum { MaxTwoAway = 5*5*5 - 3*3*3 };
  enum { MaxOneAway = 3*3*3 - 1 };
  enum { MaxOneOrTwoAway = MaxOneAway + MaxTwoAway };

  static void registerPatchMgr(PatchMgr *pmgr) {
    CkpvAccess(PatchMap_patchMgr) = pmgr;
  }

  HomePatchList *homePatchList();
  void homePatchIDList(PatchIDList &);  // expensive - for startup only
  void basePatchIDList(int pe, PatchIDList &);  // use for required proxies
  int numHomePatches(void);

  // returns the number of patches being managed 
  inline int numPatches(void) const { return nPatches; }
  inline int numPatchesOnNode(int node) { return nPatchesOnNode[node]; }
  inline int numNodesWithPatches(void) { return nNodesWithPatches; }

  // returns the number of patches in each dimension
  inline int gridsize_a(void) const { return aDim; }
  inline int gridsize_b(void) const { return bDim; }
  inline int gridsize_c(void) const { return cDim; }
  // returns the number of patches in each dimension
  inline int numaway_a(void) const { return aAway; }
  inline int numaway_b(void) const { return bAway; }
  inline int numaway_c(void) const { return cAway; }

  // returns 1 if periodic in each dimension
  inline int periodic_a(void) const { return aPeriodic; }
  inline int periodic_b(void) const { return bPeriodic; }
  inline int periodic_c(void) const { return cPeriodic; }

  // returns the origin (minimum, not center) of patch grid
  inline ScaledPosition origin(void) const {
    return ScaledPosition(aOrigin,bOrigin,cOrigin);
  }

  // returns the patch id for the given indices
  inline int pid(int aIndex, int bIndex, int cIndex);

  // returns the [abc] index for the given patch id.
  inline int index_a(int pid) const { return pid % aDim; }
  inline int index_b(int pid) const { return (pid / aDim) % bDim; }
  inline int index_c(int pid) const { return pid / (aDim*bDim); }

  // returns the min/max [abc] scaled coordinate
  inline BigReal min_a(int pid) const { return patchBounds_a[patchData[pid].aIndex*2]; }
  inline BigReal max_a(int pid) const { return patchBounds_a[patchData[pid].aIndex*2+2]; }
  inline BigReal min_b(int pid) const { return patchBounds_b[patchData[pid].bIndex*2]; }
  inline BigReal max_b(int pid) const { return patchBounds_b[patchData[pid].bIndex*2+2]; }
  inline BigReal min_c(int pid) const { return patchBounds_c[patchData[pid].cIndex*2]; }
  inline BigReal max_c(int pid) const { return patchBounds_c[patchData[pid].cIndex*2+2]; }

  // returns the center of patch scaled position
  inline ScaledPosition center(int pid) const {
    const PatchData &pd = patchData[pid];
    return ScaledPosition(patchBounds_a[pd.aIndex*2+1],
                          patchBounds_b[pd.bIndex*2+1],
                          patchBounds_c[pd.cIndex*2+1]);
  }

  // asssigns atom to patch based on position and lattice
  inline PatchID assignToPatch(Position p, const Lattice &l);

  // gives more downstream patch of pid1, pid2; handles periodicity right
  // given patches must be neighbors!!!
  inline int downstream(int pid1, int pid2);

  // returns the node where the patch currently exists.
  inline int node(int pid) const { return patchData[pid].node; }

  // returns the node where the patch's upstream proxies exist.
  inline int basenode(int pid) const { return patchData[pid].basenode; }

  // numCids(pid) returns the number of compute ids which are registered
  inline int numCids(int pid) const { return patchData[pid].numCids; }
  
  // cid(pid,i) returns the i-th compute id registered
  inline int cid(int pid, int i) const { return patchData[pid].cids[i]; }

#ifdef MEM_OPT_VERSION
  inline int numAtoms(int pid) const { return patchData[pid].numAtoms; }
  inline void setNumAtoms(int pid, int num) { patchData[pid].numAtoms = num; }

  inline int numFixedAtoms(int pid) const { return patchData[pid].numFixedAtoms; }
  inline void setNumFixedAtoms(int pid, int num) { patchData[pid].numFixedAtoms = num; }
#endif

  void assignNode(PatchID, NodeID);
  void assignBaseNode(PatchID, NodeID);
  void assignBaseNode(PatchID);

  // newCid(pid,cid) stores a compute id associated with
  // patch id pid.  Error returned when there is no room to store
  // the pid.
  void newCid(int pid, int cid);

  // oneAwayNeighbors(pid, neighbor_ids) returns the number 
  // and ids of adjacent patches.  The caller is expected to provide
  // sufficient storage for the neighbors.

  int oneAwayNeighbors(int pid, PatchID *neighbor_ids=0);

  int oneOrTwoAwayNeighbors(int pid, PatchID *neighbor_ids,
		        PatchID *downstream_ids = 0, int *transform_ids = 0);

  //LCPO
  int getPatchesInOctet(int pid, PatchID *pids, int *transform_ids = 0);

  int upstreamNeighbors(int pid, PatchID *neighbor_ids);

  int downstreamNeighbors(int pid, PatchID *neighbor_ids); 

  void printPatchMap(void);

  inline Patch *patch(PatchID pid);
  inline HomePatch *homePatch(PatchID pid);

  void registerPatch(PatchID pid, HomePatch *pptr);
  void unregisterPatch(PatchID pid, HomePatch *pptr);

  void registerPatch(PatchID pid, Patch *pptr);
  void unregisterPatch(PatchID pid, Patch *pptr);

protected:
  friend class WorkDistrib;
  int packSize(void);
  void pack(char *buf, int size);
  void unpack(char *buf);
  
  PatchMap(void);
  
private:
  struct PatchData
  {
    int node, basenode;
    short aIndex, bIndex, cIndex;
    short numCids;
    short numCidsAllocated;
    ComputeID *cids;
#ifdef MEM_OPT_VERSION
    //added to record #atoms in each patch initially
    //--Chao Mei
    unsigned short numAtoms;
    unsigned short numFixedAtoms;
#endif
  };
  int nPatches;
  int nNodesWithPatches;
  static int *nPatchesOnNode;
  static PatchData *patchData;
  static ObjectArena<ComputeID> *computeIdArena;
  BigReal *patchBounds_a;
  BigReal *patchBounds_b;
  BigReal *patchBounds_c;
  Patch **myPatch;
  HomePatch **myHomePatch;
  int aDim, bDim, cDim;
  int aAway, bAway, cAway;
  int aPeriodic, bPeriodic, cPeriodic;
  int aMaxIndex, bMaxIndex, cMaxIndex;
  BigReal aOrigin, bOrigin, cOrigin;
  BigReal aLength, bLength, cLength;

private:
  //It is used to store the atom ids that each patch has
  //we need this structure because we want to create and distribute
  //each patch one by one rather than creat all home patches at a time and then
  //send them later
  std::vector<int> *tmpPatchAtomsList;
public:
  void initTmpPatchAtomsList(){
      tmpPatchAtomsList = new std::vector<int>[nPatches];
  }
  void delTmpPatchAtomsList() {
      for(int i=0; i<nPatches; i++){
          tmpPatchAtomsList[i].clear();
      }
      delete [] tmpPatchAtomsList;
      tmpPatchAtomsList = NULL;
  }
  std::vector<int> *getTmpPatchAtomsList(){
      return tmpPatchAtomsList;
  }

};


//----------------------------------------------------------------------

inline Patch *PatchMap::patch(PatchID pid)
{
  return myPatch[pid];
}

HomePatch *PatchMap::homePatch(PatchID pid)
{
  return myHomePatch[pid];
}

#endif /* PATCHMAP_H */

