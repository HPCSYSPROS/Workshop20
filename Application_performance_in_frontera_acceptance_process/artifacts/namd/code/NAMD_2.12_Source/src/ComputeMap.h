/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMAP_H
#define COMPUTEMAP_H

#include "NamdTypes.h"
#include "ProcessorPrivate.h"
#include "ResizeArray.h"

class Compute;
class ComputeMgr;
template<class Type> class ObjectArena;

enum ComputeType
{
  computeNonbondedSelfType,
  computeNonbondedPairType,
  computeNonbondedCUDAType,
  computeNonbondedMICType,
  computeExclsType,
  computeBondsType,
  computeAnglesType,
  computeDihedralsType,
  computeImpropersType,
  computeTholeType,
  computeAnisoType,
  computeCrosstermsType,
  // JLai
  computeGromacsPairType,
  computeSelfGromacsPairType,
  // End of JLai entry
  computeSelfExclsType,
  computeSelfBondsType,
  computeSelfAnglesType,
  computeSelfDihedralsType,
  computeSelfImpropersType,
  computeSelfTholeType,
  computeSelfAnisoType,
  computeSelfCrosstermsType,
#ifdef DPMTA
  computeDPMTAType,
#endif
#ifdef DPME
  computeDPMEType,
#endif
  computePmeType,
#ifdef NAMD_CUDA
  computePmeCUDAType,
  computeNonbondedCUDA2Type,
#endif
  optPmeType,
  computeEwaldType,
  computeFullDirectType,
  computeGlobalType,
  computeExtType,
  computeQMType,
  computeGBISserType,
  computeLCPOType,
  computeFmmType,
  computeMsmSerialType,
  computeMsmMsaType,
  computeMsmType,
  computeEFieldType,
/* BEGIN gf */
  computeGridForceType,
/* END gf */
  computeStirType,
  computeSphericalBCType,
  computeCylindricalBCType,
  computeTclBCType,
  computeRestraintsType,
  computeConsForceType,
  computeConsTorqueType,
  computeErrorType
};

class ComputeMap
{
public:
  static ComputeMap *Instance();
  inline static ComputeMap *Object() { return instance; }

  void checkMap();

  ~ComputeMap(void);

  void registerCompute(ComputeID cid, Compute *c) {
    computePtrs[cid] = c;
  }

  // numComputes() returns the number of compute objects known
  // by the map.
  inline int numComputes(void) {
    return nComputes;
  }

  // node(cid) returns the node where the compute object currently exists.
  inline int node(ComputeID cid) {
    return computeData[cid].node;
  }

  inline void setNode(ComputeID cid, NodeID node) {
    computeData[cid].node = node;
  }

  // newNode(cid,node) sets up map to tell WorkDistrib to send 
  // compute to new node
  inline NodeID newNode(ComputeID cid) {
    return (computeData[cid].moveToNode);
  }

  inline void setNewNode(ComputeID cid, NodeID node) {
    computeData[cid].moveToNode = node;
  }

  // numPids(cid) returns the number of patch ids which are registered
  // with this compute object.
  int numPids(ComputeID cid);
  
  // pid(cid,i) returns the i-th patch id registered
  // with the patch.  
  int pid(ComputeID cid, int i);
  int trans(ComputeID cid, int i);

  // type(cid) returns the compute type of the given ComputeID
  ComputeType type(ComputeID cid);
  int partition(ComputeID cid);
  int numPartitions(ComputeID cid);

  inline void setNumPartitions(ComputeID cid, char numPartitions) {
    computeData[cid].numPartitions = numPartitions;
  }
  inline char newNumPartitions(ComputeID cid) {
    return (computeData[cid].newNumPartitions);
  }
  inline void setNewNumPartitions(ComputeID cid, char numPartitions) {
    computeData[cid].newNumPartitions = numPartitions;
  }

  int allocateCids();

  // storeCompute(cid,node,maxPids) tells the ComputeMap to store
  // information about the indicated patch, and allocate space
  // for up to maxPids dependents
  ComputeID storeCompute(int node,int maxPids,ComputeType type,
			 int partition=-1, int numPartitions=0);

  // newPid(cid,pid) stores the n patch ids associated with
  // compute id cid.
  void newPid(ComputeID cid, int pid, int trans = 13);

  #if defined(NAMD_MIC)
    void setDirectToDevice(const ComputeID cid, const int d);
    int directToDevice(const ComputeID cid) const;
  #endif

  ComputeID cloneCompute(ComputeID src, int partition);

  void printComputeMap(void);
  void saveComputeMap(const char *fname);
  void loadComputeMap(const char *fname);

  Compute *compute(ComputeID cid) { return computePtrs[cid]; };

  friend class ComputeMgr;

  struct PatchRec
  {
    PatchID pid;
    int trans;

    PatchRec() : pid(-1), trans(-1) { ; }
  };

  enum { numPidsAllocated=8 };

  struct ComputeData
  {
    ComputeData() { 
      node = -1; moveToNode = -1; 
      newNumPartitions = 0;
      numPids = 0;
      #if defined(NAMD_MIC)
        directToDevice = 0;
      #endif
    }
    int node;
    int moveToNode;
    ComputeType type;
    char partition;
    char numPartitions;
    char newNumPartitions;
    char numPids;
    PatchRec pids[numPidsAllocated];
    #if defined(NAMD_MIC)
      char directToDevice;
    #endif
  };
protected:
  friend class WorkDistrib;
  void pack(ComputeData *buf);
  void unpack(int n, ComputeData *buf);
  void initPtrs();
  void extendPtrs();

  ComputeMap(void);

private:
  int nComputes;
  ResizeArray<ComputeData> computeData;
  Compute **computePtrs;

  static ComputeMap *instance;
};

#endif /* COMPUTEMAP_H */

