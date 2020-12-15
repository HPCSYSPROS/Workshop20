#ifndef COMPUTEPMECUDA_H
#define COMPUTEPMECUDA_H

#include <vector>
#include <list>
#include "PmeBase.h"

#include "PatchTypes.h"       // Results
#include "Compute.h"
#include "Box.h"
#include "OwnerBox.h"
#include "ComputePmeCUDAMgr.decl.h"

#ifdef NAMD_CUDA
class HomePatch;

class ComputePmeCUDA : public Compute {
public:
  ComputePmeCUDA(ComputeID c, PatchIDList& pids);
  ComputePmeCUDA(ComputeID c, PatchID pid);
  virtual ~ComputePmeCUDA();
  void initialize();
  void atomUpdate();
  int noWork();
  void doWork();
  bool storePmeForceMsg(PmeForceMsg *msg);
private:
  struct PatchRecord {
    PatchRecord() {
      pmeForceMsg = NULL;
      patch = NULL;
      positionBox = NULL;
      avgPositionBox = NULL;
      forceBox = NULL;
    }
    // Message that contains the pointers to forces
    PmeForceMsg* pmeForceMsg;
    // Home pencil
    int homePencilY;
    int homePencilZ;
    int homePencilNode;
    // Pointer to patch
    Patch *patch;
    // Patch ID
    PatchID patchID;
    // Boxes
    Box<Patch,CompAtom> *positionBox;
    Box<Patch,CompAtom> *avgPositionBox;
    Box<Patch,Results> *forceBox;
  };

  double calcSelfEnergy(int numAtoms, CompAtom *x);
  void sendAtoms();
  void recvForces();
  void setupActivePencils();

  CmiNodeLock lock;
  int patchCounter;

  std::vector<PatchRecord> patches;

  bool selfEnergyDone;

  bool sendAtomsDone;

  PmeGrid pmeGrid;

  ComputePmeCUDAMgr *mgr;

  CProxy_ComputePmeCUDAMgr computePmeCUDAMgrProxy;

  bool atomsChanged;

};
#endif // NAMD_CUDA

#endif // COMPUTEPMECUDA_H