#include <algorithm>
#include <map>
#include <vector>
#include "NamdTypes.h"
#include "charm++.h"
#include "Patch.h"
#include "PatchMap.h"
#include "ProxyMgr.h"
#include "LJTable.h"
#include "Node.h"
#include "ObjectArena.h"
// #include "ComputeCUDAMgr.h"
#include "ReductionMgr.h"
#include "CudaComputeNonbonded.h"
#include "WorkDistrib.h"
#include "HomePatch.h"
#include "Priorities.h"
#include "ComputePmeCUDAMgr.h"
//#include "CudaUtils.h"

#include "DeviceCUDA.h"
#ifdef NAMD_CUDA
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;
#endif

#ifdef NAMD_CUDA

extern "C" void CcdCallBacksReset(void *ignored, double curWallTime);  // fix Charm++

//
// Class constructor
//
CudaComputeNonbonded::CudaComputeNonbonded(ComputeID c, int deviceID, bool doStreaming) : Compute(c),
deviceID(deviceID), doStreaming(doStreaming), nonbondedKernel(deviceID, doStreaming),
tileListKernel(deviceID, doStreaming), GBISKernel(deviceID) {

  cudaCheck(cudaSetDevice(deviceID));

	exclusionsByAtom = NULL;

  vdwTypes = NULL;
  vdwTypesSize = 0;

  exclIndexMaxDiff = NULL;
  exclIndexMaxDiffSize = 0;

  atomIndex = NULL;
  atomIndexSize = 0;

  atomStorageSize = 0;

  // Atom and charge storage
  atoms = NULL;
  atomsSize = 0;

  // Force storage
  h_forces = NULL;
  h_forcesSize = 0;
  h_forcesSlow = NULL;
  h_forcesSlowSize = 0;

  d_forces = NULL;
  d_forcesSize = 0;
  d_forcesSlow = NULL;
  d_forcesSlowSize = 0;

  // GBIS
  intRad0H = NULL;
  intRad0HSize = 0;
  intRadSH = NULL;
  intRadSHSize = 0;
  psiSumH = NULL;
  psiSumHSize = 0;
  bornRadH = NULL;
  bornRadHSize = 0;
  dEdaSumH = NULL;
  dEdaSumHSize = 0;
  dHdrPrefixH = NULL;
  dHdrPrefixHSize = 0;

  cudaPatches = NULL;

  atomsChangedIn = true;
  atomsChanged = true;
  computesChanged = true;

  forceDoneEventRecord = false;

  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->pressureProfileOn) {
    NAMD_die("CudaComputeNonbonded, pressure profile not supported");
  }

  if (simParams->GBISOn) gbisPhase = 3;

  doSkip = false;

#define CUDA_DEBUG_EVENT 171
  traceRegisterUserEvent("CUDA DEBUG", CUDA_DEBUG_EVENT);
#define CUDA_VDW_KERNEL 172
  traceRegisterUserEvent("CUDA VdW kernel", CUDA_VDW_KERNEL);
#define CUDA_GBIS1_KERNEL 173
  traceRegisterUserEvent("CUDA GBIS Phase 1 kernel", CUDA_GBIS1_KERNEL);
#define CUDA_GBIS2_KERNEL 174
  traceRegisterUserEvent("CUDA GBIS Phase 2 kernel", CUDA_GBIS2_KERNEL);
#define CUDA_GBIS3_KERNEL 175
  traceRegisterUserEvent("CUDA GBIS Phase 3 kernel", CUDA_GBIS3_KERNEL);
}

//
// Class destructor
//
CudaComputeNonbonded::~CudaComputeNonbonded() {
  cudaCheck(cudaSetDevice(deviceID));
	if (exclusionsByAtom != NULL) delete [] exclusionsByAtom;
  if (vdwTypes != NULL) deallocate_host<int>(&vdwTypes);
  if (exclIndexMaxDiff != NULL) deallocate_host<int2>(&exclIndexMaxDiff);
  if (atoms != NULL) deallocate_host<CudaAtom>(&atoms);
  if (h_forces != NULL) deallocate_host<float4>(&h_forces);
  if (h_forcesSlow != NULL) deallocate_host<float4>(&h_forcesSlow);
  if (d_forces != NULL) deallocate_device<float4>(&d_forces);
  if (d_forcesSlow != NULL) deallocate_device<float4>(&d_forcesSlow);

  // GBIS
  if (intRad0H != NULL) deallocate_host<float>(&intRad0H);
  if (intRadSH != NULL) deallocate_host<float>(&intRadSH);
  if (psiSumH != NULL) deallocate_host<GBReal>(&psiSumH);
  if (bornRadH != NULL) deallocate_host<float>(&bornRadH);
  if (dEdaSumH != NULL) deallocate_host<GBReal>(&dEdaSumH);
  if (dHdrPrefixH != NULL) deallocate_host<float>(&dHdrPrefixH);

  if (cudaPatches != NULL) deallocate_host<CudaPatchRecord>(&cudaPatches);

  if (patches.size() > 0) {
    deallocate_host<VirialEnergy>(&h_virialEnergy);
    deallocate_device<VirialEnergy>(&d_virialEnergy);
    cudaCheck(cudaStreamDestroy(stream));
    cudaCheck(cudaEventDestroy(forceDoneEvent));
    CmiDestroyLock(lock);
    delete reduction;
  }

  // NOTE: unregistering happens in [sync] -entry method
  computeMgr->sendUnregisterBoxesOnPe(pes, this);

}

void CudaComputeNonbonded::unregisterBox(int i) {
  if (patches[i].positionBox != NULL) patches[i].patch->unregisterPositionPickup(this, &patches[i].positionBox);
  if (patches[i].forceBox != NULL) patches[i].patch->unregisterForceDeposit(this, &patches[i].forceBox);
  if (patches[i].intRadBox != NULL) patches[i].patch->unregisterIntRadPickup(this, &patches[i].intRadBox);
  if (patches[i].psiSumBox != NULL) patches[i].patch->unregisterPsiSumDeposit(this, &patches[i].psiSumBox);
  if (patches[i].bornRadBox != NULL) patches[i].patch->unregisterBornRadPickup(this, &patches[i].bornRadBox);
  if (patches[i].dEdaSumBox != NULL) patches[i].patch->unregisterDEdaSumDeposit(this, &patches[i].dEdaSumBox);
  if (patches[i].dHdrPrefixBox != NULL) patches[i].patch->unregisterDHdrPrefixPickup(this, &patches[i].dHdrPrefixBox);
}

void CudaComputeNonbonded::unregisterBoxesOnPe() {
  if (rankPatches[CkMyRank()].size() == 0)
    NAMD_bug("CudaComputeNonbonded::unregisterBoxesOnPe, empty rank");
  for (int i=0;i < rankPatches[CkMyRank()].size();i++) {
    unregisterBox(rankPatches[CkMyRank()][i]);
  }
}

//
// Register inter-patch (self) compute.
// Only serialized calls allowed
//
void CudaComputeNonbonded::registerComputeSelf(ComputeID cid, PatchID pid) {
  computesChanged = true;
  addPatch(pid);
  addCompute(cid, pid, pid, 0.);
}

//
// Register pair-patch compute.
// Only serialized calls allowed
//
void CudaComputeNonbonded::registerComputePair(ComputeID cid, PatchID* pid, int* trans) {
  computesChanged = true;
  addPatch(pid[0]);
  addPatch(pid[1]);
  PatchMap* patchMap = PatchMap::Object();
  int t1 = trans[0];
  int t2 = trans[1];
  Vector offset = patchMap->center(pid[0]) - patchMap->center(pid[1]);
  offset.x += (t1%3-1) - (t2%3-1);
  offset.y += ((t1/3)%3-1) - ((t2/3)%3-1);
  offset.z += (t1/9-1) - (t2/9-1);
  addCompute(cid, pid[0], pid[1], offset);
}

//
// Add patch
//
void CudaComputeNonbonded::addPatch(PatchID pid) {
  patches.push_back(PatchRecord(pid));
}

//
// Add compute
//
void CudaComputeNonbonded::addCompute(ComputeID cid, PatchID pid1, PatchID pid2, Vector offset) {
  ComputeRecord cr;
  cr.cid = cid;
  cr.pid[0] = pid1;
  cr.pid[1] = pid2;
  cr.offset = offset;
  computes.push_back(cr);
}

//
// Update numAtoms and numFreeAtoms on a patch
//
void CudaComputeNonbonded::updatePatch(int i) {
  int numAtoms = patches[i].patch->getNumAtoms();
  int numFreeAtoms = numAtoms;
  if ( fixedAtomsOn ) {
    const CompAtomExt *aExt = patches[i].patch->getCompAtomExtInfo();
    for ( int j=0; j< numAtoms; ++j ) {
      if ( aExt[j].atomFixed ) --numFreeAtoms;
    }
  }
  patches[i].numAtoms = numAtoms;
  patches[i].numFreeAtoms = numFreeAtoms;
  cudaPatches[i].numAtoms = numAtoms;
  cudaPatches[i].numFreeAtoms = numFreeAtoms;
}

int CudaComputeNonbonded::findPid(PatchID pid) {
  for (int i=0;i < rankPatches[CkMyRank()].size();i++) {
    int j = rankPatches[CkMyRank()][i];
    if (patches[j].patchID == pid) return j;
  }
  return -1;
}

void CudaComputeNonbonded::patchReady(PatchID pid, int doneMigration, int seq) {
  if (doneMigration) {
    int i = findPid(pid);
    if (i == -1)
      NAMD_bug("CudaComputeNonbonded::patchReady, Patch ID not found");
    updatePatch(i);
  }
  CmiLock(lock);
  Compute::patchReady(pid, doneMigration, seq);
  CmiUnlock(lock);
}

void CudaComputeNonbonded::gbisP2PatchReady(PatchID pid, int seq) {
  CmiLock(lock);
  Compute::gbisP2PatchReady(pid, seq);
  CmiUnlock(lock);
}

void CudaComputeNonbonded::gbisP3PatchReady(PatchID pid, int seq) {
  CmiLock(lock);
  Compute::gbisP3PatchReady(pid, seq);
  CmiUnlock(lock);
}

void CudaComputeNonbonded::assignPatch(int i) {
  PatchMap* patchMap = PatchMap::Object();
  PatchID pid = patches[i].patchID;
  Patch* patch = patchMap->patch(pid);
  if (patch == NULL) {
    // Create ProxyPatch if none exists
    ProxyMgr::Object()->createProxy(pid);
    patch = patchMap->patch(pid);
  }
  patches[i].patch = patch;
  if (patches[i].patch == NULL) {
    NAMD_bug("CudaComputeNonbonded::assignPatch, patch not found");
  }
  patches[i].positionBox = patches[i].patch->registerPositionPickup(this);
  patches[i].forceBox    = patches[i].patch->registerForceDeposit(this);
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->GBISOn) {
    patches[i].intRadBox     = patches[i].patch->registerIntRadPickup(this);
    patches[i].psiSumBox     = patches[i].patch->registerPsiSumDeposit(this);
    patches[i].bornRadBox    = patches[i].patch->registerBornRadPickup(this);
    patches[i].dEdaSumBox    = patches[i].patch->registerDEdaSumDeposit(this);
    patches[i].dHdrPrefixBox = patches[i].patch->registerDHdrPrefixPickup(this);
  }
  // Store Pe where this patch was registered
  patches[i].pe = CkMyPe();
  //
  patches[i].isSamePhysicalNode = ( CmiPhysicalNodeID(patchMap->node(pid)) == CmiPhysicalNodeID(CkMyPe()) );
  patches[i].isSameNode = ( CkNodeOf(patchMap->node(pid)) == CkMyNode() );
}

struct pid_sortop_reverse_priority {
  bool operator() (int2 pidj, int2 pidi) {  // i and j reversed
    int ppi = PATCH_PRIORITY(pidi.x);
    int ppj = PATCH_PRIORITY(pidj.x);
    if ( ppi != ppj ) return ppi < ppj;
    return pidi.x < pidj.x;
  }
};

void CudaComputeNonbonded::assignPatchesOnPe() {
  if (rankPatches[CkMyRank()].size() == 0)
    NAMD_bug("CudaComputeNonbonded::assignPatchesOnPe, empty rank");

  // calculate priority rank of local home patch within pe
  {
    PatchMap* patchMap = PatchMap::Object();
    ResizeArray< ResizeArray<int2> > homePatchByRank(CkMyNodeSize());
    for ( int k=0; k < rankPatches[CkMyRank()].size(); ++k ) {
      int i = rankPatches[CkMyRank()][k];
      int pid = patches[i].patchID;
      int homePe = patchMap->node(pid);
      if ( CkNodeOf(homePe) == CkMyNode() ) {
        int2 pid_index;
        pid_index.x = pid;
        pid_index.y = i;
        homePatchByRank[CkRankOf(homePe)].add(pid_index);
      }
    }
    for ( int i=0; i<CkMyNodeSize(); ++i ) {
      pid_sortop_reverse_priority so;
      std::sort(homePatchByRank[i].begin(),homePatchByRank[i].end(),so);
      int masterBoost = ( CkMyRank() == i ? 2 : 0 );
      for ( int j=0; j<homePatchByRank[i].size(); ++j ) {
        int index = homePatchByRank[i][j].y;
        patches[index].reversePriorityRankInPe = j + masterBoost;
      }
    }
  }

  for (int i=0;i < rankPatches[CkMyRank()].size();i++) {
    assignPatch(rankPatches[CkMyRank()][i]);
  }
}

//
// Returns Pe of Patch ID "pid", -1 otherwise
//
// int findHomePatchPe(std::vector<PatchIDList>& rankPatchIDs, PatchID pid) {
int findHomePatchPe(PatchIDList* rankPatchIDs, PatchID pid) {
  // for (int i=0;i < rankPatchIDs.size();i++) {
  for (int i=0;i < CkMyNodeSize();i++) {
    if (rankPatchIDs[i].find(pid) != -1) return CkNodeFirst(CkMyNode()) + i;
  }
  return -1;
}

//
// Find all PEs that have Patch
//
void findProxyPatchPes(std::vector<int>& proxyPatchPes, PatchID pid) {
  proxyPatchPes.clear();
  for (int i=0;i < CkMyNodeSize();i++) {
    int pe = CkNodeFirst(CkMyNode()) + i;
    if (PatchMap::ObjectOnPe(pe)->patch(pid) != NULL)
      proxyPatchPes.push_back(pe);
  }
}

//
// Called after all computes have been registered
//
void CudaComputeNonbonded::assignPatches(ComputeMgr* computeMgrIn) {
  // Remove duplicate patches
  std::sort(patches.begin(), patches.end());
  std::vector<PatchRecord>::iterator last = std::unique(patches.begin(), patches.end());
  patches.erase(last, patches.end());
  //
  std::map<PatchID, int> pidMap;
  // Set number of patches and register boxes
  setNumPatches(patches.size());
  masterPe = CkMyPe();
  computeMgr = computeMgrIn;
  // Start patch counter
  patchesCounter = getNumPatches();
  // For each rank, list of patches
  rankPatches.resize(CkMyNodeSize());
  // For each rank, list of home patch IDs
  PatchIDList* rankHomePatchIDs = new PatchIDList[CkMyNodeSize()];
  for (int i=0;i < CkMyNodeSize();i++) {
    int pe = CkNodeFirst(CkMyNode()) + i;
    PatchMap::Object()->basePatchIDList(pe, rankHomePatchIDs[i]);
  }
  std::vector<int> proxyPatchPes;
  std::vector<int> peProxyPatchCounter(CkMyNodeSize(), 0);
  //--------------------------------------------------------
  // Build a list of PEs to avoid
  std::vector<int> pesToAvoid;
#if 0
  // Avoid other GPUs' master PEs
  for (int i=0;i < deviceCUDA->getDeviceCount();i++) {
    int pe = deviceCUDA->getMasterPeForDeviceID(i);
    if (pe != -1 && pe != masterPe) pesToAvoid.push_back(pe);
  }
  // Avoid PEs that are involved in PME
  ComputePmeCUDAMgr *computePmeCUDAMgr = ComputePmeCUDAMgr::Object();
  for (int pe=CkNodeFirst(CkMyNode());pe < CkNodeFirst(CkMyNode()) + CkMyNodeSize();pe++) {
    if (computePmeCUDAMgr->isPmePe(pe)) pesToAvoid.push_back(pe);
  }
  // Set counters of avoidable PEs to high numbers
  for (int i=0;i < pesToAvoid.size();i++) {
    int pe = pesToAvoid[i];
    peProxyPatchCounter[CkRankOf(pe)] = (1 << 20);    
  }
#endif
  // Avoid master Pe somewhat
  peProxyPatchCounter[CkRankOf(masterPe)] = 2; // patches.size();
  //--------------------------------------------------------
  for (int i=0;i < patches.size();i++) {
    PatchID pid = patches[i].patchID;
    int pe = findHomePatchPe(rankHomePatchIDs, pid);
    if (pe == -1) {
      // Patch not present on this node => try finding a ProxyPatch
      findProxyPatchPes(proxyPatchPes, pid);
      if (proxyPatchPes.size() == 0) {
        // No ProxyPatch => create one on rank that has the least ProxyPatches
        int rank = std::min_element(peProxyPatchCounter.begin(), peProxyPatchCounter.end()) - peProxyPatchCounter.begin();
        pe = CkNodeFirst(CkMyNode()) + rank;
        peProxyPatchCounter[rank]++;
      } else {
        // Choose ProxyPatch, try to avoid masterPe (current Pe) and Pes that already have a ProxyPatch,
        // this is done by finding the entry with minimum peProxyPatchCounter -value
        // Find miniumum among proxyPatchPes, i.e., find the minimum among
        // peProxyPatchCounter[CkRankOf(proxyPatchPes[j])]
        // int pppi = std::min_element(proxyPatchPes.begin(), proxyPatchPes.end(),
        //   [&](int i, int j) {return peProxyPatchCounter[CkRankOf(i)] < peProxyPatchCounter[CkRankOf(j)];})
        //   - proxyPatchPes.begin();
        // pe = proxyPatchPes[pppi];
        int minCounter = (1 << 30);
        for (int j=0;j < proxyPatchPes.size();j++) {
          if (minCounter > peProxyPatchCounter[CkRankOf(proxyPatchPes[j])]) {
            pe = proxyPatchPes[j];
            minCounter = peProxyPatchCounter[CkRankOf(pe)];
          }
        }
        if (pe == -1)
          NAMD_bug("CudaComputeNonbonded::assignPatches, Unable to choose PE with proxy patch");
        peProxyPatchCounter[CkRankOf(pe)]++;
      }
    } else if (std::find(pesToAvoid.begin(), pesToAvoid.end(), pe) != pesToAvoid.end()) {
      // Found home patch on this node, but it's on PE that should be avoided => find a new one
      int rank = std::min_element(peProxyPatchCounter.begin(), peProxyPatchCounter.end()) - peProxyPatchCounter.begin();
      pe = CkNodeFirst(CkMyNode()) + rank;
      peProxyPatchCounter[rank]++;
    }
    if (pe < CkNodeFirst(CkMyNode()) || pe >= CkNodeFirst(CkMyNode()) + CkMyNodeSize() )
      NAMD_bug("CudaComputeNonbonded::assignPatches, Invalid PE for a patch");
    rankPatches[CkRankOf(pe)].push_back(i);
    pidMap[pid] = i;
  }
  // Setup computes using pidMap
  for (int i=0;i < computes.size();i++) {
    computes[i].patchInd[0] = pidMap[computes[i].pid[0]];
    computes[i].patchInd[1] = pidMap[computes[i].pid[1]];
  }

  for (int i=0;i < CkMyNodeSize();i++) {
    if (rankPatches[i].size() > 0) pes.push_back(CkNodeFirst(CkMyNode()) + i);
  }
  computeMgr->sendAssignPatchesOnPe(pes, this);
  delete [] rankHomePatchIDs;
}

void CudaComputeNonbonded::initialize() {
  if (patches.size() > 0) {
    // Allocate CUDA version of patches
    cudaCheck(cudaSetDevice(deviceID));
    allocate_host<CudaPatchRecord>(&cudaPatches, patches.size());

    allocate_host<VirialEnergy>(&h_virialEnergy, 1);
    allocate_device<VirialEnergy>(&d_virialEnergy, 1);

#if CUDA_VERSION >= 5050
    int leastPriority, greatestPriority;
    cudaCheck(cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority));
    int priority = (doStreaming) ? leastPriority : greatestPriority;
    // int priority = greatestPriority;
    cudaCheck(cudaStreamCreateWithPriority(&stream,cudaStreamDefault, priority));
#else
    cudaCheck(cudaStreamCreate(&stream));
#endif
    cudaCheck(cudaEventCreate(&forceDoneEvent));

    buildExclusions();

    lock = CmiCreateLock();

    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  }  
}

//
// atomUpdate() can be called by any Pe
//
void CudaComputeNonbonded::atomUpdate() {
  atomsChangedIn = true;
}

//
// Compute patches[].atomStart, patches[].numAtoms, patches[].numFreeAtoms, and atomStorageSize
//
void CudaComputeNonbonded::updatePatches() {

  // Maximum number of tiles per tile list
  maxTileListLen = 0;
  int atomStart = 0;
  for (int i=0;i < patches.size();i++) {
    patches[i].atomStart = atomStart;
    cudaPatches[i].atomStart = atomStart;
    int numAtoms = patches[i].numAtoms;
    int numTiles = ((numAtoms-1)/WARPSIZE+1);
    maxTileListLen = std::max(maxTileListLen, numTiles);
    atomStart += numTiles*WARPSIZE;
  }
  atomStorageSize = atomStart;

  if (maxTileListLen >= 65536) {
    NAMD_bug("CudaComputeNonbonded::updatePatches, maximum number of tiles per tile lists (65536) blown");
  }
}

void CudaComputeNonbonded::skipPatch(int i) {
  if (CkMyPe() != patches[i].pe)
    NAMD_bug("CudaComputeNonbonded::skipPatch called on wrong Pe");
  Flags &flags = patches[i].patch->flags;
  patches[i].positionBox->skip();
  patches[i].forceBox->skip();
  if (flags.doGBIS) {
    patches[i].psiSumBox->skip();
    patches[i].intRadBox->skip();
    patches[i].bornRadBox->skip();
    patches[i].dEdaSumBox->skip();
    patches[i].dHdrPrefixBox->skip();
  }
}

void CudaComputeNonbonded::skipPatchesOnPe() {
  if (rankPatches[CkMyRank()].size() == 0)
    NAMD_bug("CudaComputeNonbonded::skipPatchesOnPe, empty rank");
  for (int i=0;i < rankPatches[CkMyRank()].size();i++) {
    skipPatch(rankPatches[CkMyRank()][i]);
  }
  bool done = false;
  CmiLock(lock);
  patchesCounter -= rankPatches[CkMyRank()].size();
  if (patchesCounter == 0) {
    patchesCounter = getNumPatches();
    done = true;
  }
  CmiUnlock(lock);
  if (done) {
    // Reduction must be done on masterPe
    computeMgr->sendFinishReductions(masterPe, this);
  }
}

void CudaComputeNonbonded::skip() {
  if (CkMyPe() != masterPe)
    NAMD_bug("CudaComputeNonbonded::skip() called on non masterPe");

  if (patches.size() == 0) return;

  doSkip = true;

  computeMgr->sendSkipPatchesOnPe(pes, this);
}

void CudaComputeNonbonded::getMaxMovementTolerance(float& maxAtomMovement, float& maxPatchTolerance) {
  if (CkMyPe() != masterPe)
    NAMD_bug("CudaComputeNonbonded::getMaxMovementTolerance() called on non masterPe");

  for (int i=0;i < patches.size();++i) {
    PatchRecord &pr = patches[i];

    float maxMove = pr.patch->flags.maxAtomMovement;
    if ( maxMove > maxAtomMovement ) maxAtomMovement = maxMove;

    float maxTol = pr.patch->flags.pairlistTolerance;
    if ( maxTol > maxPatchTolerance ) maxPatchTolerance = maxTol;
  }
}

inline void CudaComputeNonbonded::updateVdwTypesExclLoop(int first, int last, void *result, int paraNum, void *param) {
  CudaComputeNonbonded* c = (CudaComputeNonbonded *)param;
  c->updateVdwTypesExclSubset(first, last);
}

void CudaComputeNonbonded::updateVdwTypesExclSubset(int first, int last) {
  for (int i=first;i <= last;i++) {
    PatchRecord &pr = patches[i];
    int start = pr.atomStart;
    int numAtoms = pr.numAtoms;
    const CompAtom *compAtom = pr.compAtom;
    const CompAtomExt *compAtomExt = pr.patch->getCompAtomExtInfo();
    // Atoms have changed, re-do exclusions and vdw types
    int2* exclp = exclIndexMaxDiff + start;
    int* aip = atomIndex + start;
    for ( int k=0;k < numAtoms; ++k ) {
      int j = compAtomExt[k].sortOrder;
      vdwTypes[start + k] = compAtom[j].vdwType;
      aip[k] = compAtomExt[j].id;
#ifdef MEM_OPT_VERSION
      exclp[k].x = exclusionsByAtom[compAtomExt[j].exclId].y;
      exclp[k].y = exclusionsByAtom[compAtomExt[j].exclId].x;
#else // ! MEM_OPT_VERSION
      exclp[k].x = exclusionsByAtom[compAtomExt[j].id].y;
      exclp[k].y = exclusionsByAtom[compAtomExt[j].id].x;
#endif // MEM_OPT_VERSION
    }
  }
}

//
// Called every time atoms changed
//
void CudaComputeNonbonded::updateVdwTypesExcl() {
  // Re-allocate (VdwTypes, exclIndexMaxDiff) as needed
  reallocate_host<int>(&vdwTypes, &vdwTypesSize, atomStorageSize, 1.4f);
  reallocate_host<int2>(&exclIndexMaxDiff, &exclIndexMaxDiffSize, atomStorageSize, 1.4f);
  reallocate_host<int>(&atomIndex, &atomIndexSize, atomStorageSize, 1.4f);

#if CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if (useCkLoop >= 1) {
    CkLoop_Parallelize(updateVdwTypesExclLoop, 1, (void *)this, CkMyNodeSize(), 0, patches.size()-1);
  } else
#endif
  {
    updateVdwTypesExclSubset(0, patches.size()-1);
  }

  nonbondedKernel.updateVdwTypesExcl(atomStorageSize, vdwTypes, exclIndexMaxDiff, atomIndex, stream);
}

inline void CudaComputeNonbonded::copyAtomsLoop(int first, int last, void *result, int paraNum, void *param) {
  CudaComputeNonbonded* c = (CudaComputeNonbonded *)param;
  c->copyAtomsSubset(first, last);
}

void CudaComputeNonbonded::copyAtomsSubset(int first, int last) {
  for (int i=first;i <= last;++i) {
    PatchRecord &pr = patches[i];
    int numAtoms = pr.numAtoms;
    if (numAtoms > 0) {
      int start = pr.atomStart;
      const CudaAtom *src = pr.patch->getCudaAtomList();
      CudaAtom *dst = atoms + start;
      memcpy(dst, src, sizeof(CudaAtom)*numAtoms);
      // Fill the rest with the copy of the last atom
      int numAtomsAlign = ((numAtoms-1)/32+1)*32;
      CudaAtom lastAtom = src[numAtoms-1];
      for (int j=numAtoms;j < numAtomsAlign;j++) {
        dst[j] = lastAtom;
      }
    }
  }
}

void CudaComputeNonbonded::copyGBISphase(int i) {
  if (CkMyPe() != patches[i].pe)
    NAMD_bug("CudaComputeNonbonded::copyGBISphase called on wrong Pe");
  PatchRecord &pr = patches[i];
  const CompAtomExt *aExt = pr.patch->getCompAtomExtInfo();
  if (gbisPhase == 1) {
    //Copy GBIS intRadius to Host
    if (atomsChanged) {
      float *intRad0 = intRad0H + pr.atomStart;
      float *intRadS = intRadSH + pr.atomStart;
      for (int k=0;k < pr.numAtoms;++k) {
        int j = aExt[k].sortOrder;
        intRad0[k] = pr.intRad[2*j+0];
        intRadS[k] = pr.intRad[2*j+1];
      }
    }
  } else if (gbisPhase == 2) {
    float *bornRad = bornRadH + pr.atomStart;
    for ( int k=0; k < pr.numAtoms; ++k ) {
      int j = aExt[k].sortOrder;
      bornRad[k] = pr.bornRad[j];
    }
  } else if (gbisPhase == 3) {
    float *dHdrPrefix = dHdrPrefixH + pr.atomStart;
    for ( int k=0; k < pr.numAtoms; ++k ) {
      int j = aExt[k].sortOrder;
      dHdrPrefix[k] = pr.dHdrPrefix[j];
    }
  } // end phases
}

void CudaComputeNonbonded::openBox(int i) {
  if (CkMyPe() != patches[i].pe)
    NAMD_bug("CudaComputeNonbonded::openBox called on wrong Pe");
  SimParameters *simParams = Node::Object()->simParameters;
  if (!simParams->GBISOn || gbisPhase == 1) {
    patches[i].compAtom = patches[i].positionBox->open();
    copyAtomsSubset(i, i);
  }
  if (simParams->GBISOn) {
    if (gbisPhase == 1) {
      patches[i].intRad     = patches[i].intRadBox->open();
      patches[i].psiSum     = patches[i].psiSumBox->open();
    } else if (gbisPhase == 2) {
      patches[i].bornRad    = patches[i].bornRadBox->open();
      patches[i].dEdaSum    = patches[i].dEdaSumBox->open();
    } else if (gbisPhase == 3) {
      patches[i].dHdrPrefix = patches[i].dHdrPrefixBox->open();
    }
    copyGBISphase(i);
  }
}

void CudaComputeNonbonded::messageEnqueueWork() {
  if (masterPe != CkMyPe())
    NAMD_bug("CudaComputeNonbonded::messageEnqueueWork() must be called from masterPe");
  WorkDistrib::messageEnqueueWork(this);
}

void CudaComputeNonbonded::openBoxesOnPe() {
  if (rankPatches[CkMyRank()].size() == 0)
    NAMD_bug("CudaComputeNonbonded::openBoxesOnPe, empty rank");
  for (int i=0;i < rankPatches[CkMyRank()].size();i++) {
    openBox(rankPatches[CkMyRank()][i]);
  }
  bool done = false;
  CmiLock(lock);
  patchesCounter -= rankPatches[CkMyRank()].size();
  if (patchesCounter == 0) {
    patchesCounter = getNumPatches();
    done = true;
  }
  CmiUnlock(lock);
  if (done) {
    computeMgr->sendLaunchWork(masterPe, this);
  }
}

int CudaComputeNonbonded::noWork() {
  // Simply enqueu doWork on masterPe and return "no work"
  computeMgr->sendMessageEnqueueWork(masterPe, this);
  return 1;
}

void CudaComputeNonbonded::reallocateArrays() {
  cudaCheck(cudaSetDevice(deviceID));
  SimParameters *simParams = Node::Object()->simParameters;

  // Re-allocate atoms
  reallocate_host<CudaAtom>(&atoms, &atomsSize, atomStorageSize, 1.4f);

  // Re-allocate forces
  if (doStreaming) {
    reallocate_host<float4>(&h_forces, &h_forcesSize, atomStorageSize, 1.4f, cudaHostAllocMapped);
    reallocate_host<float4>(&h_forcesSlow, &h_forcesSlowSize, atomStorageSize, 1.4f, cudaHostAllocMapped);
  } else {
    reallocate_host<float4>(&h_forces, &h_forcesSize, atomStorageSize, 1.4f);
    reallocate_host<float4>(&h_forcesSlow, &h_forcesSlowSize, atomStorageSize, 1.4f);
  }
  reallocate_device<float4>(&d_forces, &d_forcesSize, atomStorageSize, 1.4f);
  reallocate_device<float4>(&d_forcesSlow, &d_forcesSlowSize, atomStorageSize, 1.4f);  

  if (simParams->GBISOn) {
    reallocate_host<float>(&intRad0H, &intRad0HSize, atomStorageSize, 1.2f);
    reallocate_host<float>(&intRadSH, &intRadSHSize, atomStorageSize, 1.2f);
    reallocate_host<GBReal>(&psiSumH, &psiSumHSize, atomStorageSize, 1.2f);
    reallocate_host<GBReal>(&dEdaSumH, &dEdaSumHSize, atomStorageSize, 1.2f);
    reallocate_host<float>(&bornRadH, &bornRadHSize, atomStorageSize, 1.2f);
    reallocate_host<float>(&dHdrPrefixH, &dHdrPrefixHSize, atomStorageSize, 1.2f);
  }
}

void CudaComputeNonbonded::doWork() {
  if (CkMyPe() != masterPe)
    NAMD_bug("CudaComputeNonbonded::doWork() called on non masterPe");

  // Read value of atomsChangedIn, which is set in atomUpdate(), and reset it.
  // atomsChangedIn can be set to true by any Pe
  // atomsChanged can only be set by masterPe
  // This use of double varibles makes sure we don't have race condition
  atomsChanged = atomsChangedIn;
  atomsChangedIn = false;

  SimParameters *simParams = Node::Object()->simParameters;

  if (patches.size() == 0) return;  // No work do to

  // Take the flags from the first patch on this Pe
  // Flags &flags = patches[rankPatches[CkMyRank()][0]].patch->flags;
  Flags &flags = patches[0].patch->flags;

  doSlow = flags.doFullElectrostatics;
  doEnergy = flags.doEnergy;
  doVirial = flags.doVirial;

  if (flags.doNonbonded) {

    if (simParams->GBISOn) {
      gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...
    }

    if (!simParams->GBISOn || gbisPhase == 1) {
      if ( computesChanged ) {
        updateComputes();
      }
      if (atomsChanged) {
        // Re-calculate patch atom numbers and storage
        updatePatches();
        reSortDone = false;
      }
      reallocateArrays();
    }

    // Open boxes on Pes and launch work to masterPe
    computeMgr->sendOpenBoxesOnPe(pes, this);

  } else {
    // No work to do, skip
    skip();
  }

}

void CudaComputeNonbonded::launchWork() {
  if (CkMyPe() != masterPe)
    NAMD_bug("CudaComputeNonbonded::launchWork() called on non masterPe");

  beforeForceCompute = CkWallTimer();

  cudaCheck(cudaSetDevice(deviceID));
  SimParameters *simParams = Node::Object()->simParameters;

  //execute only during GBIS phase 1, or if not using GBIS
  if (!simParams->GBISOn || gbisPhase == 1) {

    if ( atomsChanged || computesChanged ) {
      // Invalidate pair lists
      pairlistsValid = false;
      pairlistTolerance = 0.0f;
    }

    // Get maximum atom movement and patch tolerance
    float maxAtomMovement = 0.0f;
    float maxPatchTolerance = 0.0f;
    getMaxMovementTolerance(maxAtomMovement, maxPatchTolerance);
    // Update pair-list cutoff
    Flags &flags = patches[0].patch->flags;
    savePairlists = false;
    usePairlists = false;
    if ( flags.savePairlists ) {
      savePairlists = true;
      usePairlists = true;
    } else if ( flags.usePairlists ) {
      if ( ! pairlistsValid ||
           ( 2. * maxAtomMovement > pairlistTolerance ) ) {
        reduction->item(REDUCTION_PAIRLIST_WARNINGS) += 1;
      } else {
        usePairlists = true;
      }
    }
    if ( ! usePairlists ) {
      pairlistsValid = false;
    }
    float plcutoff = cutoff;
    if ( savePairlists ) {
      pairlistsValid = true;
      pairlistTolerance = 2. * maxPatchTolerance;
      plcutoff += pairlistTolerance;
    }
    plcutoff2 = plcutoff * plcutoff;

    // if (atomsChanged)
    //   CkPrintf("plcutoff = %f  listTolerance = %f  save = %d  use = %d\n",
    //     plcutoff, pairlistTolerance, savePairlists, usePairlists);

  } // if (!simParams->GBISOn || gbisPhase == 1)

  // Calculate PME & VdW forces
  if (!simParams->GBISOn || gbisPhase == 1) {
    doForce();
    if (doStreaming) {
      patchReadyQueue = nonbondedKernel.getPatchReadyQueue();
      patchReadyQueueLen = tileListKernel.getNumPatches();
      patchReadyQueueNext = 0;
      // Fill in empty patches [0 ... patchReadyQueueNext-1] at the top
      int numEmptyPatches = tileListKernel.getNumEmptyPatches();
      int* emptyPatches = tileListKernel.getEmptyPatches();
      for (int i=0;i < numEmptyPatches;i++) {
        patchReadyQueue[i] = emptyPatches[i];
      }
      if (patchReadyQueueLen != patches.size())
        NAMD_bug("CudaComputeNonbonded::launchWork, invalid patchReadyQueueLen");
    }
  }

  // For GBIS phase 1 at pairlist update, we must re-sort tile list
  // before calling doGBISphase1().
  if (atomsChanged && simParams->GBISOn && gbisPhase == 1) {
    // In this code path doGBISphase1() is called in forceDone()
    forceDoneSetCallback();
    return;
  }

  // GBIS Phases
  if (simParams->GBISOn) {
    if (gbisPhase == 1) {
      doGBISphase1();
    } else if (gbisPhase == 2) {
      doGBISphase2(); 
    } else if (gbisPhase == 3) {
      doGBISphase3(); 
    }
  }

  // Copy forces to host
  if (!simParams->GBISOn || gbisPhase == 3) {
    if (!doStreaming) {
      copy_DtoH<float4>(d_forces, h_forces, atomStorageSize, stream);
      if (doSlow) copy_DtoH<float4>(d_forcesSlow, h_forcesSlow, atomStorageSize, stream);
    }
  }

  if ((!simParams->GBISOn || gbisPhase == 2) && (doEnergy || doVirial)) {
    // For GBIS, energies are ready after phase 2
    nonbondedKernel.reduceVirialEnergy(tileListKernel,
      atomStorageSize, doEnergy, doVirial, doSlow, simParams->GBISOn,
      d_forces, d_forcesSlow, d_virialEnergy, stream);
    copy_DtoH<VirialEnergy>(d_virialEnergy, h_virialEnergy, 1, stream);
  }

  // Setup call back
  forceDoneSetCallback();
}

//
// GBIS Phase 1
//
void CudaComputeNonbonded::doGBISphase1() {
  cudaCheck(cudaSetDevice(deviceID));
  
  if (atomsChanged) {
    GBISKernel.updateIntRad(atomStorageSize, intRad0H, intRadSH, stream);
  }

  SimParameters *simParams = Node::Object()->simParameters;
  Lattice lattice = patches[0].patch->flags.lattice;

  GBISKernel.GBISphase1(tileListKernel, atomStorageSize,
    lattice.a().x, lattice.b().y, lattice.c().z,
    simParams->alpha_cutoff-simParams->fsMax, psiSumH, stream);
}

//
// GBIS Phase 2
//
void CudaComputeNonbonded::doGBISphase2() {
  cudaCheck(cudaSetDevice(deviceID));

  SimParameters *simParams = Node::Object()->simParameters;
  Lattice lattice = patches[0].patch->flags.lattice;

  GBISKernel.updateBornRad(atomStorageSize, bornRadH, stream);

  GBISKernel.GBISphase2(tileListKernel, atomStorageSize,
    doEnergy, doSlow,
    lattice.a().x, lattice.b().y, lattice.c().z,
    simParams->cutoff, simParams->nonbondedScaling, simParams->kappa,
    (simParams->switchingActive ? simParams->switchingDist : -1.0),
    simParams->dielectric, simParams->solvent_dielectric,
    d_forces, dEdaSumH, stream);
}

//
// GBIS Phase 3
//
void CudaComputeNonbonded::doGBISphase3() {
  cudaCheck(cudaSetDevice(deviceID));
  SimParameters *simParams = Node::Object()->simParameters;
  Lattice lattice = patches[0].patch->flags.lattice;

  if (doSlow) {
    GBISKernel.update_dHdrPrefix(atomStorageSize, dHdrPrefixH, stream);

    GBISKernel.GBISphase3(tileListKernel, atomStorageSize,
      lattice.a().x, lattice.b().y, lattice.c().z,
      simParams->alpha_cutoff-simParams->fsMax, d_forcesSlow, stream);
  }
}

//
// Calculate electrostatic & VdW forces
//
void CudaComputeNonbonded::doForce() {
  cudaCheck(cudaSetDevice(deviceID));

  Lattice lattice = patches[0].patch->flags.lattice;

  if (atomsChanged) {
    int numTileLists = calcNumTileLists();
    // Build initial tile lists and sort
    tileListKernel.buildTileLists(numTileLists, patches.size(), atomStorageSize,
      maxTileListLen,lattice.a().x, lattice.b().y, lattice.c().z,
      cudaPatches, (const float4*)atoms, plcutoff2, stream);
    // Prepare tile list for atom-based refinement
    tileListKernel.prepareTileList(stream);
  }

  if (atomsChanged) {
    // Update Vdw types and exclusion index & maxdiff
    updateVdwTypesExcl();
  }

  beforeForceCompute = CkWallTimer();

  // Calculate forces (and refine tile list if atomsChanged=true)
  nonbondedKernel.nonbondedForce(tileListKernel, atomStorageSize, atomsChanged,
    doEnergy, doVirial, doSlow,
    lattice.a().x, lattice.b().y, lattice.c().z,
    (const float4*)atoms, cutoff2, d_forces, d_forcesSlow, h_forces, h_forcesSlow,
    stream);

  if (atomsChanged) {
    tileListKernel.finishTileList(stream);
  }

  traceUserBracketEvent(CUDA_DEBUG_EVENT, beforeForceCompute, CkWallTimer());
}

//
// Count an upper estimate for the number of tile lists
//
int CudaComputeNonbonded::calcNumTileLists() {
  int numTileLists = 0;
  for (int i=0;i < computes.size();i++) {
    int pi1 = computes[i].patchInd[0];
    int numAtoms1 = patches[pi1].numAtoms;
    int numTiles1 = (numAtoms1-1)/WARPSIZE+1;
    numTileLists += numTiles1;
  }
  return numTileLists;
}

//
// Finish & submit reductions
//
void CudaComputeNonbonded::finishReductions() {

  if (CkMyPe() != masterPe)
    NAMD_bug("CudaComputeNonbonded::finishReductions() called on non masterPe");
  
  // fprintf(stderr, "%d finishReductions doSkip %d doVirial %d doEnergy %d\n", CkMyPe(), doSkip, doVirial, doEnergy);

  if (!doSkip) {

    if (doStreaming && (doVirial || doEnergy)) {
      // For streaming kernels, we must wait for virials and forces to be copied back to CPU
      if (!forceDoneEventRecord)
        NAMD_bug("CudaComputeNonbonded::finishReductions, forceDoneEvent not being recorded");
      cudaCheck(cudaEventSynchronize(forceDoneEvent));
      forceDoneEventRecord = false;
    }

    if (doVirial) {
      Tensor virialTensor;
      virialTensor.xx = h_virialEnergy->virial[0];
      virialTensor.xy = h_virialEnergy->virial[1];
      virialTensor.xz = h_virialEnergy->virial[2];
      virialTensor.yx = h_virialEnergy->virial[3];
      virialTensor.yy = h_virialEnergy->virial[4];
      virialTensor.yz = h_virialEnergy->virial[5];
      virialTensor.zx = h_virialEnergy->virial[6];
      virialTensor.zy = h_virialEnergy->virial[7];
      virialTensor.zz = h_virialEnergy->virial[8];
      // fprintf(stderr, "virialTensor %lf %lf %lf\n", virialTensor.xx, virialTensor.xy, virialTensor.xz);
      // fprintf(stderr, "virialTensor %lf %lf %lf\n", virialTensor.yx, virialTensor.yy, virialTensor.yz);
      // fprintf(stderr, "virialTensor %lf %lf %lf\n", virialTensor.zx, virialTensor.zy, virialTensor.zz);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, virialTensor);
      if (doSlow) {
        Tensor virialTensor;
        virialTensor.xx = h_virialEnergy->virialSlow[0];
        virialTensor.xy = h_virialEnergy->virialSlow[1];
        virialTensor.xz = h_virialEnergy->virialSlow[2];
        virialTensor.yx = h_virialEnergy->virialSlow[3];
        virialTensor.yy = h_virialEnergy->virialSlow[4];
        virialTensor.yz = h_virialEnergy->virialSlow[5];
        virialTensor.zx = h_virialEnergy->virialSlow[6];
        virialTensor.zy = h_virialEnergy->virialSlow[7];
        virialTensor.zz = h_virialEnergy->virialSlow[8];
        // fprintf(stderr, "virialTensor (slow) %lf %lf %lf\n", virialTensor.xx, virialTensor.xy, virialTensor.xz);
        // fprintf(stderr, "virialTensor (slow) %lf %lf %lf\n", virialTensor.yx, virialTensor.yy, virialTensor.yz);
        // fprintf(stderr, "virialTensor (slow) %lf %lf %lf\n", virialTensor.zx, virialTensor.zy, virialTensor.zz);
        ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, virialTensor);
      }
    }
    if (doEnergy) {
      // if (doSlow)
      //   printf("energyElec %lf energySlow %lf energyGBIS %lf\n", h_virialEnergy->energyElec, h_virialEnergy->energySlow, h_virialEnergy->energyGBIS);
      SimParameters *simParams = Node::Object()->simParameters;
      reduction->item(REDUCTION_LJ_ENERGY)    += h_virialEnergy->energyVdw;
      reduction->item(REDUCTION_ELECT_ENERGY) += h_virialEnergy->energyElec + ((simParams->GBISOn) ? h_virialEnergy->energyGBIS : 0.0);
      // fprintf(stderr, "energyGBIS %lf\n", h_virialEnergy->energyGBIS);
      if (doSlow) reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += h_virialEnergy->energySlow;
      // fprintf(stderr, "h_virialEnergy->energyElec %lf\n", h_virialEnergy->energyElec);
    }

    reduction->item(REDUCTION_EXCLUSION_CHECKSUM_CUDA) += tileListKernel.getNumExcluded();
  }
  reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
  reduction->submit();

  // Reset flags
  doSkip = false;
  computesChanged = false;
}

//
// Finish a single patch
//
void CudaComputeNonbonded::finishPatch(int i) {
  if (CkMyPe() != patches[i].pe)
    NAMD_bug("CudaComputeNonbonded::finishPatch called on wrong Pe");

  PatchRecord &pr = patches[i];
  pr.results = pr.forceBox->open();

  const CompAtomExt *aExt = pr.patch->getCompAtomExtInfo();
  int atomStart = pr.atomStart;
  int numAtoms = pr.numAtoms;
  if (numAtoms > 0) {
    Force *f      = pr.results->f[Results::nbond];
    Force *f_slow = pr.results->f[Results::slow];
    float4 *af      = h_forces + atomStart;
    float4 *af_slow = h_forcesSlow + atomStart;
    // float maxf = 0.0f;
    // int maxf_k;
    for ( int k=0; k<numAtoms; ++k ) {
      int j = aExt[k].sortOrder;
      f[j].x += af[k].x;
      f[j].y += af[k].y;
      f[j].z += af[k].z;
      // if (maxf < fabsf(af[k].x) || maxf < fabsf(af[k].y) || maxf < fabsf(af[k].z)) {
      //   maxf = std::max(maxf, fabsf(af[k].x));
      //   maxf = std::max(maxf, fabsf(af[k].y));
      //   maxf = std::max(maxf, fabsf(af[k].z));
      //   maxf_k = k;
      // }
      if ( doSlow ) {
        f_slow[j].x += af_slow[k].x;
        f_slow[j].y += af_slow[k].y;
        f_slow[j].z += af_slow[k].z;
      }
    }
    // if (maxf > 10000.0f) {
    //   fprintf(stderr, "%d %f %f %f\n", maxf_k, af[maxf_k].x, af[maxf_k].y, af[maxf_k].z);
    //   cudaCheck(cudaStreamSynchronize(stream));
    //   NAMD_die("maxf!");
    // }
  }

  pr.positionBox->close(&(pr.compAtom));
  pr.forceBox->close(&(pr.results));
}

//
// Finish a set of patches on this pe
//
void CudaComputeNonbonded::finishSetOfPatchesOnPe(std::vector<int>& patchSet) {
  if (patchSet.size() == 0)
    NAMD_bug("CudaComputeNonbonded::finishPatchesOnPe, empty rank");
  SimParameters *simParams = Node::Object()->simParameters;
  // Save value of gbisPhase here because it can change after the last finishGBISPhase() or finishPatch() is called
  int gbisPhaseSave = gbisPhase;
  // Close Boxes depending on Phase
  if (simParams->GBISOn) {
    for (int i=0;i < patchSet.size();i++) {
      finishGBISPhase(patchSet[i]);
    }
  }
  // Finish patches
  if (!simParams->GBISOn || gbisPhaseSave == 3) {
    for (int i=0;i < patchSet.size();i++) {
      finishPatch(patchSet[i]);
    }
  }
  bool done = false;
  CmiLock(lock);
  patchesCounter -= patchSet.size();
  if (patchesCounter == 0) {
    patchesCounter = getNumPatches();
    done = true;
  }
  CmiUnlock(lock);
  if (done) {
    // Do reductions
    if (!simParams->GBISOn || gbisPhaseSave == 3) {
      // Reduction must be done on masterPe
      computeMgr->sendFinishReductions(masterPe, this);
    }
  }
}

//
// Finish all patches that are on this pe
//
void CudaComputeNonbonded::finishPatchesOnPe() {
  finishSetOfPatchesOnPe(rankPatches[CkMyRank()]);
}

//
// Finish single patch on this pe
//
void CudaComputeNonbonded::finishPatchOnPe(int i) {
  std::vector<int> v(1, i);
  finishSetOfPatchesOnPe(v);
}

void CudaComputeNonbonded::finishPatches() {
  computeMgr->sendFinishPatchesOnPe(pes, this);
}

void CudaComputeNonbonded::finishGBISPhase(int i) {
  if (CkMyPe() != patches[i].pe)
    NAMD_bug("CudaComputeNonbonded::finishGBISPhase called on wrong Pe");
  PatchRecord &pr = patches[i];
  const CompAtomExt *aExt = pr.patch->getCompAtomExtInfo();
  int atomStart = pr.atomStart;
  if (gbisPhase == 1) {
    GBReal *psiSumMaster = psiSumH + atomStart;
    for ( int k=0; k<pr.numAtoms; ++k ) {
      int j = aExt[k].sortOrder;
      pr.psiSum[j] += psiSumMaster[k];
    }
    pr.psiSumBox->close(&(pr.psiSum));
  } else if (gbisPhase == 2) {
    GBReal *dEdaSumMaster = dEdaSumH + atomStart;
    for ( int k=0; k<pr.numAtoms; ++k ) {
      int j = aExt[k].sortOrder;
      pr.dEdaSum[j] += dEdaSumMaster[k];
    }
    pr.dEdaSumBox->close(&(pr.dEdaSum));
  } else if (gbisPhase == 3) {
    pr.intRadBox->close(&(pr.intRad)); //box 6
    pr.bornRadBox->close(&(pr.bornRad)); //box 7
    pr.dHdrPrefixBox->close(&(pr.dHdrPrefix)); //box 9
  } //end phases
}

void CudaComputeNonbonded::finishTimers() {
  SimParameters *simParams = Node::Object()->simParameters;

  if (simParams->GBISOn) {
    if (gbisPhase == 1)
      traceUserBracketEvent(CUDA_GBIS1_KERNEL, beforeForceCompute, CkWallTimer());
    if (gbisPhase == 2)
      traceUserBracketEvent(CUDA_GBIS2_KERNEL, beforeForceCompute, CkWallTimer());
    if (gbisPhase == 3)
      traceUserBracketEvent(CUDA_GBIS3_KERNEL, beforeForceCompute, CkWallTimer());
  } else {
    traceUserBracketEvent(CUDA_VDW_KERNEL, beforeForceCompute, CkWallTimer());
  }  
}

//
// Re-sort tile lists if neccessary
//
void CudaComputeNonbonded::reSortTileLists() {
  // Re-sort tile lists
  SimParameters *simParams = Node::Object()->simParameters;
  cudaCheck(cudaSetDevice(deviceID));
  tileListKernel.reSortTileLists(simParams->GBISOn, stream);
}

void CudaComputeNonbonded::forceDoneCheck(void *arg, double walltime) {
  CudaComputeNonbonded* c = (CudaComputeNonbonded *)arg;

  if (CkMyPe() != c->masterPe)
    NAMD_bug("CudaComputeNonbonded::forceDoneCheck called on non masterPe");

  SimParameters *simParams = Node::Object()->simParameters;
  cudaCheck(cudaSetDevice(c->deviceID));

  if (c->doStreaming) {
    int patchInd;
    while ( -1 != (patchInd = c->patchReadyQueue[c->patchReadyQueueNext]) ) {
      c->patchReadyQueue[c->patchReadyQueueNext] = -1;
      c->patchReadyQueueNext++;
      c->checkCount = 0;

      if ( c->patchReadyQueueNext == c->patchReadyQueueLen ) {
        c->finishTimers();
        if (c->atomsChanged && (!simParams->GBISOn || c->gbisPhase == 1) && !c->reSortDone) {
          c->reSortTileLists();
          c->reSortDone = true;
          if (simParams->GBISOn && c->gbisPhase == 1) {
            // We must do GBIS Phase 1
            c->doGBISphase1();
            c->forceDoneSetCallback();
            return;
          }
        }
      }

      // Finish patch
      int pe = c->patches[patchInd].pe;
      c->computeMgr->sendFinishPatchOnPe(pe, c, patchInd);

      // Last patch, return
      if ( c->patchReadyQueueNext == c->patchReadyQueueLen ) return;

    }
  } else {
    if (!c->forceDoneEventRecord)
      NAMD_bug("CudaComputeNonbonded::forceDoneCheck, forceDoneEvent not being recorded");
    cudaError_t err = cudaEventQuery(c->forceDoneEvent);
    if (err == cudaSuccess) {
      // Event has occurred
      c->forceDoneEventRecord = false;
      c->checkCount = 0;
      c->finishTimers();
      if (c->atomsChanged && (!simParams->GBISOn || c->gbisPhase == 1) && !c->reSortDone) {
        c->reSortTileLists();
        c->reSortDone = true;
        if (simParams->GBISOn && c->gbisPhase == 1) {
          // We must do GBIS Phase 1
          c->doGBISphase1();
          c->forceDoneSetCallback();
          return;
        }
      }
      c->finishPatches();
      return;
    } else if (err != cudaErrorNotReady) {
      // Anything else is an error
      cudaCheck(err);
      // NAMD_bug("CudaComputeNonbonded::forceDoneCheck, cudaEventQuery returned error");
    }
  }

  // if (c->checkCount % 1000 == 0)
  //   fprintf(stderr, "c->patchReadyQueueNext %d\n", c->patchReadyQueueNext);

  // Event has not occurred
  c->checkCount++;
  if (c->checkCount >= 1000000) {
    NAMD_bug("CudaComputeNonbonded::forceDoneCheck, check count exceeded");
  }

  // Call again 
  CcdCallBacksReset(0, walltime);
  CcdCallFnAfter(forceDoneCheck, arg, 0.1);
}

//
// Set call back for all the work in the stream at this point
//
void CudaComputeNonbonded::forceDoneSetCallback() {
  if (CkMyPe() != masterPe)
    NAMD_bug("CudaComputeNonbonded::forceDoneSetCallback called on non masterPe");
  beforeForceCompute = CkWallTimer();
  cudaCheck(cudaSetDevice(deviceID));
  if (!doStreaming || doVirial || doEnergy) {
    cudaCheck(cudaEventRecord(forceDoneEvent, stream));
    forceDoneEventRecord = true;
  }
  checkCount = 0;
  CcdCallBacksReset(0, CmiWallTimer());
  // Set the call back at 0.1ms
  CcdCallFnAfter(forceDoneCheck, this, 0.1);
}

struct cr_sortop_distance {
  const Lattice &l;
  cr_sortop_distance(const Lattice &lattice) : l(lattice) { }
  bool operator() (CudaComputeNonbonded::ComputeRecord i,
      CudaComputeNonbonded::ComputeRecord j) {
    Vector a = l.a();
    Vector b = l.b();
    Vector c = l.c();
    BigReal ri = (i.offset.x * a + i.offset.y * b + i.offset.z * c).length2();
    BigReal rj = (j.offset.x * a + j.offset.y * b + j.offset.z * c).length2();
    return ( ri < rj );
  }
};

static inline bool sortop_bitreverse(int a, int b) {
  if ( a == b ) return 0; 
  for ( int bit = 1; bit; bit *= 2 ) {
    if ( (a&bit) != (b&bit) ) return ((a&bit) < (b&bit));
  }
  return 0;
}

struct cr_sortop_reverse_priority {
  cr_sortop_distance &distop;
  const CudaComputeNonbonded::PatchRecord *pr;
  cr_sortop_reverse_priority(cr_sortop_distance &sod,
       const CudaComputeNonbonded::PatchRecord *patchrecs) : distop(sod), pr(patchrecs) { }
  bool pid_compare_priority(int2 pidi, int2 pidj) {
    const CudaComputeNonbonded::PatchRecord &pri = pr[pidi.y];
    const CudaComputeNonbonded::PatchRecord &prj = pr[pidj.y];
    if ( pri.isSamePhysicalNode && ! prj.isSamePhysicalNode ) return 0;
    if ( prj.isSamePhysicalNode && ! pri.isSamePhysicalNode ) return 1;
    if ( pri.isSameNode && ! prj.isSameNode ) return 0;
    if ( prj.isSameNode && ! pri.isSameNode ) return 1;
    if ( pri.isSameNode ) {  // and prj.isSameNode
      int rpri = pri.reversePriorityRankInPe;
      int rprj = prj.reversePriorityRankInPe;
      if ( rpri != rprj ) return rpri > rprj;
      return sortop_bitreverse(CkRankOf(pri.pe),CkRankOf(prj.pe));
    }
    int ppi = PATCH_PRIORITY(pidi.x);
    int ppj = PATCH_PRIORITY(pidj.x);
    if ( ppi != ppj ) return ppi < ppj;
    return pidi.x < pidj.x;
  }
  bool operator() (CudaComputeNonbonded::ComputeRecord j,
      CudaComputeNonbonded::ComputeRecord i) {  // i and j reversed
    // Choose patch i (= patch with greater priority)
    int2 pidi = pid_compare_priority(make_int2(i.pid[0], i.patchInd[0]), make_int2(i.pid[1], i.patchInd[1])) ? make_int2(i.pid[0], i.patchInd[0]) : make_int2(i.pid[1], i.patchInd[1]);
    // Choose patch j
    int2 pidj = pid_compare_priority(make_int2(j.pid[0], j.patchInd[0]), make_int2(j.pid[1], j.patchInd[1])) ? make_int2(j.pid[0], j.patchInd[0]) : make_int2(j.pid[1], j.patchInd[1]);
    if ( pidi.x != pidj.x ) return pid_compare_priority(pidi, pidj);
    return distop(i,j);
  }
};

//
// Setup computes. This is only done at the beginning and at load balancing, hence the lack of
// consideration for performance in the CPU->GPU memory copy.
//
void CudaComputeNonbonded::updateComputes() {
  cudaCheck(cudaSetDevice(deviceID));

  Lattice lattice = patches[0].patch->flags.lattice;
  cr_sortop_distance so(lattice);
  std::stable_sort(computes.begin(), computes.end(), so);

  if (doStreaming) {
    cr_sortop_reverse_priority sorp(so, patches.data());
    std::stable_sort(computes.begin(), computes.end(), sorp);
  }

  CudaComputeRecord* cudaComputes = new CudaComputeRecord[computes.size()];

  for (int i=0;i < computes.size();i++) {
    cudaComputes[i].patchInd.x = computes[i].patchInd[0];
    cudaComputes[i].patchInd.y = computes[i].patchInd[1];
    cudaComputes[i].offsetXYZ.x = computes[i].offset.x;
    cudaComputes[i].offsetXYZ.y = computes[i].offset.y;
    cudaComputes[i].offsetXYZ.z = computes[i].offset.z;
  }

  tileListKernel.updateComputes(computes.size(), cudaComputes, stream);
  cudaCheck(cudaStreamSynchronize(stream));

  delete [] cudaComputes;
}

struct exlist_sortop {
  bool operator() (int32 *li, int32 *lj) {
    return ( li[1] < lj[1] );
  }
};

//
// Builds the exclusions table. Swiped from ComputeNonbondedCUDA.C
//
void CudaComputeNonbonded::buildExclusions() {
  cudaCheck(cudaSetDevice(deviceID));

  Molecule *mol = Node::Object()->molecule;

#ifdef MEM_OPT_VERSION
  int natoms = mol->exclSigPoolSize;
#else
  int natoms = mol->numAtoms; 
#endif

	if (exclusionsByAtom != NULL) delete [] exclusionsByAtom;
  exclusionsByAtom = new int2[natoms];

  // create unique sorted lists

  ObjectArena<int32> listArena;
  ResizeArray<int32*> unique_lists;
  int32 **listsByAtom = new int32*[natoms];
  SortableResizeArray<int32> curList;
  for ( int i=0; i<natoms; ++i ) {
    curList.resize(0);
    curList.add(0);  // always excluded from self
#ifdef MEM_OPT_VERSION
    const ExclusionSignature *sig = mol->exclSigPool + i;
    int n = sig->fullExclCnt;
    for ( int j=0; j<n; ++j ) { curList.add(sig->fullOffset[j]); }
    n += 1;
#else
    const int32 *mol_list = mol->get_full_exclusions_for_atom(i);
    int n = mol_list[0] + 1;
    for ( int j=1; j<n; ++j ) {
      curList.add(mol_list[j] - i);
    }
#endif
    curList.sort();

    int j;
    for ( j=0; j<unique_lists.size(); ++j ) {
      if ( n != unique_lists[j][0] ) continue;  // no match
      int k;
      for ( k=0; k<n; ++k ) {
        if ( unique_lists[j][k+3] != curList[k] ) break;
      }
      if ( k == n ) break;  // found match
    }
    if ( j == unique_lists.size() ) {  // no match
      int32 *list = listArena.getNewArray(n+3);
      list[0] = n;
      int maxdiff = 0;
      maxdiff = -1 * curList[0];
      if ( curList[n-1] > maxdiff ) maxdiff = curList[n-1];
      list[1] = maxdiff;
      for ( int k=0; k<n; ++k ) {
        list[k+3] = curList[k];
      }
      unique_lists.add(list);
    }
    listsByAtom[i] = unique_lists[j];
  }
  // sort lists by maxdiff
  std::stable_sort(unique_lists.begin(), unique_lists.end(), exlist_sortop());
  long int totalbits = 0;
  int nlists = unique_lists.size();
  for ( int j=0; j<nlists; ++j ) {
    int32 *list = unique_lists[j];
    int maxdiff = list[1];
    list[2] = totalbits + maxdiff;
    totalbits += 2*maxdiff + 1;
  }
  for ( int i=0; i<natoms; ++i ) {
    exclusionsByAtom[i].x = listsByAtom[i][1];  // maxdiff
    exclusionsByAtom[i].y = listsByAtom[i][2];  // start
  }
  delete [] listsByAtom;

  if ( totalbits & 31 ) totalbits += ( 32 - ( totalbits & 31 ) );

  {
    long int bytesneeded = totalbits / 8;
    if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Found %d unique exclusion lists needing %ld bytes\n",
		unique_lists.size(), bytesneeded);
    }

    long int bytesavail = MAX_EXCLUSIONS * sizeof(unsigned int);
    if ( bytesneeded > bytesavail ) {
      char errmsg[512];
      sprintf(errmsg,"Found %d unique exclusion lists needing %ld bytes "
                     "but only %ld bytes can be addressed with 32-bit int.",
                     unique_lists.size(), bytesneeded, bytesavail);
      NAMD_die(errmsg);
    }
  }

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

  unsigned int *exclusion_bits = new unsigned int[totalbits/32];
  memset(exclusion_bits, 0, totalbits/8);

  long int base = 0;
  for ( int i=0; i<unique_lists.size(); ++i ) {
    base += unique_lists[i][1];
    if ( unique_lists[i][2] != (int32)base ) {
      NAMD_bug("CudaComputeNonbonded::build_exclusions base != stored");
    }
    int n = unique_lists[i][0];
    for ( int j=0; j<n; ++j ) {
      SET_EXCL(exclusion_bits,base,unique_lists[i][j+3]);
    }
    base += unique_lists[i][1] + 1;
  }

  int numExclusions = totalbits/32;

	nonbondedKernel.bindExclusions(numExclusions, exclusion_bits);

  delete [] exclusion_bits;
}

void CudaComputeNonbonded::buildTables() {
  buildForceAndEnergyTables(4096);
  buildVdwCoefTable();
}

//
// Builds the VdW Lennard-Jones coefficient table. Swiped from ComputeNonbondedCUDA.C
// NOTE: should only be called once
//
void CudaComputeNonbonded::buildVdwCoefTable() {
  cudaCheck(cudaSetDevice(deviceID));

  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  const int dim = ljTable->get_table_dim();

  // round dim up to odd multiple of 16
  int tsize = (((dim+16+31)/32)*32)-16;
  if ( tsize < dim ) NAMD_bug("CudaComputeNonbonded::buildVdwCoefTable bad tsize");

  float2 *h_vdwCoefTable = new float2[tsize*tsize];
  float2 *row = h_vdwCoefTable;
  for ( int i=0; i<dim; ++i, row += tsize ) {
    for ( int j=0; j<dim; ++j ) {
      const LJTable::TableEntry *e = ljTable->table_val(i,j);
      row[j].x = e->A * scaling;
      row[j].y = e->B * scaling;
    }
  }

	nonbondedKernel.bindVdwCoefTable(h_vdwCoefTable, tsize);

  delete [] h_vdwCoefTable;

  if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Updated CUDA LJ table with %d x %d elements.\n", dim, dim);
  }
}

//
// Builds force and energy tables. Swiped from ComputeNonbondedCUDA.C
//
void CudaComputeNonbonded::buildForceAndEnergyTables(int tableSize) {
  cudaCheck(cudaSetDevice(deviceID));

  float4* t = new float4[tableSize];
  float4* et = new float4[tableSize];  // energy table

  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  // const int r2_delta_expc = 64 * (r2_delta_exp - 127);
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  double* r2list = new double[tableSize];  // double to match cpu code
  for ( int i=1; i<tableSize; ++i ) {
    double r = ((double) tableSize) / ( (double) i + 0.5 );
    r2list[i] = r*r + r2_delta;
  }

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  for ( int i=1; i<tableSize; ++i ) {
    double r = ((double) tableSize) / ( (double) i + 0.5 );
    int table_i = (r2iilist[2*i] >> 14) + r2_delta_expc;  // table_i >= 0

    if ( r > cutoff ) {
      t[i].x = 0.;
      t[i].y = 0.;
      t[i].z = 0.;
      t[i].w = 0.;
      et[i].x = 0.;
      et[i].y = 0.;
      et[i].z = 0.;
      et[i].w = 0.;
      continue;
    }

    BigReal diffa = r2list[i] - r2_table[table_i];

    // coulomb 1/r or fast force
    // t[i].x = 1. / (r2 * r);  // -1/r * d/dr r^-1
    {
      BigReal table_a = fast_table[4*table_i];
      BigReal table_b = fast_table[4*table_i+1];
      BigReal table_c = fast_table[4*table_i+2];
      BigReal table_d = fast_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].x = 2. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].x = ener;
    }


    // pme correction for slow force
    // t[i].w = 0.;
    {
      BigReal table_a = scor_table[4*table_i];
      BigReal table_b = scor_table[4*table_i+1];
      BigReal table_c = scor_table[4*table_i+2];
      BigReal table_d = scor_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].w = 2. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].w = ener;
    }


    // vdw 1/r^6
    // t[i].y = 6. / (r8);  // -1/r * d/dr r^-6
    {
      BigReal table_a = vdwb_table[4*table_i];
      BigReal table_b = vdwb_table[4*table_i+1];
      BigReal table_c = vdwb_table[4*table_i+2];
      BigReal table_d = vdwb_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].y = 2. * -1. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].y = -1. * ener;
    }


    // vdw 1/r^12
    // t[i].z = 12e / (r8 * r4 * r2);  // -1/r * d/dr r^-12
    {
      BigReal table_a = vdwa_table[4*table_i];
      BigReal table_b = vdwa_table[4*table_i+1];
      BigReal table_c = vdwa_table[4*table_i+2];
      BigReal table_d = vdwa_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].z = 2. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].z = ener;
    }

  }

  t[0].x = 0.f;
  t[0].y = 0.f;
  t[0].z = 0.f;
  t[0].w = 0.f;
  et[0].x = et[1].x;
  et[0].y = et[1].y;
  et[0].z = et[1].z;
  et[0].w = et[1].w;

  nonbondedKernel.bindForceAndEnergyTable(tableSize, t, et);

  if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Updated CUDA force table with %d elements.\n", tableSize);
  }

  delete [] t;
  delete [] et;
  delete [] r2list;
}

#endif // NAMD_CUDA
