#include "NamdTypes.h"
#include "common.h"
#include "Node.h"
#include "ComputeCUDAMgr.h"

#include "DeviceCUDA.h"
#ifdef NAMD_CUDA
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;

//
// Class constructor
//
ComputeCUDAMgr::ComputeCUDAMgr() {
	// __sdag_init();
  numDevices = 0;
  // numNodesContributed = 0;
  // numDevicesMax = 0;
}

//
// Class constructor
//
ComputeCUDAMgr::ComputeCUDAMgr(CkMigrateMessage *) {
	// __sdag_init();
  NAMD_bug("ComputeCUDAMgr cannot be migrated");
  numDevices = 0;
  // numNodesContributed = 0;
  // numDevicesMax = 0;
}

//
// Class destructor
//
ComputeCUDAMgr::~ComputeCUDAMgr() {
  for (int i=0;i < numDevices;i++) {
    if (cudaComputeNonbondedList[i] != NULL) delete cudaComputeNonbondedList[i];
  }
}

//
// Initialize manager
// This gets called on rank 0 of each node
//
void ComputeCUDAMgr::initialize(CkQdMsg *msg) {
	if (msg != NULL) delete msg;

	numDevices = deviceCUDA->getDeviceCount();

  // Create pointers to devices
  cudaComputeNonbondedList.resize(numDevices, NULL);
}

ComputeCUDAMgr* ComputeCUDAMgr::getComputeCUDAMgr() {
  // Get pointer to ComputeCUDAMgr on this node
  CProxy_ComputeCUDAMgr computeCUDAMgrProxy = CkpvAccess(BOCclass_group).computeCUDAMgr;
  ComputeCUDAMgr* computeCUDAMgr = computeCUDAMgrProxy.ckLocalBranch();
  if (computeCUDAMgr == NULL)
    NAMD_bug("getComputeCUDAMgr, unable to locate local branch of BOC entry ComputeCUDAMgr");
  return computeCUDAMgr;
}

//
// Creates CudaComputeNonbonded object
//
CudaComputeNonbonded* ComputeCUDAMgr::createCudaComputeNonbonded(ComputeID c) {
  int deviceID = deviceCUDA->getDeviceID();
  if (cudaComputeNonbondedList.at(deviceID) != NULL)
    NAMD_bug("ComputeCUDAMgr::createCudaComputeNonbonded called twice");
  bool doStreaming = !deviceCUDA->getNoStreaming() && !Node::Object()->simParameters->GBISOn;
  cudaComputeNonbondedList[deviceID] = new CudaComputeNonbonded(c, deviceID, doStreaming);
  cudaComputeNonbondedList[deviceID]->buildTables();
  return cudaComputeNonbondedList[deviceID];
}

//
// Returns CudaComputeNonbonded for this Pe
//
CudaComputeNonbonded* ComputeCUDAMgr::getCudaComputeNonbonded() {
  // Get device ID for this Pe
  int deviceID = deviceCUDA->getDeviceID();
  CudaComputeNonbonded* p = cudaComputeNonbondedList[deviceID];
  if (p == NULL)
    NAMD_bug("ComputeCUDAMgr::getCudaComputeNonbonded(), device not created yet");
  return p;
}
#endif // NAMD_CUDA

#include "ComputeCUDAMgr.def.h"
