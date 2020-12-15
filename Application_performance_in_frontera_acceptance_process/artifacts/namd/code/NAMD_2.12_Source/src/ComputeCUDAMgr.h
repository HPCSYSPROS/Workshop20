#ifndef COMPUTECUDAMGR_H
#define COMPUTECUDAMGR_H
#include <vector>
#include "CudaUtils.h"
#include "ComputeCUDAMgr.decl.h"
#include "CudaComputeNonbonded.h"
#ifdef NAMD_CUDA

class ComputeCUDAMgr : public CBase_ComputeCUDAMgr {
public:
  // ComputeCUDAMgr_SDAG_CODE;
	ComputeCUDAMgr();
	ComputeCUDAMgr(CkMigrateMessage *);
	~ComputeCUDAMgr();
	void initialize(CkQdMsg *msg);
	void initialize_devices(CkQdMsg *msg);
  static ComputeCUDAMgr* getComputeCUDAMgr();
  CudaComputeNonbonded* createCudaComputeNonbonded(ComputeID c);
  CudaComputeNonbonded* getCudaComputeNonbonded();
private:

  // Number of CUDA devices on this node that are used in computation
  int numDevices;
  std::vector<CudaComputeNonbonded*> cudaComputeNonbondedList;
};

#else // NAMD_CUDA

class ComputeCUDAMgr : public CBase_ComputeCUDAMgr {
};

#endif // NAMD_CUDA
#endif // COMPUTECUDAMGR_H