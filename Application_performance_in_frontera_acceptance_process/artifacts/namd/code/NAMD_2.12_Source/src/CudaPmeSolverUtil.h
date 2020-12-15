#ifndef CUDAPMESOLVERUTIL_H
#define CUDAPMESOLVERUTIL_H
#include <stdio.h>
#ifdef NAMD_CUDA
#include <cuda.h>
#include <cufft.h>
#endif // NAMD_CUDA
#include "PmeSolverUtil.h"
#include "CudaUtils.h"
#include "CudaPmeSolverUtilKernel.h"

#ifdef NAMD_CUDA
void writeComplexToDisk(const float2 *d_data, const int size, const char* filename, cudaStream_t stream);
void writeHostComplexToDisk(const float2 *h_data, const int size, const char* filename);
void writeRealToDisk(const float *d_data, const int size, const char* filename, cudaStream_t stream);

#define cufftCheck(stmt) do {						\
  cufftResult err = stmt;						\
  if (err != CUFFT_SUCCESS) {						\
  	char msg[128];	\
	  sprintf(msg, "%s in file %s, function %s\n", #stmt,__FILE__,__FUNCTION__); \
	  cudaDie(msg); \
  }									\
} while(0)

//
// CUDA implementation of FFTCompute
//
class CudaFFTCompute : public FFTCompute {
private:
  cufftHandle forwardPlan, backwardPlan;
  cufftType_t forwardType, backwardType;
  int deviceID;
	cudaStream_t stream;
	void setStream();

private:
	float* allocateData(const int dataSizeRequired);
	void plan3D(int *n, int flags);
	void plan2D(int *n, int howmany, int flags);
	void plan1DX(int *n, int howmany, int flags);
	void plan1DY(int *n, int howmany, int flags);
	void plan1DZ(int *n, int howmany, int flags);
	// int ncall, plantype;

public:
	CudaFFTCompute(int deviceID, cudaStream_t stream) : deviceID(deviceID), stream(stream) {}
	~CudaFFTCompute();
	void forward();
	void backward();
};

//
// Cuda implementation of PmeKSpaceCompute class
//
class CudaPmePencilXYZ;
class CudaPmePencilZ;

class CudaPmeKSpaceCompute : public PmeKSpaceCompute {
private:
	int cuda_arch;
	int deviceID;
	cudaStream_t stream;
	// Device memory versions of (bm1, bm2, bm3)
	float *d_bm1, *d_bm2, *d_bm3;
	//float *prefac_x, *prefac_y, *prefac_z;
	struct EnergyVirial {
		double energy;
		double virial[9];
	};
	EnergyVirial* d_energyVirial;
	EnergyVirial* h_energyVirial;
	cudaEvent_t copyEnergyVirialEvent;
	bool ortho;
  // Check counter for event polling in energyAndVirialCheck()
  int checkCount;
	static void energyAndVirialCheck(void *arg, double walltime);
	CudaPmePencilXYZ* pencilXYZPtr;
	CudaPmePencilZ* pencilZPtr;
public:
	CudaPmeKSpaceCompute(PmeGrid pmeGrid, const int permutation,
		const int jblock, const int kblock, double kappa, int cuda_arch, int deviceID, cudaStream_t stream);
	~CudaPmeKSpaceCompute();
	void solve(Lattice &lattice, const bool doEnergy, const bool doVirial, float* data);
	// void waitEnergyAndVirial();
	double getEnergy();
	void getVirial(double *virial);
	void energyAndVirialSetCallback(CudaPmePencilXYZ* pencilPtr);
	void energyAndVirialSetCallback(CudaPmePencilZ* pencilPtr);
};

//
// Cuda implementation of PmeRealSpaceCompute class
//

class ComputePmeCUDADevice;

class CudaPmeRealSpaceCompute : public PmeRealSpaceCompute {
private:
	int cuda_arch;
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  bool gridTexObjActive;
  cudaTextureObject_t gridTexObj;
#endif
  int tex_data_len;
	float* tex_data;
	int deviceID;
	cudaStream_t stream;
	void setupGridTexture(float* data, int data_len);
	// Device memory for atoms
	int d_atomsCapacity;
	CudaAtom* d_atoms;
	// Device memory for patches
	// int d_patchesCapacity;
	// PatchInfo* d_patches;
	// Device memory for forces
	int d_forceCapacity;
	CudaForce* d_force;
	// // Device memory for self energy
	// double* d_selfEnergy;
  // Events
  cudaEvent_t gatherForceEvent;
  // Check counter for event polling
  int checkCount;
  // Store device pointer for event polling
  ComputePmeCUDADevice* devicePtr;
  static void cuda_gatherforce_check(void *arg, double walltime);
public:
	CudaPmeRealSpaceCompute(PmeGrid pmeGrid, const int jblock, const int kblock,
		int cuda_arch, int deviceID, cudaStream_t stream);
	~CudaPmeRealSpaceCompute();
	void copyAtoms(const int numAtoms, const CudaAtom* atoms);
	void spreadCharge(Lattice &lattice);
	void gatherForce(Lattice &lattice, CudaForce* force);
	void gatherForceSetCallback(ComputePmeCUDADevice* devicePtr_in);
	void waitGatherForceDone();
};

//
// Cuda implementation of PmeTranspose class
//
class CudaPmeTranspose : public PmeTranspose {
private:
	int deviceID;
	cudaStream_t stream;
	float2* d_data;
#ifndef P2P_ENABLE_3D
	float2* d_buffer;
#endif
	// List of device data pointers for transpose destinations on:
	// (a) this device on a different pencil (e.g. in XYZ->YZX transpose, on Y -pencil)
	// (b) different device on a different pencil
	// If NULL, use the local d_data -buffer
	std::vector<float2*> dataPtrsYZX;
	std::vector<float2*> dataPtrsZXY;

	// Batch data
	int max_nx_YZX[3];
	TransposeBatch<float2> *batchesYZX;
	int max_nx_ZXY[3];
	TransposeBatch<float2> *batchesZXY;

	void copyDataToPeerDevice(const int iblock,
		const int iblock_out, const int jblock_out, const int kblock_out,
		int deviceID_out, int permutation_out, float2* data_out);
public:
	CudaPmeTranspose(PmeGrid pmeGrid, const int permutation,
		const int jblock, const int kblock, int deviceID, cudaStream_t stream);
	~CudaPmeTranspose();
	void setDataPtrsYZX(std::vector<float2*>& dataPtrsNew, float2* data);
	void setDataPtrsZXY(std::vector<float2*>& dataPtrsNew, float2* data);
	void transposeXYZtoYZX(const float2* data);
	void transposeXYZtoZXY(const float2* data);
	// void waitTransposeDone();
	void waitStreamSynchronize();
	void copyDataDeviceToHost(const int iblock, float2* h_data, const int h_dataSize);
	void copyDataHostToDevice(const int iblock, float2* data_in, float2* data_out);
#ifndef P2P_ENABLE_3D
	void copyDataDeviceToDevice(const int iblock, float2* data_out);
	float2* getBuffer(const int iblock);
#endif
	void copyDataToPeerDeviceYZX(const int iblock, int deviceID_out, int permutation_out, float2* data_out);
	void copyDataToPeerDeviceZXY(const int iblock, int deviceID_out, int permutation_out, float2* data_out);
};
#endif // NAMD_CUDA
#endif // CUDAPMESOLVERUTIL_H