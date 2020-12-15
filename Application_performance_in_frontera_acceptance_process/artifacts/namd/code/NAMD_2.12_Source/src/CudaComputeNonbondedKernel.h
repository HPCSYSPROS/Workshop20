#ifndef CUDACOMPUTENONBONDEDKERNEL_H
#define CUDACOMPUTENONBONDEDKERNEL_H
#include "CudaUtils.h"
#include "CudaTileListKernel.h"
#ifdef NAMD_CUDA

class CudaComputeNonbondedKernel {
private:

  const int deviceID;
  const bool doStreaming;

  // Force and energy tables
  cudaArray_t forceArray;
  cudaArray_t energyArray;

#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  bool forceTableTexActive;
  cudaTextureObject_t forceTableTex;

  bool energyTableTexActive;
  cudaTextureObject_t energyTableTex;
#endif

  bool vdwCoefTableTextureBound;
#ifdef DISABLE_CUDA_TEXTURE_OBJECTS
  bool force_table_bound;
  bool energy_table_bound;
#endif

  // Exclusions
  unsigned int* overflowExclusions;
  int overflowExclusionsSize;

  int2* exclIndexMaxDiff;
  int exclIndexMaxDiffSize;

  // Atom indices
  int* atomIndex;
  int atomIndexSize;

  // VdW types
  int* vdwTypes;
  int vdwTypesSize;

  // VdW Lennard-Jones coefficient table.
  float2* vdwCoefTable;
  int vdwCoefTableSize;
  // Width (and height) of the VdW coeffient table
  int vdwCoefTableWidth;

  unsigned int* patchNumCount;
  int patchNumCountSize;

  int* patchReadyQueue;
  int patchReadyQueueSize;

public:
  CudaComputeNonbondedKernel(int deviceID, bool doStreaming);
  ~CudaComputeNonbondedKernel();

  void updateVdwTypesExcl(const int atomStorageSize, const int* h_vdwTypes,
    const int2* h_exclIndexMaxDiff, const int* h_atomIndex, cudaStream_t stream);

  void nonbondedForce(CudaTileListKernel& tlKernel,
    const int atomStorageSize, const bool doPairlist,
    const bool doEnergy, const bool doVirial, const bool doSlow,
    const float latticeX, const float latticeY, const float latticeZ,
    const float4* h_xyzq, const float cutoff2, 
    float4* d_forces, float4* d_forcesSlow,
    float4* h_forces, float4* h_forcesSlow,
    cudaStream_t stream);

  void reduceVirialEnergy(CudaTileListKernel& tlKernel,
    const int atomStorageSize, const bool doEnergy, const bool doVirial, const bool doSlow, const bool doGBIS,
    float4* d_forces, float4* d_forcesSlow,
    VirialEnergy* d_virialEnergy, cudaStream_t stream);

  void getVirialEnergy(VirialEnergy* h_virialEnergy, cudaStream_t stream);

  void bindExclusions(int numExclusions, unsigned int* exclusion_bits);
  void bindVdwCoefTable(float2* h_vdwCoefTable, int vdwCoefTableWidthIn);

#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  void bindTextureObject(int tableSize, float4* h_table, cudaArray_t& array, cudaTextureObject_t& tableTex);
#endif

  void bindForceAndEnergyTable(int tableSize,
  float4* h_forceTable,
  float4* h_energyTable);

  int* getPatchReadyQueue();
};

#endif // NAMD_CUDA
#endif // CUDACOMPUTENONBONDEDKERNEL_H