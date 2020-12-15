#include <cuda.h>
#include <cub/cub.cuh>
#include "CudaUtils.h"
#include "CudaTileListKernel.h"
#include "CudaComputeGBISKernel.h"
#define GBIS_CUDA
#include "ComputeGBIS.inl"
#include "DeviceCUDA.h"
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;

#define LARGE_FLOAT (float)(1.0e10)

#define GBISKERNEL_NUM_WARP 4

__device__ __forceinline__
void shuffleNext(float& w) {
  w = __shfl(w, (threadIdx.x+1) & (WARPSIZE-1));
}

// Generic
template <int phase> struct GBISParam {};
template <int phase> struct GBISInput {};
template <int phase> struct GBISResults {};

// ------------------------------------------------------------------------------------------------
// Phase 1 specialization

template <> struct GBISParam<1> {
  float a_cut;
};

template <> struct GBISInput<1> {
  float qi, qj;
  float intRad0j, intRadSi;
  // inp1 = intRad0
  // inp2 = intRadS
  __device__ __forceinline__ void loadI(const int i, const float* inp1, const float* inp2, const float* inp3) {
    qi = inp1[i];
    intRadSi = inp2[i];
  }
  __device__ __forceinline__ void loadJ(const int i, const float* inp1, const float* inp2, const float* inp3) {
    qj = inp2[i];
    intRad0j = inp1[i];
  }
  __device__ __forceinline__ void initQi(const GBISParam<1> param, const float q) {}
  __device__ __forceinline__ void initQj(const float q) {}
  __device__ __forceinline__ void shuffleNext() {
    qj       = __shfl(qj,       (threadIdx.x+1) & (WARPSIZE-1));
    intRad0j = __shfl(intRad0j, (threadIdx.x+1) & (WARPSIZE-1));    
  }
};

template <> struct GBISResults<1> {
  float psiSum;
  __device__ __forceinline__ void init() {psiSum = 0.0f;}
  __device__ __forceinline__ void shuffleNext() {
    psiSum = __shfl(psiSum, (threadIdx.x+1) & (WARPSIZE-1));
  }
};

// calculate H(i,j) [and not H(j,i)]
template <bool doEnergy, bool doSlow>
__device__ __forceinline__ void calcGBISPhase(const float r2, const float dx, const float dy, const float dz,
  const float cutoff2, const GBISParam<1> param, const GBISInput<1> inp,
  GBISResults<1>& resI, GBISResults<1>& resJ, float& energy) {

  if (r2 < cutoff2 && r2 > 0.01f) {
    float r_i = rsqrtf(r2);
    float r  = r2 * r_i;
    float hij;
    int dij;
    CalcH(r, r2, r_i, param.a_cut, inp.qi, inp.qj, hij, dij);
    resI.psiSum += hij;
    float hji;
    int dji;
    CalcH(r, r2, r_i, param.a_cut, inp.intRad0j, inp.intRadSi, hji, dji);
    resJ.psiSum += hji;
  }

}

__device__ __forceinline__ void writeResult(const int i, const GBISResults<1> res, float* out, float4* forces) {
  atomicAdd(&out[i], res.psiSum);
}

// ------------------------------------------------------------------------------------------------
// Phase 2

template <> struct GBISParam<2> {
  float kappa;
  float epsilon_p_i, epsilon_s_i;
  float smoothDist;
  float r_cut2, r_cut_2, r_cut_4;
  float scaling;
};

template <> struct GBISInput<2> {
  float qi, qj;
  float bornRadI, bornRadJ;
  // inp1 = bornRad
  __device__ __forceinline__ void loadI(const int i, const float* inp1, const float* inp2, const float* inp3) {
    bornRadI = inp1[i];
  }
  __device__ __forceinline__ void loadJ(const int i, const float* inp1, const float* inp2, const float* inp3) {
    bornRadJ = inp1[i];
  }
  __device__ __forceinline__ void initQi(const GBISParam<2> param, const float q) {
    qi = -q*param.scaling;    
  }
  __device__ __forceinline__ void initQj(const float q) {
    qj = q;
  }
  __device__ __forceinline__ void shuffleNext() {
    qj       = __shfl(qj,       (threadIdx.x+1) & (WARPSIZE-1));
    bornRadJ = __shfl(bornRadJ, (threadIdx.x+1) & (WARPSIZE-1));    
  }
};

template <> struct GBISResults<2> {
  float3 force;
  float dEdaSum;
  __device__ __forceinline__ void init() {force.x = 0.0f; force.y = 0.0f; force.z = 0.0f; dEdaSum = 0.0f;}
  __device__ __forceinline__ void shuffleNext() {
    force.x = __shfl(force.x, (threadIdx.x+1) & (WARPSIZE-1));
    force.y = __shfl(force.y, (threadIdx.x+1) & (WARPSIZE-1));
    force.z = __shfl(force.z, (threadIdx.x+1) & (WARPSIZE-1));
    dEdaSum = __shfl(dEdaSum, (threadIdx.x+1) & (WARPSIZE-1));
  }
};

template <bool doEnergy, bool doSlow>
__device__ __forceinline__ void calcGBISPhase(const float r2, const float dx, const float dy, const float dz,
  const float cutoff2, const GBISParam<2> param, const GBISInput<2> inp,
  GBISResults<2>& resI, GBISResults<2>& resJ, float& energyT) {

  if (r2 < cutoff2 && r2 > 0.01f) {
    float r_i = rsqrtf(r2);
    float r  = r2 * r_i;
    //float bornRadJ = sh_jBornRad[j];

    //calculate GB energy
    float qiqj = inp.qi*inp.qj;
    float aiaj = inp.bornRadI*inp.bornRadJ;
    float aiaj4 = 4*aiaj;
    float expr2aiaj4 = exp(-r2/aiaj4);
    float fij = sqrt(r2 + aiaj*expr2aiaj4);
    float f_i = 1/fij;
    float expkappa = exp(-param.kappa*fij);
    float Dij = param.epsilon_p_i - expkappa*param.epsilon_s_i;
    float gbEij = qiqj*Dij*f_i;

    //calculate energy derivatives
    float ddrfij = r*f_i*(1.f - 0.25f*expr2aiaj4);
    float ddrf_i = -ddrfij*f_i*f_i;
    float ddrDij = param.kappa*expkappa*ddrfij*param.epsilon_s_i;
    float ddrGbEij = qiqj*(ddrDij*f_i+Dij*ddrf_i);

    //NAMD smoothing function
    float scale = 1.f;
    float ddrScale = 0.f;
    float forcedEdr;
    if (param.smoothDist > 0.f) {
      scale = r2 * param.r_cut_2 - 1.f;
      scale *= scale;
      ddrScale = r*(r2-param.r_cut2)*param.r_cut_4;
      if (doEnergy) energyT += gbEij * scale;
      forcedEdr = (ddrGbEij)*scale + (gbEij)*ddrScale;
    } else {
      if (doEnergy) energyT += gbEij;
      forcedEdr = ddrGbEij;
    }

    //add dEda
    if (doSlow) {
      float dEdai = 0.5f*qiqj*f_i*f_i
                *(param.kappa*param.epsilon_s_i*expkappa-Dij*f_i)
                *(aiaj+0.25f*r2)*expr2aiaj4/inp.bornRadI*scale;//0
      resI.dEdaSum += dEdai;
      float dEdaj = 0.5f*qiqj*f_i*f_i
                *(param.kappa*param.epsilon_s_i*expkappa-Dij*f_i)
                *(aiaj+0.25f*r2)*expr2aiaj4/inp.bornRadJ*scale;//0
      resJ.dEdaSum += dEdaj;
    }

    forcedEdr *= r_i;
    float tmpx = dx*forcedEdr;
    float tmpy = dy*forcedEdr;
    float tmpz = dz*forcedEdr;
    resI.force.x += tmpx;
    resI.force.y += tmpy;
    resI.force.z += tmpz;
    resJ.force.x -= tmpx;
    resJ.force.y -= tmpy;
    resJ.force.z -= tmpz;
  }

  // GB Self Energy
  if (doEnergy) {
    if (r2 < 0.01f) {
      float fij = inp.bornRadI;//inf
      float expkappa = exp(-param.kappa*fij);//0
      float Dij = param.epsilon_p_i - expkappa*param.epsilon_s_i;
      float gbEij = inp.qi*(inp.qi / (-param.scaling) )*Dij/fij;
      energyT += 0.5f*gbEij;
    } //same atom or within cutoff
  }
}

__device__ __forceinline__ void writeResult(const int i, const GBISResults<2> res, float* out, float4* forces) {
  atomicAdd(&out[i], res.dEdaSum);
  atomicAdd(&forces[i].x, res.force.x);
  atomicAdd(&forces[i].y, res.force.y);
  atomicAdd(&forces[i].z, res.force.z);
}

// ------------------------------------------------------------------------------------------------
// Phase 3
template <> struct GBISParam<3> {
  float a_cut;
};

template <> struct GBISInput<3> {
  float qi;
  float intRadSJ, intRadJ0, intRadIS;
  float dHdrPrefixI, dHdrPrefixJ;
  // inp1 = intRad0
  // inp2 = intRadS
  // inp3 = dHdrPrefix
  __device__ __forceinline__ void loadI(const int i, const float* inp1, const float* inp2, const float* inp3) {
    qi          = inp1[i];
    intRadIS    = inp2[i];
    dHdrPrefixI = inp3[i];
  }
  __device__ __forceinline__ void loadJ(const int i, const float* inp1, const float* inp2, const float* inp3) {
    intRadJ0    = inp1[i];
    intRadSJ    = inp2[i];
    dHdrPrefixJ = inp3[i];
  }
  __device__ __forceinline__ void initQi(const GBISParam<3> param, const float q) {}
  __device__ __forceinline__ void initQj(const float q) {}
  __device__ __forceinline__ void shuffleNext() {
    intRadSJ    = __shfl(intRadSJ,    (threadIdx.x+1) & (WARPSIZE-1));
    dHdrPrefixJ = __shfl(dHdrPrefixJ, (threadIdx.x+1) & (WARPSIZE-1));
    intRadJ0    = __shfl(intRadJ0,    (threadIdx.x+1) & (WARPSIZE-1));
  }
};

template <> struct GBISResults<3> {
  float3 force;
  __device__ __forceinline__ void init() {force.x = 0.0f; force.y = 0.0f; force.z = 0.0f;}
  __device__ __forceinline__ void shuffleNext() {
    force.x = __shfl(force.x, (threadIdx.x+1) & (WARPSIZE-1));
    force.y = __shfl(force.y, (threadIdx.x+1) & (WARPSIZE-1));
    force.z = __shfl(force.z, (threadIdx.x+1) & (WARPSIZE-1));
  }
};

template <bool doEnergy, bool doSlow>
__device__ __forceinline__ void calcGBISPhase(const float r2, const float dx, const float dy, const float dz,
  const float cutoff2, const GBISParam<3> param, const GBISInput<3> inp,
  GBISResults<3>& resI, GBISResults<3>& resJ, float& energy) {

  if (r2 < cutoff2 && r2 > 0.01f) {
    float r_i = rsqrtf(r2);
    float r  = r2 * r_i;
    float dhij, dhji;
    int dij, dji;
    CalcDH(r, r2, r_i, param.a_cut, inp.qi,       inp.intRadSJ, dhij, dij);
    CalcDH(r, r2, r_i, param.a_cut, inp.intRadJ0, inp.intRadIS, dhji, dji);

    float forceAlpha = r_i*(inp.dHdrPrefixI*dhij + inp.dHdrPrefixJ*dhji);
    float tmpx = dx * forceAlpha;
    float tmpy = dy * forceAlpha;
    float tmpz = dz * forceAlpha;
    resI.force.x += tmpx;
    resI.force.y += tmpy;
    resI.force.z += tmpz;
    resJ.force.x -= tmpx;
    resJ.force.y -= tmpy;
    resJ.force.z -= tmpz;
  }

}

__device__ __forceinline__ void writeResult(const int i, const GBISResults<3> res, float* out, float4* forces) {
  atomicAdd(&forces[i].x, res.force.x);
  atomicAdd(&forces[i].y, res.force.y);
  atomicAdd(&forces[i].z, res.force.z);
}

// ------------------------------------------------------------------------------------------------

//
// GBIS kernel
//
template <bool doEnergy, bool doSlow, int phase>
__global__ void
__launch_bounds__(32*GBISKERNEL_NUM_WARP, 6)
GBIS_Kernel(const int numTileLists,
  const TileList* __restrict__ tileLists,
  const int* __restrict__ tileJatomStart,
  const PatchPairRecord* __restrict__ patchPairs,
  const float latticeX, const float latticeY, const float latticeZ,
  const float4* __restrict__ xyzq,
  const float cutoff2,
  const GBISParam<phase> param,
  const float* __restrict__ inp1,
  const float* __restrict__ inp2,
  const float* __restrict__ inp3,
  float* __restrict__ out,
  float4* __restrict__ forces,
  TileListVirialEnergy* __restrict__ virialEnergy) {

  // Single warp takes care of one list of tiles
  for (int itileList = (threadIdx.x + blockDim.x*blockIdx.x)/WARPSIZE;itileList < numTileLists;itileList += blockDim.x*gridDim.x/WARPSIZE)
  {
    TileList tmp = tileLists[itileList];
    int iatomStart = tmp.iatomStart;
    int jtileStart = tmp.jtileStart;
    int jtileEnd   = tmp.jtileEnd;
    PatchPairRecord PPStmp = patchPairs[itileList];
    int iatomSize     = PPStmp.iatomSize;
    int jatomSize     = PPStmp.jatomSize;

    float shx = tmp.offsetXYZ.x*latticeX;
    float shy = tmp.offsetXYZ.y*latticeY;
    float shz = tmp.offsetXYZ.z*latticeZ;

    // Warp index (0...warpsize-1)
    const int wid = threadIdx.x % WARPSIZE;

    // Load i-atom data (and shift coordinates)
    float4 xyzq_i = xyzq[iatomStart + wid];
    xyzq_i.x += shx;
    xyzq_i.y += shy;
    xyzq_i.z += shz;

    GBISInput<phase> inp;
    inp.loadI(iatomStart + wid, inp1, inp2, inp3);
    if (phase == 2) inp.initQi(param, xyzq_i.w);

    // Number of i loops
    int nloopi = min(iatomSize - iatomStart, WARPSIZE);

    GBISResults<phase> resI;
    resI.init();
    float energyT;
    if (doEnergy) energyT = 0.0f;

    for (int jtile=jtileStart;jtile <= jtileEnd;jtile++) {

      // Load j-atom starting index and exclusion mask
      int jatomStart = tileJatomStart[jtile];

      // Load coordinates and charge
      float4 xyzq_j  = xyzq[jatomStart + wid];

      inp.loadJ(jatomStart + wid, inp1, inp2, inp3);
      if (phase == 2) inp.initQj(xyzq_j.w);

      // Number of j loops
      int nloopj = min(jatomSize - jatomStart, WARPSIZE);

      const bool diag_tile = (iatomStart == jatomStart);
      const int modval = (diag_tile) ? 2*WARPSIZE-1 : WARPSIZE-1;
      int t = (phase != 2 && diag_tile) ? 1 : 0;
      if (phase != 2 && diag_tile) {
        inp.shuffleNext();
      }

      GBISResults<phase> resJ;
      resJ.init();

      for (;t < WARPSIZE;t++) {
        int j = (t + wid) & modval;

        float dx = __shfl(xyzq_j.x,j) - xyzq_i.x;
        float dy = __shfl(xyzq_j.y,j) - xyzq_i.y;
        float dz = __shfl(xyzq_j.z,j) - xyzq_i.z;

        float r2 = dx*dx + dy*dy + dz*dz;

        if (j < nloopj && wid < nloopi) {
          calcGBISPhase<doEnergy, doSlow>(r2, dx, dy, dz, cutoff2, param, inp, resI, resJ, energyT);
        }

        inp.shuffleNext();
        resJ.shuffleNext();
      } // t

      // Write j
      writeResult(jatomStart + wid, resJ, out, forces);

    } // jtile

    // Write i
    writeResult(iatomStart + wid, resI, out, forces);

    // Reduce energy
    if (doEnergy) {
      typedef cub::WarpReduce<double> WarpReduceDouble;
      __shared__ typename WarpReduceDouble::TempStorage tempStorage[GBISKERNEL_NUM_WARP];
      int warpId = threadIdx.x / WARPSIZE;
      volatile double energyTot = WarpReduceDouble(tempStorage[warpId]).Sum((double)energyT);
      if (wid == 0) virialEnergy[itileList].energyGBIS = energyTot;
    }

  }
}

#undef GBIS_CUDA

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

//
// Class constructor
//
CudaComputeGBISKernel::CudaComputeGBISKernel(int deviceID) : deviceID(deviceID) {
  cudaCheck(cudaSetDevice(deviceID));

  intRad0 = NULL;
  intRad0Size = 0;

  intRadS = NULL;
  intRadSSize = 0;

  psiSum = NULL;
  psiSumSize = 0;

  bornRad = NULL;
  bornRadSize = 0;

  dEdaSum = NULL;
  dEdaSumSize = 0;

  dHdrPrefix = NULL;
  dHdrPrefixSize = 0;

}

//
// Class destructor
//
CudaComputeGBISKernel::~CudaComputeGBISKernel() {
  cudaCheck(cudaSetDevice(deviceID));
  if (intRad0 != NULL) deallocate_device<float>(&intRad0);
  if (intRadS != NULL) deallocate_device<float>(&intRadS);
  if (psiSum != NULL) deallocate_device<float>(&psiSum);
  if (bornRad != NULL) deallocate_device<float>(&bornRad);
  if (dEdaSum != NULL) deallocate_device<float>(&dEdaSum);
  if (dHdrPrefix != NULL) deallocate_device<float>(&dHdrPrefix);
}

//
// Update (copy Host->Device) intRad0, intRadS
//
void CudaComputeGBISKernel::updateIntRad(const int atomStorageSize, float* intRad0H, float* intRadSH,
  cudaStream_t stream) {

  reallocate_device<float>(&intRad0, &intRad0Size, atomStorageSize, 1.2f);
  reallocate_device<float>(&intRadS, &intRadSSize, atomStorageSize, 1.2f);

  copy_HtoD<float>(intRad0H, intRad0, atomStorageSize, stream);
  copy_HtoD<float>(intRadSH, intRadS, atomStorageSize, stream);
}

//
// Update (copy Host->Device) bornRad
//
void CudaComputeGBISKernel::updateBornRad(const int atomStorageSize, float* bornRadH, cudaStream_t stream) {
  reallocate_device<float>(&bornRad, &bornRadSize, atomStorageSize, 1.2f);
  copy_HtoD<float>(bornRadH, bornRad, atomStorageSize, stream);
}

//
// Update (copy Host->Device) dHdrPrefix
//
void CudaComputeGBISKernel::update_dHdrPrefix(const int atomStorageSize, float* dHdrPrefixH, cudaStream_t stream) {
  reallocate_device<float>(&dHdrPrefix, &dHdrPrefixSize, atomStorageSize, 1.2f);
  copy_HtoD<float>(dHdrPrefixH, dHdrPrefix, atomStorageSize, stream);
}

//
// Phase 1
//
void CudaComputeGBISKernel::GBISphase1(CudaTileListKernel& tlKernel, const int atomStorageSize,
  const float latticeX, const float latticeY, const float latticeZ, const float a_cut, float* h_psiSum,
  cudaStream_t stream) {

  reallocate_device<float>(&psiSum, &psiSumSize, atomStorageSize, 1.2f);
  clear_device_array<float>(psiSum, atomStorageSize, stream);

  int nwarp = GBISKERNEL_NUM_WARP;
  int nthread = WARPSIZE*nwarp;
  int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getNumTileListsGBIS()-1)/nwarp+1);

  GBISParam<1> param;
  param.a_cut = a_cut;

  float cutoff2 = (a_cut + FS_MAX)*(a_cut + FS_MAX);

  GBIS_Kernel<false, false, 1> <<< nblock, nthread, 0, stream >>>
  (tlKernel.getNumTileListsGBIS(), tlKernel.getTileListsGBIS(), tlKernel.getTileJatomStartGBIS(),
    tlKernel.getPatchPairs(), latticeX, latticeY, latticeZ, tlKernel.get_xyzq(), cutoff2,
    param, intRad0, intRadS, NULL, psiSum, NULL, NULL);

  cudaCheck(cudaGetLastError());

  copy_DtoH<float>(psiSum, h_psiSum, atomStorageSize, stream);
}

//
// Phase 2
//
void CudaComputeGBISKernel::GBISphase2(CudaTileListKernel& tlKernel, const int atomStorageSize,
  const bool doEnergy, const bool doSlow,
  const float latticeX, const float latticeY, const float latticeZ,
  const float r_cut, const float scaling, const float kappa, const float smoothDist,
  const float epsilon_p, const float epsilon_s,
  float4* d_forces, float* h_dEdaSum, cudaStream_t stream) {

  reallocate_device<float>(&dEdaSum, &dEdaSumSize, atomStorageSize, 1.2f);
  clear_device_array<float>(dEdaSum, atomStorageSize, stream);

  if (doEnergy) {
    tlKernel.setTileListVirialEnergyGBISLength(tlKernel.getNumTileListsGBIS());
  }

  int nwarp = GBISKERNEL_NUM_WARP;
  int nthread = WARPSIZE*nwarp;

  int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getNumTileListsGBIS()-1)/nwarp+1);

  GBISParam<2> param;
  param.r_cut2      = r_cut*r_cut;
  param.r_cut_2     = 1.f / param.r_cut2;
  param.r_cut_4     = 4.f*param.r_cut_2*param.r_cut_2;
  param.epsilon_s_i = 1.f / epsilon_s;
  param.epsilon_p_i = 1.f / epsilon_p;
  param.scaling     = scaling;
  param.kappa       = kappa;
  param.smoothDist  = smoothDist;

#define CALL(DOENERGY, DOSLOW) GBIS_Kernel<DOENERGY, DOSLOW, 2> \
     <<< nblock, nthread, 0, stream >>> \
    (tlKernel.getNumTileListsGBIS(), tlKernel.getTileListsGBIS(), tlKernel.getTileJatomStartGBIS(), \
      tlKernel.getPatchPairs(), latticeX, latticeY, latticeZ, tlKernel.get_xyzq(), param.r_cut2, \
      param, bornRad, NULL, NULL, dEdaSum, d_forces, tlKernel.getTileListVirialEnergy())

  if (!doEnergy && !doSlow) CALL(false, false);
  if (!doEnergy &&  doSlow) CALL(false, true);
  if ( doEnergy && !doSlow) CALL(true,  false);
  if ( doEnergy &&  doSlow) CALL(true,  true);

  cudaCheck(cudaGetLastError());

  copy_DtoH<float>(dEdaSum, h_dEdaSum, atomStorageSize, stream);
}

//
// Phase 3
//
void CudaComputeGBISKernel::GBISphase3(CudaTileListKernel& tlKernel, const int atomStorageSize,
  const float latticeX, const float latticeY, const float latticeZ, const float a_cut,
  float4* d_forces, cudaStream_t stream) {

  int nwarp = GBISKERNEL_NUM_WARP;
  int nthread = WARPSIZE*nwarp;
  int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getNumTileListsGBIS()-1)/nwarp+1);

  GBISParam<3> param;
  param.a_cut = a_cut;

  float cutoff2 = (a_cut + FS_MAX)*(a_cut + FS_MAX);

  GBIS_Kernel<false, false, 3> <<< nblock, nthread, 0, stream >>>
  (tlKernel.getNumTileListsGBIS(), tlKernel.getTileListsGBIS(), tlKernel.getTileJatomStartGBIS(),
    tlKernel.getPatchPairs(), latticeX, latticeY, latticeZ, tlKernel.get_xyzq(), cutoff2,
    param, intRad0, intRadS, dHdrPrefix, NULL, d_forces, NULL);

  cudaCheck(cudaGetLastError());
}
