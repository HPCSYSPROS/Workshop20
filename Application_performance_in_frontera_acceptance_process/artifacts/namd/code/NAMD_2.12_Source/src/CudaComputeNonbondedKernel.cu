#include <cuda.h>
#include <cub/cub.cuh>
#include "CudaComputeNonbondedKernel.h"
#include "CudaTileListKernel.h"
#include "DeviceCUDA.h"
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;

#define OVERALLOC 1.2f

void NAMD_die(const char *);

texture<float2, 1, cudaReadModeElementType> vdwCoefTableTexture;

#ifdef DISABLE_CUDA_TEXTURE_OBJECTS
texture<float4, 1, cudaReadModeElementType> force_table;
texture<float4, 1, cudaReadModeElementType> energy_table;
#endif

#define MAX_CONST_EXCLUSIONS 2048  // cache size is 8k
__constant__ unsigned int constExclusions[MAX_CONST_EXCLUSIONS];

#define NONBONDKERNEL_NUM_WARP 4

template<bool doEnergy, bool doSlow>
__device__ __forceinline__
void calcForceEnergy(const float r2, const float qi, const float qj,
  const float dx, const float dy, const float dz,
  const int vdwtypei, const int vdwtypej, const float2* __restrict__ vdwCoefTable,
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  cudaTextureObject_t forceTableTex, cudaTextureObject_t energyTableTex,
#endif
  float3& iforce, float3& iforceSlow, float3& jforce, float3& jforceSlow,
  float& energyVdw, float& energyElec, float& energySlow) {

  int vdwIndex = vdwtypej + vdwtypei;
#if __CUDA_ARCH__ >= 350
  float2 ljab = __ldg(&vdwCoefTable[vdwIndex]);
#else
  float2 ljab = tex1Dfetch(vdwCoefTableTexture, vdwIndex);
#endif

  float rinv = rsqrtf(r2);
  float4 ei;
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  float4 fi = tex1D<float4>(forceTableTex, rinv);
  if (doEnergy) ei = tex1D<float4>(energyTableTex, rinv);
#else
  float4 fi = tex1D(force_table, rinv);
  if (doEnergy) ei = tex1D(energy_table, rinv);
#endif

  float fSlow = qi * qj;
  float f = ljab.x * fi.z + ljab.y * fi.y + fSlow * fi.x;

  if (doEnergy) {
    energyVdw  += ljab.x * ei.z + ljab.y * ei.y;
    energyElec += fSlow * ei.x;
    if (doSlow) energySlow += fSlow * ei.w;
  }
  if (doSlow) fSlow *= fi.w;

  float fx = dx * f;
  float fy = dy * f;
  float fz = dz * f;
  iforce.x += fx;
  iforce.y += fy;
  iforce.z += fz;
  jforce.x -= fx;
  jforce.y -= fy;
  jforce.z -= fz;

  if (doSlow) {
    float fxSlow = dx * fSlow;
    float fySlow = dy * fSlow;
    float fzSlow = dz * fSlow;
    iforceSlow.x += fxSlow;
    iforceSlow.y += fySlow;
    iforceSlow.z += fzSlow;
    jforceSlow.x -= fxSlow;
    jforceSlow.y -= fySlow;
    jforceSlow.z -= fzSlow;
  }
}

template<bool doSlow>
__device__ __forceinline__
void storeForces(const int pos, const float3 force, const float3 forceSlow,
  float4* __restrict__ devForces, float4* __restrict__ devForcesSlow) {
  atomicAdd(&devForces[pos].x, force.x);
  atomicAdd(&devForces[pos].y, force.y);
  atomicAdd(&devForces[pos].z, force.z);
  if (doSlow) {
    atomicAdd(&devForcesSlow[pos].x, forceSlow.x);
    atomicAdd(&devForcesSlow[pos].y, forceSlow.y);
    atomicAdd(&devForcesSlow[pos].z, forceSlow.z);
  }
}

template<bool doSlow>
__device__ __forceinline__
void storeForces(const int pos, const float3 force, const float3 forceSlow,
  float3* __restrict__ forces, float3* __restrict__ forcesSlow) {
  atomicAdd(&forces[pos].x, force.x);
  atomicAdd(&forces[pos].y, force.y);
  atomicAdd(&forces[pos].z, force.z);
  if (doSlow) {
    atomicAdd(&forcesSlow[pos].x, forceSlow.x);
    atomicAdd(&forcesSlow[pos].y, forceSlow.y);
    atomicAdd(&forcesSlow[pos].z, forceSlow.z);
  }
}

template<bool doPairlist>
__device__ __forceinline__
void shuffleNext(float& xyzq_j_w, int& vdwtypej, int& jatomIndex, int& jexclMaxdiff, int& jexclIndex) {
  xyzq_j_w = __shfl(xyzq_j_w, (threadIdx.x+1) & (WARPSIZE-1));
  vdwtypej = __shfl(vdwtypej, (threadIdx.x+1) & (WARPSIZE-1));
  if (doPairlist) {
    jatomIndex   = __shfl(jatomIndex, (threadIdx.x+1) & (WARPSIZE-1));    
    jexclIndex   = __shfl(jexclIndex, (threadIdx.x+1) & (WARPSIZE-1) );
    jexclMaxdiff = __shfl(jexclMaxdiff, (threadIdx.x+1) & (WARPSIZE-1) );
  }
}

template<bool doPairlist>
__device__ __forceinline__
void shuffleNext(float& xyzq_j_w, int& vdwtypej, int& jatomIndex) {
  xyzq_j_w = __shfl(xyzq_j_w, (threadIdx.x+1) & (WARPSIZE-1));
  vdwtypej = __shfl(vdwtypej, (threadIdx.x+1) & (WARPSIZE-1));
  if (doPairlist) {
    jatomIndex   = __shfl(jatomIndex, (threadIdx.x+1) & (WARPSIZE-1));    
  }
}

template<bool doSlow>
__device__ __forceinline__
void shuffleNext(float3& jforce, float3& jforceSlow) {
  jforce.x = __shfl(jforce.x, (threadIdx.x+1)&(WARPSIZE-1));
  jforce.y = __shfl(jforce.y, (threadIdx.x+1)&(WARPSIZE-1));
  jforce.z = __shfl(jforce.z, (threadIdx.x+1)&(WARPSIZE-1));
  if (doSlow) {
    jforceSlow.x = __shfl(jforceSlow.x, (threadIdx.x+1)&(WARPSIZE-1));
    jforceSlow.y = __shfl(jforceSlow.y, (threadIdx.x+1)&(WARPSIZE-1));
    jforceSlow.z = __shfl(jforceSlow.z, (threadIdx.x+1)&(WARPSIZE-1));
  }
}

//#define USE_NEW_EXCL_METHOD

//
// Returns the lower estimate for the distance between a bounding box and a set of atoms
//
__device__ __forceinline__ float distsq(const BoundingBox a, const float4 b) {
  float dx = max(0.0f, fabsf(a.x - b.x) - a.wx);
  float dy = max(0.0f, fabsf(a.y - b.y) - a.wy);
  float dz = max(0.0f, fabsf(a.z - b.z) - a.wz);
  float r2 = dx*dx + dy*dy + dz*dz;
  return r2;
}

#define LARGE_FLOAT (float)(1.0e10)

//
// Nonbonded force kernel
//
template <bool doEnergy, bool doVirial, bool doSlow, bool doPairlist, bool doStreaming>
__global__ void
__launch_bounds__(WARPSIZE*NONBONDKERNEL_NUM_WARP,
  doPairlist ? (10) : (doEnergy ? (10) : (10) )
  )
nonbondedForceKernel(const int start, const int numTileLists,
  const TileList* __restrict__ tileLists, TileExcl* __restrict__ tileExcls,
  const int* __restrict__ tileJatomStart,
  const int vdwCoefTableWidth, const float2* __restrict__ vdwCoefTable, const int* __restrict__ vdwTypes,
  const float latticeX, const float latticeY, const float latticeZ,
  const float4* __restrict__ xyzq, const float cutoff2,
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  cudaTextureObject_t forceTableTex, cudaTextureObject_t energyTableTex,
#endif
  // ----------
  // doPairlist
  const int atomStorageSize, const float plcutoff2, const PatchPairRecord* __restrict__ patchPairs,
  const int* __restrict__ atomIndex,
  const int2* __restrict__ exclIndexMaxDiff, const unsigned int* __restrict__ overflowExclusions,
  unsigned int* __restrict__ tileListDepth, int* __restrict__ tileListOrder,
  int* __restrict__ jtiles, TileListStat* __restrict__ tileListStat,
  const BoundingBox* __restrict__ boundingBoxes,
#ifdef USE_NEW_EXCL_METHOD
  const int* __restrict__ minmaxExclAtom,
#endif
  // ----------
  float4* __restrict__ devForces, float4* __restrict__ devForcesSlow,
  // ---- USE_STREAMING_FORCES ----
  const int numPatches,
  unsigned int* __restrict__ patchNumCount,
  const CudaPatchRecord* __restrict__ cudaPatches,
  float4* __restrict__ mapForces, float4* __restrict__ mapForcesSlow,
  int* __restrict__ mapPatchReadyQueue,
  int* __restrict__ outputOrder,
  // ------------------------------
  TileListVirialEnergy* __restrict__ virialEnergy) {

  // Single warp takes care of one list of tiles
  // for (int itileList = (threadIdx.x + blockDim.x*blockIdx.x)/WARPSIZE;itileList < numTileLists;itileList += blockDim.x*gridDim.x/WARPSIZE)
  int itileList = start + (threadIdx.x + blockDim.x*blockIdx.x)/WARPSIZE;
  if (itileList < numTileLists)
  {

    float3 iforce;
    float3 iforceSlow;
    float energyVdw, energyElec, energySlow;
    int nexcluded;
    unsigned int itileListLen;
    int2 patchInd;
    int2 patchNumList;

    // Start computation
    {
      // Warp index (0...warpsize-1)
      const int wid = threadIdx.x % WARPSIZE;

      TileList tmp = tileLists[itileList];
      int iatomStart = tmp.iatomStart;
      int jtileStart = tmp.jtileStart;
      int jtileEnd   = tmp.jtileEnd;
      patchInd     = tmp.patchInd;
      patchNumList = tmp.patchNumList;

      float shx = tmp.offsetXYZ.x*latticeX;
      float shy = tmp.offsetXYZ.y*latticeY;
      float shz = tmp.offsetXYZ.z*latticeZ;

      int iatomSize, iatomFreeSize, jatomSize, jatomFreeSize;
      if (doPairlist) {
        PatchPairRecord PPStmp = patchPairs[itileList];
        iatomSize     = PPStmp.iatomSize;
        iatomFreeSize = PPStmp.iatomFreeSize;
        jatomSize     = PPStmp.jatomSize;
        jatomFreeSize = PPStmp.jatomFreeSize;
      }

      // Write to global memory here to avoid register spilling
      if (doVirial) {
        if (wid == 0) {
          virialEnergy[itileList].shx = shx;
          virialEnergy[itileList].shy = shy;
          virialEnergy[itileList].shz = shz;
        }
      }

      // Load i-atom data (and shift coordinates)
      float4 xyzq_i = xyzq[iatomStart + wid];
      xyzq_i.x += shx;
      xyzq_i.y += shy;
      xyzq_i.z += shz;
      int vdwtypei = vdwTypes[iatomStart + wid]*vdwCoefTableWidth;

      // Load i-atom data (and shift coordinates)
      BoundingBox boundingBoxI;
      if (doPairlist) {
        boundingBoxI = boundingBoxes[iatomStart/WARPSIZE];
        boundingBoxI.x += shx;
        boundingBoxI.y += shy;
        boundingBoxI.z += shz;
      }

      // Get i-atom global index
#ifdef USE_NEW_EXCL_METHOD
      int iatomIndex, minExclAtom, maxExclAtom;
#else
      int iatomIndex;
#endif
      if (doPairlist) {
#ifdef USE_NEW_EXCL_METHOD
        iatomIndex = atomIndex[iatomStart + wid];
        int2 tmp = minmaxExclAtom[iatomStart + wid];
        minExclAtom = tmp.x;
        maxExclAtom = tmp.y;
#else
        iatomIndex = atomIndex[iatomStart + wid];
#endif
      }

      // i-forces in registers
      // float3 iforce;
      iforce.x = 0.0f;
      iforce.y = 0.0f;
      iforce.z = 0.0f;

      // float3 iforceSlow;
      if (doSlow) {
        iforceSlow.x = 0.0f;
        iforceSlow.y = 0.0f;
        iforceSlow.z = 0.0f;
      }

      // float energyVdw, energyElec, energySlow;
      if (doEnergy) {
        energyVdw = 0.0f;
        energyElec = 0.0f;
        if (doSlow) energySlow = 0.0f;
      }

      // Number of exclusions
      // NOTE: Lowest bit is used as indicator bit for tile pairs:
      //       bit 0 tile has no atoms within pairlist cutoff
      //       bit 1 tile has atoms within pairlist cutoff
      // int nexcluded;
      if (doPairlist) nexcluded = 0;

      // Number of i loops and free atoms
      int nfreei;
      if (doPairlist) {
        int nloopi = min(iatomSize - iatomStart, WARPSIZE);
        nfreei = max(iatomFreeSize - iatomStart, 0);
        if (wid >= nloopi) {
          xyzq_i.x = -LARGE_FLOAT;
          xyzq_i.y = -LARGE_FLOAT;
          xyzq_i.z = -LARGE_FLOAT;
        }
      }

      // tile list stuff
      // int itileListLen;
      // int minJatomStart;
      if (doPairlist) {
        // minJatomStart = tileJatomStart[jtileStart];
        itileListLen = 0;
      }

      // Exclusion index and maxdiff
      int iexclIndex, iexclMaxdiff;
      if (doPairlist) {
        int2 tmp = exclIndexMaxDiff[iatomStart + wid];
        iexclIndex   = tmp.x;
        iexclMaxdiff = tmp.y;
      }

      for (int jtile=jtileStart;jtile <= jtileEnd;jtile++) {

        // Load j-atom starting index and exclusion mask
        int jatomStart = tileJatomStart[jtile];

        float4 xyzq_j = xyzq[jatomStart + wid];

        // Check for early bail
        if (doPairlist) {
          float r2bb = distsq(boundingBoxI, xyzq_j);
          if (__all(r2bb > plcutoff2)) continue;
        }
        unsigned int excl = (doPairlist) ? 0 : tileExcls[jtile].excl[wid];
        int vdwtypej = vdwTypes[jatomStart + wid];

        // Get i-atom global index
        int jatomIndex;
        if (doPairlist) {
          jatomIndex = atomIndex[jatomStart + wid];
        }

        // Number of j loops and free atoms
        int nfreej;
        if (doPairlist) {
          int nloopj = min(jatomSize - jatomStart, WARPSIZE);
          nfreej = max(jatomFreeSize - jatomStart, 0);
          //if (nfreei == 0 && nfreej == 0) continue;
          if (wid >= nloopj) {
            xyzq_j.x = LARGE_FLOAT;
            xyzq_j.y = LARGE_FLOAT;
            xyzq_j.z = LARGE_FLOAT;
          }
        }

        const bool self = (iatomStart == jatomStart);
        const int modval = (self) ? 2*WARPSIZE-1 : WARPSIZE-1;

        float3 jforce;
        jforce.x = 0.0f;
        jforce.y = 0.0f;
        jforce.z = 0.0f;
        
        float3 jforceSlow;
        if (doSlow) {
          jforceSlow.x = 0.0f;
          jforceSlow.y = 0.0f;
          jforceSlow.z = 0.0f;
        }

        int t = (self) ? 1 : 0;

        if (doPairlist) {
          // Build pair list
          // NOTE: Pairlist update, we must also include the diagonal since this is used
          //       in GBIS phase 2.
          // Clear the lowest (indicator) bit
          nexcluded &= (~1);

          // For self tiles, do the diagonal term (t=0).
          // NOTE: No energies are computed here, since this self-diagonal term is only for GBIS phase 2
          if (self) {
            int j = (0 + wid) & modval;
            // NOTE: __shfl() operation can give non-sense here because j may be >= WARPSIZE.
            //       However, if (j < WARPSIZE ..) below makes sure that these non-sense
            //       results are not actually every used
            float dx = __shfl(xyzq_j.x,j) - xyzq_i.x;
            float dy = __shfl(xyzq_j.y,j) - xyzq_i.y;
            float dz = __shfl(xyzq_j.z,j) - xyzq_i.z;

            float r2 = dx*dx + dy*dy + dz*dz;

            if (j < WARPSIZE && r2 < plcutoff2) {
              // We have atom pair within the pairlist cutoff => Set indicator bit
              nexcluded |= 1;
            }
            shuffleNext<doPairlist>(xyzq_j.w, vdwtypej, jatomIndex);
          }

          for (;t < WARPSIZE;t++) {
            int j = (t + wid) & modval;

            // NOTE: __shfl() operation can give non-sense here because j may be >= WARPSIZE.
            //       However, if (j < WARPSIZE ..) below makes sure that these non-sense
            //       results are not used
            float dx = __shfl(xyzq_j.x,j) - xyzq_i.x;
            float dy = __shfl(xyzq_j.y,j) - xyzq_i.y;
            float dz = __shfl(xyzq_j.z,j) - xyzq_i.z;

            float r2 = dx*dx + dy*dy + dz*dz;

            excl >>= 1;
            if (j < WARPSIZE && r2 < plcutoff2) {
              // We have atom pair within the pairlist cutoff => Set indicator bit
              nexcluded |= 1;
              if (j < nfreej || wid < nfreei) {
                bool excluded = false;
                int indexdiff = jatomIndex - iatomIndex;
                if ( abs(indexdiff) <= iexclMaxdiff) {
                  indexdiff += iexclIndex;
                  int indexword = ((unsigned int) indexdiff) >> 5;

                  if ( indexword < MAX_CONST_EXCLUSIONS ) {
                    indexword = constExclusions[indexword];
                  } else {
                    indexword = overflowExclusions[indexword];
                  }

                  excluded = ((indexword & (1<<(indexdiff&31))) != 0);
                }
                if (excluded) nexcluded += 2;
                if (!excluded) excl |= 0x80000000;
                if (!excluded && r2 < cutoff2) {
                  calcForceEnergy<doEnergy, doSlow>(r2, xyzq_i.w, xyzq_j.w, dx, dy, dz,
                    vdwtypei, vdwtypej,
                    vdwCoefTable,
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
                    forceTableTex, energyTableTex,
#endif
                    iforce, iforceSlow, jforce, jforceSlow, energyVdw, energyElec, energySlow);
                }
              }
            }
            shuffleNext<doPairlist>(xyzq_j.w, vdwtypej, jatomIndex);
            shuffleNext<doSlow>(jforce, jforceSlow);
          } // t
        } else {
          // Just compute forces
/*          
          if (self) {
            shuffleNext<doPairlist>(xyzq_j.w, vdwtypej, jatomIndex);
            excl >>= 1;
          }
          for (;t < WARPSIZE;t++) {
            int j = (t + wid) & (WARPSIZE-1);

            float dx = __shfl(xyzq_j.x,j) - xyzq_i.x;
            float dy = __shfl(xyzq_j.y,j) - xyzq_i.y;
            float dz = __shfl(xyzq_j.z,j) - xyzq_i.z;

            float r2 = dx*dx + dy*dy + dz*dz;

            if ((excl & 1) && r2 < cutoff2) {
              calcForceEnergy<doEnergy, doSlow>(r2, xyzq_i.w, xyzq_j.w, dx, dy, dz,
                vdwtypei, vdwtypej, vdwCoefTable,
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
                forceTableTex, energyTableTex,
#endif
                iforce, iforceSlow, jforce, jforceSlow, energyVdw, energyElec, energySlow);
            } // ((excl & 1) && r2 < cutoff2)
            excl >>= 1;
            shuffleNext<doPairlist>(xyzq_j.w, vdwtypej, jatomIndex);
            shuffleNext<doSlow>(jforce, jforceSlow);
          } // t
*/
          if (self) {
            excl >>= 1;
            xyzq_j.x = __shfl(xyzq_j.x, (threadIdx.x+1) & (WARPSIZE-1));
            xyzq_j.y = __shfl(xyzq_j.y, (threadIdx.x+1) & (WARPSIZE-1));
            xyzq_j.z = __shfl(xyzq_j.z, (threadIdx.x+1) & (WARPSIZE-1));
            shuffleNext<doPairlist>(xyzq_j.w, vdwtypej, jatomIndex);
          }
          for (;t < WARPSIZE;t++) {
            if ((excl & 1)) {
              float dx = xyzq_j.x - xyzq_i.x;
              float dy = xyzq_j.y - xyzq_i.y;
              float dz = xyzq_j.z - xyzq_i.z;

              float r2 = dx*dx + dy*dy + dz*dz;

              if (r2 < cutoff2) {
                calcForceEnergy<doEnergy, doSlow>(r2, xyzq_i.w, xyzq_j.w, dx, dy, dz,
                  vdwtypei, vdwtypej,
                  vdwCoefTable,
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
                  forceTableTex, energyTableTex,
#endif
                  iforce, iforceSlow, jforce, jforceSlow, energyVdw, energyElec, energySlow);
              } // (r2 < cutoff2)
            } // (excl & 1)
            excl >>= 1;
            xyzq_j.x = __shfl(xyzq_j.x, (threadIdx.x+1) & (WARPSIZE-1));
            xyzq_j.y = __shfl(xyzq_j.y, (threadIdx.x+1) & (WARPSIZE-1));
            xyzq_j.z = __shfl(xyzq_j.z, (threadIdx.x+1) & (WARPSIZE-1));
            shuffleNext<doPairlist>(xyzq_j.w, vdwtypej, jatomIndex);
            shuffleNext<doSlow>(jforce, jforceSlow);
          } // t
        }

        // Write j-forces
        storeForces<doSlow>(jatomStart + wid, jforce, jforceSlow, devForces, devForcesSlow);

        // Write exclusions
        if (doPairlist && __any(nexcluded & 1)) {
          int anyexcl = (65536 | __any(excl));
          // Mark this jtile as non-empty:
          //  VdW:      1 if tile has atom pairs within pairlist cutoff and some these atoms interact
          //  GBIS: 65536 if tile has atom pairs within pairlist cutoff but not necessary interacting (i.e. these atoms are fixed or excluded)
          if (wid == 0) jtiles[jtile] = anyexcl;
          // Store exclusions
          tileExcls[jtile].excl[wid] = excl;
          // itileListLen:
          // lower 16 bits number of tiles with atom pairs within pairlist cutoff that interact
          // upper 16 bits number of tiles with atom pairs within pairlist cutoff (but not necessary interacting)
          itileListLen += anyexcl;
          // NOTE, this minJatomStart is only stored once for the first tile list entry
          // minJatomStart = min(minJatomStart, jatomStart);
        }

      } // jtile

      // Write i-forces
      storeForces<doSlow>(iatomStart + wid, iforce, iforceSlow, devForces, devForcesSlow);
    }
    // Done with computation

    // Save pairlist stuff
    if (doPairlist) {

      // Warp index (0...warpsize-1)
      const int wid = threadIdx.x % WARPSIZE;

      if (wid == 0) {
        // minJatomStart is in range [0 ... atomStorageSize-1]
        //int atom0 = (minJatomStart)/WARPSIZE;
        // int atom0 = 0;
        // int storageOffset = atomStorageSize/WARPSIZE;
        // int itileListLen = 0;
        // for (int jtile=jtileStart;jtile <= jtileEnd;jtile++) itileListLen += jtiles[jtile];
        // Store 0 if itileListLen == 0
        // tileListDepth[itileList] = (itileListLen > 0)*(itileListLen*storageOffset + atom0);
        tileListDepth[itileList] = itileListLen;
        tileListOrder[itileList] = itileList;
        // Number of active tilelists with tile with atom pairs within pairlist cutoff that interact
        if ((itileListLen & 65535) > 0) atomicAdd(&tileListStat->numTileLists, 1);
        // Number of active tilelists with tiles with atom pairs within pairlist cutoff (but not necessary interacting)
        if (itileListLen > 0) atomicAdd(&tileListStat->numTileListsGBIS, 1);
        // NOTE: always numTileListsGBIS >= numTileLists
      }

      typedef cub::WarpReduce<int> WarpReduceInt;
      __shared__ typename WarpReduceInt::TempStorage tempStorage[NONBONDKERNEL_NUM_WARP];
      int warpId = threadIdx.x / WARPSIZE;
      // Remove indicator bit
      nexcluded >>= 1;
      volatile int nexcludedWarp = WarpReduceInt(tempStorage[warpId]).Sum(nexcluded);
      if (wid == 0) atomicAdd(&tileListStat->numExcluded, nexcludedWarp);

    }

    if (doVirial) {
      // Warp index (0...warpsize-1)
      const int wid = threadIdx.x % WARPSIZE;

      typedef cub::WarpReduce<float> WarpReduce;
      __shared__ typename WarpReduce::TempStorage tempStorage[NONBONDKERNEL_NUM_WARP];
      int warpId = threadIdx.x / WARPSIZE;
      volatile float iforcexSum = WarpReduce(tempStorage[warpId]).Sum(iforce.x);
      volatile float iforceySum = WarpReduce(tempStorage[warpId]).Sum(iforce.y);
      volatile float iforcezSum = WarpReduce(tempStorage[warpId]).Sum(iforce.z);
      if (wid == 0) {
        virialEnergy[itileList].forcex = iforcexSum;
        virialEnergy[itileList].forcey = iforceySum;
        virialEnergy[itileList].forcez = iforcezSum;
      }

      if (doSlow) {
        iforcexSum = WarpReduce(tempStorage[warpId]).Sum(iforceSlow.x);
        iforceySum = WarpReduce(tempStorage[warpId]).Sum(iforceSlow.y);
        iforcezSum = WarpReduce(tempStorage[warpId]).Sum(iforceSlow.z);
        if (wid == 0) {
          virialEnergy[itileList].forceSlowx = iforcexSum;
          virialEnergy[itileList].forceSlowy = iforceySum;
          virialEnergy[itileList].forceSlowz = iforcezSum;
        }
      }
    }

    // Reduce energy
    if (doEnergy) {
      // NOTE: We must hand write these warp-wide reductions to avoid excess register spillage
      //       (Why does CUB suck here?)
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        energyVdw += __shfl_xor(energyVdw, i, 32);
        energyElec += __shfl_xor(energyElec, i, 32);
        if (doSlow) energySlow += __shfl_xor(energySlow, i, 32);
      }

      if (threadIdx.x % WARPSIZE == 0) {
        virialEnergy[itileList].energyVdw  = energyVdw;
        virialEnergy[itileList].energyElec = energyElec;
        if (doSlow) virialEnergy[itileList].energySlow = energySlow;
      }
    }

    if (doStreaming) {
      // Make sure devForces and devForcesSlow have been written into device memory
      // NO NEED TO SYNCHRONIZE THREADS, THIS IS WARP-LEVEL
      __threadfence();

      int patchDone[2] = {false, false};
      const int wid = threadIdx.x % WARPSIZE;
      if (wid == 0) {
        int patchCountOld0 = atomicInc(&patchNumCount[patchInd.x], (unsigned int)(patchNumList.x-1));
        patchDone[0] = (patchCountOld0 + 1 == patchNumList.x);
        if (patchInd.x != patchInd.y) {
          int patchCountOld1 = atomicInc(&patchNumCount[patchInd.y], (unsigned int)(patchNumList.y-1));
          patchDone[1] = (patchCountOld1 + 1 == patchNumList.y);
        }
      }

      patchDone[0] = __any(patchDone[0]);
      patchDone[1] = __any(patchDone[1]);

      if (patchDone[0]) {
        // Patch 1 is done, write onto host-mapped memory
        CudaPatchRecord patch = cudaPatches[patchInd.x];
        int start = patch.atomStart;
        int end   = start + patch.numAtoms;
        for (int i=start+wid;i < end;i+=WARPSIZE) {
          mapForces[i] = devForces[i];
          if (doSlow) mapForcesSlow[i] = devForcesSlow[i];
        }
      }
      if (patchDone[1]) {
        // Patch 2 is done
        CudaPatchRecord patch = cudaPatches[patchInd.y];
        int start = patch.atomStart;
        int end   = start + patch.numAtoms;
        for (int i=start+wid;i < end;i+=WARPSIZE) {
          mapForces[i] = devForces[i];
          if (doSlow) mapForcesSlow[i] = devForcesSlow[i];
        }
      }

      if (patchDone[0] || patchDone[1]) {
        // Make sure mapForces and mapForcesSlow are up-to-date
        __threadfence_system();
        // Add patch into "patchReadyQueue"
        if (wid == 0) {
          if (patchDone[0]) {
            int ind = atomicAdd(&tileListStat->patchReadyQueueCount, 1);
            // int ind = atomicInc((unsigned int *)&mapPatchReadyQueue[numPatches], numPatches-1);
            mapPatchReadyQueue[ind] = patchInd.x;
          }
          if (patchDone[1]) {
            int ind = atomicAdd(&tileListStat->patchReadyQueueCount, 1);
            // int ind = atomicInc((unsigned int *)&mapPatchReadyQueue[numPatches], numPatches-1);
            mapPatchReadyQueue[ind] = patchInd.y;
          }
        }
        // Make sure "patchReadyQueue" is visible in page-locked host memory
        __threadfence_system();
      }
    }

    if (doStreaming && outputOrder != NULL && threadIdx.x % WARPSIZE == 0) {
      int index = atomicAdd(&tileListStat->outputOrderIndex, 1);
      outputOrder[index] = itileList;
    }
  } // if (itileList < numTileLists)
}

//
// Finish up - reduce virials from nonbonded kernel
//
#define REDUCENONBONDEDVIRIALKERNEL_NUM_WARP 32
__global__ void reduceNonbondedVirialKernel(const bool doSlow,
  const int atomStorageSize,
  const float4* __restrict__ xyzq,
  const float4* __restrict__ devForces, const float4* __restrict__ devForcesSlow,
  VirialEnergy* __restrict__ virialEnergy) {

  for (int ibase = blockIdx.x*blockDim.x;ibase < atomStorageSize;ibase += blockDim.x*gridDim.x)
  {
    int i = ibase + threadIdx.x;

    // Set to zero to avoid nan*0
    float4 pos;
    pos.x = 0.0f;
    pos.y = 0.0f;
    pos.z = 0.0f;
    float4 force, forceSlow;
    force.x = 0.0f;
    force.y = 0.0f;
    force.z = 0.0f;
    forceSlow.x = 0.0f;
    forceSlow.y = 0.0f;
    forceSlow.z = 0.0f;
    if (i < atomStorageSize) {
      pos = xyzq[i];
      force = devForces[i];
      if (doSlow) forceSlow = devForcesSlow[i];
    }
    // Reduce across the entire thread block
    float vxxt = force.x*pos.x;
    float vxyt = force.x*pos.y;
    float vxzt = force.x*pos.z;
    float vyxt = force.y*pos.x;
    float vyyt = force.y*pos.y;
    float vyzt = force.y*pos.z;
    float vzxt = force.z*pos.x;
    float vzyt = force.z*pos.y;
    float vzzt = force.z*pos.z;
    // atomicAdd(&virialEnergy->virial[0], (double)vxx);
    // atomicAdd(&virialEnergy->virial[1], (double)vxy);
    // atomicAdd(&virialEnergy->virial[2], (double)vxz);
    // atomicAdd(&virialEnergy->virial[3], (double)vyx);
    // atomicAdd(&virialEnergy->virial[4], (double)vyy);
    // atomicAdd(&virialEnergy->virial[5], (double)vyz);
    // atomicAdd(&virialEnergy->virial[6], (double)vzx);
    // atomicAdd(&virialEnergy->virial[7], (double)vzy);
    // atomicAdd(&virialEnergy->virial[8], (double)vzz);

    typedef cub::BlockReduce<float, REDUCENONBONDEDVIRIALKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
    __shared__ typename BlockReduce::TempStorage tempStorage;
    volatile float vxx = BlockReduce(tempStorage).Sum(vxxt); __syncthreads();
    volatile float vxy = BlockReduce(tempStorage).Sum(vxyt); __syncthreads();
    volatile float vxz = BlockReduce(tempStorage).Sum(vxzt); __syncthreads();
    volatile float vyx = BlockReduce(tempStorage).Sum(vyxt); __syncthreads();
    volatile float vyy = BlockReduce(tempStorage).Sum(vyyt); __syncthreads();
    volatile float vyz = BlockReduce(tempStorage).Sum(vyzt); __syncthreads();
    volatile float vzx = BlockReduce(tempStorage).Sum(vzxt); __syncthreads();
    volatile float vzy = BlockReduce(tempStorage).Sum(vzyt); __syncthreads();
    volatile float vzz = BlockReduce(tempStorage).Sum(vzzt); __syncthreads();
    if (threadIdx.x == 0) {
      atomicAdd(&virialEnergy->virial[0], (double)vxx);
      atomicAdd(&virialEnergy->virial[1], (double)vxy);
      atomicAdd(&virialEnergy->virial[2], (double)vxz);
      atomicAdd(&virialEnergy->virial[3], (double)vyx);
      atomicAdd(&virialEnergy->virial[4], (double)vyy);
      atomicAdd(&virialEnergy->virial[5], (double)vyz);
      atomicAdd(&virialEnergy->virial[6], (double)vzx);
      atomicAdd(&virialEnergy->virial[7], (double)vzy);
      atomicAdd(&virialEnergy->virial[8], (double)vzz);
    }

    if (doSlow) {
      // if (isnan(forceSlow.x) || isnan(forceSlow.y) || isnan(forceSlow.z))
      float vxxSlowt = forceSlow.x*pos.x;
      float vxySlowt = forceSlow.x*pos.y;
      float vxzSlowt = forceSlow.x*pos.z;
      float vyxSlowt = forceSlow.y*pos.x;
      float vyySlowt = forceSlow.y*pos.y;
      float vyzSlowt = forceSlow.y*pos.z;
      float vzxSlowt = forceSlow.z*pos.x;
      float vzySlowt = forceSlow.z*pos.y;
      float vzzSlowt = forceSlow.z*pos.z;
      // atomicAdd(&virialEnergy->virialSlow[0], (double)vxxSlow);
      // atomicAdd(&virialEnergy->virialSlow[1], (double)vxySlow);
      // atomicAdd(&virialEnergy->virialSlow[2], (double)vxzSlow);
      // atomicAdd(&virialEnergy->virialSlow[3], (double)vyxSlow);
      // atomicAdd(&virialEnergy->virialSlow[4], (double)vyySlow);
      // atomicAdd(&virialEnergy->virialSlow[5], (double)vyzSlow);
      // atomicAdd(&virialEnergy->virialSlow[6], (double)vzxSlow);
      // atomicAdd(&virialEnergy->virialSlow[7], (double)vzySlow);
      // atomicAdd(&virialEnergy->virialSlow[8], (double)vzzSlow);
      volatile float vxxSlow = BlockReduce(tempStorage).Sum(vxxSlowt); __syncthreads();
      volatile float vxySlow = BlockReduce(tempStorage).Sum(vxySlowt); __syncthreads();
      volatile float vxzSlow = BlockReduce(tempStorage).Sum(vxzSlowt); __syncthreads();
      volatile float vyxSlow = BlockReduce(tempStorage).Sum(vyxSlowt); __syncthreads();
      volatile float vyySlow = BlockReduce(tempStorage).Sum(vyySlowt); __syncthreads();
      volatile float vyzSlow = BlockReduce(tempStorage).Sum(vyzSlowt); __syncthreads();
      volatile float vzxSlow = BlockReduce(tempStorage).Sum(vzxSlowt); __syncthreads();
      volatile float vzySlow = BlockReduce(tempStorage).Sum(vzySlowt); __syncthreads();
      volatile float vzzSlow = BlockReduce(tempStorage).Sum(vzzSlowt); __syncthreads();
      if (threadIdx.x == 0) {
        atomicAdd(&virialEnergy->virialSlow[0], (double)vxxSlow);
        atomicAdd(&virialEnergy->virialSlow[1], (double)vxySlow);
        atomicAdd(&virialEnergy->virialSlow[2], (double)vxzSlow);
        atomicAdd(&virialEnergy->virialSlow[3], (double)vyxSlow);
        atomicAdd(&virialEnergy->virialSlow[4], (double)vyySlow);
        atomicAdd(&virialEnergy->virialSlow[5], (double)vyzSlow);
        atomicAdd(&virialEnergy->virialSlow[6], (double)vzxSlow);
        atomicAdd(&virialEnergy->virialSlow[7], (double)vzySlow);
        atomicAdd(&virialEnergy->virialSlow[8], (double)vzzSlow);
      }
    }
  
  }

}

#define REDUCEVIRIALENERGYKERNEL_NUM_WARP 32
__global__ void reduceVirialEnergyKernel(
  const bool doEnergy, const bool doVirial, const bool doSlow,
  const int numTileLists,
  const TileListVirialEnergy* __restrict__ tileListVirialEnergy,
  VirialEnergy* __restrict__ virialEnergy) {

  for (int ibase = blockIdx.x*blockDim.x;ibase < numTileLists;ibase += blockDim.x*gridDim.x)
  {
    int itileList = ibase + threadIdx.x;
    TileListVirialEnergy ve;
    if (itileList < numTileLists) {
      ve = tileListVirialEnergy[itileList];
    } else {
      // Set to zero to avoid nan*0
      if (doVirial) {
        ve.shx = 0.0f;
        ve.shy = 0.0f;
        ve.shz = 0.0f;
        ve.forcex = 0.0f;
        ve.forcey = 0.0f;
        ve.forcez = 0.0f;
        ve.forceSlowx = 0.0f;
        ve.forceSlowy = 0.0f;
        ve.forceSlowz = 0.0f;
      }
      if (doEnergy) {
        ve.energyVdw = 0.0;
        ve.energyElec = 0.0;
        ve.energySlow = 0.0;
        // ve.energyGBIS = 0.0;
      }
    }

    if (doVirial) {
      typedef cub::BlockReduce<float, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
      __shared__ typename BlockReduce::TempStorage tempStorage;
      float vxxt = ve.forcex*ve.shx;
      float vxyt = ve.forcex*ve.shy;
      float vxzt = ve.forcex*ve.shz;
      float vyxt = ve.forcey*ve.shx;
      float vyyt = ve.forcey*ve.shy;
      float vyzt = ve.forcey*ve.shz;
      float vzxt = ve.forcez*ve.shx;
      float vzyt = ve.forcez*ve.shy;
      float vzzt = ve.forcez*ve.shz;
      volatile float vxx = BlockReduce(tempStorage).Sum(vxxt); __syncthreads();
      volatile float vxy = BlockReduce(tempStorage).Sum(vxyt); __syncthreads();
      volatile float vxz = BlockReduce(tempStorage).Sum(vxzt); __syncthreads();
      volatile float vyx = BlockReduce(tempStorage).Sum(vyxt); __syncthreads();
      volatile float vyy = BlockReduce(tempStorage).Sum(vyyt); __syncthreads();
      volatile float vyz = BlockReduce(tempStorage).Sum(vyzt); __syncthreads();
      volatile float vzx = BlockReduce(tempStorage).Sum(vzxt); __syncthreads();
      volatile float vzy = BlockReduce(tempStorage).Sum(vzyt); __syncthreads();
      volatile float vzz = BlockReduce(tempStorage).Sum(vzzt); __syncthreads();
      if (threadIdx.x == 0) {
        atomicAdd(&virialEnergy->virial[0], (double)vxx);
        atomicAdd(&virialEnergy->virial[1], (double)vxy);
        atomicAdd(&virialEnergy->virial[2], (double)vxz);
        atomicAdd(&virialEnergy->virial[3], (double)vyx);
        atomicAdd(&virialEnergy->virial[4], (double)vyy);
        atomicAdd(&virialEnergy->virial[5], (double)vyz);
        atomicAdd(&virialEnergy->virial[6], (double)vzx);
        atomicAdd(&virialEnergy->virial[7], (double)vzy);
        atomicAdd(&virialEnergy->virial[8], (double)vzz);
      }

      if (doSlow) {
        typedef cub::BlockReduce<float, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
        __shared__ typename BlockReduce::TempStorage tempStorage;
        float vxxt = ve.forceSlowx*ve.shx;
        float vxyt = ve.forceSlowx*ve.shy;
        float vxzt = ve.forceSlowx*ve.shz;
        float vyxt = ve.forceSlowy*ve.shx;
        float vyyt = ve.forceSlowy*ve.shy;
        float vyzt = ve.forceSlowy*ve.shz;
        float vzxt = ve.forceSlowz*ve.shx;
        float vzyt = ve.forceSlowz*ve.shy;
        float vzzt = ve.forceSlowz*ve.shz;
        volatile float vxx = BlockReduce(tempStorage).Sum(vxxt); __syncthreads();
        volatile float vxy = BlockReduce(tempStorage).Sum(vxyt); __syncthreads();
        volatile float vxz = BlockReduce(tempStorage).Sum(vxzt); __syncthreads();
        volatile float vyx = BlockReduce(tempStorage).Sum(vyxt); __syncthreads();
        volatile float vyy = BlockReduce(tempStorage).Sum(vyyt); __syncthreads();
        volatile float vyz = BlockReduce(tempStorage).Sum(vyzt); __syncthreads();
        volatile float vzx = BlockReduce(tempStorage).Sum(vzxt); __syncthreads();
        volatile float vzy = BlockReduce(tempStorage).Sum(vzyt); __syncthreads();
        volatile float vzz = BlockReduce(tempStorage).Sum(vzzt); __syncthreads();
        if (threadIdx.x == 0) {
          atomicAdd(&virialEnergy->virialSlow[0], (double)vxx);
          atomicAdd(&virialEnergy->virialSlow[1], (double)vxy);
          atomicAdd(&virialEnergy->virialSlow[2], (double)vxz);
          atomicAdd(&virialEnergy->virialSlow[3], (double)vyx);
          atomicAdd(&virialEnergy->virialSlow[4], (double)vyy);
          atomicAdd(&virialEnergy->virialSlow[5], (double)vyz);
          atomicAdd(&virialEnergy->virialSlow[6], (double)vzx);
          atomicAdd(&virialEnergy->virialSlow[7], (double)vzy);
          atomicAdd(&virialEnergy->virialSlow[8], (double)vzz);
        }
      }
    }

    if (doEnergy) {
      typedef cub::BlockReduce<double, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
      __shared__ typename BlockReduce::TempStorage tempStorage;
      volatile double energyVdw  = BlockReduce(tempStorage).Sum(ve.energyVdw); __syncthreads();
      volatile double energyElec = BlockReduce(tempStorage).Sum(ve.energyElec); __syncthreads();
      if (threadIdx.x == 0) {
          atomicAdd(&virialEnergy->energyVdw, (double)energyVdw);
          atomicAdd(&virialEnergy->energyElec, (double)energyElec);
      }
      if (doSlow) {
        volatile double energySlow = BlockReduce(tempStorage).Sum(ve.energySlow); __syncthreads();
        if (threadIdx.x == 0) atomicAdd(&virialEnergy->energySlow, (double)energySlow);
      }
      // if (doGBIS) {
      //   double energyGBIS = BlockReduce(tempStorage).Sum(ve.energyGBIS); __syncthreads();
      //   if (threadIdx.x == 0) atomicAdd(&virialEnergy->energyGBIS, (double)energyGBIS);
      // }
    }

  }

}

#define REDUCEGBISENERGYKERNEL_NUM_WARP 32
__global__ void reduceGBISEnergyKernel(const int numTileLists,
  const TileListVirialEnergy* __restrict__ tileListVirialEnergy,
  VirialEnergy* __restrict__ virialEnergy) {

  for (int ibase = blockIdx.x*blockDim.x;ibase < numTileLists;ibase += blockDim.x*gridDim.x)
  {
    int itileList = ibase + threadIdx.x;
    double energyGBISt = 0.0;
    if (itileList < numTileLists) {
      energyGBISt = tileListVirialEnergy[itileList].energyGBIS;
    }

    typedef cub::BlockReduce<double, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
    __shared__ typename BlockReduce::TempStorage tempStorage;
    volatile double energyGBIS = BlockReduce(tempStorage).Sum(energyGBISt); __syncthreads();
    if (threadIdx.x == 0) atomicAdd(&virialEnergy->energyGBIS, (double)energyGBIS);
  }

}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

CudaComputeNonbondedKernel::CudaComputeNonbondedKernel(int deviceID, bool doStreaming) : deviceID(deviceID), doStreaming(doStreaming) {
  cudaCheck(cudaSetDevice(deviceID));

  overflowExclusions = NULL;
  overflowExclusionsSize = 0;

  exclIndexMaxDiff = NULL;
  exclIndexMaxDiffSize = 0;

  atomIndex = NULL;
  atomIndexSize = 0;

  vdwTypes = NULL;
  vdwTypesSize = 0;

  vdwCoefTable = NULL;
  vdwCoefTableWidth = 0;
  vdwCoefTableSize = 0;

#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  forceTableTexActive = false;
  energyTableTexActive = false;
#endif

  vdwCoefTableTextureBound = false;
#ifdef DISABLE_CUDA_TEXTURE_OBJECTS
  force_table_bound = false;
  energy_table_bound = false;
#endif

  patchNumCount = NULL;
  patchNumCountSize = 0;

  patchReadyQueue = NULL;
  patchReadyQueueSize = 0;

}

CudaComputeNonbondedKernel::~CudaComputeNonbondedKernel() {
  cudaCheck(cudaSetDevice(deviceID));
  if (overflowExclusions != NULL) deallocate_device<unsigned int>(&overflowExclusions);
  if (exclIndexMaxDiff != NULL) deallocate_device<int2>(&exclIndexMaxDiff);
  if (atomIndex != NULL) deallocate_device<int>(&atomIndex);
  if (vdwTypes != NULL) deallocate_device<int>(&vdwTypes);
  if (vdwCoefTable != NULL) deallocate_device<float2>(&vdwCoefTable);
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  if (forceTableTexActive) cudaCheck(cudaDestroyTextureObject(forceTableTex));
  if (energyTableTexActive) cudaCheck(cudaDestroyTextureObject(energyTableTex));
#else
  if (force_table_bound) cudaCheck(cudaUnbindTexture(force_table));
  if (energy_table_bound) cudaCheck(cudaUnbindTexture(energy_table));
#endif
  if (vdwCoefTableTextureBound) cudaCheck(cudaUnbindTexture(vdwCoefTableTexture));
  if (patchNumCount != NULL) deallocate_device<unsigned int>(&patchNumCount);
  if (patchReadyQueue != NULL) deallocate_host<int>(&patchReadyQueue);
}

void CudaComputeNonbondedKernel::updateVdwTypesExcl(const int atomStorageSize, const int* h_vdwTypes,
  const int2* h_exclIndexMaxDiff, const int* h_atomIndex, cudaStream_t stream) {

  reallocate_device<int>(&vdwTypes, &vdwTypesSize, atomStorageSize, OVERALLOC);
  reallocate_device<int2>(&exclIndexMaxDiff, &exclIndexMaxDiffSize, atomStorageSize, OVERALLOC);
  reallocate_device<int>(&atomIndex, &atomIndexSize, atomStorageSize, OVERALLOC);

  copy_HtoD<int>(h_vdwTypes, vdwTypes, atomStorageSize, stream);
  copy_HtoD<int2>(h_exclIndexMaxDiff, exclIndexMaxDiff, atomStorageSize, stream);
  copy_HtoD<int>(h_atomIndex, atomIndex, atomStorageSize, stream);
}

int* CudaComputeNonbondedKernel::getPatchReadyQueue() {
  if (!doStreaming) {
    NAMD_die("CudaComputeNonbondedKernel::getPatchReadyQueue() called on non-streaming kernel");
  }
  return patchReadyQueue;
}

void CudaComputeNonbondedKernel::nonbondedForce(CudaTileListKernel& tlKernel,
  const int atomStorageSize, const bool doPairlist,
  const bool doEnergy, const bool doVirial, const bool doSlow,
  const float latticeX, const float latticeY, const float latticeZ,
  const float4* h_xyzq, const float cutoff2, 
  float4* d_forces, float4* d_forcesSlow,
  float4* h_forces, float4* h_forcesSlow,
  cudaStream_t stream) {

  if (!doPairlist) copy_HtoD<float4>(h_xyzq, tlKernel.get_xyzq(), atomStorageSize, stream);

  clear_device_array<float4>(d_forces, atomStorageSize, stream);
  if (doSlow) clear_device_array<float4>(d_forcesSlow, atomStorageSize, stream);

  tlKernel.clearTileListStat(stream);
  // clear_device_array<TileListStat>(tlKernel.getTileListStatDevPtr(), 1, stream);

  // --- streaming ----
  float4* m_forces = NULL;
  float4* m_forcesSlow = NULL;
  int* m_patchReadyQueue = NULL;
  int numPatches = 0;
  unsigned int* patchNumCountPtr = NULL;
  if (doStreaming) {
    numPatches = tlKernel.getNumPatches();
    if (reallocate_device<unsigned int>(&patchNumCount, &patchNumCountSize, numPatches)) {
      // If re-allocated, clear array
      clear_device_array<unsigned int>(patchNumCount, numPatches, stream);
    }
    patchNumCountPtr = patchNumCount;
    bool re = reallocate_host<int>(&patchReadyQueue, &patchReadyQueueSize, numPatches, cudaHostAllocMapped);
    if (re) {
      // If re-allocated, re-set to "-1"
      for (int i=0;i < numPatches;i++) patchReadyQueue[i] = -1;
    }
    cudaCheck(cudaHostGetDevicePointer(&m_patchReadyQueue, patchReadyQueue, 0));
    cudaCheck(cudaHostGetDevicePointer(&m_forces, h_forces, 0));
    cudaCheck(cudaHostGetDevicePointer(&m_forcesSlow, h_forcesSlow, 0));
  }
  // -----------------

  if (doVirial || doEnergy) {
    tlKernel.setTileListVirialEnergyLength(tlKernel.getNumTileLists());
  }

  int shMemSize = 0;

  int* outputOrderPtr = tlKernel.getOutputOrder();

  int nwarp = NONBONDKERNEL_NUM_WARP;
  int nthread = WARPSIZE*nwarp;
  int start = 0;
  while (start < tlKernel.getNumTileLists())
  {

    int nleft = tlKernel.getNumTileLists() - start;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (nleft-1)/nwarp+1);

#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
#define CALL(DOENERGY, DOVIRIAL, DOSLOW, DOPAIRLIST, DOSTREAMING) \
    nonbondedForceKernel<DOENERGY, DOVIRIAL, DOSLOW, DOPAIRLIST, DOSTREAMING> \
  <<< nblock, nthread, shMemSize, stream >>>  \
  (start, tlKernel.getNumTileLists(), tlKernel.getTileLists(), tlKernel.getTileExcls(), tlKernel.getTileJatomStart(), \
    vdwCoefTableWidth, vdwCoefTable, vdwTypes, latticeX, latticeY, latticeZ, tlKernel.get_xyzq(), cutoff2, \
    forceTableTex, energyTableTex, \
    atomStorageSize, tlKernel.get_plcutoff2(), tlKernel.getPatchPairs(), atomIndex, exclIndexMaxDiff, overflowExclusions, \
    tlKernel.getTileListDepth(), tlKernel.getTileListOrder(), tlKernel.getJtiles(), tlKernel.getTileListStatDevPtr(), \
    tlKernel.getBoundingBoxes(), d_forces, d_forcesSlow, \
    numPatches, patchNumCountPtr, tlKernel.getCudaPatches(), m_forces, m_forcesSlow, m_patchReadyQueue, \
    outputOrderPtr, tlKernel.getTileListVirialEnergy()); called=true
#else  // DISABLE_CUDA_TEXTURE_OBJECTS
#define CALL(DOENERGY, DOVIRIAL, DOSLOW, DOPAIRLIST, DOSTREAMING) \
  nonbondedForceKernel<DOENERGY, DOVIRIAL, DOSLOW, DOPAIRLIST, DOSTREAMING> \
  <<< nblock, nthread, shMemSize, stream >>>  \
  (start, tlKernel.getNumTileLists(), tlKernel.getTileLists(), tlKernel.getTileExcls(), tlKernel.getTileJatomStart(), \
    vdwCoefTableWidth, vdwCoefTable, vdwTypes, latticeX, latticeY, latticeZ, tlKernel.get_xyzq(), cutoff2, \
    atomStorageSize, tlKernel.get_plcutoff2(), tlKernel.getPatchPairs(), atomIndex, exclIndexMaxDiff, overflowExclusions, \
    tlKernel.getTileListDepth(), tlKernel.getTileListOrder(), tlKernel.getJtiles(), tlKernel.getTileListStatDevPtr(), \
    tlKernel.getBoundingBoxes(), d_forces, d_forcesSlow, \
    numPatches, patchNumCountPtr, tlKernel.getCudaPatches(), m_forces, m_forcesSlow, m_patchReadyQueue, \
    outputOrderPtr, tlKernel.getTileListVirialEnergy()); called=true
#endif // DISABLE_CUDA_TEXTURE_OBJECTS

    bool called = false;

    if (doStreaming) {
      if (!doEnergy && !doVirial && !doSlow && !doPairlist) CALL(0, 0, 0, 0, 1);
      if (!doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(0, 0, 1, 0, 1);
      if (!doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(0, 1, 0, 0, 1);
      if (!doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(0, 1, 1, 0, 1);
      if ( doEnergy && !doVirial && !doSlow && !doPairlist) CALL(1, 0, 0, 0, 1);
      if ( doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(1, 0, 1, 0, 1);
      if ( doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(1, 1, 0, 0, 1);
      if ( doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(1, 1, 1, 0, 1);

      if (!doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(0, 0, 0, 1, 1);
      if (!doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(0, 0, 1, 1, 1);
      if (!doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(0, 1, 0, 1, 1);
      if (!doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(0, 1, 1, 1, 1);
      if ( doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(1, 0, 0, 1, 1);
      if ( doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(1, 0, 1, 1, 1);
      if ( doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(1, 1, 0, 1, 1);
      if ( doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(1, 1, 1, 1, 1);
    } else {
      if (!doEnergy && !doVirial && !doSlow && !doPairlist) CALL(0, 0, 0, 0, 0);
      if (!doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(0, 0, 1, 0, 0);
      if (!doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(0, 1, 0, 0, 0);
      if (!doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(0, 1, 1, 0, 0);
      if ( doEnergy && !doVirial && !doSlow && !doPairlist) CALL(1, 0, 0, 0, 0);
      if ( doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(1, 0, 1, 0, 0);
      if ( doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(1, 1, 0, 0, 0);
      if ( doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(1, 1, 1, 0, 0);

      if (!doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(0, 0, 0, 1, 0);
      if (!doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(0, 0, 1, 1, 0);
      if (!doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(0, 1, 0, 1, 0);
      if (!doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(0, 1, 1, 1, 0);
      if ( doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(1, 0, 0, 1, 0);
      if ( doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(1, 0, 1, 1, 0);
      if ( doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(1, 1, 0, 1, 0);
      if ( doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(1, 1, 1, 1, 0);
    }

    if (!called) {
      NAMD_die("CudaComputeNonbondedKernel::nonbondedForce, none of the kernels called");
    }

#undef CALL
    cudaCheck(cudaGetLastError());

    start += nblock*nwarp;
  }

}

//
// Perform virial and energy reductions for non-bonded force calculation
//
void CudaComputeNonbondedKernel::reduceVirialEnergy(CudaTileListKernel& tlKernel,
  const int atomStorageSize, const bool doEnergy, const bool doVirial, const bool doSlow, const bool doGBIS,
  float4* d_forces, float4* d_forcesSlow,
  VirialEnergy* d_virialEnergy, cudaStream_t stream) {

  if (doEnergy || doVirial) {
    clear_device_array<VirialEnergy>(d_virialEnergy, 1, stream);
  }

  if (doVirial)
  {
    int nthread = REDUCENONBONDEDVIRIALKERNEL_NUM_WARP*WARPSIZE;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (atomStorageSize-1)/nthread+1);
    reduceNonbondedVirialKernel <<< nblock, nthread, 0, stream >>>
    (doSlow, atomStorageSize, tlKernel.get_xyzq(), d_forces, d_forcesSlow, d_virialEnergy);
    cudaCheck(cudaGetLastError());
  }

  if (doVirial || doEnergy)
  {
    int nthread = REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getTileListVirialEnergyLength()-1)/nthread+1);
    reduceVirialEnergyKernel <<< nblock, nthread, 0, stream >>>
    (doEnergy, doVirial, doSlow, tlKernel.getTileListVirialEnergyLength(), tlKernel.getTileListVirialEnergy(), d_virialEnergy);
    cudaCheck(cudaGetLastError());
  }  

  if (doGBIS && doEnergy)
  {
    int nthread = REDUCEGBISENERGYKERNEL_NUM_WARP*WARPSIZE;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getTileListVirialEnergyGBISLength()-1)/nthread+1);
    reduceGBISEnergyKernel <<< nblock, nthread, 0, stream >>>
    (tlKernel.getTileListVirialEnergyGBISLength(), tlKernel.getTileListVirialEnergy(), d_virialEnergy);
    cudaCheck(cudaGetLastError());
  }

}

void CudaComputeNonbondedKernel::bindExclusions(int numExclusions, unsigned int* exclusion_bits) {
	int nconst = ( numExclusions < MAX_CONST_EXCLUSIONS ? numExclusions : MAX_CONST_EXCLUSIONS );
	cudaCheck(cudaMemcpyToSymbol(constExclusions, exclusion_bits, nconst*sizeof(unsigned int), 0));

  reallocate_device<unsigned int>(&overflowExclusions, &overflowExclusionsSize, numExclusions);
  copy_HtoD_sync<unsigned int>(exclusion_bits, overflowExclusions, numExclusions);
}

void CudaComputeNonbondedKernel::bindVdwCoefTable(float2* h_vdwCoefTable, int vdwCoefTableWidthIn) {
  vdwCoefTableWidth = vdwCoefTableWidthIn;

  int vdwCoefTableArea = vdwCoefTableWidth*vdwCoefTableWidth;

  reallocate_device<float2>(&vdwCoefTable, &vdwCoefTableSize, vdwCoefTableArea);
  copy_HtoD_sync<float2>(h_vdwCoefTable, vdwCoefTable, vdwCoefTableArea);

  if (vdwCoefTableTextureBound) {
    cudaCheck(cudaUnbindTexture(vdwCoefTableTexture));
    vdwCoefTableTextureBound = false;
  }

  vdwCoefTableTexture.normalized = false;
  vdwCoefTableTexture.addressMode[0] = cudaAddressModeClamp;
  vdwCoefTableTexture.filterMode = cudaFilterModePoint;
  cudaCheck(cudaBindTexture((size_t*)0, vdwCoefTableTexture, vdwCoefTable, vdwCoefTableArea*sizeof(float2)));
  vdwCoefTableTextureBound = true;
}

#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
void CudaComputeNonbondedKernel::bindTextureObject(int tableSize, float4* h_table, cudaArray_t& array, cudaTextureObject_t& tableTex) {
  cudaChannelFormatDesc desc;
  desc.x = sizeof(float)*8;
  desc.y = sizeof(float)*8;
  desc.z = sizeof(float)*8;
  desc.w = sizeof(float)*8;
  desc.f = cudaChannelFormatKindFloat;
  cudaCheck(cudaMallocArray(&array, &desc, tableSize, 1));
  cudaCheck(cudaMemcpyToArray(array, 0, 0, h_table, tableSize*sizeof(float4), cudaMemcpyHostToDevice));

  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = array;

  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.filterMode = cudaFilterModeLinear;
  texDesc.normalizedCoords = 1;

  cudaCheck(cudaCreateTextureObject(&tableTex, &resDesc, &texDesc, NULL));
}
#endif

void CudaComputeNonbondedKernel::bindForceAndEnergyTable(int tableSize,
  float4* h_forceTable,
  float4* h_energyTable) {

#ifndef DISABLE_CUDA_TEXTURE_OBJECTS

  if (forceTableTexActive) {
    cudaCheck(cudaFreeArray(forceArray));
    cudaCheck(cudaDestroyTextureObject(forceTableTex));
  }

  if (energyTableTexActive) {
    cudaCheck(cudaFreeArray(energyArray));
    cudaCheck(cudaDestroyTextureObject(energyTableTex));
  }
  forceTableTex = 0;
  energyTableTex = 0;

  bindTextureObject(tableSize, h_forceTable, forceArray, forceTableTex);
  bindTextureObject(tableSize, h_energyTable, energyArray, energyTableTex);
  forceTableTexActive = true;
  energyTableTexActive = true;

#else // DISABLE_CUDA_TEXTURE_OBJECTS

  if (force_table_bound) {
   cudaCheck(cudaFreeArray(forceArray));
  }
  cudaCheck(cudaMallocArray(&forceArray, &force_table.channelDesc, tableSize, 1));

  if (energy_table_bound) {
   cudaCheck(cudaFreeArray(energyArray));
  }
  cudaCheck(cudaMallocArray(&energyArray, &energy_table.channelDesc, tableSize, 1));

  cudaCheck(cudaMemcpyToArray(forceArray, 0, 0, h_forceTable, tableSize*sizeof(float4), cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpyToArray(energyArray, 0, 0, h_energyTable, tableSize*sizeof(float4), cudaMemcpyHostToDevice));

  force_table.normalized = true;
  force_table.addressMode[0] = cudaAddressModeClamp;
  force_table.addressMode[1] = cudaAddressModeClamp;
  force_table.filterMode = cudaFilterModeLinear;

  energy_table.normalized = true;
  energy_table.addressMode[0] = cudaAddressModeClamp;
  energy_table.addressMode[1] = cudaAddressModeClamp;
  energy_table.filterMode = cudaFilterModeLinear;

  cudaCheck(cudaBindTextureToArray(force_table, forceArray));
  cudaCheck(cudaBindTextureToArray(energy_table, energyArray));

  force_table_bound = true;
  energy_table_bound = true;
#endif // DISABLE_CUDA_TEXTURE_OBJECTS
}
