// #include <map>
// #include <algorithm>
#include <cuda.h>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>
#include <cub/cub.cuh>
#include "CudaUtils.h"
#include "CudaTileListKernel.h"
#include "DeviceCUDA.h"
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;

#define OVERALLOC 1.2f

#if __CUDA_ARCH__ < 350
#define __ldg *
#endif

void NAMD_die(const char *);

//
// Calculate the number of lists that contribute to each patch
//
__global__ void calcPatchNumLists(const int numTileLists, const int numPatches,
  const TileList* __restrict__ tileLists, int* __restrict__ patchNumLists) {

  for (int i = threadIdx.x + blockIdx.x*blockDim.x;i < numTileLists;i += blockDim.x*gridDim.x)
  {
    int2 patchInd = tileLists[i].patchInd;
    atomicAdd(&patchNumLists[patchInd.x], 1);
    if (patchInd.x != patchInd.y) atomicAdd(&patchNumLists[patchInd.y], 1);
  }

}

//
// Write patchNumList back to tile list and
// Find empty patches into emptyPatches[0 ... numEmptyPatches-1]
//
__global__ void setPatchNumLists_findEmptyPatches(const int numTileLists,
  TileList* __restrict__ tileLists, const int* __restrict__ patchNumLists,
  const int numPatches, int* __restrict__ numEmptyPatches, int* __restrict__ emptyPatches) {

  for (int i = threadIdx.x + blockIdx.x*blockDim.x;i < numTileLists;i += blockDim.x*gridDim.x)
  {
    int2 patchInd = tileLists[i].patchInd;
    int2 patchNumList = make_int2(patchNumLists[patchInd.x], patchNumLists[patchInd.y]);
    tileLists[i].patchNumList = patchNumList;
  }

  for (int i = threadIdx.x + blockIdx.x*blockDim.x;i < numPatches;i += blockDim.x*gridDim.x)
  {
    if (patchNumLists[i] == 0) {
      int ind = atomicAdd(numEmptyPatches, 1);
      emptyPatches[ind] = i;
    }
  }

}

//
// Builds a sort key that removes zeros but keeps the order otherwise the same
//
__global__ void buildRemoveZerosSortKey(const int numTileLists,
  const unsigned int* __restrict__ tileListDepth, const int begin_bit, unsigned int* __restrict__ sortKey) {

  for (int itileList = threadIdx.x + blockDim.x*blockIdx.x;itileList < numTileLists;itileList += blockDim.x*gridDim.x)
  {
    int depth = (tileListDepth[itileList] >> begin_bit) & 65535;
    sortKey[itileList] = (depth == 0) ? numTileLists : itileList;
  }

}

__global__ void setupSortKey(const int numTileLists, const int maxTileListLen,
  const TileList* __restrict__ tileLists, const unsigned int* __restrict__ tileListDepth,
  const int begin_bit, const unsigned int* __restrict__ sortKeys, unsigned int* __restrict__ sortKey) {

  for (int itileList = threadIdx.x + blockDim.x*blockIdx.x;itileList < numTileLists;itileList += blockDim.x*gridDim.x)
  {
    int icompute = tileLists[itileList].icompute;
    int depth = min((tileListDepth[itileList] >> begin_bit) & 65535, maxTileListLen);
    int i = icompute*maxTileListLen + (depth - 1);
    sortKey[itileList] = (depth == 0) ? 0x7fffffff : sortKeys[i];
  }

}

template <int width>
__global__ void localSort(const int n, const int begin_bit, const int num_bit,
  unsigned int* __restrict__ keys, int* __restrict__ vals) {

  // NOTE: blockDim.x = width

  for (int base = blockDim.x*blockIdx.x;base < n;base += blockDim.x*gridDim.x)
  {
    int i = base + threadIdx.x;
    typedef cub::BlockRadixSort<unsigned int, width, 1, int> BlockRadixSort;
    __shared__ typename BlockRadixSort::TempStorage tempStorage;
    unsigned int key[1] = {(i < n) ? ((keys[i] >> begin_bit) & 65535) : 0};
    int val[1] = {(i < n) ? vals[i] : 0};
    BlockRadixSort(tempStorage).SortDescending(key, val, 0, num_bit);
    if (i < n) {
      keys[i] = key[0];
      vals[i] = val[0];
    }
    __syncthreads();
  }

}

__global__ void storeInReverse(const int numTileListsSrc, const int begin_bit,
  const int* __restrict__ outputOrder, const int* __restrict__ tileListPos,
  const int* __restrict__ tileListOrderSrc,
  const unsigned int* __restrict__ tileListDepthSrc,
  int* __restrict__ tileListOrderDst,
  unsigned int* __restrict__ tileListDepthDst) {

  for (int i = threadIdx.x + blockDim.x*blockIdx.x;i < numTileListsSrc;i += blockDim.x*gridDim.x)
  {
    int j = outputOrder[numTileListsSrc - i - 1];
    if ( ((tileListDepthSrc[j] >> begin_bit) & 65535) > 0 ) {
      int k = tileListPos[i];
      tileListDepthDst[k] = tileListDepthSrc[j];
      tileListOrderDst[k] = j; //tileListOrderSrc[j];
    }
  }
}

//
// Bit shift tileListDepth so that only lower 16 bits are used
//
__global__ void bitshiftTileListDepth(const int numTileLists, const int begin_bit,
  const int* __restrict__ outputOrder, const unsigned int* __restrict__ tileListDepthSrc,
  unsigned int* __restrict__ tileListDepthDst) {

  for (int i = threadIdx.x + blockDim.x*blockIdx.x;i < numTileLists;i+=blockDim.x*gridDim.x)
  {
    int j = outputOrder[numTileLists - i - 1];
    tileListDepthDst[i] = ((tileListDepthSrc[j] >> begin_bit) & 65535) == 0 ? 0 : 1;
  }

}

__global__ void initMinMaxListLen(const int numComputes, const int maxTileListLen,
  int2* __restrict__ minmaxListLen) {

  int2 val;
  val.x = maxTileListLen+1;
  val.y = 0;
  for (int i = threadIdx.x + blockDim.x*blockIdx.x;i < numComputes;i += blockDim.x*gridDim.x)
  {
    minmaxListLen[i] = val;
  }

}

//
// Build sortKeys[], values are in range 0 ... numTileListsDst-1
//
__global__ void buildSortKeys(const int numTileListsDst, const int maxTileListLen,
  const TileList* __restrict__ tileListsSrc,
  const int* __restrict__ tileListOrderDst,
  const unsigned int* __restrict__ tileListDepthDst,
  int2* __restrict__ minmaxListLen, unsigned int* __restrict__ sortKeys) {

  for (int i = threadIdx.x + blockDim.x*blockIdx.x;i < numTileListsDst;i += blockDim.x*gridDim.x)
  {
    int k = tileListOrderDst[i];
    int icompute = tileListsSrc[k].icompute;
    int depth    = tileListDepthDst[i] & 65535;
    // depth is in range [1 ... maxTileListLen]
    int j        = icompute*maxTileListLen + (depth-1);
    sortKeys[j] = i;
    int2 minmax = minmaxListLen[icompute];
    int2 minmaxOrig = minmax;
    if (minmax.x > depth) minmax.x = depth;
    if (minmax.y < depth) minmax.y = depth;
    if (minmax.x != minmaxOrig.x) {
      atomicMin(&minmaxListLen[icompute].x, minmax.x);
    }
    if (minmax.y != minmaxOrig.y) {
      atomicMax(&minmaxListLen[icompute].y, minmax.y);
    }
  }

}

__global__ void fillSortKeys(const int numComputes, const int maxTileListLen,
  const int2* __restrict__ minmaxListLen, unsigned int* __restrict__ sortKeys) {

  int i = (threadIdx.x + blockDim.x*blockIdx.x)/WARPSIZE;
  if (i < numComputes) {
    const int wid = threadIdx.x % WARPSIZE;
    int2 minmax = minmaxListLen[i];
    int minlen = minmax.x;
    int maxlen = minmax.y;
    // minlen, maxlen are in range [1 ... maxTileListLen]
    unsigned int minKey = sortKeys[i*maxTileListLen + minlen-1];
    unsigned int maxKey = sortKeys[i*maxTileListLen + maxlen-1];
    unsigned int aveKey = (maxKey + minKey)/2;
    for (int j=wid;j < minlen-1;j+=WARPSIZE) {
      sortKeys[i*maxTileListLen + j] = minKey;
    }
    for (int j=maxlen+wid;j < maxTileListLen;j+=WARPSIZE) {
      sortKeys[i*maxTileListLen + j] = maxKey;
    }
    for (int j=wid;j < maxTileListLen;j+=WARPSIZE) {
      if (sortKeys[i*maxTileListLen + j] == 0) {
        sortKeys[i*maxTileListLen + j] = aveKey;
      }
    }
  }

}

//
// Calculate bounding boxes for sets of WARPSIZE=32 atoms
//
#define BOUNDINGBOXKERNEL_NUM_WARP 8
__global__ void
buildBoundingBoxesKernel(const int atomStorageSize, const float4* __restrict__ xyzq,
  BoundingBox* __restrict__ boundingBoxes) {

  const int warpId = threadIdx.x / WARPSIZE;
  const int wid = threadIdx.x % WARPSIZE;

  // Loop with warp-aligned index to avoid warp-divergence
  for (int iwarp = warpId*WARPSIZE + blockIdx.x*blockDim.x;iwarp < atomStorageSize;iwarp += blockDim.x*gridDim.x) {
    // Full atom index
    const int i = iwarp + wid;
    // Bounding box index
    const int ibb = i/WARPSIZE;

    float4 xyzq_i = xyzq[min(atomStorageSize-1, i)];

    volatile float3 minxyz, maxxyz;

    typedef cub::WarpReduce<float> WarpReduce;
    __shared__ typename WarpReduce::TempStorage tempStorage[BOUNDINGBOXKERNEL_NUM_WARP];
    minxyz.x = WarpReduce(tempStorage[warpId]).Reduce(xyzq_i.x, cub::Min());
    minxyz.y = WarpReduce(tempStorage[warpId]).Reduce(xyzq_i.y, cub::Min());
    minxyz.z = WarpReduce(tempStorage[warpId]).Reduce(xyzq_i.z, cub::Min());
    maxxyz.x = WarpReduce(tempStorage[warpId]).Reduce(xyzq_i.x, cub::Max());
    maxxyz.y = WarpReduce(tempStorage[warpId]).Reduce(xyzq_i.y, cub::Max());
    maxxyz.z = WarpReduce(tempStorage[warpId]).Reduce(xyzq_i.z, cub::Max());

    if (wid == 0) {
      BoundingBox boundingBox;
      boundingBox.x = 0.5f*(minxyz.x + maxxyz.x);
      boundingBox.y = 0.5f*(minxyz.y + maxxyz.y);
      boundingBox.z = 0.5f*(minxyz.z + maxxyz.z);
      boundingBox.wx = 0.5f*(maxxyz.x - minxyz.x);
      boundingBox.wy = 0.5f*(maxxyz.y - minxyz.y);
      boundingBox.wz = 0.5f*(maxxyz.z - minxyz.z);
      boundingBoxes[ibb] = boundingBox;
    }
  }

}

//
// Returns the lower estimate for the distance between two bounding boxes
//
__device__ __forceinline__ float distsq(const BoundingBox a, const BoundingBox b) {
  float dx = max(0.0f, fabsf(a.x - b.x) - a.wx - b.wx);
  float dy = max(0.0f, fabsf(a.y - b.y) - a.wy - b.wy);
  float dz = max(0.0f, fabsf(a.z - b.z) - a.wz - b.wz);
  float r2 = dx*dx + dy*dy + dz*dz;
  return r2;
}

#if 0
//
// Performs warp-level exclusive sum on a shared memory array:
// sh_out[0 ... n-1] = exclusiveSum( sh_in[0 ... n-1] )
// sh_in and sh_out can point to same array
// Returns the total sum
//
template <typename T>
__device__ __forceinline__
int shWarpExclusiveSum(const int n, volatile T* sh_in, volatile int* sh_out) {
  const int wid = threadIdx.x % WARPSIZE;
  volatile int blockOffset = 0;
  for (int iblock=0;iblock < n;iblock += WARPSIZE) {
    // Size of the exclusive sum
    int blockLen = min(WARPSIZE, n-iblock);
    // Perform exclusive sum on sh_in[iblock ... iblock + blockLen-1]
    typedef cub::WarpScan<int> WarpScan;
    __shared__ typename WarpScan::TempStorage tempStorage;
    int data = (wid < blockLen) ? (int)sh_in[iblock + wid] : 0;
    WarpScan(tempStorage).ExclusiveSum(data, data);
    // Shift by block offset
    data += blockOffset;
    // Save last value
    int last = (int)sh_in[iblock + blockLen-1];
    // Write output
    if (wid < blockLen) sh_out[iblock + wid] = data;
    // block offset = last term of the exclusive sum + last value
    blockOffset = sh_out[iblock + blockLen-1] + last;
  }
  return blockOffset;
}
#endif

#define TILELISTKERNELNEW_NUM_WARP 4

//
// NOTE: Executed on a single thread block
//
template<int nthread>
__global__ void calcTileListPosKernel(const int numComputes,
  const CudaComputeRecord* __restrict__ computes,
  const CudaPatchRecord* __restrict__ patches,
  int* __restrict__ tilePos) {

  typedef cub::BlockScan<int, nthread> BlockScan;

  __shared__ typename BlockScan::TempStorage tempStorage;
  __shared__ int shTilePos0;

  if (threadIdx.x == nthread-1) {
    shTilePos0 = 0;
  }

  for (int base=0;base < numComputes;base+=nthread) {
    int k = base + threadIdx.x;

    int numTiles1 = (k < numComputes) ? (patches[computes[k].patchInd.x].numAtoms-1)/WARPSIZE+1 : 0;

    // Calculate positions in tile list and jtile list
    int tilePosVal;
    BlockScan(tempStorage).ExclusiveSum(numTiles1, tilePosVal);

    // Store into global memory
    if (k < numComputes) {      
      tilePos[k] = shTilePos0 + tilePosVal;
    }

    __syncthreads();
    // Store block end position
    if (threadIdx.x == nthread-1) {
      shTilePos0 += tilePosVal + numTiles1;
    }
  }
}


template<int nthread>
__global__ void updatePatchesKernel(const int numComputes,
  const int* __restrict__ tilePos,
  const CudaComputeRecord* __restrict__ computes,
  const CudaPatchRecord* __restrict__ patches,
  TileList* __restrict__ tileLists) {

  const int tid = threadIdx.x % nthread;

  // nthread threads takes care of one compute
  for (int k = (threadIdx.x + blockIdx.x*blockDim.x)/nthread;k < numComputes;k+=blockDim.x*gridDim.x/nthread)
  {
    CudaComputeRecord compute = computes[k];
    float3 offsetXYZ = compute.offsetXYZ;
    int2 patchInd = compute.patchInd;
    int numTiles1 = (patches[patchInd.x].numAtoms-1)/WARPSIZE+1;
    int itileList0 = tilePos[k];
    for (int i=tid;i < numTiles1;i+=nthread) {
      tileLists[itileList0 + i].offsetXYZ = offsetXYZ;
      tileLists[itileList0 + i].patchInd  = patchInd;
      tileLists[itileList0 + i].icompute  = k;
    }
  }

}

__host__ __device__ __forceinline__
int buildTileListsBBKernel_shmem_sizePerThread(const int maxTileListLen) {
  // Size in bytes
  int size = (
    maxTileListLen*sizeof(char)
    );
  return size;
}

__global__ void
buildTileListsBBKernel(const int numTileLists,
  TileList* __restrict__ tileLists,
  const CudaPatchRecord* __restrict__ patches,
  const int* __restrict__ tileListPos,
  const float latticeX, const float latticeY, const float latticeZ, 
  const float cutoff2, const int maxTileListLen,
  const BoundingBox* __restrict__ boundingBoxes,
  int* __restrict__ tileJatomStart,
  const int tileJatomStartSize,
  unsigned int* __restrict__ tileListDepth,
  int* __restrict__ tileListOrder,
  PatchPairRecord* __restrict__ patchPairs,
  TileListStat* __restrict__ tileListStat) {

  extern __shared__ char sh_buffer[];
  int sizePerThread = buildTileListsBBKernel_shmem_sizePerThread(maxTileListLen);
  int pos = threadIdx.x*sizePerThread;
  volatile char* sh_tile = (char*)&sh_buffer[pos];

  // Loop with warp-aligned index to avoid warp-divergence
  for (int iwarp = (threadIdx.x/WARPSIZE)*WARPSIZE + blockIdx.x*blockDim.x;iwarp < numTileLists;iwarp += blockDim.x*gridDim.x) {

    // Use one thread per tile list
    const int wid = threadIdx.x % WARPSIZE;
    const int itileList = iwarp + wid;

    int i;
    int itileListLen = 0;
    CudaPatchRecord patch1;
    CudaPatchRecord patch2;
    float3 offsetXYZ;
    int2 patchInd;
    int numTiles2;
    int icompute;

    if (itileList < numTileLists) {
      offsetXYZ = tileLists[itileList].offsetXYZ;
      patchInd  = tileLists[itileList].patchInd;
      icompute  = tileLists[itileList].icompute;
      // Get i-column
      i = itileList - tileListPos[icompute];
      float shx = offsetXYZ.x*latticeX;
      float shy = offsetXYZ.y*latticeY;
      float shz = offsetXYZ.z*latticeZ;

      // Load patches
      patch1 = patches[patchInd.x];
      patch2 = patches[patchInd.y];
      // int numTiles1 = (patch1.numAtoms-1)/WARPSIZE+1;
      numTiles2 = (patch2.numAtoms-1)/WARPSIZE+1;
      int tileStart1 = patch1.atomStart/WARPSIZE;
      int tileStart2 = patch2.atomStart/WARPSIZE;
      bool self = (tileStart1 == tileStart2);

      // Load i-atom data (and shift coordinates)
      BoundingBox boundingBoxI = boundingBoxes[i + tileStart1];
      boundingBoxI.x += shx;
      boundingBoxI.y += shy;
      boundingBoxI.z += shz;

      for (int j=0;j < numTiles2;j++) {
        sh_tile[j] = 0;
        if (!self || j >= i) {
          BoundingBox boundingBoxJ = boundingBoxes[j + tileStart2];
          float r2bb = distsq(boundingBoxI, boundingBoxJ);
          if (r2bb < cutoff2) {
            sh_tile[j] = 1;
            itileListLen++;
          }
        }
      }

      tileListDepth[itileList] = (unsigned int)itileListLen;
      tileListOrder[itileList] = itileList;
    }

    typedef cub::WarpScan<int> WarpScan;
    __shared__ typename WarpScan::TempStorage tempStorage;
    int active = (itileListLen > 0);
    int activePos;
    WarpScan(tempStorage).ExclusiveSum(active, activePos);
    int itileListPos;
    WarpScan(tempStorage).ExclusiveSum(itileListLen, itileListPos);

    int jtileStart, numJtiles;
    // Last thread in the warp knows the total number
    if (wid == WARPSIZE-1) {
      atomicAdd(&tileListStat->numTileLists, activePos + active);
      numJtiles = itileListPos + itileListLen;
      jtileStart = atomicAdd(&tileListStat->numJtiles, numJtiles);
    }
    numJtiles  = cub::ShuffleIndex(numJtiles,  WARPSIZE-1);
    jtileStart = cub::ShuffleIndex(jtileStart, WARPSIZE-1);
    if (jtileStart + numJtiles > tileJatomStartSize) {
      // tileJatomStart out of memory, exit 
      if (wid == 0) tileListStat->tilesSizeExceeded = true;
      return;
    }

    int jStart = itileListPos;
    int jEnd   = cub::ShuffleDown(itileListPos, 1);
    if (wid == WARPSIZE-1) jEnd = numJtiles;

    if (itileListLen > 0) {
      // Setup tileLists[]
      TileList TLtmp;
      TLtmp.iatomStart = patch1.atomStart + i*WARPSIZE;
      TLtmp.jtileStart = jtileStart + jStart;
      TLtmp.jtileEnd   = jtileStart + jEnd - 1;
      TLtmp.patchInd   = patchInd;
      TLtmp.offsetXYZ  = offsetXYZ;
      TLtmp.icompute   = icompute;
      // TLtmp.patchNumList.x = 0;
      // TLtmp.patchNumList.y = 0;
      tileLists[itileList] = TLtmp;
      // PatchPair
      PatchPairRecord patchPair;
      patchPair.iatomSize     = patch1.atomStart + patch1.numAtoms;
      patchPair.iatomFreeSize = patch1.atomStart + patch1.numFreeAtoms;
      patchPair.jatomSize     = patch2.atomStart + patch2.numAtoms;
      patchPair.jatomFreeSize = patch2.atomStart + patch2.numFreeAtoms;
      patchPairs[itileList] = patchPair;

      // Write tiles
      int jtile = jtileStart + jStart;
      for (int j=0;j < numTiles2;j++) {
        if (sh_tile[j]) {
          tileJatomStart[jtile] = patch2.atomStart + j*WARPSIZE;
          jtile++;
        }
      }

    }

  }

}

#define REPACKTILELISTSKERNEL_NUM_WARP 32
__global__ void
repackTileListsKernel(const int numTileLists, const int begin_bit, const int* __restrict__ tileListPos,
  const int* __restrict__ tileListOrder,
  const int* __restrict__ jtiles,
  const TileList* __restrict__ tileListsSrc, TileList* __restrict__ tileListsDst,
  const PatchPairRecord* __restrict__ patchPairsSrc, PatchPairRecord* __restrict__ patchPairsDst,
  const int* __restrict__ tileJatomStartSrc, int* __restrict__ tileJatomStartDst,
  const TileExcl* __restrict__ tileExclsSrc, TileExcl* __restrict__ tileExclsDst) {

  const int wid = threadIdx.x % WARPSIZE;

  // One warp does one tile list
  for (int i = (threadIdx.x + blockDim.x*blockIdx.x)/WARPSIZE;i < numTileLists;i+=blockDim.x*gridDim.x/WARPSIZE) 
  {
    int j = tileListOrder[i];
    int start = tileListPos[i];
    int end   = tileListPos[i+1]-1;
    if (wid == 0 && patchPairsSrc != NULL) patchPairsDst[i] = patchPairsSrc[j];
    // TileList
    int startOld   = __ldg(&tileListsSrc[j].jtileStart);
    int endOld     = __ldg(&tileListsSrc[j].jtileEnd);
    int iatomStart = __ldg(&tileListsSrc[j].iatomStart);
    float3 offsetXYZ;
    offsetXYZ.x  = __ldg(&tileListsSrc[j].offsetXYZ.x);
    offsetXYZ.y  = __ldg(&tileListsSrc[j].offsetXYZ.y);
    offsetXYZ.z  = __ldg(&tileListsSrc[j].offsetXYZ.z);
    int2 patchInd = tileListsSrc[j].patchInd;
    int icompute = __ldg(&tileListsSrc[j].icompute);
    if (wid == 0) {
      TileList tileList;
      tileList.iatomStart = iatomStart;
      tileList.offsetXYZ  = offsetXYZ;
      tileList.jtileStart = start;
      tileList.jtileEnd   = end;
      tileList.patchInd   = patchInd;
      tileList.icompute   = icompute;
      tileListsDst[i] = tileList;
    }

    if (jtiles == NULL) {
      // No jtiles, simple copy will do
      int jtile = start;
      for (int jtileOld=startOld;jtileOld <= endOld;jtileOld+=WARPSIZE,jtile+=WARPSIZE) {
        if (jtileOld + wid <= endOld) {
          tileJatomStartDst[jtile + wid] = tileJatomStartSrc[jtileOld + wid];
        }
      }
      if (tileExclsSrc != NULL) {
        int jtile = start;
        for (int jtileOld=startOld;jtileOld <= endOld;jtileOld++,jtile++) {
          tileExclsDst[jtile].excl[wid] = tileExclsSrc[jtileOld].excl[wid];
        }
      }
    } else {
      int jtile0 = start;
      for (int jtileOld=startOld;jtileOld <= endOld;jtileOld+=WARPSIZE) {
        int t = jtileOld + wid;
        int jtile = (t <= endOld) ? jtiles[t] : 0;
        jtile >>= begin_bit;
        jtile &= 65535;
        typedef cub::WarpScan<int> WarpScan;
        __shared__ typename WarpScan::TempStorage tempStorage[REPACKTILELISTSKERNEL_NUM_WARP];
        int warpId = threadIdx.x / WARPSIZE;
        int jtilePos;
        WarpScan(tempStorage[warpId]).ExclusiveSum(jtile, jtilePos);

        if (jtile) tileJatomStartDst[jtile0+jtilePos] = __ldg(&tileJatomStartSrc[t]);

        if (tileExclsSrc != NULL) {
          unsigned int b = __ballot(jtile);
          while (b != 0) {
            // k = index of thread that has data
            int k = __ffs(b) - 1;
            tileExclsDst[jtile0].excl[wid] = __ldg(&tileExclsSrc[jtileOld + k].excl[wid]);
            // remove 1 bit and advance jtile0
            b ^= ((unsigned int)1 << k);
            jtile0++;
          }
        } else {
          jtile0 += __popc(__ballot(jtile));
        }
      }
    }
  }

}

//
// NOTE: Executed on a single thread block
// oobKey = out-of-bounds key value
//
#define SORTTILELISTSKERNEL_NUM_THREAD 512
#define SORTTILELISTSKERNEL_ITEMS_PER_THREAD 22
template <typename keyT, typename valT, bool ascend>
__launch_bounds__ (SORTTILELISTSKERNEL_NUM_THREAD, 1) __global__
void sortTileListsKernel(const int numTileListsSrc, const int numTileListsDst,
  const int begin_bit, const int end_bit, const keyT oobKey,
  keyT* __restrict__ tileListDepthSrc, keyT* __restrict__ tileListDepthDst,
  valT* __restrict__ tileListOrderSrc, valT* __restrict__ tileListOrderDst) {

  typedef cub::BlockLoad<keyT*, SORTTILELISTSKERNEL_NUM_THREAD,
  SORTTILELISTSKERNEL_ITEMS_PER_THREAD, cub::BLOCK_LOAD_WARP_TRANSPOSE> BlockLoadU;

  typedef cub::BlockLoad<valT*, SORTTILELISTSKERNEL_NUM_THREAD,
  SORTTILELISTSKERNEL_ITEMS_PER_THREAD, cub::BLOCK_LOAD_WARP_TRANSPOSE> BlockLoad;

  typedef cub::BlockRadixSort<keyT, SORTTILELISTSKERNEL_NUM_THREAD,
  SORTTILELISTSKERNEL_ITEMS_PER_THREAD, valT> BlockRadixSort;

  __shared__ union {
    typename BlockLoad::TempStorage      load;
    typename BlockLoadU::TempStorage     loadU;
    typename BlockRadixSort::TempStorage sort;
  } tempStorage;

  keyT keys[SORTTILELISTSKERNEL_ITEMS_PER_THREAD];
  valT values[SORTTILELISTSKERNEL_ITEMS_PER_THREAD];

  BlockLoadU(tempStorage.loadU).Load(tileListDepthSrc, keys, numTileListsSrc, oobKey);
  __syncthreads();
  BlockLoad(tempStorage.load).Load(tileListOrderSrc, values, numTileListsSrc);
  __syncthreads();

  if (ascend)
    BlockRadixSort(tempStorage.sort).SortBlockedToStriped(keys, values, begin_bit, end_bit);
  else
    BlockRadixSort(tempStorage.sort).SortDescendingBlockedToStriped(keys, values, begin_bit, end_bit);

  cub::StoreDirectStriped<SORTTILELISTSKERNEL_NUM_THREAD>(threadIdx.x, tileListDepthDst, keys, numTileListsDst);
  cub::StoreDirectStriped<SORTTILELISTSKERNEL_NUM_THREAD>(threadIdx.x, tileListOrderDst, values, numTileListsDst);
}

__global__ void reOrderTileListDepth(const int numTileLists, const int* __restrict__ tileListOrder,
  unsigned int* __restrict__ tileListDepthSrc, unsigned int* __restrict__ tileListDepthDst) {

  for (int i = threadIdx.x + blockDim.x*blockIdx.x;i < numTileLists;i+=blockDim.x*gridDim.x)
  {
    int j = tileListOrder[i];
    tileListDepthDst[i] = tileListDepthSrc[j];
  }

} 

//
// Bit shift tileListDepth so that only lower 16 bits are used
//
__global__ void bitshiftTileListDepth(const int numTileLists, const int begin_bit,
  unsigned int* __restrict__ tileListDepth) {

  for (int itileList = threadIdx.x + blockDim.x*blockIdx.x;itileList < numTileLists;itileList+=blockDim.x*gridDim.x)
  {
    unsigned int a = tileListDepth[itileList];
    a >>= begin_bit;
    a &= 65535;
    tileListDepth[itileList] = a;
  }

}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

CudaTileListKernel::CudaTileListKernel(int deviceID, bool doStreaming) :
deviceID(deviceID), doStreaming(doStreaming) {

  cudaCheck(cudaSetDevice(deviceID));

  activeBuffer = 1;

  numPatches = 0;
  numComputes = 0;

  cudaPatches = NULL;
  cudaPatchesSize = 0;

  cudaComputes = NULL;
  cudaComputesSize = 0;

  patchNumLists = NULL;
  patchNumListsSize = 0;

  emptyPatches = NULL;
  emptyPatchesSize = 0;
  h_emptyPatches = NULL;
  h_emptyPatchesSize = 0;
  numEmptyPatches = 0;

  sortKeySrc = NULL;
  sortKeySrcSize = 0;
  sortKeyDst = NULL;
  sortKeyDstSize = 0;

  tileLists1 = NULL;
  tileLists1Size = 0;
  tileLists2 = NULL;
  tileLists2Size = 0;

  patchPairs1 = NULL;
  patchPairs1Size = 0;
  patchPairs2 = NULL;
  patchPairs2Size = 0;

  tileJatomStart1 = NULL;
  tileJatomStart1Size = 0;
  tileJatomStart2 = NULL;
  tileJatomStart2Size = 0;

  boundingBoxes = NULL;
  boundingBoxesSize = 0;

  tileListDepth1 = NULL;
  tileListDepth1Size = 0;
  tileListDepth2 = NULL;
  tileListDepth2Size = 0;

  tileListOrder1 = NULL;
  tileListOrder1Size = 0;
  tileListOrder2 = NULL;
  tileListOrder2Size = 0;

  tileExcls1 = NULL;
  tileExcls1Size = 0;
  tileExcls2 = NULL;
  tileExcls2Size = 0;

  xyzq = NULL;
  xyzqSize = 0;

  allocate_device<TileListStat>(&d_tileListStat, 1);
  allocate_host<TileListStat>(&h_tileListStat, 1);

  tileListPos = NULL;
  tileListPosSize = 0;
  tempStorage = NULL;
  tempStorageSize = 0;

  jtiles = NULL;
  jtilesSize = 0;

  tilePos = NULL;
  tilePosSize = 0;

  tileListsGBIS = NULL;
  tileListsGBISSize = 0;

  tileJatomStartGBIS = NULL;
  tileJatomStartGBISSize = 0;

  tileListVirialEnergy = NULL;
  tileListVirialEnergySize = 0;

  atomStorageSize = 0;
  numTileLists = 0;
  numTileListsGBIS = 0;
  numJtiles = 1;

  outputOrder = NULL;
  outputOrderSize = 0;
  doOutputOrder = false;

  minmaxListLen = NULL;
  minmaxListLenSize = 0;

  sortKeys = NULL;
  sortKeysSize = 0;
  sortKeys_endbit = 0;

  cudaCheck(cudaEventCreate(&tileListStatEvent));
  tileListStatEventRecord = false;
}

CudaTileListKernel::~CudaTileListKernel() {
  cudaCheck(cudaSetDevice(deviceID));
  deallocate_device<TileListStat>(&d_tileListStat);
  deallocate_host<TileListStat>(&h_tileListStat);
  //
  if (patchNumLists != NULL) deallocate_device<int>(&patchNumLists);
  if (emptyPatches != NULL) deallocate_device<int>(&emptyPatches);
  if (h_emptyPatches != NULL) deallocate_host<int>(&h_emptyPatches);
  if (sortKeySrc != NULL) deallocate_device<unsigned int>(&sortKeySrc);
  if (sortKeyDst != NULL) deallocate_device<unsigned int>(&sortKeyDst);
  //
  if (cudaPatches != NULL) deallocate_device<CudaPatchRecord>(&cudaPatches);
  if (cudaComputes != NULL) deallocate_device<CudaComputeRecord>(&cudaComputes);
  if (patchPairs1 != NULL) deallocate_device<PatchPairRecord>(&patchPairs1);
  if (patchPairs2 != NULL) deallocate_device<PatchPairRecord>(&patchPairs2);
  if (tileLists1 != NULL) deallocate_device<TileList>(&tileLists1);
  if (tileLists2 != NULL) deallocate_device<TileList>(&tileLists2);
  if (tileJatomStart1 != NULL) deallocate_device<int>(&tileJatomStart1);
  if (tileJatomStart2 != NULL) deallocate_device<int>(&tileJatomStart2);
  if (boundingBoxes != NULL) deallocate_device<BoundingBox>(&boundingBoxes);
  if (tileListDepth1 != NULL) deallocate_device<unsigned int>(&tileListDepth1);
  if (tileListDepth2 != NULL) deallocate_device<unsigned int>(&tileListDepth2);
  if (tileListOrder1 != NULL) deallocate_device<int>(&tileListOrder1);
  if (tileListOrder2 != NULL) deallocate_device<int>(&tileListOrder2);
  if (tileListPos != NULL) deallocate_device<int>(&tileListPos);
  if (tileExcls1 != NULL) deallocate_device<TileExcl>(&tileExcls1);
  if (tileExcls2 != NULL) deallocate_device<TileExcl>(&tileExcls2);
  if (tempStorage != NULL) deallocate_device<char>(&tempStorage);
  if (jtiles != NULL) deallocate_device<int>(&jtiles);
  if (tilePos != NULL) deallocate_device<int>(&tilePos);

  if (tileListsGBIS != NULL) deallocate_device<TileList>(&tileListsGBIS);
  if (tileJatomStartGBIS != NULL) deallocate_device<int>(&tileJatomStartGBIS);

  if (tileListVirialEnergy != NULL) deallocate_device<TileListVirialEnergy>(&tileListVirialEnergy);

  if (xyzq != NULL) deallocate_device<float4>(&xyzq);

  if (sortKeys != NULL) deallocate_device<unsigned int>(&sortKeys);
  if (minmaxListLen != NULL) deallocate_device<int2>(&minmaxListLen);

  cudaCheck(cudaEventDestroy(tileListStatEvent));
}

void CudaTileListKernel::prepareTileList(cudaStream_t stream) {
  clear_device_array<int>(jtiles, numJtiles, stream);
}

void CudaTileListKernel::clearTileListStat(cudaStream_t stream) {
  // clear tileListStat, for patchReadyQueueCount, which is set equal to the number of empty patches
  memset(h_tileListStat, 0, sizeof(TileListStat));
  h_tileListStat->patchReadyQueueCount = getNumEmptyPatches();
  copy_HtoD<TileListStat>(h_tileListStat, d_tileListStat, 1, stream);
}

void CudaTileListKernel::finishTileList(cudaStream_t stream) {
  copy_DtoH<TileListStat>(d_tileListStat, h_tileListStat, 1, stream);
  cudaCheck(cudaEventRecord(tileListStatEvent, stream));
  tileListStatEventRecord = true;
}

void CudaTileListKernel::updateComputes(const int numComputesIn,
  const CudaComputeRecord* h_cudaComputes, cudaStream_t stream) {

  numComputes = numComputesIn;

  reallocate_device<CudaComputeRecord>(&cudaComputes, &cudaComputesSize, numComputes);
  copy_HtoD<CudaComputeRecord>(h_cudaComputes, cudaComputes, numComputes, stream);

  if (doStreaming) doOutputOrder = true;
}

void CudaTileListKernel::writeTileList(const char* filename, const int numTileLists,
  const TileList* d_tileLists, cudaStream_t stream) {
  
  TileList* h_tileLists = new TileList[numTileLists];
  copy_DtoH<TileList>(d_tileLists, h_tileLists, numTileLists, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  FILE* handle = fopen(filename,"wt");
  for (int itileList=0;itileList < numTileLists;itileList++) {
    TileList tmp = h_tileLists[itileList];
    fprintf(handle, "%d %d %d %f %f %f %d %d %d %d\n",
      tmp.iatomStart, tmp.jtileStart, tmp.jtileEnd, tmp.offsetXYZ.x, tmp.offsetXYZ.y,
      tmp.offsetXYZ.z, tmp.patchInd.x, tmp.patchInd.y, tmp.patchNumList.x, tmp.patchNumList.y);
  }
  fclose(handle);
  delete [] h_tileLists;
}

void CudaTileListKernel::writeTileJatomStart(const char* filename, const int numJtiles,
  const int* d_tileJatomStart, cudaStream_t stream) {

  int* h_tileJatomStart = new int[numJtiles];
  copy_DtoH<int>(d_tileJatomStart, h_tileJatomStart, numJtiles, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  FILE* handle = fopen(filename,"wt");
  for (int i=0;i < numJtiles;i++) {
    fprintf(handle, "%d\n", h_tileJatomStart[i]);
  }
  fclose(handle);
  delete [] h_tileJatomStart;
}

/*
std::pair<int, int> flip_pair(const std::pair<int, int> &p)
{
    return std::pair<int, int>(p.second, p.first);
}

void CudaTileListKernel::markJtileOverlap(const int width, const int numTileLists, TileList* d_tileLists,
  const int numJtiles, int* d_tileJatomStart, cudaStream_t stream) {

  const int shCacheSize = 10;
  TileList* h_tileLists = new TileList[numTileLists];
  int* h_tileJatomStart = new int[numJtiles];
  copy_DtoH<TileList>(d_tileLists, h_tileLists, numTileLists, stream);
  copy_DtoH<int>(d_tileJatomStart, h_tileJatomStart, numJtiles, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  int ntotal = 0;
  int ncache = 0;
  for (int i=0;i < numTileLists;i+=width) {
    int jend = min(i + width, numTileLists);
    std::map<int, int> atomStartMap;
    std::map<int, int>::iterator it;
    atomStartMap.clear();
    for (int j=i;j < jend;j++) {
      TileList tmp = h_tileLists[j];
      int iatomStart = tmp.iatomStart;
      it = atomStartMap.find(iatomStart);
      if (it == atomStartMap.end()) {
        // Insert new
        atomStartMap.insert( std::pair<int, int>(iatomStart, 0) );
      } else {
        // Increase counter
        it->second--;
      }
      int jtileStart = tmp.jtileStart;
      int jtileEnd   = tmp.jtileEnd;
      ntotal += (jtileEnd - jtileStart + 1) + 1;
      for (int jtile=jtileStart;jtile <= jtileEnd;jtile++) {
        int jatomStart = h_tileJatomStart[jtile];
        it = atomStartMap.find(jatomStart);
        if (it == atomStartMap.end()) {
          // Insert new
          atomStartMap.insert( std::pair<int, int>(jatomStart, 0) );
        } else {
          // Increase counter
          it->second--;
        }
        jatomStart |= (65535 << 16);
        h_tileJatomStart[jtile] = jatomStart;
      }
      iatomStart |= (65535 << 16);
      tmp.iatomStart = iatomStart;
      h_tileLists[j] = tmp;
    }
    ncache += atomStartMap.size();
    std::multimap<int, int> imap;
    imap.clear();
    std::multimap<int, int>::iterator imap_it;
    std::transform(atomStartMap.begin(), atomStartMap.end(), std::inserter(imap, imap.begin()), flip_pair);
    if (i < 400) {
      printf("%d %d\n", ntotal, imap.size());
      for (imap_it = imap.begin();imap_it != imap.end();imap_it++) {
        if (imap_it->first != 0)
          printf("(%d %d)\n", imap_it->first, imap_it->second);
      }
    }
  }
  printf("ntotal %d ncache %d\n", ntotal, ncache);
  copy_HtoD<TileList>(h_tileLists, d_tileLists, numTileLists, stream);
  copy_HtoD<int>(h_tileJatomStart, d_tileJatomStart, numJtiles, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  delete [] h_tileLists;
  delete [] h_tileJatomStart;
}
*/

void CudaTileListKernel::buildTileLists(const int numTileListsPrev,
  const int numPatchesIn, const int atomStorageSizeIn, const int maxTileListLenIn,
  const float latticeX, const float latticeY, const float latticeZ,
  const CudaPatchRecord* h_cudaPatches, const float4* h_xyzq, const float plcutoff2In,
  cudaStream_t stream) {

  numPatches = numPatchesIn;
  atomStorageSize = atomStorageSizeIn;
  maxTileListLen = maxTileListLenIn;
  plcutoff2 = plcutoff2In;

  if (doStreaming) {
    // Re-allocate patchNumLists
    reallocate_device<int>(&patchNumLists, &patchNumListsSize, numPatches);
    reallocate_device<int>(&emptyPatches, &emptyPatchesSize, numPatches+1);
    reallocate_host<int>(&h_emptyPatches, &h_emptyPatchesSize, numPatches+1);
  }

  // Re-allocate (tileLists1, patchPairs1
  reallocate_device<TileList>(&tileLists1, &tileLists1Size, numTileListsPrev, OVERALLOC);
  reallocate_device<PatchPairRecord>(&patchPairs1, &patchPairs1Size, numTileListsPrev, OVERALLOC);

  // Copy cudaPatches to device
  reallocate_device<CudaPatchRecord>(&cudaPatches, &cudaPatchesSize, numPatches);
  copy_HtoD<CudaPatchRecord>(h_cudaPatches, cudaPatches, numPatches, stream);

  // Re-allocate temporary storage
  reallocate_device<int>(&tilePos, &tilePosSize, numComputes, OVERALLOC);
  // Calculate tile list positions (tilePos)
  {
    int nthread = 1024;
    int nblock = 1;
    calcTileListPosKernel<1024> <<< nblock, nthread, 0, stream >>>
    (numComputes, cudaComputes, cudaPatches, tilePos);
    cudaCheck(cudaGetLastError());
  }

  // Build (tileLists1.patchInd, tileLists1.offsetXYZ)
  {
    int nthread = 512;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (numComputes-1)/(nthread/32)+1);
    updatePatchesKernel<32> <<< nblock, nthread, 0, stream >>>
    (numComputes, tilePos, cudaComputes, cudaPatches, tileLists1);
    cudaCheck(cudaGetLastError());
  }

  // ---------------------------------------------------------------------------------------------


  // NOTE: tileListDepth2 and tileListOrder2 must have at least same size as
  // tileListDepth2 and tileListOrder2 since they're used in sorting
  reallocate_device<unsigned int>(&tileListDepth2, &tileListDepth2Size, numTileListsPrev + 1, OVERALLOC);
  reallocate_device<int>(&tileListOrder2, &tileListOrder2Size, numTileListsPrev, OVERALLOC);

  // Allocate with +1 to include last term in the exclusive sum
  reallocate_device<unsigned int>(&tileListDepth1, &tileListDepth1Size, numTileListsPrev + 1, OVERALLOC);

  reallocate_device<int>(&tileListOrder1, &tileListOrder1Size, numTileListsPrev, OVERALLOC);

  reallocate_device<float4>(&xyzq, &xyzqSize, atomStorageSize, OVERALLOC);

  copy_HtoD<float4>(h_xyzq, xyzq, atomStorageSize, stream);

  // Fills in boundingBoxes[0 ... numBoundingBoxes-1]
  {
    int numBoundingBoxes = atomStorageSize/WARPSIZE;
    reallocate_device<BoundingBox>(&boundingBoxes, &boundingBoxesSize, numBoundingBoxes, OVERALLOC);

    int nwarp = BOUNDINGBOXKERNEL_NUM_WARP;
    int nthread = WARPSIZE*nwarp;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (atomStorageSize-1)/nthread+1);
    buildBoundingBoxesKernel <<< nblock, nthread, 0, stream >>> (atomStorageSize, xyzq, boundingBoxes);
    cudaCheck(cudaGetLastError());
  }

  {
    int nwarp = TILELISTKERNELNEW_NUM_WARP;
    int nthread = WARPSIZE*nwarp;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsPrev-1)/nthread+1);

    int shmem_size = buildTileListsBBKernel_shmem_sizePerThread(maxTileListLen)*nthread;

    // NOTE: In the first call numJtiles = 1. buildTileListsBBKernel will return and
    //       tell the required size in h_tileListStat->numJtiles. In subsequent calls,
    //       re-allocation only happens when the size is exceeded.
    h_tileListStat->tilesSizeExceeded = true;
    int reallocCount = 0;
    while (h_tileListStat->tilesSizeExceeded) {
      reallocate_device<int>(&tileJatomStart1, &tileJatomStart1Size, numJtiles, OVERALLOC);

      clearTileListStat(stream);
      // clear_device_array<TileListStat>(d_tileListStat, 1, stream);

      buildTileListsBBKernel <<< nblock, nthread, shmem_size, stream >>>
      (numTileListsPrev, tileLists1, cudaPatches, tilePos,
        latticeX, latticeY, latticeZ,  plcutoff2, maxTileListLen,
        boundingBoxes, tileJatomStart1, tileJatomStart1Size,
        tileListDepth1, tileListOrder1, patchPairs1,
        d_tileListStat);

      cudaCheck(cudaGetLastError());

      // get (numATileLists, numJtiles, tilesSizeExceeded)
      copy_DtoH<TileListStat>(d_tileListStat, h_tileListStat, 1, stream);
      cudaCheck(cudaStreamSynchronize(stream));
      numJtiles = h_tileListStat->numJtiles;

      if (h_tileListStat->tilesSizeExceeded) {
        reallocCount++;
        if (reallocCount > 1) {
          NAMD_die("CudaTileListKernel::buildTileLists, multiple reallocations detected");
        }
      }

    }

    numTileLists = h_tileListStat->numTileLists;

    reallocate_device<int>(&jtiles, &jtilesSize, numJtiles, OVERALLOC);
  }

  // Re-allocate tileListVirialEnergy.
  // NOTE: Since numTileLists here is an upper estimate (since it's based on bounding boxes),
  //       we're quaranteed to have enough space
  reallocate_device<TileListVirialEnergy>(&tileListVirialEnergy, &tileListVirialEnergySize, numTileLists, OVERALLOC);

  reallocate_device<TileList>(&tileLists2, &tileLists2Size, numTileLists, OVERALLOC);
  reallocate_device<PatchPairRecord>(&patchPairs2, &patchPairs2Size, numTileLists, OVERALLOC);
  reallocate_device<int>(&tileJatomStart2, &tileJatomStart2Size, numJtiles, OVERALLOC);
  reallocate_device<TileExcl>(&tileExcls1, &tileExcls1Size, numJtiles, OVERALLOC);
  reallocate_device<TileExcl>(&tileExcls2, &tileExcls2Size, numJtiles, OVERALLOC);

  int numTileListsSrc = numTileListsPrev;
  int numJtilesSrc    = numJtiles;
  int numTileListsDst = numTileLists;
  int numJtilesDst    = numJtiles;

  // Sort tiles
  sortTileLists(
    false,
    0, false,
    numTileListsSrc, numJtilesSrc,
    PtrSize<TileList>(tileLists1, tileLists1Size), PtrSize<int>(tileJatomStart1, tileJatomStart1Size),
    PtrSize<unsigned int>(tileListDepth1, tileListDepth1Size), PtrSize<int>(tileListOrder1, tileListOrder1Size),
    PtrSize<PatchPairRecord>(patchPairs1, patchPairs1Size), PtrSize<TileExcl>(NULL, 0),
    numTileListsDst, numJtilesDst,
    PtrSize<TileList>(tileLists2, tileLists2Size), PtrSize<int>(tileJatomStart2, tileJatomStart2Size),
    PtrSize<unsigned int>(tileListDepth2, tileListDepth2Size), PtrSize<int>(tileListOrder2, tileListOrder2Size),
    PtrSize<PatchPairRecord>(patchPairs2, patchPairs2Size), PtrSize<TileExcl>(NULL, 0),
    stream);

  // Set active buffer to 2
  setActiveBuffer(2);

  if (doOutputOrder) reallocate_device<int>(&outputOrder, &outputOrderSize, numTileLists, OVERALLOC);
}

//
// Returns integer log2(a) rounded up
//
int ilog2(int a) {
  // if (a < 0)
  //   NAMD_die("CudaTileListKernel, ilog2: negative input value not valid");
  int k = 1;
  while (a >>= 1) k++;
  return k;
}

//
// Sort tile lists
//
void CudaTileListKernel::sortTileLists(
  const bool useJtiles,
  const int begin_bit, const bool highDepthBitsSetIn,
  // Source
  const int numTileListsSrc, const int numJtilesSrc,
  PtrSize<TileList> tileListsSrc, PtrSize<int> tileJatomStartSrc,
  PtrSize<unsigned int> tileListDepthSrc, PtrSize<int> tileListOrderSrc,
  PtrSize<PatchPairRecord> patchPairsSrc, PtrSize<TileExcl> tileExclsSrc,
  // Destination
  const int numTileListsDst, const int numJtilesDst,
  PtrSize<TileList> tileListsDst, PtrSize<int> tileJatomStartDst,
  PtrSize<unsigned int> tileListDepthDst, PtrSize<int> tileListOrderDst,
  PtrSize<PatchPairRecord> patchPairsDst, PtrSize<TileExcl> tileExclsDst,
  cudaStream_t stream) {

  bool doShiftDown = (begin_bit != 0 || highDepthBitsSetIn);

  // if (numTileListsDst == 0)
  //   NAMD_die("CudaTileListKernel::sortTileLists, numTileListsDst = 0");

  // Check that the array sizes are adequate
  if (numTileListsSrc > tileListsSrc.size || numJtilesSrc > tileJatomStartSrc.size ||
    numTileListsSrc > tileListDepthSrc.size || numTileListsSrc > tileListOrderSrc.size ||
    (patchPairsSrc.ptr != NULL && numTileListsSrc > patchPairsSrc.size) || 
    (tileExclsSrc.ptr != NULL && numJtilesSrc > tileExclsSrc.size))
    NAMD_die("CudaTileListKernel::sortTileLists, Src allocated too small");

  if (numTileListsDst > tileListsDst.size || numJtilesDst > tileJatomStartDst.size ||
    numTileListsSrc > tileListDepthDst.size || numTileListsSrc > tileListOrderDst.size ||
    (patchPairsDst.ptr != NULL && numTileListsDst > patchPairsDst.size) || 
    (tileExclsDst.ptr != NULL && numJtilesDst > tileExclsDst.size))
    NAMD_die("CudaTileListKernel::sortTileLists, Dst allocated too small");

  if (begin_bit != 0 && begin_bit != 16)
    NAMD_die("CudaTileListKernel::sortTileLists, begin_bit must be 0 or 16");

  // Number of bits needed in the sort
  int num_bit = ilog2(maxTileListLen);
  if (num_bit > 16)
    NAMD_die("CudaTileListKernel::sortTileLists, num_bit overflow");
  int end_bit = begin_bit + num_bit;

  if (doStreaming)
  {
    // ----------------------------------------------------------------------------------------
    if (doOutputOrder && useJtiles) {
      // outputOrder has been produced, put tile lists back in reverse order and produce sortKeys
      // NOTE: This is done very infrequently, typically only once when the MD run starts.

      // Calculate position from depth
      {
        // -----------------------------------------------------------------------------
        // Bit shift & mask tileListDepthDst such that only lower 16 bits are occupied
        // -----------------------------------------------------------------------------
        if (doShiftDown)
        {
          int nthread = 1024;
          int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsSrc-1)/nthread+1);
          bitshiftTileListDepth <<< nblock, nthread, 0, stream >>>
          (numTileListsSrc, begin_bit, outputOrder, tileListDepthSrc.ptr, tileListDepthDst.ptr);
          cudaCheck(cudaGetLastError());
        }

        reallocate_device<int>(&tileListPos, &tileListPosSize, numTileListsSrc, OVERALLOC);

        // --------------------------------------------------------------------
        // Compute itileList positions to store tileLists
        // ExclusiveSum(tileListDepthDst[0...numTileListsDst-1])
        // --------------------------------------------------------------------
        {
          size_t size = 0;
          cudaCheck(cub::DeviceScan::ExclusiveSum(NULL, size,
            (int *)tileListDepthDst.ptr, tileListPos, numTileListsSrc, stream));
          // Make sure tempStorage doesn't remain NULL
          if (size == 0) size = 128;
          reallocate_device<char>(&tempStorage, &tempStorageSize, size, 1.5f);
          size = tempStorageSize;
          cudaCheck(cub::DeviceScan::ExclusiveSum((void *)tempStorage, size,
            (int *)tileListDepthDst.ptr, tileListPos, numTileListsSrc, stream));
        }
      }

      // Store in reverse order from outputOrder
      {
        int nthread = 1024;
        int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsSrc-1)/nthread+1);
        storeInReverse <<< nblock, nthread, 0, stream >>>
        (numTileListsSrc, begin_bit, outputOrder, tileListPos,
          tileListOrderSrc.ptr, tileListDepthSrc.ptr,
          tileListOrderDst.ptr, tileListDepthDst.ptr);
        cudaCheck(cudaGetLastError());
      }

      // Build sortKeys
      {
        maxTileListLen_sortKeys = maxTileListLen;

        reallocate_device<unsigned int>(&sortKeys, &sortKeysSize, numComputes*maxTileListLen);
        clear_device_array<unsigned int>(sortKeys, numComputes*maxTileListLen, stream);

        // Re-allocate and initialize minmaxListLen
        {
          reallocate_device<int2>(&minmaxListLen, &minmaxListLenSize, numComputes);

          int nthread = 1024;
          int nblock = min(deviceCUDA->getMaxNumBlocks(), (numComputes-1)/nthread+1);
          initMinMaxListLen <<< nblock, nthread, 0, stream >>>
          (numComputes, maxTileListLen, minmaxListLen);
          cudaCheck(cudaGetLastError());
        }

        // Build sortKeys and calculate minmaxListLen
        {
          int nthread = 1024;
          int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nthread+1);
          buildSortKeys <<< nblock, nthread, 0, stream >>>
          (numTileListsDst, maxTileListLen, tileListsSrc.ptr, tileListOrderDst.ptr,
            tileListDepthDst.ptr, minmaxListLen, sortKeys);
          cudaCheck(cudaGetLastError());

          // Maximum value in sortKeys[] is numTileListsDst - 1
          sortKeys_endbit = ilog2(numTileListsDst);
        }

        // Fill in missing sortKeys using minmaxListLen
        {
          int nthread = 1024;
          int nwarp = nthread/WARPSIZE;
          int nblock = min(deviceCUDA->getMaxNumBlocks(), (numComputes-1)/nwarp+1);
          fillSortKeys <<< nblock, nthread, 0, stream >>>
          (numComputes, maxTileListLen, minmaxListLen, sortKeys);
          cudaCheck(cudaGetLastError());
        }

      }

      doOutputOrder = false;

    } else if (doOutputOrder) {
      // OutputOrder will be produced in next pairlist non-bond kernel.
      // This time just remove zero length lists
      // NOTE: This is done very infrequently, typically only once when the MD run starts.

      int endbit_tmp = ilog2(numTileListsSrc);

      // Remove zeros
      {
        reallocate_device<unsigned int>(&sortKeySrc, &sortKeySrcSize, numTileListsSrc, OVERALLOC);
        reallocate_device<unsigned int>(&sortKeyDst, &sortKeyDstSize, numTileListsSrc, OVERALLOC);

        int nthread = 1024;
        int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsSrc-1)/nthread+1);
        buildRemoveZerosSortKey <<< nblock, nthread, 0, stream >>> 
        (numTileListsSrc, tileListDepthSrc.ptr, begin_bit, sortKeySrc);
        cudaCheck(cudaGetLastError());
      }

      if (numTileListsSrc <= SORTTILELISTSKERNEL_NUM_THREAD*SORTTILELISTSKERNEL_ITEMS_PER_THREAD)
      {
        // Short list, sort withing a single thread block
        int nthread = SORTTILELISTSKERNEL_NUM_THREAD;
        int nblock = 1;

        unsigned int oobKey = numTileListsSrc;
        sortTileListsKernel <unsigned int, int, true> <<< nblock, nthread, 0, stream >>>
        (numTileListsSrc, numTileListsDst, 0, endbit_tmp, oobKey, sortKeySrc, sortKeyDst,
          tileListOrderSrc.ptr, tileListOrderDst.ptr);
        cudaCheck(cudaGetLastError());
      }
      else
      {
        // Long list, sort on multiple thread blocks
        size_t size = 0;
        cudaCheck(cub::DeviceRadixSort::SortPairs(NULL, size,
          sortKeySrc, sortKeyDst, tileListOrderSrc.ptr, tileListOrderDst.ptr,
          numTileListsSrc, 0, endbit_tmp, stream));
        // Make sure tempStorage doesn't remain NULL
        if (size == 0) size = 128;
        reallocate_device<char>(&tempStorage, &tempStorageSize, size, 1.5f);
        size = tempStorageSize;
        cudaCheck(cub::DeviceRadixSort::SortPairs((void *)tempStorage, size,
          sortKeySrc, sortKeyDst, tileListOrderSrc.ptr, tileListOrderDst.ptr,
          numTileListsSrc, 0, endbit_tmp, stream));
      }

      // Re-order tileListDepth using tileListOrderDst
      {
        int nthread = 1024;
        int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nthread+1);
        reOrderTileListDepth <<< nblock, nthread, 0, stream >>> 
        (numTileListsDst, tileListOrderDst.ptr,
          tileListDepthSrc.ptr, tileListDepthDst.ptr);
        cudaCheck(cudaGetLastError());
      }

    } else {
      // This is done during regular MD cycle

      if (sortKeys_endbit <= 0)
        NAMD_die("CudaTileListKernel::sortTileLists, sortKeys not produced or invalid sortKeys_endbit");

      // Setup sort keys
      {
        reallocate_device<unsigned int>(&sortKeySrc, &sortKeySrcSize, numTileListsSrc, OVERALLOC);
        reallocate_device<unsigned int>(&sortKeyDst, &sortKeyDstSize, numTileListsSrc, OVERALLOC);

        int nthread = 1024;
        int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsSrc-1)/nthread+1);
        setupSortKey <<< nblock, nthread, 0, stream >>> 
        (numTileListsSrc, maxTileListLen_sortKeys, tileListsSrc.ptr, tileListDepthSrc.ptr, begin_bit, sortKeys, sortKeySrc);
        cudaCheck(cudaGetLastError());

      }

      // Global sort
      if (numTileListsSrc <= SORTTILELISTSKERNEL_NUM_THREAD*SORTTILELISTSKERNEL_ITEMS_PER_THREAD)
      // if (false)
      {
        // Short list, sort withing a single thread block
        int nthread = SORTTILELISTSKERNEL_NUM_THREAD;
        int nblock = 1;

        unsigned int oobKey = (2 << sortKeys_endbit) - 1;
        sortTileListsKernel <unsigned int, int, true> <<< nblock, nthread, 0, stream >>>
        (numTileListsSrc, numTileListsDst, 0, sortKeys_endbit, oobKey, sortKeySrc, sortKeyDst,
          tileListOrderSrc.ptr, tileListOrderDst.ptr);
        cudaCheck(cudaGetLastError());
      }
      else
      {
        // Long list, sort on multiple thread blocks
        size_t size = 0;
        cudaCheck(cub::DeviceRadixSort::SortPairs(NULL, size,
          sortKeySrc, sortKeyDst, tileListOrderSrc.ptr, tileListOrderDst.ptr,
          numTileListsSrc, 0, sortKeys_endbit, stream));
        // Make sure tempStorage doesn't remain NULL
        if (size == 0) size = 128;
        reallocate_device<char>(&tempStorage, &tempStorageSize, size, 1.5f);
        size = tempStorageSize;
        cudaCheck(cub::DeviceRadixSort::SortPairs((void *)tempStorage, size,
          sortKeySrc, sortKeyDst, tileListOrderSrc.ptr, tileListOrderDst.ptr,
          numTileListsSrc, 0, sortKeys_endbit, stream));
      }

      // Re-order tileListDepth using tileListOrderDst
      {
        int nthread = 1024;
        int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nthread+1);
        reOrderTileListDepth <<< nblock, nthread, 0, stream >>> 
        (numTileListsDst, tileListOrderDst.ptr,
          tileListDepthSrc.ptr, tileListDepthDst.ptr);
        cudaCheck(cudaGetLastError());
      }

      // Local sort
      {
        int nthread = 32;
        int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nthread+1);
        localSort<32> <<< nblock, nthread, 0, stream >>> 
        (numTileListsDst, begin_bit, num_bit, tileListDepthDst.ptr, tileListOrderDst.ptr);
        cudaCheck(cudaGetLastError());

        // No need to shift any more
        doShiftDown = false;
      }

    }
    // ----------------------------------------------------------------------------------------

  } // (doStreaming)
  else
  {
    // --------------------------------------------------------------------
    // Sort {tileListDepthSrc, tileListOrderSrc}[0 ... numTileListsSrc-1]
    //   => {tileListDepthDst, tileListOrderDst}[0 ... numTileListsSrc-1]
    // --------------------------------------------------------------------
    if (numTileListsSrc <= SORTTILELISTSKERNEL_NUM_THREAD*SORTTILELISTSKERNEL_ITEMS_PER_THREAD)
    {
      // Short list, sort withing a single thread block
      int nthread = SORTTILELISTSKERNEL_NUM_THREAD;
      int nblock = 1;

      sortTileListsKernel<unsigned int, int, false> <<< nblock, nthread, 0, stream >>>
      (numTileListsSrc, numTileListsDst, begin_bit, end_bit, 0, tileListDepthSrc.ptr, tileListDepthDst.ptr,
        tileListOrderSrc.ptr, tileListOrderDst.ptr);
      cudaCheck(cudaGetLastError());

    }
    else
    {
      // Long list, sort on multiple thread blocks
      size_t size = 0;
      cudaCheck(cub::DeviceRadixSort::SortPairsDescending(NULL, size,
        tileListDepthSrc.ptr, tileListDepthDst.ptr, tileListOrderSrc.ptr, tileListOrderDst.ptr,
        numTileListsSrc, begin_bit, end_bit, stream));
      // Make sure tempStorage doesn't remain NULL
      if (size == 0) size = 128;
      reallocate_device<char>(&tempStorage, &tempStorageSize, size, 1.5f);
      size = tempStorageSize;
      cudaCheck(cub::DeviceRadixSort::SortPairsDescending((void *)tempStorage, size,
        tileListDepthSrc.ptr, tileListDepthDst.ptr, tileListOrderSrc.ptr, tileListOrderDst.ptr,
        numTileListsSrc, begin_bit, end_bit, stream));
    }
  }

  // -----------------------------------------------------------------------------
  // Bit shift & mask tileListDepthDst such that only lower 16 bits are occupied
  // -----------------------------------------------------------------------------
  if (doShiftDown)
  {
    int nthread = 1024;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nthread+1);
    bitshiftTileListDepth <<< nblock, nthread, 0, stream >>>
    (numTileListsDst, begin_bit, tileListDepthDst.ptr);
    cudaCheck(cudaGetLastError());
  }

  // Allocate with +1 to include last term in the exclusive sum
  reallocate_device<int>(&tileListPos, &tileListPosSize, numTileListsDst+1, OVERALLOC);

  // --------------------------------------------------------------------
  // Compute itileList positions to store tileLists
  // ExclusiveSum(tileListDepthDst[0...numTileListsDst+1])
  // NOTE: tileListDepthDst[numTileListsDst] is not accessed
  //       since this is an exclusive sum. But with this trick,
  //       tileListPos[numTileListsDst] will contain the total number
  //       of tile lists
  // --------------------------------------------------------------------
  {
    size_t size = 0;
    cudaCheck(cub::DeviceScan::ExclusiveSum(NULL, size,
      (int *)tileListDepthDst.ptr, tileListPos, numTileListsDst+1, stream));
    // Make sure tempStorage doesn't remain NULL
    if (size == 0) size = 128;
    reallocate_device<char>(&tempStorage, &tempStorageSize, size, 1.5f);
    size = tempStorageSize;
    // NOTE: Bug in CUB 1.4.1, stalls here with Geforce GTC Titan X.
    //       Tested on "manila" node at UIUC. Works OK with CUB 1.5.2
    cudaCheck(cub::DeviceScan::ExclusiveSum((void *)tempStorage, size,
      (int *)tileListDepthDst.ptr, tileListPos, numTileListsDst+1, stream));
  }

  // --------------------------------------------------------------------
  // Re-package to
  // tileListsDst[0 ... numTileListsDst-1], tileJatomStartDst[0 ... numJtilesDst-1]
  // patchPairDst[]
  // tileJatomStartDst[]
  // tileExclsDst[0 ... numJtilesDst-1]
  // --------------------------------------------------------------------
  {
    int nthread = WARPSIZE*REPACKTILELISTSKERNEL_NUM_WARP;
    int nwarp = REPACKTILELISTSKERNEL_NUM_WARP;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nwarp+1);

    repackTileListsKernel <<< nblock, nthread, 0, stream >>>
    (numTileListsDst, begin_bit, tileListPos, tileListOrderDst.ptr,
      (useJtiles) ? jtiles : NULL,
      tileListsSrc.ptr, tileListsDst.ptr,
      patchPairsSrc.ptr, patchPairsDst.ptr,
      tileJatomStartSrc.ptr, tileJatomStartDst.ptr,
      tileExclsSrc.ptr, tileExclsDst.ptr);
    cudaCheck(cudaGetLastError());
  }

  // Count the number of tileLists that contribute to each patch
  if (doStreaming)
  {
    clear_device_array<int>(patchNumLists, numPatches, stream);

    // Fill in patchNumLists[0...numPatches-1]
    int nthread = 512;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (numTileListsDst-1)/nthread+1);
    calcPatchNumLists <<< nblock, nthread, 0, stream >>> 
    (numTileListsDst, numPatches, tileListsDst.ptr, patchNumLists);
    cudaCheck(cudaGetLastError());

    // Use emptyPatches[numPatches] as the count variable
    clear_device_array<int>(&emptyPatches[numPatches], 1, stream);

    // Fill in tileListsDst[0...numTileListsDst-1].patchNumLists
    // and find empty patches into emptyPatches[0 ... numEmptyPatches - 1]
    setPatchNumLists_findEmptyPatches <<< nblock, nthread, 0, stream >>> 
    (numTileListsDst, tileListsDst.ptr, patchNumLists,
      numPatches, &emptyPatches[numPatches], emptyPatches);
    cudaCheck(cudaGetLastError());

    // // Copy emptyPatches[0 ... numPatches] to host
    copy_DtoH<int>(emptyPatches, h_emptyPatches, numPatches+1, stream);
    cudaCheck(cudaStreamSynchronize(stream));
    numEmptyPatches = h_emptyPatches[numPatches];
  }

}

//
// Re-sort tile lists after pairlist refinement. Can only be called after finishTileList() has finished copying
//
void CudaTileListKernel::reSortTileLists(const bool doGBIS, cudaStream_t stream) {
  // Store previous number of active lists
  int numTileListsPrev = numTileLists;

  // Wait for finishTileList() to stop copying
  if (!tileListStatEventRecord)
    NAMD_die("CudaTileListKernel::reSortTileLists, tileListStatEvent not recorded");
  cudaCheck(cudaEventSynchronize(tileListStatEvent));

  // Get numTileLists, numTileListsGBIS, and numExcluded
  {
    numTileLists     = h_tileListStat->numTileLists;
    numTileListsGBIS = h_tileListStat->numTileListsGBIS;
    numExcluded      = h_tileListStat->numExcluded;
  }

  // Sort {tileLists2, tileJatomStart2, tileExcl2} => {tileLists1, tileJatomStart1, tileExcl1}
  // VdW tile list in {tileLists1, tileJatomStart1, tileExcl1}
  sortTileLists(true, 0, true,
    numTileListsPrev, numJtiles,
    PtrSize<TileList>(tileLists2, tileLists2Size), PtrSize<int>(tileJatomStart2, tileJatomStart2Size),
    PtrSize<unsigned int>(tileListDepth2, tileListDepth2Size), PtrSize<int>(tileListOrder2, tileListOrder2Size),
    PtrSize<PatchPairRecord>(NULL, 0), PtrSize<TileExcl>(tileExcls2, tileExcls2Size),
    numTileLists, numJtiles,
    PtrSize<TileList>(tileLists1, tileLists1Size), PtrSize<int>(tileJatomStart1, tileJatomStart1Size),
    PtrSize<unsigned int>(tileListDepth1, tileListDepth1Size), PtrSize<int>(tileListOrder1, tileListOrder1Size),
    PtrSize<PatchPairRecord>(NULL, 0), PtrSize<TileExcl>(tileExcls1, tileExcls1Size),
    stream);

  // fprintf(stderr, "reSortTileLists, writing tile lists to disk...\n");
  // writeTileList("tileList.txt", numTileLists, tileLists1, stream);
  // writeTileJatomStart("tileJatomStart.txt", numJtiles, tileJatomStart1, stream);

  // markJtileOverlap(4, numTileLists, tileLists1, numJtiles, tileJatomStart1, stream);

  // NOTE:
  // Only {tileList1, tileJatomStart1, tileExcl1} are used from here on,
  // the rest {tileListDepth1, tileListOrder1, patchPairs1} may be re-used by the GBIS sorting

  if (doGBIS) {
    // GBIS is used => produce a second tile list
    // GBIS tile list in {tileListGBIS, tileJatomStartGBIS, patchPairs1}
    reallocate_device<TileList>(&tileListsGBIS, &tileListsGBISSize, numTileListsGBIS, OVERALLOC);
    reallocate_device<int>(&tileJatomStartGBIS, &tileJatomStartGBISSize, numJtiles, OVERALLOC);

    sortTileLists(true, 16, true,
      numTileListsPrev, numJtiles,
      PtrSize<TileList>(tileLists2, tileLists2Size), PtrSize<int>(tileJatomStart2, tileJatomStart2Size),
      PtrSize<unsigned int>(tileListDepth2, tileListDepth2Size), PtrSize<int>(tileListOrder2, tileListOrder2Size),
      PtrSize<PatchPairRecord>(patchPairs2, patchPairs2Size), PtrSize<TileExcl>(NULL, 0),
      numTileListsGBIS, numJtiles,
      PtrSize<TileList>(tileListsGBIS, tileListsGBISSize), PtrSize<int>(tileJatomStartGBIS, tileJatomStartGBISSize),
      PtrSize<unsigned int>(tileListDepth1, tileListDepth1Size), PtrSize<int>(tileListOrder1, tileListOrder1Size),
      PtrSize<PatchPairRecord>(patchPairs1, patchPairs1Size), PtrSize<TileExcl>(NULL, 0),
      stream);
  }

  // Set active buffer to be 1
  setActiveBuffer(1);

}

/*
//
// Apply outputOrder after regular (non-pairlist, non-energy) non-bonded kernel
//
void CudaTileListKernel::applyOutputOrder(cudaStream_t stream) {
  return;

  if (!doStreaming || !doOutputOrder)
    return;

  // Sort {tileList1, tileJatomStart1, tileExcl1} => {tileList2, tileJatomStart2, tileExcl2}
  // VdW tile list in {tileList1, tileJatomStart1, tileExcl1}
  sortTileLists(false, 0, true,
    numTileLists, numJtiles,
    PtrSize<TileList>(tileLists1, tileLists1Size), PtrSize<int>(tileJatomStart1, tileJatomStart1Size),
    PtrSize<unsigned int>(tileListDepth1, tileListDepth1Size), PtrSize<int>(tileListOrder1, tileListOrder1Size),
    PtrSize<PatchPairRecord>(NULL, 0), PtrSize<TileExcl>(tileExcls1, tileExcls2Size),
    numTileLists, numJtiles,
    PtrSize<TileList>(tileLists2, tileLists2Size), PtrSize<int>(tileJatomStart2, tileJatomStart2Size),
    PtrSize<unsigned int>(tileListDepth2, tileListDepth2Size), PtrSize<int>(tileListOrder2, tileListOrder2Size),
    PtrSize<PatchPairRecord>(NULL, 0), PtrSize<TileExcl>(tileExcls2, tileExcls2Size),
    stream);

  // Set active buffer to be 2
  setActiveBuffer(2);

}
*/

void CudaTileListKernel::setTileListVirialEnergyLength(int len) {
  if (len > tileListVirialEnergySize) {
    NAMD_die("CudaTileListKernel::setTileListVirialEnergyLength, size overflow");
  }
  tileListVirialEnergyLength = len;
}

void CudaTileListKernel::setTileListVirialEnergyGBISLength(int len) {
  if (len > tileListVirialEnergySize) {
    NAMD_die("CudaTileListKernel::setTileListVirialEnergyGBISLength, size overflow");
  }
  tileListVirialEnergyGBISLength = len;
}
