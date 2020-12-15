#ifndef CUDATILELISTKERNEL_H
#define CUDATILELISTKERNEL_H
#ifdef NAMD_CUDA

// Exclusion mask: bit 1 = atom pair is included, 0 = atom pair is excluded
struct TileExcl {
  unsigned int excl[32];
};

struct TileList {
  int iatomStart;
  int jtileStart;
  int jtileEnd;
  float3 offsetXYZ;
  int2 patchInd;        // Patch indices for this list
  union {  
    int2 patchNumList;    // Number of lists contributing to each patch
    // int icompute;
  };
  int icompute;
};

struct PatchPairRecord {
  int iatomSize;
  int iatomFreeSize;
  int jatomSize;
  int jatomFreeSize;
};

//
// Bounding box structure
//
struct BoundingBox {
  float x, y, z;      // Center
  float wx, wy, wz;   // Half-width
};

//
// Stripped-down CUDA version of compute record
//
struct CudaComputeRecord {
  int2 patchInd;
  float3 offsetXYZ;
};

//
// Stripped-down CUDA version of patch record
//
struct CudaPatchRecord {
  int numAtoms;
  int numFreeAtoms;
  int atomStart;
};

//
// Tile list status. Used to communicate tile list sizes between GPU and CPU
//
struct TileListStat {
  int numTileLists;
  int numTileListsGBIS;
  int numJtiles;
  int numExcluded;
  int patchReadyQueueCount;
  int outputOrderIndex;
  bool tilesSizeExceeded;
};

struct TileListVirialEnergy {
  float shx, shy, shz;
  float forcex, forcey, forcez;
  float forceSlowx, forceSlowy, forceSlowz;
  double energyVdw;
  double energyElec;
  double energySlow;
  double energyGBIS;
};

struct VirialEnergy {
  double virial[9];
  double virialSlow[9];
  double energyVdw;
  double energyElec;
  double energySlow;
  double energyGBIS;
};

class CudaTileListKernel {
private:

  template <typename T>
  struct PtrSize {
    PtrSize(T* ptr, int size) : ptr(ptr), size(size) {}
    T* ptr;
    int size;
  };

  const int deviceID;

  // Events
  cudaEvent_t tileListStatEvent;
  bool tileListStatEventRecord;

  // Pair list cutoff squared
  float plcutoff2;

  // Number of patches
  int numPatches;

  // Number of computes
  int numComputes;

  // Number of tile lists
  int numTileLists;

  // Number of tile lists for GBIS
  int numTileListsGBIS;

  // Number of tiles
  int numJtiles;

  // Maximum number of tiles per tile list
  int maxTileListLen;

  CudaPatchRecord* cudaPatches;
  int cudaPatchesSize;

  CudaComputeRecord* cudaComputes;
  int cudaComputesSize;

  // --- For Streaming ---
  const bool doStreaming;
  int* patchNumLists;
  int patchNumListsSize;

  int* emptyPatches;
  int emptyPatchesSize;
  int* h_emptyPatches;
  int h_emptyPatchesSize;
  int numEmptyPatches;

  unsigned int* sortKeySrc;
  int sortKeySrcSize;
  unsigned int* sortKeyDst;
  int sortKeyDstSize;

  int maxTileListLen_sortKeys;
  
  unsigned int* sortKeys;
  int sortKeysSize;

  int2* minmaxListLen;
  int minmaxListLenSize;

  int sortKeys_endbit;
  // ---------------------

  // Single entry pinned host and device buffers for communicating tile list status
  TileListStat* h_tileListStat;
  TileListStat* d_tileListStat;

  // Atom coordinates and charge
  float4* xyzq;
  int xyzqSize;
  // Atom coordinate storage size
  int atomStorageSize;

  // Tile lists
  TileList* tileLists1;
  int tileLists1Size;
  TileList* tileLists2;
  int tileLists2Size;
  TileList* tileListsGBIS;
  int tileListsGBISSize;

  // Pair pairs
  PatchPairRecord* patchPairs1;
  int patchPairs1Size;
  PatchPairRecord* patchPairs2;
  int patchPairs2Size;

  // j-atom start for tiles
  int* tileJatomStart1;
  int tileJatomStart1Size;
  int* tileJatomStart2;
  int tileJatomStart2Size;
  int* tileJatomStartGBIS;
  int tileJatomStartGBISSize;

  // Bounding boxes
  BoundingBox* boundingBoxes;
  int boundingBoxesSize;

  // Depth of each tile list
  unsigned int* tileListDepth1;
  int tileListDepth1Size;
  unsigned int* tileListDepth2;
  int tileListDepth2Size;

  // Tile list order
  int* tileListOrder1;
  int tileListOrder1Size;
  int* tileListOrder2;
  int tileListOrder2Size;

  // Position of each tile list = ExclusiveSum(tileListDepths)
  int* tileListPos;
  int tileListPosSize;

  // jtile occupancy and position
  int* jtiles;
  int jtilesSize;

  // Temporary buffers used in buildTileLists
  int* tilePos;
  int tilePosSize;

  // Exclusions
  TileExcl* tileExcls1;
  int tileExcls1Size;
  TileExcl* tileExcls2;
  int tileExcls2Size;

  // Temporary storage for CUB
  char* tempStorage;
  int tempStorageSize;

  // Number of exclusions detected
  int numExcluded;

  // Virials and energies for tile lists
  TileListVirialEnergy* tileListVirialEnergy;
  int tileListVirialEnergySize;

  int tileListVirialEnergyLength;
  int tileListVirialEnergyGBISLength;

  int activeBuffer;

  void setActiveBuffer(int activeBufferIn) {activeBuffer = activeBufferIn;}

  void sortTileLists(
    const bool useJtiles,
    const int begin_bit, const bool highDepthBitsSet,
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
    cudaStream_t stream);

  void writeTileList(const char* filename, const int numTileLists,
    const TileList* d_tileLists, cudaStream_t stream);
  void writeTileJatomStart(const char* filename, const int numJtiles,
    const int* d_tileJatomStart, cudaStream_t stream);
  // void markJtileOverlap(const int width, const int numTileLists, TileList* d_tileLists,
  //   const int numJtiles, int* d_tileJatomStart, cudaStream_t stream);

  int* outputOrder;
  int outputOrderSize;
  bool doOutputOrder;

public:

	CudaTileListKernel(int deviceID, bool doStreaming);
	~CudaTileListKernel();

  int getNumEmptyPatches() {return numEmptyPatches;}
  int* getEmptyPatches() {return h_emptyPatches;}

  int getNumExcluded() {return numExcluded;}

  float get_plcutoff2() {return plcutoff2;}
  int getNumTileLists() {return numTileLists;}
  int getNumTileListsGBIS() {return numTileListsGBIS;}
  int getNumJtiles() {return numJtiles;}
  BoundingBox* getBoundingBoxes() {return boundingBoxes;}
  int* getJtiles() {return jtiles;}
	float4* get_xyzq() {return xyzq;}

  TileListStat* getTileListStatDevPtr() {return d_tileListStat;}
  void clearTileListStat(cudaStream_t stream);

  int* getTileJatomStart() {return ((activeBuffer == 1) ? tileJatomStart1 : tileJatomStart2);}
  TileList* getTileLists() {
    return ((activeBuffer == 1) ? tileLists1 : tileLists2);
  }
  unsigned int* getTileListDepth() {return ((activeBuffer == 1) ? tileListDepth1 : tileListDepth2);}
  int* getTileListOrder() {return ((activeBuffer == 1) ? tileListOrder1 : tileListOrder2);}
  TileExcl* getTileExcls() {return ((activeBuffer == 1) ? tileExcls1 : tileExcls2);}
  PatchPairRecord* getPatchPairs() {return ((activeBuffer == 1) ? patchPairs1 : patchPairs2);}

  int* getTileJatomStartGBIS() {return tileJatomStartGBIS;}
  TileList* getTileListsGBIS() {return tileListsGBIS;}

  TileListVirialEnergy* getTileListVirialEnergy() {return tileListVirialEnergy;}

  CudaPatchRecord* getCudaPatches() {return cudaPatches;}

  void prepareTileList(cudaStream_t stream);
	void finishTileList(cudaStream_t stream);

  void updateComputes(const int numComputesIn,
    const CudaComputeRecord* h_cudaComputes, cudaStream_t stream);

  void buildTileLists(const int numTileListsPrev,
    const int numPatchesIn, const int atomStorageSizeIn, const int maxTileListLenIn,
    const float latticeX, const float latticeY, const float latticeZ,
    const CudaPatchRecord* h_cudaPatches, const float4* h_xyzq, const float plcutoff2In, cudaStream_t stream);

  void reSortTileLists(const bool doGBIS, cudaStream_t stream);
  // void applyOutputOrder(cudaStream_t stream);

  void setTileListVirialEnergyLength(int len);
  void setTileListVirialEnergyGBISLength(int len);
  int getTileListVirialEnergyLength() {return tileListVirialEnergyLength;}
  int getTileListVirialEnergyGBISLength() {return tileListVirialEnergyGBISLength;}

  int getNumPatches() {return numPatches;}

  int getNumComputes() {return numComputes;}
  int* getOutputOrder() {
    if (!doStreaming) return NULL;
    if (doOutputOrder) {
      return outputOrder;
    } else {
      return NULL;
    }
  }

};
#endif // NAMD_CUDA
#endif // CUDATILELISTKERNEL_H
