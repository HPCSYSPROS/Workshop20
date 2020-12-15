#ifndef PMESOLVERUTIL_H
#define PMESOLVERUTIL_H

#include "PmeBase.h"   // PmeGrid -structure
#include "NamdTypes.h" // CudaAtom, float2
#include "Lattice.h"

//
// Stores the patch origin (x, y, z) and number of atoms (w)
//
struct PatchInfo {
  int x, y, z, w;
  PatchInfo() {}
  PatchInfo(int x, int y, int z, int w) : x(x), y(y), z(z), w(w) {}
};

enum {Perm_X_Y_Z, Perm_cX_Y_Z, Perm_Y_Z_cX, Perm_Z_cX_Y};

// y = absolute grid position [0 ... pmeGrid.K2-1]
static inline int getPencilIndexY(const PmeGrid& pmeGrid, const int y) {
  return (y*pmeGrid.yBlocks + pmeGrid.yBlocks - 1)/pmeGrid.K2;
}

// z = absolute grid position [0 ... pmeGrid.K3-1]
static inline int getPencilIndexZ(const PmeGrid& pmeGrid, const int z) {
  return (z*pmeGrid.zBlocks + pmeGrid.zBlocks - 1)/pmeGrid.K3;
}

static void getPencilDim(const PmeGrid& pmeGrid, const int permutation,
  const int jblock, const int kblock,
  int& i0, int& i1, int& j0, int& j1, int& k0, int& k1) {

  int isize, jsize, ksize;
  int jblocks, kblocks;

  switch(permutation) {
    case Perm_X_Y_Z:
    isize = pmeGrid.K1;
    jsize = pmeGrid.K2;
    ksize = pmeGrid.K3;
    jblocks = pmeGrid.yBlocks;
    kblocks = pmeGrid.zBlocks;
    break;
    case Perm_cX_Y_Z:
    isize = pmeGrid.K1/2+1;
    jsize = pmeGrid.K2;
    ksize = pmeGrid.K3;
    jblocks = pmeGrid.yBlocks;
    kblocks = pmeGrid.zBlocks;
    break;
    case Perm_Y_Z_cX:
    isize = pmeGrid.K2;
    jsize = pmeGrid.K3;
    ksize = pmeGrid.K1/2+1;
    jblocks = pmeGrid.zBlocks;
    kblocks = pmeGrid.xBlocks;
    break;
    case Perm_Z_cX_Y:
    isize = pmeGrid.K3;
    jsize = pmeGrid.K1/2+1;
    ksize = pmeGrid.K2;
    jblocks = pmeGrid.xBlocks;
    kblocks = pmeGrid.yBlocks;
    break;
    default:
    NAMD_bug("getPencilDim, invalid permutation");
    break;
  }

  if (jblock < 0 || jblock >= jblocks || kblock < 0 || kblock >= kblocks)
    NAMD_bug("getPencilDim, invalid block indices");

  i0 = 0;
  i1 = isize - 1;

  j0 = jsize*jblock/jblocks;
  j1 = jsize*(jblock+1)/jblocks - 1;

  k0 = ksize*kblock/kblocks;
  k1 = ksize*(kblock+1)/kblocks - 1;
}

//
// Return block dimensions [i0...i1] x [j0...j1] x [k0...k1]
//
static void getBlockDim(const PmeGrid& pmeGrid, const int permutation,
  const int iblock, const int jblock, const int kblock,
  int& i0, int& i1, int& j0, int& j1, int& k0, int& k1) {

  getPencilDim(pmeGrid, permutation, jblock, kblock, i0, i1, j0, j1, k0, k1);

  int iblocks;

  switch(permutation) {
    case Perm_X_Y_Z:
    iblocks = pmeGrid.xBlocks;
    break;
    case Perm_cX_Y_Z:
    iblocks = pmeGrid.xBlocks;
    break;
    case Perm_Y_Z_cX:
    iblocks = pmeGrid.yBlocks;
    break;
    case Perm_Z_cX_Y:
    iblocks = pmeGrid.zBlocks;
    break;
    default:
    NAMD_bug("getBlockDim, invalid permutation");
    break;
  }

  if (iblock < 0 || iblock >= iblocks)
    NAMD_bug("getBlockDim, invalid block index");

  int isize = i1-i0+1;

  i0 = isize*iblock/iblocks;
  i1 = isize*(iblock+1)/iblocks - 1;
}

//
// Abstract FFTCompute base class for out-of-place FFTs
// NOTE: Must call init -method to initialize this class
//
class FFTCompute {
public:
  FFTCompute() {
    dataSrc = NULL;
    dataSrcSize = 0;
    dataSrcAllocated = false;
    dataDst = NULL;
    dataDstSize = 0;
    dataDstAllocated = false;
  }
  void init(
    float* dataSrc_in, int dataSrcSize_in,
    float* dataDst_in, int dataDstSize_in,
    int permutation, PmeGrid pmeGrid, 
    int pmePencilType, int jblock, int kblock, int flags) {

    if (dataSrc_in != NULL && dataSrc_in == dataDst_in)
      NAMD_bug("FFTCompute::init, only out-of-place FFTs supported");

    int permutationDst = permutation;
    if (permutation == Perm_X_Y_Z) {
      permutationDst = Perm_cX_Y_Z;
    }

    if (dataSrc_in == NULL) {
      // Sets data and dataSize
      dataSrcSize = getDataSizeRequired(permutation, pmeGrid, jblock, kblock);
      dataSrc = allocateData(dataSrcSize);
      dataSrcAllocated = true;
    } else {
      if (dataSrcSize_in < getDataSizeRequired(permutation, pmeGrid, jblock, kblock))
        NAMD_bug("FFTCompute::init, invalid dataSrcSize_in");
      dataSrcSize = dataSrcSize_in;
      dataSrc = dataSrc_in;
      dataSrcAllocated = false;
    }

    if (dataDst_in == NULL) {
      // Sets data and dataSize
      dataDstSize = getDataSizeRequired(permutationDst, pmeGrid, jblock, kblock);
      dataDst = allocateData(dataDstSize);
      dataDstAllocated = true;
    } else {
      if (dataDstSize_in < getDataSizeRequired(permutationDst, pmeGrid, jblock, kblock))
        NAMD_bug("FFTCompute::init, invalid dataDstSize_in");
      dataDstSize = dataDstSize_in;
      dataDst = dataDst_in;
      dataDstAllocated = false;
    }

    // Final sanity check
    if (dataDst == NULL || dataSrc == NULL ||
      dataDstSize < getDataSizeRequired(permutationDst, pmeGrid, jblock, kblock) ||
      dataSrcSize < getDataSizeRequired(permutation, pmeGrid, jblock, kblock))
      NAMD_bug("FFTCompute::init, error setting up data buffers");

    // Now "data" is pointer to grid data with at least size getDataSizeRequired(...)
    if (pmePencilType == 3) {
      // 3D FFT
      if (pmeGrid.xBlocks != 1 || pmeGrid.yBlocks != 1 || pmeGrid.zBlocks != 1)
        NAMD_bug("FFTCompute::init, 3D FFT requires a single pencil");
      int n[3] = {pmeGrid.K1, pmeGrid.K2, pmeGrid.K3};
      plan3D(n, flags);
    } else if (pmePencilType == 1) {
      int i0, i1, j0, j1, k0, k1;
      getPencilDim(pmeGrid, permutation, jblock, kblock, i0, i1, j0, j1, k0, k1);
      int n[1] = {i1-i0+1};
      int howmany = (j1-j0+1)*(k1-k0+1);
      if (permutation == Perm_X_Y_Z) {
        plan1DX(n, howmany, flags);
      } else if (permutation == Perm_Y_Z_cX) {
        plan1DY(n, howmany, flags);
      } else if (permutation == Perm_Z_cX_Y) {
        plan1DZ(n, howmany, flags);
      } else {
        NAMD_bug("FFTCompute::init, invalid permutation");
      }
    } else if (pmePencilType == 2) {
      // 2D FFTs of xy planes
      int i0, i1, j0, j1, k0, k1;
      getPencilDim(pmeGrid, permutation, 0, kblock, i0, i1, j0, j1, k0, k1);
      int n[2] = {pmeGrid.K1, pmeGrid.K2};
      int howmany = k1-k0+1;
      plan2D(n, howmany, flags);
    } else {
      NAMD_bug("FFTCompute::init, invalid pmePencilType");
    }       
  }
  virtual ~FFTCompute() {}
  virtual void forward()=0;
  virtual void backward()=0;
  float* getDataSrc() {return dataSrc;}
  float* getDataDst() {return dataDst;}
protected:
  int jblock, kblock;
  int isize, jsize, ksize;
  // Pointer to data, allocated here only if dataAllocated = true
  float* dataSrc;
  float* dataDst;
  // Size of data
  int dataSrcSize;
  int dataDstSize;
  // If dataAllocated = true, implementation must deallocate dataPtr
  bool dataSrcAllocated;
  bool dataDstAllocated;
private:
  // Return the size of data in floats that this FFT will require
  int getDataSizeRequired(int permutation, PmeGrid pmeGrid, 
    int jblock, int kblock) {

    this->jblock = jblock;
    this->kblock = kblock;

    int i0, i1, j0, j1, k0, k1;
    getPencilDim(pmeGrid, permutation, jblock, kblock,
      i0, i1, j0, j1, k0, k1);
    
    int size = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);

    isize = (i1-i0+1);
    jsize = (j1-j0+1);
    ksize = (k1-k0+1);

    if (permutation == Perm_X_Y_Z)
      return size;
    else
      return 2*size;
  }
  // Returns pointer to array of size dataSizeRequired
  virtual float* allocateData(const int dataSizeRequired)=0;
  virtual void plan3D(int *n, int flags)=0;
  virtual void plan2D(int *n, int howmany, int flags)=0;
  virtual void plan1DX(int *n, int howmany, int flags)=0;
  virtual void plan1DY(int *n, int howmany, int flags)=0;
  virtual void plan1DZ(int *n, int howmany, int flags)=0;
};

//
// PmeKSpaceCompute base class
//
class PmeKSpaceCompute {
protected:
  PmeGrid pmeGrid;
  //int K1_len, K2_start, K2_end, K3_start, K3_end;
  double *bm1, *bm2, *bm3;
  double kappa;
  const int permutation;
  const int jblock, kblock;
  int size1, size2, size3;
  int j0, k0;
public:
  PmeKSpaceCompute(PmeGrid pmeGrid, const int permutation,
    const int jblock, const int kblock, double kappa) : 
    pmeGrid(pmeGrid), permutation(permutation),
    jblock(jblock), kblock(kblock), kappa(kappa) {

    bm1 = new double[pmeGrid.K1];
    bm2 = new double[pmeGrid.K2];
    bm3 = new double[pmeGrid.K3];
    // Use compute_b_moduli from PmeKSpace.C
    extern void compute_b_moduli(double *bm, int K, int order);
    compute_b_moduli(bm1, pmeGrid.K1, pmeGrid.order);
    compute_b_moduli(bm2, pmeGrid.K2, pmeGrid.order);
    compute_b_moduli(bm3, pmeGrid.K3, pmeGrid.order);

    int i0, i1, j1, k1;
    getPencilDim(pmeGrid, permutation, jblock, kblock, i0, i1, j0, j1, k0, k1);
    size1 = i1-i0+1;
    size2 = j1-j0+1;
    size3 = k1-k0+1;
  }
  virtual ~PmeKSpaceCompute() {
    delete [] bm1;
    delete [] bm2;
    delete [] bm3;
  }
  virtual void solve(Lattice &lattice, const bool doEnergy, const bool doVirial, float* data)=0;
  virtual double getEnergy()=0;
  virtual void getVirial(double *virial)=0;
};

//
// Abstract PmeRealSpaceCompute base class
//
class PmeRealSpaceCompute {
protected:
  // Number of patches and atoms
  // int numPatches;
  int numAtoms;
  // Grid definition
  PmeGrid pmeGrid;
  // Grid (y, z) location
  int y0, z0;
  // Grid size in real-space
  int xsize, ysize, zsize;
  // Grid data
  int dataSize;
  float *data;
  // Pencil position 
  const int jblock, kblock;
public:
  PmeRealSpaceCompute(PmeGrid pmeGrid, const int jblock, const int kblock) : 
    pmeGrid(pmeGrid), jblock(jblock), kblock(kblock), data(NULL) {
      int x0, x1, y1, z1;
      getPencilDim(pmeGrid, Perm_X_Y_Z, jblock, kblock, x0, x1, y0, y1, z0, z1);
      xsize = x1-x0+1;
      ysize = y1-y0+1;
      zsize = z1-z0+1;
      // Allocate enough data for storing the complex data
      // dataSize = 2*(xsize/2+1)*ysize*zsize;
      // Only allocate enough data for storing the real-space data
      // Complex data is stored in FFTCompute
      dataSize = xsize*ysize*zsize;
    }
  virtual ~PmeRealSpaceCompute() {}
  // Setup patches and atoms in preparation for charge spreading and force gathering
  // virtual void setPatchesAtoms(const int numPatches, const PatchInfo* patches,
  //  const int numAtoms, const CudaAtom* atoms)=0;
  // Setup atoms in preparation for charge spreading and force gathering
  virtual void copyAtoms(const int numAtoms, const CudaAtom* atoms)=0;
  // Spread charges on grid
  virtual void spreadCharge(Lattice &lattice)=0;
  // Gather forces off the grid
  virtual void gatherForce(Lattice &lattice, CudaForce* force)=0;
  // // Calculate self energy
  // virtual double calcSelfEnergy()=0;
  float* getData() {return data;}
  int getDataSize() {return dataSize;}

  static inline double calcGridCoord(const double x, const double recip11, const int nfftx) {
    double w;
    w = x*recip11 + 2.0;
    return (double)((double)nfftx*(w - (floor(w + 0.5) - 0.5)));
  }

  static inline void calcGridCoord(const double x, const double y, const double z,
    const double recip11, const double recip22, const double recip33,
    const int nfftx, const int nffty, const int nfftz,
    double &frx, double &fry, double &frz) {
    double w;
    w = x*recip11 + 2.0;
    frx = (double)((double)nfftx*(w - (floor(w + 0.5) - 0.5)));
    w = y*recip22 + 2.0;
    fry = (double)((double)nffty*(w - (floor(w + 0.5) - 0.5)));
    w = z*recip33 + 2.0;
    frz = (double)((double)nfftz*(w - (floor(w + 0.5) - 0.5)));
  }

  static inline void calcGridCoord(const float x, const float y, const float z,
    const float recip11, const float recip22, const float recip33,
    const int nfftx, const int nffty, const int nfftz,
    float &frx, float &fry, float &frz) {
    float w;
    w = x*recip11 + 2.0f;
    frx = (float)(nfftx*(w - (floorf(w + 0.5f) - 0.5f)));
    w = y*recip22 + 2.0f;
    fry = (float)(nffty*(w - (floorf(w + 0.5f) - 0.5f)));
    w = z*recip33 + 2.0f;
    frz = (float)(nfftz*(w - (floorf(w + 0.5f) - 0.5f)));
  }

  static inline void calcGridCoord(const float x, const float y, const float z,
    const int nfftx, const int nffty, const int nfftz, float &frx, float &fry, float &frz) {
    frx = (float)(nfftx)*x;
    fry = (float)(nffty)*y;
    frz = (float)(nfftz)*z;
  }

  static inline void calcGridCoord(const double x, const double y, const double z,
    const int nfftx, const int nffty, const int nfftz, double &frx, double &fry, double &frz) {
    frx = (double)(nfftx)*x;
    fry = (double)(nffty)*y;
    frz = (double)(nfftz)*z;
  }

};

struct float2;
//
// Abstract PmeTranspose base class
//
class PmeTranspose {
protected:
  PmeGrid pmeGrid;
  const int permutation;
  const int jblock, kblock;
  int isize, jsize, ksize;
  int dataSize;
  int nblock;
  std::vector<int> pos;
public:
  PmeTranspose(PmeGrid pmeGrid, const int permutation, const int jblock, const int kblock) : 
    pmeGrid(pmeGrid), permutation(permutation), jblock(jblock), kblock(kblock) {

    int i0, i1, j0, j1, k0, k1;
    getPencilDim(pmeGrid, permutation, jblock, kblock, i0, i1, j0, j1, k0, k1);
    isize = (i1-i0+1);
    jsize = (j1-j0+1);
    ksize = (k1-k0+1);
    dataSize = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);

    switch(permutation) {
      case Perm_cX_Y_Z:
      nblock = pmeGrid.xBlocks;
      break;
      case Perm_Y_Z_cX:
      nblock = pmeGrid.yBlocks;
      break;
      case Perm_Z_cX_Y:
      nblock = pmeGrid.zBlocks;
      break;
      default:
      NAMD_bug("PmeTranspose::PmeTranspose, invalid permutation");
      break;
    }

    // Pos marks the beginning of blocks
    pos.resize(nblock+1);

    int x1;
    for (int iblock=0;iblock < nblock;iblock++) {
      // Get the dimension of the transpose block
      // We are transposing nblock of these blocks and storing them into a pencil of size
      // ny * nz * xsize at location ny*nz*x0
      int x0, y0dummy, y1dummy, z0dummy, z1dummy;
      getBlockDim(pmeGrid, permutation, iblock, jblock, kblock, x0, x1, y0dummy, y1dummy, z0dummy, z1dummy);

      pos[iblock] = x0;
    }
    // Last position begins at x1+1
    pos[nblock] = x1+1;
  }
  virtual ~PmeTranspose() {}
  virtual void transposeXYZtoYZX(const float2* data)=0;
  virtual void transposeXYZtoZXY(const float2* data)=0;
};

#endif // PMESOLVERUTIL_H