#include <stdio.h>
#include <algorithm>
#ifdef NAMD_CUDA
#include <cuda_runtime.h>
#endif
#include "ComputeNonbondedUtil.h"
#include "ComputePmeCUDAMgr.h"
#include "CudaPmeSolver.h"
#include "CudaPmeSolverUtil.h"

#ifdef NAMD_CUDA
extern "C" void CcdCallBacksReset(void *ignored, double curWallTime);  // fix Charm++

void writeComplexToDisk(const float2 *d_data, const int size, const char* filename, cudaStream_t stream) {
  fprintf(stderr, "writeComplexToDisk %d %s\n", size, filename);
  float2* h_data = new float2[size];
  copy_DtoH<float2>(d_data, h_data, size, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  FILE *handle = fopen(filename, "w");
  for (int i=0;i < size;i++)
    fprintf(handle, "%f %f\n", h_data[i].x, h_data[i].y);
  fclose(handle);
  delete [] h_data;
}

void writeHostComplexToDisk(const float2 *h_data, const int size, const char* filename) {
  FILE *handle = fopen(filename, "w");
  for (int i=0;i < size;i++)
    fprintf(handle, "%f %f\n", h_data[i].x, h_data[i].y);
  fclose(handle);
}

void writeRealToDisk(const float *d_data, const int size, const char* filename, cudaStream_t stream) {
  fprintf(stderr, "writeRealToDisk %d %s\n", size, filename);
  float* h_data = new float[size];
  copy_DtoH<float>(d_data, h_data, size, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  FILE *handle = fopen(filename, "w");
  for (int i=0;i < size;i++)
    fprintf(handle, "%f\n", h_data[i]);
  fclose(handle);
  delete [] h_data;
}

void CudaFFTCompute::plan3D(int *n, int flags) {
  cudaCheck(cudaSetDevice(deviceID));
  forwardType = CUFFT_R2C;
  backwardType = CUFFT_C2R;
  cufftCheck(cufftPlan3d(&forwardPlan, n[2], n[1], n[0], CUFFT_R2C));
  cufftCheck(cufftPlan3d(&backwardPlan, n[2], n[1], n[0], CUFFT_C2R));
  cufftCheck(cufftSetCompatibilityMode(forwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  cufftCheck(cufftSetCompatibilityMode(backwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  setStream();
  // plantype = 3;
}

void CudaFFTCompute::plan2D(int *n, int howmany, int flags) {
  cudaCheck(cudaSetDevice(deviceID));
  forwardType = CUFFT_R2C;
  backwardType = CUFFT_C2R;
  int nt[2] = {n[1], n[0]};
  cufftCheck(cufftPlanMany(&forwardPlan, 2, nt, NULL, 1, 0, NULL, 1, 0, CUFFT_R2C, howmany));
  cufftCheck(cufftPlanMany(&backwardPlan, 2, nt, NULL, 1, 0, NULL, 1, 0, CUFFT_C2R, howmany));
  cufftCheck(cufftSetCompatibilityMode(forwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  cufftCheck(cufftSetCompatibilityMode(backwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  setStream();
  // plantype = 2;
}

void CudaFFTCompute::plan1DX(int *n, int howmany, int flags) {
  cudaCheck(cudaSetDevice(deviceID));
  forwardType = CUFFT_R2C;
  backwardType = CUFFT_C2R;
  cufftCheck(cufftPlanMany(&forwardPlan, 1, n, NULL, 0, 0, NULL, 0, 0, CUFFT_R2C, howmany));
  cufftCheck(cufftPlanMany(&backwardPlan, 1, n, NULL, 0, 0, NULL, 0, 0, CUFFT_C2R, howmany));
  cufftCheck(cufftSetCompatibilityMode(forwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  cufftCheck(cufftSetCompatibilityMode(backwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  setStream();
  // plantype = 1;
}

void CudaFFTCompute::plan1DY(int *n, int howmany, int flags) {
  cudaCheck(cudaSetDevice(deviceID));
  forwardType = CUFFT_C2C;
  backwardType = CUFFT_C2C;
  cufftCheck(cufftPlanMany(&forwardPlan, 1, n, NULL, 0, 0, NULL, 0, 0, CUFFT_C2C, howmany));
  cufftCheck(cufftPlanMany(&backwardPlan, 1, n, NULL, 0, 0, NULL, 0, 0, CUFFT_C2C, howmany));
  cufftCheck(cufftSetCompatibilityMode(forwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  cufftCheck(cufftSetCompatibilityMode(backwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  setStream();
  // plantype = 1;
}

void CudaFFTCompute::plan1DZ(int *n, int howmany, int flags) {
  cudaCheck(cudaSetDevice(deviceID));
  forwardType = CUFFT_C2C;
  backwardType = CUFFT_C2C;
  cufftCheck(cufftPlanMany(&forwardPlan, 1, n, NULL, 0, 0, NULL, 0, 0, CUFFT_C2C, howmany));
  cufftCheck(cufftPlanMany(&backwardPlan, 1, n, NULL, 0, 0, NULL, 0, 0, CUFFT_C2C, howmany));
  cufftCheck(cufftSetCompatibilityMode(forwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  cufftCheck(cufftSetCompatibilityMode(backwardPlan, CUFFT_COMPATIBILITY_FFTW_PADDING));
  setStream();
  // plantype = 1;
}

CudaFFTCompute::~CudaFFTCompute() {
  cudaCheck(cudaSetDevice(deviceID));
	cufftCheck(cufftDestroy(forwardPlan));
	cufftCheck(cufftDestroy(backwardPlan));
  if (dataSrcAllocated) deallocate_device<float>(&dataSrc);
  if (dataDstAllocated) deallocate_device<float>(&dataDst);
}

float* CudaFFTCompute::allocateData(const int dataSizeRequired) {
  cudaCheck(cudaSetDevice(deviceID));
  float* tmp = NULL;
  allocate_device<float>(&tmp, dataSizeRequired);
  return tmp;
}

// int ncall = 0;

void CudaFFTCompute::forward() {
  cudaCheck(cudaSetDevice(deviceID));
  // ncall++;
  if (forwardType == CUFFT_R2C) {

    cufftCheck(cufftExecR2C(forwardPlan, (cufftReal *)dataSrc, (cufftComplex *)dataDst));

    // if (ncall == 1) {
    //   writeComplexToDisk((float2 *)dataSrc, (isize/2+1)*jsize*ksize, "dataSrc.txt", stream);
    // }

    // if (ncall == 1 && plantype == 2) {
    //   writeComplexToDisk((float2 *)data, (isize/2+1)*jsize*ksize, "data_fx_fy_z.txt", stream);
    // }

  } else if (forwardType == CUFFT_C2C) {
    // nc2cf++;
    // if (ncall == 1 && nc2cf == 1)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_y_z_fx.txt");
    // else if (ncall == 1 && nc2cf == 2)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_z_fx_fy.txt");
    cufftCheck(cufftExecC2C(forwardPlan, (cufftComplex *)dataSrc, (cufftComplex *)dataDst, CUFFT_FORWARD));
    // fprintf(stderr, "ncall %d plantype %d\n", ncall, plantype);
    // if (ncall == 1 && plantype == 1 && isize == 62) {
    //   writeComplexToDisk((float2 *)data, isize*jsize*(ksize/2+1), "data_fy_z_fx.txt", stream);
    // }
    // if (ncall == 1 && nc2cf == 1)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_fy_z_fx.txt");
    // else if (ncall == 1 && nc2cf == 2)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_fz_fx_fy.txt");
  } else {
    cudaNAMD_bug("CudaFFTCompute::forward(), unsupported FFT type");
  }
}

void CudaFFTCompute::backward() {
  cudaCheck(cudaSetDevice(deviceID));
  if (backwardType == CUFFT_C2R) {
    // if (ncall == 1) {
    //   if (plantype == 1)
    //     writeComplexToDisk((float2 *)data, 33*64*64, "data_fx_by_bz.txt");
    //   else
    //     writeComplexToDisk((float2 *)data, 33*64*64, "data_fx_fy_fz_2.txt");
    // }

    cufftCheck(cufftExecC2R(backwardPlan, (cufftComplex *)dataDst, (cufftReal *)dataSrc));

    // if (ncall == 1)
    //   if (plantype == 1)
    //     writeRealToDisk(data, 64*64*64, "data_bx_by_bz_1D.txt");
    //   else
    //     writeRealToDisk(data, 64*64*64, "data_bx_by_bz_3D.txt");
  } else if (backwardType == CUFFT_C2C) {
    // nc2cb++;
    // if (ncall == 1 && nc2cb == 1)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_fz_fx_fy_2.txt");
    // else if (ncall == 1 && nc2cb == 2)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_fy_bz_fx.txt");
    cufftCheck(cufftExecC2C(backwardPlan, (cufftComplex *)dataDst, (cufftComplex *)dataSrc, CUFFT_INVERSE));
    // if (ncall == 1 && nc2cb == 1)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_bz_fx_fy.txt");
    // else if (ncall == 1 && nc2cb == 2)
    //   writeComplexToDisk((float2 *)data, 33*64*64, "data_by_bz_fx.txt");
  } else {
    cudaNAMD_bug("CudaFFTCompute::backward(), unsupported FFT type");
  }
}

void CudaFFTCompute::setStream() {
  cudaCheck(cudaSetDevice(deviceID));
  cufftCheck(cufftSetStream(forwardPlan, stream));
  cufftCheck(cufftSetStream(backwardPlan, stream));
}

//###########################################################################
//###########################################################################
//###########################################################################

CudaPmeKSpaceCompute::CudaPmeKSpaceCompute(PmeGrid pmeGrid, const int permutation,
  const int jblock, const int kblock, double kappa, int cuda_arch, int deviceID, cudaStream_t stream) : 
  PmeKSpaceCompute(pmeGrid, permutation, jblock, kblock, kappa),
  cuda_arch(cuda_arch), deviceID(deviceID), stream(stream) {

  cudaCheck(cudaSetDevice(deviceID));

  // Copy bm1 -> prefac_x on GPU memory
  float *bm1f = new float[pmeGrid.K1];
  float *bm2f = new float[pmeGrid.K2];
  float *bm3f = new float[pmeGrid.K3];
  for (int i=0;i < pmeGrid.K1;i++) bm1f[i] = (float)bm1[i];
  for (int i=0;i < pmeGrid.K2;i++) bm2f[i] = (float)bm2[i];
  for (int i=0;i < pmeGrid.K3;i++) bm3f[i] = (float)bm3[i];
  allocate_device<float>(&d_bm1, pmeGrid.K1);
  allocate_device<float>(&d_bm2, pmeGrid.K2);
  allocate_device<float>(&d_bm3, pmeGrid.K3);
  copy_HtoD_sync<float>(bm1f, d_bm1, pmeGrid.K1);
  copy_HtoD_sync<float>(bm2f, d_bm2, pmeGrid.K2);
  copy_HtoD_sync<float>(bm3f, d_bm3, pmeGrid.K3);
  delete [] bm1f;
  delete [] bm2f;
  delete [] bm3f;
  allocate_device<EnergyVirial>(&d_energyVirial, 1);
  allocate_host<EnergyVirial>(&h_energyVirial, 1);
  // cudaCheck(cudaEventCreateWithFlags(&copyEnergyVirialEvent, cudaEventDisableTiming));
  cudaCheck(cudaEventCreate(&copyEnergyVirialEvent));
  // ncall = 0;
}

CudaPmeKSpaceCompute::~CudaPmeKSpaceCompute() {
  cudaCheck(cudaSetDevice(deviceID));
  deallocate_device<float>(&d_bm1);
  deallocate_device<float>(&d_bm2);
  deallocate_device<float>(&d_bm3);
  deallocate_device<EnergyVirial>(&d_energyVirial);
  deallocate_host<EnergyVirial>(&h_energyVirial);
  cudaCheck(cudaEventDestroy(copyEnergyVirialEvent));
}

void CudaPmeKSpaceCompute::solve(Lattice &lattice, const bool doEnergy, const bool doVirial, float* data) {
  cudaCheck(cudaSetDevice(deviceID));

  const bool doEnergyVirial = (doEnergy || doVirial);

  int nfft1, nfft2, nfft3;
  float *prefac1, *prefac2, *prefac3;
  float recip1, recip2, recip3;

  Vector recip1v = lattice.a_r();
  Vector recip2v = lattice.b_r();
  Vector recip3v = lattice.c_r();
  double recip[9] = {recip1v.x, recip1v.y, recip1v.z,
    recip2v.x, recip2v.y, recip2v.z,
    recip3v.x, recip3v.y, recip3v.z};

  if (permutation == Perm_Z_cX_Y) {
    // Z, X, Y
    nfft1 = pmeGrid.K3;
    nfft2 = pmeGrid.K1;
    nfft3 = pmeGrid.K2;
    prefac1 = d_bm3;
    prefac2 = d_bm1;
    prefac3 = d_bm2;
    recip1 = (float)recip[8];
    recip2 = (float)recip[0];
    recip3 = (float)recip[4];
  } else if (permutation == Perm_cX_Y_Z) {
    // X, Y, Z
    nfft1 = pmeGrid.K1;
    nfft2 = pmeGrid.K2;
    nfft3 = pmeGrid.K3;
    prefac1 = d_bm1;
    prefac2 = d_bm2;
    prefac3 = d_bm3;
    recip1 = (float)recip[0];
    recip2 = (float)recip[4];
    recip3 = (float)recip[8];
  } else {
    NAMD_bug("CudaPmeKSpaceCompute::solve, invalid permutation");
  }

  ortho = (recip[1] == 0.0 && recip[2] == 0.0 && recip[3] == 0.0 &&
    recip[5] == 0.0 && recip[6] == 0.0 && recip[7] == 0.0);

  if (!ortho)
    cudaNAMD_bug("CudaPmeKSpaceCompute::solve, only orthorhombic boxes are currently supported");

  // ncall++;
  // if (ncall == 1) {
  //   char filename[256];
  //   sprintf(filename,"dataf_%d_%d.txt",jblock,kblock);
  //   writeComplexToDisk((float2*)data, size1*size2*size3, filename, stream);
  // }

  // if (ncall == 1) {
  //   float2* h_data = new float2[size1*size2*size3];
  //   float2* d_data = (float2*)data;
  //   copy_DtoH<float2>(d_data, h_data, size1*size2*size3, stream);
  //   cudaCheck(cudaStreamSynchronize(stream));
  //   FILE *handle = fopen("dataf.txt", "w");
  //   for (int z=0;z < pmeGrid.K3;z++) {
  //     for (int y=0;y < pmeGrid.K2;y++) {
  //       for (int x=0;x < pmeGrid.K1/2+1;x++) {
  //         int i;
  //         if (permutation == Perm_cX_Y_Z) {
  //           i = x + y*size1 + z*size1*size2;
  //         } else {
  //           i = z + x*size1 + y*size1*size2;
  //         }
  //         fprintf(handle, "%f %f\n", h_data[i].x, h_data[i].y);
  //       }
  //     }
  //   }
  //   fclose(handle);
  //   delete [] h_data;
  // }

  // Clear energy and virial array if needed
  if (doEnergyVirial) clear_device_array<EnergyVirial>(d_energyVirial, 1, stream);

  scalar_sum_ortho(permutation == Perm_cX_Y_Z, nfft1, nfft2, nfft3, size1, size2, size3, kappa,
    recip1, recip2, recip3, prefac1, prefac2, prefac3, j0, k0, doEnergyVirial,
    &d_energyVirial->energy, d_energyVirial->virial, (float2*)data, cuda_arch, stream);

  // Copy energy and virial to host if needed
  if (doEnergyVirial) {
    copy_DtoH<EnergyVirial>(d_energyVirial, h_energyVirial, 1, stream);
    cudaCheck(cudaEventRecord(copyEnergyVirialEvent, stream));
    // cudaCheck(cudaStreamSynchronize(stream));
  }

  // if (ncall == 1) {
  //   float2* h_data = new float2[size1*size2*size3];
  //   float2* d_data = (float2*)data;
  //   copy_DtoH<float2>(d_data, h_data, size1*size2*size3, stream);
  //   cudaCheck(cudaStreamSynchronize(stream));
  //   FILE *handle = fopen("datafs.txt", "w");
  //   for (int i=0;i < size1*size2*size3;i++) {
  //     fprintf(handle, "%f %f\n", h_data[i].x, h_data[i].y);
  //   }
  //   fclose(handle);
  //   delete [] h_data;
  // }

  //fprintf(stderr, "CudaPmeKSpaceCompute::solve()\n");
}

// void CudaPmeKSpaceCompute::waitEnergyAndVirial() {
//   cudaCheck(cudaSetDevice(deviceID));
//   cudaCheck(cudaEventSynchronize(copyEnergyVirialEvent));
// }

void CudaPmeKSpaceCompute::energyAndVirialCheck(void *arg, double walltime) {
  CudaPmeKSpaceCompute* c = (CudaPmeKSpaceCompute *)arg;

  cudaError_t err = cudaEventQuery(c->copyEnergyVirialEvent);
  if (err == cudaSuccess) {
    // Event has occurred
    c->checkCount = 0;
    if (c->pencilXYZPtr != NULL)
      c->pencilXYZPtr->energyAndVirialDone();
    else if (c->pencilZPtr != NULL)
      c->pencilZPtr->energyAndVirialDone();
    else
      NAMD_bug("CudaPmeKSpaceCompute::energyAndVirialCheck, pencilXYZPtr and pencilZPtr not set");
    return;
  } else if (err == cudaErrorNotReady) {
    // Event has not occurred
    c->checkCount++;
    if (c->checkCount >= 1000000) {
      NAMD_bug("CudaPmeKSpaceCompute::energyAndVirialCheck, check count exceeded");
    }
  } else {
    // Anything else is an error
    NAMD_bug("CudaPmeKSpaceCompute::energyAndVirialCheck, cudaEventQuery returned error");
  }

  // Call again 
  CcdCallBacksReset(0, walltime);
  CcdCallFnAfter(energyAndVirialCheck, arg, 0.1);
}

void CudaPmeKSpaceCompute::energyAndVirialSetCallback(CudaPmePencilXYZ* pencilPtr) {
  cudaCheck(cudaSetDevice(deviceID));
  pencilXYZPtr = pencilPtr;
  pencilZPtr = NULL;
  checkCount = 0;
  CcdCallBacksReset(0, CmiWallTimer());
  // Set the call back at 0.1ms
  CcdCallFnAfter(energyAndVirialCheck, this, 0.1);
}

void CudaPmeKSpaceCompute::energyAndVirialSetCallback(CudaPmePencilZ* pencilPtr) {
  cudaCheck(cudaSetDevice(deviceID));
  pencilXYZPtr = NULL;
  pencilZPtr = pencilPtr;
  checkCount = 0;
  CcdCallBacksReset(0, CmiWallTimer());
  // Set the call back at 0.1ms
  CcdCallFnAfter(energyAndVirialCheck, this, 0.1);
}

double CudaPmeKSpaceCompute::getEnergy() {
  return h_energyVirial->energy;
}

void CudaPmeKSpaceCompute::getVirial(double *virial) {
  if (ortho) {
    if (permutation == Perm_Z_cX_Y) {
      // h_energyVirial->virial is storing ZZ, ZX, ZY, XX, XY, YY
      virial[0] = h_energyVirial->virial[3];
      virial[1] = h_energyVirial->virial[4];
      virial[2] = h_energyVirial->virial[1];

      virial[3] = h_energyVirial->virial[4];
      virial[4] = h_energyVirial->virial[5];
      virial[5] = h_energyVirial->virial[2];

      virial[6] = h_energyVirial->virial[1];
      virial[7] = h_energyVirial->virial[7];
      virial[8] = h_energyVirial->virial[0];
    } else if (permutation == Perm_cX_Y_Z) {
      // h_energyVirial->virial is storing XX, XY, XZ, YY, YZ, ZZ
      virial[0] = h_energyVirial->virial[0];
      virial[1] = h_energyVirial->virial[1];
      virial[2] = h_energyVirial->virial[2];

      virial[3] = h_energyVirial->virial[1];
      virial[4] = h_energyVirial->virial[3];
      virial[5] = h_energyVirial->virial[4];

      virial[6] = h_energyVirial->virial[2];
      virial[7] = h_energyVirial->virial[4];
      virial[8] = h_energyVirial->virial[5];
    }
  } else {
    NAMD_bug("CudaPmeKSpaceCompute::getVirial, only orthorhombic boxes are currently supported");
  }
}


//###########################################################################
//###########################################################################
//###########################################################################

//
// Class constructor
//
CudaPmeRealSpaceCompute::CudaPmeRealSpaceCompute(PmeGrid pmeGrid,
  const int jblock, const int kblock, int cuda_arch, int deviceID, cudaStream_t stream) : 
  PmeRealSpaceCompute(pmeGrid, jblock, kblock), cuda_arch(cuda_arch), deviceID(deviceID), stream(stream) {
  if (dataSize < xsize*ysize*zsize)
    NAMD_bug("CudaPmeRealSpaceCompute::CudaPmeRealSpaceCompute, insufficient dataSize");
  cudaCheck(cudaSetDevice(deviceID));
  d_atomsCapacity = 0;
  d_atoms = NULL;
  d_forceCapacity = 0;
  d_force = NULL;
  tex_data = NULL;
  tex_data_len = 0;
  allocate_device<float>(&data, dataSize);
  setupGridTexture(data, xsize*ysize*zsize);
  cudaCheck(cudaEventCreate(&gatherForceEvent));
}

//
// Class desctructor
//
CudaPmeRealSpaceCompute::~CudaPmeRealSpaceCompute() {
  cudaCheck(cudaSetDevice(deviceID));
  if (d_atoms != NULL) deallocate_device<CudaAtom>(&d_atoms);
  if (d_force != NULL) deallocate_device<CudaForce>(&d_force);
  // if (d_patches != NULL) deallocate_device<PatchInfo>(&d_patches);
  // deallocate_device<double>(&d_selfEnergy);
  deallocate_device<float>(&data);
  cudaCheck(cudaEventDestroy(gatherForceEvent));
}

// //
// // Copy patches and atoms to device memory
// //
// void CudaPmeRealSpaceCompute::setPatchesAtoms(const int numPatches, const PatchInfo* patches,
//   const int numAtoms, const CudaAtom* atoms) {

//   this->numPatches = numPatches;
//   this->numAtoms = numAtoms;

//   // Reallocate device arrays as neccessary
//   reallocate_device<CudaAtom>(&d_atoms, &d_atomsCapacity, numAtoms, 1.5f);
//   reallocate_device<PatchInfo>(&d_patches, &d_patchesCapacity, numPatches, 1.5f);

//   // Copy atom and patch data to device
//   copy_HtoD<CudaAtom>(atoms, d_atoms, numAtoms, stream);
//   copy_HtoD<PatchInfo>(patches, d_patches, numPatches, stream);
// }

//
// Copy atoms to device memory
//
void CudaPmeRealSpaceCompute::copyAtoms(const int numAtoms, const CudaAtom* atoms) {
  cudaCheck(cudaSetDevice(deviceID));
  this->numAtoms = numAtoms;

  // Reallocate device arrays as neccessary
  reallocate_device<CudaAtom>(&d_atoms, &d_atomsCapacity, numAtoms, 1.5f);

  // Copy atom data to device
  copy_HtoD<CudaAtom>(atoms, d_atoms, numAtoms, stream);
}

//
// Spread charges on grid
//
void CudaPmeRealSpaceCompute::spreadCharge(Lattice &lattice) {
  cudaCheck(cudaSetDevice(deviceID));
  Vector recip1v = lattice.a_r();
  Vector recip2v = lattice.b_r();
  Vector recip3v = lattice.c_r();
  double recip[9] = {recip1v.x, recip1v.y, recip1v.z,
    recip2v.x, recip2v.y, recip2v.z,
    recip3v.x, recip3v.y, recip3v.z};

  // Clear grid
  clear_device_array<float>(data, xsize*ysize*zsize, stream);

  bool ortho = (recip[1] == 0.0 && recip[2] == 0.0 && recip[3] == 0.0 &&
    recip[5] == 0.0 && recip[6] == 0.0 && recip[7] == 0.0);

  if (!ortho) {
    cudaNAMD_bug("CudaPmeRealSpaceCompute::spreadCharge, only orthorhombic boxes are currently supported");
  }

  float recip1 = (float)recip[0];
  float recip2 = (float)recip[4];
  float recip3 = (float)recip[8];

  spread_charge_ortho((const float4*)d_atoms, numAtoms,
    recip1, recip2, recip3, pmeGrid.K1, pmeGrid.K2, pmeGrid.K3, xsize, ysize, zsize,
    xsize, y0, z0, (pmeGrid.yBlocks == 1), (pmeGrid.zBlocks == 1),
    data, pmeGrid.order, stream);

  // ncall++;

  // if (ncall == 1) writeRealToDisk(data, xsize*ysize*zsize, "data.txt");
}

void CudaPmeRealSpaceCompute::cuda_gatherforce_check(void *arg, double walltime) {
  CudaPmeRealSpaceCompute* c = (CudaPmeRealSpaceCompute *)arg;
  cudaCheck(cudaSetDevice(c->deviceID));

  cudaError_t err = cudaEventQuery(c->gatherForceEvent);
  if (err == cudaSuccess) {
    // Event has occurred
    c->checkCount = 0;
//    c->deviceProxy[CkMyNode()].gatherForceDone();
    c->devicePtr->gatherForceDone();
    return;
  } else if (err == cudaErrorNotReady) {
    // Event has not occurred
    c->checkCount++;
    if (c->checkCount >= 1000000) {
      NAMD_bug("CudaPmeRealSpaceCompute::cuda_gatherforce_check, check count exceeded");
    }
  } else {
    // Anything else is an error
    NAMD_bug("CudaPmeRealSpaceCompute::cuda_gatherforce_check, cudaEventQuery returned error");
  }

  // Call again 
  CcdCallBacksReset(0, walltime);
  CcdCallFnAfter(cuda_gatherforce_check, arg, 0.1);
}

void CudaPmeRealSpaceCompute::gatherForceSetCallback(ComputePmeCUDADevice* devicePtr_in) {
  cudaCheck(cudaSetDevice(deviceID));
  devicePtr = devicePtr_in;
  checkCount = 0;
  CcdCallBacksReset(0, CmiWallTimer());
  // Set the call back at 0.1ms
  CcdCallFnAfter(cuda_gatherforce_check, this, 0.1);
}

void CudaPmeRealSpaceCompute::waitGatherForceDone() {
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaEventSynchronize(gatherForceEvent));
}

void CudaPmeRealSpaceCompute::setupGridTexture(float* data, int data_len) {
  if (cuda_arch >= 350 || (tex_data == data && tex_data_len == data_len)) return;
  tex_data = data;
  tex_data_len = data_len;
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  // Use texture objects
  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeLinear;
  resDesc.res.linear.devPtr = data;
  resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
  resDesc.res.linear.desc.x = sizeof(float)*8;
  resDesc.res.linear.sizeInBytes = data_len*sizeof(float);
  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;
  cudaCheck(cudaCreateTextureObject(&gridTexObj, &resDesc, &texDesc, NULL));
#else
  bindGridTexture(data, data_len);
#endif
}

void CudaPmeRealSpaceCompute::gatherForce(Lattice &lattice, CudaForce* force) {
  cudaCheck(cudaSetDevice(deviceID));

  // Re-allocate force array if needed
  reallocate_device<CudaForce>(&d_force, &d_forceCapacity, numAtoms, 1.5f);

  Vector recip1v = lattice.a_r();
  Vector recip2v = lattice.b_r();
  Vector recip3v = lattice.c_r();
  double recip[9] = {recip1v.x, recip1v.y, recip1v.z,
    recip2v.x, recip2v.y, recip2v.z,
    recip3v.x, recip3v.y, recip3v.z};

  bool ortho = (recip[1] == 0.0 && recip[2] == 0.0 && recip[3] == 0.0 &&
    recip[5] == 0.0 && recip[6] == 0.0 && recip[7] == 0.0);

  if (!ortho) {
    cudaNAMD_bug("CudaPmeRealSpaceCompute::gatherForce, only orthorhombic boxes are currently supported");
  }

  float recip1 = (float)recip[0];
  float recip2 = (float)recip[4];
  float recip3 = (float)recip[8];

  gather_force_ortho((const float4*)d_atoms, numAtoms,
    recip1, recip2, recip3, pmeGrid.K1, pmeGrid.K2, pmeGrid.K3,
    xsize, ysize, zsize, xsize, y0, z0, (pmeGrid.yBlocks == 1), (pmeGrid.zBlocks == 1),
    data, pmeGrid.order, (float3*)d_force, 
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
    gridTexObj,
#endif
    stream);

  copy_DtoH<CudaForce>(d_force, force, numAtoms, stream);

  cudaCheck(cudaEventRecord(gatherForceEvent, stream));
}

/*
double CudaPmeRealSpaceCompute::calcSelfEnergy() {
  double h_selfEnergy;
  clear_device_array<double>(d_selfEnergy, 1);
  calc_sum_charge_squared((const float4*)d_atoms, numAtoms, d_selfEnergy, stream);
  copy_DtoH<double>(d_selfEnergy, &h_selfEnergy, 1, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  // 1.7724538509055160273 = sqrt(pi)
  h_selfEnergy *= -ComputeNonbondedUtil::ewaldcof/1.7724538509055160273;
  return h_selfEnergy;
}
*/

//###########################################################################
//###########################################################################
//###########################################################################

CudaPmeTranspose::CudaPmeTranspose(PmeGrid pmeGrid, const int permutation,
    const int jblock, const int kblock, int deviceID, cudaStream_t stream) : 
  PmeTranspose(pmeGrid, permutation, jblock, kblock), deviceID(deviceID), stream(stream) {
  cudaCheck(cudaSetDevice(deviceID));

  allocate_device<float2>(&d_data, dataSize);
#ifndef P2P_ENABLE_3D
  allocate_device<float2>(&d_buffer, dataSize);
#endif

  // Setup data pointers to NULL, these can be overridden later on by using setDataPtrs()
  dataPtrsYZX.resize(nblock, NULL);
  dataPtrsZXY.resize(nblock, NULL);

  allocate_device< TransposeBatch<float2> >(&batchesYZX, 3*nblock);
  allocate_device< TransposeBatch<float2> >(&batchesZXY, 3*nblock);
}

CudaPmeTranspose::~CudaPmeTranspose() {
  cudaCheck(cudaSetDevice(deviceID));
  deallocate_device<float2>(&d_data);
#ifndef P2P_ENABLE_3D
  deallocate_device<float2>(&d_buffer);
#endif
  deallocate_device< TransposeBatch<float2> >(&batchesZXY);
  deallocate_device< TransposeBatch<float2> >(&batchesYZX);
}

//
// Set dataPtrsYZX
//
void CudaPmeTranspose::setDataPtrsYZX(std::vector<float2*>& dataPtrsNew, float2* data) {
  if (dataPtrsYZX.size() != dataPtrsNew.size())
    NAMD_bug("CudaPmeTranspose::setDataPtrsYZX, invalid dataPtrsNew size");
  for (int iblock=0;iblock < nblock;iblock++) {
    dataPtrsYZX[iblock] = dataPtrsNew[iblock];
  }
  // Build batched data structures
  TransposeBatch<float2> *h_batchesYZX = new TransposeBatch<float2>[3*nblock];

  for (int iperm=0;iperm < 3;iperm++) {
    int isize_out;
    if (iperm == 0) {
      // Perm_Z_cX_Y:
      // ZXY -> XYZ
      isize_out = pmeGrid.K1/2+1;
    } else if (iperm == 1) {
      // Perm_cX_Y_Z:
      // XYZ -> YZX
      isize_out = pmeGrid.K2;
    } else {
      // Perm_Y_Z_cX:
      // YZX -> ZXY
      isize_out = pmeGrid.K3;
    }

    int max_nx = 0;
    for (int iblock=0;iblock < nblock;iblock++) {

      int x0 = pos[iblock];
      int nx = pos[iblock+1] - x0;
      max_nx = std::max(max_nx, nx);

      int width_out;
      float2* data_out;
      if (dataPtrsYZX[iblock] == NULL) {
        // Local transpose, use internal buffer
        data_out = d_data + jsize*ksize*x0;
        width_out = jsize;
      } else {
        // Non-local tranpose, use buffer in dataPtr[] and the size of that buffer
        data_out = dataPtrsYZX[iblock];
        width_out = isize_out;
      }

      TransposeBatch<float2> batch;
      batch.nx        = nx;
      batch.ysize_out = width_out;
      batch.zsize_out = ksize;
      batch.data_in   = data+x0;
      batch.data_out  = data_out;

      h_batchesYZX[iperm*nblock + iblock] = batch;

    // transpose_xyz_yzx(
    //   nx, jsize, ksize,
    //   isize, jsize,
    //   width_out, ksize,
    //   data+x0, data_out, stream);
    }

    max_nx_YZX[iperm] = max_nx;
  }

  copy_HtoD< TransposeBatch<float2> >(h_batchesYZX, batchesYZX, 3*nblock, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  delete [] h_batchesYZX;
}

//
// Set dataPtrsZXY
//
void CudaPmeTranspose::setDataPtrsZXY(std::vector<float2*>& dataPtrsNew, float2* data) {
  if (dataPtrsZXY.size() != dataPtrsNew.size())
    NAMD_bug("CudaPmeTranspose::setDataPtrsZXY, invalid dataPtrsNew size");
  for (int iblock=0;iblock < nblock;iblock++) {
    dataPtrsZXY[iblock] = dataPtrsNew[iblock];
  }

  // Build batched data structures
  TransposeBatch<float2> *h_batchesZXY = new TransposeBatch<float2>[3*nblock];

  for (int iperm=0;iperm < 3;iperm++) {
    int isize_out;
    if (iperm == 0) {
      // Perm_cX_Y_Z:
      // XYZ -> ZXY
      isize_out = pmeGrid.K3;
    } else if (iperm == 1) {
      // Perm_Z_cX_Y:
      // ZXY -> YZX
      isize_out = pmeGrid.K2;
    } else {
      // Perm_Y_Z_cX:
      // YZX -> XYZ
      isize_out = pmeGrid.K1/2+1;
    }

    int max_nx = 0;
    for (int iblock=0;iblock < nblock;iblock++) {

      int x0 = pos[iblock];
      int nx = pos[iblock+1] - x0;
      max_nx = std::max(max_nx, nx);

      int width_out;
      float2* data_out;
      if (dataPtrsZXY[iblock] == NULL) {
        // Local transpose, use internal buffer
        data_out = d_data + jsize*ksize*x0;
        width_out = ksize;
      } else {
        // Non-local tranpose, use buffer in dataPtr[] and the size of that buffer
        data_out = dataPtrsZXY[iblock];
        width_out = isize_out;
      }

      TransposeBatch<float2> batch;
      batch.nx        = nx;
      batch.zsize_out = width_out;
      batch.xsize_out = nx;
      batch.data_in   = data+x0;
      batch.data_out  = data_out;

      h_batchesZXY[iperm*nblock + iblock] = batch;
    }

    max_nx_ZXY[iperm] = max_nx;
  }

  copy_HtoD< TransposeBatch<float2> >(h_batchesZXY, batchesZXY, 3*nblock, stream);
  cudaCheck(cudaStreamSynchronize(stream));
  delete [] h_batchesZXY;
}

void CudaPmeTranspose::transposeXYZtoYZX(const float2* data) {
  cudaCheck(cudaSetDevice(deviceID));

  int iperm;
  switch(permutation) {
    case Perm_Z_cX_Y:
    // ZXY -> XYZ
    iperm = 0;
    break;
    case Perm_cX_Y_Z:
    // XYZ -> YZX
    iperm = 1;
    break;
    case Perm_Y_Z_cX:
    // YZX -> ZXY
    iperm = 2;
    break;
    default:
    NAMD_bug("PmeTranspose::transposeXYZtoYZX, invalid permutation");
    break;
  }

  batchTranspose_xyz_yzx(
    nblock, batchesYZX + iperm*nblock,
    max_nx_YZX[iperm], jsize, ksize,
    isize, jsize, stream);


/*
  int isize_out;
  switch(permutation) {
    case Perm_Z_cX_Y:
    // ZXY -> XYZ
    isize_out = pmeGrid.K1/2+1;
    break;
    case Perm_cX_Y_Z:
    // XYZ -> YZX
    isize_out = pmeGrid.K2;
    break;
    case Perm_Y_Z_cX:
    // YZX -> ZXY
    isize_out = pmeGrid.K3;
    break;
    default:
    NAMD_bug("PmeTranspose::transposeXYZtoYZX, invalid permutation");
    break;
  }

  for (int iblock=0;iblock < nblock;iblock++) {

    int x0 = pos[iblock];
    int nx = pos[iblock+1] - x0;

    int width_out;
    float2* data_out;
    if (dataPtrsYZX[iblock] == NULL) {
      // Local transpose, use internal buffer
      data_out = d_data + jsize*ksize*x0;
      width_out = jsize;
    } else {
      // Non-local tranpose, use buffer in dataPtr[] and the size of that buffer
      data_out = dataPtrsYZX[iblock];
      width_out = isize_out;
    }

    transpose_xyz_yzx(
      nx, jsize, ksize,
      isize, jsize,
      width_out, ksize,
      data+x0, data_out, stream);
  }
*/
}

void CudaPmeTranspose::transposeXYZtoZXY(const float2* data) {
  cudaCheck(cudaSetDevice(deviceID));

  int iperm;
  switch(permutation) {
    case Perm_cX_Y_Z:
    // XYZ -> ZXY
    iperm = 0;
    break;
    case Perm_Z_cX_Y:
    // ZXY -> YZX
    iperm = 1;
    break;
    case Perm_Y_Z_cX:
    // YZX -> XYZ
    iperm = 2;
    break;
    default:
    NAMD_bug("PmeTranspose::transposeXYZtoZXY, invalid permutation");
    break;
  }

  batchTranspose_xyz_zxy(
    nblock, batchesZXY + iperm*nblock,
    max_nx_ZXY[iperm], jsize, ksize,
    isize, jsize, stream);

/*
  int isize_out;
  switch(permutation) {
    case Perm_cX_Y_Z:
    // XYZ -> ZXY
    isize_out = pmeGrid.K3;
    break;
    case Perm_Z_cX_Y:
    // ZXY -> YZX
    isize_out = pmeGrid.K2;
    break;
    case Perm_Y_Z_cX:
    // YZX -> XYZ
    isize_out = pmeGrid.K1/2+1;
    break;
    default:
    NAMD_bug("PmeTranspose::transposeXYZtoZXY, invalid permutation");
    break;
  }

  for (int iblock=0;iblock < nblock;iblock++) {

    int x0 = pos[iblock];
    int nx = pos[iblock+1] - x0;

    int width_out;
    float2* data_out;
    if (dataPtrsZXY[iblock] == NULL) {
      // Local transpose, use internal buffer
      data_out = d_data + jsize*ksize*x0;
      width_out = ksize;
    } else {
      // Non-local tranpose, use buffer in dataPtr[] and the size of that buffer
      data_out = dataPtrsZXY[iblock];
      width_out = isize_out;
    }

    transpose_xyz_zxy(
      nx, jsize, ksize,
      isize, jsize,
      width_out, nx,
      data+x0, data_out, stream);
  }
*/
}

void CudaPmeTranspose::waitStreamSynchronize() {
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaStreamSynchronize(stream));
}

void CudaPmeTranspose::copyDataDeviceToHost(const int iblock, float2* h_data, const int h_dataSize) {
  cudaCheck(cudaSetDevice(deviceID));

  if (iblock >= nblock)
    NAMD_bug("CudaPmeTranspose::copyDataDeviceToHost, block index exceeds number of blocks");

  int x0 = pos[iblock];
  int nx = pos[iblock+1] - x0;

  int copySize  = jsize*ksize*nx;
  int copyStart = jsize*ksize*x0;

  if (copyStart + copySize > dataSize)
    NAMD_bug("CudaPmeTranspose::copyDataDeviceToHost, dataSize exceeded");

  if (copySize > h_dataSize) 
    NAMD_bug("CudaPmeTranspose::copyDataDeviceToHost, h_dataSize exceeded");

  copy_DtoH<float2>(d_data+copyStart, h_data, copySize, stream);
}

void CudaPmeTranspose::copyDataHostToDevice(const int iblock, float2* data_in, float2* data_out) {
  cudaCheck(cudaSetDevice(deviceID));

  if (iblock >= nblock)
    NAMD_bug("CudaPmeTranspose::copyDataHostToDevice, block index exceeds number of blocks");

  // Determine block size = how much we're copying
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, permutation, iblock, jblock, kblock, i0, i1, j0, j1, k0, k1);
  int ni = i1-i0+1;
  int nj = j1-j0+1;
  int nk = k1-k0+1;

  copy3D_HtoD<float2>(data_in, data_out,
    0, 0, 0,
    ni, nj,
    i0, 0, 0,
    isize, jsize,
    ni, nj, nk, stream);
}

#ifndef P2P_ENABLE_3D
//
// Copy from temporary buffer to final buffer
//
void CudaPmeTranspose::copyDataDeviceToDevice(const int iblock, float2* data_out) {
  cudaCheck(cudaSetDevice(deviceID));

  if (iblock >= nblock)
    NAMD_bug("CudaPmeTranspose::copyDataDeviceToDevice, block index exceeds number of blocks");

  // Determine block size = how much we're copying
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, permutation, iblock, jblock, kblock, i0, i1, j0, j1, k0, k1);
  int ni = i1-i0+1;
  int nj = j1-j0+1;
  int nk = k1-k0+1;

  float2* data_in = d_buffer + i0*nj*nk;

  copy3D_DtoD<float2>(data_in, data_out,
    0, 0, 0,
    ni, nj,
    i0, 0, 0,
    isize, jsize,
    ni, nj, nk, stream);
}

//
// Return temporary buffer for block "iblock"
//
float2* CudaPmeTranspose::getBuffer(const int iblock) {
  if (iblock >= nblock)
    NAMD_bug("CudaPmeTranspose::getBuffer, block index exceeds number of blocks");

  // Determine block size = how much we're copying
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, permutation, iblock, jblock, kblock, i0, i1, j0, j1, k0, k1);
  int ni = i1-i0+1;
  int nj = j1-j0+1;
  int nk = k1-k0+1;

  return d_buffer + i0*nj*nk;
}
#endif

void CudaPmeTranspose::copyDataToPeerDeviceYZX(const int iblock, int deviceID_out, int permutation_out,
  float2* data_out) {

  int iblock_out = jblock;
  int jblock_out = kblock;
  int kblock_out = iblock;

  copyDataToPeerDevice(iblock, iblock_out, jblock_out, kblock_out, deviceID_out, permutation_out, data_out);
}

void CudaPmeTranspose::copyDataToPeerDeviceZXY(const int iblock, int deviceID_out, int permutation_out,
  float2* data_out) {

  int iblock_out = kblock;
  int jblock_out = iblock;
  int kblock_out = jblock;

  copyDataToPeerDevice(iblock, iblock_out, jblock_out, kblock_out, deviceID_out, permutation_out, data_out);
}

void CudaPmeTranspose::copyDataToPeerDevice(const int iblock,
  const int iblock_out, const int jblock_out, const int kblock_out,
  int deviceID_out, int permutation_out, float2* data_out) {

  cudaCheck(cudaSetDevice(deviceID));

  // Determine block size = how much we're copying
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, permutation_out, iblock_out, jblock_out, kblock_out, i0, i1, j0, j1, k0, k1);
  int ni = i1-i0+1;
  int nj = j1-j0+1;
  int nk = k1-k0+1;

  getPencilDim(pmeGrid, permutation_out, jblock_out, kblock_out, i0, i1, j0, j1, k0, k1);
  int isize_out = i1-i0+1;
  int jsize_out = j1-j0+1;

  int x0 = pos[iblock];
  float2* data_in = d_data + jsize*ksize*x0;

#ifndef P2P_ENABLE_3D
  // Copy into temporary peer device buffer
  copy_PeerDtoD<float2>(deviceID, deviceID_out, data_in, data_out, ni*nj*nk, stream);
#else
  copy3D_PeerDtoD<float2>(deviceID, deviceID_out, data_in, data_out,
    0, 0, 0,
    ni, nj,
    0, 0, 0,
    isize_out, jsize_out,
    ni, nj, nk, stream);
#endif

}

#endif // NAMD_CUDA
