#include <stdio.h>
#include "CudaUtils.h"

#ifdef NAMD_CUDA
#include <cuda.h>

__global__ void read_CUDA_ARCH_kernel(int *cuda_arch) {
  if (threadIdx.x == 0) {

#if __CUDA_ARCH__ == 100
    *cuda_arch = 100;
#elif __CUDA_ARCH__ == 110
    *cuda_arch = 110;
#elif __CUDA_ARCH__ == 120
    *cuda_arch = 120;
#elif __CUDA_ARCH__ == 130
    *cuda_arch = 130;
#elif __CUDA_ARCH__ == 200
    *cuda_arch = 200;
#elif __CUDA_ARCH__ == 210
    *cuda_arch = 210;
#elif __CUDA_ARCH__ == 300
    *cuda_arch = 300;
#elif __CUDA_ARCH__ == 350
    *cuda_arch = 350;
#elif __CUDA_ARCH__ == 500
    *cuda_arch = 500;
#else
    *cuda_arch = 500;
#endif

  }
}

//
// Reads the value of __CUDA_ARCH__ from device code
//
int read_CUDA_ARCH() {
  int *d_cuda_arch;
  int h_cuda_arch;
  allocate_device<int>(&d_cuda_arch, 1);
  
  read_CUDA_ARCH_kernel <<< 1, 1 >>> (d_cuda_arch);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    char str[1024];
    sprintf(str, "Error executing CUDA kernel read_CUDA_ARCH_kernel in file %s\nError string: %s\nPossible cause: Device compute capability is less than the compute capability the code was compiled for.\n",
      __FILE__,cudaGetErrorString(err));
    cudaDie(str);
  }
  cudaCheck(cudaDeviceSynchronize());

  copy_DtoH_sync<int>(d_cuda_arch, &h_cuda_arch, 1);
  deallocate_device<int>(&d_cuda_arch);

  return h_cuda_arch;
}

#endif // NAMD_CUDA