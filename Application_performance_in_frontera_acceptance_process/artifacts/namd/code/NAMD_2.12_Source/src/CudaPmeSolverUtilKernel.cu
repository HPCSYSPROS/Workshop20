#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <stdio.h>
#ifdef NAMD_CUDA
#include <cuda.h>
#endif
#include "CudaUtils.h"
#include "CudaPmeSolverUtilKernel.h"

#ifdef NAMD_CUDA
// CCELEC is 1/ (4 pi eps ) in AKMA units, conversion from SI
// units: CCELEC = e*e*Na / (4*pi*eps*1Kcal*1A)
//
//      parameter :: CCELEC=332.0636D0 ! old value of dubious origin
//      parameter :: CCELEC=331.843D0  ! value from 1986-1987 CRC Handbook
//                                   ! of Chemistry and Physics
//  real(chm_real), parameter ::  &
//       CCELEC_amber    = 332.0522173D0, &
//       CCELEC_charmm   = 332.0716D0   , &
//       CCELEC_discover = 332.054D0    , &
//       CCELEC_namd     = 332.0636D0   
//const double ccelec = 332.0636;
//const double half_ccelec = 0.5*ccelec;
//const float ccelec_float = 332.0636f;

/*
// Structure into which virials are stored
// NOTE: This structure is only used for computing addresses
struct Virial_t {
  double sforce_dp[27][3];
  long long int sforce_fp[27][3];
  double virmat[9];
  // Energies start here ...
};
*/

// Local structure for scalar_sum -function for energy and virial reductions
struct RecipVirial_t {
  double energy;
  double virial[6];
};

//
// Performs scalar sum on data(nfft1, nfft2, nfft3)
// T = float or double
// T2 = float2 or double2
//
template <typename T, typename T2, bool calc_energy_virial, bool orderXYZ>
__global__ void scalar_sum_ortho_kernel(const int nfft1, const int nfft2, const int nfft3,
          const int size1, const int size2, const int size3,
          const int nf1, const int nf2, const int nf3,
          const T recip11, const T recip22, const T recip33,
          const T* prefac1, const T* prefac2, const T* prefac3,
          const T fac, const T piv_inv,
          const int k2_00, const int k3_00,
          T2* data,
          double* __restrict__ energy_recip,
          double* __restrict__ virial) {
  // Shared memory required: sizeof(T)*(nfft1 + nfft2 + nfft3)
  extern __shared__ T sh_prefac[];

  // Create pointers to shared memory
  T* sh_prefac1 = (T *)&sh_prefac[0];
  T* sh_prefac2 = (T *)&sh_prefac[nfft1];
  T* sh_prefac3 = (T *)&sh_prefac[nfft1 + nfft2];

  // Calculate start position (k1, k2, k3) for each thread
  unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int k3 = tid/(size1*size2);
  tid -= k3*size1*size2;
  int k2 = tid/size1;
  int k1 = tid - k2*size1;

  // Starting position in data
  int pos = k1 + (k2 + k3*size2)*size1;

  // Move (k2, k3) to the global coordinate (k1_00 = 0 since this is the pencil direction)
  k2 += k2_00;
  k3 += k3_00;

  // Calculate limits w.r.t. global coordinates
  const int lim2 = size2 + k2_00;
  const int lim3 = size3 + k3_00;

  // Calculate increments (k1_inc, k2_inc, k3_inc)
  int tot_inc = blockDim.x*gridDim.x;
  const int k3_inc = tot_inc/(size1*size2);
  tot_inc -= k3_inc*size1*size2;
  const int k2_inc = tot_inc/size1;
  const int k1_inc = tot_inc - k2_inc*size1;

  // Set data[0] = 0 for the global (0,0,0)
  if (k1 == 0 && k2 == 0 && k3 == 0) {
    T2 zero;
    zero.x = (T)0;
    zero.y = (T)0;
    data[0] = zero;
    // Increment position
    k1 += k1_inc;
    pos += k1_inc;
    if (k1 >= size1) {
      k1 -= size1;
      k2++;
    }
    k2 += k2_inc;
    pos += k2_inc*size1;
    if (k2 >= lim2) {
      k2 -= size2;
      k3++;
    }
    k3 += k3_inc;
    pos += k3_inc*size1*size2;
  }

  // Load prefac data into shared memory
  {
    int t = threadIdx.x;
    while (t < nfft1) {
      sh_prefac1[t] = prefac1[t];
      t += blockDim.x;
    }
    t = threadIdx.x;
    while (t < nfft2) {
      sh_prefac2[t] = prefac2[t];
      t += blockDim.x;
    }
    t = threadIdx.x;
    while (t < nfft3) {
      sh_prefac3[t] = prefac3[t];
      t += blockDim.x;
    }
  }
  __syncthreads();

  double energy = 0.0;
  double virial0 = 0.0;
  double virial1 = 0.0;
  double virial2 = 0.0;
  double virial3 = 0.0;
  double virial4 = 0.0;
  double virial5 = 0.0;

  while (k3 < lim3) {

    T2 q = data[pos];

    int m1 = k1;
    int m2 = k2;
    int m3 = k3;
    if (k1 >= nf1) m1 -= nfft1;
    if (k2 >= nf2) m2 -= nfft2;
    if (k3 >= nf3) m3 -= nfft3;

    T mhat1 = recip11*m1;
    T mhat2 = recip22*m2;
    T mhat3 = recip33*m3;

    T msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3;
    T msq_inv = (T)1.0/msq;

    // NOTE: check if it's faster to pre-calculate exp()
    T eterm = exp(-fac*msq)*piv_inv*sh_prefac1[k1]*sh_prefac2[k2]*sh_prefac3[k3]*msq_inv;

    if (calc_energy_virial) {
      T tmp1  = eterm*(q.x*q.x + q.y*q.y);
      T vterm  = ((T)2)*(fac + msq_inv);
      T tmp2   = tmp1*vterm;

      energy += (double)tmp1;
      virial0 += (double)(tmp1*(vterm*mhat1*mhat1 - ((T)1)));
      virial1 += (double)(tmp2*mhat1*mhat2);
      virial2 += (double)(tmp2*mhat1*mhat3);
      virial3 += (double)(tmp1*(vterm*mhat2*mhat2 - ((T)1)));
      virial4 += (double)(tmp2*mhat2*mhat3);
      virial5 += (double)(tmp1*(vterm*mhat3*mhat3 - ((T)1)));

      // The following is put into a separate if {} -block to avoid divergence within warp and
      // save registers
      if ((orderXYZ && k1 >= 1 && k1 < nfft1) ||
        (!orderXYZ && k2 >= 1 && k2 < nfft2)) {

        int k1s, k2s, k3s;
        if (orderXYZ) {
          k1s = nfft1 - (k1+1) + 1;
          k2s = ((nfft2-(k2+1)+1) % nfft2);
          k3s = ((nfft3-(k3+1)+1) % nfft3);
        } else {
          k1s = ((nfft1-(k1+1)+1) % nfft1);
          k2s = nfft2 - (k2+1) + 1;
          k3s = ((nfft3-(k3+1)+1) % nfft3);          
        }

        int m1s = k1s;
        int m2s = k2s;
        int m3s = k3s;

        if (k1s >= nf1) m1s -= nfft1;
        if (k2s >= nf2) m2s -= nfft2;
        if (k3s >= nf3) m3s -= nfft3;

        T mhat1s = recip11*m1s;
        T mhat2s = recip22*m2s;
        T mhat3s = recip33*m3s;

        T msqs = mhat1s*mhat1s + mhat2s*mhat2s + mhat3s*mhat3s;
        T msqs_inv = ((T)1)/msqs;

        T eterms = exp(-fac*msqs)*piv_inv*sh_prefac1[k1s]*sh_prefac2[k2s]*sh_prefac3[k3s]*msqs_inv;

        T tmp1s  = eterms*(q.x*q.x + q.y*q.y);
        T vterms  = ((T)2)*(fac + msqs_inv);
        T tmp2s   = tmp1s*vterms;

        energy += (double)tmp1s;
        virial0 += (double)(tmp1s*(vterms*mhat1s*mhat1s - ((T)1)));
        virial1 += (double)(tmp2s*mhat1s*mhat2s);
        virial2 += (double)(tmp2s*mhat1s*mhat3s);
        virial3 += (double)(tmp1s*(vterms*mhat2s*mhat2s - ((T)1)));
        virial4 += (double)(tmp2s*mhat2s*mhat3s);
        virial5 += (double)(tmp1s*(vterms*mhat3s*mhat3s - ((T)1)));
      }

    }

    q.x *= eterm;
    q.y *= eterm;
    data[pos] = q;
    
    // Increment position
    k1 += k1_inc;
    pos += k1_inc;
    if (k1 >= size1) {
      k1 -= size1;
      k2++;
    }
    k2 += k2_inc;
    pos += k2_inc*size1;
    if (k2 >= lim2) {
      k2 -= size2;
      k3++;
    }
    k3 += k3_inc;
    pos += k3_inc*size2*size1;
  }

  // Reduce energy and virial
  if (calc_energy_virial) {
#if __CUDA_ARCH__ < 300
    // Requires blockDim.x*sizeof(RecipVirial_t) amount of shared memory
    volatile RecipVirial_t* sh_ev = (RecipVirial_t *)sh_prefac;
    // NOTE: this __syncthreads() is needed because we're using a single shared memory buffer
    __syncthreads();
    sh_ev[threadIdx.x].energy  = energy;
    sh_ev[threadIdx.x].virial[0] = virial0;
    sh_ev[threadIdx.x].virial[1] = virial1;
    sh_ev[threadIdx.x].virial[2] = virial2;
    sh_ev[threadIdx.x].virial[3] = virial3;
    sh_ev[threadIdx.x].virial[4] = virial4;
    sh_ev[threadIdx.x].virial[5] = virial5;
    __syncthreads();
#endif
#if __CUDA_ARCH__ < 300
    for (int d=1;d < blockDim.x;d *= 2) {
      int t = threadIdx.x + d;
      double energy_val = (t < blockDim.x) ? sh_ev[t].energy : 0.0;
      double virial0_val = (t < blockDim.x) ? sh_ev[t].virial[0] : 0.0;
      double virial1_val = (t < blockDim.x) ? sh_ev[t].virial[1] : 0.0;
      double virial2_val = (t < blockDim.x) ? sh_ev[t].virial[2] : 0.0;
      double virial3_val = (t < blockDim.x) ? sh_ev[t].virial[3] : 0.0;
      double virial4_val = (t < blockDim.x) ? sh_ev[t].virial[4] : 0.0;
      double virial5_val = (t < blockDim.x) ? sh_ev[t].virial[5] : 0.0;
      __syncthreads();
      sh_ev[threadIdx.x].energy += energy_val;
      sh_ev[threadIdx.x].virial[0] += virial0_val;
      sh_ev[threadIdx.x].virial[1] += virial1_val;
      sh_ev[threadIdx.x].virial[2] += virial2_val;
      sh_ev[threadIdx.x].virial[3] += virial3_val;
      sh_ev[threadIdx.x].virial[4] += virial4_val;
      sh_ev[threadIdx.x].virial[5] += virial5_val;
      __syncthreads();
    }
#else
    const int tid = threadIdx.x & (warpSize-1);
    const int base = (threadIdx.x/warpSize);
    volatile RecipVirial_t* sh_ev = (RecipVirial_t *)sh_prefac;
    // Reduce within warps
    for (int d=warpSize/2;d >= 1;d /= 2) {
      energy += __hiloint2double(__shfl(__double2hiint(energy), tid+d),
         __shfl(__double2loint(energy), tid+d));
      virial0 += __hiloint2double(__shfl(__double2hiint(virial0), tid+d),
          __shfl(__double2loint(virial0), tid+d));
      virial1 += __hiloint2double(__shfl(__double2hiint(virial1), tid+d),
          __shfl(__double2loint(virial1), tid+d));
      virial2 += __hiloint2double(__shfl(__double2hiint(virial2), tid+d),
          __shfl(__double2loint(virial2), tid+d));
      virial3 += __hiloint2double(__shfl(__double2hiint(virial3), tid+d),
          __shfl(__double2loint(virial3), tid+d));
      virial4 += __hiloint2double(__shfl(__double2hiint(virial4), tid+d),
          __shfl(__double2loint(virial4), tid+d));
      virial5 += __hiloint2double(__shfl(__double2hiint(virial5), tid+d),
          __shfl(__double2loint(virial5), tid+d));
    }
    // Reduce between warps
    // NOTE: this __syncthreads() is needed because we're using a single shared memory buffer
    __syncthreads();
    if (tid == 0) {
      sh_ev[base].energy = energy;
      sh_ev[base].virial[0] = virial0;
      sh_ev[base].virial[1] = virial1;
      sh_ev[base].virial[2] = virial2;
      sh_ev[base].virial[3] = virial3;
      sh_ev[base].virial[4] = virial4;
      sh_ev[base].virial[5] = virial5;
    }
    __syncthreads();
    if (base == 0) {
      energy = (tid < blockDim.x/warpSize) ? sh_ev[tid].energy : 0.0;
      virial0 = (tid < blockDim.x/warpSize) ? sh_ev[tid].virial[0] : 0.0;
      virial1 = (tid < blockDim.x/warpSize) ? sh_ev[tid].virial[1] : 0.0;
      virial2 = (tid < blockDim.x/warpSize) ? sh_ev[tid].virial[2] : 0.0;
      virial3 = (tid < blockDim.x/warpSize) ? sh_ev[tid].virial[3] : 0.0;
      virial4 = (tid < blockDim.x/warpSize) ? sh_ev[tid].virial[4] : 0.0;
      virial5 = (tid < blockDim.x/warpSize) ? sh_ev[tid].virial[5] : 0.0;
      for (int d=warpSize/2;d >= 1;d /= 2) {
  energy += __hiloint2double(__shfl(__double2hiint(energy), tid+d),
           __shfl(__double2loint(energy), tid+d));
  virial0 += __hiloint2double(__shfl(__double2hiint(virial0), tid+d),
            __shfl(__double2loint(virial0), tid+d));
  virial1 += __hiloint2double(__shfl(__double2hiint(virial1), tid+d),
            __shfl(__double2loint(virial1), tid+d));
  virial2 += __hiloint2double(__shfl(__double2hiint(virial2), tid+d),
            __shfl(__double2loint(virial2), tid+d));
  virial3 += __hiloint2double(__shfl(__double2hiint(virial3), tid+d),
            __shfl(__double2loint(virial3), tid+d));
  virial4 += __hiloint2double(__shfl(__double2hiint(virial4), tid+d),
            __shfl(__double2loint(virial4), tid+d));
  virial5 += __hiloint2double(__shfl(__double2hiint(virial5), tid+d),
            __shfl(__double2loint(virial5), tid+d));
      }
    }
    
#endif

    if (threadIdx.x == 0) {
#if __CUDA_ARCH__ < 300
      energy = sh_ev[0].energy;
      virial0 = sh_ev[0].virial[0];
      virial1 = sh_ev[0].virial[1];
      virial2 = sh_ev[0].virial[2];
      virial3 = sh_ev[0].virial[3];
      virial4 = sh_ev[0].virial[4];
      virial5 = sh_ev[0].virial[5];
#endif
      atomicAdd(energy_recip, energy*0.5);
      virial0 *= -0.5;
      virial1 *= -0.5;
      virial2 *= -0.5;
      virial3 *= -0.5;
      virial4 *= -0.5;
      virial5 *= -0.5;
      atomicAdd(&virial[0], virial0);
      atomicAdd(&virial[1], virial1);
      atomicAdd(&virial[2], virial2);
      atomicAdd(&virial[3], virial3);
      atomicAdd(&virial[4], virial4);
      atomicAdd(&virial[5], virial5);
    }

  }

}

#ifdef DISABLE_CUDA_TEXTURE_OBJECTS
texture<float, 1, cudaReadModeElementType> gridTexRef;
#endif

template <typename T>
__forceinline__ __device__ void write_grid(const float val, const int ind,
             T* data) {
  atomicAdd(&data[ind], (T)val);
}

//
// General version for any order
//
template <typename T, int order>
__forceinline__ __device__ void calc_one_theta(const T w, T *theta) {

  theta[order-1] = ((T)0);
  theta[1] = w;
  theta[0] = ((T)1) - w;

#pragma unroll
  for (int k=3;k <= order-1;k++) {
    T div = ((T)1) / (T)(k-1);
    theta[k-1] = div*w*theta[k-2];
#pragma unroll
    for (int j=1;j <= k-2;j++) {
      theta[k-j-1] = div*((w+j)*theta[k-j-2] + (k-j-w)*theta[k-j-1]);
    }
    theta[0] = div*(((T)1) - w)*theta[0];
  }
      
  //--- one more recursion
  T div = ((T)1) / (T)(order-1);
  theta[order-1] = div*w*theta[order-2];
#pragma unroll
  for (int j=1;j <= order-2;j++) {
    theta[order-j-1] = div*((w+j)*theta[order-j-2] + (order-j-w)*theta[order-j-1]);
  }
    
  theta[0] = div*(((T)1) - w)*theta[0];
}

//
// Calculate theta and dtheta for general order bspline
//
template <typename T, typename T3, int order>
__forceinline__ __device__ void calc_theta_dtheta(T wx, T wy, T wz, T3 *theta, T3 *dtheta) {

  theta[order-1].x = ((T)0);
  theta[order-1].y = ((T)0);
  theta[order-1].z = ((T)0);
  theta[1].x = wx;
  theta[1].y = wy;
  theta[1].z = wz;
  theta[0].x = ((T)1) - wx;
  theta[0].y = ((T)1) - wy;
  theta[0].z = ((T)1) - wz;

#pragma unroll
  for (int k=3;k <= order-1;k++) {
    T div = ((T)1) / (T)(k-1);
    theta[k-1].x = div*wx*theta[k-2].x;
    theta[k-1].y = div*wy*theta[k-2].y;
    theta[k-1].z = div*wz*theta[k-2].z;
#pragma unroll
    for (int j=1;j <= k-2;j++) {
      theta[k-j-1].x = div*((wx + j)*theta[k-j-2].x + (k-j-wx)*theta[k-j-1].x);
      theta[k-j-1].y = div*((wy + j)*theta[k-j-2].y + (k-j-wy)*theta[k-j-1].y);
      theta[k-j-1].z = div*((wz + j)*theta[k-j-2].z + (k-j-wz)*theta[k-j-1].z);
    }
    theta[0].x = div*(((T)1) - wx)*theta[0].x;
    theta[0].y = div*(((T)1) - wy)*theta[0].y;
    theta[0].z = div*(((T)1) - wz)*theta[0].z;
  }

  //--- perform standard b-spline differentiation
  dtheta[0].x = -theta[0].x;
  dtheta[0].y = -theta[0].y;
  dtheta[0].z = -theta[0].z;
#pragma unroll
  for (int j=2;j <= order;j++) {
    dtheta[j-1].x = theta[j-2].x - theta[j-1].x;
    dtheta[j-1].y = theta[j-2].y - theta[j-1].y;
    dtheta[j-1].z = theta[j-2].z - theta[j-1].z;
  }
      
  //--- one more recursion
  T div = ((T)1) / (T)(order-1);
  theta[order-1].x = div*wx*theta[order-2].x;
  theta[order-1].y = div*wy*theta[order-2].y;
  theta[order-1].z = div*wz*theta[order-2].z;
#pragma unroll
  for (int j=1;j <= order-2;j++) {
    theta[order-j-1].x = div*((wx + j)*theta[order-j-2].x + (order-j-wx)*theta[order-j-1].x);
    theta[order-j-1].y = div*((wy + j)*theta[order-j-2].y + (order-j-wy)*theta[order-j-1].y);
    theta[order-j-1].z = div*((wz + j)*theta[order-j-2].z + (order-j-wz)*theta[order-j-1].z);
  }
    
  theta[0].x = div*(((T)1) - wx)*theta[0].x;
  theta[0].y = div*(((T)1) - wy)*theta[0].y;
  theta[0].z = div*(((T)1) - wz)*theta[0].z;
}

//
// Spreads the charge on the grid. Calculates theta and dtheta on the fly
// blockDim.x                   = Number of atoms each block loads
// blockDim.y*blockDim.x/order3 = Number of atoms we spread at once
//
template <typename AT, int order, bool periodicY, bool periodicZ>
__global__ void
spread_charge_ortho_kernel(const float4 *xyzq, const int ncoord,
          const float recip11, const float recip22, const float recip33,
          const int nfftx, const int nffty, const int nfftz,
          const int xsize, const int ysize, const int zsize,
          const int xdim, const int y00, const int z00,
          AT* data) {

  // Shared memory use:
  // order = 4: 1920 bytes
  // order = 6: 2688 bytes
  // order = 8: 3456 bytes
  __shared__ int sh_ix[32];
  __shared__ int sh_iy[32];
  __shared__ int sh_iz[32];
  __shared__ float sh_thetax[order*32];
  __shared__ float sh_thetay[order*32];
  __shared__ float sh_thetaz[order*32];

  // Process atoms pos to pos_end-1 (blockDim.x)
  const unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int pos_end = min((blockIdx.x+1)*blockDim.x, ncoord);

  if (pos < pos_end && threadIdx.y == 0) {

    float4 xyzqi = xyzq[pos];
    float x = xyzqi.x;
    float y = xyzqi.y;
    float z = xyzqi.z;
    float q = xyzqi.w;

    // float w;
    // w = x*recip11 + 2.0f;
    // float frx = (float)(nfftx*(w - (floorf(w + 0.5f) - 0.5f)));
    // w = y*recip22 + 2.0f;
    // float fry = (float)(nffty*(w - (floorf(w + 0.5f) - 0.5f)));
    // w = z*recip33 + 2.0f;
    // float frz = (float)(nfftz*(w - (floorf(w + 0.5f) - 0.5f)));

    float frx = ((float)nfftx)*x;
    float fry = ((float)nffty)*y;
    float frz = ((float)nfftz)*z;

    int frxi = (int)frx;
    int fryi = (int)fry;
    int frzi = (int)frz;

    float wx = frx - (float)frxi;
    float wy = fry - (float)fryi;
    float wz = frz - (float)frzi;

    if (!periodicY && y00 == 0 && fryi >= ysize) fryi -= nffty;
    if (!periodicZ && z00 == 0 && frzi >= zsize) frzi -= nfftz;

    sh_ix[threadIdx.x] = frxi;
    sh_iy[threadIdx.x] = fryi - y00;
    sh_iz[threadIdx.x] = frzi - z00;

    float theta[order];

    calc_one_theta<float, order>(wx, theta);
#pragma unroll
    for (int i=0;i < order;i++) sh_thetax[threadIdx.x*order + i] = q*theta[i];

    calc_one_theta<float, order>(wy, theta);
#pragma unroll
    for (int i=0;i < order;i++) sh_thetay[threadIdx.x*order + i] = theta[i];

    calc_one_theta<float, order>(wz, theta);
#pragma unroll
    for (int i=0;i < order;i++) sh_thetaz[threadIdx.x*order + i] = theta[i];

  }

  __syncthreads();

  // Grid point location, values of (ix0, iy0, iz0) are in range 0..order-1
  // NOTE: Only tid=0...order*order*order-1 do any computation
  const int order3 = ((order*order*order-1)/warpSize + 1)*warpSize;
  const int tid = (threadIdx.x + threadIdx.y*blockDim.x) % order3;   // 0...order3-1
  const int x0 = tid % order;
  const int y0 = (tid / order) % order;
  const int z0 = tid / (order*order);

  // Loop over atoms pos..pos_end-1
  int iadd = blockDim.x*blockDim.y/order3;
  int i = (threadIdx.x + threadIdx.y*blockDim.x)/order3;
  int iend = pos_end - blockIdx.x*blockDim.x;
  for (;i < iend;i += iadd) {
    int x = sh_ix[i] + x0;
    int y = sh_iy[i] + y0;
    int z = sh_iz[i] + z0;
      
    if (x >= nfftx) x -= nfftx;

    if (periodicY  && (y >= nffty)) y -= nffty;
    if (!periodicY && (y < 0 || y >= ysize)) continue;

    if (periodicZ  && (z >= nfftz)) z -= nfftz;
    if (!periodicZ && (z < 0 || z >= zsize)) continue;
      
    // Get position on the grid
    int ind = x + xdim*(y + ysize*(z));
    
    // Here we unroll the 6x6x6 loop with 216 threads.
    // NOTE: We use 7*32=224 threads to do this
    // Calculate interpolated charge value and store it to global memory

    if (tid < order*order*order)
      write_grid<AT>(sh_thetax[i*order+x0]*sh_thetay[i*order+y0]*sh_thetaz[i*order+z0], ind, data);

    // if (tid < order*order*order)
    //   write_grid<AT>(1.0f, ind, data);

    // if (tid == 0)
    //   write_grid<AT>(1.0f, ind, data);
  }

}

/*
#define PATCHSPLIT 1

//
// Spreads the charge on the grid. Calculates theta and dtheta on the fly
// blockDim.x                   = Number of atoms each block loads
// blockDim.y*blockDim.x/order3 = Number of atoms we spread at once
//
template <typename AT, int order>
__global__ void
spread_charge_ortho_kernel(const int4 *patch, const int npatch,
  const float4 *xyzq, const int ncoord,
  const float recip11, const float recip22, const float recip33,
  const int nfftx, const int nffty, const int nfftz,
  const int xsize, const int ysize,
  AT* data) {

  // Shared memory use:
  // order = 4: 1920 bytes
  // order = 6: 2688 bytes
  // order = 8: 3456 bytes
  __shared__ int sh_ix[32];
  __shared__ int sh_iy[32];
  __shared__ int sh_iz[32];
  __shared__ float sh_thetax[order*32];
  __shared__ float sh_thetay[order*32];
  __shared__ float sh_thetaz[order*32];

  // blockIdx.x = patch index
  int ipatch = blockIdx.x / PATCHSPLIT;
  int fpatch = blockIdx.x % PATCHSPLIT;
  if (ipatch >= npatch) return;

  int4 patchval = patch[ipatch];
  int pos = (ipatch == 0) ? 0 : patch[ipatch-1].w;
  const int ox = patchval.x;
  const int oy = patchval.y;
  const int oz = patchval.z;
  int pos_end = patchval.w;
  int n = pos_end - pos;
  pos += n*fpatch/PATCHSPLIT;
  pos_end = pos + n*(fpatch+1)/PATCHSPLIT;

  // This block takes care of atoms in the range [pos, pos_end-1]
  for (;pos < pos_end;pos += blockDim.x) {

    __syncthreads();

    if (pos + threadIdx.x < pos_end && threadIdx.y == 0) {

      float4 xyzqi = xyzq[pos + threadIdx.x];
      float x = xyzqi.x;
      float y = xyzqi.y;
      float z = xyzqi.z;
      float q = xyzqi.w;

      float w;

      w = x*recip11 + 2.0f;
      float frx = (float)(nfftx*(w - (floorf(w + 0.5f) - 0.5f)));
      w = y*recip22 + 2.0f;
      float fry = (float)(nffty*(w - (floorf(w + 0.5f) - 0.5f)));
      w = z*recip33 + 2.0f;
      float frz = (float)(nfftz*(w - (floorf(w + 0.5f) - 0.5f)));

      int frxi = (int)frx;
      int fryi = (int)fry;
      int frzi = (int)frz;

      sh_ix[threadIdx.x] = frxi + ox;
      sh_iy[threadIdx.x] = fryi + oy;
      sh_iz[threadIdx.x] = frzi + oz;

      float wx = frx - (float)frxi;
      float wy = fry - (float)fryi;
      float wz = frz - (float)frzi;

      float theta[order];

      calc_one_theta<float, order>(wx, theta);
  #pragma unroll
      for (int i=0;i < order;i++) sh_thetax[threadIdx.x*order + i] = q*theta[i];

      calc_one_theta<float, order>(wy, theta);
  #pragma unroll
      for (int i=0;i < order;i++) sh_thetay[threadIdx.x*order + i] = theta[i];

      calc_one_theta<float, order>(wz, theta);
  #pragma unroll
      for (int i=0;i < order;i++) sh_thetaz[threadIdx.x*order + i] = theta[i];

    }

    __syncthreads();

    // Grid point location, values of (ix0, iy0, iz0) are in range 0..order-1
    // NOTE: Only tid=0...order*order*order-1 do any computation
    const int order3 = ((order*order*order-1)/warpSize + 1)*warpSize;
    const int tid = (threadIdx.x + threadIdx.y*blockDim.x) % order3;   // 0...order3-1
    const int x0 = tid % order;
    const int y0 = (tid / order) % order;
    const int z0 = tid / (order*order);

    // Loop over atoms
    // iadd = number of atoms to compute per iteration
    int iadd = blockDim.x*blockDim.y/order3;
    int i = (threadIdx.x + threadIdx.y*blockDim.x)/order3;
    int iend = min(blockDim.x, pos_end-pos);
    for (;i < iend;i += iadd) {
      int x = sh_ix[i] + x0;
      int y = sh_iy[i] + y0;
      int z = sh_iz[i] + z0;
        
      if (x >= nfftx) x -= nfftx;
      if (y >= nffty) y -= nffty;
      if (z >= nfftz) z -= nfftz;
        
      // Get position on the grid
      int ind = x + xsize*(y + ysize*(z));
        
      // Here we unroll the 6x6x6 loop with 216 threads.
      // NOTE: We use 7*32=224 threads to do this
      // Calculate interpolated charge value and store it to global memory
      if (tid < order*order*order)
        write_grid<AT>(sh_thetax[i*order+x0]*sh_thetay[i*order+y0]*sh_thetaz[i*order+z0], ind, data);
    }
  }
}
*/

//-----------------------------------------------------------------------------------------
// Generic version can not be used
template <typename T> __forceinline__ __device__
void gather_force_store(const float fx, const float fy, const float fz,
      const int stride, const int pos,
      T* force) {
}

// Template specialization for "float"
template <> __forceinline__ __device__
void gather_force_store<float>(const float fx, const float fy, const float fz, 
             const int stride, const int pos, 
             float* force) {
  // Store into non-strided float XYZ array
  force[pos]          = fx;
  force[pos+stride]   = fy;
  force[pos+stride*2] = fz;
}

// Template specialization for "float3"
template <> __forceinline__ __device__
void gather_force_store<float3>(const float fx, const float fy, const float fz, 
        const int stride, const int pos, 
        float3* force) {
  // Store into non-strided "float3" array
  force[pos].x = fx;
  force[pos].y = fy;
  force[pos].z = fz;
}
//-----------------------------------------------------------------------------------------

// Per atom data structure for the gather_force -kernels
template <typename T, int order>
struct gather_t {
  int ix;
  int iy;
  int iz;
  T charge;
  T thetax[order];
  T thetay[order];
  T thetaz[order];
  T dthetax[order];
  T dthetay[order];
  T dthetaz[order];
  float f1;
  float f2;
  float f3;
};

//
// Gathers forces from the grid
// blockDim.x            = Number of atoms each block loads
// blockDim.x*blockDim.y = Total number of threads per block
//
template <typename CT, typename FT, int order, bool periodicY, bool periodicZ>
__global__ void gather_force_ortho(const float4 *xyzq, const int ncoord,
              const int nfftx, const int nffty, const int nfftz,
              const int xsize, const int ysize, const int zsize,
              const int xdim, const int y00, const int z00,
              const float recip1, const float recip2, const float recip3,
              const float* data,      // NOTE: data is used for loads when __CUDA_ARCH__ >= 350
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
              const cudaTextureObject_t gridTexObj,
#endif
              const int stride,
              FT *force) {

  const int tid = threadIdx.x + threadIdx.y*blockDim.x; // 0...63

  // Shared memory
  __shared__ gather_t<CT, order> shmem[32];
#if __CUDA_ARCH__ < 300
  __shared__ float3 shred_buf[32*2];
  volatile float3 *shred = &shred_buf[(tid/8)*8];
#endif

  const int pos = blockIdx.x*blockDim.x + threadIdx.x;
  const int pos_end = min((blockIdx.x+1)*blockDim.x, ncoord);

  // Load atom data into shared memory
  if (pos < pos_end && threadIdx.y == 0) {

    float4 xyzqi = xyzq[pos];
    float x = xyzqi.x;
    float y = xyzqi.y;
    float z = xyzqi.z;
    float q = xyzqi.w;

    // float w;
    // w = x*recip1 + 2.0f;
    // float frx = (float)(nfftx*(w - (floorf(w + 0.5f) - 0.5f)));
    // w = y*recip2 + 2.0f;
    // float fry = (float)(nffty*(w - (floorf(w + 0.5f) - 0.5f)));
    // w = z*recip3 + 2.0f;
    // float frz = (float)(nfftz*(w - (floorf(w + 0.5f) - 0.5f)));

    float frx = ((float)nfftx)*x;
    float fry = ((float)nffty)*y;
    float frz = ((float)nfftz)*z;

    int frxi = (int)frx;
    int fryi = (int)fry;
    int frzi = (int)frz;

    float wx = frx - (float)frxi;
    float wy = fry - (float)fryi;
    float wz = frz - (float)frzi;

    if (!periodicY && y00 == 0 && fryi >= ysize) fryi -= nffty;
    if (!periodicZ && z00 == 0 && frzi >= zsize) frzi -= nfftz;

    shmem[threadIdx.x].ix = frxi;
    shmem[threadIdx.x].iy = fryi - y00;
    shmem[threadIdx.x].iz = frzi - z00;
    shmem[threadIdx.x].charge = q;

    float3 theta_tmp[order];
    float3 dtheta_tmp[order];
    calc_theta_dtheta<float, float3, order>(wx, wy, wz, theta_tmp, dtheta_tmp);
    
#pragma unroll
    for (int i=0;i < order;i++) shmem[threadIdx.x].thetax[i] = theta_tmp[i].x;

#pragma unroll
    for (int i=0;i < order;i++) shmem[threadIdx.x].thetay[i] = theta_tmp[i].y;

#pragma unroll
    for (int i=0;i < order;i++) shmem[threadIdx.x].thetaz[i] = theta_tmp[i].z;

#pragma unroll
    for (int i=0;i < order;i++) shmem[threadIdx.x].dthetax[i] = dtheta_tmp[i].x;

#pragma unroll
    for (int i=0;i < order;i++) shmem[threadIdx.x].dthetay[i] = dtheta_tmp[i].y;

#pragma unroll
    for (int i=0;i < order;i++) shmem[threadIdx.x].dthetaz[i] = dtheta_tmp[i].z;

  }
  __syncthreads();

  // We divide the order x order x order cube into 8 sub-cubes.
  // These sub-cubes are taken care by a single thread
  // The size of the sub-cubes is:
  // order=4 : 2x2x2
  // order=6 : 3x3x3
  // order=8 : 4x4x4
  const int nsc = (order == 4) ? 2 : ((order == 6) ? 3 : 4);
  // Calculate the starting index on the sub-cube for this thread
  // tid = 0...63
  const int t = (tid % 8);         // sub-cube index (0...7)
  // t = (tx0 + ty0*2 + tz0*4)/nsc
  // (tx0, ty0, tz0) gives the starting index of the 3x3x3 sub-cube
  const int tz0 = (t / 4)*nsc;
  const int ty0 = ((t / 2) % 2)*nsc;
  const int tx0 = (t % 2)*nsc;

  //
  // Calculate forces for 32 atoms. We have 32*2 = 64 threads
  // Loop is iterated 4 times:
  //                         (iterations)
  // Threads 0...7   = atoms 0, 8,  16, 24
  // Threads 8...15  = atoms 1, 9,  17, 25
  // Threads 16...31 = atoms 2, 10, 18, 26
  //                ...
  // Threads 56...63 = atoms 7, 15, 23, 31
  //

  int base = tid/8;
  const int base_end = pos_end - blockIdx.x*blockDim.x;
  while (base < base_end) {

    float f1 = 0.0f;
    float f2 = 0.0f;
    float f3 = 0.0f;
    int ix0 = shmem[base].ix;
    int iy0 = shmem[base].iy;
    int iz0 = shmem[base].iz;

    // Each thread calculates a nsc x nsc x nsc sub-cube
#pragma unroll
    for (int i=0;i < nsc*nsc*nsc;i++) {
      int tz = tz0 + (i/(nsc*nsc));
      int ty = ty0 + ((i/nsc) % nsc);
      int tx = tx0 + (i % nsc);

      int ix = ix0 + tx;
      int iy = iy0 + ty;
      int iz = iz0 + tz;
      if (ix >= nfftx) ix -= nfftx;

      if (periodicY  && (iy >= nffty)) iy -= nffty;
      if (!periodicY && (iy < 0 || iy >= ysize)) continue;

      if (periodicZ  && (iz >= nfftz)) iz -= nfftz;
      if (!periodicZ && (iz < 0 || iz >= zsize)) continue;

      int ind = ix + (iy + iz*ysize)*xdim;

#if __CUDA_ARCH__ >= 350
      float q0 = __ldg(&data[ind]);
#else
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
      float q0 = tex1Dfetch<float>(gridTexObj, ind);
#else
      float q0 = tex1Dfetch(gridTexRef, ind);
#endif
#endif
      float thx0 = shmem[base].thetax[tx];
      float thy0 = shmem[base].thetay[ty];
      float thz0 = shmem[base].thetaz[tz];
      float dthx0 = shmem[base].dthetax[tx];
      float dthy0 = shmem[base].dthetay[ty];
      float dthz0 = shmem[base].dthetaz[tz];
      f1 += dthx0 * thy0 * thz0 * q0;
      f2 += thx0 * dthy0 * thz0 * q0;
      f3 += thx0 * thy0 * dthz0 * q0;
    }

    //-------------------------

    // Reduce
#if __CUDA_ARCH__ >= 300
    const int i = threadIdx.x & 7;

    f1 += __shfl(f1, i+4, 8);
    f2 += __shfl(f2, i+4, 8);
    f3 += __shfl(f3, i+4, 8);

    f1 += __shfl(f1, i+2, 8);
    f2 += __shfl(f2, i+2, 8);
    f3 += __shfl(f3, i+2, 8);

    f1 += __shfl(f1, i+1, 8);
    f2 += __shfl(f2, i+1, 8);
    f3 += __shfl(f3, i+1, 8);

    if (i == 0) {
      shmem[base].f1 = f1;
      shmem[base].f2 = f2;
      shmem[base].f3 = f3;
    }

#else
    const int i = threadIdx.x & 7;
    shred[i].x = f1;
    shred[i].y = f2;
    shred[i].z = f3;

    if (i < 4) {
      shred[i].x += shred[i+4].x;
      shred[i].y += shred[i+4].y;
      shred[i].z += shred[i+4].z;
    }

    if (i < 2) {
      shred[i].x += shred[i+2].x;
      shred[i].y += shred[i+2].y;
      shred[i].z += shred[i+2].z;
    }

    if (i == 0) {
      shmem[base].f1 = shred[0].x + shred[1].x;
      shmem[base].f2 = shred[0].y + shred[1].y;
      shmem[base].f3 = shred[0].z + shred[1].z;
    }
#endif

    base += 8;
  }

  // Write forces
  __syncthreads();
  if (pos < pos_end && threadIdx.y == 0) {
    float f1 = shmem[threadIdx.x].f1;
    float f2 = shmem[threadIdx.x].f2;
    float f3 = shmem[threadIdx.x].f3;
    float q = -shmem[threadIdx.x].charge; //*ccelec_float;
    float fx = q*recip1*f1*nfftx;
    float fy = q*recip2*f2*nffty;
    float fz = q*recip3*f3*nfftz;
    gather_force_store<FT>(fx, fy, fz, stride, pos, force);
  }

}

/*
//
// Gathers forces from the grid
// blockDim.x            = Number of atoms each block loads
// blockDim.x*blockDim.y = Total number of threads per block
//
template <typename CT, typename FT, int order>
__global__ void gather_force_ortho(const int4* patch, const int npatch,
  const float4 *xyzq, const int ncoord,
  const int nfftx, const int nffty, const int nfftz,
  const int xsize, const int ysize,
  const float recip1, const float recip2, const float recip3,
  const float* data,      // NOTE: data is used for loads when __CUDA_ARCH__ >= 350
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  const cudaTextureObject_t gridTexObj,
#endif
  const int stride,
  FT *force) {

  const int tid = threadIdx.x + threadIdx.y*blockDim.x; // 0...63

  // Shared memory
  __shared__ gather_t<CT, order> shmem[32];
#if __CUDA_ARCH__ < 300
  __shared__ float3 shred_buf[32*2];
  volatile float3 *shred = &shred_buf[(tid/8)*8];
#endif

  // blockIdx.x = patch index
  if (blockIdx.x >= npatch) return;

  int4 patchval = patch[blockIdx.x];
  int pos = (blockIdx.x == 0) ? 0 : patch[blockIdx.x-1].w;
  const int ox = patchval.x;
  const int oy = patchval.y;
  const int oz = patchval.z;
  const int pos_end = patchval.w;

  // This block takes care of atoms in the range [pos, pos_end-1]
  for (;pos < pos_end;pos += blockDim.x) {

    __syncthreads();

    if (pos + threadIdx.x < pos_end && threadIdx.y == 0) {

      float4 xyzqi = xyzq[pos];
      float x = xyzqi.x;
      float y = xyzqi.y;
      float z = xyzqi.z;
      float q = xyzqi.w;

      float w;

      w = x*recip1 + 2.0f;
      float frx = (float)(nfftx*(w - (floorf(w + 0.5f) - 0.5f)));

      w = y*recip2 + 2.0f;
      float fry = (float)(nffty*(w - (floorf(w + 0.5f) - 0.5f)));

      w = z*recip3 + 2.0f;
      float frz = (float)(nfftz*(w - (floorf(w + 0.5f) - 0.5f)));

      int frxi = (int)frx;
      int fryi = (int)fry;
      int frzi = (int)frz;

      shmem[threadIdx.x].ix = frxi + ox;
      shmem[threadIdx.x].iy = fryi + oy;
      shmem[threadIdx.x].iz = frzi + oz;
      shmem[threadIdx.x].charge = q;

      float wx = frx - (float)frxi;
      float wy = fry - (float)fryi;
      float wz = frz - (float)frzi;

      float3 theta_tmp[order];
      float3 dtheta_tmp[order];
      calc_theta_dtheta<float, float3, order>(wx, wy, wz, theta_tmp, dtheta_tmp);
      
#pragma unroll
      for (int i=0;i < order;i++) shmem[threadIdx.x].thetax[i] = theta_tmp[i].x;

#pragma unroll
      for (int i=0;i < order;i++) shmem[threadIdx.x].thetay[i] = theta_tmp[i].y;

#pragma unroll
      for (int i=0;i < order;i++) shmem[threadIdx.x].thetaz[i] = theta_tmp[i].z;

#pragma unroll
      for (int i=0;i < order;i++) shmem[threadIdx.x].dthetax[i] = dtheta_tmp[i].x;

#pragma unroll
      for (int i=0;i < order;i++) shmem[threadIdx.x].dthetay[i] = dtheta_tmp[i].y;

#pragma unroll
      for (int i=0;i < order;i++) shmem[threadIdx.x].dthetaz[i] = dtheta_tmp[i].z;

    }
    __syncthreads();

    // We divide the order x order x order cube into 8 sub-cubes.
    // These sub-cubes are taken care by a single thread
    // The size of the sub-cubes is:
    // order=4 : 2x2x2
    // order=6 : 3x3x3
    // order=8 : 4x4x4
    const int nsc = (order == 4) ? 2 : ((order == 6) ? 3 : 4);
    // Calculate the starting index on the sub-cube for this thread
    // tid = 0...63
    const int t = (tid % 8);         // sub-cube index (0...7)
    // t = (tx0 + ty0*2 + tz0*4)/nsc
    // (tx0, ty0, tz0) gives the starting index of the 3x3x3 sub-cube
    const int tz0 = (t / 4)*nsc;
    const int ty0 = ((t / 2) % 2)*nsc;
    const int tx0 = (t % 2)*nsc;

    //
    // Calculate forces for 32 atoms. We have 32*2 = 64 threads
    // Loop is iterated 4 times:
    //                         (iterations)
    // Threads 0...7   = atoms 0, 8,  16, 24
    // Threads 8...15  = atoms 1, 9,  17, 25
    // Threads 16...31 = atoms 2, 10, 18, 26
    //                ...
    // Threads 56...63 = atoms 7, 15, 23, 31
    //

    int base = tid/8;
    //const int base_end = pos_end - blockIdx.x*blockDim.x;
    const int base_end = min(blockDim.x, pos_end-pos);
    while (base < base_end) {

      float f1 = 0.0f;
      float f2 = 0.0f;
      float f3 = 0.0f;
      int ix0 = shmem[base].ix;
      int iy0 = shmem[base].iy;
      int iz0 = shmem[base].iz;

      // Each thread calculates a nsc x nsc x nsc sub-cube
#pragma unroll
      for (int i=0;i < nsc*nsc*nsc;i++) {
        int tz = tz0 + (i/(nsc*nsc));
        int ty = ty0 + ((i/nsc) % nsc);
        int tx = tx0 + (i % nsc);

        int ix = ix0 + tx;
        int iy = iy0 + ty;
        int iz = iz0 + tz;
        if (ix >= nfftx) ix -= nfftx;
        if (iy >= nffty) iy -= nffty;
        if (iz >= nfftz) iz -= nfftz;
#if __CUDA_ARCH__ >= 350
        float q0 = __ldg(&data[ix + (iy + iz*ysize)*xsize]);
#else
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        float q0 = tex1Dfetch<float>(gridTexObj, ix + (iy + iz*ysize)*xsize);
#else
        float q0 = tex1Dfetch(gridTexRef, ix + (iy + iz*ysize)*xsize);
#endif
#endif
        float thx0 = shmem[base].thetax[tx];
        float thy0 = shmem[base].thetay[ty];
        float thz0 = shmem[base].thetaz[tz];
        float dthx0 = shmem[base].dthetax[tx];
        float dthy0 = shmem[base].dthetay[ty];
        float dthz0 = shmem[base].dthetaz[tz];
        f1 += dthx0 * thy0 * thz0 * q0;
        f2 += thx0 * dthy0 * thz0 * q0;
        f3 += thx0 * thy0 * dthz0 * q0;
      }

      //-------------------------

      // Reduce
#if __CUDA_ARCH__ >= 300
      const int i = threadIdx.x & 7;

      f1 += __shfl(f1, i+4, 8);
      f2 += __shfl(f2, i+4, 8);
      f3 += __shfl(f3, i+4, 8);

      f1 += __shfl(f1, i+2, 8);
      f2 += __shfl(f2, i+2, 8);
      f3 += __shfl(f3, i+2, 8);

      f1 += __shfl(f1, i+1, 8);
      f2 += __shfl(f2, i+1, 8);
      f3 += __shfl(f3, i+1, 8);

      if (i == 0) {
        shmem[base].f1 = f1;
        shmem[base].f2 = f2;
        shmem[base].f3 = f3;
      }

#else
      const int i = threadIdx.x & 7;
      shred[i].x = f1;
      shred[i].y = f2;
      shred[i].z = f3;

      if (i < 4) {
        shred[i].x += shred[i+4].x;
        shred[i].y += shred[i+4].y;
        shred[i].z += shred[i+4].z;
      }

      if (i < 2) {
        shred[i].x += shred[i+2].x;
        shred[i].y += shred[i+2].y;
        shred[i].z += shred[i+2].z;
      }

      if (i == 0) {
        shmem[base].f1 = shred[0].x + shred[1].x;
        shmem[base].f2 = shred[0].y + shred[1].y;
        shmem[base].f3 = shred[0].z + shred[1].z;
      }
#endif

      base += 8;
    }

    // Write forces
    __syncthreads();
    if (pos < pos_end && threadIdx.y == 0) {
      float f1 = shmem[threadIdx.x].f1;
      float f2 = shmem[threadIdx.x].f2;
      float f3 = shmem[threadIdx.x].f3;
      float q = -shmem[threadIdx.x].charge;//*ccelec_float;
      float fx = q*recip1*f1*nfftx;
      float fy = q*recip2*f2*nffty;
      float fz = q*recip3*f3*nfftz;
      gather_force_store<FT>(fx, fy, fz, stride, pos, force);
    }
  }
}
*/

/*
//
// Calculates sum of squared charge. Used in calculation of self energy
//
__global__ void calc_sum_charge_squared_kernel(const int ncoord, const float4* xyzq,
          double* __restrict__ sum_charge_squared) {
  // Shared memory
  // Required space: blockDim.x*sizeof(double)
  extern __shared__ double sh_q2[];

  int i = threadIdx.x + blockIdx.x*blockDim.x;
  float q = 0.0f;
  if (i < ncoord) q = xyzq[i].w;
  sh_q2[threadIdx.x] = q*q;
  __syncthreads();
  for(int d=1;d < blockDim.x;d *= 2) {
    int t = threadIdx.x + d;
    double q2_val = (t < blockDim.x) ? sh_q2[t] : 0.0;
    __syncthreads();
    sh_q2[threadIdx.x] += q2_val;
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    atomicAdd(sum_charge_squared, sh_q2[0]);
  }

}
*/

const int TILEDIM = 32;
const int TILEROWS = 8;

template <typename T>
__device__ __forceinline__
void transpose_xyz_yzx_device(
  const int x_in, const int y_in, const int z_in,
  const int x_out, const int y_out,
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int ysize_out, const int zsize_out,
  const T* data_in, T* data_out) {

  // Shared memory
  __shared__ T tile[TILEDIM][TILEDIM+1];

  // Read (x,y) data_in into tile (shared memory)
  for (int j=0;j < TILEDIM;j += TILEROWS)
    if ((x_in < nx) && (y_in + j < ny) && (z_in < nz))
      tile[threadIdx.y + j][threadIdx.x] = data_in[x_in + (y_in + j + z_in*ysize_in)*xsize_in];

  __syncthreads();

  // Write (y,x) tile into data_out
  const int z_out = z_in;
  for (int j=0;j < TILEDIM;j += TILEROWS)
    if ((x_out + j < nx) && (y_out < ny) && (z_out < nz))
      data_out[y_out + (z_out + (x_out+j)*zsize_out)*ysize_out] = tile[threadIdx.x][threadIdx.y + j];
}

//
// Transposes a 3d matrix out-of-place: data_in(x, y, z) -> data_out(y, z, x)
//
template <typename T>
__global__ void transpose_xyz_yzx_kernel(
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int ysize_out, const int zsize_out,
  const T* data_in, T* data_out) {

  int x_in = blockIdx.x * TILEDIM + threadIdx.x;
  int y_in = blockIdx.y * TILEDIM + threadIdx.y;
  int z_in = blockIdx.z           + threadIdx.z;

  int x_out = blockIdx.x * TILEDIM + threadIdx.y;
  int y_out = blockIdx.y * TILEDIM + threadIdx.x;

  transpose_xyz_yzx_device<T>(
    x_in, y_in, z_in,
    x_out, y_out,
    nx, ny, nz,
    xsize_in, ysize_in,
    ysize_out, zsize_out,
    data_in, data_out);

/*
  // Shared memory
  __shared__ T tile[TILEDIM][TILEDIM+1];

  int x = blockIdx.x * TILEDIM + threadIdx.x;
  int y = blockIdx.y * TILEDIM + threadIdx.y;
  int z = blockIdx.z           + threadIdx.z;

  // Read (x,y) data_in into tile (shared memory)
  for (int j=0;j < TILEDIM;j += TILEROWS)
    if ((x < nx) && (y + j < ny) && (z < nz))
      tile[threadIdx.y + j][threadIdx.x] = data_in[x + (y + j + z*ysize_in)*xsize_in];

  __syncthreads();

  // Write (y,x) tile into data_out
  x = blockIdx.x * TILEDIM + threadIdx.y;
  y = blockIdx.y * TILEDIM + threadIdx.x;
  for (int j=0;j < TILEDIM;j += TILEROWS)
    if ((x + j < nx) && (y < ny) && (z < nz))
      data_out[y + (z + (x+j)*zsize_out)*ysize_out] = tile[threadIdx.x][threadIdx.y + j];
*/
}

//
// Transposes a batch of 3d matrices out-of-place: data_in(x, y, z) -> data_out(y, z, x)
// Batch index bi is encoded in blockIdx.z, where 
// blockIdx.z = 0...nz-1 are for batch 1
// blockIdx.z = nz...2*nz-1 are for batch 2
// ...
// gridDim.z = nz*numBatches
//
template <typename T>
__global__ void batchTranspose_xyz_yzx_kernel(
  const TransposeBatch<T>* batches,
  const int ny, const int nz, 
  const int xsize_in, const int ysize_in) {

  int x_in = blockIdx.x * TILEDIM + threadIdx.x;
  int y_in = blockIdx.y * TILEDIM + threadIdx.y;
  int z_in = (blockIdx.z % nz)    + threadIdx.z;

  int x_out = blockIdx.x * TILEDIM + threadIdx.y;
  int y_out = blockIdx.y * TILEDIM + threadIdx.x;

  int bi = blockIdx.z/nz;

  TransposeBatch<T> batch = batches[bi];
  int nx        = batch.nx;
  int ysize_out = batch.ysize_out;
  int zsize_out = batch.zsize_out;
  T* data_in    = batch.data_in;
  T* data_out   = batch.data_out;

  transpose_xyz_yzx_device<T>(
    x_in, y_in, z_in,
    x_out, y_out,
    nx, ny, nz,
    xsize_in, ysize_in,
    ysize_out, zsize_out,
    data_in, data_out);

}

/*
//
// Transposes a 3d matrix out-of-place: data_in(x, y, z) -> data_out(y, z, x)
//
template <typename T>
__forceinline__ __device__
void transpose_xyz_yzx_dev(
  const int blockz,
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int xsize_out, const int ysize_out,
  const T* data_in, T* data_out) {

  // Shared memory
  __shared__ T tile[TILEDIM][TILEDIM+1];

  int x = blockIdx.x * TILEDIM + threadIdx.x;
  int y = blockIdx.y * TILEDIM + threadIdx.y;
  // int z = blockIdx.z           + threadIdx.z;
  int z = blockz               + threadIdx.z;

  // Read (x,y) data_in into tile (shared memory)
  for (int j=0;j < TILEDIM;j += TILEROWS)
    if ((x < nx) && (y + j < ny) && (z < nz))
      tile[threadIdx.y + j][threadIdx.x] = data_in[x + (y + j + z*ysize_in)*xsize_in];

  __syncthreads();

  // Write (y,x) tile into data_out
  x = blockIdx.x * TILEDIM + threadIdx.y;
  y = blockIdx.y * TILEDIM + threadIdx.x;
  for (int j=0;j < TILEDIM;j += TILEROWS)
    if ((x + j < nx) && (y < ny) && (z < nz))
      data_out[y + (z + (x+j)*ysize_out)*xsize_out] = tile[threadIdx.x][threadIdx.y + j];

}

//
// Transposes a 3d matrix out-of-place: data_in(x, y, z) -> data_out(y, z, x)
// (nx, ny, nz)                     = size of the transposed volume
// (xsize_in, ysize_in, zsize_in)   = size of the input data
// into nblock memory blocks
//
template <typename T>
__global__ void transpose_xyz_yzx_kernel(
  const int nblock,
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int xsize_out, const int ysize_out,
  const T* data_in, T* data_out) {

  const int iblock = blockIdx.z/nz;

  if (iblock < nblock) {
    transpose_xyz_yzx_dev(blockIdx.z % nz, nx, ny, nz,
      xsize_in, ysize_in, xsize_out, ysize_out,
      data_in, data_out);
  }

}
*/

template <typename T>
__device__ __forceinline__
void transpose_xyz_zxy_device(
  const int x_in, const int y_in, const int z_in,
  const int x_out, const int z_out,
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int zsize_out, const int xsize_out,
  const T* data_in, T* data_out) {

  // Shared memory
  __shared__ T tile[TILEDIM][TILEDIM+1];

  // Read (x,z) data_in into tile (shared memory)
  for (int k=0;k < TILEDIM;k += TILEROWS)
    if ((x_in < nx) && (y_in < ny) && (z_in + k < nz))
      tile[threadIdx.y + k][threadIdx.x] = data_in[x_in + (y_in + (z_in + k)*ysize_in)*xsize_in];

  __syncthreads();

  // Write (z,x) tile into data_out
  const int y_out = y_in;
  for (int k=0;k < TILEDIM;k += TILEROWS)
    if ((x_out + k < nx) && (y_out < ny) && (z_out < nz))
      data_out[z_out + (x_out + k + y_out*xsize_out)*zsize_out] = tile[threadIdx.x][threadIdx.y + k];
}

//
// Transposes a 3d matrix out-of-place: data_in(x, y, z) -> data_out(z, x, y)
//
template <typename T>
__global__ void transpose_xyz_zxy_kernel(
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int zsize_out, const int xsize_out,
  const T* data_in, T* data_out) {

  int x_in = blockIdx.x * TILEDIM + threadIdx.x;
  int y_in = blockIdx.z           + threadIdx.z;
  int z_in = blockIdx.y * TILEDIM + threadIdx.y;

  int x_out = blockIdx.x * TILEDIM + threadIdx.y;
  int z_out = blockIdx.y * TILEDIM + threadIdx.x;

  transpose_xyz_zxy_device<T>(
    x_in, y_in, z_in, x_out, z_out,
    nx, ny, nz,
    xsize_in, ysize_in,
    zsize_out, xsize_out,
    data_in, data_out);

}

//
// Transposes a batch of 3d matrices out-of-place: data_in(x, y, z) -> data_out(z, x, y)
// Batch index bi is encoded in blockIdx.z, where 
// blockIdx.z = 0...ny-1 are for batch 1
// blockIdx.z = ny...2*ny-1 are for batch 2
// ...
// gridDim.z = ny*numBatches
//
template <typename T>
__global__ void batchTranspose_xyz_zxy_kernel(
  const TransposeBatch<T>* batches,
  const int ny, const int nz, 
  const int xsize_in, const int ysize_in) {

  int x_in = blockIdx.x * TILEDIM + threadIdx.x;
  int y_in = (blockIdx.z % ny)    + threadIdx.z;
  int z_in = blockIdx.y * TILEDIM + threadIdx.y;

  int x_out = blockIdx.x * TILEDIM + threadIdx.y;
  int z_out = blockIdx.y * TILEDIM + threadIdx.x;

  int bi = blockIdx.z/ny;

  TransposeBatch<T> batch = batches[bi];
  int nx        = batch.nx;
  int zsize_out = batch.zsize_out;
  int xsize_out = batch.xsize_out;
  T* data_in    = batch.data_in;
  T* data_out   = batch.data_out;

  transpose_xyz_zxy_device<T>(
    x_in, y_in, z_in, x_out, z_out,
    nx, ny, nz,
    xsize_in, ysize_in,
    zsize_out, xsize_out,
    data_in, data_out);

}

//#######################################################################################
//#######################################################################################
//#######################################################################################

void spread_charge_ortho(const float4 *atoms, const int numAtoms,
  const float recip11, const float recip22, const float recip33,
  const int nfftx, const int nffty, const int nfftz,
  const int xsize, const int ysize, const int zsize,
  const int xdim, const int y00, const int z00, 
  const bool periodicY, const bool periodicZ,
  float* data, const int order, cudaStream_t stream) {

  dim3 nthread, nblock;

  switch(order) {
  case 4:
    nthread.x = 32;
    nthread.y = 4;
    nthread.z = 1;
    nblock.x = (numAtoms - 1)/nthread.x + 1;
    nblock.y = 1;
    nblock.z = 1;
    if (periodicY && periodicZ)
      spread_charge_ortho_kernel<float, 4, true, true> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else if (periodicY)
      spread_charge_ortho_kernel<float, 4, true, false> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else if (periodicZ)
      spread_charge_ortho_kernel<float, 4, false, true> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else
      spread_charge_ortho_kernel<float, 4, false, false> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    break;

  case 6:
    nthread.x = 32;
    nthread.y = 7;
    nthread.z = 1;
    nblock.x = (numAtoms - 1)/nthread.x + 1;
    nblock.y = 1;
    nblock.z = 1;
    if (periodicY && periodicZ)
      spread_charge_ortho_kernel<float, 6, true, true> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else if (periodicY)
      spread_charge_ortho_kernel<float, 6, true, false> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else if (periodicZ)
      spread_charge_ortho_kernel<float, 6, false, true> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else
      spread_charge_ortho_kernel<float, 6, false, false> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    break;

  case 8:
    nthread.x = 32;
    nthread.y = 16;
    nthread.z = 1;
    nblock.x = (numAtoms - 1)/nthread.x + 1;
    nblock.y = 1;
    nblock.z = 1;
    if (periodicY && periodicZ)
      spread_charge_ortho_kernel<float, 8, true, true> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else if (periodicY)
      spread_charge_ortho_kernel<float, 8, true, false> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else if (periodicZ)
      spread_charge_ortho_kernel<float, 8, false, true> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    else
      spread_charge_ortho_kernel<float, 8, false, false> <<< nblock, nthread, 0, stream >>>
        (atoms, numAtoms, recip11, recip22, recip33,
         nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00, data);
    break;

  default:
    char str[128];
    sprintf(str, "spread_charge_ortho, order %d not implemented",order);
    cudaNAMD_bug(str);
  }
  cudaCheck(cudaGetLastError());

}

void scalar_sum_ortho(const bool orderXYZ, const int nfft1, const int nfft2, const int nfft3,
  const int size1, const int size2, const int size3, const double kappa,
  const float recip11, const float recip22, const float recip33,
  const float* prefac1, const float* prefac2, const float* prefac3,
  const int k2_00, const int k3_00,
  const bool doEnergyVirial, double* energy, double* virial, float2* data,
  const int cuda_arch, cudaStream_t stream) {

  // Best performance:
  // cuda_arch = 200:
  // energy & virial & (C2075 | K40c) & 512x14: 102.7 (C2075) | 70.4 (K240c)
  // C2075 & 768x12: 27.4
  int nthread = 512;
  int nblock = 10;

  if (cuda_arch < 300) {
    if (doEnergyVirial) {
      nthread = 512;
      nblock = 14;
    } else {
      nthread = 768;
      nblock = 12;
    }
  } else {
    if (doEnergyVirial) {
      nthread = 1024;
      nblock = 14;
    } else {
      nthread = 1024;
      nblock = 14;
    }
  }


  int shmem_size = sizeof(float)*(nfft1 + nfft2 + nfft3);
  if (doEnergyVirial) {
    if (cuda_arch < 300) {
      shmem_size = max(shmem_size, (int)(nthread*sizeof(RecipVirial_t)));
    } else {
      const int warpSize = 32;      
      shmem_size = max(shmem_size, (int)((nthread/warpSize)*sizeof(RecipVirial_t)));
    }
  }

  double inv_vol = recip11*recip22*recip33;
  float piv_inv = (float)(inv_vol/M_PI);
  float fac = (float)(M_PI*M_PI/(kappa*kappa));

  int nf1 = nfft1/2 + (nfft1 % 2);
  int nf2 = nfft2/2 + (nfft2 % 2);
  int nf3 = nfft3/2 + (nfft3 % 2);

  if (doEnergyVirial) {
    if (orderXYZ) {
      scalar_sum_ortho_kernel<float, float2, true, true> <<< nblock, nthread, shmem_size, stream >>>
      (nfft1, nfft2, nfft3, size1, size2, size3,
        nf1, nf2, nf3, recip11, recip22, recip33,
        prefac1, prefac2, prefac3,
        fac, piv_inv, k2_00, k3_00, data, energy, virial);
    } else {
      scalar_sum_ortho_kernel<float, float2, true, false> <<< nblock, nthread, shmem_size, stream >>>
      (nfft1, nfft2, nfft3, size1, size2, size3,
        nf1, nf2, nf3, recip11, recip22, recip33,
        prefac1, prefac2, prefac3,
        fac, piv_inv, k2_00, k3_00, data, energy, virial);
    }
  } else {
    if (orderXYZ) {
      scalar_sum_ortho_kernel<float, float2, false, true> <<< nblock, nthread, shmem_size, stream >>>
      (nfft1, nfft2, nfft3, size1, size2, size3,
        nf1, nf2, nf3, recip11, recip22, recip33,
        prefac1, prefac2, prefac3,
        fac, piv_inv, k2_00, k3_00, data, NULL, NULL);
    } else {
      scalar_sum_ortho_kernel<float, float2, false, false> <<< nblock, nthread, shmem_size, stream >>>
      (nfft1, nfft2, nfft3, size1, size2, size3,
        nf1, nf2, nf3, recip11, recip22, recip33,
        prefac1, prefac2, prefac3,
        fac, piv_inv, k2_00, k3_00, data, NULL, NULL);
    }
  }
  cudaCheck(cudaGetLastError());

}

void gather_force_ortho(const float4 *atoms, const int numAtoms,
  const float recip11, const float recip22, const float recip33,
  const int nfftx, const int nffty, const int nfftz,
  const int xsize, const int ysize, const int zsize,
  const int xdim, const int y00, const int z00, 
  const bool periodicY, const bool periodicZ,
  const float* data, const int order, float3* force, 
#ifndef DISABLE_CUDA_TEXTURE_OBJECTS
  const cudaTextureObject_t gridTexObj,
#endif
  cudaStream_t stream) {

  dim3 nthread(32, 2, 1);
  dim3 nblock((numAtoms - 1)/nthread.x + 1, 1, 1);
  // dim3 nblock(npatch, 1, 1);

  switch(order) {
    case 4:
    if (periodicY && periodicZ)
      gather_force_ortho<float, float3, 4, true, true> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else if (periodicY)
      gather_force_ortho<float, float3, 4, true, false> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else if (periodicZ)
      gather_force_ortho<float, float3, 4, false, true> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else
      gather_force_ortho<float, float3, 4, false, false> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    break;

    case 6:
    if (periodicY && periodicZ)
      gather_force_ortho<float, float3, 6, true, true> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else if (periodicY)
      gather_force_ortho<float, float3, 6, true, false> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else if (periodicZ)
      gather_force_ortho<float, float3, 6, false, true> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else
      gather_force_ortho<float, float3, 6, false, false> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    break;
 
    case 8:
    if (periodicY && periodicZ)
      gather_force_ortho<float, float3, 8, true, true> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else if (periodicY)
      gather_force_ortho<float, float3, 8, true, false> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else if (periodicZ)
      gather_force_ortho<float, float3, 8, false, true> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    else
      gather_force_ortho<float, float3, 8, false, false> <<< nblock, nthread, 0, stream >>>
      (atoms, numAtoms, nfftx, nffty, nfftz, xsize, ysize, zsize, xdim, y00, z00,
        recip11, recip22, recip33, data,
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
        gridTexObj,
  #endif
        1, force);
    break;

    default:
    char str[128];
    sprintf(str, "gather_force_ortho, order %d not implemented",order);
    cudaNAMD_bug(str);
  }
  cudaCheck(cudaGetLastError());

}

#ifdef DISABLE_CUDA_TEXTURE_OBJECTS
void bindGridTexture(float* data, int data_len) {
  gridTexRef.normalized = 0;
  gridTexRef.filterMode = cudaFilterModePoint;
  gridTexRef.addressMode[0] = cudaAddressModeClamp;
  gridTexRef.channelDesc.x = sizeof(float)*8;
  gridTexRef.channelDesc.y = 0;
  gridTexRef.channelDesc.z = 0;
  gridTexRef.channelDesc.w = 0;
  gridTexRef.channelDesc.f = cudaChannelFormatKindFloat;
  cudaCheck(cudaBindTexture(NULL, gridTexRef, data, data_len*sizeof(float)));
}
#endif

/*
void calc_sum_charge_squared(const float4 *atoms, const int numAtoms, double* sum_charge_squared,
  cudaStream_t stream) {
  
  int nthread = 256;
  int nblock = (numAtoms-1)/nthread+1;
  int shmem_size = nthread*sizeof(double);
  calc_sum_charge_squared_kernel<<< nblock, nthread, shmem_size, stream >>>
    (numAtoms, atoms, sum_charge_squared);
  cudaCheck(cudaGetLastError());

}
*/

//
// Transpose
//
void transpose_xyz_yzx(
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int ysize_out, const int zsize_out,
  const float2* data_in, float2* data_out, cudaStream_t stream) {

  dim3 numthread(TILEDIM, TILEROWS, 1);
  dim3 numblock((nx-1)/TILEDIM+1, (ny-1)/TILEDIM+1, nz);

  transpose_xyz_yzx_kernel<float2> <<< numblock, numthread, 0, stream >>>
  (nx, ny, nz, xsize_in, ysize_in,
    ysize_out, zsize_out,
    data_in, data_out);

  cudaCheck(cudaGetLastError());
}

//
// Batched transpose
//
void batchTranspose_xyz_yzx(
  const int numBatches, TransposeBatch<float2>* batches, 
  const int max_nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in, cudaStream_t stream) {

  dim3 numthread(TILEDIM, TILEROWS, 1);
  dim3 numblock((max_nx-1)/TILEDIM+1, (ny-1)/TILEDIM+1, nz*numBatches);

  batchTranspose_xyz_yzx_kernel<float2> <<< numblock, numthread, 0, stream >>>
  (batches, ny, nz, xsize_in, ysize_in);

  cudaCheck(cudaGetLastError());
}

//
// Transpose
//
void transpose_xyz_zxy(
  const int nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  const int zsize_out, const int xsize_out,
  const float2* data_in, float2* data_out, cudaStream_t stream) {

  dim3 numthread(TILEDIM, TILEROWS, 1);
  dim3 numblock((nx-1)/TILEDIM+1, (nz-1)/TILEDIM+1, ny);

  transpose_xyz_zxy_kernel<float2> <<< numblock, numthread, 0, stream >>>
  (nx, ny, nz, xsize_in, ysize_in,
    zsize_out, xsize_out,
    data_in, data_out);

  cudaCheck(cudaGetLastError());
}

//
// Batched transpose
//
void batchTranspose_xyz_zxy(
  const int numBatches, TransposeBatch<float2>* batches, 
  const int max_nx, const int ny, const int nz,
  const int xsize_in, const int ysize_in,
  cudaStream_t stream) {

  dim3 numthread(TILEDIM, TILEROWS, 1);
  dim3 numblock((max_nx-1)/TILEDIM+1, (nz-1)/TILEDIM+1, ny*numBatches);

  batchTranspose_xyz_zxy_kernel<float2> <<< numblock, numthread, 0, stream >>>
  (batches, ny, nz, xsize_in, ysize_in);

  cudaCheck(cudaGetLastError());
}

#endif // NAMD_CUDA
