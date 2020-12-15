
#include "ComputePmeCUDAKernel.h"

void NAMD_die(const char *);

#ifdef NAMD_CUDA

#include <cuda.h>

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 200
// allow linking, must prevent execution elsewhere
__device__ void atomicAdd(float *, float) { }
#endif

#if 0
#include <stdio.h>
#endif

// must run with at least order**2 threads
#define warps_per_block 8
#define atoms_per_warp 4
#define atoms_per_block (atoms_per_warp * warps_per_block)


static const float order_4_coeffs[4][4] =
{{1, -3, 3, -1}, {0, 3, -6, 3}, {0, 3, 0, -3}, {0, 1, 4, 1}};
// divide by 3! = 6

static const float order_6_coeffs[6][6] =
{{1, -5, 10, -10, 5, -1}, {0, 5, -20, 30, -20, 5}, {0, 10, -20, 0, 
  20, -10}, {0, 10, 20, -60, 20, 10}, {0, 5, 50, 0, -50, -5}, {0, 1, 
  26, 66, 26, 1}};
// divide by 5! = 120

static const float order_8_coeffs[8][8] =
{{1, -7, 21, -35, 35, -21, 7, -1}, {0, 7, -42, 105, -140, 105, -42, 
  7}, {0, 21, -84, 105, 0, -105, 84, -21}, {0, 35, 0, -315, 560, -315,
   0, 35}, {0, 35, 280, -665, 0, 665, -280, -35}, {0, 21, 504, 
  315, -1680, 315, 504, 21}, {0, 7, 392, 1715, 
  0, -1715, -392, -7}, {0, 1, 120, 1191, 2416, 1191, 120, 1}};
// divide by 7! = 5040

void cuda_init_bspline_coeffs(float **c, float **dc, int order) {
  float *coeffs = new float[order*order];
  float *dcoeffs = new float[order*order];
  double divisor;
  static const float *scoeffs;
  switch ( order ) {
  case 4:
    scoeffs = &order_4_coeffs[0][0];
    divisor = 6;
    break;
  case 6:
    scoeffs = &order_6_coeffs[0][0];
    divisor = 120;
    break;
  case 8:
    scoeffs = &order_8_coeffs[0][0];
    divisor = 5040;
    break;
  default:
    NAMD_die("unsupported PMEInterpOrder");
  }
  double sum = 0;
  for ( int i=0, p=order-1; i<order; ++i,--p ) {
    for ( int j=0; j<order; ++j ) {
      double c = scoeffs[i*order+(order-1-j)];  // reverse order
      sum += c;
      c /= divisor;
      coeffs[i*order+j] = c;
      dcoeffs[i*order+j] = (double)p * c;
      // printf("%d %d %f %f\n", i, j, c, (double)p*c);
    }
  }
  // printf("checksum: %f %f\n", sum, divisor);
  if ( sum != divisor )
    NAMD_die("cuda_init_bspline_coeffs static data checksum error");
  cudaMalloc((void**) c, order*order*sizeof(float));
  cudaMalloc((void**) dc, order*order*sizeof(float));
  cudaMemcpy(*c, coeffs, order*order*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(*dc, dcoeffs, order*order*sizeof(float), cudaMemcpyHostToDevice);
  delete [] coeffs;
  delete [] dcoeffs;
}

#define BSPLINE_DEFS \
  __shared__ union { \
    float a2d[order][order]; \
    float a1d[order*order]; \
  } bspline_coeffs; \
  __shared__ volatile union { \
    float a2d[atoms_per_warp][7]; \
    float a1d[atoms_per_warp*7]; \
  } atoms[warps_per_block]; \
  if ( threadIdx.x < order*order ) { \
    bspline_coeffs.a1d[threadIdx.x] = coeffs[threadIdx.x]; \
  } \
  __syncthreads();

// simplify warp-synchronous operations
#define WARP (threadIdx.x>>5)
#define THREAD (threadIdx.x&31)
#define FX bspline_factors[threadIdx.x>>5][0]
#define FY bspline_factors[threadIdx.x>>5][1]
#define FZ bspline_factors[threadIdx.x>>5][2]
#define DFX bspline_dfactors[threadIdx.x>>5][0]
#define DFY bspline_dfactors[threadIdx.x>>5][1]
#define DFZ bspline_dfactors[threadIdx.x>>5][2]

#define FX2 bspline_2factors[threadIdx.x>>5][0]
#define FY2 bspline_2factors[threadIdx.x>>5][1]
#define FZ2 bspline_2factors[threadIdx.x>>5][2]

#define ATOM atoms[threadIdx.x>>5].a2d[i_atom]
#define AX atoms[threadIdx.x>>5].a2d[i_atom][0]
#define AY atoms[threadIdx.x>>5].a2d[i_atom][1]
#define AZ atoms[threadIdx.x>>5].a2d[i_atom][2]
#define AQ atoms[threadIdx.x>>5].a2d[i_atom][3]
#define AI atoms[threadIdx.x>>5].a2d[i_atom][4]
#define AJ atoms[threadIdx.x>>5].a2d[i_atom][5]
#define AK atoms[threadIdx.x>>5].a2d[i_atom][6]



template <int order>
__global__ void cuda_pme_charges_batched_dev(
  const float * __restrict__ coeffs, \
  float * const * __restrict__ q_arr, int * __restrict__ f_arr, int * __restrict__ fz_arr, \
  float ** __restrict__ a_data_ptr, int* n_atoms_ptr, \
  int* K1_ptr, int* K2_ptr, int* K3_ptr
) {

  int patchIndex = blockIdx.y;
  int K1 = K1_ptr[patchIndex];
  int K2 = K2_ptr[patchIndex];
  int K3 = K3_ptr[patchIndex];
  int n_atoms = n_atoms_ptr[patchIndex];

  BSPLINE_DEFS
  __shared__ volatile float bspline_factors[warps_per_block][3][order];
  int atoms_this_warp = atoms_per_warp;
  {
    int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp;
    if ( aidx + atoms_per_warp > n_atoms ) {
      atoms_this_warp = n_atoms - aidx;  // may be negative
      if ( atoms_this_warp < 0 ) atoms_this_warp = 0;
    }
    if ( THREAD < atoms_this_warp*7 ) {
      float* a_data = a_data_ptr[patchIndex];
      atoms[WARP].a1d[THREAD] = a_data[aidx*7+THREAD];
    }
  }
  for ( int i_atom=0; i_atom < atoms_this_warp; ++i_atom ) {
    float aq=AQ;
    int ai=(int)AI;
    int aj=(int)AJ;
    int ak=(int)AK;

    if ( THREAD < 3 * order ) {
      const float af = ATOM[THREAD/order];  // x, y, or z
      float f = bspline_coeffs.a2d[0][THREAD % order];
      for ( int i=1; i<order; ++i ) {
        f = af * f + bspline_coeffs.a2d[i][THREAD % order];
      }
      bspline_factors[WARP][THREAD/order][THREAD % order] = f;
    }
    for ( int i=THREAD; i < order*order; i+=32 ) {
      int ti = i/order;
      int tj = i%order;
      int gti = ti + ai;  if ( gti >= K1 ) gti -= K1;
      int gtj = tj + aj;  if ( gtj >= K2 ) gtj -= K2;
      f_arr[gti * K2 + gtj] = 1;
    }
    if ( THREAD < order ) {
      int gtk = ak;  if ( gtk >= K3 ) gtk -= K3;
      gtk += THREAD;  // padded to increase coalescing (maybe, but reduces re-use)

      fz_arr[gtk] = 1;
    }

    for ( int i=THREAD; i < order*order*order; i+=32 ) {
      int ti = i/(order*order);
      int tj = (i/order)%order;
      int tk = i%order;
      float val = aq * FX[ti] * FY[tj] * FZ[tk];
      
      int gti = ti + ai;  if ( gti >= K1 ) gti -= K1;
      int gtj = tj + aj;  if ( gtj >= K2 ) gtj -= K2;
      int gtk = ak;  if ( gtk >= K3 ) gtk -= K3;
      gtk += tk;  // padded to increase coalescing (maybe, but reduces re-use)
      float *dest = q_arr[gti * K2 + gtj];  // could load as constant
      atomicAdd(dest+gtk,val);
    }
  } // i_atom
}


template <int order>
__global__ void cuda_pme_charges_dev(
  const float * __restrict__ coeffs, \
  float * const * __restrict__ q_arr, int * __restrict__ f_arr, int * __restrict__ fz_arr, \
  float * __restrict__ a_data, int n_atoms, \
  int K1, int K2, int K3
) {
  BSPLINE_DEFS
  __shared__ volatile float bspline_factors[warps_per_block][3][order];
  int atoms_this_warp = atoms_per_warp;
  {
    int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp;
    if ( aidx + atoms_per_warp > n_atoms ) {
      atoms_this_warp = n_atoms - aidx;  // may be negative
      if ( atoms_this_warp < 0 ) atoms_this_warp = 0;
    }
    if ( THREAD < atoms_this_warp*7 ) {
      atoms[WARP].a1d[THREAD] = a_data[aidx*7+THREAD];
    }
  }
  for ( int i_atom=0; i_atom < atoms_this_warp; ++i_atom ) {
    float aq=AQ;
    int ai=(int)AI;
    int aj=(int)AJ;
    int ak=(int)AK;
#if 0
if ( THREAD == 0 && (
  AI != (int)AI || AJ != (int)AJ || AK != (int)AK ||
  AQ < -30.f || AQ > 30.f || AX < 0.f || AX > 1.f ||
  AY < 0.f || AY > 1.f || AZ < 0.f || AZ > 1.f )
) {
  int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp + i_atom;
  printf("%d/%d %f %f %f %f %f %f %f\n", aidx, n_atoms, AI, AJ, AK, AQ, AX, AY, AZ );
  continue;
}
#endif
    if ( THREAD < 3 * order ) {
      const float af = ATOM[THREAD/order];  // x, y, or z
      float f = bspline_coeffs.a2d[0][THREAD % order];
      for ( int i=1; i<order; ++i ) {
        f = af * f + bspline_coeffs.a2d[i][THREAD % order];
      }
      bspline_factors[WARP][THREAD/order][THREAD % order] = f;
    }
    for ( int i=THREAD; i < order*order; i+=32 ) {
      int ti = i/order;
      int tj = i%order;
      int gti = ti + ai;  if ( gti >= K1 ) gti -= K1;
      int gtj = tj + aj;  if ( gtj >= K2 ) gtj -= K2;
#if 0
if ( gti < 0 || gtj < 0 || gti >= K1 || gtj >= K2 ) {
  int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp + i_atom;
  printf("%d/%d %d %d %d %f %f %f %f %f %f %f\n", aidx, n_atoms, i, gti, gtj, AI, AJ, AK, AQ, AX, AY, AZ);
} else
#endif
      f_arr[gti * K2 + gtj] = 1;
    }
    if ( THREAD < order ) {
      int gtk = ak;  if ( gtk >= K3 ) gtk -= K3;
      gtk += THREAD;  // padded to increase coalescing (maybe, but reduces re-use)
#if 0
if ( gtk < 0 || gtk >= K3+order-1 ) {
  int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp + i_atom;
  printf("%d/%d %d %d %f %f %f %f %f %f %f\n", aidx, n_atoms, THREAD, gtk, AI, AJ, AK, AQ, AX, AY, AZ);
} else
#endif
      fz_arr[gtk] = 1;
    }
    for ( int i=THREAD; i < order*order*order; i+=32 ) {
      int ti = i/(order*order);
      int tj = (i/order)%order;
      int tk = i%order;
      float val = aq * FX[ti] * FY[tj] * FZ[tk];
      int gti = ti + ai;  if ( gti >= K1 ) gti -= K1;
      int gtj = tj + aj;  if ( gtj >= K2 ) gtj -= K2;
      int gtk = ak;  if ( gtk >= K3 ) gtk -= K3;
      gtk += tk;  // padded to increase coalescing (maybe, but reduces re-use)
      float *dest = q_arr[gti * K2 + gtj];  // could load as constant
#if 0
if ( gti < 0 || gtj < 0 || gtk < 0 || gti >= K1 || gtj >= K2 || gtk >= K3+order-1 || dest == 0 ) {
  int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp + i_atom;
  printf("%d/%d %d %d %d %d %f %f %f %ld %f %f %f %f\n", aidx, n_atoms, i, gti, gtj, gtk, AI, AJ, AK, dest, AQ, AX, AY, AZ);
} else
#endif
       atomicAdd(dest+gtk,val);
    }
  } // i_atom
}

template <int order>
__global__ void cuda_pme_forces_dev(
  const float * __restrict__ coeffs, \
  float * const * __restrict__ q_arr, \
  float * const * __restrict__ afn_g, \
  /* float * __restrict__ a_data,
  float * __restrict__ f_data, int n_atoms, */
  int K1, int K2, int K3
) {
  __shared__ float *afn_s[3];
  if ( threadIdx.x < 3 ) afn_s[threadIdx.x] = afn_g[3*blockIdx.y+threadIdx.x];
  BSPLINE_DEFS
  __shared__ volatile union {
    float a2d[atoms_per_warp][3];
    float a1d[atoms_per_warp*3];
  } force_output[warps_per_block];
  __shared__ float2 bspline_2factors[warps_per_block][3][order];
  __shared__ volatile union {
    float a2d[32][3];
    float a1d[32*3];
  } force_reduction[warps_per_block];
  int atoms_this_warp = atoms_per_warp;
  {
    const int n_atoms = afn_s[2] - afn_s[1];
    int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp;
    if ( aidx + atoms_per_warp > n_atoms ) {
      atoms_this_warp = n_atoms - aidx;  // may be negative
      if ( atoms_this_warp < 0 ) atoms_this_warp = 0;
    }
    if ( THREAD < atoms_this_warp*7 ) {
      atoms[WARP].a1d[THREAD] = afn_s[0][aidx*7+THREAD];
    }
  }
  for ( int i_atom=0; i_atom < atoms_this_warp; ++i_atom ) {
    float aq=AQ;
    int ai=(int)AI;
    int aj=(int)AJ;
    int ak=(int)AK;
    if ( THREAD < 3 * order ) {
      const float af = ATOM[THREAD/order];  // x, y, or z
      float df = 0.f;
      float c = bspline_coeffs.a2d[0][THREAD % order];
      float f = c;
      for ( int i=1; i<order; ++i ) {
        df = af * df - (order-i) * c;
        c = bspline_coeffs.a2d[i][THREAD % order];
        f = af * f + c;
      }
      bspline_2factors[WARP][THREAD/order][THREAD % order] = make_float2(f,df);
    }
    __threadfence_block();
    float force_x = 0.f;
    float force_y = 0.f;
    float force_z = 0.f;

    for ( int i=THREAD; i < order*order*order; i+=32 ) {
      int ti = i/(order*order);
      int tj = (i/order)%order;
      int tk = i%order;

      const float2 fx2 = FX2[ti];
      const float2 fy2 = FY2[tj];
      const float2 fz2 = FZ2[tk];

      const float fx = fx2.x;
      const float fy = fy2.x;
      const float fz = fz2.x;
      const float dfx = fx2.y;
      const float dfy = fy2.y;
      const float dfz = fz2.y;

      float dx = K1 * aq * dfx * fy * fz;
      float dy = K2 * aq * fx * dfy * fz;
      float dz = K3 * aq * fx * fy * dfz;
    
      int gti = ti + ai;  if ( gti >= K1 ) gti -= K1;
      int gtj = tj + aj;  if ( gtj >= K2 ) gtj -= K2;
      int gtk = ak;  if ( gtk >= K3 ) gtk -= K3;
      
      gtk += tk;  // padded to increase coalescing (maybe, but reduces re-use)
      const float * __restrict__ src = q_arr[gti * K2 + gtj];  // could load as constant
      float pot = src[gtk];
      force_x += dx * pot;
      force_y += dy * pot;
      force_z += dz * pot;
    }
    force_reduction[WARP].a2d[THREAD][0] = force_x;
    force_reduction[WARP].a2d[THREAD][1] = force_y;
    force_reduction[WARP].a2d[THREAD][2] = force_z;
    if ( THREAD < 24 ) {
      force_reduction[WARP].a1d[THREAD] += force_reduction[WARP].a1d[THREAD + 24]
                                         + force_reduction[WARP].a1d[THREAD + 48]
                                         + force_reduction[WARP].a1d[THREAD + 72];
    }
    if ( THREAD < 6 ) {
      force_reduction[WARP].a1d[THREAD] += force_reduction[WARP].a1d[THREAD + 6]
                                         + force_reduction[WARP].a1d[THREAD + 12]
                                         + force_reduction[WARP].a1d[THREAD + 18];
    }
    if ( THREAD < 3 ) {
      force_output[WARP].a2d[i_atom][THREAD] = force_reduction[WARP].a1d[THREAD]
                                             + force_reduction[WARP].a1d[THREAD + 3];
    }
  } // i_atom
  if ( THREAD < atoms_this_warp*3 ) {
    int aidx = blockIdx.x * atoms_per_block + WARP * atoms_per_warp;
    afn_s[1][aidx*3+THREAD] = force_output[WARP].a1d[THREAD];
  }
}

CUDA_PME_CHARGES_PROTOTYPE {
  int nblocks = (n_atoms + atoms_per_block - 1) / atoms_per_block;
  if ( ! nblocks ) return;

#define CALL(ORDER) if ( order == ORDER ) \
                      cuda_pme_charges_dev<ORDER><<<nblocks,32*warps_per_block,0,stream>>>( \
                        coeffs, q_arr, f_arr, fz_arr, a_data, n_atoms, K1, K2, K3)
  CALL(4); 
  else CALL(6);
  else CALL(8);
  else NAMD_die("unsupported PMEInterpOrder");
#undef CALL
}

CUDA_PME_CHARGES_BATCHED_PROTOTYPE {
  int nblocksX = (n_max_atoms + atoms_per_block - 1) / atoms_per_block;
  if ( (! nblocksX) || (! numPatches) ) return;
  dim3 gridSize;
  gridSize.x = nblocksX;
  gridSize.y = numPatches;
  gridSize.z = 1;
#define CALL(ORDER) if ( order == ORDER ) \
                      cuda_pme_charges_batched_dev<ORDER><<<gridSize,32*warps_per_block,0,stream>>>( \
                        coeffs, q_arr, f_arr, fz_arr, a_data_ptr, n_atoms_ptr, K1_ptr, K2_ptr, K3_ptr)
  CALL(4);
  else CALL(6);
  else CALL(8);
  else NAMD_die("unsupported PMEInterpOrder");
#undef CALL
}


#ifdef CUDA_VERSION
#if CUDA_VERSION < 4020
void dummy() { }
#define cudaFuncSetSharedMemConfig(X,Y) dummy()
#endif
#endif


CUDA_PME_FORCES_PROTOTYPE {
  int nblocks = (maxn + atoms_per_block - 1) / atoms_per_block;
  if ( (! nblocks) || (! dimy) ) return;

#define CALL(ORDER) if ( order == ORDER ) \
                      cudaFuncSetSharedMemConfig(cuda_pme_forces_dev<ORDER>,cudaSharedMemBankSizeEightByte), \
                      cuda_pme_forces_dev<ORDER><<<dim3(nblocks,dimy),32*warps_per_block,0,stream>>>( \
                        coeffs, q_arr, afn, /* a_data, f_data, n_atoms, */ K1, K2, K3)
  CALL(4); 
  else CALL(6);
  else CALL(8);
  else NAMD_die("unsupported PMEInterpOrder");
#undef CALL
}

#endif // NAMD_CUDA

