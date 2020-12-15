
#include "CudaUtils.h"
#include "ComputeNonbondedCUDAKernel.h"
#include <stdio.h>

#ifdef NAMD_CUDA

#ifdef WIN32
#define __thread __declspec(thread)
#endif

texture<unsigned int, 1, cudaReadModeElementType> tex_exclusions;
static __thread int exclusions_size;
static __thread unsigned int *exclusions;

__constant__ unsigned int const_exclusions[MAX_CONST_EXCLUSIONS];
static __thread unsigned int *overflow_exclusions;

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

void cuda_bind_exclusions(const unsigned int *t, int n) {
  exclusions_size = n;
  static __thread int exclusions_alloc;
  if ( exclusions && exclusions_alloc < exclusions_size ) {
    cudaFree(exclusions);
    cuda_errcheck("freeing exclusions");
    cudaFree(overflow_exclusions);
    cuda_errcheck("freeing overflow_exclusions");
    exclusions = 0;
  }
  if ( ! exclusions ) {
    exclusions_alloc = exclusions_size;
    cudaMalloc((void**) &exclusions, n*sizeof(unsigned int));
    cuda_errcheck("malloc exclusions");
    cudaMalloc((void**) &overflow_exclusions, n*sizeof(unsigned int));
    cuda_errcheck("malloc overflow_exclusions");
  }
  cudaMemcpy(exclusions, t, n*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy exclusions");
  tex_exclusions.normalized = false;
  tex_exclusions.addressMode[0] = cudaAddressModeClamp;
  tex_exclusions.filterMode = cudaFilterModePoint;
  cudaBindTexture(NULL, tex_exclusions, exclusions, n*sizeof(unsigned int));
  cuda_errcheck("binding exclusions to texture");

  cudaMemcpy(overflow_exclusions, t,
		n*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to overflow_exclusions");
  int nconst = ( n < MAX_CONST_EXCLUSIONS ? n : MAX_CONST_EXCLUSIONS );
  cudaMemcpyToSymbol(const_exclusions, t, nconst*sizeof(unsigned int), 0);
  cuda_errcheck("memcpy to const_exclusions");
}


texture<float2, 1, cudaReadModeElementType> lj_table;
static __thread int lj_table_size;

void cuda_bind_lj_table(const float2 *t, int _lj_table_size) {
    static __thread float2 *ct;
    static __thread int lj_table_alloc;
    lj_table_size = _lj_table_size;
    if ( ct && lj_table_alloc < lj_table_size ) {
      cudaFree(ct);
      cuda_errcheck("freeing lj table");
      ct = 0;
    }
    if ( ! ct ) {
      lj_table_alloc = lj_table_size;
      cudaMalloc((void**) &ct, lj_table_size*lj_table_size*sizeof(float2));
      cuda_errcheck("allocating lj table");
    }
    cudaMemcpy(ct, t, lj_table_size*lj_table_size*sizeof(float2),
                                            cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to lj table");

    lj_table.normalized = false;
    lj_table.addressMode[0] = cudaAddressModeClamp;
    lj_table.filterMode = cudaFilterModePoint;

    cudaBindTexture((size_t*)0, lj_table, ct,
        lj_table_size*lj_table_size*sizeof(float2));
    cuda_errcheck("binding lj table to texture");
}


texture<float4, 1, cudaReadModeElementType> force_table;
texture<float4, 1, cudaReadModeElementType> energy_table;

void cuda_bind_force_table(const float4 *t, const float4 *et) {
    static __thread cudaArray *ct;
    static __thread cudaArray *ect;
    if ( ! ct ) {
      cudaMallocArray(&ct, &force_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating force table");
    }
    if ( ! ect ) {
      cudaMallocArray(&ect, &energy_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating energy table");
    }
    cudaMemcpyToArray(ct, 0, 0, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    // cudaMemcpy(ct, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to force table");
    cudaMemcpyToArray(ect, 0, 0, et, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to energy table");

    force_table.normalized = true;
    force_table.addressMode[0] = cudaAddressModeClamp;
    force_table.addressMode[1] = cudaAddressModeClamp;
    force_table.filterMode = cudaFilterModeLinear;

    energy_table.normalized = true;
    energy_table.addressMode[0] = cudaAddressModeClamp;
    energy_table.addressMode[1] = cudaAddressModeClamp;
    energy_table.filterMode = cudaFilterModeLinear;

    cudaBindTextureToArray(force_table, ct);
    cuda_errcheck("binding force table to texture");

    cudaBindTextureToArray(energy_table, ect);
    cuda_errcheck("binding energy table to texture");
}

static __thread int num_patches;
static __thread int num_virials;
static __thread int num_atoms;
// Size of the device array followed by the array pointer
static __thread int patch_pairs_size;
static __thread patch_pair* patch_pairs;
static __thread int atom_params_size;
static __thread atom_param* atom_params;
static __thread int vdw_types_size;
static __thread int* vdw_types;
static __thread int atoms_size;
static __thread atom* atoms;

static __thread int tmpforces_size;
static __thread float4* tmpforces;
static __thread int slow_tmpforces_size;
static __thread float4* slow_tmpforces;

static __thread int tmpvirials_size;
static __thread float* tmpvirials;
static __thread int slow_tmpvirials_size;
static __thread float* slow_tmpvirials;

static __thread int global_counters_size;
static __thread unsigned int* global_counters;
static __thread int plist_size;
static __thread unsigned int* plist;
static __thread int exclmasks_size;
static __thread exclmask* exclmasks;
// Device pointers to the page-locked host arrays (provided by ComputeNonbondedCUDA -class)
static __thread float4* forces;
static __thread float4* slow_forces;
static __thread int* force_ready_queue;
static __thread float* virials;
static __thread float* slow_virials;
static __thread int* block_order;

//GBIS arrays
static __thread int intRad0D_size;
static __thread float *intRad0D;

static __thread int intRadSD_size;
static __thread float *intRadSD;

static __thread GBReal *psiSumD;        // host-mapped memory

static __thread int tmp_psiSumD_size;
static __thread GBReal *tmp_psiSumD;

static __thread int bornRadD_size;
static __thread float *bornRadD;

static __thread GBReal *dEdaSumD;       // host-mapped memory

static __thread int tmp_dEdaSumD_size;
static __thread GBReal *tmp_dEdaSumD;

static __thread int dHdrPrefixD_size;
static __thread float *dHdrPrefixD;

static __thread int GBIS_P1_counters_size;
static __thread unsigned int *GBIS_P1_counters;

static __thread int GBIS_P2_counters_size;
static __thread unsigned int *GBIS_P2_counters;

static __thread int GBIS_P3_counters_size;
static __thread unsigned int *GBIS_P3_counters;

static __thread float *energy_gbis;     // host-mapped memory

static __thread int tmp_energy_gbis_size;
static __thread float *tmp_energy_gbis;

__thread int max_grid_size;

__thread cudaStream_t stream;
__thread cudaStream_t stream2;
 
void cuda_init() {
  patch_pairs_size = 0;
  patch_pairs = NULL;

  atom_params_size = 0;
  atom_params = NULL;

  vdw_types_size = 0;
  vdw_types = NULL;

  atoms_size = 0;
  atoms = NULL;  

  tmpforces_size = 0;
  tmpforces = NULL;

  slow_tmpforces_size = 0;
  slow_tmpforces = NULL;

  tmpvirials_size = 0;
  tmpvirials = NULL;

  slow_tmpvirials_size = 0;
  slow_tmpvirials = NULL;

  global_counters_size = 0;
  global_counters = NULL;

  plist_size = 0;
  plist = NULL;

  exclmasks_size = 0;
  exclmasks = NULL;

  forces = NULL;
  slow_forces = NULL;

  force_ready_queue = NULL;

  exclusions_size = 0;
  exclusions = NULL;

  // --------------------
  // For GBIS
  // --------------------
  intRad0D_size = 0;
  intRad0D = NULL;

  intRadSD_size = 0;
  intRadSD = NULL;

  psiSumD = NULL;        // host-mapped memory

  tmp_psiSumD_size = 0;
  tmp_psiSumD = NULL;

  bornRadD_size = 0;
  bornRadD = NULL;

  dEdaSumD = NULL;       // host-mapped memory

  tmp_dEdaSumD_size = 0;
  tmp_dEdaSumD = NULL;

  dHdrPrefixD_size = 0;
  dHdrPrefixD = NULL;

  GBIS_P1_counters_size = 0;
  GBIS_P1_counters = NULL;

  GBIS_P2_counters_size = 0;
  GBIS_P2_counters = NULL;

  GBIS_P3_counters_size = 0;
  GBIS_P3_counters = NULL;

  energy_gbis = NULL;     // host-mapped memory

  tmp_energy_gbis_size = 0;
  tmp_energy_gbis = NULL;

  int dev;
  cudaGetDevice(&dev);
  cuda_errcheck("cudaGetDevice");
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  cuda_errcheck("cudaGetDeviceProperties");
  max_grid_size = deviceProp.maxGridSize[1];
}

void cuda_bind_patch_pairs(patch_pair *h_patch_pairs, int npatch_pairs,
			   int npatches, int natoms, int plist_len, 
			   int nexclmask) {
  num_patches = npatches;
  num_virials = npatches;
  num_atoms = natoms;
  reallocate_device<patch_pair>(&patch_pairs, &patch_pairs_size, npatch_pairs, 1.2f);
  reallocate_device<atom>(&atoms, &atoms_size, num_atoms, 1.2f);
  reallocate_device<atom_param>(&atom_params, &atom_params_size, num_atoms, 1.2f);
  reallocate_device<int>(&vdw_types, &vdw_types_size, num_atoms, 1.2f);
  reallocate_device<unsigned int>(&global_counters, &global_counters_size, num_patches+2, 1.2f);
  reallocate_device<float4>(&tmpforces, &tmpforces_size, num_atoms, 1.2f);
  reallocate_device<float4>(&slow_tmpforces, &slow_tmpforces_size, num_atoms, 1.2f);
  reallocate_device<unsigned int>(&plist, &plist_size, plist_len, 1.2f);
  reallocate_device<exclmask>(&exclmasks, &exclmasks_size, nexclmask, 1.2f);
  reallocate_device<float>(&tmpvirials, &tmpvirials_size, num_patches*16, 1.2f);
  reallocate_device<float>(&slow_tmpvirials, &slow_tmpvirials_size, num_patches*16, 1.2f);

  // For GBIS
  reallocate_device<unsigned int>(&GBIS_P1_counters, &GBIS_P1_counters_size, num_patches, 1.2f);
  reallocate_device<unsigned int>(&GBIS_P2_counters, &GBIS_P2_counters_size, num_patches, 1.2f);
  reallocate_device<unsigned int>(&GBIS_P3_counters, &GBIS_P3_counters_size, num_patches, 1.2f);
  reallocate_device<float>(&intRad0D, &intRad0D_size, num_atoms, 1.2f);
  reallocate_device<float>(&intRadSD, &intRadSD_size, num_atoms, 1.2f);
  reallocate_device<GBReal>(&tmp_psiSumD, &tmp_psiSumD_size, num_atoms, 1.2f);
  reallocate_device<float>(&bornRadD, &bornRadD_size, num_atoms, 1.2f);
  reallocate_device<GBReal>(&tmp_dEdaSumD, &tmp_dEdaSumD_size, num_atoms, 1.2f);
  reallocate_device<float>(&dHdrPrefixD, &dHdrPrefixD_size, num_atoms, 1.2f);
  reallocate_device<float>(&tmp_energy_gbis, &tmp_energy_gbis_size, num_patches, 1.2f);

  cudaMemcpy(patch_pairs, h_patch_pairs, npatch_pairs*sizeof(patch_pair), cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to patch_pairs");

  cudaMemset(global_counters, 0, (num_patches+2)*sizeof(unsigned int));
  cuda_errcheck("memset global_counters");

  cudaMemset(GBIS_P1_counters, 0, num_patches*sizeof(unsigned int));
  cuda_errcheck("memset GBIS_P1_counters");

  cudaMemset(GBIS_P2_counters, 0, num_patches*sizeof(unsigned int));
  cuda_errcheck("memset GBIS_P2_counters");

  cudaMemset(GBIS_P3_counters, 0, num_patches*sizeof(unsigned int));
  cuda_errcheck("memset GBIS_P3_counters");

  cudaMemset(tmpforces, 0, num_atoms*sizeof(float4));
  cuda_errcheck("memset tmpforces");

  cudaMemset(tmpvirials, 0, num_patches*sizeof(float)*16);
  cuda_errcheck("memset tmpvirials");

  cudaMemset(slow_tmpforces, 0, num_atoms*sizeof(float4));
  cuda_errcheck("memset slow_tmpforces");

  cudaMemset(slow_tmpvirials, 0, num_patches*sizeof(float)*16);
  cuda_errcheck("memset slow_tmpvirials");

}

void cuda_bind_atom_params(const atom_param *t) {
  cudaMemcpyAsync(atom_params, t, num_atoms * sizeof(atom_param),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to atom_params");
}

void cuda_bind_vdw_types(const int *t) {
  cudaMemcpyAsync(vdw_types, t, num_atoms * sizeof(int),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to vdw_types");
}

void cuda_bind_atoms(const atom *a) {
  cuda_errcheck("before memcpy to atoms");
  cudaMemcpyAsync(atoms, a, num_atoms * sizeof(atom),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to atoms");
}

void cuda_bind_forces(float4 *f, float4 *f_slow) {
  cudaHostGetDevicePointer(&forces, f, 0);
  cuda_errcheck("cudaHostGetDevicePointer forces");
  cudaHostGetDevicePointer(&slow_forces, f_slow, 0);
  cuda_errcheck("cudaHostGetDevicePointer slow_forces");
}

void cuda_bind_virials(float *v, int *queue, int *blockorder) {
  cudaHostGetDevicePointer(&virials, v, 0);
  cuda_errcheck("cudaHostGetDevicePointer virials");
  slow_virials = virials + num_virials*16;
  cudaHostGetDevicePointer(&force_ready_queue, queue, 0);
  cuda_errcheck("cudaHostGetDevicePointer force_ready_queue");
  cudaHostGetDevicePointer(&block_order, blockorder, 0);
  cuda_errcheck("cudaHostGetDevicePointer block_order");
}

//GBIS bindings
void cuda_bind_GBIS_energy(float *e) {
  cudaHostGetDevicePointer(&energy_gbis, e, 0);
  cuda_errcheck("cudaHostGetDevicePointer energy_gbis");
}
void cuda_bind_GBIS_intRad(float *intRad0H, float *intRadSH) {
  cudaMemcpyAsync(intRad0D, intRad0H, num_atoms * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(intRadSD, intRadSH, num_atoms * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to intRad");
}

void cuda_bind_GBIS_psiSum(GBReal *psiSumH) {
  cudaHostGetDevicePointer(&psiSumD, psiSumH, 0);
  cuda_errcheck("cudaHostGetDevicePointer psiSum");
}

void cuda_bind_GBIS_bornRad(float *bornRadH) {
  cudaMemcpyAsync(bornRadD, bornRadH, num_atoms * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to bornRad");
}

void cuda_bind_GBIS_dEdaSum(GBReal *dEdaSumH) {
  cudaHostGetDevicePointer(&dEdaSumD, dEdaSumH, 0);
  cuda_errcheck("cudaHostGetDevicePointer dEdaSum");
}

void cuda_bind_GBIS_dHdrPrefix(float *dHdrPrefixH) {
  cudaMemcpyAsync(dHdrPrefixD, dHdrPrefixH, num_atoms * sizeof(float),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to dHdrPrefix");
}
// end GBIS methods

#if 0
void cuda_load_forces(float4 *f, float4 *f_slow, int begin, int count) {
  // printf("load forces %d %d %d\n",begin,count,num_atoms);
  cudaMemcpyAsync(f+begin, forces+begin, count * sizeof(float4),
				cudaMemcpyDeviceToHost, stream);
  if ( f_slow ) {
    cudaMemcpyAsync(f_slow+begin, slow_forces+begin, count * sizeof(float4),
				cudaMemcpyDeviceToHost, stream);
  }
  cuda_errcheck("memcpy from forces");
}

void cuda_load_virials(float *v, int doSlow) {
  int count = force_lists_size;
  if ( doSlow ) count *= 2;
  cudaMemcpyAsync(v, virials, count * 16*sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cuda_errcheck("memcpy from virials");
}
#endif

#if 0
__host__ __device__ static int3 patch_coords_from_id(
        dim3 PATCH_GRID, int id) {

  return make_int3( id % PATCH_GRID.x,
                ( id / PATCH_GRID.x ) % PATCH_GRID.y,
                id / ( PATCH_GRID.x * PATCH_GRID.y ) );
}

__host__ __device__ static int patch_id_from_coords(
        dim3 PATCH_GRID, int3 coords) {

  // handles periodic boundaries
  int x = (coords.x + 4 * PATCH_GRID.x) % PATCH_GRID.x;
  int y = (coords.y + 4 * PATCH_GRID.y) % PATCH_GRID.y;
  int z = (coords.z + 4 * PATCH_GRID.z) % PATCH_GRID.z;

  return ( z * PATCH_GRID.y + y ) * PATCH_GRID.x + x;
}

__host__ __device__ static int3 patch_offset_from_neighbor(int neighbor) {

  // int3 coords = patch_coords_from_id(make_uint3(3,3,3), 13 + neighbor);
  int3 coords = patch_coords_from_id(make_uint3(3,3,3), neighbor);
  return make_int3(coords.x - 1, coords.y - 1, coords.z - 1);

}
#endif
 
#define BLOCK_SIZE 128
#define SHARED_SIZE 32

#define MAKE_PAIRLIST
#define DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef MAKE_PAIRLIST
#define DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"

void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc,
			   float cutoff2, float plcutoff2,
			   int cbegin, int ccount, int ctotal,
			   int doSlow, int doEnergy, int usePairlists, int savePairlists,
			   int doStreaming, int saveOrder, cudaStream_t &strm) {

 if ( ccount ) {
   if ( usePairlists ) {
     if ( ! savePairlists ) plcutoff2 = 0.;
   } else {
     plcutoff2 = cutoff2;
   }

   cudaMemsetAsync(tmpforces, 0, num_atoms*sizeof(float4), strm);
   cudaMemsetAsync(tmpvirials, 0, num_patches*sizeof(float)*16, strm);
   if ( doSlow ) {
     cudaMemsetAsync(slow_tmpforces, 0, num_atoms*sizeof(float4), strm);
     cudaMemsetAsync(slow_tmpvirials, 0, num_patches*sizeof(float)*16, strm);
   }

   int grid_dim = max_grid_size;  // maximum allowed
   for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
     if ( grid_dim > ccount - cstart ) grid_dim = ccount - cstart;

     dim3 nthread3(WARPSIZE, NUM_WARP, 1);

#define CALL(X) X<<< grid_dim, nthread3, 0, strm >>>			\
       (patch_pairs, atoms, atom_params, vdw_types, plist,		\
	tmpforces, (doSlow?slow_tmpforces:NULL),			\
	forces, (doSlow?slow_forces:NULL),				\
	tmpvirials, (doSlow?slow_tmpvirials:NULL),			\
	virials, (doSlow?slow_virials:NULL),				\
	global_counters, (doStreaming?force_ready_queue:NULL),		\
	overflow_exclusions, num_patches,				\
	cbegin+cstart, ctotal, (saveOrder?block_order:NULL),		\
	exclmasks, lj_table_size,					\
	lata, latb, latc, cutoff2, plcutoff2, doSlow)
//end definition

     if ( doEnergy ) {
       if ( doSlow ) {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_slow_energy_pairlist);
         else CALL(dev_nonbonded_slow_energy);
       } else {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_energy_pairlist);
         else CALL(dev_nonbonded_energy);
       }
     } else {
       if ( doSlow ) {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_slow_pairlist);
         else CALL(dev_nonbonded_slow);
       } else {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_pairlist);
         else CALL(dev_nonbonded);
       }
     }

     cuda_errcheck("dev_nonbonded");
   }
 }

}

//import GBIS Kernel definitions
#include "ComputeGBISCUDAKernel.h"

//////////////////////////////////////////
//  GBIS P1
//////////////////////////////////////////
void cuda_GBIS_P1(
	int cbegin,
  int ccount,
  int pbegin,
  int pcount,
  float a_cut,
  float rho_0,
  float3 lata,
  float3 latb,
  float3 latc,
  cudaStream_t &strm
) {

 if ( ccount ) {
  cudaMemsetAsync(tmp_psiSumD, 0, num_atoms*sizeof(GBReal), strm);

  int grid_dim = max_grid_size;  // maximum allowed
  for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
    if (grid_dim > ccount - cstart) {
      grid_dim = ccount - cstart;
    }

    dim3 nthread3(WARPSIZE, NUM_WARP, 1);
    GBIS_P1_Kernel<<<grid_dim, nthread3, 0, strm>>>(
      patch_pairs+cbegin+cstart,
      atoms,
      intRad0D,
      intRadSD,
      tmp_psiSumD,
      psiSumD,
      a_cut,
      rho_0,
      lata,
      latb,
      latc,
      GBIS_P1_counters 
      );
    cuda_errcheck("dev_GBIS_P1");
  } // end for
 }
} // end GBIS P1

//////////////////////////////////////////
//  GBIS P2
//////////////////////////////////////////
void cuda_GBIS_P2(
	int cbegin,
  int ccount,
  int pbegin,
  int pcount,
  float a_cut,
  float r_cut,
  float scaling,
  float kappa,
  float smoothDist,
  float epsilon_p,
  float epsilon_s,
  float3 lata,
  float3 latb,
  float3 latc,
  int doEnergy,
  int doFullElec,
  cudaStream_t &strm
) {

 if ( ccount ) {
  cudaMemsetAsync(tmp_dEdaSumD, 0, num_atoms*sizeof(GBReal), strm);
  cudaMemsetAsync(tmp_energy_gbis, 0, num_patches*sizeof(float), strm);

  int grid_dim = max_grid_size;  // maximum allowed
  for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
    if (grid_dim > ccount - cstart)
      grid_dim = ccount - cstart;

    dim3 nthread3(WARPSIZE, NUM_WARP, 1);
    GBIS_P2_Kernel<<<grid_dim, nthread3, 0, strm>>>(
      patch_pairs+cbegin+cstart,
      atoms,
      bornRadD,
      tmp_dEdaSumD,
      dEdaSumD,
      a_cut,
      r_cut,
      scaling,
      kappa,
      smoothDist,
      epsilon_p,
      epsilon_s,
      lata,
      latb,
      latc,
      doEnergy,
      doFullElec,
      tmpforces,
      forces,
      tmp_energy_gbis,
      energy_gbis,
      GBIS_P2_counters 
      );
    cuda_errcheck("dev_GBIS_P2");
  } // end for
 }
} // end P2

//////////////////////////////////////////
//  GBIS P3
//////////////////////////////////////////
void cuda_GBIS_P3(
	int cbegin,
  int ccount,
  int pbegin,
  int pcount,
  float a_cut,
  float rho_0,
  float scaling,
  float3 lata,
  float3 latb,
  float3 latc,
  cudaStream_t &strm
) {
  int grid_dim = max_grid_size;  // maximum allowed
  for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
    if (grid_dim > ccount - cstart)
      grid_dim = ccount - cstart;

    dim3 nthread3(WARPSIZE, NUM_WARP, 1);
    GBIS_P3_Kernel<<<grid_dim, nthread3, 0, strm>>>(
      patch_pairs+cbegin+cstart,
      atoms,
      intRad0D,
      intRadSD,
      dHdrPrefixD,
      a_cut,
      rho_0,
      scaling,
      lata,
      latb,
      latc,
      slow_tmpforces,
      slow_forces,
      GBIS_P3_counters 
      );
    cuda_errcheck("dev_GBIS_P3");
  }
}

#if 0
int cuda_stream_finished() {
  return ( cudaStreamQuery(stream) == cudaSuccess );
}
#endif


#else  // NAMD_CUDA

// for make depends
#include "ComputeNonbondedCUDAKernelBase.h"
#include "ComputeGBISCUDAKernel.h"

#endif  // NAMD_CUDA

