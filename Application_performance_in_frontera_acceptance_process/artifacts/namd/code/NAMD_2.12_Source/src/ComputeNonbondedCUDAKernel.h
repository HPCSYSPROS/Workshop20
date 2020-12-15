#ifdef NAMD_CUDA
//this type defined in multiple files
typedef float GBReal;

void cuda_errcheck(const char *msg);

#ifndef __CUDACC__
#undef __align__(X)
#define __align__(X)
#endif


// Number of warps per block for Non-bonded CUDA kernel
#define NUM_WARP 4
#define WARPSIZE 32

// Exclusion mask: bit 1 = atom pair is excluded
struct exclmask {
  unsigned int excl[32];
};

struct __align__(16) patch_pair {
  float3 offset;
  int patch1_start;      // Coordinate/force start for this patch
  int patch1_size;       // Size of the patch
  int patch2_start;
  int patch2_size;
  int patch1_ind;        // Patch index
  int patch2_ind;
  int patch1_num_pairs;  // Number of pairs that involve this patch
  int patch2_num_pairs;
  union {
    bool patch_done[2];      // After-GPU-computation shared memory temporary storage
    struct {
      int plist_start;       // Pair list start
      int plist_size;        // Pair list size
    };
  };
  int exclmask_start;    // Exclusion mask start
  int patch1_free_size;  // Size of the free atoms in patch
  int patch2_free_size;  // Size of the free atoms in patch
//  int pad1, pad2;
};

#define PATCH_PAIR_SIZE (sizeof(patch_pair)/4)

struct __align__(16) atom {  // must be multiple of 16!
  float3 position;
  float charge;
};

struct __align__(16) atom_param {  // must be multiple of 16!
  int vdw_type;
  int index;
  int excl_index;
  int excl_maxdiff;  // maxdiff == 0 -> only excluded from self
};

#define COPY_ATOM( DEST, SOURCE ) { \
  DEST.position.x = SOURCE.position.x; \
  DEST.position.y = SOURCE.position.y; \
  DEST.position.z = SOURCE.position.z; \
  DEST.charge = SOURCE.charge; \
  }

#define COPY_PARAM( DEST, SOURCE ) { \
  DEST.sqrt_epsilon = SOURCE.sqrt_epsilon; \
  DEST.half_sigma = SOURCE.half_sigma; \
  DEST.index = SOURCE.index; \
  DEST.excl_index = SOURCE.excl_index; \
  DEST.excl_maxdiff = SOURCE.excl_maxdiff; \
  }

#define COPY_ATOM_TO_SHARED( ATOM, PARAM, SHARED ) { \
    COPY_ATOM( SHARED, ATOM ) \
    COPY_PARAM( SHARED, PARAM ) \
  }

#define COPY_ATOM_FROM_SHARED( ATOM, PARAM, SHARED ) { \
    COPY_ATOM( ATOM, SHARED ) \
    COPY_PARAM( PARAM, SHARED ) \
  }

// 2^11 ints * 2^5 bits = 2^16 bits = range of unsigned short excl_index
// 2^27 ints * 2^5 bits = 2^32 bits = range of unsigned int excl_index
#define MAX_EXCLUSIONS (1<<27)
#define MAX_CONST_EXCLUSIONS 2048  // cache size is 8k

void cuda_bind_exclusions(const unsigned int *t, int n);

void cuda_bind_lj_table(const float2 *t, int _lj_table_size);

// #define FORCE_TABLE_SIZE 512
// maximum size of CUDA array 1D texture reference is 2^13 = 8192
// #define FORCE_TABLE_SIZE 8192
// CUDA docs lie, older devices can only handle 4096
#define FORCE_TABLE_SIZE 4096

void cuda_bind_force_table(const float4 *t, const float4 *et);

void cuda_init();

void cuda_bind_patch_pairs(patch_pair *h_patch_pairs, int npatch_pairs,
			   int npatches, int natoms, int nexclmask, int plist_len);

void cuda_bind_atom_params(const atom_param *t);
void cuda_bind_vdw_types(const int *t);

void cuda_bind_atoms(const atom *a);

void cuda_bind_forces(float4 *f, float4 *f_slow);

void cuda_bind_virials(float *v, int *queue, int *blockorder);

void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc,
			   float cutoff2, float plcutoff2,
			   int cbegin, int ccount, int ctotal,
			   int doSlow, int doEnergy, int usePairlists, int savePairlists,
			   int doStreaming, int saveOrder, cudaStream_t &strm);

//GBIS methods
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
  );
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
  );
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
  );

void cuda_bind_GBIS_intRad(float *intRad0H, float *intRadSH);
void cuda_bind_GBIS_energy(float *energy_gbis);
void cuda_bind_GBIS_psiSum(GBReal *psiSumH);
void cuda_bind_GBIS_bornRad(float *bornRadH);
void cuda_bind_GBIS_dEdaSum(GBReal *dEdaSumH);
void cuda_bind_GBIS_dHdrPrefix(float *dHdrPrefixH);

//end GBIS methods

int cuda_stream_finished();

#endif  // NAMD_CUDA

