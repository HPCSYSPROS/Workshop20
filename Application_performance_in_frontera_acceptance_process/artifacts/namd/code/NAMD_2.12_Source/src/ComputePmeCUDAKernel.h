#ifdef NAMD_CUDA

#ifndef __CUDACC__
#undef __align__(X)
#define __align__(X)
#endif

void cuda_init_bspline_coeffs(float **c, float **dc, int order);

#define CUDA_PME_CHARGES_PROTOTYPE \
void cuda_pme_charges( \
  const float *coeffs, \
  float * const *q_arr, int *f_arr, int *fz_arr, \
  float *a_data, int n_atoms, \
  int K1, int K2, int K3, \
  int order, cudaStream_t stream)

CUDA_PME_CHARGES_PROTOTYPE;


#define CUDA_PME_CHARGES_BATCHED_PROTOTYPE \
void cuda_pme_charges_batched( \
  const float *coeffs, \
  float * const *q_arr, int *f_arr, int *fz_arr, \
  float **a_data_ptr, int* n_atoms_ptr, \
  int* K1_ptr, int* K2_ptr, int* K3_ptr, \
  int order, int numPatches, int n_max_atoms, cudaStream_t stream)

CUDA_PME_CHARGES_BATCHED_PROTOTYPE;


#define CUDA_PME_FORCES_PROTOTYPE \
void cuda_pme_forces( \
  const float *coeffs, \
  float * const *q_arr, \
  float * const *afn, int dimy, int maxn, \
  /* float *a_data, float *f_data, int n_atoms, */ \
  int K1, int K2, int K3, \
  int order, cudaStream_t stream)

CUDA_PME_FORCES_PROTOTYPE;

#endif // NAMD_CUDA

