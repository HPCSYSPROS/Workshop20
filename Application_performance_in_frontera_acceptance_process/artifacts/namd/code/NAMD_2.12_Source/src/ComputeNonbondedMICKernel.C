// NOTE: See ComputeNonbondedMICKernel.h for a general description of NAMD on MIC.

// Only compile the contents of this file if this build supports MIC.
#ifdef NAMD_MIC


#include <stdlib.h>
#include <stdio.h>
#include <offload.h>
#include "ComputeNonbondedMICKernel.h"
#include <assert.h>
#include <math.h>

// DMK - DEBUG
#if MIC_TRACK_DEVICE_MEM_USAGE != 0
  #include <sys/types.h>
  #include <unistd.h>
#endif

// Setup __ASSUME_ALIGNED macro to take the appropriate action based on macro flags
#if (CHECK_ASSUME_ALIGNED != 0)
  #define __ASSUME_ALIGNED(v) assert(((unsigned long long int)(v)) % 64 == 0); __assume_aligned((v), MIC_ALIGN);
#else
  #define __ASSUME_ALIGNED(v) __assume_aligned((v), MIC_ALIGN);
#endif

// Setup __ASSUME macro to take the appropriate action based on macro flags
#if (CHECK_ASSUME != 0)
  #define __ASSUME(s) assert(s); __assume(s);
#else
  #define __ASSUME(s) __assume(s);
#endif

// Setup RESTRICT macro
#define RESTRICT __restrict

#ifdef WIN32
  #define __thread __declspec(thread)
#endif


////////////////////////////////////////////////////////////////////////////////
// Global data used during simulation (tables, constants, etc. that are setup
//   during startup but essentially read-only and constant throughout the
//   steady-state simulation.

__attribute__((target(mic))) double * device__table_four = NULL;
__attribute__((target(mic))) float * device__table_four_float = NULL;
__attribute__((target(mic))) int device__table_four_n_16 = 0;

__attribute__((target(mic))) double * device__lj_table = NULL;
__attribute__((target(mic))) float * device__lj_table_float = NULL;
__attribute__((target(mic))) int device__lj_table_dim = 0;
__attribute__((target(mic))) int device__lj_table_size = 0;

__attribute__((target(mic))) unsigned int * device__exclusion_bits = NULL;
__attribute__((target(mic))) long int device__exclusion_bits_size = 0;

__attribute__((target(mic))) mic_constants * device__constants = NULL;

__attribute__((target(mic))) double * device__table_four_copy = NULL;
__attribute__((target(mic))) float * device__table_four_float_copy = NULL;
__attribute__((target(mic))) double * device__lj_table_copy = NULL;
__attribute__((target(mic))) float * device__lj_table_float_copy = NULL;
__attribute__((target(mic))) unsigned int * device__exclusion_bits_copy = NULL;
__attribute__((target(mic))) mic_constants * device__constants_copy = NULL;

__attribute__((target(mic))) patch_pair* device__patch_pairs_copy = NULL;
__attribute__((target(mic))) force_list* device__force_lists_copy = NULL;
__attribute__((target(mic))) atom* device__atoms_copy = NULL;
__attribute__((target(mic))) atom_param* device__atom_params_copy = NULL;
__attribute__((target(mic))) int device__patch_pairs_copy_size = 0;
__attribute__((target(mic))) int device__force_lists_copy_size = 0;
__attribute__((target(mic))) int device__atoms_copy_size = 0;
__attribute__((target(mic))) int device__atom_params_copy_size = 0;


////////////////////////////////////////////////////////////////////////////////
// Device variables which exist both on the host and the MIC device and/or are
//   used to manage moving data betwen the host and the MIC device.

// DMK - DEBUG - PE and node info for printing debug output
__thread int host__pe = -1;
__thread int host__node = -1;
__attribute__((target(mic))) int device__pe = -1;
__attribute__((target(mic))) int device__node = -1;

__thread int singleKernelFlag = 0;

__thread patch_pair * host__patch_pairs = NULL;
__thread int host__patch_pairs_size = 0;
__thread int host__patch_pairs_bufSize = 0;
__attribute__((target(mic))) const patch_pair * device__patch_pairs = NULL;
__attribute__((target(mic))) int device__patch_pairs_size = 0;

__thread force_list * host__force_lists = NULL;
__thread int host__force_lists_size = 0;
__thread int host__force_lists_bufSize = 0;
__attribute__((target(mic))) force_list * device__force_lists = NULL;
__attribute__((target(mic))) int device__force_lists_size = 0;

__attribute__((target(mic))) uintptr_t device__pairlists = 0;  // NOTE: Don't use a global in case the device is shared between multiple threads, so each thread's patch pair lists has its own pairlist pointer
__attribute__((target(mic))) int device__pairlists_alloc_size = 0;

__attribute__((target(mic))) int** device__pl_array = NULL;
__attribute__((target(mic))) int* device__pl_size = NULL;
__attribute__((target(mic))) double** device__r2_array = NULL;

__thread atom * host__atoms = NULL;
__thread atom_param * host__atom_params = NULL;
__thread double4 * host__forces = NULL;
__thread double4 * host__slow_forces = NULL;
__thread int host__atoms_size = 0;
__thread int host__atoms_bufSize = 0;
__attribute__((target(mic))) atom * device__atoms = NULL;
__attribute__((target(mic))) atom_param * device__atom_params = NULL;
__attribute__((target(mic))) double4 * device__forces = NULL;
__attribute__((target(mic))) double4 * device__slow_forces = NULL;
__attribute__((target(mic))) int device__atoms_size = 0;

__thread size_t host__force_buffers_req_size = 0;
__attribute__((target(mic))) size_t device__force_buffers_req_size = 0;
__attribute__((target(mic))) size_t device__force_buffers_alloc_size = 0;
__attribute__((target(mic))) double4 * device__force_buffers = NULL;
__attribute__((target(mic))) double4 * device__slow_force_buffers = NULL;

__attribute__((target(mic))) mic_position3_t device__lata;
__attribute__((target(mic))) mic_position3_t device__latb;
__attribute__((target(mic))) mic_position3_t device__latc;

__thread mic_kernel_data * host__kernel_data = NULL;

__attribute__((target(mic))) int device__numOMPThreads = -1;

__thread int tag_atom_params;
__thread int tag_remote_kernel;
__thread int tag_local_kernel;

__thread int patch_pairs_copySize;
__thread int force_lists_copySize;
__thread int atom_params_copySize;

#if MIC_DEVICE_FPRINTF != 0
__attribute__((target(mic))) FILE * device__fout = NULL;
#endif

// DMK - DEBUG - Based on the number of kernels invoked, track the timestep number.
//__thread int host__timestep = 0;
__attribute__((target(mic))) int device__timestep = 0;


// DMK - TRACING / TIMING
#include <sys/time.h>
__declspec(target(mic))
double getCurrentTime() {
  timeval now;
  gettimeofday(&now, NULL);
  return (double)(now.tv_sec + (now.tv_usec * 1.0e-6));
}

// DMK - TRACING
#if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
  #if MIC_DEVICE_TRACING_DETAILED != 0
    __thread double * host__device_times_computes = NULL;
    __thread double * host__device_times_patches = NULL;
    __attribute__((target(mic))) double * device__device_times_computes = NULL;
    __attribute__((target(mic))) double * device__device_times_patches = NULL;
  #endif
  __thread double * host__device_times_start = NULL;
  __attribute__((target(mic))) double * device__device_times_start = NULL;
#endif


// DMK - TODO : A device version of the die function, which does not have access
//   to the host info.  Only to be used within the kernel itself.  Perhaps add
//   some startup code to send host/PE# and device# to device and then target
//   mic_die to the device as well, removing the need for this separate function.
__attribute__((target(mic)))
void mic_dev_die(const char * const str) {
  #ifdef __MIC__
    const char * const loc = "on device";
  #else
    const char * const loc = "on host";
  #endif
  if (str != NULL) {
    printf("[MIC_DIE] :: \"%s\" (%s)\n", str, loc);
  } else {
    printf("[MIC_DIE] :: mic_dev_die called (%s)\n", loc);
  }
  fflush(NULL);
  abort();
}


__attribute__((target(mic)))
void mic_print_config() {

  #if MIC_PRINT_CONFIG != 0

    // DMK - TODO | FIXME : Create a mechanism, so that it is only printed once if
    //   there are multiple MIC devices being used.
  
    printf("device :: MULTIPLE_THREADS  (%d)\n", MULTIPLE_THREADS);
    printf("device :: # OpenMP Threads : %d\n", device__numOMPThreads);

    printf("device :: MIC_HANDCODE_FORCE                     (%d)\n", MIC_HANDCODE_FORCE);
    printf("device ::   MIC_HANDCODE_FORCE_PFDIST            (%d)\n", MIC_HANDCODE_FORCE_PFDIST);
    printf("device ::   MIC_HANDCODE_FORCE_USEGATHER_NBTBL   (%d)\n", MIC_HANDCODE_FORCE_USEGATHER_NBTBL);
    printf("device ::   MIC_HANDCODE_FORCE_USEGATHER         (%d)\n", MIC_HANDCODE_FORCE_USEGATHER);
    printf("device ::   MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT  (%d)\n", MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT);
    printf("device ::   MIC_HANDCODE_FORCE_COMBINE_FORCES    (%d)\n", MIC_HANDCODE_FORCE_COMBINE_FORCES);

    printf("device :: MIC_HANDCODE_FORCE_SINGLE            (%d)\n", MIC_HANDCODE_FORCE_SINGLE);
    printf("device :: MIC_HANDCODE_FORCE_CALCR2TABLE       (%d)\n", MIC_HANDCODE_FORCE_CALCR2TABLE);
    printf("device :: MIC_HANDCODE_FORCE_SOA_VS_AOS        (%d)\n", MIC_HANDCODE_FORCE_SOA_VS_AOS);

    printf("device :: MIC_SORT_ATOMS  (%d)\n", MIC_SORT_ATOMS);
    printf("device :: MIC_SORT_LISTS  (%d)\n", MIC_SORT_LISTS);

    printf("device :: MIC_HANDCODE_PLGEN    (%d)\n", MIC_HANDCODE_PLGEN);
    printf("device :: MIC_TILE_PLGEN        (%d)\n", MIC_TILE_PLGEN);
    printf("device :: MIC_CONDITION_NORMAL  (%d)\n", MIC_CONDITION_NORMAL);
    printf("device :: MIC_PAD_PLGEN         (%d)\n", MIC_PAD_PLGEN);

    printf("device :: MIC_SORT_COMPUTES  (%d)\n", MIC_SORT_COMPUTES);

    printf("device :: MIC_SYNC_INPUT  (%d)\n", MIC_SYNC_INPUT);

    printf("device :: REFINE_PAIRLISTS            (%d)\n", REFINE_PAIRLISTS);
    printf("device ::   REFINE_PAIRLIST_HANDCODE  (%d)\n", REFINE_PAIRLIST_HANDCODE);
    printf("device ::   REFINE_PAIRLISTS_XYZ      (%d)\n", REFINE_PAIRLISTS_XYZ);

    printf("device :: MIC_FULL_CHECK          (%d)\n", MIC_FULL_CHECK);
    printf("device :: MIC_EXCL_CHECKSUM_FULL  (%d)\n", MIC_EXCL_CHECKSUM_FULL);

    printf("device :: MIC_ALIGN  (%d)\n", MIC_ALIGN);

    fflush(NULL);

  #endif // MIC_PRINT_CONFIG
}


// Function to initialize a given MIC device by initializing various device
//   variables, arrays, etc.
// Input:
//  - pe : the PE number for the host core associated with the device
//  - deviceNum : the device number to be initialized
void mic_init_device(const int pe, const int node, const int deviceNum) {

  // Record the PE and node info associated with the given device
  host__pe = pe;
  host__node = node;

  // Initialize kernel data structures
  host__kernel_data = (mic_kernel_data*)(_MM_MALLOC_WRAPPER(2 * sizeof(mic_kernel_data), 64, "mic_kernel_data"));
  __ASSERT(host__kernel_data != NULL);
  mic_kernel_data * kernel_data = host__kernel_data;

  // Initialize flags for offload buffers that are not copied every timestep
  patch_pairs_copySize = 0;
  force_lists_copySize = 0;
  atom_params_copySize = 0;

  // Allocate data for the host__device_times_start buffer (constant size throughout simulation)
  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    host__device_times_start = (double*)_MM_MALLOC_WRAPPER(10 * sizeof(double), 64, "host__device_times_start");
    __ASSERT(host__device_times_start != NULL);
    double * device_times_start = host__device_times_start;
    #define DEVICE_TIMES_CLAUSE  nocopy(device_times_start[0:10] : alloc_if(1) free_if(0) align(64)) \
                                 nocopy(device__device_times_start)
  #else
    #define DEVICE_TIMES_CLAUSE
  #endif

  // Initialize the device itself via an offload pragma section
  #pragma offload target(mic:deviceNum) \
    in(pe) in(node) nocopy(device__pe) nocopy(device__node) \
    in(kernel_data[0:2] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
    nocopy(device__pl_array) nocopy(device__pl_size) nocopy(device__r2_array) \
    nocopy(device__numOMPThreads) \
    DEVICE_TIMES_CLAUSE \
    DEVICE_FPRINTF_CLAUSE
  {
    __MUST_BE_MIC;
    __FULL_CHECK(__ASSERT(node >= 0 && pe >= 0));

    device__pe = pe;
    device__node = node;

    #if MIC_DEVICE_FPRINTF != 0
    {
      char filename[128] = { 0 };
      sprintf(filename, "/tmp/namd_deviceDebugInfo.%d", device__pe);
      if (device__node <= 0) {
        printf("[MIC-DEBUG] :: Generating debug output to device file for MICs.\n");
        fflush(NULL);
      }
      device__fout = fopen(filename, "w");
    }
    #endif

    DEVICE_FPRINTF("Device on PE %d (node: %d) initializing...\n", device__pe, device__node);

    // Get the number of threads available to the code
    #if MULTIPLE_THREADS != 0
      device__numOMPThreads = omp_get_max_threads();
    #else
      device__numOMPThreads = 1;
    #endif

    // Initialize the r2 arrays (scratch buffers used by computes to refine pairlists
    //   each timestep, when pairlist refinement is enabled)
    #if REFINE_PAIRLISTS != 0
      device__pl_array = (int**)(_MM_MALLOC_WRAPPER(device__numOMPThreads * sizeof(int*), 64, "device__pl_array"));       __ASSERT(device__pl_array != NULL);
      device__pl_size = (int*)(_MM_MALLOC_WRAPPER(device__numOMPThreads * sizeof(int), 64, "device__pl_size"));           __ASSERT(device__pl_size != NULL);
      device__r2_array = (double**)(_MM_MALLOC_WRAPPER(device__numOMPThreads * sizeof(double*), 64, "device__r2_array");  __ASSERT(device__r2_array != NULL);
      for (int i = 0; i < device__numOMPThreads; i++) {
        device__pl_array[i] = NULL;
        device__pl_size[i] = 0;
        device__r2_array[i] = NULL;
      }
    #endif

    // Copy the pointer to the device_times_start buffer
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
      device__device_times_start = device_times_start;
    #endif
    
    // DMK : NOTE | TODO | FIXME - This will print the configuration of the MICs, but the condition
    //   assumes that there will be at least one MIC on the first node (node 0).  Further, it assumes
    //   that all MICs will be the same (i.e. only prints on one device).  If there are cases where
    //   these assumptions do not hold, we need to expand this code to cover those cases.  For now,
    //   leaving it this way because it reduces the amount of output during the run (espeically when
    //   scaling to multiple nodes and/or multiple MICs per node) and is likely the case (for now).
    if (node == 0 && deviceNum == 0) { mic_print_config(); }

  } // end pragma offload

  #undef DEVICE_TIMES_CLAUSE
}


// Function that returns 0 if executed on the host, 1 if executed on the MIC device.  Used to
//   indicate if a target is available or not during application startup.
// Input: N/A
// Output:
//  - 0 if executed on host, 1 if executed on device
__attribute__((target(mic))) int mic_check_internal() {
    int retval;
    #ifdef __MIC__
        retval = 1;
    #else
        retval = 0;
    #endif
    return retval;
}


// Function to check that the given device is available for offloading.
// Input:
//  - dev : The device number (0+) to check
// Output:
//  - 1 if device is available (will call mic_dev_die if not available to prevent further application progress)
int mic_check(int deviceNum) {
  int dev_ok = 0;
  #pragma offload target(mic:deviceNum) inout(dev_ok)
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_check :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif
    dev_ok = mic_check_internal();
  }
  return dev_ok;
}


// DMK - NOTE : The following is debug code used for verifying data structures.  Should not be used in
//   production runs.
#if MIC_DATA_STRUCT_VERIFY != 0


__declspec(target(mic))
void _verify_parallel_memcpy(void * const dst,
                             const void * const src,
                             const size_t size
                            ) {
  const char* const src_c = (const char* const)src;
  char* const dst_c = (char* const)dst;
  #pragma omp parallel for schedule(static, 4096)
  for (int i = 0; i < size; i++) { dst_c[i] = src_c[i]; }
}

template<class T>
__declspec(target(mic))
void* _verify_remalloc(int &copySize, const int size, T* &ptr) {
  if (copySize < size) {
    copySize = (int)(1.2f * size);  // Add some extra buffer room
    copySize = (copySize + 1023) & (~1023);  // Round up to multiple of 1024
    _MM_FREE_WRAPPER(ptr);
    ptr = (T*)_MM_MALLOC_WRAPPER(copySize * sizeof(T), 64, "_verify_remalloc");
    __FULL_CHECK(assert(ptr != NULL));
  }
  return ptr;
}


__declspec(target(mic))
int _verify_pairlist(const int pe, const int timestep, const int isRemote, const int phase,
                     const int i_upper, const int j_upper,
                     const int * const pairlist_base,
                     const int ppI, const char * const plStr,
                     const patch_pair &pp, const int plTypeIndex
                    ) {

  const atom_param * const pExt_0 = device__atom_params + pp.patch1_atom_start;
  const atom_param * const pExt_1 = device__atom_params + pp.patch2_atom_start;

  if (pairlist_base == NULL) {
    if (timestep > 0) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: pairlist_base == NULL !!!\n", pe, timestep, isRemote, phase); return 1; }
  } else if (i_upper <= 0) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: i_upper <= 0 !!!\n", pe, timestep, isRemote, phase); return 1;
  } else if (j_upper <= 0) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: j_upper <= 0 !!!\n", pe, timestep, isRemote, phase); return 1;
  } else {

    int currentI = 0;
    int inPadding = 0;
    int pairlist_size = pairlist_base[1] - 2;
    int pairlist_allocSize = pairlist_base[0];
    const int * const pairlist = pairlist_base + 2;

    if (pairlist_size + 2 > pairlist_allocSize) {
      printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: pairlist_size + 2 > pairlist_allocSize !!!\n", pe, timestep, isRemote, phase); return 1;
    }

    for (int k = 0; k < pairlist_size; k++) {

      int somethingWrong = 0;
      int i = (pairlist[k] >> 16) & 0xFFFF;
      int j =  pairlist[k]        & 0xFFFF;

      // Verify the interaction is in the correct list
      if (j != 0xFFFF) { // If this is a valid entry, not padding
        const int maxDiff = pExt_1[j].excl_maxdiff;
        const int indexDiff = pExt_0[i].index - pExt_1[j].index;
        if (indexDiff < -1 * maxDiff || indexDiff > maxDiff) { // Not in patten, so must be NORM
          if (plTypeIndex != PL_NORM_INDEX) {
            printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: interaction outside maxDiff range in non-NORM pairlist !!!"
                   " - ppI:%d, i:%d, j:%d, indexDiff:%d, maxDiff:%d, plTypeIndex:%d\n",
                   pe, timestep, isRemote, phase,
                   ppI, i, j, indexDiff, maxDiff, plTypeIndex
                  );
            somethingWrong = 1;
	  }
        } else { // In pattern, so verify against bits
          const int offset = (2 * indexDiff) + pExt_1[j].excl_index;
          const int offset_major = offset / (sizeof(unsigned int) * 8);
          const int offset_minor = offset % (sizeof(unsigned int) * 8);
          const int exclFlags = ((device__exclusion_bits[offset_major]) >> offset_minor) & 0x03;
          if (exclFlags == 0x3) { 
            printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: invalid exclusion bit pattern detected !!!"
                   " - ppI:%d, i:%d, j:%d, indexDiff:%d, maxDiff:%d, plTypeIndex:%d\n",
                   pe, timestep, isRemote, phase,
                   ppI, i, j, indexDiff, maxDiff, plTypeIndex
                  );
            somethingWrong = 1;
	  } else if (exclFlags != plTypeIndex) { 
            printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: plTypeIndex/exclution_bits mismatch !!!"
                   " - i:%d, j:%d, indexDiff:%d, maxDiff:%d, plTypeIndex:%d, exclFlags:%d\n",
                   pe, timestep, isRemote, phase,
                   i, j, indexDiff, maxDiff, plTypeIndex, exclFlags
                  );
            somethingWrong = 1;
	  }
        }
      }

      // Check for a new "i" value
      if (i != currentI) {
        if (i < currentI) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: i < currentI !!!\n", pe, timestep, isRemote, phase); somethingWrong = 1; }
        if (k % 16 != 0) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: i changed when k % 16 != 0 !!!\n", pe, timestep, isRemote, phase); somethingWrong = 1; }
        currentI = i;
        inPadding = 0;  // New i means we are transitioning into a non-padding region
      }

      // If we are in a padding region...
      if (inPadding) {
        if (i < 0 || i >= i_upper) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: invalid i value detected in padded region!!!\n", pe, timestep, isRemote, phase); somethingWrong = 1; }
        if (j != 0xFFFF) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: invalid j value detected in padded region !!!\n", pe, timestep, isRemote, phase); somethingWrong = 1; }

      // Otherwise, we are in a non-padding region...
      } else {
        if (i < 0 || i >= i_upper) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: invalid i value detected !!!\n", pe, timestep, isRemote, phase); somethingWrong = 1; }
        if (j == 0xFFFF) {  // Detect transition into a padded region
          inPadding = 1;
	} else {
          if (j < 0 || j >= j_upper) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: invalid j value detected !!!\n", pe, timestep, isRemote, phase); somethingWrong = 1; }
	}
      }

      if (somethingWrong) {
        char buf[2048];
        int lo = k - 3; if (lo < 0) { lo = 0; }
        int hi = k + 3; if (hi >= pairlist_size) { hi = pairlist_size; }
        char *str = buf;
        str += sprintf(str, "[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d ::   pl.%s[%d:%d->%d,%d:%d] = { ... ", pe, timestep, isRemote, phase, plStr, ppI, lo, hi, k, pairlist_size);
        for (int _k = lo; _k < hi; _k++) { str += sprintf(str, "%s0x%08x%s ", ((_k==k)?("<<"):("")), pairlist[_k], ((_k==k)?(">>"):(""))); }
        str += sprintf(str, "}, i_upper:%d, j_upper:%d\n", i_upper, j_upper);
        printf(buf);
      }
    }
  }

  return 0;
}


__declspec(target(mic))
int _verify_forces(const int pe, const int timestep, const int isRemote, const int phase,
                   const int index, const double4 * const forceBuffers, const double4 * const forces,
                   const char * const typeStr
                  ) {

  force_list &fl = device__force_lists[index];

  // Verify the indexes / sizes
  if (fl.force_list_start < 0 || fl.force_list_start >= device__force_buffers_req_size) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__force_lists[%d].force_list_start invalid (%d, size:%lu) !!!\n", pe, timestep, isRemote, phase, index, fl.force_list_start, device__force_buffers_req_size); return 1; }
  if (fl.force_output_start < 0 || fl.force_output_start >= device__atoms_size) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__force_lists[%d].force_output_start invalid (%d, size:%d) !!!\n", pe, timestep, isRemote, phase, index, fl.force_output_start, device__atoms_size); return 1; }
  if (fl.patch_stride < 0) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__force_lists[%d].patch_stride(%d) < 0 !!!\n", pe, timestep, isRemote, phase, index, fl.patch_stride); return 1; }
  if (fl.force_list_size < 0) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__force_lists[%d].force_list_size(%d) < 0 !!!\n", pe, timestep, isRemote, phase, index, fl.force_list_size); return 1; }
  if (fl.force_list_start + (4 * fl.patch_stride * fl.force_list_size) > device__force_buffers_req_size) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: fl[%d].force_list_start(%d) + (4 * fl.patch_stride(%d) * fl.force_list_size(%d)) > device__force_buffers_req_size(%lu) !!!\n",
           pe, timestep, isRemote, phase, index, fl.force_list_start, fl.patch_stride, fl.force_list_size, device__force_buffers_req_size);
    return 1;
  }

  // Check sub-lists for large forces
  double* f = (double*)(forceBuffers + fl.force_list_start);
  for (int i = 0; i < 4 * fl.patch_stride * fl.force_list_size; i++) {
    if (fabsf(f[i]) > 2.5e2) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: Large %sforce detected (1) !!! f[%d]:%le\n", pe, timestep, isRemote, phase, typeStr, i, f[i]); } // NOTE: Want to print them all, so don't return here
  }

  // Check final list for large forces
  f = (double*)(forces + fl.force_output_start);
  for (int i = 0; i < 4 * fl.patch_stride; i++) {
    if (fabsf(f[i]) > 2.5e2) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: Large %sforce detected (2) !!! f[%d]:%le\n", pe, timestep, isRemote, phase, typeStr, i, f[i]); } // NOTE: Want to print them all, so don't return here
  }

  return 0;
}

__declspec(target(mic))
int _verify_buffers_match(const int pe, const int timestep, const int isRemote, const int phase,
                          const char * const buf, const char * const golden, const unsigned int bufSize,
                          const char * const str
                         ) {
  if (buf == NULL) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: buf == NULL (%s) !!!\n", pe, timestep, isRemote, phase, str); return 1;
  } else if (golden == NULL) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: golden == NULL (%s) !!!\n", pe, timestep, isRemote, phase, str); return 1;
  } else if (bufSize < 0) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: bufSize <= 0 (%s) !!!\n", pe, timestep, isRemote, phase, str); return 1;
  } else {
    int mismatchCount = 0;
    #pragma omp parallel for reduction(+:mismatchCount)
    for (int i = 0; i < bufSize; i++) {
      #if 1
        if (buf[i] != golden[i]) {
          printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: buf[%d/%d]:0x%02x != golden[%d]:0x%02x (%s) !!!\n",
                 pe, timestep, isRemote, phase, i, bufSize, buf[i], i, golden[i], str
                );
          mismatchCount++;
	}
      #else
        if (start < 0) {
          if (buf[i] != golden[i]) {
            printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: buf[%d/%d]:0x%02x != golden[%d]:0x%02x (%s) !!!\n",
                   pe, timestep, isRemote, phase, i, bufSize, buf[i], i, golden[i], str
                  );
            start = i;
            mismatchFound = 1;
	  }
        } else {
          if (buf[i] == golden[i] || i == bufSize - 1) {
            if (i < bufSize - 1 && buf[i+1] != golden[i+1]) { // NOTE: Ignore single byte matches found within consecutive mismatches
              printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: buf[%d/%d]:0x%02x != golden[%d]:0x%02x (%s) !!!\n",
                     pe, timestep, isRemote, phase, i, bufSize, buf[i], i, golden[i], str
                    );
            } else {
              printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: buf mismatch range %d -> %d (%s) !!!\n",
                     pe, timestep, isRemote, phase, start, ((buf[i] != golden[i] && i == bufSize - 1) ? (i+1) : (i)), str
                    );
              start = -1;
            }
          } else if (buf[i] != golden[i] && i - start < 16) {
            printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: buf[%d/%d]:0x%02x != golden[%d]:0x%02x (%s) !!!\n",
                   pe, timestep, isRemote, phase, i, bufSize, buf[i], i, golden[i], str
                  );
          }
        }
      #endif
    }
    if (mismatchCount > 0) { return 1; }
  }
  return 0;
}

__declspec(target(mic))
void _verify_tables(const int pe, const int timestep, const int isRemote, const int phase) {

  // Verify that certain tables/buffers haven't been modified
  if (_verify_buffers_match(pe, timestep, isRemote, phase, (char*)(device__exclusion_bits), (char*)(device__exclusion_bits_copy), device__exclusion_bits_size * sizeof(unsigned int), "exclusion_bits")) { return; }
  if (_verify_buffers_match(pe, timestep, isRemote, phase, (char*)(device__constants), (char*)(device__constants_copy), sizeof(mic_constants), "constants")) { return; }
  #if MIC_HANDCODE_FORCE_SINGLE != 0
    if (_verify_buffers_match(pe, timestep, isRemote, phase, (char*)(device__table_four_float), (char*)(device__table_four_float_copy), 61 * device__table_four_n_16 * sizeof(float), "table_four(float)")) { return; }
    if (_verify_buffers_match(pe, timestep, isRemote, phase, (char*)(device__lj_table_float), (char*)(device__lj_table_float_copy), device__lj_table_size / sizeof(double) * sizeof(float), "lj_table(float)")) { return; }
  #else
    if (_verify_buffers_match(pe, timestep, isRemote, phase, (char*)(device__table_four), (char*)(device__table_four_copy), 61 * device__table_four_n_16 * sizeof(double), "table_four(double)")) { return; }
    if (_verify_buffers_match(pe, timestep, isRemote, phase, (char*)(device__lj_table), (char*)(device__lj_table_copy), device__lj_table_size * sizeof(double), "lj_table(double)")) { return; }
  #endif

  // Check to see if the constants struct has been placed inside the exclusion_bits buffer
  // NOTE: This is a specific check related to a specific error that was occuring (leaving in for now)
  char *p1 = (char*)(device__constants);
  char *p2 = (char*)(device__exclusion_bits);
  int pDiff = p1 - p2;
  if (p1 >= p2 && p1 < p2 + (device__exclusion_bits_size * sizeof(unsigned int))) {
    printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: !!!! OVERLAP OF EXCLUSIONS AND CONSTANTS DETECTED !!!! - constants:%p(%ld), exclusions:%p(%ld)\n",
           pe, timestep, isRemote, phase, device__constants, sizeof(device__constants), device__exclusion_bits, device__exclusion_bits_size * sizeof(unsigned int)
          );
  }

  fflush(NULL);
}


__declspec(target(mic))
void _verify_data_structures(const int pe, const int timestep, const int isRemote, const int phase,
                             const int check_lo, const int check_hi, const int doSlow
                            ) {

  // Verify that the basic arrays are allocated
  if (device__atoms == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__atoms == NULL !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__atom_params == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__atom_params == NULL !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__forces == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__forces == NULL !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__slow_forces == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__slow_forces == NULL !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__patch_pairs == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__patch_pairs == NULL !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__force_lists == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__force_lists == NULL !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__pairlists == NULL) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__pairlists == NULL !!!\n", pe, timestep, isRemote, phase); return; }

  // Verify sizes (as much as possible)
  if (device__patch_pairs_size <= 0) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__patch_pairs_size <= 0 !!!\n", pe, timestep, isRemote, phase); return; }
  if (device__force_lists_size <= 0) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__force_lists_size <= 0 !!!\n", pe, timestep, isRemote, phase); return; }

  // Verify the tables
  _verify_tables(pe, timestep, isRemote, phase);

  // For each patch pair
  #if MULTIPLE_THREADS != 0
    #pragma omp parallel for
  #endif
  for (int k = check_lo; k < check_hi; k++) {

    const patch_pair &pp = device__patch_pairs[k];

    // Check indexing
    if (pp.patch1_atom_start < 0 || pp.patch1_atom_start >= device__atoms_size) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__patch_pairs[%d].patch1_atom_start invalid (%d, size:%d) !!!\n", pe, timestep, isRemote, phase, k, pp.patch1_atom_start, device__atoms_size); }
    if (pp.patch2_atom_start < 0 || pp.patch2_atom_start >= device__atoms_size) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__patch_pairs[%d].patch2_atom_start invalid (%d, size:%d) !!!\n", pe, timestep, isRemote, phase, k, pp.patch2_atom_start, device__atoms_size); }
    if (pp.patch1_force_start < 0 || pp.patch1_force_start >= device__force_buffers_req_size) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__patch_pairs[%d].patch1_force_start invalid (%d, size:%lu) !!!\n", pe, timestep, isRemote, phase, k, pp.patch1_force_start, device__force_buffers_req_size); }
    if (pp.patch2_force_start < 0 || pp.patch2_force_start >= device__force_buffers_req_size) { printf("[VERIFICATION-ERROR] :: PE:%d.%d.%d.%d :: device__patch_pairs[%d].patch2_force_start invalid (%d, size:%lu) !!!\n", pe, timestep, isRemote, phase, k, pp.patch2_force_start, device__force_buffers_req_size); }

    // Verify the pairlists
    int i_upper = pp.patch1_size;
    int j_upper = pp.patch2_size;
    const int ** const pairlists = (const int ** const)device__pairlists;
    _verify_pairlist(pe, timestep, isRemote, phase, i_upper, j_upper, pairlists[NUM_PAIRLIST_TYPES * k + PL_NORM_INDEX], k, "NORM", pp, PL_NORM_INDEX);
    _verify_pairlist(pe, timestep, isRemote, phase, i_upper, j_upper, pairlists[NUM_PAIRLIST_TYPES * k + PL_EXCL_INDEX], k, "EXCL", pp, PL_EXCL_INDEX);
    _verify_pairlist(pe, timestep, isRemote, phase, i_upper, j_upper, pairlists[NUM_PAIRLIST_TYPES * k +  PL_MOD_INDEX], k,  "MOD", pp,  PL_MOD_INDEX);
  }

  // Verify forces if this is the local kernel (end of timestep)
  if (isRemote == 0) {
    #if MULTIPLE_THREADS != 0
      #pragma omp parallel for
    #endif
    for (int k = 0; k < device__force_lists_size; k++) {
      _verify_forces(pe, timestep, isRemote, phase, k, device__force_buffers, device__forces, "");
      if (doSlow) { _verify_forces(pe, timestep, isRemote, phase, k, device__slow_force_buffers, device__slow_forces, "slow "); }
    }
  }

  fflush(NULL);
}

#endif // MIC_DATA_STRUCT_VERIFY



// Function called during startup to push the table_four data (lookup table used
//   during force computation) that has been calculated on the host down to the
//   given device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - table_four : pointer to the table data itself
//   - table_four_n : dimension of the table_four table (used for indexing and
//       size calculations)
// Output: N/A
void mic_bind_table_four(const int deviceNum,
                         const double *table_four,
                         const int table_four_n,
                         const int table_four_n_16
                        ) {

  // Verify parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));
  __FULL_CHECK(__ASSERT(table_four != NULL));
  __FULL_CHECK(__ASSERT(table_four_n > 0));
  __FULL_CHECK(__ASSERT(table_four_n_16 > 0 && table_four_n_16 >= table_four_n));
  __FULL_CHECK(__ASSERT(table_four_n_16 % 16 == 0));

  // Copy the table pointer and dimension information into the device variables.
  //   Note that there are actually several sub-tables within the table_four
  //   table, each a multiple of table_four_n elements in length.
  int numTableFourElems = 61 * table_four_n_16;  // WARNING !!! Must match ComputeNonbondedUtil.C (const 61) !!!

  // Transfer the table_four data and dimension data to the given device
  #pragma offload target(mic:deviceNum) \
    in(table_four_n_16) nocopy(device__table_four_n_16) \
    in(numTableFourElems) \
    in(table_four[0:numTableFourElems] : alloc_if(1) free_if(1)) \
    nocopy(device__table_four) nocopy(device__table_four_float) \
    nocopy(device__table_four_copy) nocopy(device__table_four_float_copy)
  {
    __MUST_BE_MIC;

    // Copy and verify the table_four_n_16 value
    device__table_four_n_16 = table_four_n_16;
    __FULL_CHECK(__ASSERT(device__table_four_n_16 > 0));

    // Allocate memory on the device for the table_four table.  Depending on the precision being used
    //   for the table, either create a 64-bit or 32-bit copy of the data on the device.
    __FULL_CHECK(__ASSERT(device__table_four == NULL && device__table_four_float == NULL));  // Check for double allocation
    #if MIC_HANDCODE_FORCE_SINGLE != 0

      // Allocate and copy data into a 32-bit version of the table
      device__table_four_float = (float*)(_MM_MALLOC_WRAPPER(sizeof(float) * numTableFourElems, 64, "device__table_four_float"));
      __ASSERT(device__table_four_float != NULL);
      for (int i = 0; i < numTableFourElems; i++) {
        device__table_four_float[i] = (float)(table_four[i]);
      }

      // Create a second copy for verification purposes
      #if MIC_DATA_STRUCT_VERIFY != 0
        device__table_four_float_copy = (float*)(_MM_MALLOC_WRAPPER(sizeof(float) * numTableFourElems, 64, "device__table_four_float_copy"));
        __ASSERT(device__table_four_float_copy != NULL);
        memcpy(device__table_four_float_copy, device__table_four_float, numTableFourElems * sizeof(float));
      #endif

    #else // MIC_HANDCODE_FORCE_SINGLE != 0

      // Allocate and copy data into a 64-bit version of the table
      device__table_four = (double*)(_MM_MALLOC_WRAPPER(sizeof(double) * numTableFourElems, 64, "device__table_four"));
      __ASSERT(device__table_four != NULL);
      for (int i = 0; i < numTableFourElems; i++) {
        device__table_four[i] = table_four[i];
      }
      #if MIC_DATA_STRUCT_VERIFY != 0
        device__table_four_copy = (double*)(_MM_MALLOC_WRAPPER(sizeof(double) * numTableFourElems, 64, "device__table_four_copy"));
        __ASSERT(device__table_four_copy != NULL);
        memcpy(device__table_four_copy, device__table_four, numTableFourElems * sizeof(double));
      #endif

    #endif // MIC_HANDCODE_FORCE_SINGLE != 0

  } // end pragma offload
}


// Function called during startup to push the lj_table data (lookup table used
//   during force computation) that has been calculated on the host down to the
//   given device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - lj_table : pointer to the table data itself (64-bit version on the host)
//   - lj_table_dim : dimension of the lj_table (used for indexing)
//   - lj_table_size : size of the lj_table
// Output: N/A
void mic_bind_lj_table(const int deviceNum,
                       const char * lj_table,
                       const int lj_table_dim,
                       const int lj_table_size
                      ) {

  // Verify Parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));
  __FULL_CHECK(__ASSERT(lj_table != NULL));
  __FULL_CHECK(__ASSERT(lj_table_dim > 0));
  __FULL_CHECK(__ASSERT(lj_table_size > 0));

  // Transfer the table data, dimension, and size information to the given device.
  #pragma offload target(mic:deviceNum) \
    in(lj_table_dim) nocopy(device__lj_table_dim) \
    in(lj_table_size) nocopy(device__lj_table_size) \
    in(lj_table[0:lj_table_size] : alloc_if(1) free_if(1)) \
    nocopy(device__lj_table) nocopy(device__lj_table_float) \
    nocopy(device__lj_table_copy) nocopy(device__lj_table_float_copy) \
    nocopy(device__pe)
  {
    __MUST_BE_MIC;

    // Copy and verify the LJ table size and dimension info
    device__lj_table_dim = lj_table_dim;
    device__lj_table_size = lj_table_size;
    __FULL_CHECK(__ASSERT(device__lj_table_dim > 0 && device__lj_table_size > 0));

    // Allocate memory on the device for the LJ table.  Depending on the precision being used
    //   for the table, either create a 64-bit or 32-bit copy of the data on the device.
    __FULL_CHECK(__ASSERT(device__lj_table == NULL && device__lj_table_float == NULL));  // Check for double allocation
    int numElements = device__lj_table_size * sizeof(char) / sizeof(double);
    #if MIC_HANDCODE_FORCE_SINGLE != 0

      // Allocate and copy the data into a 32-bit version of the table
      device__lj_table_float = (float*)(_MM_MALLOC_WRAPPER(sizeof(float) * numElements, 64, "device__lj_table_float"));
      __ASSERT(device__lj_table_float != NULL);
      for (int i = 0; i < numElements; i++) {
        device__lj_table_float[i] = (float)(((double*)lj_table)[i]);
      }

      // Create a copy of the table for verification purposes later
      #if MIC_DATA_STRUCT_VERIFY != 0
        device__lj_table_float_copy = (float*)(_MM_MALLOC_WRAPPER(sizeof(float) * numElements, 64, "device__lj_table_float_copy"));
        __ASSERT(device__lj_table_float_copy != NULL);
        memcpy(device__lj_table_float_copy, device__lj_table_float, sizeof(float) * numElements);
      #endif

    #else // MIC_HANDCODE_FORCE_SINGLE != 0

      // Allocate and copy the data into a 64-bit version of the table
      device__lj_table = (double*)(_MM_MALLOC_WRAPPER(sizeof(double) * numElements, 64, "device__lj_table"));
      __ASSERT(device__lj_table_float != NULL);
      for (int i = 0; i < numElements; i++) {
        device__lj_table[i] = ((double*)lj_table)[i];
      }

      // Create a copy of the table for verificiation purposes later
      #if MIC_DATA_STRUCT_VERIFY != 0
        device__lj_table_copy = (double*)(_MM_MALLOC_WRAPPER(sizeof(double) * numElements, 64, "device__lj_table_copy"));
        __ASSERT(device__lj_table_copy != NULL);
        memcpy(device__lj_table_copy, device__lj_table, sizeof(double) * numElements);
      #endif

    #endif // MIC_HANDCODE_FORCE_SINGLE != 0

  } // end pragma offload
}


// Function called during startup to push the exclusion data (lookup table used
//   during pairlist generation) that has been calculated on the host down to the
//   given device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - exclusion_bits : pointer to the exclusion data itself
//   - exclusion_bits_size : size of the exclsion data
// Output: N/A
void mic_bind_exclusions(const int deviceNum,
                         unsigned int *exclusion_bits,
                         const long int exclusion_bits_size
                        ) {

  // Verify Parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));
  __FULL_CHECK(__ASSERT(exclusion_bits != NULL));
  __FULL_CHECK(__ASSERT(exclusion_bits_size > 0));

  // Transfer the exclusion data and size information down to the given device.
  #pragma offload target(mic:deviceNum) \
    in(exclusion_bits_size) nocopy(device__exclusion_bits_size) \
    in(exclusion_bits[0:exclusion_bits_size] : alloc_if(1) free_if(1)) \
    nocopy(device__exclusion_bits) nocopy(device__exclusion_bits_copy)
  {
    __MUST_BE_MIC;

    // Copy and verify the table size info
    device__exclusion_bits_size = exclusion_bits_size;
    __FULL_CHECK(__ASSERT(device__exclusion_bits_size > 0));

    // Create a copy of the exclusion bits buffer on the device
    __ASSERT(device__exclusion_bits == NULL);  // Check for double allocate
    device__exclusion_bits = (unsigned int *)_MM_MALLOC_WRAPPER(sizeof(unsigned int) * device__exclusion_bits_size, 64, "device__exclusion_bits");
    __ASSERT(device__exclusion_bits != NULL);
    memcpy(device__exclusion_bits, exclusion_bits, sizeof(unsigned int) * device__exclusion_bits_size);

    // Create a copy of the table for verification purposes later
    #if MIC_DATA_STRUCT_VERIFY != 0
      device__exclusion_bits_copy = (unsigned int*)(_MM_MALLOC_WRAPPER(device__exclusion_bits_size * sizeof(unsigned int), 64, " device__exclusion_bits_copy"));
      __ASSERT(device__exclusion_bits_copy != NULL);
      memcpy(device__exclusion_bits_copy, device__exclusion_bits, sizeof(unsigned int) * device__exclusion_bits_size);
    #else
      device__exclusion_bits_copy = NULL;
    #endif

  } // end pragma offload
}


// Function called during startup to push (and calculate) several variables that
//   are constant and used throughout the simulation.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   ... : remaining parameters are simply simulation constants
// Output: N/A
void mic_bind_constants(const int deviceNum,
                        const double cutoff2,
                        const double dielectric_1,
                        const double scaling,
                        const double scale14,
                        const double r2_delta,
                        const int r2_delta_exp,
                        const int commOnly
                       ) {
  
  // Verify Parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));

  // Create a mic_constants data structure and copy (or calculate) the
  //   various constants that are to be pushed to the device in to this
  //   data structure.
  mic_constants constants;
  constants.cutoff2 = cutoff2;
  constants.dielectric_1 = dielectric_1;
  constants.scaling = scaling;
  constants.scale14 = scale14;
  constants.modf_mod = 1.0 - scale14;
  constants.r2_delta = r2_delta;
  constants.r2_delta_exp = r2_delta_exp;
  constants.r2_delta_expc = 64 * (r2_delta_exp - 1023);
  constants.commOnly = commOnly;
  constants.singleKernelFlag = singleKernelFlag;

  // Transfer the constants to the given device
  #pragma offload target(mic:deviceNum) \
    in(constants) nocopy(device__constants) \
    nocopy(device__pe) nocopy(device__exclusion_bits) \
    nocopy(device__exclusion_bits_copy) nocopy(device__exclusion_bits_size) \
    nocopy(device__table_four_copy) nocopy(device__table_four_float_copy) \
    nocopy(device__lj_table_copy) nocopy(device__lj_table_float_copy) \
    nocopy(device__constants_copy)
  {
    __MUST_BE_MIC;

    // Create a copy of the constant info on the device
    __ASSERT(device__constants == NULL);  // Check for double allocation
    device__constants = (mic_constants*)_MM_MALLOC_WRAPPER(sizeof(mic_constants), 64, "device__constants");
    __ASSERT(device__constants != NULL);
    memcpy(device__constants, &constants, sizeof(mic_constants));

    // Correct r2_delta_expc for use with single precision if using mixed precision is
    //   being used on the device
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      device__constants->r2_delta_expc = 64 * (device__constants->r2_delta_exp - 127);
    #endif

    // Copy the data for verification purposes later
    #if MIC_DATA_STRUCT_VERIFY != 0
      device__constants_copy = (mic_constants*)_MM_MALLOC_WRAPPER(sizeof(mic_constants), 64, "device__constants_copy");
      __ASSERT(device__constants_copy != NULL);
      memcpy(device__constants_copy, device__constants, sizeof(mic_constants));
    #endif

  } // end pragma offload
}


// Function called after patch pairs have been modified on the host to push those
//   modifications to the MIC device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - patch_pairs : array containing individual patch_pair records (data to copy)
//   - patch_pairs_size : the length (in elements) of the patch_pairs array (valid/used elements)
//   - patch_pairs_bufSize : the actual length (in elements) of the patch_pairs array (allocated)
// Output: N/A
void mic_bind_patch_pairs_only(const int deviceNum,
                               const patch_pair *patch_pairs,
                               const int patch_pairs_size,
                               const int patch_pairs_bufSize
			      ) {

  // NOTE: This function does not actually transfer the patch pairs, it only (re)allocates
  //   the device buffer(s) associated with them.  See mic_nonbonded_compute() for transfer.

  // Verify parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));
  __FULL_CHECK(__ASSERT(patch_pairs != NULL));
  __FULL_CHECK(__ASSERT(patch_pairs_size > 0));
  __FULL_CHECK(__ASSERT(patch_pairs_bufSize > 0));
  __FULL_CHECK(__ASSERT(patch_pairs_bufSize >= patch_pairs_size));

  // Check if the buffer currently allocated on the device is too small.  If so,
  //   reallocate the buffer to be larger.
  if (patch_pairs_bufSize > host__patch_pairs_bufSize) {  // If buffer is too small, reallocate

    // Setup clause required to free device_times_computes buffer on the device
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
      double * device_times_computes = host__device_times_computes;
      #define TIMES_COMPUTES_FREE_CLAUSE  nocopy(device_times_computes : alloc_if(0) free_if(1))
    #else
      #define TIMES_COMPUTES_FREE_CLAUSE
    #endif

    // Free old buffers
    if (host__patch_pairs != NULL) {

      const patch_pair * _patch_pairs = host__patch_pairs;
      #pragma offload target(mic:deviceNum) \
        nocopy(_patch_pairs : alloc_if(0) free_if(1)) \
        TIMES_COMPUTES_FREE_CLAUSE
      { __MUST_BE_MIC; }
    }

    // (Re)allocate memory for device_times_computes buffer on the host, and setup clause for device allocation
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
      if (host__device_times_computes != NULL) { _MM_FREE_WRAPPER(host__device_times_computes); }
      host__device_times_computes = (double*)_MM_MALLOC_WRAPPER(2 * patch_pairs_bufSize * sizeof(double), 64, "host__device_times_computes");
      __ASSERT(host__device_times_computes != NULL);
      device_times_computes = host__device_times_computes;
      #define TIMES_COMPUTES_ALLOC_CLAUSE  nocopy(device_times_computes[0:2*patch_pairs_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
    #else
      #define TIMES_COMPUTES_ALLOC_CLAUSE
    #endif

    // Record the new host-side pointer and buffer size
    host__patch_pairs = (patch_pair*)patch_pairs;
    host__patch_pairs_bufSize = patch_pairs_bufSize;

    // Allocate new memory
    #pragma offload target(mic:deviceNum) \
      in(patch_pairs_size) nocopy(device__patch_pairs_size) \
      in(patch_pairs_bufSize) \
      nocopy(patch_pairs[0:patch_pairs_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      nocopy(device__pairlists) nocopy(device__pairlists_alloc_size) \
      nocopy(device__numOMPThreads) nocopy(device__pe) \
      TIMES_COMPUTES_ALLOC_CLAUSE
    {
      __MUST_BE_MIC;

      // Copy and verify the number of patch pairs
      device__patch_pairs_size = patch_pairs_size;
      __ASSERT(device__patch_pairs_size > 0);

      // Make sure there are enough pairlists pointers (one per pairlist type per patch_pair)
      const int numPairlists = NUM_PAIRLIST_TYPES * device__patch_pairs_size;
      if (numPairlists > device__pairlists_alloc_size) {
        int **new_pairlists = (int**)(_MM_MALLOC_WRAPPER(numPairlists * sizeof(int*), 64, "device__pairlists"));
        int initStart = 0;
        if (device__pairlists != 0) {
          int **old_pairlists = (int**)(device__pairlists);
          memcpy((void*)new_pairlists, (void*)old_pairlists, sizeof(int*) * device__pairlists_alloc_size);
	  _MM_FREE_WRAPPER(old_pairlists);
          initStart = device__pairlists_alloc_size;
        }
        for (int i = initStart; i < numPairlists; i++) { new_pairlists[i] = NULL; }
        device__pairlists = (uintptr_t)new_pairlists;
        device__pairlists_alloc_size = numPairlists;
      }

    } // end pragma offload

    #undef TIMES_COMPUTES_FREE_CLAUSE
    #undef TIMES_COMPUTES_ALLOC_CLAUSE
  }

  // Record the current 'used' length of the array and flag the array for data transfer
  host__patch_pairs_size = patch_pairs_size;
  patch_pairs_copySize = patch_pairs_size;
}


// Function called after force lists have been modified on the host to push those
//   modifications to the MIC device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - force_lists : array containing individual force_list records (data to copy)
//   - force_lists_size : the length (in elements) of the force_lists array (valid/used elements)
//   - force_lists_bufSize : the actual length (in elements) of the force_lists array (allocated)
// Output: N/A
void mic_bind_force_lists_only(const int deviceNum,
                               force_list *force_lists,
                               const int force_lists_size,
                               const int force_lists_bufSize
                              ) {

  // NOTE: This function does not actually transfer the force lists, it only (re)allocates
  //   the device buffer(s) associated with them.  See mic_nonbonded_compute() for transfer.

  // Verify parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));
  __FULL_CHECK(__ASSERT(force_lists != NULL));
  __FULL_CHECK(__ASSERT(force_lists_size > 0));
  __FULL_CHECK(__ASSERT(force_lists_bufSize > 0));
  __FULL_CHECK(__ASSERT(force_lists_bufSize >= force_lists_size));

  // Check if the buffer currently allocated on the device is too small.  If so, reallocate the buffer.
  if (force_lists_bufSize > host__force_lists_bufSize) {

    // Setup clause needed to free the device_times_patches buffer on the device
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
      double * device_times_patches = host__device_times_patches;
      #define TIMES_PATCHES_FREE_CLAUSE  nocopy(device_times_patches : alloc_if(0) free_if(1))
    #else
      #define TIMES_PATCHES_FREE_CLAUSE
    #endif

    // Free the old buffer
    if (host__force_lists != NULL) {
      const force_list * _force_lists = host__force_lists;
      #pragma offload target(mic:deviceNum) \
        nocopy(_force_lists : alloc_if(0) free_if(1)) \
        TIMES_PATCHES_FREE_CLAUSE
      { __MUST_BE_MIC; }
    }

    // (Re)allocate the device_times_patches buffer on the host and setup the clause required to allocate on the device
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
      if (host__device_times_patches != NULL) { _MM_FREE_WRAPPER(host__device_times_patches); }
      host__device_times_patches = (double*)_MM_MALLOC_WRAPPER(2 * force_lists_bufSize * sizeof(double), 64, "host__device_times_patches");
      __ASSERT(host__device_times_patches != NULL);
      device_times_patches = host__device_times_patches;
      #define TIMES_PATCHES_ALLOC_CLAUSE  nocopy(device_times_patches[0:2*force_lists_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
    #else
      #define TIMES_PATCHES_ALLOC_CLAUSE
    #endif

    // Record teh new host-side pointer and buffer size
    host__force_lists = force_lists;
    host__force_lists_bufSize = force_lists_bufSize;

    // Allocate the new memory
    #pragma offload target(mic:deviceNum) \
      in(force_lists_size) nocopy(device__force_lists_size) \
      nocopy(force_lists[0:force_lists_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      TIMES_PATCHES_ALLOC_CLAUSE
    {
      __MUST_BE_MIC;

      // Copy and verify the number of force lists
      device__force_lists_size = force_lists_size;
      __ASSERT(device__force_lists_size > 0);
    }

    #undef TIMES_PATCHES_FREE_CLAUSE
    #undef TIMES_PATCHES_ALLOC_CLAUSE

  } // end pragma offload

  // Record the current 'used' length of the array and flag the array for data transfer
  host__force_lists_size = force_lists_size;
  force_lists_copySize = force_lists_size;
}


// Function called after the atom list has been modified on the host to push those
//   modifications to the MIC device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - atoms : array containing individual atom records (data to copy)
//   - atom_params : array containing individual atom_param records (data to copy)
//   - forces : array that will contain output data from the device (data to copy back)
//   - slow_forces : array that will contain output data from the device (data to copy back)
//   - atoms_size : the length (in elements) of the specified arrays (valid/used elements)
//   - atoms_bufSize : the actual length (in elements) of the specified arrays (allocated)
// Output: N/A
void mic_bind_atoms_only(const int deviceNum,
                         atom *atoms,
                         atom_param *atom_params,
                         double4 *forces,
                         double4 *slow_forces,
                         const int atoms_size,
                         const int atoms_bufSize
			) {

  // NOTE: This function does not actually transfer the atom/force data, it only (re)allocates
  //   the device buffer(s) associated with them.  See mic_nonbonded_compute() for transfer.

  // Verify parameters
  __FULL_CHECK(__ASSERT(deviceNum >= 0));
  __FULL_CHECK(__ASSERT(atoms != NULL));
  __FULL_CHECK(__ASSERT(atom_params != NULL));
  __FULL_CHECK(__ASSERT(forces != NULL));
  __FULL_CHECK(__ASSERT(slow_forces != NULL));
  __FULL_CHECK(__ASSERT(atoms_size > 0));
  __FULL_CHECK(__ASSERT(atoms_bufSize > 0));
  __FULL_CHECK(__ASSERT(atoms_bufSize >= atoms_size));

  // Check if the buffers are large enough (reallocate if not)
  if (atoms_bufSize > host__atoms_bufSize) {

    // Free the old buffers
    if (host__atoms != NULL) {
      const atom * _atoms = host__atoms;
      const atom_param * _atom_params = host__atom_params;
      const double4 * _forces = host__forces;
      const double4 * _slow_forces = host__slow_forces;

      #pragma offload target(mic:deviceNum) \
        nocopy(_atoms : alloc_if(0) free_if(1)) \
        nocopy(_atom_params : alloc_if(0) free_if(1)) \
        nocopy(_forces : alloc_if(0) free_if(1)) \
        nocopy(_slow_forces : alloc_if(0) free_if(1))
      { __MUST_BE_MIC; }
    }

    // Copy the new buffer and size info
    host__atoms = atoms;
    host__atom_params = atom_params;
    host__forces = forces;
    host__slow_forces = slow_forces;
    host__atoms_size = atoms_size;
    host__atoms_bufSize = atoms_bufSize;

    // Allocate the new buffers
    #pragma offload target(mic:deviceNum) \
      in(atoms_size : into(device__atoms_size)) \
      nocopy(atoms[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      nocopy(atom_params[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      nocopy(forces[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      nocopy(slow_forces[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
    { __MUST_BE_MIC; }
  }

  // Signal the copy
  host__atoms_size = atoms_size;
  atom_params_copySize = atoms_size;
}


// Function called to allocate memory on the MIC device related to the individual and private
//   force buffers for each of the compute objects.  Note that this function simply allocates
//   the memory on the device, as there is no host-side equivalent buffer (i.e. device only).
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - force_buffers_size : the amount of memory to be allocated on the device for the buffer
void mic_bind_force_buffers_only(const int deviceNum,
                                 const size_t force_buffers_size
                                ) {

  // Check if the requested size is larger than the already allocated size.  If so, flag the buffer
  //   as needing to be reallocated by increasing the request size passed to the device (occurs later).
  if (host__force_buffers_req_size < force_buffers_size) {
    host__force_buffers_req_size = (int)(force_buffers_size * 1.1f);
    host__force_buffers_req_size = (host__force_buffers_req_size + 4095) & (~4095);  // Round up to 4K
  }
}


#if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0

void mic_submit_patch_data(const int deviceNum,
                           void * const hostPtr,
                           void* &hostPtr_prev,
                           const int transferBytes,
                           const int allocBytes_host,
                           int &allocBytes_device,
                           uint64_t &devicePtr,
                           void* &signal
                          ) {

  // NOTE: This function is called once per patch per timestep by one of the host threads
  // TODO : Since signals are thread specific, cannot use them here since the master does
  //   the checking for signal completions.  Figure out a way to make this asynchronous,
  //   but just do synchronous transfer for now to see what happens.

  // Check to see if the device's copy of the buffer needs to be reallocated
  // NOTE: If we are in a position where the buffer might grow in size, the contents
  //   need to be rewritten in full (i.e. atom migration has occured), so don't worry
  //   about the previous contents of the buffer
  int allocFlag = 0;
  if (hostPtr != hostPtr_prev || allocBytes_host > allocBytes_device) {

    // Free the old device buffer if it has been allocated previously
    if (devicePtr != 0) {
      char* __hostPtr_prev = (char*)hostPtr_prev;
      #pragma offload target(mic:deviceNum) \
        out(__hostPtr_prev[0:0] : alloc_if(0) free_if(1))
      { }
    }

    // Flag that allocation is needed
    allocFlag = 1;
  }

  // Transfer the data (allocating if needed)
  char* __hostPtr = (char*)hostPtr;
  if (allocFlag != 0) {

    uint64_t __devicePtr = 0;
    allocBytes_device = allocBytes_host;
    #pragma offload target(mic:deviceNum) \
      in(__hostPtr[0:allocBytes_host] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      out(__devicePtr)
    { __devicePtr = (uint64_t)__hostPtr; }

    // Record the pointer for use later
    devicePtr = __devicePtr;
    signal = NULL; // Synchronous, so no need for the host to wait

  } else {

    #pragma offload_transfer target(mic:deviceNum) \
      in(__hostPtr[0:transferBytes] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
      signal(hostPtr)
    { }

    signal = hostPtr;  // Asynchronous, so force the host to wait on the signal later
  }

  // Take not of the host buffer that was transfered
  hostPtr_prev = hostPtr;
}

#endif  // MIC_SUBMIT_ATOMS_ON_ARRIVAL


// Kernel bodies : Include CopmuteNonbondedMICKernelBase.h to declare each
//   version of the force computation functions (calc_pair, calc_self, etc.).
//   Each time the header file is included, different macros are set/unset
//   to control the exact contents of 'that version of the kernel.'

#define NBPAIR 1
#define NBSELF 2

#define NBTYPE NBPAIR
  #include "ComputeNonbondedMICKernelBase.h"
  #define CALCENERGY
    #include "ComputeNonbondedMICKernelBase.h"
  #undef CALCENERGY
  #define FULLELECT
    #include "ComputeNonbondedMICKernelBase.h"
    #define CALCENERGY
      #include "ComputeNonbondedMICKernelBase.h"
    #undef CALCENERGY
  #undef FULLELECT
#undef NBTYPE

#define NBTYPE NBSELF
  #include "ComputeNonbondedMICKernelBase.h"
  #define CALCENERGY
    #include "ComputeNonbondedMICKernelBase.h"
  #undef CALCENERGY
  #define FULLELECT
    #include "ComputeNonbondedMICKernelBase.h"
    #define CALCENERGY
      #include "ComputeNonbondedMICKernelBase.h"
    #undef CALCENERGY
  #undef FULLELECT
#undef NBTYPE


// This function is the main "computation" function, called each timestep to
//   initiate computation on the MIC device.  A signal is setup as part of
//   the offload pragma, which is periodically checked by the
//   mic_check_remote_kernel_complete and mic_check_local_kernel_complete
//   functions (which are called periodically by the Charm++ runtime system),
//   so progress can continue once the device has finished its computation.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - isRemote : flag indicating a 'remote' (1) or 'local' (0) kernel invocation (currently, must be 'local')
//   - lata : lattice information ('a' along with latb and latc)
//   - latb : lattice information ('b' along with lata and latc)
//   - latc : lattice information ('c' along with lata and latb)
//   - doSlow : flag indicating if slow forces should be calculated in this kernel invocation
//   - doEnergy : flag indiciating if energy information should be calculated in this kernel invocation
//   - usePairlists : flag indicating if pairlists should be used within this kernel invocation (currently, must be)
//   - savePairlists : flag indicating if pairlists should be saved across kernel invocations (currently, must be)
//   - atomsChanged : flag indicating whether or not the atoms (and related) buffers have changed (length, not content)
// Output: N/A
void mic_nonbonded_forces(const int deviceNum,
                          const int isRemote,
                          const int numLocalAtoms,
                          const int numLocalComputes,
                          const int numLocalPatches,
                          const mic_position3_t lata,
                          const mic_position3_t latb,
                          const mic_position3_t latc,
                          const int doSlow,
                          const int doEnergy,
                          const int usePairlists,
                          const int savePairlists,
                          const int atomsChanged
                         ) {

  const int numAtoms = host__atoms_size;
  const int numPatchPairs = host__patch_pairs_size;
  const int numForceLists = host__force_lists_size;

  // Get and set the tag to use
  int *tag_kernel = ((isRemote) ? (&tag_remote_kernel) : (&tag_local_kernel));
  *tag_kernel = 1;

  // Setup the kernel data structures
  mic_kernel_data * remote_kernel_data = host__kernel_data + 1; //host__remote_kernel_data;
  mic_kernel_data * local_kernel_data = host__kernel_data; //host__local_kernel_data;
  mic_kernel_data * kernel_data = ((isRemote) ? (remote_kernel_data) : (local_kernel_data));
  mic_kernel_data * _kernel_data = host__kernel_data;
  #if MIC_SYNC_INPUT != 0
    // NOTE: With MIC_SYNC_INPUT, all input is pushed before either of the "local" or "remote" kernels
    //   are executed, so setup both kernel_data structs as the "remote" kernel is being issued
    if (singleKernelFlag != 0 || isRemote) {
      remote_kernel_data->isRemote = 1;
      local_kernel_data->isRemote = 0;
      remote_kernel_data->numLocalAtoms = local_kernel_data->numLocalAtoms = numLocalAtoms;
      remote_kernel_data->numLocalComputes = local_kernel_data->numLocalComputes = numLocalComputes;
      remote_kernel_data->numLocalPatches = local_kernel_data->numLocalPatches = numLocalPatches;
      remote_kernel_data->doSlow = local_kernel_data->doSlow = doSlow;
      remote_kernel_data->doEnergy = local_kernel_data->doEnergy = doEnergy;
      remote_kernel_data->usePairlists = local_kernel_data->usePairlists = usePairlists;
      remote_kernel_data->savePairlists = local_kernel_data->savePairlists = savePairlists;
      remote_kernel_data->lata.x = local_kernel_data->lata.x = lata.x;
      remote_kernel_data->lata.y = local_kernel_data->lata.y = lata.y;
      remote_kernel_data->lata.z = local_kernel_data->lata.z = lata.z;
      remote_kernel_data->latb.x = local_kernel_data->latb.x = latb.x;
      remote_kernel_data->latb.y = local_kernel_data->latb.y = latb.y;
      remote_kernel_data->latb.z = local_kernel_data->latb.z = latb.z;
      remote_kernel_data->latc.x = local_kernel_data->latc.x = latc.x;
      remote_kernel_data->latc.y = local_kernel_data->latc.y = latc.y;
      remote_kernel_data->latc.z = local_kernel_data->latc.z = latc.z;
      remote_kernel_data->numAtoms = local_kernel_data->numAtoms = numAtoms;
      remote_kernel_data->numPatchPairs = local_kernel_data->numPatchPairs = numPatchPairs;
      remote_kernel_data->numForceLists = local_kernel_data->numForceLists = numForceLists;
      remote_kernel_data->forceBuffersReqSize = local_kernel_data->forceBuffersReqSize = host__force_buffers_req_size;
    }
  #else
    kernel_data->isRemote = isRemote;
    kernel_data->numLocalAtoms = numLocalAtoms;
    kernel_data->numLocalComputes = numLocalComputes;
    kernel_data->numLocalPatches = numLocalPatches;
    kernel_data->doSlow = doSlow;
    kernel_data->doEnergy = doEnergy;
    kernel_data->usePairlists = usePairlists;
    kernel_data->savePairlists = savePairlists;
    kernel_data->lata.x = lata.x;
    kernel_data->lata.y = lata.y;
    kernel_data->lata.z = lata.z;
    kernel_data->latb.x = latb.x;
    kernel_data->latb.y = latb.y;
    kernel_data->latb.z = latb.z;
    kernel_data->latc.x = latc.x;
    kernel_data->latc.y = latc.y;
    kernel_data->latc.z = latc.z;
    kernel_data->numAtoms = numAtoms;
    kernel_data->numPatchPairs = numPatchPairs;
    kernel_data->numForceLists = numForceLists;
    kernel_data->forceBuffersReqSize = host__force_buffers_req_size;
  #endif

  // For buffers that are only periodically copied/updated, get the
  //   flags that indicate whether or not those buffers should be
  //   transfered during this kernel invocation, and then reset them
  const int slowForcesNumAtoms = ((doSlow) ? (numAtoms) : (0));
  const int _patch_pairs_copySize = patch_pairs_copySize;  // from __thread variable to stack variable for use within offload pragma
  const int _force_lists_copySize = force_lists_copySize;  // from __thread variable to stack variable for use within offload pragma
  const int _atom_params_copySize = atom_params_copySize;  // from __thread variable to stack variable for use within offload pragma
  patch_pairs_copySize = 0;
  force_lists_copySize = 0;
  atom_params_copySize = 0;

  // Calculate the start/size of any data buffers that need to be copied to the device this timestep
  // NOTE: Some "remote" computes still need "local patch" data, so copy all patch data before/during "remote."
  int toCopySize_atoms = (singleKernelFlag != 0 || isRemote) ? (numAtoms) : (0);  // The remote kernel will copy in all atoms
  int toCopySize_atom_params = (singleKernelFlag != 0 || isRemote) ? (_atom_params_copySize) : (0);  // If there are atom_params to copy, the remote kernel will copy in all atom_params
  int toCopySize_patch_pairs = (singleKernelFlag != 0 || isRemote) ? (_patch_pairs_copySize) : (0);  // This could potentially be split between kernels (but is relatively rare, so don't worry about it for now)
  int toCopySize_force_lists = (singleKernelFlag != 0 || isRemote) ? (_force_lists_copySize) : (0);  // This could potentially be split between kernels (but is relatively rare, so don't worry about it for now)
  int toCopyStart_forces = (singleKernelFlag != 0) ? (0) : ((isRemote) ? (numLocalAtoms) : (0));
  int toCopySize_forces = (singleKernelFlag != 0) ? (numAtoms) : ((isRemote) ? (numAtoms - numLocalAtoms) : (numLocalAtoms));
  int toCopyStart_slow_forces = (doSlow) ? (toCopyStart_forces) : (0);
  int toCopySize_slow_forces = (doSlow) ? (toCopySize_forces) : (0);

  // If tracing data is to be recorded, setup local variables to the host buffers that will
  //   receive the data and setup the clauses that will be used to transfer the data
  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    double * device_times_start = host__device_times_start;
    int toCopySize_device_times_start = ((isRemote) ? (0) : (10));
    #if MIC_DEVICE_TRACING_DETAILED != 0
      double * device_times_computes = host__device_times_computes;
      double * device_times_patches = host__device_times_patches;
      int toCopySize_device_times_computes = ((isRemote) ? (0) : (2 * numPatchPairs));
      int toCopySize_device_times_patches = ((isRemote) ? (0) : (2 * numForceLists));
      #if MIC_SYNC_OUTPUT != 0
        #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED \
          in(device_times_computes[0:0] : alloc_if(0) free_if(0)) \
          in(device_times_patches[0:0] : alloc_if(0) free_if(0)) \
          nocopy(device__device_times_computes) nocopy(device__device_times_patches)
      #else
        #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED \
          out(device_times_computes[0:toCopySize_device_times_computes] : alloc_if(0) free_if(0)) \
          out(device_times_patches[0:toCopySize_device_times_patches] : alloc_if(0) free_if(0)) \
          nocopy(device__device_times_computes) nocopy(device__device_times_patches)
      #endif
    #else
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED
    #endif
    #if MIC_SYNC_OUTPUT != 0
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
        MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED \
        in(device_times_start[0:0] : alloc_if(0) free_if(0)) \
        nocopy(device__device_times_start)
    #else
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
        MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED \
        out(device_times_start[0:toCopySize_device_times_start] : alloc_if(0) free_if(0)) \
        nocopy(device__device_times_start)
    #endif
  #else
    #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE
  #endif

  // Setup local pointers and offload clauses for transfering atom input data
  atom * atoms = host__atoms;
  atom_param * atom_params = host__atom_params;
  double4 * forces = host__forces;
  double4 * slow_forces = host__slow_forces;
  #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
    #define MIC_DEVICE_ATOMS_CLAUSE
  #else
    #define MIC_DEVICE_ATOMS_CLAUSE \
      in(atoms[0:toCopySize_atoms] : alloc_if(0) free_if(0)) \
      in(atom_params[0:toCopySize_atom_params] : alloc_if(0) free_if(0))
  #endif

  // If kernel data transfer stats are being collected, calculate and print the amount of data
  //   being transfered, separated by type : in, out, or inout.
  // NOTE: Only include data required for production runs (e.g. not tracing data).
  #if MIC_KERNEL_DATA_TRANSFER_STATS != 0
    int transfer_in = sizeof(isRemote)
      + (3 * sizeof(lata)) // for lata, latb, and latc
      + sizeof(size_t) // for device__force_buffers_req_size
      + (sizeof(atom) * toCopySize_atoms)
      + (sizeof(atom_param) * toCopySize_atom_params)
      + (sizeof(patch_pair) * toCopySize_patch_pairs)
      + (sizeof(force_list) * toCopySize_force_lists);
    int transfer_inout = sizeof(mic_kernel_data);
    int transfer_out = sizeof(double4) * (toCopySize_forces + toCopySize_slow_forces)
      + sizeof(int); // for device__timestep
    #if MIC_DEBUG > 0
      MICP("timestep:%06d - transfering %d bytes (isRemote:%d, in:%d, inout:%d, out:%d)\n", \
           host__timestep, transfer_in + transfer_inout + transfer_out, \
           isRemote, transfer_in, transfer_inout, transfer_out \
          ); MICPF;
    #else
      printf("[MIC-DEBUG] :: PE %04d-%05d-%1d :: transfering %d bytes (in:%d, inout:%d, out:%d)\n",
             host__pe, host__timestep, isRemote,
             transfer_in + transfer_inout + transfer_out, transfer_in, transfer_inout, transfer_out
            );
    #endif
  #endif

  // Setup input clauses required for this kernel.
  // NOTE: If MIC_SYNC_INPUT is being used, use an offload pragma to first transfer the data to
  //   the device and use nocopy clauses for the actual kernel offload pragmas.
  patch_pair * patch_pairs = host__patch_pairs;
  force_list * force_lists = host__force_lists;
  #if MIC_SYNC_INPUT != 0

    // Push the data to the device
    if (singleKernelFlag != 0 || isRemote != 0) {
      MICP("<< pushing input data >>\n"); MICPF;
      #if MIC_TRACING != 0
        double input_start = CmiWallTimer();
      #endif
      #pragma offload_transfer target(mic:deviceNum) \
        in(patch_pairs[0:toCopySize_patch_pairs] : alloc_if(0) free_if(0)) \
        in(force_lists[0:toCopySize_force_lists] : alloc_if(0) free_if(0)) \
        in(_kernel_data[0:2] : alloc_if(0) free_if(0)) \
        MIC_DEVICE_ATOMS_CLAUSE
      { }
      #if MIC_TRACING != 0
        double input_end = CmiWallTimer();
        MIC_TRACING_RECORD(MIC_EVENT_SYNC_INPUT_PRAGMA, input_start, input_end);
      #endif
    }

    // Setup the clauses for the actual kernel offload pragma (so that data is not transfered)
    #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
      #define TMP_SYNC_ATOMS_CLAUSE
    #else
      #define TMP_SYNC_ATOMS_CLAUSE \
        in(atoms[0:0] : alloc_if(0) free_if(0)) \
        in(atom_params[0:0] : alloc_if(0) free_if(0))
    #endif

    #if MIC_SYNC_OUTPUT != 0
      #define TMP_SYNC_KERNEL_DATA_CLAUSE \
        in(_kernel_data[0:0] : alloc_if(0) free_if(0))
    #else
      #define TMP_SYNC_KERNEL_DATA_CLAUSE \
        out(_kernel_data[isRemote:1] : alloc_if(0) free_if(0))
    #endif

    #define MIC_DEVICE_SYNC_INPUT_CLAUSE \
      nocopy(device__lata) nocopy(device__latb) nocopy(device__latc) \
      nocopy(device__atoms_size) \
      in(patch_pairs[0:0] : alloc_if(0) free_if(0)) nocopy(device__patch_pairs_size) \
      in(force_lists[0:0] : alloc_if(0) free_if(0)) nocopy(device__force_lists_size) \
      TMP_SYNC_KERNEL_DATA_CLAUSE \
      TMP_SYNC_ATOMS_CLAUSE

  #else

    #if MIC_SYNC_OUTPUT != 0
      #define TMP_SYNC_KERNEL_DATA_CLAUSE \
        in(_kernel_data[isRemote:1] : alloc_if(0) free_if(0))
    #else
      #define TMP_SYNC_KERNEL_DATA_CLAUSE \
        inout(_kernel_data[isRemote:1] : alloc_if(0) free_if(0))
    #endif

    // Setup the clauses for the actual kernel offload pragma (so that data is transfered)
    #define MIC_DEVICE_SYNC_INPUT_CLAUSE \
      nocopy(device__lata) nocopy(device__latb) nocopy(device__latc) \
      in(patch_pairs[0:toCopySize_patch_pairs] : alloc_if(0) free_if(0)) \
      in(force_lists[0:toCopySize_force_lists] : alloc_if(0) free_if(0)) \
      nocopy(device__atoms_size) nocopy(device__patch_pairs_size) nocopy(device__force_lists_size) \
      TMP_SYNC_KERNEL_DATA_CLAUSE \
      MIC_DEVICE_ATOMS_CLAUSE

  #endif

  MICP("<< issuing %s kernel >>\n", (isRemote) ? ("remote") : ("local")); MICPF;

  #if MIC_SYNC_OUTPUT != 0
    #define MIC_DEVICE_SYNC_OUTPUT_CLAUSE \
      in(forces[0:0] : alloc_if(0) free_if(0)) \
      in(slow_forces[0:0] : alloc_if(0) free_if(0)) \
      nocopy(device__timestep)
  #else
    #define MIC_DEVICE_SYNC_OUTPUT_CLAUSE \
      out(forces[toCopyStart_forces:toCopySize_forces] : alloc_if(0) free_if(0)) \
      out(slow_forces[toCopyStart_slow_forces:toCopySize_slow_forces] : alloc_if(0) free_if(0)) \
      out(device__timestep : into(host__timestep))
  #endif

  #if MIC_DATA_STRUCT_VERIFY != 0
    #define MIC_DEVICE_DATA_STRUCT_VERIFY_CLAUSE \
      nocopy(device__exclusion_bits_copy) nocopy(device__constants_copy) \
      nocopy(device__table_four_copy) nocopy(device__table_four_float_copy) \
      nocopy(device__lj_table_copy) nocopy(device__lj_table_float_copy) \
      nocopy(device__patch_pairs_copy) nocopy(device__force_lists_copy) \
      nocopy(device__atoms_copy) nocopy(device__atom_params_copy) \
      nocopy(device__patch_pairs_copy_size) nocopy(device__force_lists_copy_size) \
      nocopy(device__atoms_copy_size) nocopy(device__atom_params_copy_size)
  #else
    #define MIC_DEVICE_DATA_STRUCT_VERIFY_CLAUSE
  #endif

  #if REFINE_PAIRLISTS != 0
    #define MIC_DEVICE_REFINE_PAIRLISTS_CLAUSE \
      nocopy(device__r2_array) \
      nocopy(device__pl_array) \
      nocopy(device__pl_size)
  #else
    #define MIC_DEVICE_REFINE_PAIRLISTS_CLAUSE
  #endif

  #if MIC_TRACING != 0
    double pragma_start = CmiWallTimer();
  #endif

  // Trigger computation on the MIC device (asynchronous offload pragma that does the kernel work)
  #pragma offload target(mic:deviceNum) \
    in(isRemote) \
    MIC_DEVICE_SYNC_INPUT_CLAUSE \
    MIC_DEVICE_SYNC_OUTPUT_CLAUSE \
    MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
    nocopy(device__force_buffers) nocopy(device__slow_force_buffers) \
    nocopy(device__force_buffers_alloc_size) nocopy(device__force_buffers_req_size) \
    nocopy(device__pairlists) nocopy(device__pairlists_alloc_size) \
    nocopy(device__table_four) nocopy(device__table_four_float) nocopy(device__table_four_n_16) \
    nocopy(device__lj_table) nocopy(device__lj_table_float) nocopy(device__lj_table_dim) \
    nocopy(device__exclusion_bits) \
    nocopy(device__constants) \
    nocopy(device__atoms) nocopy(device__atom_params) \
    nocopy(device__forces) nocopy(device__slow_forces) \
    nocopy(device__patch_pairs) nocopy(device__force_lists) \
    nocopy(device__numOMPThreads) \
    MIC_DEVICE_REFINE_PAIRLISTS_CLAUSE \
    MIC_DEVICE_DATA_STRUCT_VERIFY_CLAUSE \
    nocopy(device__pe) nocopy(device__node) \
    DEVICE_FPRINTF_CLAUSE
  {
//    DEVICE_FPRINTF_CLAUSE \
//    signal(tag_kernel)
    __MUST_BE_MIC;
    __FULL_CHECK(__ASSERT(isRemote == 0 || isRemote == 1));
    mic_kernel_data * kernel_data = _kernel_data + isRemote;
    __FULL_CHECK(__ASSERT(kernel_data != NULL));
    __FULL_CHECK(__ASSERT(device__numOMPThreads > 0));

    #if (MIC_DEVICE_FPRINTF != 0) && (MIC_DEVICE_FPRINTF_REOPEN_FREQ > 0)
      if ((device__timestep != 0 && device__timestep % MIC_DEVICE_FPRINTF_REOPEN_FREQ == 0) && (kernel_data->isRemote == 0)) {
        char filename[128] = { 0 };
        sprintf(filename, "/tmp/namd_deviceDebugInfo.%d", device__pe);
        if (device__node <= 0) {
          printf("[MIC-DEBUG] :: Reopening debug output files on MIC devices (timestep: %d).\n", device__timestep);
          fflush(NULL);
        }
        fclose(device__fout);
        device__fout = fopen(filename, "w");
        DEVICE_FPRINTF("Device file on PE %d (node: %d) - reopen (timestep: %d)...\n", device__pe, device__node, device__timestep);
      }
    #endif

    DEVICE_FPRINTF("%d %s ", device__timestep, (kernel_data->isRemote != 0 ? "R" : "L"));

    // DMK - TODO | FIXME - For now, copy values over into device__xxx variables
    //   so the older version of the code will still work (remove them eventually)
    device__patch_pairs = patch_pairs;
    device__force_lists = force_lists;
    device__atoms = atoms;
    device__atom_params = atom_params;
    device__forces = forces;
    device__slow_forces = slow_forces;
    __ASSERT(device__patch_pairs != NULL);
    __ASSERT(device__force_lists != NULL);
    __ASSERT(device__atoms != NULL);
    __ASSERT(device__atom_params != NULL);
    __ASSERT(device__forces != NULL);
    __ASSERT(device__slow_forces != NULL);
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
      device__device_times_start = device_times_start;
      __ASSERT(device__device_times_start);
      #if MIC_DEVICE_TRACING_DETAILED != 0
        device__device_times_computes = device_times_computes;
        device__device_times_patches = device_times_patches;
        __ASSERT(device__device_times_computes);
        __ASSERT(device__device_times_patches);
      #endif
    #endif
    device__lata.x = kernel_data->lata.x;
    device__lata.y = kernel_data->lata.y;
    device__lata.z = kernel_data->lata.z;
    device__latb.x = kernel_data->latb.x;
    device__latb.y = kernel_data->latb.y;
    device__latb.z = kernel_data->latb.z;
    device__latc.x = kernel_data->latc.x;
    device__latc.y = kernel_data->latc.y;
    device__latc.z = kernel_data->latc.z;
    device__atoms_size = kernel_data->numAtoms;
    device__patch_pairs_size = kernel_data->numPatchPairs;
    device__force_lists_size = kernel_data->numForceLists;
    device__force_buffers_req_size = kernel_data->forceBuffersReqSize;

    // TRACING - Record the start time of the kernel
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      __ASSERT(device__device_times_start != NULL);
      device__device_times_start[((isRemote)?(0):(1))] = getCurrentTime();
    #endif

    // If data structure verification is being performed, create copies of the buffers that can
    //   change throughout the simulation at the start of the timestep on the device, and verify
    //   the lookup tables before any work is performed.
    #if MIC_DATA_STRUCT_VERIFY != 0

      _verify_tables(device__pe, device__timestep, isRemote, -1);

      _verify_remalloc<patch_pair>(device__patch_pairs_copy_size, device__patch_pairs_size, device__patch_pairs_copy);
      _verify_remalloc<force_list>(device__force_lists_copy_size, device__force_lists_size, device__force_lists_copy);
      _verify_remalloc<atom>(device__atoms_copy_size, device__atoms_size, device__atoms_copy);
      _verify_remalloc<atom_param>(device__atom_params_copy_size, device__atoms_size, device__atom_params_copy);
      __ASSERT(device__patch_pairs_copy != NULL && device__force_lists_copy != NULL);
      __ASSERT(device__atoms_copy != NULL && device__atom_params_copy != NULL);

      _verify_parallel_memcpy(device__patch_pairs_copy, device__patch_pairs, sizeof(patch_pair) * device__patch_pairs_size);
      _verify_parallel_memcpy(device__force_lists_copy, device__force_lists, sizeof(force_list) * device__force_lists_size);
      _verify_parallel_memcpy(device__atoms_copy, device__atoms, sizeof(atom) * device__atoms_size);
      _verify_parallel_memcpy(device__atom_params_copy, device__atom_params, sizeof(atom_param) * device__atoms_size);
    #endif

    // Make sure there is enough memory allocated for the force buffers (reallocate if not)
    if (device__force_buffers_req_size > device__force_buffers_alloc_size) {
      if (device__force_buffers != NULL) { _MM_FREE_WRAPPER(device__force_buffers); }
      if (device__slow_force_buffers != NULL) { _MM_FREE_WRAPPER(device__slow_force_buffers); }
      device__force_buffers_alloc_size = device__force_buffers_req_size;
      size_t reqSize = sizeof(double) * 4 * device__force_buffers_alloc_size;
      device__force_buffers = (double4*)(_MM_MALLOC_WRAPPER(reqSize, 64, "device__force_buffers"));
      device__slow_force_buffers = (double4*)(_MM_MALLOC_WRAPPER(reqSize, 64, "device__slow_force_buffers"));
      __ASSERT(device__force_buffers != NULL && device__slow_force_buffers != NULL);
    }

    // Declare and initialize the overall virial summation/accumulation variables that will
    //   eventually be passed back to the host
    double virial_xx = 0.0;
    double virial_xy = 0.0;
    double virial_xz = 0.0;
    double virial_yy = 0.0;
    double virial_yz = 0.0;
    double virial_zz = 0.0;
    double fullElectVirial_xx = 0.0;
    double fullElectVirial_xy = 0.0;
    double fullElectVirial_xz = 0.0;
    double fullElectVirial_yy = 0.0;
    double fullElectVirial_yz = 0.0;
    double fullElectVirial_zz = 0.0;
    double vdwEnergy = 0.0;
    double electEnergy = 0.0;
    double fullElectEnergy = 0.0;
    #if MIC_EXCL_CHECKSUM_FULL != 0
      int exclusionSum = 0;
    #endif

    // Calculate the lower and upper bounds for the patch pairs that will be computed within this kernel
    const int ppI_start = (isRemote) ? (kernel_data->numLocalComputes) : (0);
    const int ppI_end = (device__constants->singleKernelFlag != 0 || isRemote) ? (device__patch_pairs_size) : (kernel_data->numLocalComputes);
    __FULL_CHECK(__ASSERT(ppI_start >= 0 && ppI_start <= device__patch_pairs_size));
    __FULL_CHECK(__ASSERT(ppI_end >= 0 && ppI_end <= device__patch_pairs_size));
    __FULL_CHECK(__ASSERT(ppI_start <= ppI_end));

    // TRACING - Record the start time for the patch pairs (computes) parallel region
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(2):(3))] = getCurrentTime();
    #endif

    DEVICE_FPRINTF(" C(%d,%d)", ppI_start, ppI_end);

    // Process each patch pair (compute)
    #pragma novector
    #if MULTIPLE_THREADS != 0
      #if MIC_EXCL_CHECKSUM_FULL != 0
        #define EXCL_CHECKSUM_CLAUSE  reduction(+ : exclusionSum)
      #else
        #define EXCL_CHECKSUM_CLAUSE
      #endif
      #pragma omp parallel for schedule(dynamic, 1) \
        reduction(+ : vdwEnergy, electEnergy, fullElectEnergy) \
        reduction(+ : virial_xx, virial_xy, virial_xz, virial_yy, virial_yz, virial_zz) \
        reduction(+ : fullElectVirial_xx, fullElectVirial_xy, fullElectVirial_xz, fullElectVirial_yy, fullElectVirial_yz, fullElectVirial_zz ) \
        EXCL_CHECKSUM_CLAUSE
      #undef EXCL_CHECKSUM_CLAUSE
    #endif
    for (int ppI = ppI_start; ppI < ppI_end; ppI++) {

      DEVICE_FPRINTF(" C");

      // TRACING - Record the start of each patch pair (compute)
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        __FULL_CHECK(__ASSERT(device__device_times_computes != NULL));
        device__device_times_computes[ppI * 2] = getCurrentTime();
      #endif

      // Create and populate the mic_params data structure (passed into the kernel function)
      mic_params params;

      // Initialize the pointer to the patch_pairs structure for this compute
      params.pp = (patch_pair*)(device__patch_pairs + ppI);
      __FULL_CHECK(__ASSERT(params.pp != NULL));

      #if MIC_SORT_LISTS != 0
        const int abSwapFlag = ((params.pp->patch1_size < params.pp->patch2_size) ? (1) : (0));
        #define ABSWAP(t, f)  ((abSwapFlag)?(t):(f))
      #else
        #define ABSWAP(t, f)  (f)
      #endif

      // DMK - DEBUG - Set a few values used for debugging
      params.ppI = ppI;
      params.p1 = ABSWAP(params.pp->p2, params.pp->p1);
      params.p2 = ABSWAP(params.pp->p1, params.pp->p2);
      params.pe = device__pe;
      params.timestep = device__timestep;

      // If pairlist refinement is enabled, initialize the related scratch buffer pointers
      // In some cases, buffers used by the compute kernel functions (calc_self, etc.) only
      //   require scratch buffers (i.e. not used for input or output).  As such, we can
      //   allocate these buffers on a one-per-thread basis rather than a one-per-compute
      //   basis to save some memory (and get memory location reuse for any given thread).
      //   These fields relate to the pointer and sizes for those buffers.
      #if REFINE_PAIRLISTS != 0
        const int threadIndex = omp_get_thread_num();
        params.plArray_ptr = device__pl_array + threadIndex;
        params.plSize_ptr = device__pl_size + threadIndex;
        params.r2Array_ptr = device__r2_array + threadIndex;
      #endif

      // Initialize pointers to the various lookup tables, constants, etc. that are used
      //   from within the compute functions
      #if MIC_HANDCODE_FORCE_SINGLE != 0
        params.table_four_base_ptr = (void*)device__table_four_float;
        params.lj_table_base_ptr = (void*)device__lj_table_float;
      #else
        params.table_four_base_ptr = (void*)device__table_four;
        params.lj_table_base_ptr = (void*)device__lj_table;
      #endif
      params.table_four_n_16 = device__table_four_n_16;
      params.lj_table_dim = device__lj_table_dim;
      params.exclusion_bits = device__exclusion_bits;
      params.constants = device__constants;
      __FULL_CHECK(__ASSERT(params.table_four_base_ptr != NULL));
      __FULL_CHECK(__ASSERT(params.lj_table_base_ptr != NULL));
      __FULL_CHECK(__ASSERT(params.table_four_n_16 % 16 == 0));
      __FULL_CHECK(__ASSERT(params.exclusion_bits != NULL));
      __FULL_CHECK(__ASSERT(params.constants != NULL));

      // Setup pairlist flags and pairlist buffer (an array of pointers, one for each pairlist,
      //   for each compute)
      params.usePairlists = kernel_data->usePairlists;
      params.savePairlists = kernel_data->savePairlists;
      params.pairlists_ptr = &(((int**)device__pairlists)[NUM_PAIRLIST_TYPES * ppI]);
      __FULL_CHECK(__ASSERT(device__pairlists != NULL));
      __FULL_CHECK(__ASSERT(params.pairlists_ptr != NULL));

      // Setup the sizes of the atom lists and force lists
      int n0 = ABSWAP(params.pp->patch2_size, params.pp->patch1_size);
      int n1 = ABSWAP(params.pp->patch1_size, params.pp->patch2_size);
      int n0_16 = (n0 + 15) & (~15); // Round up to nearest multiple of 16
      int n1_16 = (n1 + 15) & (~15);
      params.numAtoms[0] = n0;
      params.numAtoms[1] = n1;
      params.numAtoms_16[0] = n0_16;
      params.numAtoms_16[1] = n1_16;
      __FULL_CHECK(__ASSERT(params.numAtoms[0] >= 0));
      __FULL_CHECK(__ASSERT(params.numAtoms[1] >= 0));

      // Setup the pointers to the input particle data for this compute
      // NOTE: All of the particle data is sent in two chunks (atoms and atom_params, where
      //   atoms changes every timestep but atom_params only changes periodically).  Each
      //   of these buffers contains several sub-arrays (one for each field) in an
      //   structure-of-arrays (SOA) format.  This code, along with the force code
      //   below it, is simply creating the array pointers for each "field's array."
      #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
        mic_position_t* pBase0 = (mic_position_t*)(params.pp->patch1_atomDataPtr);
        mic_position_t* pBase1 = (mic_position_t*)(params.pp->patch2_atomDataPtr);
        int *pExtBase0 = (int*)(pBase0 + n0_16);
        int *pExtBase1 = (int*)(pBase1 + n1_16);
      #else
        __FULL_CHECK(__ASSERT(device__atoms != NULL));
        __FULL_CHECK(__ASSERT(device__atom_params != NULL));

        __FULL_CHECK(__ASSERT(params.pp->patch1_atom_start >= 0));
        __FULL_CHECK(__ASSERT((params.pp->patch1_atom_start < device__atoms_size) || (params.pp->patch1_atom_start == device__atoms_size && n0 == 0)));
        __FULL_CHECK(__ASSERT(params.pp->patch2_atom_start >= 0));
        __FULL_CHECK(__ASSERT((params.pp->patch2_atom_start < device__atoms_size) || (params.pp->patch2_atom_start == device__atoms_size && n1 == 0)));

        mic_position_t* pBase0 = (mic_position_t*)(device__atoms + ABSWAP(params.pp->patch2_atom_start, params.pp->patch1_atom_start));
        mic_position_t* pBase1 = (mic_position_t*)(device__atoms + ABSWAP(params.pp->patch1_atom_start, params.pp->patch2_atom_start));
        int* pExtBase0 = (int*)(device__atom_params + ABSWAP(params.pp->patch2_atom_start, params.pp->patch1_atom_start));
        int* pExtBase1 = (int*)(device__atom_params + ABSWAP(params.pp->patch1_atom_start, params.pp->patch2_atom_start));
      #endif
      __FULL_CHECK(__ASSERT(pBase0 != NULL));
      __FULL_CHECK(__ASSERT(pBase1 != NULL));
      __FULL_CHECK(__ASSERT(pExtBase0 != NULL));
      __FULL_CHECK(__ASSERT(pExtBase1 != NULL));
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        params.p[0] = (atom*)pBase0;
        params.p[1] = (atom*)pBase1;
        params.pExt[0] = (atom_param*)pExtBase0;
        params.pExt[1] = (atom_param*)pExtBase1;
      #else
        params.p_x[0] = pBase0 + (0 * n0_16); params.p_x[1] = pBase1 + (0 * n1_16);
        params.p_y[0] = pBase0 + (1 * n0_16); params.p_y[1] = pBase1 + (1 * n1_16);
        params.p_z[0] = pBase0 + (2 * n0_16); params.p_z[1] = pBase1 + (2 * n1_16);
        params.p_q[0] = pBase0 + (3 * n0_16); params.p_q[1] = pBase1 + (3 * n1_16);
        params.pExt_vdwType[0] = pExtBase0 + (0 * n0_16);
        params.pExt_index[0] = pExtBase0 + (1 * n0_16);
        params.pExt_exclIndex[0] = pExtBase0 + (2 * n0_16);
        params.pExt_exclMaxDiff[0] = pExtBase0 + (3 * n0_16);
        params.pExt_vdwType[1] = pExtBase1 + (0 * n1_16);
        params.pExt_index[1] = pExtBase1 + (1 * n1_16);
        params.pExt_exclIndex[1] = pExtBase1 + (2 * n1_16);
        params.pExt_exclMaxDiff[1] = pExtBase1 + (3 * n1_16);
      #endif

      // Setup the pointers to the output force data for this compute
      // NOTE: Forces are output every timestep, but slow forces (fullf) are output only
      //   during some timesteps.
      __FULL_CHECK(__ASSERT(device__force_lists != NULL));
      __FULL_CHECK(__ASSERT(device__force_buffers != NULL));
      __FULL_CHECK(__ASSERT(device__slow_force_buffers != NULL));
      __FULL_CHECK(__ASSERT(params.pp->patch1_force_list_index < device__force_lists_size));
      __FULL_CHECK(__ASSERT(params.pp->patch2_force_list_index < device__force_lists_size));
      __FULL_CHECK(__ASSERT(params.pp->patch1_force_start < device__force_buffers_req_size));
      __FULL_CHECK(__ASSERT(params.pp->patch2_force_start < device__force_buffers_req_size));
      params.doSlow = kernel_data->doSlow;  // Flag indicating if slow forces should be calculate this timestep or not
      params.fl[0] = (force_list*)(device__force_lists + ABSWAP(params.pp->patch2_force_list_index, params.pp->patch1_force_list_index));
      params.fl[1] = (force_list*)(device__force_lists + ABSWAP(params.pp->patch1_force_list_index, params.pp->patch2_force_list_index));
      double *ffBase0 = (double*)(device__force_buffers + ABSWAP(params.pp->patch2_force_start, params.pp->patch1_force_start));
      double *ffBase1 = (double*)(device__force_buffers + ABSWAP(params.pp->patch1_force_start, params.pp->patch2_force_start));
      double *fullfBase0 = (double*)(device__slow_force_buffers + ABSWAP(params.pp->patch2_force_start, params.pp->patch1_force_start));
      double *fullfBase1 = (double*)(device__slow_force_buffers + ABSWAP(params.pp->patch1_force_start, params.pp->patch2_force_start));
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        params.ff[0] = (double4*)ffBase0;
        params.ff[1] = (double4*)ffBase1;
        params.fullf[0] = (double4*)fullfBase0;
        params.fullf[1] = (double4*)fullfBase1;
        __ASSERT(params.ff[0] != NULL);
        __ASSERT(params.ff[1] != NULL);
        __ASSERT(params.fullf[0] != NULL);
        __ASSERT(params.fullf[1] != NULL);
      #else
        params.ff_x[0] = ffBase0 + (0 * n0_16); params.ff_x[1] = ffBase1 + (0 * n1_16);
        params.ff_y[0] = ffBase0 + (1 * n0_16); params.ff_y[1] = ffBase1 + (1 * n1_16);
        params.ff_z[0] = ffBase0 + (2 * n0_16); params.ff_z[1] = ffBase1 + (2 * n1_16);
        params.ff_w[0] = ffBase0 + (3 * n0_16); params.ff_w[1] = ffBase1 + (3 * n1_16);
        params.fullf_x[0] = fullfBase0 + (0 * n0_16); params.fullf_x[1] = fullfBase1 + (0 * n1_16);
        params.fullf_y[0] = fullfBase0 + (1 * n0_16); params.fullf_y[1] = fullfBase1 + (1 * n1_16);
        params.fullf_z[0] = fullfBase0 + (2 * n0_16); params.fullf_z[1] = fullfBase1 + (2 * n1_16);
        params.fullf_w[0] = fullfBase0 + (3 * n0_16); params.fullf_w[1] = fullfBase1 + (3 * n1_16);
        __ASSERT(params.ff_x[0] != NULL); __ASSERT(params.ff_x[1] != NULL);
        __ASSERT(params.ff_y[0] != NULL); __ASSERT(params.ff_y[1] != NULL);
        __ASSERT(params.ff_z[0] != NULL); __ASSERT(params.ff_z[1] != NULL);
        __ASSERT(params.ff_w[0] != NULL); __ASSERT(params.ff_w[1] != NULL);
        __ASSERT(params.fullf_x[0] != NULL); __ASSERT(params.fullf_x[1] != NULL);
        __ASSERT(params.fullf_y[0] != NULL); __ASSERT(params.fullf_y[1] != NULL);
        __ASSERT(params.fullf_z[0] != NULL); __ASSERT(params.fullf_z[1] != NULL);
        __ASSERT(params.fullf_w[0] != NULL); __ASSERT(params.fullf_w[1] != NULL);
      #endif

      // Create the offsets for the first list of particles.
      // NOTE: In this version of the nonbonded code, the positions of the atoms are stored
      //   as offsets within the patch in which they are located.  These offsets represent
      //   the offsets between the two patches being interacted (including periodic boundaries).
      params.offset.x = params.pp->offset.x * device__lata.x
                      + params.pp->offset.y * device__latb.x
                      + params.pp->offset.z * device__latc.x;
      params.offset.y = params.pp->offset.x * device__lata.y
                      + params.pp->offset.y * device__latb.y
                      + params.pp->offset.z * device__latc.y;
      params.offset.z = params.pp->offset.x * device__lata.z
                      + params.pp->offset.y * device__latb.z
                      + params.pp->offset.z * device__latc.z;
      #if MIC_SORT_LISTS != 0
        params.offset.x *= ABSWAP(-1.0f, 1.0f);
        params.offset.y *= ABSWAP(-1.0f, 1.0f);
        params.offset.z *= ABSWAP(-1.0f, 1.0f);
      #endif

      // DMK - DEBUG - Record the center point for each of the input patches
      float patch1_center_x = params.pp->patch1_center_x * device__lata.x
                             + params.pp->patch1_center_y * device__latb.x
                             + params.pp->patch1_center_z * device__latc.x;
      float patch1_center_y = params.pp->patch1_center_x * device__lata.y
                             + params.pp->patch1_center_y * device__latb.y
                             + params.pp->patch1_center_z * device__latc.y;
      float patch1_center_z = params.pp->patch1_center_x * device__lata.z
                             + params.pp->patch1_center_y * device__latb.z
                             + params.pp->patch1_center_z * device__latc.z;
      float patch2_center_x = params.pp->patch2_center_x * device__lata.x
                             + params.pp->patch2_center_y * device__latb.x
                             + params.pp->patch2_center_z * device__latc.x;
      float patch2_center_y = params.pp->patch2_center_x * device__lata.y
                             + params.pp->patch2_center_y * device__latb.y
                             + params.pp->patch2_center_z * device__latc.y;
      float patch2_center_z = params.pp->patch2_center_x * device__lata.z
                             + params.pp->patch2_center_y * device__latb.z
                             + params.pp->patch2_center_z * device__latc.z;
      params.patch1_center_x = ABSWAP(patch2_center_x, patch1_center_x);
      params.patch2_center_x = ABSWAP(patch1_center_x, patch2_center_x);
      params.patch1_center_y = ABSWAP(patch2_center_y, patch1_center_y);
      params.patch2_center_y = ABSWAP(patch1_center_y, patch2_center_y);
      params.patch1_center_z = ABSWAP(patch2_center_z, patch1_center_z);
      params.patch2_center_z = ABSWAP(patch1_center_z, patch2_center_z);

      // Initialize the virial accumulators for this compute (zero them out)
      params.virial_xx = 0;
      params.virial_xy = 0;
      params.virial_xz = 0;
      params.virial_yy = 0;
      params.virial_yz = 0;
      params.virial_zz = 0;
      params.fullElectVirial_xx = 0;
      params.fullElectVirial_xy = 0;
      params.fullElectVirial_xz = 0;
      params.fullElectVirial_yy = 0;
      params.fullElectVirial_yz = 0;
      params.fullElectVirial_zz = 0;
      params.vdwEnergy = 0.0;
      params.electEnergy = 0.0;
      params.fullElectEnergy = 0.0;
      #if MIC_EXCL_CHECKSUM_FULL != 0
        params.exclusionSum = 0;
      #endif

      // Select the version of the kernel to call based on the timestep's requirements
      //   and what type of compute this compute is
      int isSelf = (params.pp->patch1_force_start == params.pp->patch2_force_start);
                    // NOTE: Many ways to check this (arbitrary test used here)
      int selfBit = ((isSelf) ? (0x01) : (0x00));
      int doSlowBit = ((kernel_data->doSlow) ? (0x02) : (0x00));
      int doEnergyBit = ((kernel_data->doEnergy) ? (0x04) : (0x00));
      int kernelSelect = selfBit | doSlowBit | doEnergyBit;
      DEVICE_FPRINTF("%d", kernelSelect);
      switch (kernelSelect) {
        case 0x00: calc_pair(params); break;
        case 0x01: calc_self(params); break;
        case 0x02: calc_pair_fullelect(params); break;
        case 0x03: calc_self_fullelect(params); break;
        case 0x04: calc_pair_energy(params); break;
        case 0x05: calc_self_energy(params); break;
        case 0x06: calc_pair_energy_fullelect(params); break;
        case 0x07: calc_self_energy_fullelect(params); break;
        default:
          mic_dev_die("!!! INVALID KERNEL SELECTION ON MIC DEVICE !!!\n");
          break;
      } // end switch (kernelSelect)

      // Contribute this compute's virial summations into the overall virial summation
      //   that will be passed back to the host
      #if MIC_SORT_LISTS != 0
        if (abSwapFlag) {
          virial_xx -= params.virial_xx;
          virial_xy -= params.virial_xy;
          virial_xz -= params.virial_xz;
          virial_yy -= params.virial_yy;
          virial_yz -= params.virial_yz;
          virial_zz -= params.virial_zz;
          fullElectVirial_xx -= params.fullElectVirial_xx;
          fullElectVirial_xy -= params.fullElectVirial_xy;
          fullElectVirial_xz -= params.fullElectVirial_xz;
          fullElectVirial_yy -= params.fullElectVirial_yy;
          fullElectVirial_yz -= params.fullElectVirial_yz;
          fullElectVirial_zz -= params.fullElectVirial_zz;
          vdwEnergy -= params.vdwEnergy;
          electEnergy -= params.electEnergy;
          fullElectEnergy -= params.fullElectEnergy;
        } else {
      #endif
          virial_xx += params.virial_xx;
          virial_xy += params.virial_xy;
          virial_xz += params.virial_xz;
          virial_yy += params.virial_yy;
          virial_yz += params.virial_yz;
          virial_zz += params.virial_zz;
          fullElectVirial_xx += params.fullElectVirial_xx;
          fullElectVirial_xy += params.fullElectVirial_xy;
          fullElectVirial_xz += params.fullElectVirial_xz;
          fullElectVirial_yy += params.fullElectVirial_yy;
          fullElectVirial_yz += params.fullElectVirial_yz;
          fullElectVirial_zz += params.fullElectVirial_zz;
          vdwEnergy += params.vdwEnergy;
          electEnergy += params.electEnergy;
          fullElectEnergy += params.fullElectEnergy;
      #if MIC_SORT_LISTS != 0
	}
      #endif

      #if MIC_EXCL_CHECKSUM_FULL != 0
        exclusionSum += params.exclusionSum;
      #endif

      #undef ABSWAP

      // TRACING - Record the end time for the compute
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        device__device_times_computes[ppI * 2 + 1] = getCurrentTime();
      #endif

    } // end parallel for (ppI < device__patch_pairs_size)

    // TRACING - Record the end time for all the computes
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(4):(5))] = getCurrentTime();
    #endif

    // Store the reduced virial values into the kernel_data structure to be passed
    //   back to the host core
    kernel_data->virial_xx = virial_xx;
    kernel_data->virial_xy = virial_xy;
    kernel_data->virial_xz = virial_xz;
    kernel_data->virial_yy = virial_yy;
    kernel_data->virial_yz = virial_yz;
    kernel_data->virial_zz = virial_zz;
    kernel_data->fullElectVirial_xx = fullElectVirial_xx;
    kernel_data->fullElectVirial_xy = fullElectVirial_xy;
    kernel_data->fullElectVirial_xz = fullElectVirial_xz;
    kernel_data->fullElectVirial_yy = fullElectVirial_yy;
    kernel_data->fullElectVirial_yz = fullElectVirial_yz;
    kernel_data->fullElectVirial_zz = fullElectVirial_zz;
    kernel_data->vdwEnergy = vdwEnergy;
    kernel_data->electEnergy = electEnergy;
    kernel_data->fullElectEnergy = fullElectEnergy;
    #if MIC_EXCL_CHECKSUM_FULL != 0
      kernel_data->exclusionSum = exclusionSum;
    #endif

    // Calculate the upper and lower bounds for the force lists to be processed by this kernel
    int flI_start = (isRemote) ? (kernel_data->numLocalPatches) : (0);
    int flI_end = (device__constants->singleKernelFlag != 0 || isRemote) ? (device__force_lists_size) : (kernel_data->numLocalPatches);
    __FULL_CHECK(__ASSERT(flI_start >= 0 && flI_start <= device__force_lists_size));
    __FULL_CHECK(__ASSERT(flI_end >= 0 && flI_start <= device__force_lists_size));
    __FULL_CHECK(__ASSERT(flI_start <= flI_end));

    // Because there are fewer patches than computes, create more parallelism by sub-dividing the
    //   the force lists (patches), causing each patch to be processed by multiple threads
    const int numThreads = device__numOMPThreads; //omp_get_max_threads();
    const int numForceLoopIters = flI_end - flI_start;
    int numForceLoopSplits = 1;
    if ((numForceLoopIters > 0) && ((2 * numThreads) > numForceLoopIters)) {
      numForceLoopSplits = (2 * numThreads) / numForceLoopIters;  // NOTE: 2 * numThreads to break it up more (smaller chunks of work)
      if (numForceLoopSplits < 2) { numForceLoopSplits = 2; }  // At least split in half
      if (numForceLoopSplits > 16) { numForceLoopSplits = 16; }  // Don't split any single patch too fine (threading overhead costs)
    }
    DEVICE_FPRINTF(" F(%d,%d:%d)", flI_start, flI_end, numForceLoopSplits);
    flI_start *= numForceLoopSplits;
    flI_end *= numForceLoopSplits;

    // TRACING - Record the start time for the processing of all force lists (patches)
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(6):(7))] = getCurrentTime();
    #endif

    // Process each of the force lists (patches)
    #pragma novector
    #if MULTIPLE_THREADS != 0
      #pragma omp parallel for schedule(dynamic)
    #endif
    for (int _flI = flI_start; _flI < flI_end; _flI++) {

      // From _flI, calculate the flI (patch index) and flIPart (thread per patch) values
      int flI = _flI / numForceLoopSplits;
      __FULL_CHECK(__ASSERT(flI >= 0 && flI < device__force_lists_size));
      int flIPart = _flI % numForceLoopSplits;
      
      // TRACING - Record the start of each force reduction
      // DMK - FIXME | NOTE : THIS ONLY RECORDS ONE TASK/ITERATION OF EACH SPLIT SET, SO
      //   PROJECTIONS WILL ONLY SHOW SOME OF THE FORCE "TASKS" !!! (the problem is
      //   numForceLoopSplits can differ between the "remote" and "local" kernels).
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        // DMK - NOTE : Don't select the same part each time (e.g. only flIPart == 0)
        if (flIPart == flI % numForceLoopSplits) {
          device__device_times_patches[flI * 2] = getCurrentTime();
	}
      #endif

      // Get the output force list pointer and force array length
      const force_list &fl = device__force_lists[flI];
      const int f_len = fl.patch_stride * 4; // NOTE : number of individual doubles
      __ASSUME(f_len % 16 == 0);

      // Calculate this "task's" portion of the patch object's work (flIPartILo -> flIPartIHi)
      int flIPartJLo = (int)(((float)(f_len)) * (((float)(flIPart    )) / ((float)(numForceLoopSplits))));
      int flIPartJHi = (int)(((float)(f_len)) * (((float)(flIPart + 1)) / ((float)(numForceLoopSplits))));
      // NOTE: Force flIPartJLo and flIPartJHi to cacheline boundaries to avoid false sharing
      flIPartJLo = (flIPartJLo + 7) & (~7);
      flIPartJHi = (flIPartJHi + 7) & (~7);
      if (flIPartJHi > f_len) { flIPartJHi = f_len; }
      __ASSUME(flIPartJLo % 8 == 0);
      __ASSUME(flIPartJHi % 8 == 0);  // NOTE: true after capping at f_len since f_len % 16 == 0

      // Sum the force output for each compute contributing to this force list (patch)
      {
        // Setup the pointer to the final force array that will be passed back up to the host
        __FULL_CHECK(__ASSERT(device__forces != NULL));
        __FULL_CHECK(__ASSERT(fl.force_output_start >= 0));
        __FULL_CHECK(__ASSERT((fl.force_output_start < device__atoms_size) || (fl.force_output_start == device__atoms_size && f_len == 0)));
        double * RESTRICT fDst = (double*)(device__forces + fl.force_output_start);
        __FULL_CHECK(__ASSERT(fDst != NULL));
        __ASSUME_ALIGNED(fDst);

        // Setup the pointer to the list of arrays, each with output from one of
        //   the compute objects that contributed to this patch's force data
        __FULL_CHECK(__ASSERT(device__force_buffers != NULL));
        __FULL_CHECK(__ASSERT(fl.force_list_start >= 0));
        __FULL_CHECK(__ASSERT((fl.force_list_start < device__force_buffers_req_size) || (fl.force_list_start == device__force_buffers_req_size && f_len == 0)));
        double *fSrcBase = (double*)(device__force_buffers + fl.force_list_start);
        __FULL_CHECK(__ASSERT(fSrcBase != NULL));
        __ASSUME_ALIGNED(fSrcBase);

        // Accumulate the forces from the various computes contributing to this patch
        memset(fDst + flIPartJLo, 0, sizeof(double) * (flIPartJHi - flIPartJLo));
        for (int i = 0; i < fl.force_list_size; i++) {
          #pragma simd
          for (int j = flIPartJLo; j < flIPartJHi; j++) {
            fDst[j] += fSrcBase[i * f_len + j];
          }
        }
      }

      // Sum the output for each contributing compute
      if (kernel_data->doSlow) {

        // Setup the pointer to the final force array that will be passed back up to the host
        __FULL_CHECK(__ASSERT(device__slow_forces != NULL));
        double * RESTRICT fDst = (double*)(device__slow_forces + fl.force_output_start);
        __FULL_CHECK(__ASSERT(fDst != NULL));
        __ASSUME_ALIGNED(fDst);

        // Setup the pointer to the list of arrays, each with output from one of
        //   the compute objects that contributed to this patch's force data
        __FULL_CHECK(__ASSERT(device__slow_force_buffers != NULL));
        double *fSrcBase = (double*)(device__slow_force_buffers + fl.force_list_start);
        __FULL_CHECK(__ASSERT(fSrcBase != NULL));
        __ASSUME_ALIGNED(fSrcBase);

        // Accumulate the forces from the various computes contributing to this patch
        memset(fDst + flIPartJLo, 0, sizeof(double) * (flIPartJHi - flIPartJLo));
        for (int i = 0; i < fl.force_list_size; i++) {
          #pragma simd
          for (int j = flIPartJLo; j < flIPartJHi; j++) {
            fDst[j] += fSrcBase[i * f_len + j];
          }
        }
      }
     
      // TRACING - Record the end time for each force list (patch)
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        device__device_times_patches[flI * 2 + 1] = getCurrentTime();
      #endif

    } // end parallel for (flI < device__force_lists_size)

    // TRACING - Record the end time for all the force lists (patches)
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(8):(9))] = getCurrentTime();
    #endif

    // DMK - DEBUG
    #if MIC_DATA_STRUCT_VERIFY != 0
      _verify_data_structures(device__pe, device__timestep, isRemote, 2, ppI_start, ppI_end, kernel_data->doSlow);
      _verify_buffers_match(device__pe, device__timestep, isRemote, 2, (char*)device__patch_pairs, (char*)device__patch_pairs_copy, sizeof(patch_pairs) * device__patch_pairs_size, "device__patch_pairs_copy");
      _verify_buffers_match(device__pe, device__timestep, isRemote, 2, (char*)device__force_lists, (char*)device__force_lists_copy, sizeof(force_lists) * device__force_lists_size, "device__force_lists_copy");
      _verify_buffers_match(device__pe, device__timestep, isRemote, 2, (char*)device__atoms, (char*)device__atoms_copy, sizeof(atom) * device__atoms_size, "device__atoms_copy");
      _verify_buffers_match(device__pe, device__timestep, isRemote, 2, (char*)device__atom_params, (char*)device__atom_params_copy, sizeof(atom_param) * device__atoms_size, "device__atom_params_copy");
    #endif

    // DMK - DEBUG
    #if MIC_TRACK_DEVICE_MEM_USAGE != 0
      if (device__timestep % 100 == 0 && isRemote == 0) {

        uint64_t plTotalMem = sizeof(void*) * device__pairlists_alloc_size;
        uint64_t plUsedMem = plTotalMem;
        uint64_t plPeakUsedMem = plTotalMem;
        if (device__pairlists != NULL) {
          for (int k = 0; k < device__pairlists_alloc_size; k++) {
            int * plPtr = ((int**)device__pairlists)[k];
            if (plPtr != NULL) {
              plTotalMem += plPtr[0] * sizeof(int) + 64;
              plUsedMem += plPtr[1] * sizeof(int);
              plPeakUsedMem += plPtr[-1] * sizeof(int);
            }
	  }
	}

        uint64_t fbsTotalMem = device__force_buffers_alloc_size * 4 * sizeof(double);
        MemInfo mi;
        readMemInfo(&mi);

        printf("[MIC-MEMINFO] :: PE %d :: MIC mem usage (MB) -- ts: %d - "
               "f: %.3lf, t: %.3lf, c: %.3lf, a: %.3lf, i: %.3lf - "
               "vS: %.3lf, vP: %.3lf - "
               "plTM: %.3lf, plUM: %.3lf, plUPM: %.3lf, fbsTM(x2): %.3lf\n",
               ((int)(device__pe)),
               ((int)(device__timestep)),
               ((double)(mi.memFree)) * 1.0e-6,
               ((double)(mi.memTotal)) * 1.0e-6,
               ((double)(mi.cached)) * 1.0e-6,
               ((double)(mi.active)) * 1.0e-6,
               ((double)(mi.inactive)) * 1.0e-6,
               ((double)(mi.vmSize)) * 1.0e-6,
               ((double)(mi.vmPeak)) * 1.0e-6,
               ((double)(plTotalMem)) * 1.0e-6,
               ((double)(plUsedMem)) * 1.0e-6,
               ((double)(plPeakUsedMem)) * 1.0e-6,
               ((double)(fbsTotalMem)) * 1.0e-6
              ); fflush(NULL);
      }
    #endif

    // DMK - DEBUG - Increment the timestep counter on the device
    if (kernel_data->isRemote == 0) { device__timestep++; }

    DEVICE_FPRINTF(" |\n");

  } // end pragma offload

  #if MIC_TRACING != 0
    double pragma_end = CmiWallTimer();
    MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_PRAGMA, pragma_start, pragma_end);
  #endif

  #if MIC_SYNC_INPUT != 0
    #undef TMP_SYNC_ATOMS_CLAUSE
  #endif
  #undef TMP_SYNC_KERNEL_DATA_CLAUSE
  #undef MIC_DEVICE_SYNC_INPUT_CLAUSE

  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)  
    #undef MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED
  #endif
  #undef MIC_DEVICE_TIMING_PRAGMA_CLAUSE

  #undef MIC_DEVICE_DATA_STRUCT_VERIFY_CLAUSE
  #undef MIC_DEVICE_REFINE_PAIRLISTS_CLAUSE
}

void mic_transfer_output(const int deviceNum,
                         const int isRemote,
                         const int numLocalAtoms,
                         const int doSlow
                        ) {

  const int numAtoms = host__atoms_size;

  mic_kernel_data * kernel_data = host__kernel_data;

  int toCopyStart_forces = (singleKernelFlag != 0) ? (0) : ((isRemote) ? (numLocalAtoms) : (0));
  int toCopySize_forces = (singleKernelFlag != 0) ? (numAtoms) : ((isRemote) ? (numAtoms - numLocalAtoms) : (numLocalAtoms));
  int toCopyStart_slow_forces = (doSlow) ? (toCopyStart_forces) : (0);
  int toCopySize_slow_forces = (doSlow) ? (toCopySize_forces) : (0);

  double4 * forces = host__forces;
  double4 * slow_forces = host__slow_forces;

  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    const int numPatchPairs = host__patch_pairs_size;
    const int numForceLists = host__force_lists_size;
    double * device_times_start = host__device_times_start;
    int toCopySize_device_times_start = ((isRemote) ? (0) : (10));
    #if MIC_DEVICE_TRACING_DETAILED != 0
      double * device_times_computes = host__device_times_computes;
      double * device_times_patches = host__device_times_patches;
      int toCopySize_device_times_computes = ((isRemote) ? (0) : (2 * numPatchPairs));
      int toCopySize_device_times_patches = ((isRemote) ? (0) : (2 * numForceLists));
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED \
        out(device_times_computes[0:toCopySize_device_times_computes] : alloc_if(0) free_if(0)) \
        out(device_times_patches[0:toCopySize_device_times_patches] : alloc_if(0) free_if(0))
    #else
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED
    #endif
    #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
      MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED \
      out(device_times_start[0:toCopySize_device_times_start] : alloc_if(0) free_if(0))
  #else
    #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE
  #endif

  #pragma offload_transfer target(mic:deviceNum) \
    out(forces[toCopyStart_forces:toCopySize_forces] : alloc_if(0) free_if(0)) \
    out(slow_forces[toCopyStart_slow_forces:toCopySize_slow_forces] : alloc_if(0) free_if(0)) \
    out(kernel_data[isRemote:1] : alloc_if(0) free_if(0)) \
    MIC_DEVICE_TIMING_PRAGMA_CLAUSE
  { }
  //  out(device__timestep : into(host__timestep))

  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    #undef MIC_DEVICE_TIMING_PRAGMA_CLAUSE_DETAILED
  #endif
  #undef MIC_DEVICE_TIMING_PRAGMA_CLAUSE
}

int mic_check_remote_kernel_complete(const int deviceNum) {
  // DMK - NOTE : Disable these warnings for now since the device code is only called once,
  //   but the state variable goes from 0 -> 1 -> 2 -> 0, which checks at both 1 and 2
  //if (!tag_remote_kernel) { printf("WARNING :: mic_check_remote_kernel_complete :: called when kernel not active.\n"); }
  // if (_Offload_signaled(deviceNum, &tag_remote_kernel)) {
    tag_remote_kernel = 0;
    return 1;
  // }
  // return 0;
}

int mic_check_local_kernel_complete(const int deviceNum) {
  // DMK - NOTE : Disable these warnings for now since the device code is only called once,
  //   but the state variable goes from 0 -> 1 -> 2 -> 0, which checks at both 1 and 2
  //if (!tag_local_kernel) { printf("WARNING :: mic_check_local_kernel_complete :: called when kernel not active.\n"); }
  // if (_Offload_signaled(deviceNum, &tag_local_kernel)) {
    tag_local_kernel = 0;
    return 1;
  // }
  // return 0;
}


void mic_free_device(const int deviceNum) {

  mic_kernel_data * kernel_data = host__kernel_data;

  // Cleanup kernel data (for buffers allocated via offload pragmas, use "free_if(1)")
  #pragma offload target(mic:deviceNum) \
    nocopy(kernel_data : alloc_if(0) free_if(1)) \
    nocopy(device__table_four) nocopy(device__table_four_float) \
    nocopy(device__lj_table) nocopy(device__lj_table_float) \
    nocopy(device__exclusion_bits) nocopy(device__exclusion_bits_copy) \
    nocopy(device__constants) \
    nocopy(device__force_buffers) nocopy(device__slow_force_buffers) \
    nocopy(device__pairlists) nocopy(device__pairlists_alloc_size) \
    nocopy(device__pl_array) nocopy(device__pl_size) nocopy(device__r2_array) \
    nocopy(device__patch_pairs : alloc_if(0) free_if(1)) \
    nocopy(device__force_lists : alloc_if(0) free_if(1)) \
    nocopy(device__numOMPThreads) \
    DEVICE_FPRINTF_CLAUSE
  {
    // Cleanup pairlist memory
    if (device__pairlists != NULL) {
      int **pairlists = (int**)(device__pairlists);
      for (int i = 0; i < device__pairlists_alloc_size; i++) {
        if (pairlists[i] != NULL) { _MM_FREE_WRAPPER(pairlists[i]); }
      }
      _MM_FREE_WRAPPER(pairlists);
      device__pairlists = 0;
      device__pairlists_alloc_size = 0;
    }

    // Cleanup refined pairlist memory
    #if REFINE_PAIRLISTS != 0
      for (int i = 0; i < device__numOMPThreads; i++) {
        if (device__pl_array[i] != NULL) { _MM_FREE_WRAPPER(device__pl_array[i]); }
        if (device__r2_array[i] != NULL) { _MM_FREE_WRAPPER(device__r2_array[i]); }
      }
      _MM_FREE_WRAPPER(device__pl_array); device__pl_array = NULL;
      _MM_FREE_WRAPPER(device__r2_array); device__r2_array = NULL;
      _MM_FREE_WRAPPER(device__pl_size); device__pl_size = NULL;
    #endif

    // Clean up tables
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      if (device__table_four != NULL) { _MM_FREE_WRAPPER(device__table_four); device__table_four = NULL; }
      if (device__table_four_float != NULL) { _MM_FREE_WRAPPER(device__table_four_float); device__table_four_float = NULL; }
      if (device__lj_table != NULL) { _MM_FREE_WRAPPER(device__lj_table); device__lj_table = NULL; }
      if (device__lj_table_float != NULL) { _MM_FREE_WRAPPER(device__lj_table_float); device__lj_table_float = NULL; }
    #endif
    if (device__exclusion_bits != NULL) { _MM_FREE_WRAPPER(device__exclusion_bits); device__exclusion_bits = NULL; }
    if (device__exclusion_bits_copy != NULL) { _MM_FREE_WRAPPER(device__exclusion_bits_copy); device__exclusion_bits_copy = NULL; }
    if (device__constants != NULL) { _MM_FREE_WRAPPER(device__constants); device__constants = NULL; }

    // Clean up intermediate force buffers
    if (device__force_buffers != NULL) { _MM_FREE_WRAPPER(device__force_buffers); device__force_buffers = NULL; }
    if (device__slow_force_buffers != NULL) { _MM_FREE_WRAPPER(device__slow_force_buffers); device__slow_force_buffers = NULL; }

    #if MIC_DEVICE_FPRINTF != 0
      fclose(device__fout);
    #endif
  }

  if (host__kernel_data != NULL) { _MM_FREE_WRAPPER(host__kernel_data); host__kernel_data = NULL; }
}


// DMK - DEBUG
#if MIC_TRACK_DEVICE_MEM_USAGE != 0

__declspec(target(mic))
void printMemInfo(int device__pe, int device__timestep, MemInfo * mi) {
  printf("[MIC-MEMINFO] :: PE %d :: MIC mem usage (MB) -- ts: %d - "
         "f: %.3lf, t: %.3lf, c: %.3lf, a: %.3lf, i: %.3lf - "
         "vS: %.3lf, vP: %.3lf\n",
         ((int)(device__pe)),
         ((int)(device__timestep)),
         ((double)(mi->memFree)) * 1.0e-6,
         ((double)(mi->memTotal)) * 1.0e-6,
         ((double)(mi->cached)) * 1.0e-6,
         ((double)(mi->active)) * 1.0e-6,
         ((double)(mi->inactive)) * 1.0e-6,
         ((double)(mi->vmSize)) * 1.0e-6,
         ((double)(mi->vmPeak)) * 1.0e-6
        ); fflush(NULL);
}

__declspec(target(mic))
void readMemInfo_processLine(MemInfo * memInfo, char * n, char * v, char * u) {
  assert(n != NULL && v != NULL);

  size_t * loc = NULL;
  if (strcmp(n, "MemTotal:") == 0) { loc = &(memInfo->memTotal); }
  if (strcmp(n, "MemFree:") == 0) { loc = &(memInfo->memFree); }
  if (strcmp(n, "Cached:") == 0) { loc = &(memInfo->cached); }
  if (strcmp(n, "Active:") == 0) { loc = &(memInfo->active); }
  if (strcmp(n, "Inactive:") == 0) { loc = &(memInfo->inactive); }
  if (strcmp(n, "VmSize:") == 0) { loc = &(memInfo->vmSize); }
  if (strcmp(n, "VmPeak:") == 0) { loc = &(memInfo->vmPeak); }
  if (loc == NULL) { return; }

  *loc = (size_t)(atoi(v));
  if (u != NULL) {
    if (strcmp(u, "kB") == 0) { (*loc) *= 1024; }
  }
}

__declspec(target(mic))
bool readMemInfo(MemInfo * memInfo) {
  if (memInfo == NULL) { return false; }
  memset(memInfo, 0, sizeof(MemInfo));

  FILE * procMem = fopen("/proc/meminfo", "r");
  if (procMem == NULL) {
    printf("[WARNING] :: Unable to read /proc/meminfo\n");
    return false;
  }
  char line[256], n[128], v[128], u[128];
  while (!feof(procMem)) {
    int numFilled = fscanf(procMem, "%[^\n]\n", line);
    if (numFilled == 0) { printf("empty line\n"); break; }
    numFilled = sscanf(line, "%s%s%s", n, v, u);
    if (numFilled == 3) {
      readMemInfo_processLine(memInfo, n, v, u);
    } else if (numFilled == 2) { // No unit or prev line had name only
      int vLen = strlen(v);
      if (v[vLen-1] == ':') { // Prev was name only (correct)
        memcpy(n, v, vLen + 1);
      }
      readMemInfo_processLine(memInfo, n, v, NULL);
    }
  }
  fclose(procMem);

  pid_t pid = getpid();
  char filename[256];
  sprintf(filename, "/proc/%d/status", (int)pid);
  procMem = fopen(filename, "r");
  if (procMem == NULL) {
    printf("[WARNING] :: Unable to read %s\n", filename);
    return false;
  }
  while (!feof(procMem)) {
    int numFilled = fscanf(procMem, "%[^\n]\n", line);
    if (numFilled == 0) { printf("empty line\n"); break; }
    numFilled = sscanf(line, "%s%s%s", n, v, u);
    if (numFilled == 3) {
      readMemInfo_processLine(memInfo, n, v, u);
    } else if (numFilled == 2) { // No unit
      readMemInfo_processLine(memInfo, n, v, NULL);
    }
  }
  fclose(procMem);

  return true;
}

#endif  // MIC_TRACK_DEVICE_MEM_USAGE != 0


__declspec(target(mic))
void* _mm_malloc_withPrint(size_t s, int a, char * l) {
  void * ptr = _mm_malloc(s, a);
  printf("[MIC-MALLOC] :: _mm_malloc - l: \"%s\", s: %lu (%lu), ptr: %p\n", l, s, sizeof(s), ptr);
  return ptr;
}

__declspec(target(mic))
void _mm_free_withPrint(void * ptr) {
  printf("[MIC-MALLOC] :: _mm_free - ptr: %p\n", ptr);
  _mm_free(ptr);
}


#else  // NAMD_MIC

#include "ComputeNonbondedMICKernel.h"
#include "ComputeNonbondedMICKernelBase.h"

#endif  // NAMD_MIC
