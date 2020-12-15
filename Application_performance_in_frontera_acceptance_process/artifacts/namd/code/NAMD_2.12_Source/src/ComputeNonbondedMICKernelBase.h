
#ifdef NAMD_MIC

////////////////////////////////////////////////////////////////////////////////
///// Begin Macro Definition Region

#ifdef ARCH_POWERPC
  #include <builtins.h>
#endif

#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
  #include <emmintrin.h>  // We're using SSE2 intrinsics
  #if defined(__INTEL_COMPILER)
    #define __align(X) __declspec(align(X) )
  #elif defined(__GNUC__) || defined(__PGI)
    #define __align(X)  __attribute__((aligned(X) ))
  #else
    #define __align(X) __declspec(align(X) )
  #endif
#endif

// DMK - NOTE : Removing the includes for now, since we are using the "raw data"
//   from the host data structures (i.e. packing them up at startup and passing
//   that raw data to the card).  As data structures are included on the device
//   as well (if that happens), then some of these could be included.
//#ifdef DEFINITION // (
//  #include "LJTable.h"
//  #include "Molecule.h"
//  #include "ComputeNonbondedUtil.h"
//#endif // )
//#include "Parameters.h"
//#if NAMD_ComputeNonbonded_SortAtoms != 0
//  #include "PatchMap.h"
//#endif

// determining class name
#undef NAME
#undef CLASS
#undef CLASSNAME
#define NAME CLASSNAME(calc)

#undef PAIR
#if NBTYPE == NBPAIR
  #define PAIR(X) X
  #define CLASS ComputeNonbondedPair
  #define CLASSNAME(X) ENERGYNAME( X ## _pair )
#else
  #define PAIR(X)
#endif

#undef SELF
#if NBTYPE == NBSELF
  #define SELF(X) X
  #define CLASS ComputeNonbondedSelf
  #define CLASSNAME(X) ENERGYNAME( X ## _self )
#else
  #define SELF(X)
#endif

#undef ENERGYNAME
#undef ENERGY
#undef NOENERGY
#ifdef CALCENERGY
  #define ENERGY(X) X
  #define NOENERGY(X)
  #define ENERGYNAME(X) SLOWONLYNAME( X ## _energy )
#else
  #define ENERGY(X)
  #define NOENERGY(X) X
  #define ENERGYNAME(X) SLOWONLYNAME( X )
#endif

#undef SLOWONLYNAME
#undef FAST
#ifdef SLOWONLY
  #define FAST(X)
  #define SLOWONLYNAME(X) MERGEELECTNAME( X ## _slow )
#else
  #define FAST(X) X
  #define SLOWONLYNAME(X) MERGEELECTNAME( X )
#endif

#undef MERGEELECTNAME
#undef SHORT
#undef NOSHORT
#ifdef MERGEELECT
  #define SHORT(X)
  #define NOSHORT(X) X
  #define MERGEELECTNAME(X) FULLELECTNAME( X ## _merge )
#else
  #define SHORT(X) X
  #define NOSHORT(X)
  #define MERGEELECTNAME(X) FULLELECTNAME( X )
#endif

#undef FULLELECTNAME
#undef FULL
#undef NOFULL
#ifdef FULLELECT
  #define FULLELECTNAME(X) TABENERGYNAME( X ## _fullelect )
  #define FULL(X) X
  #define NOFULL(X)
#else
  #define FULLELECTNAME(X) TABENERGYNAME( X )
  #define FULL(X)
  #define NOFULL(X) X
#endif

#undef TABENERGYNAME
#undef TABENERGY
#undef NOTABENERGY
#ifdef TABENERGYFLAG
  #define TABENERGYNAME(X) FEPNAME( X ## _tabener )
  #define TABENERGY(X) X
  #define NOTABENERGY(X)
#else
  #define TABENERGYNAME(X) FEPNAME( X )
  #define TABENERGY(X)
  #define NOTABENERGY(X) X
#endif

// Here are what these macros stand for:
// FEP/NOT_FEP: FEP free energy perturbation is active/inactive 
//      (does NOT use LAM)
// LES: locally-enhanced sampling is active
// LAM: scale the potential by a factor lambda (except FEP)
// INT: measure interaction energies
// PPROF: pressure profiling

#undef FEPNAME
#undef FEP
#undef LES
#undef INT
#undef PPROF
#undef LAM
#undef CUDA
#undef MIC
#undef ALCH
#undef TI
#undef GO
#define FEPNAME(X) LAST( X )
#define FEP(X)
#define ALCHPAIR(X)
#define NOT_ALCHPAIR(X) X
#define LES(X)
#define INT(X)
#define PPROF(X)
#define LAM(X)
#define CUDA(X)
#define MIC(X)
#define ALCH(X)
#define TI(X)
#define GO(X)
#ifdef FEPFLAG
  #undef FEPNAME
  #undef FEP
  #undef ALCH
  #define FEPNAME(X) LAST( X ## _fep )
  #define FEP(X) X
  #define ALCH(X) X
#endif
#ifdef TIFLAG
  #undef FEPNAME
  #undef TI
  #undef ALCH
  #define FEPNAME(X) LAST( X ## _ti )
  #define TI(X) X
  #define ALCH(X) X
#endif
#ifdef LESFLAG
  #undef FEPNAME
  #undef LES
  #undef LAM
  #define FEPNAME(X) LAST( X ## _les )
  #define LES(X) X
  #define LAM(X) X
#endif
#ifdef INTFLAG
  #undef FEPNAME
  #undef INT
  #define FEPNAME(X) LAST( X ## _int )
  #define INT(X) X
#endif
#ifdef PPROFFLAG
  #undef FEPNAME
  #undef INT
  #undef PPROF
  #define FEPNAME(X) LAST( X ## _pprof )
  #define INT(X) X
  #define PPROF(X) X
#endif
#ifdef GOFORCES
  #undef FEPNAME
  #undef GO
  #define FEPNAME(X) LAST( X ## _go )
  #define GO(X) X
#endif
#ifdef NAMD_CUDA
  #undef CUDA
  #define CUDA(X) X
#endif
#ifdef NAMD_MIC
  #undef MIC
  #define MIC(X) X
#endif

#define LAST(X) X

// see if things are really messed up
SELF( PAIR( foo bar ) )
LES( FEP( foo bar ) )
LES( INT( foo bar ) )
FEP( INT( foo bar ) )
LAM( INT( foo bar ) )
FEP( NOENERGY( foo bar ) )
ENERGY( NOENERGY( foo bar ) )
TABENERGY(NOTABENERGY( foo bar ) )

#define ALLOCA_ALIGNED(p, s, a) { \
  char *ptr = alloca((s) + (a)); \
  int64 a1 = (a) - 1; \
  int64 ptr64 = ((int64)ptr) + a1); \
  (p) = (void*)(ptr64 - (ptr64 & a1)); \
}

#define __MIC_PAD_PLGEN_CTRL (MIC_PAD_PLGEN)

///// End Macro Definition Region
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Declare the compute function

__attribute__((target(mic))) void NAME (mic_params &params) {

  ///// Setup Various Arrays/Values Required by this Compute Function /////

  #if MIC_HANDCODE_FORCE_SINGLE != 0
    #define CALC_TYPE float
  #else
    #define CALC_TYPE double
  #endif

  // Setup nonbonded table pointers (using the main table pointer, create pointers for each
  //   of the non-overlapping sub-tables)
  const int table_four_n_16 = params.table_four_n_16;
  __ASSUME(table_four_n_16 % 16 == 0);
  __ASSERT(params.table_four_base_ptr != NULL);
  const CALC_TYPE * RESTRICT const table_noshort = (CALC_TYPE*)(params.table_four_base_ptr)                       ; __ASSUME_ALIGNED(table_noshort);
  const CALC_TYPE * RESTRICT const   table_short = (CALC_TYPE*)(params.table_four_base_ptr) + 16 * table_four_n_16; __ASSUME_ALIGNED(table_short);
  const CALC_TYPE * RESTRICT const    slow_table = (CALC_TYPE*)(params.table_four_base_ptr) + 32 * table_four_n_16; __ASSUME_ALIGNED(slow_table);
  //const CALC_TYPE * RESTRICT const    fast_table = (CALC_TYPE*)(params.table_four_base_ptr) + 36 * table_four_n_16; __ASSUME_ALIGNED(fast_table);
  //const CALC_TYPE * RESTRICT const    scor_table = (CALC_TYPE*)(params.table_four_base_ptr) + 40 * table_four_n_16; __ASSUME_ALIGNED(scor_table);
  //const CALC_TYPE * RESTRICT const    corr_table = (CALC_TYPE*)(params.table_four_base_ptr) + 44 * table_four_n_16; __ASSUME_ALIGNED(corr_table);
  //const CALC_TYPE * RESTRICT const    full_table = (CALC_TYPE*)(params.table_four_base_ptr) + 48 * table_four_n_16; __ASSUME_ALIGNED(full_table);
  //const CALC_TYPE * RESTRICT const    vdwa_table = (CALC_TYPE*)(params.table_four_base_ptr) + 52 * table_four_n_16; __ASSUME_ALIGNED(vdwa_table);
  //const CALC_TYPE * RESTRICT const    vdwb_table = (CALC_TYPE*)(params.table_four_base_ptr) + 56 * table_four_n_16; __ASSUME_ALIGNED(vdwb_table);
  //const CALC_TYPE * RESTRICT const      r2_table = (CALC_TYPE*)(params.table_four_base_ptr) + 60 * table_four_n_16; __ASSUME_ALIGNED(r2_table);
  // DMK - NOTE : Other table pointers will be useful as this code grows in scope, so including them all now

  // Setup LJ table pointers
  const CALC_TYPE * RESTRICT const lj_table_base_ptr = (CALC_TYPE*)(params.lj_table_base_ptr);
  __ASSERT(lj_table_base_ptr != NULL);
  const int lj_table_dim = params.lj_table_dim;
  __ASSUME_ALIGNED(lj_table_base_ptr);

  // Constants
  const CALC_TYPE r2_delta = params.constants->r2_delta;
  const int r2_delta_exp = params.constants->r2_delta_exp;
  const int r2_delta_expc = params.constants->r2_delta_expc;
  const CALC_TYPE scaling = params.constants->scaling;
  const CALC_TYPE modf_mod = params.constants->modf_mod;

  // Cutoff values
  const CALC_TYPE plcutoff = params.pp->plcutoff;
  const CALC_TYPE plcutoff2 = plcutoff * plcutoff;
  const CALC_TYPE plcutoff2_delta = plcutoff2 + r2_delta;
  const CALC_TYPE cutoff2 = params.constants->cutoff2;
  const CALC_TYPE cutoff2_delta = cutoff2 + r2_delta;

  // Number of atoms in each list of atoms
  const int i_upper = params.numAtoms[0];
  const int j_upper = params.numAtoms[1];
  const int i_upper_16 = params.numAtoms_16[0];
  const int j_upper_16 = params.numAtoms_16[1];
  __ASSERT(i_upper >= 0); __ASSERT(j_upper >= 0);
  __ASSERT(i_upper_16 >= 0); __ASSERT(j_upper_16 >= 0);
  __ASSUME(i_upper_16 % 16 == 0);
  __ASSUME(j_upper_16 % 16 == 0);

  // Setup pointers to atom input arrays
  #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
    atom * RESTRICT p_0 = params.p[0]; __ASSERT(p_0 != NULL);  __ASSUME_ALIGNED(p_0);
    atom_param * RESTRICT pExt_0 = params.pExt[0]; __ASSERT(pExt_0 != NULL);  __ASSUME_ALIGNED(pExt_0);
    #if (0 SELF(+1))
      #define p_1  p_0
      #define pExt_1  pExt_0
    #else
      atom * RESTRICT p_1 = params.p[1]; __ASSERT(p_1 != NULL);  __ASSUME_ALIGNED(p_1);
      atom_param * RESTRICT pExt_1 = params.pExt[1]; __ASSERT(pExt_1 != NULL);  __ASSUME_ALIGNED(pExt_1);
    #endif
  #else
    mic_position_t * RESTRICT p_0_x = params.p_x[0];
    mic_position_t * RESTRICT p_0_y = params.p_y[0];
    mic_position_t * RESTRICT p_0_z = params.p_z[0];
    mic_position_t * RESTRICT p_0_q = params.p_q[0];
    int * RESTRICT pExt_0_vdwType = params.pExt_vdwType[0];
    int * RESTRICT pExt_0_index = params.pExt_index[0];
    int * RESTRICT pExt_0_exclIndex = params.pExt_exclIndex[0];
    int * RESTRICT pExt_0_exclMaxDiff = params.pExt_exclMaxDiff[0];
    #if (0 SELF(+1))
      #define p_1_x  p_0_x
      #define p_1_y  p_0_y
      #define p_1_z  p_0_z
      #define p_1_q  p_0_q
      #define p_1_vdwType  p_0_vdwType
      #define p_1_index  p_0_index
      #define p_1_exclIndex  p_0_exclIndex
      #define p_1_exclMaxDiff  p_0_exclMaxDiff
    #else
      mic_position_t * RESTRICT p_1_x = params.p_x[1];
      mic_position_t * RESTRICT p_1_y = params.p_y[1];
      mic_position_t * RESTRICT p_1_z = params.p_z[1];
      mic_position_t * RESTRICT p_1_q = params.p_q[1];
      int * RESTRICT pExt_1_vdwType = params.pExt_vdwType[1];
      int * RESTRICT pExt_1_index = params.pExt_index[1];
      int * RESTRICT pExt_1_exclIndex = params.pExt_exclIndex[1];
      int * RESTRICT pExt_1_exclMaxDiff = params.pExt_exclMaxDiff[1];
    #endif
    __ASSERT(p_0_x != NULL); __ASSERT(p_0_y != NULL); __ASSERT(p_0_z != NULL); __ASSERT(p_0_q != NULL);
    __ASSERT(p_1_x != NULL); __ASSERT(p_1_y != NULL); __ASSERT(p_1_z != NULL); __ASSERT(p_1_q != NULL);
    __ASSERT(pExt_0_vdwType != NULL); __ASSERT(pExt_0_index != NULL);
    __ASSERT(pExt_1_vdwType != NULL); __ASSERT(pExt_1_index != NULL);
    __ASSERT(pExt_0_exclIndex != NULL); __ASSERT(pExt_0_exclMaxDiff != NULL);
    __ASSERT(pExt_1_exclIndex != NULL); __ASSERT(pExt_1_exclMaxDiff != NULL);
    __ASSUME_ALIGNED(p_0_x); __ASSUME_ALIGNED(p_1_x);
    __ASSUME_ALIGNED(p_0_y); __ASSUME_ALIGNED(p_1_y);
    __ASSUME_ALIGNED(p_0_z); __ASSUME_ALIGNED(p_1_z);
    __ASSUME_ALIGNED(p_0_q); __ASSUME_ALIGNED(p_1_q);
    __ASSUME_ALIGNED(pExt_0_vdwType); __ASSUME_ALIGNED(pExt_0_index);
    __ASSUME_ALIGNED(pExt_1_vdwType); __ASSUME_ALIGNED(pExt_1_index);
    __ASSUME_ALIGNED(pExt_0_exclIndex); __ASSUME_ALIGNED(pExt_0_exclMaxDiff);
    __ASSUME_ALIGNED(pExt_1_exclIndex); __ASSUME_ALIGNED(pExt_1_exclMaxDiff);
  #endif

  // Setup pointers to force output arrays and clear those arrays (init to zero)
  #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
    double4 * RESTRICT f_0 = params.ff[0];
    #if (0 SELF(+1))
      #define f_1  f_0
    #else
      double4 * RESTRICT f_1 = params.ff[1];
    #endif
    __ASSERT(f_0 != NULL); __ASSUME_ALIGNED(f_0);
    __ASSERT(f_1 != NULL); __ASSUME_ALIGNED(f_1);
    memset(f_0, 0, sizeof(double4) * i_upper_16);
    PAIR( memset(f_1, 0, sizeof(double4) * j_upper_16); )
    #if (0 FULL(+1))
      double4 * RESTRICT fullf_0 = params.fullf[0];
      #if (0 SELF(+1))
        #define fullf_1  fullf_0
      #else
        double4 * RESTRICT fullf_1 = params.fullf[1];
      #endif
      __ASSERT(fullf_0 != NULL); __ASSUME_ALIGNED(fullf_0);
      __ASSERT(fullf_1 != NULL); __ASSUME_ALIGNED(fullf_1);
      memset(fullf_0, 0, sizeof(double4) * i_upper_16);
      PAIR( memset(fullf_1, 0, sizeof(double4) * j_upper_16); )
    #endif
  #else
    double * RESTRICT f_0_x = params.ff_x[0];
    double * RESTRICT f_0_y = params.ff_y[0];
    double * RESTRICT f_0_z = params.ff_z[0];
    double * RESTRICT f_0_w = params.ff_w[0];
    #if (0 SELF(+1))
      #define f_1_x  f_0_x
      #define f_1_y  f_0_y
      #define f_1_z  f_0_z
      #define f_1_w  f_0_w
    #else
      double * RESTRICT f_1_x = params.ff_x[1];
      double * RESTRICT f_1_y = params.ff_y[1];
      double * RESTRICT f_1_z = params.ff_z[1];
      double * RESTRICT f_1_w = params.ff_w[1];
    #endif
    __ASSERT(f_0_x != NULL); __ASSERT(f_0_y != NULL); __ASSERT(f_0_z != NULL); __ASSERT(f_0_w != NULL);
    __ASSERT(f_1_x != NULL); __ASSERT(f_1_y != NULL); __ASSERT(f_1_z != NULL); __ASSERT(f_1_w != NULL);
    __ASSUME_ALIGNED(f_0_x); __ASSUME_ALIGNED(f_1_x);
    __ASSUME_ALIGNED(f_0_y); __ASSUME_ALIGNED(f_1_y);
    __ASSUME_ALIGNED(f_0_z); __ASSUME_ALIGNED(f_1_z);
    __ASSUME_ALIGNED(f_0_w); __ASSUME_ALIGNED(f_1_w);
    memset(f_0_x, 0, 4 * sizeof(double) * i_upper_16);
    PAIR( memset(f_1_x, 0, 4 * sizeof(double) * j_upper_16); )
    #if (0 FULL(+1))
      double * RESTRICT fullf_0_x = params.fullf_x[0];
      double * RESTRICT fullf_0_y = params.fullf_y[0];
      double * RESTRICT fullf_0_z = params.fullf_z[0];
      double * RESTRICT fullf_0_w = params.fullf_w[0];
      #if (0 SELF(+1))
        #define fullf_1_x  fullf_0_x
        #define fullf_1_y  fullf_0_y
        #define fullf_1_z  fullf_0_z
        #define fullf_1_w  fullf_0_w
      #else
        double * RESTRICT fullf_1_x = params.fullf_x[1];
        double * RESTRICT fullf_1_y = params.fullf_y[1];
        double * RESTRICT fullf_1_z = params.fullf_z[1];
        double * RESTRICT fullf_1_w = params.fullf_w[1];
      #endif
      __ASSERT(fullf_0_x != NULL); __ASSERT(fullf_0_y != NULL); __ASSERT(fullf_0_z != NULL); __ASSERT(fullf_0_w != NULL);
      __ASSERT(fullf_1_x != NULL); __ASSERT(fullf_1_y != NULL); __ASSERT(fullf_1_z != NULL); __ASSERT(fullf_1_w != NULL);
      __ASSUME_ALIGNED(fullf_0_x); __ASSUME_ALIGNED(fullf_1_x);
      __ASSUME_ALIGNED(fullf_0_y); __ASSUME_ALIGNED(fullf_1_y);
      __ASSUME_ALIGNED(fullf_0_z); __ASSUME_ALIGNED(fullf_1_z);
      __ASSUME_ALIGNED(fullf_0_w); __ASSUME_ALIGNED(fullf_1_w);
      memset(fullf_0_x, 0, 4 * sizeof(double) * i_upper_16);
      PAIR( memset(fullf_1_x, 0, 4 * sizeof(double) * j_upper_16); )
    #endif
  #endif

  // If the commOnly flag has been set, return now (after having zero'd the output forces)
  if (params.constants->commOnly) { return; }

  // If either atom list is size zero, return now (after having zero'd the output forces)
  if (i_upper <= 0 || j_upper <= 0) { return; }

  CALC_TYPE offset_x = (CALC_TYPE)(params.offset.x);
  CALC_TYPE offset_y = (CALC_TYPE)(params.offset.y);
  CALC_TYPE offset_z = (CALC_TYPE)(params.offset.z);

  ///// Generate the Pairlists /////

  // Grab the pairlist pointers for the various pairlists that will be used
  int * RESTRICT pairlist_norm = params.pairlists_ptr[PL_NORM_INDEX];
  int * RESTRICT pairlist_mod  = params.pairlists_ptr[PL_MOD_INDEX];
  int * RESTRICT pairlist_excl = params.pairlists_ptr[PL_EXCL_INDEX];
  #if MIC_CONDITION_NORMAL != 0
    int * RESTRICT pairlist_normhi = params.pairlists_ptr[PL_NORMHI_INDEX];
  #endif

  // NOTE : The first 2 integers are used for sizing information.  The actual
  //   pairlist entries (i.e. what we want to be 64 byte aligned) start at
  //   element 2 (thus the +/- 14 shifts in the macros below).
  #define PAIRLIST_ALLOC_CHECK(pl, init_size) { \
    if ((pl) == NULL) { \
      const int plInitSize = (init_size) + (16 * (MIC_HANDCODE_FORCE_PFDIST + 1)); \
      if (device__timestep == 0) { /* NOTE: Avoid all the prints on timestep 0 */ \
        (pl) = (int*)(_mm_malloc(plInitSize * sizeof(int) + 64, 64)); \
      } else { \
        (pl) = (int*)_MM_MALLOC_WRAPPER(plInitSize * sizeof(int) + 64, 64, "pairlist_alloc_check"); \
      } \
      __ASSERT((pl) != NULL); \
      (pl) += 14; \
      (pl)[0] = plInitSize; \
      (pl)[1] = 2; \
    } \
  }

  #define PAIRLIST_GROW_CHECK(pl, offset, mul) { \
    if ((offset) + j_upper + 8 + MIC_PREFETCH_DISTANCE + (16 * (MIC_HANDCODE_FORCE_PFDIST + 1)) > (pl)[0]) { \
      int newPlAllocSize = (int)(((pl)[0]) * (mul) + j_upper + 64); /* NOTE: add at least j_upper (max that can be added before next check) in case pl[0] is small to begin with */ \
      newPlAllocSize = (~2047) & (newPlAllocSize + 2047); /* NOTE: round up to 2K, add at least 2K */ \
      int* RESTRICT newPl = (int*)_MM_MALLOC_WRAPPER(newPlAllocSize * sizeof(int), 64, "pairlist_grow_check"); \
      __ASSERT((newPl) != NULL); \
      newPl += 14; \
      memcpy(newPl, (pl), (offset) * sizeof(int)); \
      _MM_FREE_WRAPPER((pl) - 14); \
      (pl) = newPl; \
      (pl)[0] = newPlAllocSize; \
    } \
  }

  DEVICE_FPRINTF("P");

  // Regenerate the pairlists, if need be
  if ((!(params.usePairlists)) || (params.savePairlists)) {

    // For the sake of getting a good estimate of the virtual memory require for pairlists,
    //   do an actual count of the interactions within this compute (part).
    // NOTE: This will only occur during timestep 0.  Later allocations grow the buffer by
    //   a fixed factor.

    if (pairlist_norm == NULL) { // NOTE: Only checking one pointer, but all will be allocated together
      // NOTE: This only occurs once when the compute is run for the first time

      int plCount_norm = 16;
      int plCount_mod = 16;
      int plCount_excl = 16;
      #if MIC_CONDITION_NORMAL != 0
        int plCount_normhi = 16;
      #endif

      int numParts = params.pp->numParts;
      int part = params.pp->part;
      for (int i = part; i < i_upper; i += numParts) {

        // Load position information for the current "i" atom
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          CALC_TYPE p_i_x = p_0[i].x + offset_x;
          CALC_TYPE p_i_y = p_0[i].y + offset_y;
          CALC_TYPE p_i_z = p_0[i].z + offset_z;
        #else
          CALC_TYPE p_i_x = p_0_x[i] + offset_x;
          CALC_TYPE p_i_y = p_0_y[i] + offset_y;
          CALC_TYPE p_i_z = p_0_z[i] + offset_z;
        #endif
          
        for (int j = 0 SELF(+i+1); j < j_upper; j++) {

          // Load the "j" atom's position, calculate/check the distance (squared)
          //   between the "i" and "j" atoms
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            CALC_TYPE p_d_x = p_i_x - p_1[j].x;
            CALC_TYPE p_d_y = p_i_y - p_1[j].y;
            CALC_TYPE p_d_z = p_i_z - p_1[j].z;
          #else
            CALC_TYPE p_d_x = p_i_x - p_1_x[j];
            CALC_TYPE p_d_y = p_i_y - p_1_y[j];
            CALC_TYPE p_d_z = p_i_z - p_1_z[j];
          #endif
          CALC_TYPE r2 = (p_d_x * p_d_x) + (p_d_y * p_d_y) + (p_d_z * p_d_z);
          if (r2 <= plcutoff2) {

            // Check the exclusion bits to set the modified and excluded flags
            // NOTE: This code MUST match the exclusion lists generation code
            //   on ComputeNonbondedMIC.C, which generates the flags.  In
            //   particular, the order of the 2 bits must be correct, per the
            //   pairlist format listed below (i.e. isModified -> MSB).
            int exclFlags = 0x00;
            #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
              const int indexDiff = pExt_0[i].index - pExt_1[j].index;
              const int maxDiff = pExt_1[j].excl_maxdiff;
            #else
              const int indexDiff = pExt_0_index[i] - pExt_1_index[j];
              const int maxDiff = pExt_1_exclMaxDiff[j];
            #endif
            if (indexDiff >= -1 * maxDiff && indexDiff <= maxDiff) {
              #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                const int offset = (2 * indexDiff) + pExt_1[j].excl_index;
              #else
                const int offset = (2 * indexDiff) + pExt_1_exclIndex[j];
              #endif
              const int offset_major = offset / (sizeof(unsigned int) * 8);
              const int offset_minor = offset % (sizeof(unsigned int) * 8); // NOTE: Reverse indexing direction relative to offset_major
              exclFlags = ((params.exclusion_bits[offset_major]) >> offset_minor) & 0x03;
	    }

            // Create the pairlist entry value (i and j value) and store it in to the
            //   appropriate pairlist based on the type of interaction (exclFlags value)
            if (exclFlags == 0) {
              #if MIC_CONDITION_NORMAL != 0
                if (r2 <= normhi_split) { plCount_norm++; }
                else { plCount_normhi++; }
              #else
                plCount_norm++;
              #endif
            } else if (exclFlags == 1) {
              plCount_excl++;
            } else if (exclFlags == 2) {
              plCount_mod++;
            }
              
          } // end if (r2 <= plcutoff2)
        } // end for (j < j_upper)
      } // end for (i < i_upper)

      // If padding the pairlists, add some extra room for padding ('vector width * num sub-lists' total)
      #if __MIC_PAD_PLGEN_CTRL != 0
        plCount_norm += (16 * i_upper);
        plCount_mod += (16 * i_upper);
        plCount_excl += (16 * i_upper);
        #if MIC_CONDITION_NORMAL != 0
          plCount_normhi += (16 * i_upper);
        #endif
      #endif

      // Add a constant percent and round up to a multiple of 1024
      const float plPercent = 1.4;
      const int plRound = 2048;
      plCount_norm = (int)(plCount_norm * plPercent); plCount_norm = (~(plRound-1)) & (plCount_norm + (plRound-1));
      plCount_mod = (int)(plCount_mod * plPercent); plCount_mod = (~(plRound-1)) & (plCount_mod + (plRound-1));
      plCount_excl = (int)(plCount_excl * plPercent); plCount_excl = (~(plRound-1)) & (plCount_excl + (plRound-1));
      #if MIC_CONDITION_NORMAL != 0
        plCount_normhi = (int)(plCount_normhi * plPercent); plCount_normhi = (~(plRound-1)) & (plCount_normhi + (plRound-1));
        plCount_norm += plCount_normhi;
      #endif

      // DMK - DEBUG - printing allocations, but reduce output by skipping initial pairlist allocations in timestep 0 (lots of them)
      // NOTE: This check is temporary for debugging, remove (keep wrapper version) when no longer needed
      if (device__timestep == 0) { /* NOTE: Avoid all the allocation prints during timestep 0 */
        pairlist_norm = (int*)(_mm_malloc(plCount_norm * sizeof(int), 64)); __ASSERT(pairlist_norm != NULL);
        pairlist_mod = (int*)(_mm_malloc(plCount_mod * sizeof(int), 64)); __ASSERT(pairlist_mod != NULL);
        pairlist_excl = (int*)(_mm_malloc(plCount_excl * sizeof(int), 64)); __ASSERT(pairlist_excl != NULL);
        #if MIC_CONDITION_NORMAL != 0
          pairlist_normhi = (int*)(_mm_malloc(plCount_normhi * sizeof(int), 64)); __ASSERT(pairlist_normhi != NULL);
        #endif
      } else {
        pairlist_norm = (int*)(_MM_MALLOC_WRAPPER(plCount_norm * sizeof(int), 64, "pairlist_norm")); __ASSERT(pairlist_norm != NULL);
        pairlist_mod = (int*)(_MM_MALLOC_WRAPPER(plCount_mod * sizeof(int), 64, "pairlist_mod")); __ASSERT(pairlist_mod != NULL);
        pairlist_excl = (int*)(_MM_MALLOC_WRAPPER(plCount_excl * sizeof(int), 64, "pairlist_excl")); __ASSERT(pairlist_excl != NULL);
        #if MIC_CONDITION_NORMAL != 0
          pairlist_normhi = (int*)(_MM_MALLOC_WRAPPER(plCount_normhi * sizeof(int), 64, "pairlist_normhi")); __ASSERT(pairlist_normhi != NULL);
        #endif
      }

      pairlist_norm += 14; pairlist_norm[0] = plCount_norm; pairlist_norm[1] = 2;
      pairlist_mod += 14; pairlist_mod[0] = plCount_mod; pairlist_mod[1] = 2;
      pairlist_excl += 14; pairlist_excl[0] = plCount_excl; pairlist_excl[1] = 2;
      #if MIC_CONDITION_NORMAL != 0
        pairlist_normhi += 14; pairlist_normhi[0] = plCount_normhi; pairlist_normhi[1] = 2;
      #endif

      // DMK - DEBUG - Pairlist memory stats
      #if MIC_TRACK_DEVICE_MEM_USAGE != 0
        pairlist_norm[-1] = 0;
        pairlist_mod[-1] = 0;
        pairlist_excl[-1] = 0;
        #if MIC_CONDITION_NORMAL != 0
          pairlist_normhi[-1] = 0;
        #endif
      #endif

    } // end if (pairlist_norm == NULL)

    // Create the pairlists from scratch.  The first two elements in each pairlist
    //   are special values (element 0 is the allocated size of the memory buffer
    //   and element 1 is the length of the pairlist itself, including the these
    //   first two elements).
    int plOffset_norm = 2;  // Current length of the normal pairlist
    int plOffset_mod = 2;   // Current length of the modified pairlist
    int plOffset_excl = 2;  // Current length of the excluded pairlist
    #if MIC_CONDITION_NORMAL != 0
      double normhi_split = (0.25 * cutoff2) + (0.75 * plcutoff2);  // Weighted average of cutoff2 and plcutoff2
      int plOffset_normhi = 2;  // Current length of the "normal high" pairlist (entires
                                //   that belong in the normal pairlist, but that are
                                //   further than the normhi_split distance, such that
                                //   cutoff2 <= normhi_split <= plcutoff2).
    #endif

    // If tiling of the pairlist is enabled...
    #if MIC_TILE_PLGEN != 0

      // Calculate the loop tiling size
      int plTileSize = (i_upper > j_upper) ? (i_upper) : (j_upper); //(i_upper + j_upper) / 2;  // Avg. of list sizes ...
      plTileSize /= MIC_TILE_PLGEN;              // ... divided by MIC_TILE_PLGEN and ...
      plTileSize = (plTileSize + 15) & (~15);    // ... rounded up to multiple of 16
      __ASSUME(plTileSize % 16 == 0);

      // The "i" and "j" tile loops
      for (int _i = 0; _i < i_upper SELF(-1); _i += plTileSize) {
        for (int _j = 0 SELF(+ _i); _j < j_upper; _j += plTileSize) {

          // Calculate the "i" loop bounds for the current tile
          int i_lo = _i;
          int i_hi = _i + plTileSize;
          i_hi = (i_hi > i_upper SELF(-1)) ? (i_upper SELF(-1)) : (i_hi);

          // The "i" loop, iterating over the first list of atoms
          for (int i = i_lo; i < i_hi; i++) {

            // If the pairlist is not long enough, grow the pairlist
            PAIRLIST_GROW_CHECK(pairlist_norm, plOffset_norm, 1.4);
            PAIRLIST_GROW_CHECK(pairlist_mod , plOffset_mod,  1.2);
            PAIRLIST_GROW_CHECK(pairlist_excl, plOffset_excl, 1.2);
            #if MIC_CONDITION_NORMAL != 0
              PAIRLIST_GROW_CHECK(pairlist_normhi, plOffset_normhi, 1.3);
            #endif

	    // Load position information for the current "i" atom
            #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
              CALC_TYPE p_i_x = p_0[i].x + offset_x;
              CALC_TYPE p_i_y = p_0[i].y + offset_y;
              CALC_TYPE p_i_z = p_0[i].z + offset_z;
            #else
              CALC_TYPE p_i_x = p_0_x[i] + offset_x; //params.offset.x;
              CALC_TYPE p_i_y = p_0_y[i] + offset_y; //params.offset.y;
              CALC_TYPE p_i_z = p_0_z[i] + offset_z; //params.offset.z;
            #endif

            // Calculate the "j" loop bounds for the current tile
            int j_lo = PAIR(_j) SELF((_i == _j) ? (i+1) : (_j));
            int j_hi = _j + plTileSize;
            j_hi = (j_hi > j_upper) ? (j_upper) : (j_hi);

            #if MIC_HANDCODE_PLGEN != 0

              __m512 p_i_x_vec = _mm512_set_1to16_ps(p_i_x);
              __m512 p_i_y_vec = _mm512_set_1to16_ps(p_i_y);
              __m512 p_i_z_vec = _mm512_set_1to16_ps(p_i_z);
              __m512i pExt_i_index_vec = _mm512_set_1to16_epi32(pExt_0_index[i]);

              __m512i i_shifted_vec = _mm512_slli_epi32(_mm512_set_1to16_epi32(i), 16); 
              __m512i j_vec = _mm512_set_16to16_epi32(15, 14, 13, 12, 11, 10,  9,  8,
                                                       7,  6,  5,  4,  3,  2,  1,  0);
              // For self computes, round (i+1) down to nearest multiple of 16 and then
              //   add that value to j_vec (skipping iterations where all j values
              //   are < i+1)
              PAIR( j_vec = _mm512_add_epi32(j_vec, _mm512_set_1to16_epi32(j_lo)); )
              #if (0 SELF(+1))
                const int j_lo_16 = (j_lo) & (~0x0F);
                j_vec = _mm512_add_epi32(j_vec, _mm512_set_1to16_epi32(j_lo_16));
                const int j_lo_m1 = j_lo - 1;
              #endif

              // The "j" loop, iterating over the second list of atoms
              #pragma novector
              #pragma loop count(35)
	      for (int j = PAIR(j_lo) SELF(j_lo_16); j < j_hi; j += 16) {
      
                // Create the active_mask
                __mmask16 active_mask = _mm512_cmplt_epi32_mask(j_vec, _mm512_set_1to16_epi32(j_hi));
                #if (0 SELF(+1))
                  active_mask = _mm512_kand(active_mask, _mm512_cmpgt_epi32_mask(j_vec, _mm512_set_1to16_epi32(j_lo_m1)));
                #endif

                // Create the pairlist entry values for these would-be interactions before advancing
                //   the j_vec values (i in upper 16 bits and j in lower 16 bits)
                __m512i plEntry_vec = _mm512_or_epi32(i_shifted_vec, j_vec);

                // Advance j_vec counter
                j_vec = _mm512_add_epi32(j_vec, _mm512_set_1to16_epi32(16));

                // Load the positions for the "j" atoms
                #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                  __m512i p_j_tmp0a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j     )); // Load  first set of 4 atoms
                  __m512i p_j_tmp1a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j +  4)); // Load second set of 4 atoms
                  __m512i p_j_tmp2a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j +  8)); // Load  third set of 4 atoms
                  __m512i p_j_tmp3a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j + 12)); // Load fourth set of 4 atoms
                  __mmask16 k_2x2_0 = _mm512_int2mask(0xAAAA);
                  __mmask16 k_2x2_1 = _mm512_int2mask(0x5555);
                  __m512i p_j_tmp0b_vec = _mm512_mask_swizzle_epi32(p_j_tmp0a_vec, k_2x2_0, p_j_tmp1a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i p_j_tmp1b_vec = _mm512_mask_swizzle_epi32(p_j_tmp1a_vec, k_2x2_1, p_j_tmp0a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i p_j_tmp2b_vec = _mm512_mask_swizzle_epi32(p_j_tmp2a_vec, k_2x2_0, p_j_tmp3a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i p_j_tmp3b_vec = _mm512_mask_swizzle_epi32(p_j_tmp3a_vec, k_2x2_1, p_j_tmp2a_vec, _MM_SWIZ_REG_CDAB);
                  __mmask16 k_4x4_0 = _mm512_int2mask(0xCCCC);
                  __mmask16 k_4x4_1 = _mm512_int2mask(0x3333);
                  __m512i p_j_tmp0c_vec = _mm512_mask_swizzle_epi32(p_j_tmp0b_vec, k_4x4_0, p_j_tmp2b_vec, _MM_SWIZ_REG_BADC);
                  __m512i p_j_tmp1c_vec = _mm512_mask_swizzle_epi32(p_j_tmp1b_vec, k_4x4_0, p_j_tmp3b_vec, _MM_SWIZ_REG_BADC);
                  __m512i p_j_tmp2c_vec = _mm512_mask_swizzle_epi32(p_j_tmp2b_vec, k_4x4_1, p_j_tmp0b_vec, _MM_SWIZ_REG_BADC);
                  __m512i p_j_perm_pattern = _mm512_set_1to16_epi32(15, 11, 7, 3,  14, 10, 6, 2,  13, 9, 5, 1,  12, 8, 4, 0);
                  __m512 p_j_x_vec = _mm512_castsi512_ps(_mm512_permutevar_epi32(p_j_perm_pattern, p_j_tmp0c_vec));
                  __m512 p_j_y_vec = _mm512_castsi512_ps(_mm512_permutevar_epi32(p_j_perm_pattern, p_j_tmp1c_vec));
                  __m512 p_j_z_vec = _mm512_castsi512_ps(_mm512_permutevar_epi32(p_j_perm_pattern, p_j_tmp2c_vec));
                #else
                  __m512 p_j_x_vec = _mm512_mask_load_ps(_mm512_setzero_ps(), active_mask, p_1_x + j);
                  __m512 p_j_y_vec = _mm512_mask_load_ps(_mm512_setzero_ps(), active_mask, p_1_y + j);
                  __m512 p_j_z_vec = _mm512_mask_load_ps(_mm512_setzero_ps(), active_mask, p_1_z + j);
                #endif

                // Calculate the distance between "i" atom and these "j" atoms (in double)
                __m512 p_ij_x_vec = _mm512_sub_ps(p_i_x_vec, p_j_x_vec);
                __m512 p_ij_y_vec = _mm512_sub_ps(p_i_y_vec, p_j_y_vec);
                __m512 p_ij_z_vec = _mm512_sub_ps(p_i_z_vec, p_j_z_vec);
                __m512 r2_vec = _mm512_mul_ps(p_ij_x_vec, p_ij_x_vec);
                r2_vec = _mm512_add_ps(r2_vec, _mm512_mul_ps(p_ij_y_vec, p_ij_y_vec));
                r2_vec = _mm512_add_ps(r2_vec, _mm512_mul_ps(p_ij_z_vec, p_ij_z_vec));

                // Do cutoff distance check.  If there are no particles within cutoff, move on
                //   to the next iteration.
                __mmask16 cutoff_mask = _mm512_mask_cmple_ps_mask(active_mask, r2_vec, _mm512_set_1to16_ps((float)(plcutoff2)));

                if (_mm512_kortestz(cutoff_mask, cutoff_mask)) { continue; }

                #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                  __m512i pExt_j_tmp0a_vec = _mm512_load_epi32(pExt_1 + j     ); // Load  first set of 4 atoms
                  __m512i pExt_j_tmp1a_vec = _mm512_load_epi32(pExt_1 + j +  4); // Load second set of 4 atoms
                  __m512i pExt_j_tmp2a_vec = _mm512_load_epi32(pExt_1 + j +  8); // Load  third set of 4 atoms
                  __m512i pExt_j_tmp3a_vec = _mm512_load_epi32(pExt_1 + j + 12); // Load fourth set of 4 atoms
                  __m512i pExt_j_tmp0b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp0a_vec, k_2x2_0, pExt_j_tmp1a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i pExt_j_tmp1b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp1a_vec, k_2x2_1, pExt_j_tmp0a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i pExt_j_tmp2b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp2a_vec, k_2x2_0, pExt_j_tmp3a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i pExt_j_tmp3b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp3a_vec, k_2x2_1, pExt_j_tmp2a_vec, _MM_SWIZ_REG_CDAB);
                  __m512i pExt_j_tmp1c_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp1b_vec, k_4x4_0, pExt_j_tmp3b_vec, _MM_SWIZ_REG_BADC);
                  __m512i pExt_j_tmp2c_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp2b_vec, k_4x4_1, pExt_j_tmp0b_vec, _MM_SWIZ_REG_BADC);
                  __m512i pExt_j_tmp3c_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp3b_vec, k_4x4_1, pExt_j_tmp1b_vec, _MM_SWIZ_REG_BADC);
                  __m512i pExt_j_perm_pattern = _mm512_set_16to16_epi32(15, 11, 7, 3,  14, 10, 6, 2,  13, 9, 5, 1,  12, 8, 4, 0);
                  __m512i pExt_j_index_vec     = _mm512_permutevar_epi32(p_j_perm_pattern, pExt_j_tmp1c_vec);
                  __m512i pExt_j_exclIndex_vec = _mm512_permutevar_epi32(p_j_perm_pattern, pExt_j_tmp2c_vec);
                  __m512i pExt_j_maxDiff_vec   = _mm512_permutevar_epi32(p_j_perm_pattern, pExt_j_tmp3c_vec);
                #else
                  __m512i pExt_j_index_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_1_index + j);
                  __m512i pExt_j_maxDiff_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_1_exclMaxDiff + j);
                  __m512i pExt_j_exclIndex_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_1_exclIndex + j);
                #endif

                // Check the index difference vs the maxDiff value to see if there the exclusion bits
                //   for this particular pair of atoms needs to be loaded.
                __m512i indexDiff_vec = _mm512_mask_sub_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_i_index_vec, pExt_j_index_vec);
                __mmask16 indexRange_min_mask = _mm512_mask_cmpge_epi32_mask(cutoff_mask, indexDiff_vec, _mm512_sub_epi32(_mm512_setzero_epi32(), pExt_j_maxDiff_vec)); // NOTE: indexDiff >= -1 * maxDiff
                __mmask16 indexRange_max_mask = _mm512_mask_cmple_epi32_mask(cutoff_mask, indexDiff_vec, pExt_j_maxDiff_vec);
                __mmask16 indexRange_mask = _mm512_kand(indexRange_min_mask, indexRange_max_mask);

                // Calculate offsets into the exclusion flags
                __m512i offset_vec = _mm512_mask_add_epi32(_mm512_setzero_epi32(), indexRange_mask, _mm512_slli_epi32(indexDiff_vec, 1), pExt_j_exclIndex_vec);
                __m512i offset_major_vec = _mm512_srli_epi32(offset_vec, 5); // NOTE : offset / 32
                __m512i offset_minor_vec = _mm512_and_epi32(_mm512_set_1to16_epi32(0x001f), offset_vec); // NOTE : offset % 32

                // Gather exclFlags using offset_major values and then extra 2-bit fields using offset_minor.
                #if 0
                  int offset_major[16] __attribute__((aligned(64)));
                  int exclFlags[16] __attribute__((aligned(64)));
                  _mm512_store_epi32(offset_major, offset_major_vec);
                  exclFlags[ 0] = params.exclusion_bits[offset_major[ 0]];
                  exclFlags[ 1] = params.exclusion_bits[offset_major[ 1]];
                  exclFlags[ 2] = params.exclusion_bits[offset_major[ 2]];
                  exclFlags[ 3] = params.exclusion_bits[offset_major[ 3]];
                  exclFlags[ 4] = params.exclusion_bits[offset_major[ 4]];
                  exclFlags[ 5] = params.exclusion_bits[offset_major[ 5]];
                  exclFlags[ 6] = params.exclusion_bits[offset_major[ 6]];
                  exclFlags[ 7] = params.exclusion_bits[offset_major[ 7]];
                  exclFlags[ 8] = params.exclusion_bits[offset_major[ 8]];
                  exclFlags[ 9] = params.exclusion_bits[offset_major[ 9]];
                  exclFlags[10] = params.exclusion_bits[offset_major[10]];
                  exclFlags[11] = params.exclusion_bits[offset_major[11]];
                  exclFlags[12] = params.exclusion_bits[offset_major[12]];
                  exclFlags[13] = params.exclusion_bits[offset_major[13]];
                  exclFlags[14] = params.exclusion_bits[offset_major[14]];
                  exclFlags[15] = params.exclusion_bits[offset_major[15]];
                  __m512i exclFlags_vec = _mm512_load_epi32(exclFlags);
                #else
                  __m512i exclFlags_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), cutoff_mask, offset_major_vec, params.exclusion_bits, _MM_SCALE_4);
                #endif
                exclFlags_vec = _mm512_mask_srlv_epi32(_mm512_setzero_epi32(), indexRange_mask, exclFlags_vec, offset_minor_vec);
                exclFlags_vec = _mm512_and_epi32(exclFlags_vec, _mm512_set_1to16_epi32(0x03));  // NOTE : Mask out all but 2 LSBs

                // Create masks for each type of interaction (normal, modified, excluded) and
                //   store the generated pairlist entry values into the pairlists
                __mmask16 norm_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, exclFlags_vec, _mm512_setzero_epi32());
                #if MIC_CONDITION_NORMAL != 0
	          __mmask16 normhi_mask = _mm512_mask_cmplt_ps_mask(norm_mask, _mm512_set_1to16_ps(normhi_split), r2_vec);
                  norm_mask = _mm512_kxor(norm_mask, normhi_mask); // Unset any bits that were set in normhi
                #endif            
                __mmask16 excl_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, exclFlags_vec, _mm512_set_1to16_epi32(1));
                __mmask16  mod_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, exclFlags_vec, _mm512_set_1to16_epi32(2));
                _mm512_mask_packstorelo_epi32(pairlist_norm + plOffset_norm     , norm_mask, plEntry_vec);
                _mm512_mask_packstorehi_epi32(pairlist_norm + plOffset_norm + 16, norm_mask, plEntry_vec);
                #if MIC_CONDITION_NORMAL != 0
                  _mm512_mask_packstorelo_epi32(pairlist_normhi + plOffset_normhi     , normhi_mask, plEntry_vec);
                  _mm512_mask_packstorehi_epi32(pairlist_normhi + plOffset_normhi + 16, normhi_mask, plEntry_vec);
                #endif
                _mm512_mask_packstorelo_epi32(pairlist_excl + plOffset_excl     , excl_mask, plEntry_vec);
                _mm512_mask_packstorehi_epi32(pairlist_excl + plOffset_excl + 16, excl_mask, plEntry_vec);
                _mm512_mask_packstorelo_epi32(pairlist_mod  + plOffset_mod      ,  mod_mask, plEntry_vec);
                _mm512_mask_packstorehi_epi32(pairlist_mod  + plOffset_mod  + 16,  mod_mask, plEntry_vec);
                __m512i one_vec = _mm512_set_1to16_epi32(1);
                plOffset_norm += _mm512_mask_reduce_add_epi32(norm_mask, one_vec);
                #if MIC_CONDITION_NORMAL != 0
                  plOffset_normhi += _mm512_mask_reduce_add_epi32(normhi_mask, one_vec);
                #endif
                plOffset_excl += _mm512_mask_reduce_add_epi32(excl_mask, one_vec);
                plOffset_mod  += _mm512_mask_reduce_add_epi32( mod_mask, one_vec);
              }

            #else

              // The "j" loop, iterating over the second list of atoms
              for (int j = j_lo; j < j_hi; j++) {

                // Load the "j" atom's position, calculate/check the distance (squared)
                //   between the "i" and "j" atoms
                #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                  CALC_TYPE p_d_x = p_i_x - p_1[j].x;
                  CALC_TYPE p_d_y = p_i_y - p_1[j].y;
                  CALC_TYPE p_d_z = p_i_z - p_1[j].z;
                #else
                  CALC_TYPE p_d_x = p_i_x - p_1_x[j];
                  CALC_TYPE p_d_y = p_i_y - p_1_y[j];
                  CALC_TYPE p_d_z = p_i_z - p_1_z[j];
                #endif
                CALC_TYPE r2 = (p_d_x * p_d_x) + (p_d_y * p_d_y) + (p_d_z * p_d_z);
                if (r2 <= plcutoff2) {

                  // Check the exclusion bits to set the modified and excluded flags
                  // NOTE: This code MUST match the exclusion lists generation code
                  //   on ComputeNonbondedMIC.C, which generates the flags.  In
                  //   particular, the order of the 2 bits must be correct, per the
                  //   pairlist format listed below (i.e. isModified -> MSB).
                  int exclFlags = 0x00;
                  #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                    const int indexDiff = pExt_0[i].index - pExt_1[i].index;
                    const int maxDiff = pExt_1[j].excl_maxdiff;
                  #else
                    const int indexDiff = pExt_0_index[i] - pExt_1_index[j];
                    const int maxDiff = pExt_1_exclMaxDiff[j];
                  #endif
                  if (indexDiff >= -1 * maxDiff && indexDiff <= maxDiff) {
                    #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                      const int offset = (2 * indexDiff) + pExt_1[j].excl_index;
                    #else
                      const int offset = (2 * indexDiff) + pExt_1_exclIndex[j];
                    #endif
                    const int offset_major = offset / (sizeof(unsigned int) * 8);
                    const int offset_minor = offset % (sizeof(unsigned int) * 8); // NOTE: Reverse indexing direction relative to offset_major
                    exclFlags = ((params.exclusion_bits[offset_major]) >> offset_minor) & 0x03;
                  }

                  // Create the pairlist entry value and store it
                  const int plEntry = ((i & 0xFFFF) << 16) | (j & 0xFFFF);
                  if (exclFlags == 0) {
                    #if MIC_CONDITION_NORMAL != 0
                      if (r2 <= normhi_split) {
                        pairlist_norm[plOffset_norm++] = plEntry;
                      } else {
                        pairlist_normhi[plOffset_normhi++] = plEntry;
                      }
                    #else
                      pairlist_norm[plOffset_norm++] = plEntry;
                    #endif
                  } else if (exclFlags == 1) {
                    pairlist_excl[plOffset_excl++] = plEntry;
                  } else if (exclFlags == 2) {
                    pairlist_mod[plOffset_mod++] = plEntry;
                  }

                } // end if (r2 <= plcutoff2
	      } // end for (j < j_hi)

            #endif

            #if __MIC_PAD_PLGEN_CTRL != 0
              #if MIC_HANDCODE_FORCE_SINGLE != 0
                const int padLen = 16;
              #else
                const int padLen = 8;
              #endif
              const int padValue = (i << 16) | 0xFFFF;
              // NOTE: xxx % 8 != 2 because pairlist offset includes first two ints (pairlist size and alloc size)
              while (plOffset_norm % padLen != 2) { pairlist_norm[plOffset_norm++] = padValue; }
              while (plOffset_mod  % padLen != 2) { pairlist_mod [plOffset_mod++ ] = padValue; }
              while (plOffset_excl % padLen != 2) { pairlist_excl[plOffset_excl++] = padValue; }
              #if MIC_CONDITION_NORMAL != 0
                while (plOffset_normhi % padLen != 2) { pairlist_normhi[plOffset_normhi++] = padValue; }
              #endif 
           #endif

          } // end for (i < i_hi)
        } // end (_j < j_upper; _j += plTileSize)
      } // end (_i < i_upper; _i += plTileSize)

    #else  // MIC_TILE_PLGEN

      // Do the distance checks for all possible pairs of atoms between the
      //   two input atom arrays
      // The "i" loop, iterating over the first list of atoms

      int numParts = params.pp->numParts;
      int part = params.pp->part;

      // Generate plGenILo and plGenHi values (the portion of the "i" atoms that will be processed by this "part").
      //   Note that this is done differently for self computes, because of the unequal work per "i" atom.  For
      //   pair computes, the "i" iteration space is just divided up evenly.
      #if 0
      int plGenILo = 0;
      int plGenIHi = i_upper SELF(-1);
      if (numParts > 1) {

        #if (0 SELF(+1))

          // For self computes where the iteration space forms a triangle, divide the rows in the
          //   iteration space so more "shorter" rows are grouped together while fewer "longer" rows
          //   are grouped together.
          float totalArea = ((i_upper) * (i_upper - 1)) / 2.0f;
          float areaPerPart = totalArea / numParts;

          // NOTE : We want to divide the triangular iteration space (0 <= i < i_upper - 1, i < j < i_upper).
          //   Since we know the area per part and the area of the iteration space "before" some arbitrary
          //   index X (i.e. 0 <= i < X), we can calculate indexes for each part.
          //     A = X(X-1)/2 where A = area "before" X
          //     solving for X we get X^2 - X - 2A = 0, or via the quadradic equation, X = (1 +- sqrt(1+8A))/2
          // NOTE: This calculation might be a little messy, so double check the border cases
          plGenILo = ((part==0)?(0) : ((int)((1.0f+sqrtf(1+8*(areaPerPart*part)))/2.0f)) );
          plGenIHi = ((part==numParts-1)?(i_upper-1) : ((int)((1.0f+sqrtf(1+8*(areaPerPart*(part+1))))/2.0f)) );
          // Reverse the indexes since this calculation assumes i=0 has little work i=i_upper-1 has lots of
          //   work (i.e. upper triangular versus lower triangular)
          int plGenTmp = plGenILo;
          plGenILo = (i_upper - 1) - plGenIHi;
          plGenIHi = (i_upper - 1) - plGenTmp;

        #else  // SELF

          // For pair computes where the iteration space forms a square, divide the rows in the
          //   iteration space evenly since they all have the same area
          plGenILo = (int)(((float)(i_upper)) * ((float)(part    )) / ((float)(numParts)));
          plGenIHi = (int)(((float)(i_upper)) * ((float)(part + 1)) / ((float)(numParts)));

        #endif  // SELF

        // Round up plGenILo and plGenIHi to a multiple of 16 so there is no false sharing
        //   on force output cachelines
        plGenILo = (plGenILo + 15) & (~15);
        plGenIHi = (plGenIHi + 15) & (~15);
        if (plGenIHi > i_upper SELF(-1)) { plGenIHi = i_upper SELF(-1); }
      }

      #pragma loop_count (300)
      for (int i = plGenILo; i < plGenIHi; i++) {

      #else

      #pragma loop_count (300)
      for (int i = part; i < i_upper SELF(-1); i += numParts) {

      #endif

        // If the pairlist is not long enough, grow the pairlist
        PAIRLIST_GROW_CHECK(pairlist_norm, plOffset_norm, 50);
        PAIRLIST_GROW_CHECK(pairlist_mod , plOffset_mod,   8);
        PAIRLIST_GROW_CHECK(pairlist_excl, plOffset_excl,  8);
        #if MIC_CONDITION_NORMAL != 0
          PAIRLIST_GROW_CHECK(pairlist_normhi, plOffset_normhi, 15);
        #endif

        // Load position information for the current "i" atom
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          CALC_TYPE p_i_x = p_0[i].x + offset_x;
          CALC_TYPE p_i_y = p_0[i].y + offset_y;
          CALC_TYPE p_i_z = p_0[i].z + offset_z;
        #else
          CALC_TYPE p_i_x = p_0_x[i] + offset_x;
          CALC_TYPE p_i_y = p_0_y[i] + offset_y;
          CALC_TYPE p_i_z = p_0_z[i] + offset_z;
        #endif

        #if MIC_HANDCODE_PLGEN != 0

	  // Setup loop constants ("i" atom values) and iterator vector.
          __m512 p_i_x_vec = _mm512_set_1to16_ps(p_i_x);
          __m512 p_i_y_vec = _mm512_set_1to16_ps(p_i_y);
          __m512 p_i_z_vec = _mm512_set_1to16_ps(p_i_z);
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            __m512i pExt_i_index_vec = _mm512_set_1to16_epi32(pExt_0[i].index);
          #else
            __m512i pExt_i_index_vec = _mm512_set_1to16_epi32(pExt_0_index[i]);
          #endif

          __m512i i_shifted_vec = _mm512_slli_epi32(_mm512_set_1to16_epi32(i), 16); 
          __m512i j_vec = _mm512_set_16to16_epi32(15, 14, 13, 12, 11, 10,  9,  8,
                                                   7,  6,  5,  4,  3,  2,  1,  0);

          // For self computes, round (i+1) down to nearest multiple of 16 and then
          //   add that value to j_vec (skipping iterations where all j values
          //   are < i+1)
          #if (0 SELF(+1))
            const int ip1_16 = (i+1) & (~0x0F);
            j_vec = _mm512_add_epi32(j_vec, _mm512_set_1to16_epi32(ip1_16));
          #endif

          // The "j" loop, iterating over the second list of atoms
          #pragma novector
          #pragma loop count(35)
          for (int j = 0 SELF(+ ip1_16); j < j_upper; j += 16) {

            // Create the active_mask
            __mmask16 active_mask = _mm512_cmplt_epi32_mask(j_vec, _mm512_set_1to16_epi32(j_upper));
            #if (0 SELF(+1))
              active_mask = _mm512_kand(active_mask, _mm512_cmpgt_epi32_mask(j_vec, _mm512_set_1to16_epi32(i)));
            #endif

            // Create the pairlist entry values for these would-be interactions before advancing
            //   the j_vec values (i in upper 16 bits and j in lower 16 bits)
            __m512i plEntry_vec = _mm512_or_epi32(i_shifted_vec, j_vec);

            // Advance j_vec counter
            j_vec = _mm512_add_epi32(j_vec, _mm512_set_1to16_epi32(16));

            // Load the positions for the "j" atoms
            #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
              __m512i p_j_tmp0a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j     )); // Load  first set of 4 atoms
              __m512i p_j_tmp1a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j +  4)); // Load second set of 4 atoms
              __m512i p_j_tmp2a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j +  8)); // Load  third set of 4 atoms
              __m512i p_j_tmp3a_vec = _mm512_castps_si512(_mm512_load_ps(p_1 + j + 12)); // Load fourth set of 4 atoms
              __mmask16 k_2x2_0 = _mm512_int2mask(0xAAAA);
              __mmask16 k_2x2_1 = _mm512_int2mask(0x5555);
              __m512i p_j_tmp0b_vec = _mm512_mask_swizzle_epi32(p_j_tmp0a_vec, k_2x2_0, p_j_tmp1a_vec, _MM_SWIZ_REG_CDAB);
              __m512i p_j_tmp1b_vec = _mm512_mask_swizzle_epi32(p_j_tmp1a_vec, k_2x2_1, p_j_tmp0a_vec, _MM_SWIZ_REG_CDAB);
              __m512i p_j_tmp2b_vec = _mm512_mask_swizzle_epi32(p_j_tmp2a_vec, k_2x2_0, p_j_tmp3a_vec, _MM_SWIZ_REG_CDAB);
              __m512i p_j_tmp3b_vec = _mm512_mask_swizzle_epi32(p_j_tmp3a_vec, k_2x2_1, p_j_tmp2a_vec, _MM_SWIZ_REG_CDAB);
              __mmask16 k_4x4_0 = _mm512_int2mask(0xCCCC);
              __mmask16 k_4x4_1 = _mm512_int2mask(0x3333);
              __m512i p_j_tmp0c_vec = _mm512_mask_swizzle_epi32(p_j_tmp0b_vec, k_4x4_0, p_j_tmp2b_vec, _MM_SWIZ_REG_BADC);
              __m512i p_j_tmp1c_vec = _mm512_mask_swizzle_epi32(p_j_tmp1b_vec, k_4x4_0, p_j_tmp3b_vec, _MM_SWIZ_REG_BADC);
              __m512i p_j_tmp2c_vec = _mm512_mask_swizzle_epi32(p_j_tmp2b_vec, k_4x4_1, p_j_tmp0b_vec, _MM_SWIZ_REG_BADC);
              __m512i p_j_perm_pattern = _mm512_set_16to16_epi32(15, 11, 7, 3,  14, 10, 6, 2,  13, 9, 5, 1,  12, 8, 4, 0);
              __m512 p_j_x_vec = _mm512_castsi512_ps(_mm512_permutevar_epi32(p_j_perm_pattern, p_j_tmp0c_vec));
              __m512 p_j_y_vec = _mm512_castsi512_ps(_mm512_permutevar_epi32(p_j_perm_pattern, p_j_tmp1c_vec));
              __m512 p_j_z_vec = _mm512_castsi512_ps(_mm512_permutevar_epi32(p_j_perm_pattern, p_j_tmp2c_vec));
            #else
              __m512 p_j_x_vec = _mm512_mask_load_ps(_mm512_setzero_ps(), active_mask, p_1_x + j);
              __m512 p_j_y_vec = _mm512_mask_load_ps(_mm512_setzero_ps(), active_mask, p_1_y + j);
              __m512 p_j_z_vec = _mm512_mask_load_ps(_mm512_setzero_ps(), active_mask, p_1_z + j);
            #endif

            // Calculate the distance between "i" atom and these "j" atoms (in double)
            __m512 p_ij_x_vec = _mm512_sub_ps(p_i_x_vec, p_j_x_vec);
            __m512 p_ij_y_vec = _mm512_sub_ps(p_i_y_vec, p_j_y_vec);
            __m512 p_ij_z_vec = _mm512_sub_ps(p_i_z_vec, p_j_z_vec);
            __m512 r2_vec = _mm512_mul_ps(p_ij_x_vec, p_ij_x_vec);
            r2_vec = _mm512_add_ps(r2_vec, _mm512_mul_ps(p_ij_y_vec, p_ij_y_vec));
            r2_vec = _mm512_add_ps(r2_vec, _mm512_mul_ps(p_ij_z_vec, p_ij_z_vec));

            // Do cutoff distance check.  If there are no particles within cutoff, move on
            //   to the next iteration.
            __mmask16 cutoff_mask = _mm512_mask_cmple_ps_mask(active_mask, r2_vec, _mm512_set_1to16_ps((float)(plcutoff2)));

            // If nothing passed the cutoff check, then move on to the next set of atoms (iteration)
            if (_mm512_kortestz(cutoff_mask, cutoff_mask)) { continue; }

            #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
              __m512i pExt_j_tmp0a_vec = _mm512_load_epi32(pExt_1 + j     ); // Load  first set of 4 atoms
              __m512i pExt_j_tmp1a_vec = _mm512_load_epi32(pExt_1 + j +  4); // Load second set of 4 atoms
              __m512i pExt_j_tmp2a_vec = _mm512_load_epi32(pExt_1 + j +  8); // Load  third set of 4 atoms
              __m512i pExt_j_tmp3a_vec = _mm512_load_epi32(pExt_1 + j + 12); // Load fourth set of 4 atoms
              __m512i pExt_j_tmp0b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp0a_vec, k_2x2_0, pExt_j_tmp1a_vec, _MM_SWIZ_REG_CDAB);
              __m512i pExt_j_tmp1b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp1a_vec, k_2x2_1, pExt_j_tmp0a_vec, _MM_SWIZ_REG_CDAB);
              __m512i pExt_j_tmp2b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp2a_vec, k_2x2_0, pExt_j_tmp3a_vec, _MM_SWIZ_REG_CDAB);
              __m512i pExt_j_tmp3b_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp3a_vec, k_2x2_1, pExt_j_tmp2a_vec, _MM_SWIZ_REG_CDAB);
              __m512i pExt_j_tmp1c_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp1b_vec, k_4x4_0, pExt_j_tmp3b_vec, _MM_SWIZ_REG_BADC);
              __m512i pExt_j_tmp2c_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp2b_vec, k_4x4_1, pExt_j_tmp0b_vec, _MM_SWIZ_REG_BADC);
              __m512i pExt_j_tmp3c_vec = _mm512_mask_swizzle_epi32(pExt_j_tmp3b_vec, k_4x4_1, pExt_j_tmp1b_vec, _MM_SWIZ_REG_BADC);
              __m512i pExt_j_perm_pattern = _mm512_set_16to16_epi32(15, 11, 7, 3,  14, 10, 6, 2,  13, 9, 5, 1,  12, 8, 4, 0);
              __m512i pExt_j_index_vec     = _mm512_permutevar_epi32(p_j_perm_pattern, pExt_j_tmp1c_vec);
              __m512i pExt_j_exclIndex_vec = _mm512_permutevar_epi32(p_j_perm_pattern, pExt_j_tmp2c_vec);
              __m512i pExt_j_maxDiff_vec   = _mm512_permutevar_epi32(p_j_perm_pattern, pExt_j_tmp3c_vec);
            #else
              __m512i pExt_j_index_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_1_index + j);
              __m512i pExt_j_maxDiff_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_1_exclMaxDiff + j);
              __m512i pExt_j_exclIndex_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_1_exclIndex + j);
            #endif

            // Check the index difference vs the maxDiff value to see if there the exclusion bits
            //   for this particular pair of atoms needs to be loaded.
            __m512i indexDiff_vec = _mm512_mask_sub_epi32(_mm512_setzero_epi32(), cutoff_mask, pExt_i_index_vec, pExt_j_index_vec);
            __mmask16 indexRange_min_mask = _mm512_mask_cmpge_epi32_mask(cutoff_mask, indexDiff_vec, _mm512_sub_epi32(_mm512_setzero_epi32(), pExt_j_maxDiff_vec)); // NOTE: indexDiff >= -1 * maxDiff
            __mmask16 indexRange_max_mask = _mm512_mask_cmple_epi32_mask(cutoff_mask, indexDiff_vec, pExt_j_maxDiff_vec);
            __mmask16 indexRange_mask = _mm512_kand(indexRange_min_mask, indexRange_max_mask);

            // Calculate offsets into the exclusion flags
            __m512i offset_vec = _mm512_mask_add_epi32(_mm512_setzero_epi32(), indexRange_mask, _mm512_slli_epi32(indexDiff_vec, 1), pExt_j_exclIndex_vec);
            __m512i offset_major_vec = _mm512_srli_epi32(offset_vec, 5); // NOTE : offset / 32
            __m512i offset_minor_vec = _mm512_and_epi32(_mm512_set_1to16_epi32(0x001f), offset_vec); // NOTE : offset % 32

            // Gather exclFlags using offset_major values and then extra 2-bit fields using offset_minor.
            __m512i exclFlags_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), cutoff_mask, offset_major_vec, params.exclusion_bits, _MM_SCALE_4);
            exclFlags_vec = _mm512_mask_srlv_epi32(_mm512_setzero_epi32(), indexRange_mask, exclFlags_vec, offset_minor_vec);
            exclFlags_vec = _mm512_and_epi32(exclFlags_vec, _mm512_set_1to16_epi32(0x03));  // NOTE : Mask out all but 2 LSBs

            // Create masks for each type of interaction (normal, modified, excluded) and
            //   store the generated pairlist entry values into the pairlists
            __mmask16 norm_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, exclFlags_vec, _mm512_setzero_epi32());
            #if MIC_CONDITION_NORMAL != 0
	      __mmask16 normhi_mask = _mm512_mask_cmplt_ps_mask(norm_mask, _mm512_set_1to16_ps(normhi_split), r2_vec);
              norm_mask = _mm512_kxor(norm_mask, normhi_mask); // Unset any bits that were set in normhi
            #endif            
            __mmask16 excl_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, exclFlags_vec, _mm512_set_1to16_epi32(1));
            __mmask16  mod_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, exclFlags_vec, _mm512_set_1to16_epi32(2));
            _mm512_mask_packstorelo_epi32(pairlist_norm + plOffset_norm     , norm_mask, plEntry_vec);
            _mm512_mask_packstorehi_epi32(pairlist_norm + plOffset_norm + 16, norm_mask, plEntry_vec);
            #if MIC_CONDITION_NORMAL != 0
              _mm512_mask_packstorelo_epi32(pairlist_normhi + plOffset_normhi     , normhi_mask, plEntry_vec);
              _mm512_mask_packstorehi_epi32(pairlist_normhi + plOffset_normhi + 16, normhi_mask, plEntry_vec);
            #endif
            _mm512_mask_packstorelo_epi32(pairlist_excl + plOffset_excl     , excl_mask, plEntry_vec);
            _mm512_mask_packstorehi_epi32(pairlist_excl + plOffset_excl + 16, excl_mask, plEntry_vec);
            _mm512_mask_packstorelo_epi32(pairlist_mod  + plOffset_mod      ,  mod_mask, plEntry_vec);
            _mm512_mask_packstorehi_epi32(pairlist_mod  + plOffset_mod  + 16,  mod_mask, plEntry_vec);
            __m512i one_vec = _mm512_set_1to16_epi32(1);

            // Move the offsets forward by the number of atoms added to each list
            plOffset_norm += _mm512_mask_reduce_add_epi32(norm_mask, one_vec);
            #if MIC_CONDITION_NORMAL != 0
              plOffset_normhi += _mm512_mask_reduce_add_epi32(normhi_mask, one_vec);
            #endif
            plOffset_excl += _mm512_mask_reduce_add_epi32(excl_mask, one_vec);
            plOffset_mod  += _mm512_mask_reduce_add_epi32( mod_mask, one_vec);
          }

        #else // if MIC_HANDCODE_PLGEN != 0

          // The "j" loop, iterating over the second list of atoms
          #if (0 PAIR(+1))
            #pragma loop_count (30)
          #elif (0 SELF(+1))
            #pragma loop_count (300)
          #endif
          for (int j = 0 SELF(+i+1); j < j_upper; j++) {

            // Load the "j" atom's position, calculate/check the distance (squared)
            //   between the "i" and "j" atoms
            #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
              CALC_TYPE p_d_x = p_i_x - p_1[j].x;
              CALC_TYPE p_d_y = p_i_y - p_1[j].y;
              CALC_TYPE p_d_z = p_i_z - p_1[j].z;
            #else
              CALC_TYPE p_d_x = p_i_x - p_1_x[j];
              CALC_TYPE p_d_y = p_i_y - p_1_y[j];
              CALC_TYPE p_d_z = p_i_z - p_1_z[j];
            #endif
            CALC_TYPE r2 = (p_d_x * p_d_x) + (p_d_y * p_d_y) + (p_d_z * p_d_z);
            if (r2 <= plcutoff2) {

              // Check the exclusion bits to set the modified and excluded flags
              // NOTE: This code MUST match the exclusion lists generation code
              //   on ComputeNonbondedMIC.C, which generates the flags.  In
              //   particular, the order of the 2 bits must be correct, per the
              //   pairlist format listed below (i.e. isModified -> MSB).
              int exclFlags = 0x00;
              #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                const int indexDiff = pExt_0[i].index - pExt_1[j].index;
                const int maxDiff = pExt_1[j].excl_maxdiff;
              #else
                const int indexDiff = pExt_0_index[i] - pExt_1_index[j];
                const int maxDiff = pExt_1_exclMaxDiff[j];
              #endif
              if (indexDiff >= -1 * maxDiff && indexDiff <= maxDiff) {
                #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                  const int offset = (2 * indexDiff) + pExt_1[j].excl_index;
                #else
                  const int offset = (2 * indexDiff) + pExt_1_exclIndex[j];
                #endif
                const int offset_major = offset / (sizeof(unsigned int) * 8);
                const int offset_minor = offset % (sizeof(unsigned int) * 8); // NOTE: Reverse indexing direction relative to offset_major
                exclFlags = ((params.exclusion_bits[offset_major]) >> offset_minor) & 0x03;
	      }

              // Create the pairlist entry value (i and j value) and store it in to the
              //   appropriate pairlist based on the type of interaction (exclFlags value)
              const int plEntry = ((i & 0xFFFF) << 16) | (j & 0xFFFF);
              if (exclFlags == 0) {
                #if MIC_CONDITION_NORMAL != 0
                  if (r2 <= normhi_split) {
                    pairlist_norm[plOffset_norm++] = plEntry;
                  } else {
                    pairlist_normhi[plOffset_normhi++] = plEntry;
                  }
                #else
                  pairlist_norm[plOffset_norm++] = plEntry;
                #endif
              } else if (exclFlags == 1) {
                pairlist_excl[plOffset_excl++] = plEntry;
              } else if (exclFlags == 2) {
                pairlist_mod[plOffset_mod++] = plEntry;
              }

            } // end if (r2 < plcutoff2)
          } // end for (j < j_upper)

        #endif // if MIC_HANDCODE_PLGEN != 0

        //#if MIC_PAD_PLGEN != 0
        #if __MIC_PAD_PLGEN_CTRL != 0
          //#if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
          #if MIC_HANDCODE_FORCE_SINGLE != 0
            const int padLen = 16;
          #else
            const int padLen = 8;
          #endif
          const int padValue = (i << 16) | 0xFFFF;
          // NOTE: xxx % padLen != 2 because offset includes first two ints (pairlist sizes), so 0 entries <-> offset 2, 16 entries <-> offset 18, etc.
          //   Alignment is setup so that the entries (minus the first two ints) are cacheline aligned.
          while (plOffset_norm % padLen != 2) { pairlist_norm[plOffset_norm++] = padValue; }
          while (plOffset_mod  % padLen != 2) { pairlist_mod [plOffset_mod++ ] = padValue; }
          while (plOffset_excl % padLen != 2) { pairlist_excl[plOffset_excl++] = padValue; }
          #if MIC_CONDITION_NORMAL != 0
            while (plOffset_normhi % padLen != 2) { pairlist_normhi[plOffset_normhi++] = padValue; }
          #endif
        #endif

      } // end for (i < i_upper)

    #endif  // if MIC_TILE_PLGEN != 0

    // If we are conditioning the normal pairlist, place the contents of pairlist_normhi
    //   at the end of pairlist_norm
    #if MIC_CONDITION_NORMAL != 0

      // Make sure there is enough room in pairlist_norm for the contents of pairlist_normhi
      if (plOffset_norm + plOffset_normhi > pairlist_norm[0]) {
        int newPlAllocSize = pairlist_norm[0] + 2 * plOffset_normhi + 8;
        int* RESTRICT newPl = (int*)_MM_MALLOC_WRAPPER(newPlAllocSize * sizeof(int) + 64, 64, "pairlist_norm grow for conditioning");
        __ASSERT(newPl != NULL);
        newPl += 14; // Align pairlist entries, which start at index 2
        memcpy(newPl, pairlist_norm, plOffset_norm * sizeof(int));
        _MM_FREE_WRAPPER(pairlist_norm - 14);
        pairlist_norm = newPl;
        pairlist_norm[0] = newPlAllocSize;
      }

      // Copy the contents of pairlist_normhi into pairlist_norm
      for (int i = 0; i < plOffset_normhi - 2; i++) {
        pairlist_norm[plOffset_norm + i] = pairlist_normhi[2 + i];
      }
      plOffset_norm += plOffset_normhi - 2;
      plOffset_normhi = 2;

    #endif

    // Store the current offsets (pairlist lengths) in to the pairlist data structure itself
    pairlist_norm[1] = plOffset_norm; // NOTE: The size includes the initial 2 ints
    pairlist_mod[1] = plOffset_mod;
    pairlist_excl[1] = plOffset_excl;
    #if MIC_CONDITION_NORMAL != 0
      pairlist_normhi[1] = plOffset_normhi;
    #endif

    // DMK - DEBUG - Pairlist memory size info
    #if MIC_TRACK_DEVICE_MEM_USAGE != 0
      if (plOffset_norm > pairlist_norm[-1]) { pairlist_norm[-1] = plOffset_norm; } // NOTE: Because of alignment of pairlist payload, known to have 14 ints allocated prior to pairlist start
      if (plOffset_mod > pairlist_mod[-1]) { pairlist_mod[-1] = plOffset_mod; } // NOTE: Because of alignment of pairlist payload, known to have 14 ints allocated prior to pairlist start
      if (plOffset_excl > pairlist_excl[-1]) { pairlist_excl[-1] = plOffset_excl; } // NOTE: Because of alignment of pairlist payload, known to have 14 ints allocated prior to pairlist start
      #if MIC_CONDITION_NORMAL != 0
        if (plOffset_normhi > pairlist_normhi[-1]) { pairlist_normhi[-1] = plOffset_normhi; } // NOTE: Because of alignment of pairlist payload, known to have 14 ints allocated prior to pairlist start
      #endif
    #endif

    // If prefetch is enabled, pad-out some extra (valid) entries in the pairlist
    //   that can be safely prefetched (i.e. without segfaulting), even though they will have
    //   no actual effect (i.e. forces generated related to them)
    #if MIC_PREFETCH_DISTANCE > 0
      for (int pfI = 0; pfI < MIC_PREFETCH_DISTANCE; pfI++) {
        pairlist_norm[plOffset_norm + pfI] = 0;
        pairlist_mod [plOffset_mod  + pfI] = 0;
        pairlist_excl[plOffset_excl + pfI] = 0;
      }
    #endif

    // If we need to pad the pairlists with valid entries, create "zero" entries at the end
    #if MIC_HANDCODE_FORCE != 0 && MIC_HANDCODE_FORCE_PFDIST > 0
      for (int pfI = 0; pfI < 16 * MIC_HANDCODE_FORCE_PFDIST; pfI++) {
        #if MIC_HANDCODE_FORCE_SINGLE != 0
          pairlist_norm[plOffset_norm + pfI] = -1;
          pairlist_mod [plOffset_mod  + pfI] = -1;
          pairlist_excl[plOffset_excl + pfI] = -1;
        #else
          pairlist_norm[plOffset_norm + pfI] = 0;
          pairlist_mod [plOffset_mod  + pfI] = 0;
          pairlist_excl[plOffset_excl + pfI] = 0;
        #endif
      }
    #endif

    // Save the pointer values back into the pairlist pointer array, so the buffers can
    //   be reused across timesteps, in case the pointers were changed here
    params.pairlists_ptr[PL_NORM_INDEX] = pairlist_norm;
    params.pairlists_ptr[ PL_MOD_INDEX] = pairlist_mod;
    params.pairlists_ptr[PL_EXCL_INDEX] = pairlist_excl;
    #if MIC_CONDITION_NORMAL != 0
      params.pairlists_ptr[PL_NORMHI_INDEX] = pairlist_normhi;
    #endif

  } // end if (params.savePairlists)

  #undef PAIRLIST_ALLOC_CHECK
  #undef PAIRLIST_GROW_CHECK

  // Declare the scalars (or vector-width-long arrays if using force splitting) related
  //   to accumulating virial and energy output, initializing these values to zero
  #if (0 FAST(+1))
    #if (0 ENERGY(+1))
      double vdwEnergy = 0;
      SHORT( double electEnergy = 0; )
    #endif
    #if (0 SHORT(+1))
      double virial_xx = 0;
      double virial_xy = 0;
      double virial_xz = 0;
      double virial_yy = 0;
      double virial_yz = 0;
      double virial_zz = 0;
    #endif
  #endif
  #if (0 FULL(+1))
    #if (0 ENERGY(+1))
      double fullElectEnergy = 0;
    #endif
    double fullElectVirial_xx = 0;
    double fullElectVirial_xy = 0;
    double fullElectVirial_xz = 0;
    double fullElectVirial_yy = 0;
    double fullElectVirial_yz = 0;
    double fullElectVirial_zz = 0;
  #endif

  // NORMAL LOOP
  DEVICE_FPRINTF("N");
  #define PAIRLIST pairlist_norm
  #define NORMAL(X) X
  #define EXCLUDED(X)
  #define MODIFIED(X)
    #include "ComputeNonbondedMICKernelBase2.h"
  #undef PAIRLIST
  #undef NORMAL
  #undef EXCLUDED
  #undef MODIFIED

  // MODIFIED LOOP
  DEVICE_FPRINTF("M");
  #define PAIRLIST pairlist_mod
  #define NORMAL(X)
  #define EXCLUDED(X)
  #define MODIFIED(X) X
    #include "ComputeNonbondedMICKernelBase2.h"
  #undef PAIRLIST
  #undef NORMAL
  #undef EXCLUDED
  #undef MODIFIED

  // EXCLUDED LOOP
  #if defined(FULLELECT)

    DEVICE_FPRINTF("E");
    #define PAIRLIST pairlist_excl
    #undef FAST
    #define FAST(X)
    #define NORMAL(X)
    #define EXCLUDED(X) X
    #define MODIFIED(X)
      #include "ComputeNonbondedMICKernelBase2.h"
    #undef PAIRLIST
    #undef FAST
    #ifdef SLOWONLY
      #define FAST(X)
    #else
      #define FAST(X) X
    #endif
    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED

  #else  // defined(FULLELECT)

    // If we are not executing the excluded loop, then we need to count exclusive
    //   interactions per atom.  Contribute the number of entries in the EXCLUDED
    //   pairlist (valid only if padding pairlists)
    DEVICE_FPRINTF("e");
    #if MIC_EXCL_CHECKSUM_FULL != 0
    {
      #if __MIC_PAD_PLGEN_CTRL != 0
        // NOTE: If using padding, the pairlist will have 'invalid' entries mixed in with the
        //   valid entries, so we need to scan through and only count the valid entries
        const int * const pairlist_excl_base = pairlist_excl + 2;
        const int pairlist_excl_len = pairlist_excl[1] - 2;
        int exclSum = 0;
        __ASSUME_ALIGNED(pairlist_excl_base);
        #pragma simd reduction(+ : exclSum)
        for (int plI = 0; plI < pairlist_excl_len; plI++) {
          if ((pairlist_excl_base[plI] & 0xFFFF) != 0xFFFF) {
            exclSum += 1;
	  }
        }
        params.exclusionSum += exclSum;
      #else
        params.exclusionSum += pairlist_excl[1] - 2; // NOTE: Size includes first two elements (alloc size and size), don't count those
      #endif
    }
    #endif

  #endif  // defined(FULLELECT)

  // If this is a self compute, do the virial calculation here
  #if (0 SHORT( FAST( SELF(+1))))
    for (int i = 0; i < i_upper; i++) {
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        virial_xx += f_0[i].x * (p_0[i].x + params.patch1_center_x);
        virial_xy += f_0[i].x * (p_0[i].y + params.patch1_center_y);
        virial_xz += f_0[i].x * (p_0[i].z + params.patch1_center_z);
        virial_yy += f_0[i].y * (p_0[i].y + params.patch1_center_y);
        virial_yz += f_0[i].y * (p_0[i].z + params.patch1_center_z);
        virial_zz += f_0[i].z * (p_0[i].z + params.patch1_center_z);
      #else
        virial_xx += f_0_x[i] * (p_0_x[i] + params.patch1_center_x);
        virial_xy += f_0_x[i] * (p_0_y[i] + params.patch1_center_y);
        virial_xz += f_0_x[i] * (p_0_z[i] + params.patch1_center_z);
        virial_yy += f_0_y[i] * (p_0_y[i] + params.patch1_center_y);
        virial_yz += f_0_y[i] * (p_0_z[i] + params.patch1_center_z);
        virial_zz += f_0_z[i] * (p_0_z[i] + params.patch1_center_z);
      #endif
    }
  #endif
  #if (0 FULL( SELF(+1)))
    for (int i = 0; i < i_upper; i++) {
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        fullElectVirial_xx += fullf_0[i].x * (p_0[i].x + params.patch1_center_x);
        fullElectVirial_xy += fullf_0[i].x * (p_0[i].y + params.patch1_center_y);
        fullElectVirial_xz += fullf_0[i].x * (p_0[i].z + params.patch1_center_z);
        fullElectVirial_yy += fullf_0[i].y * (p_0[i].y + params.patch1_center_y);
        fullElectVirial_yz += fullf_0[i].y * (p_0[i].z + params.patch1_center_z);
        fullElectVirial_zz += fullf_0[i].z * (p_0[i].z + params.patch1_center_z);
      #else
        fullElectVirial_xx += fullf_0_x[i] * (p_0_x[i] + params.patch1_center_x);
        fullElectVirial_xy += fullf_0_x[i] * (p_0_y[i] + params.patch1_center_y);
        fullElectVirial_xz += fullf_0_x[i] * (p_0_z[i] + params.patch1_center_z);
        fullElectVirial_yy += fullf_0_y[i] * (p_0_y[i] + params.patch1_center_y);
        fullElectVirial_yz += fullf_0_y[i] * (p_0_z[i] + params.patch1_center_z);
        fullElectVirial_zz += fullf_0_z[i] * (p_0_z[i] + params.patch1_center_z);
      #endif
    }
  #endif

  FAST( SHORT( params.virial_xx = virial_xx; ) )
  FAST( SHORT( params.virial_xy = virial_xy; ) )
  FAST( SHORT( params.virial_xz = virial_xz; ) )
  FAST( SHORT( params.virial_yy = virial_yy; ) )
  FAST( SHORT( params.virial_yz = virial_yz; ) )
  FAST( SHORT( params.virial_zz = virial_zz; ) )
  FULL( params.fullElectVirial_xx = fullElectVirial_xx; )
  FULL( params.fullElectVirial_xy = fullElectVirial_xy; )
  FULL( params.fullElectVirial_xz = fullElectVirial_xz; )
  FULL( params.fullElectVirial_yy = fullElectVirial_yy; )
  FULL( params.fullElectVirial_yz = fullElectVirial_yz; )
  FULL( params.fullElectVirial_zz = fullElectVirial_zz; )
  FAST( ENERGY( params.vdwEnergy = vdwEnergy; ) )
  FAST( ENERGY( SHORT( params.electEnergy = electEnergy; ) ) )
  FULL( ENERGY( params.fullElectEnergy = fullElectEnergy; ) )

  #undef CALC_TYPE
}


// Undefine atom and force macros for SELF computes
#if (0 SELF(+1))
  #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
    #undef p_1
    #undef pExt_1
    #undef f_1
    #if (0 FULL(+1))
      #undef fullf_1
    #endif
  #else
    #undef p_1_x
    #undef p_1_y
    #undef p_1_z
    #undef p_1_q
    #undef p_1_vdwType
    #undef p_1_index
    #undef p_1_exclIndex
    #undef p_1_exclMaxDiff
    #undef f_1_x
    #undef f_1_y
    #undef f_1_z
    #undef f_1_w
    #if (0 FULL(+1))
      #undef fullf_1_x
      #undef fullf_1_y
      #undef fullf_1_z
      #undef fullf_1_w
    #endif
  #endif
#endif

#undef __MIC_PAD_PLGEN_CTRL

#else  // NAMD_MIC

#include "ComputeNonbondedMICKernelBase2.h"

#endif  // NAMD_MIC

