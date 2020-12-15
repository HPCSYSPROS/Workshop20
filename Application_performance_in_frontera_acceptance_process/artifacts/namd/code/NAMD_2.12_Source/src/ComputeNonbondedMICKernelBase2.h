
#ifdef NAMD_MIC

{ // Start a block so variables can be declared outsude the force loop

  // Validate the control macros for this file, by ensuring that certain invalid
  //   combinations result in compiler errors, preventing those combinations
  EXCLUDED( FAST( foo bar ) )
  EXCLUDED( MODIFIED( foo bar ) )
  EXCLUDED( NORMAL( foo bar ) )
  NORMAL( MODIFIED( foo bar ) )
  ALCHPAIR( NOT_ALCHPAIR( foo bar ) )

  // If pairlist refinement is enabled, then refine the pairlists by doing the exact
  //   distance checks for all the atom pairs in the current pairlist as a separate
  //   loop before the actual force computation loop (removing the distance check
  //   from the force computation loop and ensuring that all iterations of the force
  //   computation loop are within cutoff).
  #if REFINE_PAIRLISTS != 0

    int plSize = *(params.plSize_ptr);
    int * RESTRICT plArray = *(params.plArray_ptr);
    double * RESTRICT r2Array = *(params.r2Array_ptr);

    // Grow the scratch arrays as required
    #if REFINE_PAIRLISTS_XYZ != 0
      const int plMinSize = ((PAIRLIST[1] - 2) + 7) & (~(7));  // Roundup to multiple of 16
      const int plMallocSize = 4 * plMinSize;
    #else
      const int plMinSize = PAIRLIST[1] - 2;
      const int plMallocSize = plMinSize;
    #endif
    if (plMallocSize > plSize) {
      if (plArray != NULL) { _mm_free(plArray); }
      if (r2Array != NULL) { _mm_free(r2Array); }
      const int newPlSize = plMallocSize * 1.2; // r2
      plArray = (   int*)_mm_malloc(newPlSize * sizeof(   int), 64); __ASSERT(plArray != NULL);
      r2Array = (double*)_mm_malloc(newPlSize * sizeof(double), 64); __ASSERT(r2Array != NULL);
      *(params.plArray_ptr) = plArray;
      *(params.r2Array_ptr) = r2Array;
      *(params.plSize_ptr) = newPlSize;
    }

    __ASSUME_ALIGNED(r2Array); // NOTE: plArray handled below (common code)
    #if REFINE_PAIRLISTS_XYZ != 0
      double * RESTRICT p_ij_x_array = r2Array + (1 * plMinSize); __ASSUME_ALIGNED(p_ij_x_array);
      double * RESTRICT p_ij_y_array = r2Array + (2 * plMinSize); __ASSUME_ALIGNED(p_ij_y_array);
      double * RESTRICT p_ij_z_array = r2Array + (3 * plMinSize); __ASSUME_ALIGNED(p_ij_z_array);
    #endif

    #if REFINE_PAIRLIST_HANDCODE != 0
    {

      __m512i plI_vec = _mm512_set_16to16_epi32(0, 0, 0, 0, 0, 0, 0, 0,
                                                7, 6, 5, 4, 3, 2, 1, 0);
      __m512i ij_mask_vec = _mm512_set_1to16_epi32(0x0000FFFF);
      __m512i ij_store_perm_pattern = _mm512_set_16to16_epi32( 9,  7,  9,  6,  9,  5,  9,  4,
                                                               9,  3,  9,  2,  9,  1,  9,  0  );
      int offset = 0;
      plSize = PAIRLIST[1] - 2;
      int * RESTRICT orig_plArray = PAIRLIST + 2;

      #pragma loop count (100)
      #pragma novector
      for (int plI = 0; plI < plMinSize; plI += 8) {

        // Create the active_mask
        __mmask16 active_mask = _mm512_mask_cmplt_epi32_mask(_mm512_int2mask(0x00FF), plI_vec, _mm512_set_1to16_epi32(plSize));

        // Load the i and j values from the pairlist array
        // NOTE: The "hi" part of this "unaligned load" should never be required since plArray is actually
        //   aligned (unaligned load being used because of 32-bit indexes vs 64-bit forces).
        __m512i ij_vec = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), active_mask, orig_plArray + plI     );
                ij_vec = _mm512_mask_loadunpackhi_epi32(                ij_vec, active_mask, orig_plArray + plI + 16);
        __m512i i_vec = _mm512_and_epi32(_mm512_mask_srli_epi32(_mm512_setzero_epi32(), active_mask, ij_vec, 16), ij_mask_vec);
        __m512i j_vec = _mm512_mask_and_epi32(_mm512_setzero_epi32(), active_mask, ij_vec, ij_mask_vec);

        // NOTE: Make these the same size as pointers so that they can be directly loaded in to
        //   the index register that can be used with a base pointer register.
        uintptr_t tmpI[8] __attribute__((aligned(64)));
        uintptr_t tmpJ[8] __attribute__((aligned(64)));
        _mm512_store_epi32(tmpI, _mm512_permutevar_epi32(ij_store_perm_pattern, i_vec));
        _mm512_store_epi32(tmpJ, _mm512_permutevar_epi32(ij_store_perm_pattern, j_vec));

        // Increment the vectorized loop counter
        plI_vec = _mm512_mask_add_epi32(plI_vec, _mm512_int2mask(0x00FF), plI_vec, _mm512_set_1to16_epi32(8));

        // Load position/charge data for the i and j atoms
        __m512d p_i_x_vec, p_i_y_vec, p_i_z_vec, p_j_x_vec, p_j_y_vec, p_j_z_vec;
        {
          __m512d zero_vec = _mm512_setzero_pd();
          p_i_x_vec = p_i_y_vec = p_i_z_vec = zero_vec;
          p_j_x_vec = p_j_y_vec = p_j_z_vec = zero_vec;

          // NOTE: If 1 asm statement per 1 vcvtps2pd instrunction with each using the same mask (_k), the
          //   compiler resorts to using multiple masks instead of just one (using a single statment to avoid).
          #define __VCVTPS2PD_LOAD_ELEM(i) \
          { __mmask16 _k = _mm512_int2mask(0x01 << (i)); uintptr_t _i = tmpI[(i)]; uintptr_t _j = tmpJ[(i)]; \
            asm("vcvtps2pd (%6,%12,4){{1to8}}, %0{{%14}}\n" \
                "vcvtps2pd (%7,%12,4){{1to8}}, %1{{%14}}\n" \
                "vcvtps2pd (%8,%12,4){{1to8}}, %2{{%14}}\n" \
                "vcvtps2pd (%9,%13,4){{1to8}}, %3{{%14}}\n" \
                "vcvtps2pd (%10,%13,4){{1to8}}, %4{{%14}}\n" \
                "vcvtps2pd (%11,%13,4){{1to8}}, %5{{%14}}\n" \
                : "+v"(p_i_x_vec), "+v"(p_i_y_vec), "+v"(p_i_z_vec), \
                  "+v"(p_j_x_vec), "+v"(p_j_y_vec), "+v"(p_j_z_vec)  \
                : "r"(p_0_x), "r"(p_0_y), "r"(p_0_z), \
                  "r"(p_1_x), "r"(p_1_y), "r"(p_1_z), \
                  "r"(_i), "r"(_j), "Yk"(_k) \
               ); \
          }
          __VCVTPS2PD_LOAD_ELEM(0);
          __VCVTPS2PD_LOAD_ELEM(1);
          __VCVTPS2PD_LOAD_ELEM(2);
          __VCVTPS2PD_LOAD_ELEM(3);
          __VCVTPS2PD_LOAD_ELEM(4);
          __VCVTPS2PD_LOAD_ELEM(5);
          __VCVTPS2PD_LOAD_ELEM(6);
          __VCVTPS2PD_LOAD_ELEM(7);
          #undef __VCVTPS2PD_LOAD_ELEM
        }

        p_i_x_vec = _mm512_add_pd(p_i_x_vec, _mm512_set_1to8_pd(params.offset.x));
        p_i_y_vec = _mm512_add_pd(p_i_y_vec, _mm512_set_1to8_pd(params.offset.y));
        p_i_z_vec = _mm512_add_pd(p_i_z_vec, _mm512_set_1to8_pd(params.offset.z));
        __m512d p_ij_x_vec = _mm512_sub_pd(p_i_x_vec, p_j_x_vec);
        __m512d p_ij_y_vec = _mm512_sub_pd(p_i_y_vec, p_j_y_vec);
        __m512d p_ij_z_vec = _mm512_sub_pd(p_i_z_vec, p_j_z_vec);

        __m512d r2_vec = _mm512_add_pd(_mm512_mul_pd(p_ij_x_vec, p_ij_x_vec), _mm512_set_1to8_pd(r2_delta));
        r2_vec = _mm512_add_pd(_mm512_mul_pd(p_ij_y_vec, p_ij_y_vec), r2_vec);
        r2_vec = _mm512_add_pd(_mm512_mul_pd(p_ij_z_vec, p_ij_z_vec), r2_vec);

        __mmask16 cutoff_mask = _mm512_mask_cmplt_pd_mask(active_mask, r2_vec, _mm512_set_1to8_pd(cutoff2_delta));

        _mm512_mask_packstorelo_epi32(plArray + offset     , cutoff_mask, ij_vec);
        _mm512_mask_packstorehi_epi32(plArray + offset + 16, cutoff_mask, ij_vec);
        _mm512_mask_packstorelo_pd(r2Array + offset    , cutoff_mask, r2_vec);
        _mm512_mask_packstorehi_pd(r2Array + offset + 8, cutoff_mask, r2_vec);
        #if REFINE_PAIRLISTS_XYZ != 0
          _mm512_mask_packstorelo_pd(p_ij_x_array + offset    , cutoff_mask, p_ij_x_vec);
          _mm512_mask_packstorehi_pd(p_ij_x_array + offset + 8, cutoff_mask, p_ij_x_vec);
          _mm512_mask_packstorelo_pd(p_ij_y_array + offset    , cutoff_mask, p_ij_y_vec);
          _mm512_mask_packstorehi_pd(p_ij_y_array + offset + 8, cutoff_mask, p_ij_y_vec);
          _mm512_mask_packstorelo_pd(p_ij_z_array + offset    , cutoff_mask, p_ij_z_vec);
          _mm512_mask_packstorehi_pd(p_ij_z_array + offset + 8, cutoff_mask, p_ij_z_vec);
        #endif

        offset += _mm512_mask_reduce_add_epi32(cutoff_mask, _mm512_set_1to16_epi32(1));
      }

      plSize = offset;

    }
    #else

      // Refine the saved pairlist
      const int savedPlSize = PAIRLIST[1];
      plSize = 0;
      // DMK - NOTE | TODO | FIXME : If just plSize is used, the compiler doesn't recognize
      //   the compress idiom, so used r2Size separately, and once the compiler recognizes
      //   this, fix and/or report it.  For now, use two index variables so compress happens.
      int r2Size = 0;
      #if REFINE_PAIRLISTS_XYZ != 0
        int ijXSize = 0;
        int ijYSize = 0;
        int ijZSize = 0;
      #endif
      for (int plI = 2; plI < savedPlSize; plI++) {
        const int plEntry = PAIRLIST[plI];
        const int i = (plEntry >> 16) & 0xFFFF;
        const int j = (plEntry      ) & 0xFFFF;
        double p_ij_x = (p_0_x[i] + params.offset.x) - p_1_x[j];
        double p_ij_y = (p_0_y[i] + params.offset.y) - p_1_y[j];
        double p_ij_z = (p_0_z[i] + params.offset.z) - p_1_z[j];
        double r2_ij_delta = r2_delta + (p_ij_x * p_ij_x) + (p_ij_y * p_ij_y) + (p_ij_z * p_ij_z);
        if (r2_ij_delta < cutoff2_delta) {
          plArray[plSize++] = plEntry;
          r2Array[r2Size++] = r2_ij_delta;
          //plSize++;
          #if REFINE_PAIRLISTS_XYZ != 0
            p_ij_x_array[ijXSize++] = p_ij_x;
            p_ij_y_array[ijYSize++] = p_ij_y;
            p_ij_z_array[ijZSize++] = p_ij_z;
          #endif
        }
      }

    #endif // MIC_HANDCODE_FORCE check

  #else  // REFINE_PAIRLISTS != 0

    const int plSize = PAIRLIST[1] - 2;
    const int * RESTRICT plArray = PAIRLIST + 2;

  #endif  // REFINE_PAIRLISTS != 0

  // NOTE : The above code defines plArray and plSize, which are the pointer and the
  //   length of the pairlist to be processed by the force computation loop, whether
  //   or not pairlist refinement is being used.  Neither the pointer or the size should
  //   include the first two elements (allocated size and actual size) of the pairlists
  //   that are saved across iterations.
  __ASSERT(plArray != NULL);
  __ASSERT(plSize >= 0);
  __ASSUME_ALIGNED(plArray);

  // Include code for the force computation loop itself
  #if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
    #include "ComputeNonbondedMICKernelBase2_handcode_single.h"
  #else
    #include "ComputeNonbondedMICKernelBase2_scalar.h"
  #endif

} // End block

#else  // NAMD_MIC

#include "ComputeNonbondedMICKernelBase2_handcode_single.h"
#include "ComputeNonbondedMICKernelBase2_scalar.h"

#endif  // NAMD_MIC

