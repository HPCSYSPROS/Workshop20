#ifdef NAMD_MIC

  #define GATHER_PS_I32_OFFSET(v, p, i, o) \
  { \
    __mmask16 k0 = _mm512_int2mask(0x0008); \
    __mmask16 k1 = _mm512_int2mask(0x0080); \
    __mmask16 k2 = _mm512_int2mask(0x0800); \
    __mmask16 k3 = _mm512_int2mask(0x8000); \
    (v) = _mm512_mask_loadunpacklo_ps(                      _mm512_setzero_ps(), k3, &((p)[(i)[13] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k2, &((p)[(i)[ 9] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k1, &((p)[(i)[ 5] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k0, &((p)[(i)[ 1] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(_mm512_swizzle_ps((v), _MM_SWIZ_REG_CDAB), k3, &((p)[(i)[12] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k2, &((p)[(i)[ 8] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k1, &((p)[(i)[ 4] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k0, &((p)[(i)[ 0] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(_mm512_swizzle_ps((v), _MM_SWIZ_REG_BADC), k3, &((p)[(i)[14] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k2, &((p)[(i)[10] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k1, &((p)[(i)[ 6] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k0, &((p)[(i)[ 2] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(_mm512_swizzle_ps((v), _MM_SWIZ_REG_CDAB), k3, &((p)[(i)[15] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k2, &((p)[(i)[11] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k1, &((p)[(i)[ 7] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_ps(                                      (v), k0, &((p)[(i)[ 3] + (o)])); \
  }
  #define GATHER_PS_I32(v, p, i)  GATHER_PS_I32_OFFSET((v), (p), (i), 0)


  #define GATHER_EPI32_I32_OFFSET(v, p, i, o) \
  { \
    __mmask16 k0 = _mm512_int2mask(0x0008); \
    __mmask16 k1 = _mm512_int2mask(0x0080); \
    __mmask16 k2 = _mm512_int2mask(0x0800); \
    __mmask16 k3 = _mm512_int2mask(0x8000); \
    (v) = _mm512_mask_loadunpacklo_epi32(                      _mm512_setzero_epi32(), k3, &((p)[(i)[13] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k2, &((p)[(i)[ 9] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k1, &((p)[(i)[ 5] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k0, &((p)[(i)[ 1] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(_mm512_swizzle_epi32((v), _MM_SWIZ_REG_CDAB), k3, &((p)[(i)[12] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k2, &((p)[(i)[ 8] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k1, &((p)[(i)[ 4] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k0, &((p)[(i)[ 0] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(_mm512_swizzle_epi32((v), _MM_SWIZ_REG_BADC), k3, &((p)[(i)[14] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k2, &((p)[(i)[10] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k1, &((p)[(i)[ 6] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k0, &((p)[(i)[ 2] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(_mm512_swizzle_epi32((v), _MM_SWIZ_REG_CDAB), k3, &((p)[(i)[15] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k2, &((p)[(i)[11] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k1, &((p)[(i)[ 7] + (o)])); \
    (v) = _mm512_mask_loadunpacklo_epi32(                                         (v), k0, &((p)[(i)[ 3] + (o)])); \
  }
  #define GATHER_EPI32_I32(v, p, i)  GATHER_EPI32_I32_OFFSET((v), (p), (i), 0)


  #define SCATTER_INC_PD_I32_STEP(v, p, i, k, idx, pattern) \
    t = _mm512_mask_loadunpacklo_pd(t, (k), &((p)[(i)[(idx)]])); \
    t = _mm512_mask_add_pd(t, (k), _mm512_swizzle_pd((v), (pattern), t); \
    _mm512_mask_packstorelo_pd(&((p)[(i)[(idx)]]), (k), t);

  #define SCATTER_INC_PD_I32(v, p, i) \
  { \
    __mmask16 k0 = _mm512_int2mask(0x02); \
    __mmask16 k1 = _mm512_int2mask(0x20); \
    __m512d t = _mm512_setzero_pd(); \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  0, _MM_SWIZ_REG_CDAB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  1, _MM_SWIZ_REG_DCBA) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  2, _MM_SWIZ_REG_DACB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  3, _MM_SWIZ_REG_BADC) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  4, _MM_SWIZ_REG_CDAB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  5, _MM_SWIZ_REG_DCBA) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  6, _MM_SWIZ_REG_DACB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  7, _MM_SWIZ_REG_BADC) \
  }


  #define SCATTER_INC_PS_I32_STEP(v, p, i, k, idx, pattern) \
    t = _mm512_mask_loadunpacklo_ps(t, (k), &((p)[(i)[(idx)]])); \
    t = _mm512_mask_add_ps(t, (k), _mm512_swizzle_ps((v), (pattern), t); \
    _mm512_mask_packstorelo_ps(&((p)[(i)[(idx)]]), (k), t);
    
  #define SCATTER_INC_PS_I32(v, p, i) \
  { \
    __mmask16 k0 = _mm512_int2mask(0x0008); \
    __mmask16 k1 = _mm512_int2mask(0x0080); \
    __mmask16 k2 = _mm512_int2mask(0x0800); \
    __mmask16 k3 = _mm512_int2mask(0x8000); \
    __m512 t = _mm512_setzero_ps(); \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  0, _MM_SWIZ_REG_CDAB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  1, _MM_SWIZ_REG_DCBA) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  2, _MM_SWIZ_REG_DACB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k0,  3, _MM_SWIZ_REG_BADC) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  4, _MM_SWIZ_REG_CDAB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  5, _MM_SWIZ_REG_DCBA) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  6, _MM_SWIZ_REG_DACB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k1,  7, _MM_SWIZ_REG_BADC) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k2,  8, _MM_SWIZ_REG_CDAB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k2,  9, _MM_SWIZ_REG_DCBA) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k2, 10, _MM_SWIZ_REG_DACB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k2, 11, _MM_SWIZ_REG_BADC) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k3, 12, _MM_SWIZ_REG_CDAB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k3, 13, _MM_SWIZ_REG_DCBA) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k3, 14, _MM_SWIZ_REG_DACB) \
    SCATTER_INC_PS_I32_STEP((v), (p), (i), k3, 15, _MM_SWIZ_REG_BADC) \
  }

  // DMK - NOTE - The instructions for these macros are they way they are so
  //   values are converted to double before accumulating (e.g. rather than
  //   sum the floats, convert, and then accumulate).
  #define CONTRIB_ADD_PS2PD(pd_vec, ps_vec) \
  { __m512 ps_tmp_vec = (ps_vec); \
    (pd_vec) = _mm512_add_pd((pd_vec), _mm512_cvtpslo_pd(ps_tmp_vec)); \
    (pd_vec) = _mm512_add_pd((pd_vec), _mm512_cvtpslo_pd(_mm512_permute4f128_ps(ps_tmp_vec, _MM_PERM_BADC))); \
  }
  #define CONTRIB_SUB_PS2PD(pd_vec, ps_vec) \
  { __m512 ps_tmp_vec = (ps_vec); \
    (pd_vec) = _mm512_sub_pd((pd_vec), _mm512_cvtpslo_pd(ps_tmp_vec)); \
    (pd_vec) = _mm512_sub_pd((pd_vec), _mm512_cvtpslo_pd(_mm512_permute4f128_ps(ps_tmp_vec, _MM_PERM_BADC))); \
  }


  #define APPLY_FORCES_PS2PD_STEP_ADD(f_x, f_y, f_z, i, v_x, v_y, v_z, p, k) \
    tx = _mm512_mask_loadunpacklo_pd(tx, (k), (f_x) + (i)); \
    ty = _mm512_mask_loadunpacklo_pd(ty, (k), (f_y) + (i)); \
    tz = _mm512_mask_loadunpacklo_pd(tz, (k), (f_z) + (i)); \
    tx = _mm512_mask_add_pd(tx, (k), tx, _mm512_swizzle_pd((v_x), (p))); \
    ty = _mm512_mask_add_pd(ty, (k), ty, _mm512_swizzle_pd((v_y), (p))); \
    tz = _mm512_mask_add_pd(tz, (k), tz, _mm512_swizzle_pd((v_z), (p))); \
    _mm512_mask_packstorelo_pd((f_x) + (i), (k), tx); \
    _mm512_mask_packstorelo_pd((f_y) + (i), (k), ty); \
    _mm512_mask_packstorelo_pd((f_z) + (i), (k), tz);

  #define APPLY_FORCES_PS2PD_STEP_SUB(f_x, f_y, f_z, i, v_x, v_y, v_z, p, k) \
    tx = _mm512_mask_loadunpacklo_pd(tx, (k), (f_x) + (i)); \
    ty = _mm512_mask_loadunpacklo_pd(ty, (k), (f_y) + (i)); \
    tz = _mm512_mask_loadunpacklo_pd(tz, (k), (f_z) + (i)); \
    tx = _mm512_mask_sub_pd(tx, (k), tx, _mm512_swizzle_pd((v_x), (p))); \
    ty = _mm512_mask_sub_pd(ty, (k), ty, _mm512_swizzle_pd((v_y), (p))); \
    tz = _mm512_mask_sub_pd(tz, (k), tz, _mm512_swizzle_pd((v_z), (p))); \
    _mm512_mask_packstorelo_pd((f_x) + (i), (k), tx); \
    _mm512_mask_packstorelo_pd((f_y) + (i), (k), ty); \
    _mm512_mask_packstorelo_pd((f_z) + (i), (k), tz);

  //#if MIC_PAD_PLGEN != 0
  #if __MIC_PAD_PLGEN_CTRL != 0

    #if MIC_HANDCODE_FORCE_COMBINE_FORCES != 0

      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0

        #define APPLY_FORCES_PS2PD_JSTEP(j_lo, j_hi, v, fv) \
        { \
          FAST(SHORT(  v_tmp = _mm512_mask_loadunpacklo_pd( v_tmp, k_lo,     f_1 + tmpJ32[(j_lo)]); )) \
          FULL(       fv_tmp = _mm512_mask_loadunpacklo_pd(fv_tmp, k_lo, fullf_1 + tmpJ32[(j_lo)]);  ) \
          FAST(SHORT(  v_tmp = _mm512_sub_pd( v_tmp, ( v)); )) \
          FULL(       fv_tmp = _mm512_sub_pd(fv_tmp, (fv));  ) \
          FAST(SHORT(          _mm512_mask_packstorelo_pd(    f_1 + tmpJ32[(j_lo)], k_lo,  v_tmp); )) \
          FULL(                _mm512_mask_packstorelo_pd(fullf_1 + tmpJ32[(j_lo)], k_lo, fv_tmp);  ) \
          FAST(SHORT(  v_tmp = _mm512_mask_loadunpacklo_pd( v_tmp, k_hi,     f_1 + tmpJ32[(j_hi)]); )) \
          FULL(       fv_tmp = _mm512_mask_loadunpacklo_pd(fv_tmp, k_hi, fullf_1 + tmpJ32[(j_hi)]);  ) \
          FAST(SHORT(  v_tmp = _mm512_sub_pd( v_tmp, ( v)); )) \
          FULL(       fv_tmp = _mm512_sub_pd(fv_tmp, (fv));  ) \
          FAST(SHORT(          _mm512_mask_packstorelo_pd(    f_1 + tmpJ32[(j_hi)], k_hi,  v_tmp); )) \
          FULL(                _mm512_mask_packstorelo_pd(fullf_1 + tmpJ32[(j_hi)], k_hi, fv_tmp);  ) \
	}

        #define APPLY_FORCES_PS2PD \
          { \
          /* Set 'w' values (exclusion checksum) to -1.0f when counting so 'sub' operation for x,y,z,w actually adds 1.0f */ \
          FAST(SHORT( __m512     tmp_w_vec = _mm512_mask_mov_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_set_1to16_ps(0.0f MODIFIED(- 1.0f) EXCLUDED(- 1.0f))); )) \
          FULL(       __m512 fulltmp_w_vec = _mm512_setzero_ps();  ) \
          /* Transpose the values so that each set of 4 lanes has x,y,z,w values for a given interaction. */ \
          /*   NOTE: This rearranges the values via a 4x4 transpose.                                      */ \
	  __mmask16 k_2x2_0 = _mm512_int2mask(0xAAAA); \
          __mmask16 k_2x2_1 = _mm512_int2mask(0x5555); \
          FAST(SHORT( __m512i tmp_a0 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(    tmp_x_vec), k_2x2_0, _mm512_castps_si512(    tmp_y_vec), _MM_SWIZ_REG_CDAB); )) \
          FAST(SHORT( __m512i tmp_a1 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(    tmp_y_vec), k_2x2_1, _mm512_castps_si512(    tmp_x_vec), _MM_SWIZ_REG_CDAB); )) \
          FAST(SHORT( __m512i tmp_a2 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(    tmp_z_vec), k_2x2_0, _mm512_castps_si512(    tmp_w_vec), _MM_SWIZ_REG_CDAB); )) \
          FAST(SHORT( __m512i tmp_a3 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(    tmp_w_vec), k_2x2_1, _mm512_castps_si512(    tmp_z_vec), _MM_SWIZ_REG_CDAB); )) \
          FULL(       __m512i tmp_b0 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(fulltmp_x_vec), k_2x2_0, _mm512_castps_si512(fulltmp_y_vec), _MM_SWIZ_REG_CDAB);  ) \
          FULL(       __m512i tmp_b1 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(fulltmp_y_vec), k_2x2_1, _mm512_castps_si512(fulltmp_x_vec), _MM_SWIZ_REG_CDAB);  ) \
          FULL(       __m512i tmp_b2 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(fulltmp_z_vec), k_2x2_0, _mm512_castps_si512(fulltmp_w_vec), _MM_SWIZ_REG_CDAB);  ) \
          FULL(       __m512i tmp_b3 = _mm512_mask_swizzle_epi32(_mm512_castps_si512(fulltmp_w_vec), k_2x2_1, _mm512_castps_si512(fulltmp_z_vec), _MM_SWIZ_REG_CDAB);  ) \
	  __mmask16 k_4x4_0 = _mm512_int2mask(0xCCCC); \
          __mmask16 k_4x4_1 = _mm512_int2mask(0x3333); \
          FAST(SHORT( __m512     tmp_0_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_a0, k_4x4_0, tmp_a2, _MM_SWIZ_REG_BADC)); )) \
          FAST(SHORT( __m512     tmp_1_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_a1, k_4x4_0, tmp_a3, _MM_SWIZ_REG_BADC)); )) \
          FAST(SHORT( __m512     tmp_2_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_a2, k_4x4_1, tmp_a0, _MM_SWIZ_REG_BADC)); )) \
          FAST(SHORT( __m512     tmp_3_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_a3, k_4x4_1, tmp_a1, _MM_SWIZ_REG_BADC)); )) \
          FULL(       __m512 fulltmp_0_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b0, k_4x4_0, tmp_b2, _MM_SWIZ_REG_BADC));  ) \
          FULL(       __m512 fulltmp_1_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b1, k_4x4_0, tmp_b3, _MM_SWIZ_REG_BADC));  ) \
          FULL(       __m512 fulltmp_2_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b2, k_4x4_1, tmp_b0, _MM_SWIZ_REG_BADC));  ) \
          FULL(       __m512 fulltmp_3_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b3, k_4x4_1, tmp_b1, _MM_SWIZ_REG_BADC));  ) \
          /* Convert the floats to doubles.  NOTE: This must be done prior to any math because of the difference in magnitudes. */ \
          FAST(SHORT( __m512d  v_0_lo = _mm512_cvtpslo_pd(                           tmp_0_vec                ); )) \
          FAST(SHORT( __m512d  v_0_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_0_vec, _MM_PERM_BADC)); )) \
          FAST(SHORT( __m512d  v_1_lo = _mm512_cvtpslo_pd(                           tmp_1_vec                ); )) \
          FAST(SHORT( __m512d  v_1_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_1_vec, _MM_PERM_BADC)); )) \
          FAST(SHORT( __m512d  v_2_lo = _mm512_cvtpslo_pd(                           tmp_2_vec                ); )) \
          FAST(SHORT( __m512d  v_2_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_2_vec, _MM_PERM_BADC)); )) \
          FAST(SHORT( __m512d  v_3_lo = _mm512_cvtpslo_pd(                           tmp_3_vec                ); )) \
          FAST(SHORT( __m512d  v_3_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_3_vec, _MM_PERM_BADC)); )) \
          FULL(       __m512d fv_0_lo = _mm512_cvtpslo_pd(                       fulltmp_0_vec                );  ) \
          FULL(       __m512d fv_0_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_0_vec, _MM_PERM_BADC));  ) \
          FULL(       __m512d fv_1_lo = _mm512_cvtpslo_pd(                       fulltmp_1_vec                );  ) \
          FULL(       __m512d fv_1_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_1_vec, _MM_PERM_BADC));  ) \
          FULL(       __m512d fv_2_lo = _mm512_cvtpslo_pd(                       fulltmp_2_vec                );  ) \
          FULL(       __m512d fv_2_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_2_vec, _MM_PERM_BADC));  ) \
          FULL(       __m512d fv_3_lo = _mm512_cvtpslo_pd(                       fulltmp_3_vec                );  ) \
          FULL(       __m512d fv_3_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_3_vec, _MM_PERM_BADC));  ) \
	  /* Apply the forces to the 'j' atoms. */ \
          __mmask16 k_lo = _mm512_int2mask(0x000F); \
          __mmask16 k_hi = _mm512_int2mask(0x00F0); \
          __m512d  v_tmp = _mm512_setzero_pd(); \
          __m512d fv_tmp = _mm512_setzero_pd(); \
          APPLY_FORCES_PS2PD_JSTEP( 0,  4, v_0_lo, fv_0_lo); \
          APPLY_FORCES_PS2PD_JSTEP( 1,  5, v_1_lo, fv_1_lo); \
          APPLY_FORCES_PS2PD_JSTEP( 2,  6, v_2_lo, fv_2_lo); \
          APPLY_FORCES_PS2PD_JSTEP( 3,  7, v_3_lo, fv_3_lo); \
          APPLY_FORCES_PS2PD_JSTEP( 8, 12, v_0_hi, fv_0_hi); \
          APPLY_FORCES_PS2PD_JSTEP( 9, 13, v_1_hi, fv_1_hi); \
          APPLY_FORCES_PS2PD_JSTEP(10, 14, v_2_hi, fv_2_hi); \
          APPLY_FORCES_PS2PD_JSTEP(11, 15, v_3_hi, fv_3_hi); \
          /* Reduce (add) the forces into a single force vector and apply it to the 'i' atom. */ \
          FAST(SHORT( __m512d  v_tmp0_xyzw = _mm512_add_pd( v_0_lo,  v_0_hi); )) \
          FULL(       __m512d fv_tmp0_xyzw = _mm512_add_pd(fv_0_lo, fv_0_hi);  ) \
          FAST(SHORT( __m512d  v_tmp1_xyzw = _mm512_add_pd( v_1_lo,  v_1_hi); )) \
	  FULL(       __m512d fv_tmp1_xyzw = _mm512_add_pd(fv_1_lo, fv_1_hi);  ) \
          FAST(SHORT( __m512d  v_tmp2_xyzw = _mm512_add_pd( v_2_lo,  v_2_hi); )) \
          FULL(       __m512d fv_tmp2_xyzw = _mm512_add_pd(fv_2_lo, fv_2_hi);  ) \
          FAST(SHORT( __m512d  v_tmp3_xyzw = _mm512_add_pd( v_3_lo,  v_3_hi); )) \
          FULL(       __m512d fv_tmp3_xyzw = _mm512_add_pd(fv_3_lo, fv_3_hi);  ) \
          FAST(SHORT( __m512d  v_tmp4_xyzw = _mm512_add_pd( v_tmp0_xyzw,  v_tmp1_xyzw); )) \
          FULL(       __m512d fv_tmp4_xyzw = _mm512_add_pd(fv_tmp0_xyzw, fv_tmp1_xyzw);  ) \
          FAST(SHORT( __m512d  v_tmp5_xyzw = _mm512_add_pd( v_tmp2_xyzw,  v_tmp3_xyzw); )) \
          FULL(       __m512d fv_tmp5_xyzw = _mm512_add_pd(fv_tmp2_xyzw, fv_tmp3_xyzw);  ) \
          FAST(SHORT( __m512d  v_tmp6_xyzw = _mm512_add_pd( v_tmp4_xyzw,  v_tmp5_xyzw); )) \
          FULL(       __m512d fv_tmp6_xyzw = _mm512_add_pd(fv_tmp4_xyzw, fv_tmp5_xyzw);  ) \
          FAST(SHORT( __m512d  v_xyzw = _mm512_add_pd( v_tmp6_xyzw, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512( v_tmp6_xyzw), _MM_PERM_BADC))); )) \
          FULL(       __m512d fv_xyzw = _mm512_add_pd(fv_tmp6_xyzw, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(fv_tmp6_xyzw), _MM_PERM_BADC)));  ) \
          FAST(SHORT(MODIFIED(  v_xyzw = _mm512_mask_sub_pd( v_xyzw, _mm512_int2mask(0x88), _mm512_setzero_pd(),  v_xyzw); ))) /* NOTE: w *= -1.0 for add operation below */ \
          FAST(SHORT(EXCLUDED(  v_xyzw = _mm512_mask_sub_pd( v_xyzw, _mm512_int2mask(0x88), _mm512_setzero_pd(),  v_xyzw); ))) \
          FULL(MODIFIED(       fv_xyzw = _mm512_mask_sub_pd(fv_xyzw, _mm512_int2mask(0x88), _mm512_setzero_pd(), fv_xyzw);  )) \
          FULL(EXCLUDED(       fv_xyzw = _mm512_mask_sub_pd(fv_xyzw, _mm512_int2mask(0x88), _mm512_setzero_pd(), fv_xyzw);  )) \
          int iTest_index = _mm_tzcnt_32(_mm512_mask2int(cutoff_mask)); \
          uintptr_t iTest_i = tmpI32[iTest_index]; \
          FAST(SHORT( __m512d  i_xyzw = _mm512_mask_loadunpacklo_pd(_mm512_setzero_pd(), k_lo,     f_0 + iTest_i); )) \
          FULL(       __m512d fi_xyzw = _mm512_mask_loadunpacklo_pd(_mm512_setzero_pd(), k_lo, fullf_0 + iTest_i);  ) \
          FAST(SHORT(          i_xyzw = _mm512_mask_add_pd( i_xyzw, k_lo,  i_xyzw,  v_xyzw); )) \
          FULL(               fi_xyzw = _mm512_mask_add_pd(fi_xyzw, k_lo, fi_xyzw, fv_xyzw);  ) \
          FAST(SHORT(                   _mm512_mask_packstorelo_pd(    f_0 + iTest_i, k_lo,  i_xyzw); )) \
          FULL(                         _mm512_mask_packstorelo_pd(fullf_0 + iTest_i, k_lo, fi_xyzw);  ) \
          }

      #else // if MIC_HANDCODE_FORC_SOA_VS_AOS != 0

        #define APPLY_FORCES_PS2PD \
        { \
          FAST(SHORT( __m512d  v_x_lo = _mm512_cvtpslo_pd(                           tmp_x_vec                ); )) \
          FAST(SHORT( __m512d  v_x_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_x_vec, _MM_PERM_BADC)); )) \
          FAST(SHORT( __m512d  v_y_lo = _mm512_cvtpslo_pd(                           tmp_y_vec                ); )) \
          FAST(SHORT( __m512d  v_y_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_y_vec, _MM_PERM_BADC)); )) \
          FAST(SHORT( __m512d  v_z_lo = _mm512_cvtpslo_pd(                           tmp_z_vec                ); )) \
          FAST(SHORT( __m512d  v_z_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_z_vec, _MM_PERM_BADC)); )) \
          FULL(       __m512d fv_x_lo = _mm512_cvtpslo_pd(                       fulltmp_x_vec                );  ) \
          FULL(       __m512d fv_x_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_x_vec, _MM_PERM_BADC));  ) \
          FULL(       __m512d fv_y_lo = _mm512_cvtpslo_pd(                       fulltmp_y_vec                );  ) \
          FULL(       __m512d fv_y_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_y_vec, _MM_PERM_BADC));  ) \
          FULL(       __m512d fv_z_lo = _mm512_cvtpslo_pd(                       fulltmp_z_vec                );  ) \
          FULL(       __m512d fv_z_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_z_vec, _MM_PERM_BADC));  ) \
          int iTest_index = _mm_tzcnt_32(_mm512_mask2int(cutoff_mask)); \
          uintptr_t iTest_i = tmpI32[iTest_index]; \
          __mmask16 cutoff_lo = _mm512_kand(cutoff_mask, _mm512_int2mask(0xFF)); \
          __mmask16 cutoff_hi = _mm512_kmerge2l1h(cutoff_mask, _mm512_kxor(cutoff_mask, cutoff_mask)); \
          FAST(SHORT(     f_0_x[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd( v_x_hi,  v_x_lo)); )) \
          FAST(SHORT(     f_0_y[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd( v_y_hi,  v_y_lo)); )) \
          FAST(SHORT(     f_0_z[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd( v_z_hi,  v_z_lo)); )) \
          FULL(       fullf_0_x[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd(fv_x_hi, fv_x_lo));  ) \
          FULL(       fullf_0_y[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd(fv_y_hi, fv_y_lo));  ) \
          FULL(       fullf_0_z[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd(fv_z_hi, fv_z_lo));  ) \
          __m512i j_lo_vec = _mm512_mask_mov_epi32(_mm512_setzero_epi32(), cutoff_lo, j_vec); \
          __m512i j_hi_vec = _mm512_mask_permute4f128_epi32(_mm512_setzero_epi32(), cutoff_hi, j_vec, _MM_PERM_BADC); \
          FAST(SHORT( __m512d  f_j_x_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec,     f_1_x, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_y_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec,     f_1_y, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_z_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec,     f_1_z, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_x_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec,     f_1_x, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_y_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec,     f_1_y, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_z_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec,     f_1_z, _MM_SCALE_8); )) \
          FULL(       __m512d ff_j_x_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, fullf_1_x, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_y_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, fullf_1_y, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_z_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, fullf_1_z, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_x_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, fullf_1_x, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_y_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, fullf_1_y, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_z_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, fullf_1_z, _MM_SCALE_8);  ) \
          FAST(SHORT(  f_j_x_lo_vec = _mm512_mask_sub_pd( f_j_x_lo_vec, cutoff_lo,  f_j_x_lo_vec,  v_x_lo); )) \
          FAST(SHORT(  f_j_y_lo_vec = _mm512_mask_sub_pd( f_j_y_lo_vec, cutoff_lo,  f_j_y_lo_vec,  v_y_lo); )) \
          FAST(SHORT(  f_j_z_lo_vec = _mm512_mask_sub_pd( f_j_z_lo_vec, cutoff_lo,  f_j_z_lo_vec,  v_z_lo); )) \
          FAST(SHORT(  f_j_x_hi_vec = _mm512_mask_sub_pd( f_j_x_hi_vec, cutoff_hi,  f_j_x_hi_vec,  v_x_hi); )) \
          FAST(SHORT(  f_j_y_hi_vec = _mm512_mask_sub_pd( f_j_y_hi_vec, cutoff_hi,  f_j_y_hi_vec,  v_y_hi); )) \
          FAST(SHORT(  f_j_z_hi_vec = _mm512_mask_sub_pd( f_j_z_hi_vec, cutoff_hi,  f_j_z_hi_vec,  v_z_hi); )) \
          FULL(       ff_j_x_lo_vec = _mm512_mask_sub_pd(ff_j_x_lo_vec, cutoff_lo, ff_j_x_lo_vec, fv_x_lo);  ) \
          FULL(       ff_j_y_lo_vec = _mm512_mask_sub_pd(ff_j_y_lo_vec, cutoff_lo, ff_j_y_lo_vec, fv_y_lo);  ) \
          FULL(       ff_j_z_lo_vec = _mm512_mask_sub_pd(ff_j_z_lo_vec, cutoff_lo, ff_j_z_lo_vec, fv_z_lo);  ) \
          FULL(       ff_j_x_hi_vec = _mm512_mask_sub_pd(ff_j_x_hi_vec, cutoff_hi, ff_j_x_hi_vec, fv_x_hi);  ) \
          FULL(       ff_j_y_hi_vec = _mm512_mask_sub_pd(ff_j_y_hi_vec, cutoff_hi, ff_j_y_hi_vec, fv_y_hi);  ) \
          FULL(       ff_j_z_hi_vec = _mm512_mask_sub_pd(ff_j_z_hi_vec, cutoff_hi, ff_j_z_hi_vec, fv_z_hi);  ) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_x, cutoff_lo, j_lo_vec,  f_j_x_lo_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_y, cutoff_lo, j_lo_vec,  f_j_y_lo_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_z, cutoff_lo, j_lo_vec,  f_j_z_lo_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_x, cutoff_hi, j_hi_vec,  f_j_x_hi_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_y, cutoff_hi, j_hi_vec,  f_j_y_hi_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_z, cutoff_hi, j_hi_vec,  f_j_z_hi_vec, _MM_SCALE_8); )) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_x, cutoff_lo, j_lo_vec, ff_j_x_lo_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_y, cutoff_lo, j_lo_vec, ff_j_y_lo_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_z, cutoff_lo, j_lo_vec, ff_j_z_lo_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_x, cutoff_hi, j_hi_vec, ff_j_x_hi_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_y, cutoff_hi, j_hi_vec, ff_j_y_hi_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_z, cutoff_hi, j_hi_vec, ff_j_z_hi_vec, _MM_SCALE_8);  ) \
        }

      #endif // MIC_HANDCODE_FORCE_SOA_VS_AOS != 0

    #else // MIC_HANDCODE_FORCE_COMBINE_FORCES != 0

      #define APPLY_FORCES_PS2PD(v_x, v_y, v_z, f_i_x, f_i_y, f_i_z, f_j_x, f_j_y, f_j_z, i, j) \
      { \
        __m512d v_x_lo = _mm512_cvtpslo_pd(v_x); \
        __m512d v_x_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps((v_x), _MM_PERM_BADC)); \
        __m512d v_y_lo = _mm512_cvtpslo_pd(v_y); \
        __m512d v_y_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps((v_y), _MM_PERM_BADC)); \
        __m512d v_z_lo = _mm512_cvtpslo_pd(v_z); \
        __m512d v_z_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps((v_z), _MM_PERM_BADC)); \
        __mmask16 cutoff_lo = _mm512_kand(cutoff_mask, _mm512_int2mask(0xFF)); \
        __mmask16 cutoff_hi = _mm512_kmerge2l1h(cutoff_mask, _mm512_kxor(cutoff_mask, cutoff_mask)); \
        __m512i j_lo_vec = _mm512_mask_mov_epi32(_mm512_setzero_epi32(), cutoff_lo, j_vec); \
        __m512i j_hi_vec = _mm512_mask_permute4f128_epi32(_mm512_setzero_epi32(), cutoff_hi, j_vec, _MM_PERM_BADC); \
        int iTest_index = _mm_tzcnt_32(_mm512_mask2int(cutoff_mask)); \
        uintptr_t iTest_i = tmpI32[iTest_index]; \
        f_i_x[iTest_i] += _mm512_mask_reduce_add_pd(cutoff_hi, v_x_hi) + _mm512_mask_reduce_add_pd(cutoff_lo, v_x_lo); \
        f_i_y[iTest_i] += _mm512_mask_reduce_add_pd(cutoff_hi, v_y_hi) + _mm512_mask_reduce_add_pd(cutoff_lo, v_y_lo); \
        f_i_z[iTest_i] += _mm512_mask_reduce_add_pd(cutoff_hi, v_z_hi) + _mm512_mask_reduce_add_pd(cutoff_lo, v_z_lo); \
        __m512d f_j_x_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, (f_j_x), _MM_SCALE_8); \
        __m512d f_j_y_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, (f_j_y), _MM_SCALE_8); \
        __m512d f_j_z_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, (f_j_z), _MM_SCALE_8); \
        __m512d f_j_x_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, (f_j_x), _MM_SCALE_8); \
        __m512d f_j_y_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, (f_j_y), _MM_SCALE_8); \
        __m512d f_j_z_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, (f_j_z), _MM_SCALE_8); \
        f_j_x_lo_vec = _mm512_mask_sub_pd(_mm512_setzero_pd(), cutoff_lo, f_j_x_lo_vec, v_x_lo); \
        f_j_y_lo_vec = _mm512_mask_sub_pd(_mm512_setzero_pd(), cutoff_lo, f_j_y_lo_vec, v_y_lo); \
        f_j_z_lo_vec = _mm512_mask_sub_pd(_mm512_setzero_pd(), cutoff_lo, f_j_z_lo_vec, v_z_lo); \
        f_j_x_hi_vec = _mm512_mask_sub_pd(_mm512_setzero_pd(), cutoff_hi, f_j_x_hi_vec, v_x_hi); \
        f_j_y_hi_vec = _mm512_mask_sub_pd(_mm512_setzero_pd(), cutoff_hi, f_j_y_hi_vec, v_y_hi); \
        f_j_z_hi_vec = _mm512_mask_sub_pd(_mm512_setzero_pd(), cutoff_hi, f_j_z_hi_vec, v_z_hi); \
        _mm512_mask_i32loscatter_pd((f_j_x), cutoff_lo, j_lo_vec, f_j_x_lo_vec, _MM_SCALE_8); \
        _mm512_mask_i32loscatter_pd((f_j_y), cutoff_lo, j_lo_vec, f_j_y_lo_vec, _MM_SCALE_8); \
        _mm512_mask_i32loscatter_pd((f_j_z), cutoff_lo, j_lo_vec, f_j_z_lo_vec, _MM_SCALE_8); \
        _mm512_mask_i32loscatter_pd((f_j_x), cutoff_hi, j_hi_vec, f_j_x_hi_vec, _MM_SCALE_8); \
        _mm512_mask_i32loscatter_pd((f_j_y), cutoff_hi, j_hi_vec, f_j_y_hi_vec, _MM_SCALE_8); \
        _mm512_mask_i32loscatter_pd((f_j_z), cutoff_hi, j_hi_vec, f_j_z_hi_vec, _MM_SCALE_8); \
      }

    #endif // MIC_HANDCODE_FORCE_COMBINE_FORCES != 0

  #else // MIC_PAD_PLGEN != 0

    #if MIC_HANDCODE_FORCE_COMBINE_FORCES != 0

      #define APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x, v_y, v_z, fv_x, fv_y, fv_z, i, p, k) \
        FAST(SHORT(  tx = _mm512_mask_loadunpacklo_pd( tx, (k),     f_0_x + (i)); )) \
        FAST(SHORT(  ty = _mm512_mask_loadunpacklo_pd( ty, (k),     f_0_y + (i)); )) \
        FAST(SHORT(  tz = _mm512_mask_loadunpacklo_pd( tz, (k),     f_0_z + (i)); )) \
        FULL(       ftx = _mm512_mask_loadunpacklo_pd(ftx, (k), fullf_0_x + (i));  ) \
        FULL(       fty = _mm512_mask_loadunpacklo_pd(fty, (k), fullf_0_y + (i));  ) \
        FULL(       ftz = _mm512_mask_loadunpacklo_pd(ftz, (k), fullf_0_z + (i));  ) \
        FAST(SHORT(  tx = _mm512_mask_add_pd( tx, (k),  tx, _mm512_swizzle_pd(( v_x), (p))); )) \
        FAST(SHORT(  ty = _mm512_mask_add_pd( ty, (k),  ty, _mm512_swizzle_pd(( v_y), (p))); )) \
        FAST(SHORT(  tz = _mm512_mask_add_pd( tz, (k),  tz, _mm512_swizzle_pd(( v_z), (p))); )) \
        FULL(       ftx = _mm512_mask_add_pd(ftx, (k), ftx, _mm512_swizzle_pd((fv_x), (p)));  ) \
        FULL(       fty = _mm512_mask_add_pd(fty, (k), fty, _mm512_swizzle_pd((fv_y), (p)));  ) \
        FULL(       ftz = _mm512_mask_add_pd(ftz, (k), ftz, _mm512_swizzle_pd((fv_z), (p)));  ) \
        FAST(SHORT( _mm512_mask_packstorelo_pd(    f_0_x + (i), (k),  tx); )) \
        FAST(SHORT( _mm512_mask_packstorelo_pd(    f_0_y + (i), (k),  ty); )) \
        FAST(SHORT( _mm512_mask_packstorelo_pd(    f_0_z + (i), (k),  tz); )) \
        FULL(       _mm512_mask_packstorelo_pd(fullf_0_x + (i), (k), ftx);  ) \
        FULL(       _mm512_mask_packstorelo_pd(fullf_0_y + (i), (k), fty);  ) \
        FULL(       _mm512_mask_packstorelo_pd(fullf_0_z + (i), (k), ftz);  )

      #define APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x, v_y, v_z, fv_x, fv_y, fv_z, i, p, k) \
        FAST(SHORT(  tx = _mm512_mask_loadunpacklo_pd( tx, (k),     f_1_x + (i)); )) \
        FAST(SHORT(  ty = _mm512_mask_loadunpacklo_pd( ty, (k),     f_1_y + (i)); )) \
        FAST(SHORT(  tz = _mm512_mask_loadunpacklo_pd( tz, (k),     f_1_z + (i)); )) \
        FULL(       ftx = _mm512_mask_loadunpacklo_pd(ftx, (k), fullf_1_x + (i));  ) \
        FULL(       fty = _mm512_mask_loadunpacklo_pd(fty, (k), fullf_1_y + (i));  ) \
        FULL(       ftz = _mm512_mask_loadunpacklo_pd(ftz, (k), fullf_1_z + (i));  ) \
        FAST(SHORT(  tx = _mm512_mask_sub_pd( tx, (k),  tx, _mm512_swizzle_pd(( v_x), (p))); )) \
        FAST(SHORT(  ty = _mm512_mask_sub_pd( ty, (k),  ty, _mm512_swizzle_pd(( v_y), (p))); )) \
        FAST(SHORT(  tz = _mm512_mask_sub_pd( tz, (k),  tz, _mm512_swizzle_pd(( v_z), (p))); )) \
        FULL(       ftx = _mm512_mask_sub_pd(ftx, (k), ftx, _mm512_swizzle_pd((fv_x), (p)));  ) \
        FULL(       fty = _mm512_mask_sub_pd(fty, (k), fty, _mm512_swizzle_pd((fv_y), (p)));  ) \
        FULL(       ftz = _mm512_mask_sub_pd(ftz, (k), ftz, _mm512_swizzle_pd((fv_z), (p)));  ) \
        FAST(SHORT( _mm512_mask_packstorelo_pd(    f_1_x + (i), (k),  tx); )) \
        FAST(SHORT( _mm512_mask_packstorelo_pd(    f_1_y + (i), (k),  ty); )) \
        FAST(SHORT( _mm512_mask_packstorelo_pd(    f_1_z + (i), (k),  tz); )) \
        FULL(       _mm512_mask_packstorelo_pd(fullf_1_x + (i), (k), ftx);  ) \
        FULL(       _mm512_mask_packstorelo_pd(fullf_1_y + (i), (k), fty);  ) \
        FULL(       _mm512_mask_packstorelo_pd(fullf_1_z + (i), (k), ftz);  )

      #define APPLY_FORCES_PS2PD \
      { \
        __mmask16 k_lo = _mm512_int2mask(0x02); \
        __mmask16 k_hi = _mm512_int2mask(0x20); \
        FAST(SHORT( __m512d  tx = _mm512_setzero_pd(); )) \
        FAST(SHORT( __m512d  ty = _mm512_setzero_pd(); )) \
        FAST(SHORT( __m512d  tz = _mm512_setzero_pd(); )) \
        FULL(       __m512d ftx = _mm512_setzero_pd();  ) \
        FULL(       __m512d fty = _mm512_setzero_pd();  ) \
        FULL(       __m512d ftz = _mm512_setzero_pd();  ) \
        FAST(SHORT( __m512d  v_x_lo = _mm512_cvtpslo_pd(                           tmp_x_vec                ); )) \
        FAST(SHORT( __m512d  v_x_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_x_vec, _MM_PERM_BADC)); )) \
        FAST(SHORT( __m512d  v_y_lo = _mm512_cvtpslo_pd(                           tmp_y_vec                ); )) \
        FAST(SHORT( __m512d  v_y_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_y_vec, _MM_PERM_BADC)); )) \
        FAST(SHORT( __m512d  v_z_lo = _mm512_cvtpslo_pd(                           tmp_z_vec                ); )) \
        FAST(SHORT( __m512d  v_z_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(    tmp_z_vec, _MM_PERM_BADC)); )) \
        FULL(       __m512d fv_x_lo = _mm512_cvtpslo_pd(                       fulltmp_x_vec                );  ) \
        FULL(       __m512d fv_x_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_x_vec, _MM_PERM_BADC));  ) \
        FULL(       __m512d fv_y_lo = _mm512_cvtpslo_pd(                       fulltmp_y_vec                );  ) \
        FULL(       __m512d fv_y_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_y_vec, _MM_PERM_BADC));  ) \
        FULL(       __m512d fv_z_lo = _mm512_cvtpslo_pd(                       fulltmp_z_vec                );  ) \
        FULL(       __m512d fv_z_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(fulltmp_z_vec, _MM_PERM_BADC));  ) \
        int iTest_index = _mm_tzcnt_32(_mm512_mask2int(cutoff_mask)); \
        uintptr_t iTest_i = tmpI32[iTest_index]; \
        __mmask16 iTest_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, _mm512_set_1to16_epi32(iTest_i), i_vec); \
        iTest_mask = _mm512_kxor(iTest_mask, cutoff_mask); \
        if (_mm512_kortestz(iTest_mask, iTest_mask)) { \
          __mmask16 cutoff_lo = _mm512_kand(cutoff_mask, _mm512_int2mask(0xFF)); \
          __mmask16 cutoff_hi = _mm512_kmerge2l1h(cutoff_mask, _mm512_kxor(cutoff_mask, cutoff_mask)); \
          FAST(SHORT(     f_0_x[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd( v_x_hi,  v_x_lo)); )) \
          FAST(SHORT(     f_0_y[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd( v_y_hi,  v_y_lo)); )) \
          FAST(SHORT(     f_0_z[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd( v_z_hi,  v_z_lo)); )) \
          FULL(       fullf_0_x[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd(fv_x_hi, fv_x_lo));  ) \
          FULL(       fullf_0_y[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd(fv_y_hi, fv_y_lo));  ) \
          FULL(       fullf_0_z[iTest_i] += _mm512_reduce_add_pd(_mm512_add_pd(fv_z_hi, fv_z_lo));  ) \
          __m512i j_lo_vec = _mm512_mask_mov_epi32(_mm512_setzero_epi32(), cutoff_lo, j_vec); \
          __m512i j_hi_vec = _mm512_mask_permute4f128_epi32(_mm512_setzero_epi32(), cutoff_hi, j_vec, _MM_PERM_BADC); \
          FAST(SHORT( __m512d  f_j_x_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec,     f_1_x, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_y_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec,     f_1_y, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_z_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec,     f_1_z, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_x_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec,     f_1_x, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_y_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec,     f_1_y, _MM_SCALE_8); )) \
          FAST(SHORT( __m512d  f_j_z_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec,     f_1_z, _MM_SCALE_8); )) \
          FULL(       __m512d ff_j_x_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, fullf_1_x, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_y_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, fullf_1_y, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_z_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, fullf_1_z, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_x_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, fullf_1_x, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_y_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, fullf_1_y, _MM_SCALE_8);  ) \
          FULL(       __m512d ff_j_z_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, fullf_1_z, _MM_SCALE_8);  ) \
          FAST(SHORT(  f_j_x_lo_vec = _mm512_mask_sub_pd( f_j_x_lo_vec, cutoff_lo,  f_j_x_lo_vec,  v_x_lo); )) \
          FAST(SHORT(  f_j_y_lo_vec = _mm512_mask_sub_pd( f_j_y_lo_vec, cutoff_lo,  f_j_y_lo_vec,  v_y_lo); )) \
          FAST(SHORT(  f_j_z_lo_vec = _mm512_mask_sub_pd( f_j_z_lo_vec, cutoff_lo,  f_j_z_lo_vec,  v_z_lo); )) \
          FAST(SHORT(  f_j_x_hi_vec = _mm512_mask_sub_pd( f_j_x_hi_vec, cutoff_hi,  f_j_x_hi_vec,  v_x_hi); )) \
          FAST(SHORT(  f_j_y_hi_vec = _mm512_mask_sub_pd( f_j_y_hi_vec, cutoff_hi,  f_j_y_hi_vec,  v_y_hi); )) \
          FAST(SHORT(  f_j_z_hi_vec = _mm512_mask_sub_pd( f_j_z_hi_vec, cutoff_hi,  f_j_z_hi_vec,  v_z_hi); )) \
          FULL(       ff_j_x_lo_vec = _mm512_mask_sub_pd(ff_j_x_lo_vec, cutoff_lo, ff_j_x_lo_vec, fv_x_lo);  ) \
          FULL(       ff_j_y_lo_vec = _mm512_mask_sub_pd(ff_j_y_lo_vec, cutoff_lo, ff_j_y_lo_vec, fv_y_lo);  ) \
          FULL(       ff_j_z_lo_vec = _mm512_mask_sub_pd(ff_j_z_lo_vec, cutoff_lo, ff_j_z_lo_vec, fv_z_lo);  ) \
          FULL(       ff_j_x_hi_vec = _mm512_mask_sub_pd(ff_j_x_hi_vec, cutoff_hi, ff_j_x_hi_vec, fv_x_hi);  ) \
          FULL(       ff_j_y_hi_vec = _mm512_mask_sub_pd(ff_j_y_hi_vec, cutoff_hi, ff_j_y_hi_vec, fv_y_hi);  ) \
          FULL(       ff_j_z_hi_vec = _mm512_mask_sub_pd(ff_j_z_hi_vec, cutoff_hi, ff_j_z_hi_vec, fv_z_hi);  ) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_x, cutoff_lo, j_lo_vec,  f_j_x_lo_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_y, cutoff_lo, j_lo_vec,  f_j_y_lo_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_z, cutoff_lo, j_lo_vec,  f_j_z_lo_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_x, cutoff_hi, j_hi_vec,  f_j_x_hi_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_y, cutoff_hi, j_hi_vec,  f_j_y_hi_vec, _MM_SCALE_8); )) \
          FAST(SHORT( _mm512_mask_i32loscatter_pd(    f_1_z, cutoff_hi, j_hi_vec,  f_j_z_hi_vec, _MM_SCALE_8); )) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_x, cutoff_lo, j_lo_vec, ff_j_x_lo_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_y, cutoff_lo, j_lo_vec, ff_j_y_lo_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_z, cutoff_lo, j_lo_vec, ff_j_z_lo_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_x, cutoff_hi, j_hi_vec, ff_j_x_hi_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_y, cutoff_hi, j_hi_vec, ff_j_y_hi_vec, _MM_SCALE_8);  ) \
          FULL(       _mm512_mask_i32loscatter_pd(fullf_1_z, cutoff_hi, j_hi_vec, ff_j_z_hi_vec, _MM_SCALE_8);  ) \
        } else { \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 0], _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 1], _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 2], _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 3], _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 4], _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 5], _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 6], _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpI32[ 7], _MM_SWIZ_REG_BADC, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[ 8], _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[ 9], _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[10], _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[11], _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[12], _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[13], _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[14], _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpI32[15], _MM_SWIZ_REG_BADC, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 0], _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 1], _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 2], _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 3], _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 4], _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 5], _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 6], _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_lo, v_y_lo, v_z_lo, fv_x_lo, fv_y_lo, fv_z_lo, tmpJ32[ 7], _MM_SWIZ_REG_BADC, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[ 8], _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[ 9], _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[10], _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[11], _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[12], _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[13], _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[14], _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB_COMBO(v_x_hi, v_y_hi, v_z_hi, fv_x_hi, fv_y_hi, fv_z_hi, tmpJ32[15], _MM_SWIZ_REG_BADC, k_hi); \
        } \
      }
      
    #else // MIC_HANDCODE_FORCE_COMBINE_FORCES != 0

      #define APPLY_FORCES_PS2PD(v_x, v_y, v_z, f_i_x, f_i_y, f_i_z, f_j_x, f_j_y, f_j_z, i, j) \
      { \
        __mmask16 k_lo = _mm512_int2mask(0x02); \
        __mmask16 k_hi = _mm512_int2mask(0x20); \
        __m512d tx = _mm512_setzero_pd(); \
        __m512d ty = _mm512_setzero_pd(); \
        __m512d tz = _mm512_setzero_pd(); \
        __m512d v_x_lo = _mm512_cvtpslo_pd(v_x); \
        __m512d v_x_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps((v_x), _MM_PERM_BADC)); \
        __m512d v_y_lo = _mm512_cvtpslo_pd(v_y); \
        __m512d v_y_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps((v_y), _MM_PERM_BADC)); \
        __m512d v_z_lo = _mm512_cvtpslo_pd(v_z); \
        __m512d v_z_hi = _mm512_cvtpslo_pd(_mm512_permute4f128_ps((v_z), _MM_PERM_BADC)); \
        int iTest_index = _mm_tzcnt_32(_mm512_mask2int(cutoff_mask)); \
        uintptr_t iTest_i = tmpI32[iTest_index]; \
        __mmask16 iTest_mask = _mm512_mask_cmpeq_epi32_mask(cutoff_mask, _mm512_set_1to16_epi32(iTest_i), i_vec); \
        iTest_mask = _mm512_kxor(iTest_mask, cutoff_mask); \
        if (_mm512_kortestz(iTest_mask, iTest_mask)) { \
          __mmask16 cutoff_lo = _mm512_kand(cutoff_mask, _mm512_int2mask(0xFF)); \
          __mmask16 cutoff_hi = _mm512_kmerge2l1h(cutoff_mask, _mm512_kxor(cutoff_mask, cutoff_mask)); \
          f_i_x[iTest_i] += _mm512_mask_reduce_add_pd(cutoff_hi, v_x_hi) + _mm512_mask_reduce_add_pd(cutoff_lo, v_x_lo); \
          f_i_y[iTest_i] += _mm512_mask_reduce_add_pd(cutoff_hi, v_y_hi) + _mm512_mask_reduce_add_pd(cutoff_lo, v_y_lo); \
          f_i_z[iTest_i] += _mm512_mask_reduce_add_pd(cutoff_hi, v_z_hi) + _mm512_mask_reduce_add_pd(cutoff_lo, v_z_lo); \
          __m512i j_lo_vec = _mm512_mask_mov_epi32(_mm512_setzero_epi32(), cutoff_lo, j_vec); \
          __m512i j_hi_vec = _mm512_mask_permute4f128_epi32(_mm512_setzero_epi32(), cutoff_hi, j_vec, _MM_PERM_BADC); \
          __m512d f_j_x_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, (f_j_x), _MM_SCALE_8); \
          __m512d f_j_y_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, (f_j_y), _MM_SCALE_8); \
          __m512d f_j_z_lo_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_lo, j_lo_vec, (f_j_z), _MM_SCALE_8); \
          __m512d f_j_x_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, (f_j_x), _MM_SCALE_8); \
          __m512d f_j_y_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, (f_j_y), _MM_SCALE_8); \
          __m512d f_j_z_hi_vec = _mm512_mask_i32logather_pd(_mm512_setzero_pd(), cutoff_hi, j_hi_vec, (f_j_z), _MM_SCALE_8); \
          f_j_x_lo_vec = _mm512_mask_sub_pd(f_j_x_lo_vec, cutoff_lo, f_j_x_lo_vec, v_x_lo); \
          f_j_y_lo_vec = _mm512_mask_sub_pd(f_j_y_lo_vec, cutoff_lo, f_j_y_lo_vec, v_y_lo); \
          f_j_z_lo_vec = _mm512_mask_sub_pd(f_j_z_lo_vec, cutoff_lo, f_j_z_lo_vec, v_z_lo); \
          f_j_x_hi_vec = _mm512_mask_sub_pd(f_j_x_hi_vec, cutoff_hi, f_j_x_hi_vec, v_x_hi); \
          f_j_y_hi_vec = _mm512_mask_sub_pd(f_j_y_hi_vec, cutoff_hi, f_j_y_hi_vec, v_y_hi); \
          f_j_z_hi_vec = _mm512_mask_sub_pd(f_j_z_hi_vec, cutoff_hi, f_j_z_hi_vec, v_z_hi); \
          _mm512_mask_i32loscatter_pd((f_j_x), cutoff_lo, j_lo_vec, f_j_x_lo_vec, _MM_SCALE_8); \
          _mm512_mask_i32loscatter_pd((f_j_y), cutoff_lo, j_lo_vec, f_j_y_lo_vec, _MM_SCALE_8); \
          _mm512_mask_i32loscatter_pd((f_j_z), cutoff_lo, j_lo_vec, f_j_z_lo_vec, _MM_SCALE_8); \
          _mm512_mask_i32loscatter_pd((f_j_x), cutoff_hi, j_hi_vec, f_j_x_hi_vec, _MM_SCALE_8); \
          _mm512_mask_i32loscatter_pd((f_j_y), cutoff_hi, j_hi_vec, f_j_y_hi_vec, _MM_SCALE_8); \
          _mm512_mask_i32loscatter_pd((f_j_z), cutoff_hi, j_hi_vec, f_j_z_hi_vec, _MM_SCALE_8); \
        } else { \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 0], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 1], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 2], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 3], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 4], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 5], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 6], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 7], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_BADC, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 8], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[ 9], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[10], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[11], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[12], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[13], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[14], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_ADD(f_i_x, f_i_y, f_i_z, (i)[15], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_BADC, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 0], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 1], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 2], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 3], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 4], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 5], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 6], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 7], v_x_lo, v_y_lo, v_z_lo, _MM_SWIZ_REG_BADC, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 8], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_CDAB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[ 9], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DCBA, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[10], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DACB, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[11], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_BADC, k_lo); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[12], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_CDAB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[13], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DCBA, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[14], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_DACB, k_hi); \
          APPLY_FORCES_PS2PD_STEP_SUB(f_j_x, f_j_y, f_j_z, (j)[15], v_x_hi, v_y_hi, v_z_hi, _MM_SWIZ_REG_BADC, k_hi); \
        } \
      }

    #endif // MIC_HANDCODE_FORCE_COMBINE_FORCES != 0
  #endif // MIC_PAD_PLGEN != 0

  __m512i plI_vec = _mm512_set_16to16_epi32(15, 14, 13, 12, 11, 10,  9,  8,
                                             7,  6,  5,  4,  3,  2,  1,  0);
  const int plSize_16 = (plSize + 15) & (~15); // Round up to a multiple of 16

  FAST(ENERGY( __m512d vdwEnergy_vec = _mm512_setzero_pd(); ))
  FAST(SHORT(ENERGY( __m512d electEnergy_vec = _mm512_setzero_pd(); )))
  FULL(ENERGY( __m512d fullElectEnergy_vec = _mm512_setzero_pd(); ))
  #if (0 PAIR(FAST(SHORT(+1))))
    __m512d virial_xx_vec = _mm512_setzero_pd();
    __m512d virial_xy_vec = _mm512_setzero_pd();
    __m512d virial_xz_vec = _mm512_setzero_pd();
    __m512d virial_yy_vec = _mm512_setzero_pd();
    __m512d virial_yz_vec = _mm512_setzero_pd();
    __m512d virial_zz_vec = _mm512_setzero_pd();
  #endif
  #if (0 PAIR(FULL(+1)))
    __m512d fullElectVirial_xx_vec = _mm512_setzero_pd();
    __m512d fullElectVirial_xy_vec = _mm512_setzero_pd();
    __m512d fullElectVirial_xz_vec = _mm512_setzero_pd();
    __m512d fullElectVirial_yy_vec = _mm512_setzero_pd();
    __m512d fullElectVirial_yz_vec = _mm512_setzero_pd();
    __m512d fullElectVirial_zz_vec = _mm512_setzero_pd();
  #endif

  __m512i ij_mask_vec = _mm512_set_1to16_epi32(0x0000FFFF);
  __m512i ij_store_perm_pattern = _mm512_set_16to16_epi32( 9,  7,  9,  6,  9,  5,  9,  4,
                                                           9,  3,  9,  2,  9,  1,  9,  0  );

  #if (MIC_EXCL_CHECKSUM_FULL != 0) && (0 EXCLUDED(+1) MODIFIED(+1))
    __m512i exclusionSum_vec = _mm512_setzero_epi32();
  #endif

  #if (0 NORMAL(+1))
    #if (0 PAIR(+1))
      #pragma loop count (500)
    #else
      #pragma loop count (7000)
    #endif
  #else
    #if (0 PAIR(+1))
      #pragma loop count (2)
    #else
      #pragma loop count (30)
    #endif
  #endif
  #pragma prefetch plArray:_MM_HINT_NTA
  #pragma novector
  for (int plI = 0; plI < plSize_16; plI += 16) {

    // Create the active_mask
    __mmask16 active_mask = _mm512_cmplt_epi32_mask(plI_vec, _mm512_set_1to16_epi32(plSize));

    #if MIC_HANDCODE_FORCE_PFDIST != 0
      __m512i future_ij_vec = _mm512_load_epi32(plArray + plI + (16 * MIC_HANDCODE_FORCE_PFDIST));
      __mmask16 future_ij_mask = _mm512_cmpneq_epi32_mask(future_ij_vec, _mm512_set_1to16_epi32(-1));
      __m512i future_j_vec = _mm512_and_epi32(future_ij_vec, ij_mask_vec);
      _mm512_mask_prefetch_i32gather_ps(future_j_vec, future_ij_mask, p_1_x, _MM_SCALE_4, _MM_HINT_T0);
      _mm512_mask_prefetch_i32gather_ps(future_j_vec, future_ij_mask, p_1_y, _MM_SCALE_4, _MM_HINT_T0);
      _mm512_mask_prefetch_i32gather_ps(future_j_vec, future_ij_mask, p_1_z, _MM_SCALE_4, _MM_HINT_T0);
    #endif

    // Load the i and j values from the pairlist array and modify the active_mask accordingly
    __m512i ij_vec = _mm512_mask_load_epi32(_mm512_setzero_epi32(), active_mask, plArray + plI);
    #if __MIC_PAD_PLGEN_CTRL != 0
      __m512i ij_mask_vec = _mm512_set_1to16_epi32(0xFFFF);
      __m512i ij_lo_vec = _mm512_and_epi32(ij_vec, ij_mask_vec);
      active_mask = _mm512_mask_cmpneq_epi32_mask(active_mask, ij_lo_vec, ij_mask_vec);
    #endif
    __m512i i_vec = _mm512_and_epi32(_mm512_mask_srli_epi32(_mm512_setzero_epi32(), active_mask, ij_vec, 16), ij_mask_vec);
    __m512i j_vec = _mm512_mask_and_epi32(_mm512_setzero_epi32(), active_mask, ij_vec, ij_mask_vec);

    // Store the index out to an array of scalar values so they can be individually accessed by lane
    uintptr_t tmpI32[16] __attribute__((aligned(64)));
    uintptr_t tmpJ32[16] __attribute__((aligned(64)));
    _mm512_store_epi64(tmpI32    , _mm512_mask_permutevar_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0x5555), ij_store_perm_pattern,                           i_vec                ));
    _mm512_store_epi64(tmpI32 + 8, _mm512_mask_permutevar_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0x5555), ij_store_perm_pattern, _mm512_permute4f128_epi32(i_vec, _MM_PERM_BADC)));
    _mm512_store_epi64(tmpJ32    , _mm512_mask_permutevar_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0x5555), ij_store_perm_pattern,                           j_vec                ));
    _mm512_store_epi64(tmpJ32 + 8, _mm512_mask_permutevar_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0x5555), ij_store_perm_pattern, _mm512_permute4f128_epi32(j_vec, _MM_PERM_BADC)));

    // Increment the vectorized loop counter
    plI_vec = _mm512_add_epi32(plI_vec, _mm512_set_1to16_epi32(16));

    // Load position/charge data for the i and j atoms
    __m512 p_i_x_vec, p_i_y_vec, p_i_z_vec, p_i_q_vec, p_j_x_vec, p_j_y_vec, p_j_z_vec, p_j_q_vec;
    #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
      __m512i pExt_i_vdwType_vec, pExt_j_vdwType_vec;
    #endif
    {
      __m512 zero_vec = _mm512_setzero_ps();

      int iTest_index = _mm_tzcnt_32(_mm512_mask2int(active_mask)); // Index of lease significant 1 bit in active mask (NOTE: active_mask should be not all zeros)
      unsigned int iTest_i = tmpI32[iTest_index];  // NOTE: iTest_i = a valid "i" value
      #if __MIC_PAD_PLGEN_CTRL != 0
      #else
      __mmask16 iTest_mask = _mm512_mask_cmpeq_epi32_mask(active_mask, _mm512_set_1to16_epi32(iTest_i), i_vec);
      iTest_mask =_mm512_kxor(iTest_mask, active_mask);
      if (_mm512_kortestz(iTest_mask, iTest_mask)) {
      #endif
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          p_i_x_vec = _mm512_set_1to16_ps(p_0[iTest_i].x);
          p_i_y_vec = _mm512_set_1to16_ps(p_0[iTest_i].y);
          p_i_z_vec = _mm512_set_1to16_ps(p_0[iTest_i].z);
          p_i_q_vec = _mm512_set_1to16_ps(p_0[iTest_i].charge);
          {
            __mmask16 k0 = _mm512_int2mask(0x000F);
            __mmask16 k1 = _mm512_int2mask(0x00F0);
            __mmask16 k2 = _mm512_int2mask(0x0F00);
            __mmask16 k3 = _mm512_int2mask(0xF000);
            __m512i tmp_a0 = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), k0, p_1 + tmpJ32[ 0]);
                    tmp_a0 = _mm512_mask_loadunpacklo_epi32(                tmp_a0, k1, p_1 + tmpJ32[ 4]);
                    tmp_a0 = _mm512_mask_loadunpacklo_epi32(                tmp_a0, k2, p_1 + tmpJ32[ 8]);
                    tmp_a0 = _mm512_mask_loadunpacklo_epi32(                tmp_a0, k3, p_1 + tmpJ32[12]);
            __m512i tmp_a1 = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), k0, p_1 + tmpJ32[ 1]);
                    tmp_a1 = _mm512_mask_loadunpacklo_epi32(                tmp_a1, k1, p_1 + tmpJ32[ 5]);
                    tmp_a1 = _mm512_mask_loadunpacklo_epi32(                tmp_a1, k2, p_1 + tmpJ32[ 9]);
                    tmp_a1 = _mm512_mask_loadunpacklo_epi32(                tmp_a1, k3, p_1 + tmpJ32[13]);
            __m512i tmp_a2 = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), k0, p_1 + tmpJ32[ 2]);
                    tmp_a2 = _mm512_mask_loadunpacklo_epi32(                tmp_a2, k1, p_1 + tmpJ32[ 6]);
                    tmp_a2 = _mm512_mask_loadunpacklo_epi32(                tmp_a2, k2, p_1 + tmpJ32[10]);
                    tmp_a2 = _mm512_mask_loadunpacklo_epi32(                tmp_a2, k3, p_1 + tmpJ32[14]);
            __m512i tmp_a3 = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), k0, p_1 + tmpJ32[ 3]);
                    tmp_a3 = _mm512_mask_loadunpacklo_epi32(                tmp_a3, k1, p_1 + tmpJ32[ 7]);
                    tmp_a3 = _mm512_mask_loadunpacklo_epi32(                tmp_a3, k2, p_1 + tmpJ32[11]);
                    tmp_a3 = _mm512_mask_loadunpacklo_epi32(                tmp_a3, k3, p_1 + tmpJ32[15]);
            __mmask16 k_2x2_0 = _mm512_int2mask(0xAAAA);
            __mmask16 k_2x2_1 = _mm512_int2mask(0x5555);
            __m512i tmp_b0 = _mm512_mask_swizzle_epi32(tmp_a0, k_2x2_0, tmp_a1, _MM_SWIZ_REG_CDAB);
            __m512i tmp_b1 = _mm512_mask_swizzle_epi32(tmp_a1, k_2x2_1, tmp_a0, _MM_SWIZ_REG_CDAB);
            __m512i tmp_b2 = _mm512_mask_swizzle_epi32(tmp_a2, k_2x2_0, tmp_a3, _MM_SWIZ_REG_CDAB);
            __m512i tmp_b3 = _mm512_mask_swizzle_epi32(tmp_a3, k_2x2_1, tmp_a2, _MM_SWIZ_REG_CDAB);
            __mmask16 k_4x4_0 = _mm512_int2mask(0xCCCC);
            __mmask16 k_4x4_1 = _mm512_int2mask(0x3333);
            p_j_x_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b0, k_4x4_0, tmp_b2, _MM_SWIZ_REG_BADC));
            p_j_y_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b1, k_4x4_0, tmp_b3, _MM_SWIZ_REG_BADC));
            p_j_z_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b2, k_4x4_1, tmp_b0, _MM_SWIZ_REG_BADC));
            p_j_q_vec = _mm512_castsi512_ps(_mm512_mask_swizzle_epi32(tmp_b3, k_4x4_1, tmp_b1, _MM_SWIZ_REG_BADC));
          }
        #else
          p_i_x_vec = _mm512_set_1to16_ps(p_0_x[iTest_i]);
          p_i_y_vec = _mm512_set_1to16_ps(p_0_y[iTest_i]);
          p_i_z_vec = _mm512_set_1to16_ps(p_0_z[iTest_i]);
          p_i_q_vec = _mm512_set_1to16_ps(p_0_q[iTest_i]);
          #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
            pExt_i_vdwType_vec = _mm512_set_1to16_epi32(pExt_0_vdwType[iTest_i]);
          #endif
          #if MIC_HANDCODE_FORCE_USEGATHER != 0
            p_j_x_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_x, _MM_SCALE_4);
            p_j_y_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_y, _MM_SCALE_4);
            p_j_z_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_z, _MM_SCALE_4);
            p_j_q_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_q, _MM_SCALE_4);
            #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
              pExt_j_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), active_mask, j_vec, pExt_1_vdwType, _MM_SCALE_4);
            #endif
          #else
            GATHER_PS_I32(p_j_x_vec, p_1_x, tmpJ32);
            GATHER_PS_I32(p_j_y_vec, p_1_y, tmpJ32);
            GATHER_PS_I32(p_j_z_vec, p_1_z, tmpJ32);
            GATHER_PS_I32(p_j_q_vec, p_1_q, tmpJ32);
            #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
              GATHER_EPI32_I32(pExt_j_vdwType_vec, pExt_1_vdwType, tmpJ32);
            #endif
          #endif
        #endif // MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
      #if __MIC_PAD_PLGEN_CTRL != 0
      #else
      } else {
        #if MIC_HANDCODE_FORCE_USEGATHER != 0
          p_i_x_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, i_vec, p_0_x, _MM_SCALE_4);
          p_i_y_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, i_vec, p_0_y, _MM_SCALE_4);
          p_i_z_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, i_vec, p_0_z, _MM_SCALE_4);
          p_i_q_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, i_vec, p_0_q, _MM_SCALE_4);
          #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
            pExt_i_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), active_mask, i_vec, pExt_0_vdwType, _MM_SCALE_4);
          #endif
          p_j_x_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_x, _MM_SCALE_4);
          p_j_y_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_y, _MM_SCALE_4);
          p_j_z_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_z, _MM_SCALE_4);
          p_j_q_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), active_mask, j_vec, p_1_q, _MM_SCALE_4);
          #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
            pExt_j_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), active_mask, j_vec, pExt_1_vdwType, _MM_SCALE_4);
          #endif
        #else
          GATHER_PS_I32(p_i_x_vec, p_0_x, tmpI32);
          GATHER_PS_I32(p_i_y_vec, p_0_y, tmpI32);
          GATHER_PS_I32(p_i_z_vec, p_0_z, tmpI32);
          GATHER_PS_I32(p_i_q_vec, p_0_q, tmpI32);
          #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
            GATHER_EPI32_I32(pExt_i_vdwType_vec, pExt_0_vdwType, tmpI32);
          #endif
         GATHER_PS_I32(p_j_x_vec, p_1_x, tmpJ32);
          GATHER_PS_I32(p_j_y_vec, p_1_y, tmpJ32);
          GATHER_PS_I32(p_j_z_vec, p_1_z, tmpJ32);
          GATHER_PS_I32(p_j_q_vec, p_1_q, tmpJ32);
          #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT != 0
            GATHER_EPI32_I32(pExt_j_vdwType_vec, pExt_1_vdwType, tmpJ32);
          #endif
        #endif
      }
      #endif
    }

    // Calculate the delta-x, delta-y, delta-z, and radius-squared
    p_i_x_vec = _mm512_add_ps(p_i_x_vec, _mm512_set_1to16_ps(offset_x));
    p_i_y_vec = _mm512_add_ps(p_i_y_vec, _mm512_set_1to16_ps(offset_y));
    p_i_z_vec = _mm512_add_ps(p_i_z_vec, _mm512_set_1to16_ps(offset_z));
    __m512 p_ij_x_vec = _mm512_sub_ps(p_i_x_vec, p_j_x_vec);
    __m512 p_ij_y_vec = _mm512_sub_ps(p_i_y_vec, p_j_y_vec);
    __m512 p_ij_z_vec = _mm512_sub_ps(p_i_z_vec, p_j_z_vec);
    __m512 r2_vec = _mm512_add_ps(_mm512_mul_ps(p_ij_x_vec, p_ij_x_vec), _mm512_set_1to16_ps(r2_delta));
    r2_vec = _mm512_add_ps(_mm512_mul_ps(p_ij_y_vec, p_ij_y_vec), r2_vec);
    r2_vec = _mm512_add_ps(_mm512_mul_ps(p_ij_z_vec, p_ij_z_vec), r2_vec);

    // Do the 'within cutoff' check
    // DMK - NOTE - When intrinsics are used, the compiler generates extra instructions to clear all but the least significant
    //   8 bits of the mask register storing the result of the comparison.
    __mmask16 cutoff_mask = _mm512_mask_cmplt_ps_mask(active_mask, r2_vec, _mm512_set_1to16_ps(cutoff2_delta));
    if (_mm512_kortestz(cutoff_mask, cutoff_mask)) { continue; }  // If the mask is completely unset, move on to the next vector

    #if (MIC_EXCL_CHECKSUM_FULL != 0) && (0 EXCLUDED(+1) MODIFIED(+1))
      exclusionSum_vec = _mm512_mask_add_epi32(exclusionSum_vec, cutoff_mask, exclusionSum_vec, _mm512_set_1to16_epi32(1));
    #endif

    // Calculate kqq = p_0_q[i] * p_1_q[j]
    __m512 kqq_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, p_i_q_vec, p_j_q_vec);

    // Calculate table_i = (r2 >> 46) + r2_delta_expc
    __m512i r2i_vec = _mm512_castps_si512(r2_vec);
    __m512i table_i_vec = _mm512_srli_epi32(r2i_vec, 17);
    table_i_vec = _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_i_vec, _mm512_set_1to16_epi32(r2_delta_expc));

    // Calculate diffa = r2 - r2_table[table_i]
    __m512 r2_table_vec;
    {
      // From ComputeNonbondedUtil.C                    Simplified:
      //   r2_base = r2_delta * (1 << (i/64))             r2_base = r2_delta * (1 << (i/64))
      //   r2_del = r2_base / 64.0;                       r2_del = r2_base / 64.0;
      //   r2 = r2_base - r2_delta + r2_del * (i%64)      r2_table[i] = r2_base - r2_delta + r2_del * (i%64) + r2_delta;
      //   r2_table[i] = r2 + r2_delta;                               = r2_base + r2_del * (i%64)
      // NOTE: For i = 0, r2_table[0] = r2_delta + (r2_delta / 64) * 0 = r2_delta, so there no need
      //   to special case if table_i = 0 then r2_table[0] = r2_delta (see ComputeNonbondedUtil.C:606)
      __m512 r2_delta_vec = _mm512_set_1to16_ps(r2_delta);
      __m512 t0_vec = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_sllv_epi32(_mm512_set_1to16_epi32(1), _mm512_srli_epi32(table_i_vec, 6)), _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
      __m512 r2_base_vec = _mm512_mul_ps(r2_delta_vec, t0_vec); // NOTE: r2_delta * (1 << (i/64))
      __m512 r2_del_vec = _mm512_mul_ps(r2_base_vec, _mm512_set_1to16_ps(0.015625f)); // NOTE: r2_base / 64
      __m512i t1_vec = _mm512_and_epi32(table_i_vec, _mm512_set_1to16_epi32(0x3F)); // NOTE: (i%64)
      __m512 t2_vec = _mm512_cvtfxpnt_round_adjustepi32_ps(t1_vec, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE); // NOTE: (float)(i%64)
      r2_table_vec = _mm512_add_ps(_mm512_mul_ps(r2_del_vec, t2_vec), r2_base_vec);
    }
    __m512 diffa_vec = _mm512_sub_ps(r2_vec, r2_table_vec);

    // Load LJ A and B values
    #if (0 FAST(+1))

      // Load lj_pars from the table and multiply by scaling
      __m512 A_vec, B_vec;
          __m512i pExt_i_vdwType_vec;  // NOTE : Moved outside block so could be printed in large force detection code below
          __m512i pExt_j_vdwType_vec;
      {
        // Load VDW types for the "i" and "j" atoms
        #if MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT == 0
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            // DMK - NOTE : This code assumes that vdw_type is the first of four members of atom_param
            pExt_i_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), cutoff_mask, _mm512_slli_epi32(i_vec, 2), pExt_0, _MM_SCALE_4);
	    pExt_j_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), cutoff_mask, _mm512_slli_epi32(j_vec, 2), pExt_1, _MM_SCALE_4);
          #else
            #if MIC_HANDCODE_FORCE_USEGATHER != 0
              pExt_i_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), cutoff_mask, i_vec, pExt_0_vdwType, _MM_SCALE_4);
              pExt_j_vdwType_vec = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), cutoff_mask, j_vec, pExt_1_vdwType, _MM_SCALE_4);
            #else
              GATHER_EPI32_I32(pExt_i_vdwType_vec, pExt_0_vdwType, tmpI32);
              GATHER_EPI32_I32(pExt_j_vdwType_vec, pExt_1_vdwType, tmpJ32);
            #endif
          #endif
        #endif

        // Caclulate offsets into LJ table
        __m512i t0_vec = _mm512_fmadd_epi32(pExt_i_vdwType_vec, _mm512_set_1to16_epi32(lj_table_dim), pExt_j_vdwType_vec);
        __m512i ljOffset_vec = _mm512_slli_epi32(t0_vec, 2); // NOTE: *4
        MODIFIED( ljOffset_vec = _mm512_add_epi32(ljOffset_vec, _mm512_set_1to16_epi32(2)); )

        // Gather A and B values
        #if MIC_HANDCODE_FORCE_USEGATHER != 0
          __m512i ljOffset_p1_vec = _mm512_add_epi32(ljOffset_vec, _mm512_set_1to16_epi32(1));
          A_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask,    ljOffset_vec, lj_table_base_ptr, _MM_SCALE_4);
	  B_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, ljOffset_p1_vec, lj_table_base_ptr, _MM_SCALE_4);
        #else
          unsigned int ljOffset[16] __attribute__((aligned(64)));
          _mm512_store_epi32(ljOffset, ljOffset_vec);
          GATHER_PS_I32       (A_vec, lj_table_base_ptr, ljOffset   );
          GATHER_PS_I32_OFFSET(B_vec, lj_table_base_ptr, ljOffset, 1);
        #endif

        // Scale the A and B values
        __m512 scaling_vec = _mm512_set_1to16_ps(scaling);
        A_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, A_vec, scaling_vec);
        B_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, B_vec, scaling_vec);
      }
    #endif // FAST

    // Load the table_four values
    // NOTE : In this loop, each vector lane is calculating an interaction between two particles, i and j.
    //   As part of this calculation, a table lookup is done (this code) where each table entry has 16
    //   values. Rather than do 16 gathers, extract the index for each interaction, do vector loads on
    //   the table entries (i.e. 1 table entry = 16 values = 2 512-bit vector registers), and then transpose
    //   those values in-register to create the vectorized registers that can be used by the code below
    //   (i.e. 1 vector register per each of the 16 fields).  Once the values have been loaded, only
    //   swizzles and permutes are required to do this, so it shouldn't be very expensive to do.
    const CALC_TYPE * const table_four_base = SHORT(table_short) NOSHORT(table_noshort);
    __m512 table_four_i_0_vec,  table_four_i_1_vec,  table_four_i_2_vec,  table_four_i_3_vec,
           table_four_i_4_vec,  table_four_i_5_vec,  table_four_i_6_vec,  table_four_i_7_vec,
           table_four_i_8_vec,  table_four_i_9_vec,  table_four_i_10_vec, table_four_i_11_vec,
           table_four_i_12_vec, table_four_i_13_vec, table_four_i_14_vec, table_four_i_15_vec;
    {
      #if MIC_HANDCODE_FORCE_USEGATHER_NBTBL != 0

        #if 0

          __m512i table_four_i_offsets = _mm512_slli_epi32(table_i_vec, 4); // NOTE: table_i * 16 (16 elements per table entry)
          table_four_i_0_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask,                                                            table_four_i_offsets                             , table_four_base, _MM_SCALE_4);
          table_four_i_1_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 1)), table_four_base, _MM_SCALE_4);
          table_four_i_2_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 2)), table_four_base, _MM_SCALE_4);
          table_four_i_3_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 3)), table_four_base, _MM_SCALE_4);
          table_four_i_4_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 4)), table_four_base, _MM_SCALE_4);
          table_four_i_5_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 5)), table_four_base, _MM_SCALE_4);
          table_four_i_6_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 6)), table_four_base, _MM_SCALE_4);
          table_four_i_7_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 7)), table_four_base, _MM_SCALE_4);
          table_four_i_8_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 8)), table_four_base, _MM_SCALE_4);
          table_four_i_9_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32( 9)), table_four_base, _MM_SCALE_4);
          table_four_i_10_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32(10)), table_four_base, _MM_SCALE_4);
          table_four_i_11_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32(11)), table_four_base, _MM_SCALE_4);
          table_four_i_12_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32(12)), table_four_base, _MM_SCALE_4);
          table_four_i_13_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32(13)), table_four_base, _MM_SCALE_4);
          table_four_i_14_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32(14)), table_four_base, _MM_SCALE_4);
          table_four_i_15_vec = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, _mm512_mask_add_epi32(_mm512_setzero_epi32(), cutoff_mask, table_four_i_offsets, _mm512_set_1to16_epi32(15)), table_four_base, _MM_SCALE_4);

        #else

          __m512i table_four_i_offsets = _mm512_slli_epi32(table_i_vec, 4); // NOTE: table_i * 16 (16 elements per table entry)
          table_four_i_0_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base, _MM_SCALE_4);
          table_four_i_1_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 1, _MM_SCALE_4);
          table_four_i_2_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 2, _MM_SCALE_4);
          table_four_i_3_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 3, _MM_SCALE_4);
          table_four_i_4_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 4, _MM_SCALE_4);
          table_four_i_5_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 5, _MM_SCALE_4);
          table_four_i_6_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 6, _MM_SCALE_4);
          table_four_i_7_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 7, _MM_SCALE_4);
          table_four_i_8_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 8, _MM_SCALE_4);
          table_four_i_9_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 9, _MM_SCALE_4);
          table_four_i_10_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 10, _MM_SCALE_4);
          table_four_i_11_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 11, _MM_SCALE_4);
          table_four_i_12_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 12, _MM_SCALE_4);
          table_four_i_13_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 13, _MM_SCALE_4);
          table_four_i_14_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 14, _MM_SCALE_4);
          table_four_i_15_vec  = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), cutoff_mask, table_four_i_offsets, table_four_base + 15, _MM_SCALE_4);

        #endif

      #else

      unsigned int table_four_i_offsets[16] __attribute__((aligned(64)));
      _mm512_store_epi32(table_four_i_offsets, _mm512_slli_epi32(table_i_vec, 4)); // NOTE: table_i * 16 (16 elements per table entry)

      // Load and transpose 256 values (16 x 16 matrix)
      __m512i tmpA_00 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 0]])));
      __m512i tmpA_01 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 1]])));
      __m512i tmpA_02 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 2]])));
      __m512i tmpA_03 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 3]])));
      __m512i tmpA_04 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 4]])));
      __m512i tmpA_05 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 5]])));
      __m512i tmpA_06 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 6]])));
      __m512i tmpA_07 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 7]])));
      __m512i tmpA_08 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 8]])));
      __m512i tmpA_09 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[ 9]])));
      __m512i tmpA_10 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[10]])));
      __m512i tmpA_11 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[11]])));
      __m512i tmpA_12 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[12]])));
      __m512i tmpA_13 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[13]])));
      __m512i tmpA_14 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[14]])));
      __m512i tmpA_15 = _mm512_castps_si512(_mm512_load_ps(&(table_four_base[table_four_i_offsets[15]])));
      __mmask16 k_2x2_0 = _mm512_int2mask(0xAAAA);
      __mmask16 k_2x2_1 = _mm512_int2mask(0x5555);
      __m512i tmpB_00 = _mm512_mask_swizzle_epi32(tmpA_00, k_2x2_0, tmpA_01, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_01 = _mm512_mask_swizzle_epi32(tmpA_01, k_2x2_1, tmpA_00, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_02 = _mm512_mask_swizzle_epi32(tmpA_02, k_2x2_0, tmpA_03, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_03 = _mm512_mask_swizzle_epi32(tmpA_03, k_2x2_1, tmpA_02, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_04 = _mm512_mask_swizzle_epi32(tmpA_04, k_2x2_0, tmpA_05, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_05 = _mm512_mask_swizzle_epi32(tmpA_05, k_2x2_1, tmpA_04, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_06 = _mm512_mask_swizzle_epi32(tmpA_06, k_2x2_0, tmpA_07, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_07 = _mm512_mask_swizzle_epi32(tmpA_07, k_2x2_1, tmpA_06, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_08 = _mm512_mask_swizzle_epi32(tmpA_08, k_2x2_0, tmpA_09, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_09 = _mm512_mask_swizzle_epi32(tmpA_09, k_2x2_1, tmpA_08, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_10 = _mm512_mask_swizzle_epi32(tmpA_10, k_2x2_0, tmpA_11, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_11 = _mm512_mask_swizzle_epi32(tmpA_11, k_2x2_1, tmpA_10, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_12 = _mm512_mask_swizzle_epi32(tmpA_12, k_2x2_0, tmpA_13, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_13 = _mm512_mask_swizzle_epi32(tmpA_13, k_2x2_1, tmpA_12, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_14 = _mm512_mask_swizzle_epi32(tmpA_14, k_2x2_0, tmpA_15, _MM_SWIZ_REG_CDAB);
      __m512i tmpB_15 = _mm512_mask_swizzle_epi32(tmpA_15, k_2x2_1, tmpA_14, _MM_SWIZ_REG_CDAB);
      __mmask16 k_4x4_0 = _mm512_int2mask(0xCCCC);
      __mmask16 k_4x4_1 = _mm512_int2mask(0x3333);
      __m512i tmpC_00 = _mm512_mask_swizzle_epi32(tmpB_00, k_4x4_0, tmpB_02, _MM_SWIZ_REG_BADC);
      __m512i tmpC_01 = _mm512_mask_swizzle_epi32(tmpB_01, k_4x4_0, tmpB_03, _MM_SWIZ_REG_BADC);
      __m512i tmpC_02 = _mm512_mask_swizzle_epi32(tmpB_02, k_4x4_1, tmpB_00, _MM_SWIZ_REG_BADC);
      __m512i tmpC_03 = _mm512_mask_swizzle_epi32(tmpB_03, k_4x4_1, tmpB_01, _MM_SWIZ_REG_BADC);
      __m512i tmpC_04 = _mm512_mask_swizzle_epi32(tmpB_04, k_4x4_0, tmpB_06, _MM_SWIZ_REG_BADC);
      __m512i tmpC_05 = _mm512_mask_swizzle_epi32(tmpB_05, k_4x4_0, tmpB_07, _MM_SWIZ_REG_BADC);
      __m512i tmpC_06 = _mm512_mask_swizzle_epi32(tmpB_06, k_4x4_1, tmpB_04, _MM_SWIZ_REG_BADC);
      __m512i tmpC_07 = _mm512_mask_swizzle_epi32(tmpB_07, k_4x4_1, tmpB_05, _MM_SWIZ_REG_BADC);
      __m512i tmpC_08 = _mm512_mask_swizzle_epi32(tmpB_08, k_4x4_0, tmpB_10, _MM_SWIZ_REG_BADC);
      __m512i tmpC_09 = _mm512_mask_swizzle_epi32(tmpB_09, k_4x4_0, tmpB_11, _MM_SWIZ_REG_BADC);
      __m512i tmpC_10 = _mm512_mask_swizzle_epi32(tmpB_10, k_4x4_1, tmpB_08, _MM_SWIZ_REG_BADC);
      __m512i tmpC_11 = _mm512_mask_swizzle_epi32(tmpB_11, k_4x4_1, tmpB_09, _MM_SWIZ_REG_BADC);
      __m512i tmpC_12 = _mm512_mask_swizzle_epi32(tmpB_12, k_4x4_0, tmpB_14, _MM_SWIZ_REG_BADC);
      __m512i tmpC_13 = _mm512_mask_swizzle_epi32(tmpB_13, k_4x4_0, tmpB_15, _MM_SWIZ_REG_BADC);
      __m512i tmpC_14 = _mm512_mask_swizzle_epi32(tmpB_14, k_4x4_1, tmpB_12, _MM_SWIZ_REG_BADC);
      __m512i tmpC_15 = _mm512_mask_swizzle_epi32(tmpB_15, k_4x4_1, tmpB_13, _MM_SWIZ_REG_BADC);
      __mmask16 k_8x8_0 = _mm512_int2mask(0xF0F0);
      __mmask16 k_8x8_1 = _mm512_int2mask(0x0F0F);
      __m512i tmpD_00 = _mm512_mask_permute4f128_epi32(tmpC_00, k_8x8_0, tmpC_04, _MM_PERM_CDAB);
      __m512i tmpD_01 = _mm512_mask_permute4f128_epi32(tmpC_01, k_8x8_0, tmpC_05, _MM_PERM_CDAB);
      __m512i tmpD_02 = _mm512_mask_permute4f128_epi32(tmpC_02, k_8x8_0, tmpC_06, _MM_PERM_CDAB);
      __m512i tmpD_03 = _mm512_mask_permute4f128_epi32(tmpC_03, k_8x8_0, tmpC_07, _MM_PERM_CDAB);
      __m512i tmpD_04 = _mm512_mask_permute4f128_epi32(tmpC_04, k_8x8_1, tmpC_00, _MM_PERM_CDAB);
      __m512i tmpD_05 = _mm512_mask_permute4f128_epi32(tmpC_05, k_8x8_1, tmpC_01, _MM_PERM_CDAB);
      __m512i tmpD_06 = _mm512_mask_permute4f128_epi32(tmpC_06, k_8x8_1, tmpC_02, _MM_PERM_CDAB);
      __m512i tmpD_07 = _mm512_mask_permute4f128_epi32(tmpC_07, k_8x8_1, tmpC_03, _MM_PERM_CDAB);
      __m512i tmpD_08 = _mm512_mask_permute4f128_epi32(tmpC_08, k_8x8_0, tmpC_12, _MM_PERM_CDAB);
      __m512i tmpD_09 = _mm512_mask_permute4f128_epi32(tmpC_09, k_8x8_0, tmpC_13, _MM_PERM_CDAB);
      __m512i tmpD_10 = _mm512_mask_permute4f128_epi32(tmpC_10, k_8x8_0, tmpC_14, _MM_PERM_CDAB);
      __m512i tmpD_11 = _mm512_mask_permute4f128_epi32(tmpC_11, k_8x8_0, tmpC_15, _MM_PERM_CDAB);
      __m512i tmpD_12 = _mm512_mask_permute4f128_epi32(tmpC_12, k_8x8_1, tmpC_08, _MM_PERM_CDAB);
      __m512i tmpD_13 = _mm512_mask_permute4f128_epi32(tmpC_13, k_8x8_1, tmpC_09, _MM_PERM_CDAB);
      __m512i tmpD_14 = _mm512_mask_permute4f128_epi32(tmpC_14, k_8x8_1, tmpC_10, _MM_PERM_CDAB);
      __m512i tmpD_15 = _mm512_mask_permute4f128_epi32(tmpC_15, k_8x8_1, tmpC_11, _MM_PERM_CDAB);
      __mmask16 k_16x16_0 = _mm512_int2mask(0xFF00);
      __mmask16 k_16x16_1 = _mm512_int2mask(0x00FF);
      table_four_i_0_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_00, k_16x16_0, tmpD_08, _MM_PERM_BADC));
      table_four_i_1_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_01, k_16x16_0, tmpD_09, _MM_PERM_BADC));
      table_four_i_2_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_02, k_16x16_0, tmpD_10, _MM_PERM_BADC));
      table_four_i_3_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_03, k_16x16_0, tmpD_11, _MM_PERM_BADC));
      table_four_i_4_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_04, k_16x16_0, tmpD_12, _MM_PERM_BADC));
      table_four_i_5_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_05, k_16x16_0, tmpD_13, _MM_PERM_BADC));
      table_four_i_6_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_06, k_16x16_0, tmpD_14, _MM_PERM_BADC));
      table_four_i_7_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_07, k_16x16_0, tmpD_15, _MM_PERM_BADC));
      table_four_i_8_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_08, k_16x16_1, tmpD_00, _MM_PERM_BADC));
      table_four_i_9_vec  = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_09, k_16x16_1, tmpD_01, _MM_PERM_BADC));
      table_four_i_10_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_10, k_16x16_1, tmpD_02, _MM_PERM_BADC));
      table_four_i_11_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_11, k_16x16_1, tmpD_03, _MM_PERM_BADC));
      table_four_i_12_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_12, k_16x16_1, tmpD_04, _MM_PERM_BADC));
      table_four_i_13_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_13, k_16x16_1, tmpD_05, _MM_PERM_BADC));
      table_four_i_14_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_14, k_16x16_1, tmpD_06, _MM_PERM_BADC));
      table_four_i_15_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(tmpD_15, k_16x16_1, tmpD_07, _MM_PERM_BADC));

      #endif
    }

    #if (0 FAST(+1))

      __m512 vdw_d_vec = _mm512_sub_ps(_mm512_mul_ps(A_vec, table_four_i_0_vec), _mm512_mul_ps(B_vec, table_four_i_4_vec));
      __m512 vdw_c_vec = _mm512_sub_ps(_mm512_mul_ps(A_vec, table_four_i_1_vec), _mm512_mul_ps(B_vec, table_four_i_5_vec));
      __m512 vdw_b_vec = _mm512_sub_ps(_mm512_mul_ps(A_vec, table_four_i_2_vec), _mm512_mul_ps(B_vec, table_four_i_6_vec));
      __m512 vdw_a_vec = _mm512_sub_ps(_mm512_mul_ps(A_vec, table_four_i_3_vec), _mm512_mul_ps(B_vec, table_four_i_7_vec));

      #if (0 ENERGY(+1))
        __m512 vdw_val_tmp_0_vec = _mm512_mul_ps(vdw_d_vec, _mm512_set_1to16_ps(1.0f/6.0f));
        __m512 vdw_val_tmp_1_vec = _mm512_mul_ps(vdw_c_vec, _mm512_set_1to16_ps(1.0f/4.0f));
        __m512 vdw_val_tmp_2_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, vdw_val_tmp_0_vec), vdw_val_tmp_1_vec);
        __m512 vdw_val_tmp_3_vec = _mm512_mul_ps(vdw_b_vec, _mm512_set_1to16_ps(1.0f/2.0f));
        __m512 vdw_val_tmp_4_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, vdw_val_tmp_2_vec), vdw_val_tmp_3_vec);
        __m512 vdw_val_vec       = _mm512_add_ps(_mm512_mul_ps(diffa_vec, vdw_val_tmp_4_vec), vdw_a_vec);
        CONTRIB_SUB_PS2PD(vdwEnergy_vec, vdw_val_vec);
      #endif

      #if (0 SHORT(+1))

        #if (0 NORMAL(+1))
          __m512 fast_d_vec = _mm512_mul_ps(kqq_vec, table_four_i_8_vec);
          __m512 fast_c_vec = _mm512_mul_ps(kqq_vec, table_four_i_9_vec);
          __m512 fast_b_vec = _mm512_mul_ps(kqq_vec, table_four_i_10_vec);
          __m512 fast_a_vec = _mm512_mul_ps(kqq_vec, table_four_i_11_vec);
        #endif
        #if (0 MODIFIED(+1))
          __m512 modfckqq_vec = _mm512_mul_ps(_mm512_set_1to16_ps(1.0f - modf_mod), kqq_vec);
          __m512 fast_d_vec = _mm512_mul_ps(modfckqq_vec, table_four_i_8_vec);
          __m512 fast_c_vec = _mm512_mul_ps(modfckqq_vec, table_four_i_9_vec);
          __m512 fast_b_vec = _mm512_mul_ps(modfckqq_vec, table_four_i_10_vec);
          __m512 fast_a_vec = _mm512_mul_ps(modfckqq_vec, table_four_i_11_vec);
        #endif

        #if (0 ENERGY(+1))
          __m512 fast_val_tmp_0_vec = _mm512_mul_ps(fast_d_vec, _mm512_set_1to16_ps(1.0f/6.0f));
          __m512 fast_val_tmp_1_vec = _mm512_mul_ps(fast_c_vec, _mm512_set_1to16_ps(1.0f/4.0f));
          __m512 fast_val_tmp_2_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, fast_val_tmp_0_vec), fast_val_tmp_1_vec);
          __m512 fast_val_tmp_3_vec = _mm512_mul_ps(fast_b_vec, _mm512_set_1to16_ps(1.0f/2.0f));
          __m512 fast_val_tmp_4_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, fast_val_tmp_2_vec), fast_val_tmp_3_vec);
          __m512 fast_val_vec       = _mm512_add_ps(_mm512_mul_ps(diffa_vec, fast_val_tmp_4_vec), fast_a_vec);
          CONTRIB_SUB_PS2PD(electEnergy_vec, fast_val_vec);
        #endif

        #if (0 NOT_ALCHPAIR(+1))
          fast_d_vec = _mm512_add_ps(fast_d_vec, vdw_d_vec);
          fast_c_vec = _mm512_add_ps(fast_c_vec, vdw_c_vec);
          fast_b_vec = _mm512_add_ps(fast_b_vec, vdw_b_vec);
          fast_a_vec = _mm512_add_ps(fast_a_vec, vdw_a_vec);
        #endif

	__m512 fast_dir_tmp_0_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, fast_d_vec), fast_c_vec);
        __m512 fast_dir_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, fast_dir_tmp_0_vec), fast_b_vec);
        __m512 force_r_vec = fast_dir_vec;  // NOTE: No-op left in as placeholder for when LAM is added

        const __m512 tmp_x_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, force_r_vec, p_ij_x_vec);
        PAIR( CONTRIB_ADD_PS2PD(virial_xx_vec, _mm512_mul_ps(tmp_x_vec, p_ij_x_vec)); )
        PAIR( CONTRIB_ADD_PS2PD(virial_xy_vec, _mm512_mul_ps(tmp_x_vec, p_ij_y_vec)); )
        PAIR( CONTRIB_ADD_PS2PD(virial_xz_vec, _mm512_mul_ps(tmp_x_vec, p_ij_z_vec)); )

        const __m512 tmp_y_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, force_r_vec, p_ij_y_vec);
        PAIR( CONTRIB_ADD_PS2PD(virial_yy_vec, _mm512_mul_ps(tmp_y_vec, p_ij_y_vec)); )
        PAIR( CONTRIB_ADD_PS2PD(virial_yz_vec, _mm512_mul_ps(tmp_y_vec, p_ij_z_vec)); )

        const __m512 tmp_z_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, force_r_vec, p_ij_z_vec);
        PAIR( CONTRIB_ADD_PS2PD(virial_zz_vec, _mm512_mul_ps(tmp_z_vec, p_ij_z_vec)); )

        #if MIC_HANDCODE_FORCE_COMBINE_FORCES == 0
          APPLY_FORCES_PS2PD(tmp_x_vec, tmp_y_vec, tmp_z_vec, f_0_x, f_0_y, f_0_z, f_1_x, f_1_y, f_1_z, tmpI32, tmpJ32);
        #endif

      #endif // SHORT
    #endif // FAST

    #if (0 FULL(+1))

      // Load the slow_table values, if need be
      #if (0 SHORT( EXCLUDED(+1) MODIFIED(+1) ))

        __m512 slow_i_0_vec, slow_i_1_vec, slow_i_2_vec, slow_i_3_vec;
	{
          // Create the indexes
          unsigned int slow_i_offsets[16] __attribute__((aligned(64)));
          _mm512_store_epi32(slow_i_offsets, _mm512_slli_epi32(table_i_vec, 2)); // table_i * 4

          // Load the values (64 values, 16x4 matrix, so we only care about elements 0-3)
          // DMK - NOTE : Since the slow table is aligned, all four elements being loaded by the loadunpacklo
          //   will be grouped on a single cacheline, so no need to also include a loadunpackhi instruction
          #define __LOAD_SLOW(dstA, offset) \
            __m512i (dstA) = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0x000F), slow_table + (offset));
          __LOAD_SLOW(slow_i_tmp00_vec, slow_i_offsets[ 0]);
          __LOAD_SLOW(slow_i_tmp01_vec, slow_i_offsets[ 1]);
          __LOAD_SLOW(slow_i_tmp02_vec, slow_i_offsets[ 2]);
          __LOAD_SLOW(slow_i_tmp03_vec, slow_i_offsets[ 3]);
          __LOAD_SLOW(slow_i_tmp04_vec, slow_i_offsets[ 4]);
          __LOAD_SLOW(slow_i_tmp05_vec, slow_i_offsets[ 5]);
          __LOAD_SLOW(slow_i_tmp06_vec, slow_i_offsets[ 6]);
          __LOAD_SLOW(slow_i_tmp07_vec, slow_i_offsets[ 7]);
          __LOAD_SLOW(slow_i_tmp08_vec, slow_i_offsets[ 8]);
          __LOAD_SLOW(slow_i_tmp09_vec, slow_i_offsets[ 9]);
          __LOAD_SLOW(slow_i_tmp10_vec, slow_i_offsets[10]);
          __LOAD_SLOW(slow_i_tmp11_vec, slow_i_offsets[11]);
          __LOAD_SLOW(slow_i_tmp12_vec, slow_i_offsets[12]);
          __LOAD_SLOW(slow_i_tmp13_vec, slow_i_offsets[13]);
          __LOAD_SLOW(slow_i_tmp14_vec, slow_i_offsets[14]);
          __LOAD_SLOW(slow_i_tmp15_vec, slow_i_offsets[15]);
          #undef __LOAD_SLOW

          // Transpose the values
          // NOTE : Transpose each 2x2 sub-matrix, but only elements 0-4 have data we care about
          #define __TRANSPOSE_2x2__(dstA, dstB, srcA, srcB) \
            __m512i dstA = _mm512_mask_swizzle_epi32((srcA), _mm512_int2mask(0xAAAA), (srcB), _MM_SWIZ_REG_CDAB); \
            __m512i dstB = _mm512_mask_swizzle_epi32((srcB), _mm512_int2mask(0x5555), (srcA), _MM_SWIZ_REG_CDAB);
          __TRANSPOSE_2x2__(slow_2x2_tmp00_vec, slow_2x2_tmp01_vec, slow_i_tmp00_vec, slow_i_tmp01_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp02_vec, slow_2x2_tmp03_vec, slow_i_tmp02_vec, slow_i_tmp03_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp04_vec, slow_2x2_tmp05_vec, slow_i_tmp04_vec, slow_i_tmp05_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp06_vec, slow_2x2_tmp07_vec, slow_i_tmp06_vec, slow_i_tmp07_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp08_vec, slow_2x2_tmp09_vec, slow_i_tmp08_vec, slow_i_tmp09_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp10_vec, slow_2x2_tmp11_vec, slow_i_tmp10_vec, slow_i_tmp11_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp12_vec, slow_2x2_tmp13_vec, slow_i_tmp12_vec, slow_i_tmp13_vec);
          __TRANSPOSE_2x2__(slow_2x2_tmp14_vec, slow_2x2_tmp15_vec, slow_i_tmp14_vec, slow_i_tmp15_vec);
          #undef __TRANSPOSE_2x2__
          #define __TRANSPOSE_4x4__(dstA, dstB, dstC, dstD, srcA, srcB, srcC, srcD) \
            __m512i dstA = _mm512_mask_swizzle_epi32((srcA), _mm512_int2mask(0xCCCC), (srcC), _MM_SWIZ_REG_BADC); \
            __m512i dstB = _mm512_mask_swizzle_epi32((srcB), _mm512_int2mask(0xCCCC), (srcD), _MM_SWIZ_REG_BADC); \
            __m512i dstC = _mm512_mask_swizzle_epi32((srcC), _mm512_int2mask(0x3333), (srcA), _MM_SWIZ_REG_BADC); \
            __m512i dstD = _mm512_mask_swizzle_epi32((srcD), _mm512_int2mask(0x3333), (srcB), _MM_SWIZ_REG_BADC);
          __TRANSPOSE_4x4__(slow_4x4_tmp00_vec, slow_4x4_tmp01_vec, slow_4x4_tmp02_vec, slow_4x4_tmp03_vec, slow_2x2_tmp00_vec, slow_2x2_tmp01_vec, slow_2x2_tmp02_vec, slow_2x2_tmp03_vec);
          __TRANSPOSE_4x4__(slow_4x4_tmp04_vec, slow_4x4_tmp05_vec, slow_4x4_tmp06_vec, slow_4x4_tmp07_vec, slow_2x2_tmp04_vec, slow_2x2_tmp05_vec, slow_2x2_tmp06_vec, slow_2x2_tmp07_vec);
          __TRANSPOSE_4x4__(slow_4x4_tmp08_vec, slow_4x4_tmp09_vec, slow_4x4_tmp10_vec, slow_4x4_tmp11_vec, slow_2x2_tmp08_vec, slow_2x2_tmp09_vec, slow_2x2_tmp10_vec, slow_2x2_tmp11_vec);
          __TRANSPOSE_4x4__(slow_4x4_tmp12_vec, slow_4x4_tmp13_vec, slow_4x4_tmp14_vec, slow_4x4_tmp15_vec, slow_2x2_tmp12_vec, slow_2x2_tmp13_vec, slow_2x2_tmp14_vec, slow_2x2_tmp15_vec);
          #undef __TRANSPOSE_2x2__
          __m512i slow_4x8_tmp00_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp00_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp04_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp01_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp01_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp05_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp02_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp02_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp06_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp03_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp03_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp07_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp04_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp08_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp12_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp05_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp09_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp13_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp06_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp10_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp14_vec, _MM_PERM_CDAB);
          __m512i slow_4x8_tmp07_vec = _mm512_mask_permute4f128_epi32(slow_4x4_tmp11_vec, _mm512_int2mask(0x00F0), slow_4x4_tmp15_vec, _MM_PERM_CDAB);
          slow_i_0_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(slow_4x8_tmp00_vec, _mm512_int2mask(0xFF00), slow_4x8_tmp04_vec, _MM_PERM_BADC));
          slow_i_1_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(slow_4x8_tmp01_vec, _mm512_int2mask(0xFF00), slow_4x8_tmp05_vec, _MM_PERM_BADC));
          slow_i_2_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(slow_4x8_tmp02_vec, _mm512_int2mask(0xFF00), slow_4x8_tmp06_vec, _MM_PERM_BADC));
          slow_i_3_vec = _mm512_castsi512_ps(_mm512_mask_permute4f128_epi32(slow_4x8_tmp03_vec, _mm512_int2mask(0xFF00), slow_4x8_tmp07_vec, _MM_PERM_BADC));
	}

      #endif // SHORT && (EXCLUDED || MODIFIED)

      __m512 slow_d_vec = NOSHORT( table_four_i_8_vec) SHORT(table_four_i_12_vec);
      __m512 slow_c_vec = NOSHORT( table_four_i_9_vec) SHORT(table_four_i_13_vec);
      __m512 slow_b_vec = NOSHORT(table_four_i_10_vec) SHORT(table_four_i_14_vec);
      __m512 slow_a_vec = NOSHORT(table_four_i_11_vec) SHORT(table_four_i_15_vec);

      #if (0 EXCLUDED(+1))
        #if (0 SHORT(+1))
          slow_a_vec = _mm512_add_ps(                                         slow_i_3_vec , slow_a_vec);
          slow_b_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(2.0f), slow_i_2_vec), slow_b_vec);
          slow_c_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(4.0f), slow_i_1_vec), slow_c_vec);
          slow_d_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(6.0f), slow_i_0_vec), slow_d_vec);
        #endif
        #if (0 NOSHORT(+1))
          slow_d_vec = _mm512_sub_ps(slow_d_vec, table_four_i_12_vec);
          slow_c_vec = _mm512_sub_ps(slow_c_vec, table_four_i_13_vec);
          slow_b_vec = _mm512_sub_ps(slow_b_vec, table_four_i_14_vec);
          slow_a_vec = _mm512_sub_ps(slow_a_vec, table_four_i_15_vec);
        #endif
      #endif // EXCLUDED
      #if (0 MODIFIED(+1))
        #if (0 SHORT(+1))
          slow_a_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(       modf_mod), slow_i_3_vec), slow_a_vec);
          slow_b_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(2.0f * modf_mod), slow_i_2_vec), slow_b_vec);
          slow_c_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(4.0f * modf_mod), slow_i_1_vec), slow_c_vec);
          slow_d_vec = _mm512_add_ps(_mm512_mul_ps(_mm512_set_1to16_ps(6.0f * modf_mod), slow_i_0_vec), slow_d_vec);
        #endif
        #if (0 NOSHORT(+1))
          __m512 modf_mod_vec = _mm512_set_1to16_ps(modf_mod);
          slow_d_vec = _mm512_sub_ps(slow_d_vec, _mm512_mul_ps(modf_mod, table_four_i_12_vec));
          slow_c_vec = _mm512_sub_ps(slow_c_vec, _mm512_mul_ps(modf_mod, table_four_i_13_vec));
          slow_b_vec = _mm512_sub_ps(slow_b_vec, _mm512_mul_ps(modf_mod, table_four_i_14_vec));
          slow_a_vec = _mm512_sub_ps(slow_a_vec, _mm512_mul_ps(modf_mod, table_four_i_15_vec));
        #endif
      #endif // MODIFIED
      slow_d_vec = _mm512_mul_ps(slow_d_vec, kqq_vec);
      slow_c_vec = _mm512_mul_ps(slow_c_vec, kqq_vec);
      slow_b_vec = _mm512_mul_ps(slow_b_vec, kqq_vec);
      slow_a_vec = _mm512_mul_ps(slow_a_vec, kqq_vec);

      #if (0 ENERGY(+1))
        __m512 slow_val_tmp_0_vec = _mm512_mul_ps(slow_d_vec, _mm512_set_1to16_ps(1.0f/6.0f));
        __m512 slow_val_tmp_1_vec = _mm512_mul_ps(slow_c_vec, _mm512_set_1to16_ps(1.0f/4.0f));
        __m512 slow_val_tmp_2_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, slow_val_tmp_0_vec), slow_val_tmp_1_vec);
        __m512 slow_val_tmp_3_vec = _mm512_mul_ps(slow_b_vec, _mm512_set_1to16_ps(1.0f/2.0f));
        __m512 slow_val_tmp_4_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, slow_val_tmp_2_vec), slow_val_tmp_3_vec);
        __m512 slow_val_vec       = _mm512_add_ps(_mm512_mul_ps(diffa_vec, slow_val_tmp_4_vec), slow_a_vec);
        CONTRIB_SUB_PS2PD(fullElectEnergy_vec, slow_val_vec);
      #endif

      #if (0 NOT_ALCHPAIR(FAST(NOSHORT(+1))))
        slow_d_vec = _mm512_add_ps(vdw_d_vec, slow_d_vec);
        slow_c_vec = _mm512_add_ps(vdw_c_vec, slow_c_vec);
        slow_b_vec = _mm512_add_ps(vdw_b_vec, slow_b_vec);
        slow_a_vec = _mm512_add_ps(vdw_a_vec, slow_a_vec);
      #endif

      __m512 slow_dir_tmp_0_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, slow_d_vec), slow_c_vec);
      __m512 slow_dir_vec = _mm512_add_ps(_mm512_mul_ps(diffa_vec, slow_dir_tmp_0_vec), slow_b_vec);
      __m512 fullforce_r_vec = slow_dir_vec;  // NOTE: No-op left in as a placeholder for when LAM is added

      const __m512 fulltmp_x_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, fullforce_r_vec, p_ij_x_vec);
      PAIR( CONTRIB_ADD_PS2PD(fullElectVirial_xx_vec, _mm512_mul_ps(fulltmp_x_vec, p_ij_x_vec)); )
      PAIR( CONTRIB_ADD_PS2PD(fullElectVirial_xy_vec, _mm512_mul_ps(fulltmp_x_vec, p_ij_y_vec)); )
      PAIR( CONTRIB_ADD_PS2PD(fullElectVirial_xz_vec, _mm512_mul_ps(fulltmp_x_vec, p_ij_z_vec)); )

      const __m512 fulltmp_y_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, fullforce_r_vec, p_ij_y_vec);
      PAIR( CONTRIB_ADD_PS2PD(fullElectVirial_yy_vec, _mm512_mul_ps(fulltmp_y_vec, p_ij_y_vec)); )
      PAIR( CONTRIB_ADD_PS2PD(fullElectVirial_yz_vec, _mm512_mul_ps(fulltmp_y_vec, p_ij_z_vec)); )

      const __m512 fulltmp_z_vec = _mm512_mask_mul_ps(_mm512_setzero_ps(), cutoff_mask, fullforce_r_vec, p_ij_z_vec);
      PAIR( CONTRIB_ADD_PS2PD(fullElectVirial_zz_vec, _mm512_mul_ps(fulltmp_z_vec, p_ij_z_vec)); )

      #if MIC_HANDCODE_FORCE_COMBINE_FORCES == 0
        APPLY_FORCES_PS2PD(fulltmp_x_vec, fulltmp_y_vec, fulltmp_z_vec, fullf_0_x, fullf_0_y, fullf_0_z, fullf_1_x, fullf_1_y, fullf_1_z, tmpI32, tmpJ32);
      #endif

    #endif // FULL

    #if MIC_HANDCODE_FORCE_COMBINE_FORCES != 0
      APPLY_FORCES_PS2PD;
    #endif

    // DMK - DEBUG - LARGE FORCE DETECTION
    #if (MIC_DATA_STRUCT_VERIFY != 0) && (MIC_DATA_STRUCT_VERIFY_KERNELS != 0)
      const float lfd_const = 3.2e2f; // Set to 5.0e1 to print a few per timestep for ApoA1
      __mmask16 lfd_mask = _mm512_int2mask(0x0000);
      #define LFD_CHECK(v)  _mm512_mask_cmplt_ps_mask(cutoff_mask, _mm512_set_1to16_ps(lfd_const), _mm512_mask_mul_ps(v, _mm512_cmplt_ps_mask(v, _mm512_setzero_ps()), v, _mm512_set_1to16_ps(-1.0f)))
      #if (0 FAST(SHORT(+1)))
        lfd_mask = _mm512_kor(lfd_mask, LFD_CHECK(tmp_x_vec));
        lfd_mask = _mm512_kor(lfd_mask, LFD_CHECK(tmp_y_vec));
        lfd_mask = _mm512_kor(lfd_mask, LFD_CHECK(tmp_z_vec));
      #endif
      #if (0 FULL(+1))
        lfd_mask = _mm512_kor(lfd_mask, LFD_CHECK(fulltmp_x_vec));
        lfd_mask = _mm512_kor(lfd_mask, LFD_CHECK(fulltmp_y_vec));
        lfd_mask = _mm512_kor(lfd_mask, LFD_CHECK(fulltmp_z_vec));
      #endif
      #undef LFD_CHECK
      #pragma unroll(16)
      for (int _k = 0; _k < 16; _k++) {
        lfd_mask = _mm512_kor(lfd_mask, _mm512_int2mask(((pExt_0[tmpI32[_k]].index == 53839 && pExt_1[tmpJ32[_k]].index == 53827) ? (0x1) : (0x0)) << _k));
      }
      if (_mm512_mask2int(lfd_mask) != 0) {

        float tmp_force[16] __attribute__((aligned(64)));
        float tmp_fullforce[16] __attribute__((aligned(64)));
        float tmp_ijx[16] __attribute__((aligned(64)));
        float tmp_ijy[16] __attribute__((aligned(64)));
        float tmp_ijz[16] __attribute__((aligned(64)));
        float tmp_x[16] __attribute__((aligned(64)));
        float tmp_y[16] __attribute__((aligned(64)));
        float tmp_z[16] __attribute__((aligned(64)));
        float fulltmp_x[16] __attribute__((aligned(64)));
        float fulltmp_y[16] __attribute__((aligned(64)));
        float fulltmp_z[16] __attribute__((aligned(64)));
        float pi_x[16] __attribute__((aligned(64)));
        float pi_y[16] __attribute__((aligned(64)));
        float pi_z[16] __attribute__((aligned(64)));
        float pi_q[16] __attribute__((aligned(64)));
        float pj_x[16] __attribute__((aligned(64)));
        float pj_y[16] __attribute__((aligned(64)));
        float pj_z[16] __attribute__((aligned(64)));
        float pj_q[16] __attribute__((aligned(64)));
        float r2[16] __attribute__((aligned(64)));
        int table_i[16] __attribute__((aligned(64)));
        float kqq[16] __attribute__((aligned(64)));
        float diffa[16] __attribute__((aligned(64)));
        float vdw_a[16] __attribute__((aligned(64)));
        float vdw_b[16] __attribute__((aligned(64)));
        float vdw_c[16] __attribute__((aligned(64)));
        float vdw_d[16] __attribute__((aligned(64)));
        float fast_a[16] __attribute__((aligned(64)));
        float fast_b[16] __attribute__((aligned(64)));
        float fast_c[16] __attribute__((aligned(64)));
        float fast_d[16] __attribute__((aligned(64)));
        float A[16] __attribute__((aligned(64)));
        float B[16] __attribute__((aligned(64)));
        float table_four_i_0[16] __attribute__((aligned(64)));
        float table_four_i_1[16] __attribute__((aligned(64)));
        float table_four_i_2[16] __attribute__((aligned(64)));
        float table_four_i_3[16] __attribute__((aligned(64)));
        float table_four_i_4[16] __attribute__((aligned(64)));
        float table_four_i_5[16] __attribute__((aligned(64)));
        float table_four_i_6[16] __attribute__((aligned(64)));
        float table_four_i_7[16] __attribute__((aligned(64)));
        int i_vdwType[16] __attribute__((aligned(64)));
        int j_vdwType[16] __attribute__((aligned(64)));

        _mm512_store_ps(pi_x, p_i_x_vec);
        _mm512_store_ps(pi_y, p_i_y_vec);
        _mm512_store_ps(pi_z, p_i_z_vec);
        _mm512_store_ps(pi_q, p_i_q_vec);
        _mm512_store_ps(pj_x, p_j_x_vec);
        _mm512_store_ps(pj_y, p_j_y_vec);
        _mm512_store_ps(pj_z, p_j_z_vec);
        _mm512_store_ps(pj_q, p_j_q_vec);
        _mm512_store_ps(r2, r2_vec);
        _mm512_store_epi32(table_i, table_i_vec);
        _mm512_store_ps(kqq, kqq_vec);
        _mm512_store_ps(diffa, diffa_vec);
        _mm512_store_ps(table_four_i_0, table_four_i_0_vec);
        _mm512_store_ps(table_four_i_1, table_four_i_1_vec);
        _mm512_store_ps(table_four_i_2, table_four_i_2_vec);
        _mm512_store_ps(table_four_i_3, table_four_i_3_vec);
        _mm512_store_ps(table_four_i_4, table_four_i_4_vec);
        _mm512_store_ps(table_four_i_5, table_four_i_5_vec);
        _mm512_store_ps(table_four_i_6, table_four_i_6_vec);
        _mm512_store_ps(table_four_i_7, table_four_i_7_vec);
        #if (0 FAST(+1))
          _mm512_store_ps(vdw_a, vdw_a_vec);
          _mm512_store_ps(vdw_b, vdw_b_vec);
          _mm512_store_ps(vdw_c, vdw_c_vec);
          _mm512_store_ps(vdw_d, vdw_d_vec);
          _mm512_store_ps(A, A_vec);
          _mm512_store_ps(B, B_vec);
          _mm512_store_epi64(i_vdwType, pExt_i_vdwType_vec);
          _mm512_store_epi64(j_vdwType, pExt_j_vdwType_vec);
        #else
          _mm512_store_ps(vdw_a, _mm512_setzero_ps());
          _mm512_store_ps(vdw_b, _mm512_setzero_ps());
          _mm512_store_ps(vdw_c, _mm512_setzero_ps());
          _mm512_store_ps(vdw_d, _mm512_setzero_ps());
          _mm512_store_ps(A, _mm512_setzero_ps());
          _mm512_store_ps(B, _mm512_setzero_ps());
          _mm512_store_epi64(i_vdwType, _mm512_setzero_epi32());
          _mm512_store_epi64(j_vdwType, _mm512_setzero_epi32());
        #endif
        #if (0 FAST(SHORT(+1)))
          _mm512_store_ps(tmp_force, force_r_vec);
          _mm512_store_ps(tmp_x, tmp_x_vec);
          _mm512_store_ps(tmp_y, tmp_y_vec);
          _mm512_store_ps(tmp_z, tmp_z_vec);
          _mm512_store_ps(fast_a, fast_a_vec);
          _mm512_store_ps(fast_b, fast_b_vec);
          _mm512_store_ps(fast_c, fast_c_vec);
          _mm512_store_ps(fast_d, fast_d_vec);
        #else
          _mm512_store_ps(tmp_force, _mm512_setzero_ps());
          _mm512_store_ps(tmp_x, _mm512_setzero_ps());
          _mm512_store_ps(tmp_y, _mm512_setzero_ps());
          _mm512_store_ps(tmp_z, _mm512_setzero_ps());
          _mm512_store_ps(fast_a, _mm512_setzero_ps());
          _mm512_store_ps(fast_b, _mm512_setzero_ps());
          _mm512_store_ps(fast_c, _mm512_setzero_ps());
          _mm512_store_ps(fast_d, _mm512_setzero_ps());
        #endif
        #if (0 FULL(+1))
          _mm512_store_ps(tmp_fullforce, fullforce_r_vec);
          _mm512_store_ps(fulltmp_x, fulltmp_x_vec);
          _mm512_store_ps(fulltmp_y, fulltmp_y_vec);
          _mm512_store_ps(fulltmp_z, fulltmp_z_vec);
        #else
          _mm512_store_ps(tmp_fullforce, _mm512_setzero_ps());
          _mm512_store_ps(fulltmp_x, _mm512_setzero_ps());
          _mm512_store_ps(fulltmp_y, _mm512_setzero_ps());
          _mm512_store_ps(fulltmp_z, _mm512_setzero_ps());
        #endif
        _mm512_store_ps(tmp_ijx, p_ij_x_vec);
        _mm512_store_ps(tmp_ijy, p_ij_y_vec);
        _mm512_store_ps(tmp_ijz, p_ij_z_vec);

        #define PRNT(idx) \
          if (_mm512_mask2int(lfd_mask) & (1 << (idx))) { \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :  tmp x:%+.3e, y:%+.3e, z:%+.3e  :  fulltmp x:%+.3e, y:%+.3e, z:%+.3e\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   tmp_x[idx], tmp_y[idx], tmp_z[idx], \
                   fulltmp_x[idx], fulltmp_y[idx], fulltmp_z[idx] \
		  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    calc_" SELF("self") PAIR("pair") FULL("_fullelect") ENERGY("_energy") "." NORMAL("NORM") MODIFIED("MOD") EXCLUDED("EXCL") "\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]) \
	          ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    i:%ld, j:%ld, i_index:%d, j_index:%d\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   tmpI32[idx], tmpJ32[idx], pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    p_i x:%+.3e, y:%+.3e, z:%+.3e, q:%+.3e  :  p_j x:%+.3e, y:%+.3e, z:%+.3e, q:%+.3e  :  r2 %+.3e\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   pi_x[idx], pi_y[idx], pi_z[idx], pi_q[idx], \
                   pj_x[idx], pj_y[idx], pj_z[idx], pj_q[idx], \
                   r2[idx] \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    offset x:%+.3e, y:%+.3e, z:%+.3e\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   offset_x, offset_y, offset_z \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    forces f&s:%+.3e, full:%+.3e  :  p_ij x:%+.3e, y:%+.3e, z:%+.3e\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   tmp_force[idx], tmp_fullforce[idx], \
                   tmp_ijx[idx], tmp_ijy[idx], tmp_ijz[idx] \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    table_i:%d, kqq::%+.3e, diffa:%+.3e\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   table_i[idx], kqq[idx], diffa[idx] \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    vdw a:%+.3e, b:%+.3e, c:%+.3e, d:%+.3e  :  fast a:%+.3e, b:%+.3e, c:%+.3e, d:%+.3e, \n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   vdw_a[idx], vdw_b[idx], vdw_c[idx], vdw_d[idx], \
                   fast_a[idx], fast_b[idx], fast_c[idx], fast_d[idx] \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    A:%+.3e, B:%+.3e  :  vdwType i:%d, j:%d\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   A[idx], B[idx], i_vdwType[idx], j_vdwType[idx] \
                  ); \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    table_four_i 0:%+.3e, 1:%+.3e, 2:%+.3e, 3:%+.3e, 4:%+.3e, 5:%+.3e, 6:%+.3e, 7:%+.3e\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   table_four_i_0[idx], table_four_i_1[idx], table_four_i_2[idx], table_four_i_3[idx], table_four_i_4[idx], table_four_i_5[idx], table_four_i_6[idx], table_four_i_7[idx] \
                  ); \
            int indexDiff = pExt_0[tmpI32[idx]].index - pExt_1[tmpJ32[idx]].index; \
            int maxDiff = pExt_1[tmpJ32[idx]].excl_maxdiff; \
            int offset = -1, offset_major = -1, offset_minor = -1, flags = 0; \
            if (indexDiff >= -1 * maxDiff && indexDiff <= maxDiff) { \
              offset = (2 * indexDiff) + pExt_1[tmpJ32[idx]].excl_index; \
              offset_major = offset / (sizeof(unsigned int) * 8); \
              offset_minor = offset % (sizeof(unsigned int) * 8); \
              flags = ((params.exclusion_bits[offset_major]) >> offset_minor) & 0x03; \
            } \
            printf("_lfd :: %05d %06d : %8d, %8d : %06d, %3d, %3d, %3d, %3d  :    excl: exclIndex:%d, indexDiff:%d, maxDiff:%d, offset:%d, flags:%d\n", \
                   params.pe, params.timestep, \
                   pExt_0[tmpI32[idx]].index, pExt_1[tmpJ32[idx]].index, \
                   params.ppI, params.p1, params.p2, (int)(tmpI32[idx]), (int)(tmpJ32[idx]), \
                   pExt_1[tmpJ32[idx]].excl_index, indexDiff, maxDiff, offset, flags \ 
                  ); \
          }
        PRNT( 0); PRNT( 1); PRNT( 2); PRNT( 3); PRNT( 4); PRNT( 5); PRNT( 6); PRNT( 7);
        PRNT( 8); PRNT( 9); PRNT(10); PRNT(11); PRNT(12); PRNT(13); PRNT(14); PRNT(15);
        #undef PRNT
      } // end if (_mm512_mask2int(lfd_mask) != 0)
    #endif  // (MIC_DATA_STRUCT_VERIFY != 0) && (MIC_DATA_STRUCT_VERIFY_KERNELS != 0)
  }

  FAST(ENERGY( vdwEnergy += _mm512_reduce_add_pd(vdwEnergy_vec); ))
  FAST(SHORT(ENERGY( electEnergy += _mm512_reduce_add_pd(electEnergy_vec); )))
  FULL(ENERGY( fullElectEnergy += _mm512_reduce_add_pd(fullElectEnergy_vec); ))
  #if (0 PAIR(FAST(SHORT(+1))))
    virial_xx += _mm512_reduce_add_pd(virial_xx_vec);
    virial_xy += _mm512_reduce_add_pd(virial_xy_vec);
    virial_xz += _mm512_reduce_add_pd(virial_xz_vec);
    virial_yy += _mm512_reduce_add_pd(virial_yy_vec);
    virial_yz += _mm512_reduce_add_pd(virial_yz_vec);
    virial_zz += _mm512_reduce_add_pd(virial_zz_vec);
  #endif
  #if (0 PAIR(FULL(+1)))
    fullElectVirial_xx += _mm512_reduce_add_pd(fullElectVirial_xx_vec);
    fullElectVirial_xy += _mm512_reduce_add_pd(fullElectVirial_xy_vec);
    fullElectVirial_xz += _mm512_reduce_add_pd(fullElectVirial_xz_vec);
    fullElectVirial_yy += _mm512_reduce_add_pd(fullElectVirial_yy_vec);
    fullElectVirial_yz += _mm512_reduce_add_pd(fullElectVirial_yz_vec);
    fullElectVirial_zz += _mm512_reduce_add_pd(fullElectVirial_zz_vec);
  #endif

  #if (MIC_EXCL_CHECKSUM_FULL != 0) && (0 EXCLUDED(+1) MODIFIED(+1))
    params.exclusionSum += _mm512_reduce_add_epi32(exclusionSum_vec);
  #endif

  #undef GATHER_PS_I32_OFFSET
  #undef GATHER_PS_I32
  #undef GATHER_EPI32_I32_OFFSET
  #undef GATHER_EPI32_I32
  #undef SCATTER_INC_PD_I32_STEP
  #undef SCATTER_INC_PD_I32
  #undef SCATTER_INC_PS_I32_STEP
  #undef SCATTER_INC_PS_I32
  #undef CONTRIB_ADD_PS2PD
  #undef CONTRIB_SUB_PS2PD
  #undef APPLY_FORCES_PS2PD_STEP_ADD
  #undef APPLY_FORCES_PS2PD_STEP_SUB
  #undef APPLY_FORCES_PS2PD_STEP_ADD_COMBO
  #undef APPLY_FORCES_PS2PD_STEP_SUB_COMBO
  #undef APPLY_FORCES_PS2PD_JSTEP
  #undef APPLY_FORCES_PS2PD

#endif  // NAMD_MIC
