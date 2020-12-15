#ifdef NAMD_MIC

  // For each entry in the pairlist...

  // Auto-vectorize via pairlist padding
  #if __MIC_PAD_PLGEN_CTRL != 0

    // Set the number of elements/lanes per vector unit width for the data type that will be used
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      const int _plI_fs_outer_step = 16;  // 32-bit
    #else
      const int _plI_fs_outer_step = 8;  // 64-bit
    #endif

    // Create an "outer" loop that iterates over the the entire loop, stepping by the
    //   number of lanes in the vector units
    #pragma novector
    for (int _plI_fs_outer = 0; _plI_fs_outer < plSize; _plI_fs_outer += _plI_fs_outer_step) {

      // Preload i value here (use broadcast)...
      const int i = (plArray[_plI_fs_outer] >> 16) & 0xFFFF;

      // Preload x,y,z,q values here
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0[i].x) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0[i].y) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0[i].z) + ((CALC_TYPE)params.offset.z);
        const CALC_TYPE p_i_q = (CALC_TYPE)(p_0[i].charge);
        const int p_i_vdwType = pExt_0[i].vdw_type;
      #else
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0_x[i]) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0_y[i]) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0_z[i]) + ((CALC_TYPE)params.offset.z);
        const CALC_TYPE p_i_q = (CALC_TYPE)(p_0_q[i]);
        const int p_i_vdwType = pExt_0_vdwType[i];
      #endif

      // Create variables to hold the force contributions for the given "i" atom in the "inner" loop below
      double tmp_x_i_sum = 0.0;
      double tmp_y_i_sum = 0.0;
      double tmp_z_i_sum = 0.0;
      double tmp_w_i_sum = 0.0;
      double fulltmp_x_i_sum = 0.0;
      double fulltmp_y_i_sum = 0.0;
      double fulltmp_z_i_sum = 0.0;

      #if MIC_EXCL_CHECKSUM_FULL != 0
        int exclusionSum = 0;
        #define EXCL_CHECKSUM_CLAUSE  reduction(+ : exclusionSum)
      #else
        #define EXCL_CHECKSUM_CLAUSE
      #endif

      // Create an "inner" loop with one iteration per vector unit lane
      #pragma simd vectorlength(16) \
                   reduction(+ : tmp_x_i_sum, tmp_y_i_sum, tmp_z_i_sum, tmp_w_i_sum, \
                                 fulltmp_x_i_sum, fulltmp_y_i_sum, fulltmp_z_i_sum ) \
                   EXCL_CHECKSUM_CLAUSE
      for (int _plI_fs_inner = 0; _plI_fs_inner < _plI_fs_outer_step; _plI_fs_inner++) {
        const int plI = _plI_fs_outer + _plI_fs_inner;
        if ((plArray[plI] & 0xFFFF) != 0xFFFF) {

  // Scalar version of the code
  #else

    // DMK - NOTE : These loop_count values are loose, lower-bound guesses on my part (TODO : verify & refine)
    #if (0 PAIR(+1))
      #pragma loop_count (1000)
    #elif (0 SELF(+1))
      #pragma loop_count (10000)
    #endif
    for (int plI = 0; plI < plSize; plI++) {  

  #endif

    // Load the particle indicies
    const int ij = plArray[plI];
    #if __MIC_PAD_PLGEN_CTRL != 0
      // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
    #else
      const int i = (ij >> 16) & 0xFFFF;
    #endif
    const int j = (ij      ) & 0xFFFF;

    // TODO | FIXME - Spread these out throughout the loop body (if possible) and
    //   change based on AoS versus SoA
    #if MIC_PREFETCH_DISTANCE > 0
      const int pfIJ = plArray[plI + MIC_PREFETCH_DISTANCE];
      const int pfI = (pfIJ >> 16) & 0xFFFF;
      const int pfJ = (pfIJ      ) & 0xFFFF;
      _mm_prefetch((char*)(p_0_x + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_0_y + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_0_z + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_0_q + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_0_x + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_0_y + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_0_z + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_x + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_y + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_z + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_q + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_1_x + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_1_y + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_1_z + pfJ), MIC_PREFETCH_HINT);
    #endif

    // Load atom information
    #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
      #if __MIC_PAD_PLGEN_CTRL != 0
        // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
      #else
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0[i].x) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0[i].y) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0[i].z) + ((CALC_TYPE)params.offset.z);
      #endif
      const CALC_TYPE p_j_x = (CALC_TYPE)(p_1[j].x);  // Neighboring gather to be optimized
      const CALC_TYPE p_j_y = (CALC_TYPE)(p_1[j].y);
      const CALC_TYPE p_j_z = (CALC_TYPE)(p_1[j].z);
    #else
      #if __MIC_PAD_PLGEN_CTRL != 0
        // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
      #else
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0_x[i]) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0_y[i]) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0_z[i]) + ((CALC_TYPE)params.offset.z);
      #endif
      const CALC_TYPE p_j_x = (CALC_TYPE)(p_1_x[j]);
      const CALC_TYPE p_j_y = (CALC_TYPE)(p_1_y[j]);
      const CALC_TYPE p_j_z = (CALC_TYPE)(p_1_z[j]);
    #endif

    // Load position deltas and r2
    CALC_TYPE p_ij_x = p_i_x - p_j_x;
    CALC_TYPE p_ij_y = p_i_y - p_j_y;
    CALC_TYPE p_ij_z = p_i_z - p_j_z;

    #if REFINE_PAIRLISTS != 0
    CALC_TYPE r2 = (CALC_TYPE)(r2Array[plI]);
    #else
    CALC_TYPE r2 = (p_ij_x * p_ij_x) + (p_ij_y * p_ij_y) + (p_ij_z * p_ij_z) + r2_delta;
    if (r2 < cutoff2_delta) {
    #endif

      #if (MIC_EXCL_CHECKSUM_FULL != 0) && (0 EXCLUDED(+1) MODIFIED(+1))
        #if __MIC_PAD_PLGEN_CTRL != 0
          exclusionSum += 1;
        #else
          params.exclusionSum += 1;
        #endif
      #endif

      // Calculate the table_i value (table index)
      #if MIC_HANDCODE_FORCE_SINGLE != 0
        const unsigned int table_i = ((int)((__intel_castf32_u32(r2)) >> 17)) + r2_delta_expc;
      #else
        const unsigned int table_i = ((int)((__intel_castf64_u64(r2)) >> 46)) + r2_delta_expc;
      #endif

      #if MIC_HANDCODE_FORCE_CALCR2TABLE != 0
        // From ComputeNonbondedUtil.C                    Simplified:
        //   r2_base = r2_delta * (1 << (i/64))             r2_base = r2_delta * (1 << (i/64))
        //   r2_del = r2_base / 64.0;                       r2_del = r2_base / 64.0;
        //   r2 = r2_base - r2_delta + r2_del * (i%64)      r2_table[i] = r2_base - r2_delta + r2_del * (i%64) + r2_delta;
        //   r2_table[i] = r2 + r2_delta;                               = r2_base + r2_del * (i%64)
        // NOTE: For i = 0, r2_table[0] = r2_delta + (r2_delta / 64) * 0 = r2_delta, so there no need
        //   to special case if table_i = 0 then r2_table[0] = r2_delta (see ComputeNonbondedUtil.C:606)
        CALC_TYPE r2_base = r2_delta * (1 << (table_i >> 6)); // avoid original divide (table_i / 64)
        CALC_TYPE r2_del = r2_base * ((CALC_TYPE)0.015625f);  // avoid original divide (r2_base / 64)
        CALC_TYPE r2_table_i = r2_base + r2_del * (table_i & 0x3F); //(table_i % 64);  // NOTE: removing '+ r2_delta - r2_delta'
      #else
        CALC_TYPE r2_table_i = r2_table[table_i];
      #endif
      CALC_TYPE diffa = r2 - r2_table_i;
      const CALC_TYPE * const table_four_ptr = SHORT(table_short) NOSHORT(table_noshort);
      const int table_four_idx = 16 * table_i;

      // NOTE : These charge values are already scaled by
      //   'sqrt(COULOMB * scaling * dielectric_1).'  See HomePatch.C.
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        #if __MIC_PAD_PLGEN_CTRL != 0
          // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
        #else
          const CALC_TYPE p_i_q = (CALC_TYPE)(p_0[i].charge);
        #endif
        const CALC_TYPE p_j_q = (CALC_TYPE)(p_1[j].charge);
      #else
        #if __MIC_PAD_PLGEN_CTRL != 0
          // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
        #else
          const CALC_TYPE p_i_q = (CALC_TYPE)(p_0_q[i]);
        #endif
        const CALC_TYPE p_j_q = (CALC_TYPE)(p_1_q[j]);
      #endif
      CALC_TYPE kqq = p_i_q * p_j_q;

      #if (0 FAST(+1))

        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
          #else
            int p_i_vdwType = pExt_0[i].vdw_type;
          #endif
          int p_j_vdwType = pExt_1[j].vdw_type;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
          #else
            int p_i_vdwType = pExt_0_vdwType[i];
          #endif
          int p_j_vdwType = pExt_1_vdwType[j];
        #endif

        // Lookup A and B values in the LJ table
        const int lj_pars_offset = (4 * (p_i_vdwType * lj_table_dim + p_j_vdwType)) MODIFIED(+ 2);
        CALC_TYPE A = scaling * lj_table_base_ptr[lj_pars_offset    ];
        CALC_TYPE B = scaling * lj_table_base_ptr[lj_pars_offset + 1];

        // 16x16 AoS table lookup with transpose
        //CALC_TYPE vdw_d = A * table_four_ptr[table_four_idx + 0] - B * table_four_ptr[table_four_idx + 2];
        //CALC_TYPE vdw_c = A * table_four_ptr[table_four_idx + 1] - B * table_four_ptr[table_four_idx + 3];
        //CALC_TYPE vdw_b = A * table_four_ptr[table_four_idx + 4] - B * table_four_ptr[table_four_idx + 6];
        //CALC_TYPE vdw_a = A * table_four_ptr[table_four_idx + 5] - B * table_four_ptr[table_four_idx + 7];
        CALC_TYPE vdw_d = A * table_four_ptr[table_four_idx + 0] - B * table_four_ptr[table_four_idx + 4];
        CALC_TYPE vdw_c = A * table_four_ptr[table_four_idx + 1] - B * table_four_ptr[table_four_idx + 5];
        CALC_TYPE vdw_b = A * table_four_ptr[table_four_idx + 2] - B * table_four_ptr[table_four_idx + 6];
        CALC_TYPE vdw_a = A * table_four_ptr[table_four_idx + 3] - B * table_four_ptr[table_four_idx + 7];

        #if (0 ENERGY(+1))
          CALC_TYPE vdw_val = ((diffa * vdw_d * (1/6.0) + vdw_c * (1/4.0)) * diffa + vdw_b * (1/2.0)) * diffa + vdw_a;
          vdwEnergy -= vdw_val;
          // DMK - TODO | FIXME : Apply vdw_val to FEP(vdwEnergy_s)
        #endif

        #if (0 SHORT(+1))

          #if (0 NORMAL(+1))
            CALC_TYPE fast_d = kqq * table_four_ptr[table_four_idx +  8];
            CALC_TYPE fast_c = kqq * table_four_ptr[table_four_idx +  9];
            CALC_TYPE fast_b = kqq * table_four_ptr[table_four_idx + 10];
            CALC_TYPE fast_a = kqq * table_four_ptr[table_four_idx + 11];
          #endif
          #if (0 MODIFIED(+1))
            CALC_TYPE modfckqq = (1.0 - modf_mod) * kqq;
            CALC_TYPE fast_d = modfckqq * table_four_ptr[table_four_idx +  8];
            CALC_TYPE fast_c = modfckqq * table_four_ptr[table_four_idx +  9];
            CALC_TYPE fast_b = modfckqq * table_four_ptr[table_four_idx + 10];
            CALC_TYPE fast_a = modfckqq * table_four_ptr[table_four_idx + 11];
          #endif

          #if (0 ENERGY(+1))
            CALC_TYPE fast_val = ((diffa * fast_d * (1/6.0) + fast_c * (1/4.0)) * diffa + fast_b * (1/2.0)) * diffa + fast_a;
            #if (0 NOT_ALCHPAIR(+1))
              electEnergy -= fast_val;
              // DMK - TODO | FIXME : Apply fast_val to FEP(electEnergy_s)
            #endif
          #endif

          #if (0 NOT_ALCHPAIR(+1))
            fast_d += vdw_d;
            fast_c += vdw_c;
            fast_b += vdw_b;
            fast_a += vdw_a;
          #endif

          CALC_TYPE fast_dir = (fast_d * diffa + fast_c) * diffa + fast_b;
          CALC_TYPE force_r = fast_dir;

          CALC_TYPE tmp_x = force_r * p_ij_x;
          PAIR( virial_xx += tmp_x * p_ij_x; )
          PAIR( virial_xy += tmp_x * p_ij_y; )
          PAIR( virial_xz += tmp_x * p_ij_z; )
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_x_i_sum += tmp_x;
            #else
              f_0[i].x += tmp_x;
            #endif
            f_1[j].x -= tmp_x;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_x_i_sum += tmp_x;
            #else
              f_0_x[i] += tmp_x;
            #endif
            f_1_x[j] -= tmp_x;
          #endif

          CALC_TYPE tmp_y = force_r * p_ij_y;
          PAIR( virial_yy += tmp_y * p_ij_y; )
          PAIR( virial_yz += tmp_y * p_ij_z; )
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_y_i_sum += tmp_y;
            #else
              f_0[i].y += tmp_y;   /// Move out after inner loop
            #endif
            f_1[j].y -= tmp_y;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_y_i_sum += tmp_y;
            #else
              f_0_y[i] += tmp_y;
            #endif
            f_1_y[j] -= tmp_y;
          #endif

          CALC_TYPE tmp_z = force_r * p_ij_z;
          PAIR( virial_zz += tmp_z * p_ij_z; )
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_z_i_sum += tmp_z;
            #else
              f_0[i].z += tmp_z;
            #endif
            f_1[j].z -= tmp_z;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_z_i_sum += tmp_z;
            #else
              f_0_z[i] += tmp_z;
            #endif
            f_1_z[j] -= tmp_z;
          #endif

        #endif // SHORT
      #endif // FAST

      #if (0 FULL(+1))

        CALC_TYPE slow_d = table_four_ptr[table_four_idx +  8 SHORT(+ 4)];
        CALC_TYPE slow_c = table_four_ptr[table_four_idx +  9 SHORT(+ 4)];
        CALC_TYPE slow_b = table_four_ptr[table_four_idx + 10 SHORT(+ 4)];
        CALC_TYPE slow_a = table_four_ptr[table_four_idx + 11 SHORT(+ 4)];

        #if (0 SHORT( EXCLUDED(+1) MODIFIED(+1) ))
          const int slow_idx = 4 * table_i;
        #endif
        #if (0 EXCLUDED(+1))
          #if (0 SHORT(+1))
	    //slow_a += 1.0 * slow_table[slow_idx + 0];  // AoS transpose (4 members)
            //slow_b += 2.0 * slow_table[slow_idx + 1];
            //slow_c += 4.0 * slow_table[slow_idx + 2];
            //slow_d += 6.0 * slow_table[slow_idx + 3];
            slow_a += 1.0 * slow_table[slow_idx + 3];  // AoS transpose (4 members)
            slow_b += 2.0 * slow_table[slow_idx + 2];
            slow_c += 4.0 * slow_table[slow_idx + 1];
            slow_d += 6.0 * slow_table[slow_idx + 0];
          #endif
          #if (0 NOSHORT(+1))
            slow_d -= table_four_ptr[table_four_idx + 12];
            slow_c -= table_four_ptr[table_four_idx + 13];
            slow_b -= table_four_ptr[table_four_idx + 14];
            slow_a -= table_four_ptr[table_four_idx + 15];
          #endif
        #endif
        #if (0 MODIFIED(+1))
          #if (0 SHORT(+1))
            //slow_a += 1.0 * modf_mod * slow_table[slow_idx + 0];
            //slow_b += 2.0 * modf_mod * slow_table[slow_idx + 1];
            //slow_c += 4.0 * modf_mod * slow_table[slow_idx + 2];
            //slow_d += 6.0 * modf_mod * slow_table[slow_idx + 3];
            slow_a += 1.0 * modf_mod * slow_table[slow_idx + 3];
            slow_b += 2.0 * modf_mod * slow_table[slow_idx + 2];
            slow_c += 4.0 * modf_mod * slow_table[slow_idx + 1];
            slow_d += 6.0 * modf_mod * slow_table[slow_idx + 0];
          #endif
          #if (0 NOSHORT(+1))
            slow_d -= modf_mod * table_four_ptr[table_four_idx + 12];
            slow_c -= modf_mod * table_four_ptr[table_four_idx + 13];
            slow_b -= modf_mod * table_four_ptr[table_four_idx + 14];
            slow_a -= modf_mod * table_four_ptr[table_four_idx + 15];
          #endif
        #endif
        slow_d *= kqq;
        slow_c *= kqq;
        slow_b *= kqq;
        slow_a *= kqq;

        #if (0 ENERGY(+1))
          CALC_TYPE slow_val = ((diffa * slow_d * (1/6.0) + slow_c * (1/4.0)) * diffa + slow_b * (1/2.0)) * diffa + slow_a;
          #if (0 NOT_ALCHPAIR(+1))
            fullElectEnergy -= slow_val;
            // DMK - TODO | FIXME : Apply slow_val to FEP(fullElectEnergy_s)
          #endif
        #endif

        #if (0 NOT_ALCHPAIR(FAST(NOSHORT(+1))))
          slow_d += vdw_d;
          slow_c += vdw_c;
          slow_b += vdw_b;
          slow_a += vdw_a;
        #endif

        CALC_TYPE slow_dir = (diffa * slow_d + slow_c) * diffa + slow_b;
        CALC_TYPE fullforce_r = slow_dir;

        CALC_TYPE fulltmp_x = fullforce_r * p_ij_x;
        PAIR( fullElectVirial_xx += fulltmp_x * p_ij_x; )
        PAIR( fullElectVirial_xy += fulltmp_x * p_ij_y; )
        PAIR( fullElectVirial_xz += fulltmp_x * p_ij_z; )
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_x_i_sum += fulltmp_x;
          #else
            fullf_0[i].x += fulltmp_x;
          #endif
          fullf_1[j].x -= fulltmp_x;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_x_i_sum += fulltmp_x;
          #else
            fullf_0_x[i] += fulltmp_x;
          #endif
          fullf_1_x[j] -= fulltmp_x;
        #endif

        CALC_TYPE fulltmp_y = fullforce_r * p_ij_y;
        PAIR( fullElectVirial_yy += fulltmp_y * p_ij_y; )
        PAIR( fullElectVirial_yz += fulltmp_y * p_ij_z; )
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_y_i_sum += fulltmp_y;
          #else
            fullf_0[i].y += fulltmp_y;
          #endif
          fullf_1[j].y -= fulltmp_y;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_y_i_sum += fulltmp_y;
          #else
            fullf_0_y[i] += fulltmp_y;
          #endif
          fullf_1_y[j] -= fulltmp_y;
        #endif

        CALC_TYPE fulltmp_z = fullforce_r * p_ij_z;
        PAIR( fullElectVirial_zz += fulltmp_z * p_ij_z; )
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_z_i_sum += fulltmp_z;
          #else
            fullf_0[i].z += fulltmp_z;
          #endif
          fullf_1[j].z -= fulltmp_z;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_z_i_sum += fulltmp_z;
          #else
            fullf_0_z[i] += fulltmp_z;
          #endif
          fullf_1_z[j] -= fulltmp_z;
        #endif

      #endif // FULL

    #if REFINE_PAIRLISTS == 0
    } // end if (r2 < cutoff2_delta)
    #endif

  // End of loops auto-vectorized via pairlist padding
  #if __MIC_PAD_PLGEN_CTRL != 0

      } // end if
    } // end for

    #if MIC_EXCL_CHECKSUM_FULL != 0
      params.exclusionSum += exclusionSum;
    #endif
    #undef EXCL_CHECKSUM_CLAUSE

    #if (0 FAST(SHORT(+1)))
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        f_0[i].x += tmp_x_i_sum;
        f_0[i].y += tmp_y_i_sum;
        f_0[i].z += tmp_z_i_sum;
      #else
        f_0_x[i] += tmp_x_i_sum;
        f_0_y[i] += tmp_y_i_sum;
        f_0_z[i] += tmp_z_i_sum;
      #endif
    #endif

    #if (0 FULL(+1))
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        fullf_0[i].x += fulltmp_x_i_sum;
        fullf_0[i].y += fulltmp_y_i_sum;
        fullf_0[i].z += fulltmp_z_i_sum;
      #else
        fullf_0_x[i] += fulltmp_x_i_sum;
        fullf_0_y[i] += fulltmp_y_i_sum;
        fullf_0_z[i] += fulltmp_z_i_sum;
      #endif
    #endif

  } // end for


  // End of scalar loop
  #else

  } // end pairlist-loop

  #endif

#endif  // NAMD_MIC
