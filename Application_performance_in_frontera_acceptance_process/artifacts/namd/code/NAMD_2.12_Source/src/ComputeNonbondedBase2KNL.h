/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef KNL_MAKE_DEPENDS_INCLUDE

#if 0
#if (VDW_SWITCH_MODE == VDW_SWITCH_MODE_FORCE) || (VDW_SWITCH_MODE == VDW_SWITCH_MODE_MARTINI)
// REMOVE WHEN r2list NO LONGER NEEDED!
#pragma ivdep
for (k=0; k<npairi; ++k) {      
  r2list[k] = r2list_f[k] + r2_delta;
}
#endif
#endif

EXCLUDED( FAST( foo bar ) )
EXCLUDED( MODIFIED( foo bar ) )
EXCLUDED( NORMAL( foo bar ) )
NORMAL( MODIFIED( foo bar ) )
ALCHPAIR( NOT_ALCHPAIR( foo bar ) )

EXCLUDED( foo bar )
MODIFIED( foo bar )
ALCHPAIR( foo bar )
TABENERGY( foo bar )
NOFAST( foo bar )

#pragma ivdep

#if ( FULL( EXCLUDED( SHORT( 1+ ) ) ) 0 ) 
// avoid bug in Intel 15.0 compiler
#pragma novector
#else
#ifdef PRAGMA_SIMD
#ifndef TABENERGYFLAG
#if __INTEL_COMPILER_BUILD_DATE == 20160721
#warning disabled pragma simd on innner loop due to compiler segfault
#else
#pragma simd assert SHORT(FAST(reduction(+:f_i_x,f_i_y,f_i_z)) ENERGY(FAST(reduction(+:vdwEnergy) SHORT(reduction(+:electEnergy))))) \
             FULL(reduction(+:fullf_i_x,fullf_i_y,fullf_i_z) ENERGY(reduction(+:fullElectEnergy)))
#endif
#endif
#pragma loop_count avg=100
#else // PRAGMA_SIMD
#pragma loop_count avg=4
#endif // PRAGMA_SIMD
#endif
    for (k=0; k<npairi; ++k) {      

      const float r2 = r2list_f[k];
      const float r_1 = 1.f / sqrtf(r2);
      const float r_2 = r_1 * r_1;
      const float knl_table_r_1 = r_1 > 1.f ? 1.f : r_1;
      const float knl_table_f = (KNL_TABLE_SIZE-2) * knl_table_r_1;
      const int knl_table_i = knl_table_f;
      const float knl_diff = knl_table_f - knl_table_i;

      const int j = pairlisti[k];
      //register const CompAtom *p_j = p_1 + j;
#define p_j (p_1+j)
#define pFlt_j (pFlt_1+j)

#if (VDW_SWITCH_MODE == VDW_SWITCH_MODE_FORCE) || (VDW_SWITCH_MODE == VDW_SWITCH_MODE_MARTINI)
#if 0
      int table_i = (r2iilist[2*k] >> 14) + r2_delta_expc;  // table_i >= 0 
      
      float diffa = r2list[k] - r2_table[table_i];
      //const BigReal* const table_four_i = table_four + 16*table_i;
#define table_four_i (table_four + 16*table_i)
#endif
#endif

      //const LJTable::TableEntry * lj_pars = 
      //        lj_row + 2 * p_j->vdwType MODIFIED(+ 1);
      const int lj_index = 2 * pFlt_j->vdwType MODIFIED(+ 1);
#define lj_pars (lj_row+lj_index)
      
#if ( SHORT( 1+ ) 0 ) 
      //Force *f_j = f_1 + j;
#define f_j (f_1+j)
#endif
	
#if ( FULL( 1+ ) 0 )
      //Force *fullf_j = fullf_1 + j;
#define fullf_j (fullf_1+j)
#endif

      float kqq = kq_i_f * p_j->charge;

      LES( float lambda_pair = lambda_table_i[p_j->partition]; )

      register const float p_ij_x = xlist[k];
      register const float p_ij_y = ylist[k];
      register const float p_ij_z = zlist[k];

      const float A = scaling_f * lj_pars->A;
      const float B = scaling_f * lj_pars->B;

#if VDW_SWITCH_MODE == VDW_SWITCH_MODE_FORCE
      { int vdw_switch_mode_force; }  // for preprocessor debugging only
      float vdw_b = 0.f;
      {
        const float r_6 = r_2 * r_2 * r_2;
        float vdwa_energy, vdwb_energy, vdwa_gradient, vdwb_gradient;
        // from Steinbach & Brooks, JCC 15, pgs 667-683, 1994, eqns 10-13
        if ( r2 > switchOn2_f ) {
          const float tmpa = r_6 - cutoff_6_f;
          vdwa_energy = k_vdwa_f * tmpa * tmpa;
          const float tmpb = r_1 * r_2 - cutoff_3_f;
          vdwb_energy = k_vdwb_f * tmpb * tmpb;
          vdwa_gradient = -6.f * k_vdwa_f * tmpa * r_2 * r_6;
          vdwb_gradient = -3.f * k_vdwb_f * tmpb * r_2 * r_2 * r_1;
        } else {
          const float r_12 = r_6 * r_6;
          vdwa_energy = r_12 + v_vdwa_f;
          vdwb_energy = r_6 + v_vdwb_f;
          vdwa_gradient = -6.f * r_2 * r_12;
          vdwb_gradient = -3.f * r_2 * r_6;
        }
        vdw_b = -2.f * ( A * vdwa_gradient - B * vdwb_gradient );
        ENERGY(
          vdwEnergy += A * vdwa_energy - B * vdwb_energy;
        )
      }
#elif VDW_SWITCH_MODE == VDW_SWITCH_MODE_MARTINI
      { int vdw_switch_mode_martini; }  // for preprocessor debugging only
      float vdw_b = 0.f;
      {
        const float r = r2 * r_1;
        const float r12 = (r-switchOn_f)*(r-switchOn_f);
        const float r13 = (r-switchOn_f)*(r-switchOn_f)*(r-switchOn_f);

        ENERGY(
          const float LJshifttempA = -(1.f/3.f)*A12_f*r13 - (1.f/4.f)*B12_f*r12*r12 - C12_f;
          const float LJshifttempB = -(1.f/3.f)*A6_f*r13 - (1.f/4.f)*B6_f*r12*r12 - C6_f;
          const float shiftValA = ( r > switchOn_f ? LJshifttempA : -C12_f);
          const float shiftValB = ( r > switchOn_f ? LJshifttempB : -C6_f);
        )

        const float LJdshifttempA = -A12_f*r12 - B12_f*r13;
        const float LJdshifttempB = -A6_f*r12 - B6_f*r13;
        const float dshiftValA = ( r > switchOn_f ? LJdshifttempA*0.5f*r_1 : 0 );
        const float dshiftValB = ( r > switchOn_f ? LJdshifttempB*0.5f*r_1 : 0 );

        const float r_6 = r_2 * r_2 * r_2;
        const float r_12 = r_6 * r_6;

        ENERGY(
          const float vdwa_energy = r_12 + shiftValA;
          const float vdwb_energy = r_6 + shiftValB;
        )

        const float vdwa_gradient = -6.f * r_2 * r_12 + dshiftValA ;
        const float vdwb_gradient = -3.f * r_2 * r_6 + dshiftValB;

        vdw_b = -2.f * ( A * vdwa_gradient - B * vdwb_gradient );
        ENERGY(
          vdwEnergy += A * vdwa_energy - B * vdwb_energy;
        )
      }
#elif VDW_SWITCH_MODE == VDW_SWITCH_MODE_ENERGY
      { int vdw_switch_mode_energy; }  // for preprocessor debugging only
      float vdw_b = 0.f;
      {
        const float r_6 = r_2 * r_2 * r_2;
        const float r_12 = r_6 * r_6;
        const float c2 = cutoff2_f-r2;
        const float c4 = c2*(c3_f-2.f*c2);
        const float switchVal =         // used for Lennard-Jones
                        ( r2 > switchOn2_f ? c2*c4*c1_f : 1.f );
        const float dSwitchVal =        // d switchVal / d r2
                        ( r2 > switchOn2_f ? 2.f*c1_f*(c2*c2-c4) : 0.f );
        const float vdwa_gradient = ( dSwitchVal - 6.f * switchVal * r_2 ) * r_12;
        const float vdwb_gradient = ( dSwitchVal - 3.f * switchVal * r_2 ) * r_6;
        vdw_b = -2.f * ( A * vdwa_gradient - B * vdwb_gradient );
        ENERGY(
          vdwEnergy += switchVal * ( A * r_12 - B * r_6 );
        )
      }
#else
#error VDW_SWITCH_MODE not recognized
#endif  // VDW_SWITCH_MODE

#if ( SHORT(1+) 0 ) // Short-range electrostatics

      NORMAL(
      float fast_b = kqq * ( knl_fast_grad_table[knl_table_i] * (1.f-knl_diff) +
                             knl_fast_grad_table[knl_table_i+1] * knl_diff );
      )

      {
      ENERGY(
        float fast_val = kqq * ( knl_fast_ener_table[knl_table_i] * (1.f-knl_diff) +
                                 knl_fast_ener_table[knl_table_i+1] * knl_diff );
        electEnergy -=  LAM(lambda_pair *) fast_val;
      ) //ENERGY
      }

      // Combined short-range electrostatics and VdW force:
        fast_b += vdw_b;

      float fast_dir = fast_b;

      float force_r =  LAM(lambda_pair *) fast_dir;
          
      BigReal tmp_x = force_r * p_ij_x;
      f_i_x += tmp_x;
      f_j->x -= tmp_x;

      BigReal tmp_y = force_r * p_ij_y;
      f_i_y += tmp_y;
      f_j->y -= tmp_y;
      
      BigReal tmp_z = force_r * p_ij_z;
      f_i_z += tmp_z;
      f_j->z -= tmp_z;

#endif // SHORT

#if ( FULL( 1+ ) 0 )
  #if ( SHORT( 1+ ) 0 )
      float slow_b = kqq * ( knl_scor_grad_table[knl_table_i] * (1.f-knl_diff) +
                             knl_scor_grad_table[knl_table_i+1] * knl_diff );
      ENERGY(
        float slow_val = kqq * ( knl_scor_ener_table[knl_table_i] * (1.f-knl_diff) +
                                 knl_scor_ener_table[knl_table_i+1] * knl_diff );
      )
  #else
      float slow_b = kqq * ( knl_corr_grad_table[knl_table_i] * (1.f-knl_diff) +
                             knl_corr_grad_table[knl_table_i+1] * knl_diff );
      ENERGY(
        float slow_val = kqq * ( knl_corr_ener_table[knl_table_i] * (1.f-knl_diff) +
                                 knl_corr_ener_table[knl_table_i+1] * knl_diff );
      )
  #endif

      ENERGY(
        fullElectEnergy -= LAM(lambda_pair *) slow_val;
      ) // ENERGY
          
#if     (NOSHORT(1+) 0)
        slow_b += vdw_b;
#endif

      register float slow_dir = slow_b;
      float fullforce_r = slow_dir LAM(* lambda_pair);
          
      {
      BigReal ftmp_x = fullforce_r * p_ij_x;
      fullf_i_x += ftmp_x;
      fullf_j->x -= ftmp_x;
      BigReal ftmp_y = fullforce_r * p_ij_y;
      fullf_i_y += ftmp_y;
      fullf_j->y -= ftmp_y;
      BigReal ftmp_z = fullforce_r * p_ij_z;
      fullf_i_z += ftmp_z;
      fullf_j->z -= ftmp_z;
      }
#endif //FULL

   } // for pairlist

#undef p_j
#undef lj_pars
#undef table_four_i
#undef slow_i
#undef f_j
#undef fullf_j

#endif // KNL_MAKE_DEPENDS_INCLUDE

