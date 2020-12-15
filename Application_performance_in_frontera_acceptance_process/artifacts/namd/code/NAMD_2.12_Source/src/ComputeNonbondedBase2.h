/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

EXCLUDED( FAST( foo bar ) )
EXCLUDED( MODIFIED( foo bar ) )
EXCLUDED( NORMAL( foo bar ) )
NORMAL( MODIFIED( foo bar ) )
ALCHPAIR( NOT_ALCHPAIR( foo bar ) )

ALCHPAIR(
  // get alchemical nonbonded scaling parameters (once per pairlist)
  myLambda = ALCH1(lambdaUp) ALCH2(lambdaDown);
  FEP(myLambda2 = ALCH1(lambda2Up) ALCH2(lambda2Down);)
  myElecLambda =  ALCH1(elecLambdaUp) ALCH2(elecLambdaDown); 
  FEP(myElecLambda2 = ALCH1(elecLambda2Up) ALCH2(elecLambda2Down);)
  myVdwLambda =  ALCH1(vdwLambdaUp) ALCH2(vdwLambdaDown); 
  FEP(myVdwLambda2 = ALCH1(vdwLambda2Up) ALCH2(vdwLambda2Down);) 
  myVdwShift =  ALCH1(vdwShiftUp) ALCH2(vdwShiftDown); 
  FEP(myVdwShift2 =  ALCH1(vdwShift2Up) ALCH2(vdwShift2Down);) 
)

#ifdef  A2_QPX
#if ( SHORT(1+) 0 )
NORMAL(kq_iv    = vec_splats(kq_i); )
MODIFIED(kq_iv    = vec_splats((1.0-modf_mod) *kq_i); )
#endif

#if ( FULL( 1+ ) 0 )
EXCLUDED(
	 SHORT(
	       full_cnst = (vector4double)(6., 4., 2., 1.);
	       ) 
	 NOSHORT(
		 full_cnst = (vector4double)(1., 1., 1., 1.);
		 )
	 )
MODIFIED(
	 SHORT(
	       full_cnst = (vector4double)(6., 4., 2., 1.);
	       full_cnst = vec_mul (full_cnst, vec_splats(modf_mod));          
	       ) 
	 NOSHORT(
		 full_cnst = vec_splats(modf_mod);
		 )
	 )      
#endif
#endif

#ifdef ARCH_POWERPC
     __alignx(64, table_four);
     __alignx(32, p_1);
#pragma unroll(1)
#pragma ibm independent_loop
#endif

#ifndef ARCH_POWERPC
#pragma ivdep
#endif


  
#if ( FULL( EXCLUDED( SHORT( 1+ ) ) ) 0 ) 
// avoid bug in Intel 15.0 compiler
#pragma novector
#else
#ifdef PRAGMA_SIMD
#ifndef TABENERGYFLAG
#ifndef GOFORCES
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
      TABENERGY(
      const int numtypes = simParams->tableNumTypes;
      const float table_spacing = simParams->tableSpacing;
      const int npertype = (int) (mynearbyint(simParams->tableMaxDist / simParams->tableSpacing) + 1);
      )

      int table_i = (r2iilist[2*k] >> 14) + r2_delta_expc;  // table_i >= 0 
      const int j = pairlisti[k];
      //register const CompAtom *p_j = p_1 + j;
#define p_j (p_1+j)
#ifdef  A2_QPX
      register double *p_j_d = (double *) p_j;
#endif
      // const CompAtomExt *pExt_j = pExt_1 + j;
      
      BigReal diffa = r2list[k] - r2_table[table_i];
      //const BigReal* const table_four_i = table_four + 16*table_i;
#define table_four_i (table_four + 16*table_i)

#if  ( FAST( 1 + ) TABENERGY( 1 + ) 0 ) // FAST or TABENERGY
      //const LJTable::TableEntry * lj_pars = 
      //        lj_row + 2 * p_j->vdwType MODIFIED(+ 1);
      const int lj_index = 2 * p_j->vdwType MODIFIED(+ 1);
#define lj_pars (lj_row+lj_index)
#ifdef  A2_QPX
      double *lj_pars_d = (double *) lj_pars;
#endif
#endif
      
      TABENERGY(
      register const int tabtype = -1 - ( lj_pars->A < 0 ? lj_pars->A : 0 );
      )
      
#if ( SHORT( FAST( 1+ ) ) 0 ) 
      //Force *f_j = f_1 + j;
#define f_j (f_1+j)
#endif
	
#if ( FULL( 1+ ) 0 )
      //Force *fullf_j = fullf_1 + j;
#define fullf_j (fullf_1+j)
#endif

      //Power PC aliasing and alignment constraints
#ifdef ARCH_POWERPC      
#if ( FULL( 1+ ) 0 )
#pragma disjoint (*table_four, *fullf_1)
#pragma disjoint (*p_1,          *fullf_1)
#ifdef  A2_QPX
#pragma disjoint (*p_j_d,        *fullf_1)
#endif
#pragma disjoint (*r2_table,     *fullf_1)
#pragma disjoint (*r2list,       *fullf_1)
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*f_1    ,      *fullf_1)
#pragma disjoint (*fullf_1,      *f_1)
#endif   //Short + fast
#endif   //Full

#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*table_four, *f_1)
#pragma disjoint (*p_1,          *f_1)
#pragma disjoint (*r2_table,     *f_1)
#pragma disjoint (*r2list,       *f_1)
#pragma disjoint (*lj_row,      *f_1)
#ifdef  A2_QPX
#pragma disjoint (*p_j_d,        *f_1)
#endif
#endif //Short + Fast

      __alignx(64, table_four_i);
      FAST (
      __alignx(32, lj_pars);
      )
      __alignx(32, p_j);
#endif   //ARCH_POWERPC

      /*
      BigReal modf = 0.0;
      int atom2 = p_j->id;
      register char excl_flag = ( (atom2 >= excl_min && atom2 <= excl_max) ?
					excl_flags[atom2-excl_min] : 0 );
      if ( excl_flag ) { ++exclChecksum; }
      SELF( if ( j < j_hgroup ) { excl_flag = EXCHCK_FULL; } )
      if ( excl_flag ) {
	if ( excl_flag == EXCHCK_FULL ) {
	  lj_pars = lj_null_pars;
	  modf = 1.0;
	} else {
	  ++lj_pars;
	  modf = modf_mod;
	}
      }
      */

      BigReal kqq = kq_i * p_j->charge;

      
#ifdef  A2_QPX
      float * cg = (float *)&p_j->charge;
#if ( FULL( 1+ ) 0 )
#pragma disjoint (*cg, *fullf_1)
#endif   //Full

#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*cg, *f_1)
#endif   //Short + fast
#endif

      LES( BigReal lambda_pair = lambda_table_i[p_j->partition]; )

#ifndef  A2_QPX
      register const BigReal p_ij_x = p_i_x - p_j->position.x;
      register const BigReal p_ij_y = p_i_y - p_j->position.y;
      register const BigReal p_ij_z = p_i_z - p_j->position.z;
#else
      vector4double charge_v = vec_lds(0, cg);
      vector4double kqqv = vec_mul(kq_iv, charge_v );
      vector4double p_ij_v = vec_sub(p_i_v, vec_ld (0, p_j_d));
#define p_ij_x     vec_extract(p_i_v, 0)
#define p_ij_y     vec_extract(p_i_v, 1)
#define p_ij_z     vec_extract(p_i_v, 2)
#endif  

#if ( FAST(1+) 0 )
      const BigReal A = scaling * lj_pars->A;
      const BigReal B = scaling * lj_pars->B;
#ifndef  A2_QPX
      BigReal vdw_d = A * table_four_i[0] - B * table_four_i[4];
      BigReal vdw_c = A * table_four_i[1] - B * table_four_i[5];
      BigReal vdw_b = A * table_four_i[2] - B * table_four_i[6];
      BigReal vdw_a = A * table_four_i[3] - B * table_four_i[7];
#else
      const vector4double Av = vec_mul(scalingv, vec_lds(0, lj_pars_d));
      const vector4double Bv = vec_mul(scalingv, vec_lds(8, lj_pars_d));
      vector4double vdw_v = vec_msub( Av, vec_ld(0, (BigReal*)table_four_i), vec_mul(Bv, vec_ld(4*sizeof(BigReal), (BigReal*)table_four_i)) );
#define   vdw_d  vec_extract(vdw_v, 0)
#define   vdw_c  vec_extract(vdw_v, 1)
#define   vdw_b  vec_extract(vdw_v, 2)
#define   vdw_a  vec_extract(vdw_v, 3)
#endif

      ALCHPAIR (
        // Alchemical free energy calculation
        // Pairlists are separated so that lambda-coupled pairs are handled
        // independently from normal nonbonded (inside ALCHPAIR macro).
        // The separation-shifted van der Waals potential and a shifted 
        // electrostatics potential for decoupling are calculated explicitly.
        // Would be faster with lookup tables but because only a small minority
        // of nonbonded pairs are lambda-coupled the impact is minimal. 
        // Explicit calculation also makes things easier to modify.
        
        const BigReal r2 = r2list[k] - r2_delta;

        // These are now inline functions (in ComputeNonbondedFep.C) to 
        // tidy the code

        FEP(fep_vdw_forceandenergies(A,B,r2,myVdwShift,myVdwShift2,switchdist2,
          cutoff2, myVdwLambda, myVdwLambda2, Fep_WCA_repuOn, Fep_WCA_dispOn,
          Fep_Wham, WCA_rcut1, WCA_rcut2, WCA_rcut3, switchfactor,
          vdwForceSwitching, &alch_vdw_energy, &alch_vdw_force, 
          &alch_vdw_energy_2, &alch_vdw_energy_2_Left);)
        TI(ti_vdw_force_energy_dUdl(A,B,r2,myVdwShift,switchdist2,
          cutoff2, myVdwLambda, alchVdwShiftCoeff, switchfactor,
          vdwForceSwitching, &alch_vdw_energy, &alch_vdw_force, 
          &alch_vdw_dUdl);)
      )
      
	//NOT_ALCHPAIR(
	//TABENERGY(
#if (NOT_ALCHPAIR(1+) 0)
#if (TABENERGY(1+) 0)
        if (tabtype >= 0) {
          register BigReal r1;
          r1 = sqrt(p_ij_x*p_ij_x + p_ij_y*p_ij_y + p_ij_z*p_ij_z);

          //CkPrintf("%i %i %f %f %i\n", npertype, tabtype, r1, table_spacing, (int) (mynearbyint(r1 / table_spacing)));
          register int eneraddress;
          eneraddress = 2 * ((npertype * tabtype) + ((int) mynearbyint(r1 / table_spacing)));
          //CkPrintf("Using distance bin %i for distance %f\n", eneraddress, r1);
#ifndef A2_QPX
	  vdw_d = 0.;
	  vdw_c = 0.;
	  vdw_b = table_ener[eneraddress + 1] / r1;
	  vdw_a = (-1/2.) * diffa * vdw_b;
#else
	  vec_insert(0.,                               vdw_v, 0);
	  vec_insert(0.,                               vdw_v, 1);
	  vec_insert(table_ener[eneraddress + 1] / r1, vdw_v, 2);
	  vec_insert((-1/2.) * diffa * vdw_b,          vdw_v, 3);
#endif
	  ENERGY(
            register BigReal vdw_val = table_ener[eneraddress];
            //CkPrintf("Found vdw energy of %f\n", vdw_val);
            vdwEnergy += LAM(lambda_pair *) vdw_val;
            FEP( vdwEnergy_s += d_lambda_pair * vdw_val; )
            FEP( vdwEnergy_s_Left += d_lambda_pair * vdw_val; )
          )
        }  else {
	  //)
#endif
      ENERGY(
        register BigReal vdw_val =
          ( ( diffa * vdw_d * (1/6.)+ vdw_c * (1/4.)) * diffa + vdw_b *(1/2.)) * diffa + vdw_a;

        vdwEnergy -= LAM(lambda_pair *) vdw_val;

        FEP(vdwEnergy_s -= vdw_val;)
        FEP(vdwEnergy_s_Left -= vdw_val;)
      )
	//TABENERGY( } ) /* endif (tabtype >= 0) */
#if (TABENERGY (1+) 0)
	}
#endif
      //) // NOT_ALCHPAIR
#endif

      ALCHPAIR(
        ENERGY(vdwEnergy   += alch_vdw_energy;)
        FEP(vdwEnergy_s += alch_vdw_energy_2;)
        FEP(vdwEnergy_s_Left += alch_vdw_energy_2_Left;)
        TI(ALCH1(vdwEnergy_ti_1) ALCH2(vdwEnergy_ti_2) += alch_vdw_dUdl;)
      ) // ALCHPAIR
      
#endif // FAST

#if ( FAST(1+) 0 )
     INT( 
      register BigReal vdw_dir;
      vdw_dir = ( diffa * vdw_d + vdw_c ) * diffa + vdw_b;
      //BigReal force_r =  LAM(lambda_pair *) vdw_dir;
      reduction[pairVDWForceIndex_X] += force_sign * vdw_dir * p_ij_x;
      reduction[pairVDWForceIndex_Y] += force_sign * vdw_dir * p_ij_y;
      reduction[pairVDWForceIndex_Z] += force_sign * vdw_dir * p_ij_z;
      )

#if ( SHORT(1+) 0 ) // Short-range electrostatics

#ifndef  A2_QPX
      NORMAL(
      BigReal fast_d = kqq * table_four_i[8];
      BigReal fast_c = kqq * table_four_i[9];
      BigReal fast_b = kqq * table_four_i[10];
      BigReal fast_a = kqq * table_four_i[11];
      )
      MODIFIED(
      BigReal modfckqq = (1.0-modf_mod) * kqq;
      BigReal fast_d = modfckqq * table_four_i[8];
      BigReal fast_c = modfckqq * table_four_i[9];
      BigReal fast_b = modfckqq * table_four_i[10];
      BigReal fast_a = modfckqq * table_four_i[11];
      )
#else
      vector4double fastv = vec_mul(kqqv, vec_ld(8 * sizeof(BigReal), (BigReal*)table_four_i));
#define fast_d   vec_extract(fastv, 0)
#define fast_c   vec_extract(fastv, 1)
#define fast_b   vec_extract(fastv, 2)
#define fast_a   vec_extract(fastv, 3)
#endif
    
      {
      ENERGY(
      register BigReal fast_val =
	( ( diffa * fast_d * (1/6.)+ fast_c * (1/4.)) * diffa + fast_b *(1/2.)) * diffa + fast_a;
      NOT_ALCHPAIR (
        electEnergy -=  LAM(lambda_pair *) fast_val;
        FEP(electEnergy_s -= fast_val;)
      )
      ) //ENERGY
      ALCHPAIR(
        ENERGY(electEnergy   -= myElecLambda * fast_val;)
        if( Fep_Wham ) {
        	if( Fep_ElecOn ) {
            FEP(electEnergy_s -= (myElecLambda * fast_val + fast_val);)
        	}
        	else {	// repulsive or dispersive, no constribution from charges
            FEP(electEnergy_s -= (myElecLambda * fast_val);)
        	}
        }
        else{ // the default (orginal) FEP
          FEP(electEnergy_s -= myElecLambda2 * fast_val;)
        } 
        TI(
          NOENERGY(register BigReal fast_val = 
            ( ( diffa * fast_d * (1/6.)+ fast_c * (1/4.)) * diffa + fast_b *(1/2.)) * diffa + fast_a;)
          ALCH1(electEnergy_ti_1) ALCH2(electEnergy_ti_2)  -= fast_val;
        )
      )

      INT(
      register BigReal fast_dir =
      ( diffa * fast_d + fast_c ) * diffa + fast_b;
      // force_r -= -1.0 * LAM(lambda_pair *) fast_dir;
      reduction[pairElectForceIndex_X] +=  force_sign * fast_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] +=  force_sign * fast_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] +=  force_sign * fast_dir * p_ij_z;
      )
      }


     /*****  JE - Go  *****/
     // Now Go energy should appear in VDW place -- put vdw_b back into place
#if    ( NORMAL (1+) 0)
#if    ( GO (1+) 0)


// JLai
#ifndef CODE_REDUNDANT
#define CODE_REDUNDANT 0
#endif
#if CODE_REDUNDANT
       if (ComputeNonbondedUtil::goGroPair) {
	 // Explicit goGroPair calculation; only calculates goGroPair if goGroPair is turned on
	 //
	 // get_gro_force has an internal checklist that sees if atom_i and atom_j are
	 // in the explicit pairlist.  This is done because there is no guarantee that a 
	 // processor will have atom_i and atom_j so we cannot loop over the explict atom pairs.
	 // We can only loop over all pairs.
	 //
	 // NOTE: It does not look like fast_b is not normalized by the r vector.
	 //
	 // JLai
	 BigReal groLJe = 0.0;
	 BigReal groGausse = 0.0;
	 const CompAtomExt *pExt_z = pExt_1 + j;
             BigReal groForce = mol->get_gro_force2(p_ij_x, p_ij_y, p_ij_z,pExt_i.id,pExt_z->id,&groLJe,&groGausse);
	     NAMD_die("Failsafe.  This line should never be reached\n");
#ifndef A2_QPX
             fast_b += groForce;
#else 
             vec_insert(fast_b + groForce, fastv, 2);
#endif
	 ENERGY(
	     NOT_ALCHPAIR (
		 // JLai
		 groLJEnergy += groLJe;
		 groGaussEnergy += groGausse;
		 )
	     ) //ENERGY
       }
#endif 
       BigReal goNative = 0;
       BigReal goNonnative = 0;
       BigReal goForce = 0;
       register const CompAtomExt *pExt_j = pExt_1 + j;
       if (ComputeNonbondedUtil::goMethod == 2) {
	 goForce = mol->get_go_force2(p_ij_x, p_ij_y, p_ij_z, pExt_i.id, pExt_j->id,&goNative,&goNonnative);
       } else {
	 //  Ported by JLai -- JE - added (
	 const BigReal r2go = square(p_ij_x, p_ij_y, p_ij_z);
	 const BigReal rgo = sqrt(r2go);
       
	 if (ComputeNonbondedUtil::goMethod == 1) {
	   goForce = mol->get_go_force(rgo, pExt_i.id, pExt_j->id, &goNative, &goNonnative);
	 } else if (ComputeNonbondedUtil::goMethod == 3) {  
	   goForce = mol->get_go_force_new(rgo, pExt_i.id, pExt_j->id, &goNative, &goNonnative);
	 } else {
	   NAMD_die("I SHOULDN'T BE HERE.  DYING MELODRAMATICALLY.\n");
	 }
       }
       
#ifndef A2_QPX
       fast_b += goForce;
#else 
       vec_insert(fast_b + goForce, fastv, 2);
#endif
       {
       ENERGY(
	 NOT_ALCHPAIR (
		       // JLai
		       goEnergyNative +=  goNative;
		       goEnergyNonnative += goNonnative;
	 )
       ) //ENERGY                                                                                                                                             	   
	 INT(
	   reduction[pairVDWForceIndex_X] +=  force_sign * goForce * p_ij_x;
	   reduction[pairVDWForceIndex_Y] +=  force_sign * goForce * p_ij_y;
	   reduction[pairVDWForceIndex_Z] +=  force_sign * goForce * p_ij_z;
	 )
       }
       // End of INT 

       //DebugM(3,"rgo:" << rgo << ", pExt_i.id:" << pExt_i.id << ", pExt_j->id:" << pExt_j->id << \
	 //      ", goForce:" << goForce << ", fast_b:" << fast_b << std::endl);
#endif       //     ) // End of GO macro 
       /*****  JE - End Go  *****/
       // End of port JL
#endif      //) // End of Normal MACRO

      // Combined short-range electrostatics and VdW force:
#if    (  NOT_ALCHPAIR(1+) 0)
#ifndef A2_QPX
        fast_d += vdw_d;
        fast_c += vdw_c;
        fast_b += vdw_b;
        fast_a += vdw_a;  // not used!
#else
	fastv = vec_add(fastv, vdw_v);
#endif
#endif

      register BigReal fast_dir =
                  (diffa * fast_d + fast_c) * diffa + fast_b;

      BigReal force_r =  LAM(lambda_pair *) fast_dir;
      ALCHPAIR(
        force_r *= myElecLambda; 
        force_r += alch_vdw_force;
        // special ALCH forces already multiplied by relevant lambda
      )
          
#ifndef NAMD_CUDA
#ifndef  A2_QPX
      register BigReal tmp_x = force_r * p_ij_x;
      f_i_x += tmp_x;
      f_j->x -= tmp_x;

      register BigReal tmp_y = force_r * p_ij_y;
      f_i_y += tmp_y;
      f_j->y -= tmp_y;
      
      register BigReal tmp_z = force_r * p_ij_z;
      f_i_z += tmp_z;
      f_j->z -= tmp_z;
#else
      vector4double force_rv = vec_splats (force_r);
      vector4double tmp_v = vec_mul(force_rv, p_ij_v);
      f_i_v = vec_add(f_i_v, tmp_v);

#define tmp_x   vec_extract(tmp_v, 0)
#define tmp_y   vec_extract(tmp_v, 1)
#define tmp_z   vec_extract(tmp_v, 2)

      f_j->x -= tmp_x;
      f_j->y -= tmp_y;
      f_j->z -= tmp_z;
#endif

      PPROF(
        const BigReal p_j_z = p_j->position.z;
        int n2 = (int)floor((p_j_z-pressureProfileMin)*invThickness);
        pp_clamp(n2, pressureProfileSlabs);
        int p_j_partition = p_j->partition;

        pp_reduction(pressureProfileSlabs, n1, n2, 
                     p_i_partition, p_j_partition, pressureProfileAtomTypes,
                     tmp_x*p_ij_x, tmp_y * p_ij_y, tmp_z*p_ij_z,
                     pressureProfileReduction);
      )
#endif

#endif // SHORT
#endif // FAST

#if ( FULL (EXCLUDED( SHORT ( 1+ ) ) ) 0 ) 
      //const BigReal* const slow_i = slow_table + 4*table_i;
#define slow_i (slow_table + 4*table_i)

#ifdef ARCH_POWERPC  //Alignment and aliasing constraints
      __alignx (32, slow_table);
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*slow_table, *f_1)
#endif
#pragma disjoint (*slow_table, *fullf_1)
#endif  //ARCH_POWERPC

#endif //FULL 


#if ( FULL (MODIFIED( SHORT ( 1+ ) ) ) 0 ) 
      //const BigReal* const slow_i = slow_table + 4*table_i;
#define slow_i (slow_table + 4*table_i)

#ifdef ARCH_POWERPC //Alignment and aliasing constraints
      __alignx (32, slow_table);
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*slow_table, *f_1)
#endif
#pragma disjoint (*slow_table, *fullf_1)
#endif //ARCH_POWERPC

#endif //FULL
      
#if ( FULL( 1+ ) 0 )
#ifndef  A2_QPX
      BigReal slow_d = table_four_i[8 SHORT(+ 4)];
      BigReal slow_c = table_four_i[9 SHORT(+ 4)];
      BigReal slow_b = table_four_i[10 SHORT(+ 4)];
      BigReal slow_a = table_four_i[11 SHORT(+ 4)];
      EXCLUDED(
      SHORT(
      slow_a +=    slow_i[3];
      slow_b += 2.*slow_i[2];
      slow_c += 4.*slow_i[1];
      slow_d += 6.*slow_i[0];
      )
      NOSHORT(
      slow_d -= table_four_i[12];
      slow_c -= table_four_i[13];
      slow_b -= table_four_i[14];
      slow_a -= table_four_i[15];
      )
      )
      MODIFIED(
      SHORT(
      slow_a +=    modf_mod * slow_i[3];
      slow_b += 2.*modf_mod * slow_i[2];
      slow_c += 4.*modf_mod * slow_i[1];
      slow_d += 6.*modf_mod * slow_i[0];
      )
      NOSHORT(
      slow_d -= modf_mod * table_four_i[12];
      slow_c -= modf_mod * table_four_i[13];
      slow_b -= modf_mod * table_four_i[14];
      slow_a -= modf_mod * table_four_i[15];
      )
      )
      slow_d *= kqq;
      slow_c *= kqq;
      slow_b *= kqq;
      slow_a *= kqq;
#else
      vector4double slow_v = vec_ld((8 SHORT(+ 4)) * sizeof(BigReal), (BigReal*)table_four_i);
      EXCLUDED(
	       SHORT(
		     slow_v = vec_madd(full_cnst, vec_ld(0, (BigReal*)slow_i), slow_v);
		     )
	       NOSHORT(
		       slow_v = vec_sub(slow_v, vec_ld(12*sizeof(BigReal), (BigReal*)table_four_i));
		       )
	       );
      MODIFIED(
	       SHORT(
		     slow_v = vec_madd(full_cnst,  vec_ld(0, (BigReal*)slow_i), slow_v);
		     )
	       NOSHORT(
		       slow_v = vec_nmsub(full_cnst, vec_ld(12*sizeof(BigReal), (BigReal*)table_four_i), slow_v);
		       )
	       );
      slow_v = vec_mul (slow_v, vec_splats(kqq));

#define slow_d   vec_extract(slow_v, 0)
#define slow_c   vec_extract(slow_v, 1)
#define slow_b   vec_extract(slow_v, 2)
#define slow_a   vec_extract(slow_v, 3)

#endif

      ENERGY(
      register BigReal slow_val =
        ( ( diffa * slow_d *(1/6.)+ slow_c * (1/4.)) * diffa + slow_b *(1/2.)) * diffa + slow_a;
      
      NOT_ALCHPAIR (
        fullElectEnergy -= LAM(lambda_pair *) slow_val;
        FEP(fullElectEnergy_s -= slow_val;) 
      )
      ) // ENERGY
          
      ALCHPAIR(
        ENERGY(fullElectEnergy   -= myElecLambda * slow_val;)
        if( Fep_Wham ) {
        	if(Fep_ElecOn)	{
        		FEP(fullElectEnergy_s -= (myElecLambda * slow_val + slow_val);)
        	}
        	else	{
            FEP(fullElectEnergy_s -= myElecLambda * slow_val;)
        	}
        }
        else{ // orignal FEP
          FEP(fullElectEnergy_s -= myElecLambda2 * slow_val;)
        }
        TI(
          NOENERGY(register BigReal slow_val =
	        ( ( diffa * slow_d *(1/6.)+ slow_c * (1/4.)) * diffa + slow_b *(1/2.)) * diffa + slow_a;)
          ALCH1(fullElectEnergy_ti_1) ALCH2(fullElectEnergy_ti_2) -= slow_val;
        )
      )

      INT( {
      register BigReal slow_dir =
	( diffa * slow_d + slow_c ) * diffa + slow_b;
      reduction[pairElectForceIndex_X] += force_sign * slow_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] += force_sign * slow_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] += force_sign * slow_dir * p_ij_z;
      } )


#if     (NOT_ALCHPAIR (1+) 0)
#if     (FAST(1+) 0)
#if     (NOSHORT(1+) 0)
#ifndef A2_QPX
        slow_d += vdw_d;
        slow_c += vdw_c;
        slow_b += vdw_b;
        slow_a += vdw_a; // unused!
#else
	slow_v = vec_add (slow_v, vdw_v);
#endif
#endif
#endif
#endif

      register BigReal slow_dir = (diffa * slow_d + slow_c) * diffa + slow_b;
      BigReal fullforce_r = slow_dir LAM(* lambda_pair);
      ALCHPAIR (
        fullforce_r *= myElecLambda;
        FAST( NOSHORT(
          fullforce_r += alch_vdw_force;
        ))
      )
          
#ifndef NAMD_CUDA
      {
#ifndef  A2_QPX
      register BigReal ftmp_x = fullforce_r * p_ij_x;
      fullf_i_x += ftmp_x;
      fullf_j->x -= ftmp_x;
      register BigReal ftmp_y = fullforce_r * p_ij_y;
      fullf_i_y += ftmp_y;
      fullf_j->y -= ftmp_y;
      register BigReal ftmp_z = fullforce_r * p_ij_z;
      fullf_i_z += ftmp_z;
      fullf_j->z -= ftmp_z;
#else
      vector4double fforce_rv = vec_splats (fullforce_r);
      vector4double ftmp_v = vec_mul(fforce_rv, p_ij_v);
      fullf_i_v = vec_add(fullf_i_v, ftmp_v);
      
#define ftmp_x  vec_extract(ftmp_v, 0)
#define ftmp_y  vec_extract(ftmp_v, 1)
#define ftmp_z  vec_extract(ftmp_v, 2)

      fullf_j->x -= ftmp_x;
      fullf_j->y -= ftmp_y;
      fullf_j->z -= ftmp_z;

#endif

      PPROF(
        const BigReal p_j_z = p_j->position.z;
        int n2 = (int)floor((p_j_z-pressureProfileMin)*invThickness);
        pp_clamp(n2, pressureProfileSlabs);
        int p_j_partition = p_j->partition;

        pp_reduction(pressureProfileSlabs, n1, n2, 
                     p_i_partition, p_j_partition, pressureProfileAtomTypes,
                     ftmp_x*p_ij_x, ftmp_y * p_ij_y, ftmp_z*p_ij_z,
                     pressureProfileReduction);

      )
      }
#endif
#endif //FULL

   } // for pairlist

#undef p_j
#undef lj_pars
#undef table_four_i
#undef slow_i
#undef f_j
#undef fullf_j

