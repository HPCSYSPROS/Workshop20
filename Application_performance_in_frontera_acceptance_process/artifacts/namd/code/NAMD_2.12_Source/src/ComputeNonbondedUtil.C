/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifdef WIN32
extern "C" {
  double erfc(double);
}
#endif

#include "InfoStream.h"
#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"
#include "Parameters.h"
#include "MsmMacros.h"
#include <stdio.h>

#ifdef NAMD_CUDA
  void send_build_cuda_force_table();
#endif

#ifdef NAMD_MIC
  extern void send_build_mic_force_table();
#endif

Bool		ComputeNonbondedUtil::commOnly;
Bool		ComputeNonbondedUtil::fixedAtomsOn;
Bool            ComputeNonbondedUtil::qmForcesOn;
BigReal         ComputeNonbondedUtil::cutoff;
BigReal         ComputeNonbondedUtil::cutoff2;
float           ComputeNonbondedUtil::cutoff2_f;
BigReal         ComputeNonbondedUtil::dielectric_1;
const LJTable*  ComputeNonbondedUtil::ljTable = 0;
const Molecule* ComputeNonbondedUtil::mol;
BigReal		ComputeNonbondedUtil::r2_delta;
BigReal		ComputeNonbondedUtil::r2_delta_1;
int		ComputeNonbondedUtil::r2_delta_exp;
int		ComputeNonbondedUtil::rowsize;
int		ComputeNonbondedUtil::columnsize;
BigReal*	ComputeNonbondedUtil::table_alloc = 0;
BigReal*	ComputeNonbondedUtil::table_ener = 0;
BigReal*	ComputeNonbondedUtil::table_short;
BigReal*	ComputeNonbondedUtil::table_noshort;
BigReal*	ComputeNonbondedUtil::fast_table;
BigReal*	ComputeNonbondedUtil::scor_table;
BigReal*	ComputeNonbondedUtil::slow_table;
BigReal*	ComputeNonbondedUtil::corr_table;
BigReal*	ComputeNonbondedUtil::full_table;
BigReal*	ComputeNonbondedUtil::vdwa_table;
BigReal*	ComputeNonbondedUtil::vdwb_table;
BigReal*	ComputeNonbondedUtil::r2_table;
#if defined(NAMD_MIC)
  BigReal*      ComputeNonbondedUtil::mic_table_base_ptr;
  int           ComputeNonbondedUtil::mic_table_n;
  int           ComputeNonbondedUtil::mic_table_n_16;
#endif
#ifdef NAMD_KNL
float*          ComputeNonbondedUtil::knl_table_alloc;
float*          ComputeNonbondedUtil::knl_fast_ener_table;
float*          ComputeNonbondedUtil::knl_fast_grad_table;
float*          ComputeNonbondedUtil::knl_scor_ener_table;
float*          ComputeNonbondedUtil::knl_scor_grad_table;
float*          ComputeNonbondedUtil::knl_slow_ener_table;
float*          ComputeNonbondedUtil::knl_slow_grad_table;
float*          ComputeNonbondedUtil::knl_corr_ener_table;
float*          ComputeNonbondedUtil::knl_corr_grad_table;
float*          ComputeNonbondedUtil::knl_full_ener_table;
float*          ComputeNonbondedUtil::knl_full_grad_table;
#endif
BigReal         ComputeNonbondedUtil::scaling;
BigReal         ComputeNonbondedUtil::scale14;
BigReal         ComputeNonbondedUtil::switchOn;
BigReal         ComputeNonbondedUtil::switchOn_1;
BigReal         ComputeNonbondedUtil::switchOn2;
BigReal         ComputeNonbondedUtil::v_vdwa;
BigReal         ComputeNonbondedUtil::v_vdwb;
BigReal         ComputeNonbondedUtil::k_vdwa;
BigReal         ComputeNonbondedUtil::k_vdwb;
BigReal         ComputeNonbondedUtil::cutoff_3;
BigReal         ComputeNonbondedUtil::cutoff_6;
float           ComputeNonbondedUtil::v_vdwa_f;
float           ComputeNonbondedUtil::v_vdwb_f;
float           ComputeNonbondedUtil::k_vdwa_f;
float           ComputeNonbondedUtil::k_vdwb_f;
float           ComputeNonbondedUtil::cutoff_3_f;
float           ComputeNonbondedUtil::cutoff_6_f;
float           ComputeNonbondedUtil::switchOn_f;
float           ComputeNonbondedUtil::A6_f;
float           ComputeNonbondedUtil::B6_f;
float           ComputeNonbondedUtil::C6_f;
float           ComputeNonbondedUtil::A12_f;
float           ComputeNonbondedUtil::B12_f;
float           ComputeNonbondedUtil::C12_f;
BigReal         ComputeNonbondedUtil::c0;
BigReal         ComputeNonbondedUtil::c1;
BigReal         ComputeNonbondedUtil::c3;
BigReal         ComputeNonbondedUtil::c5;
BigReal         ComputeNonbondedUtil::c6;
BigReal         ComputeNonbondedUtil::c7;
BigReal         ComputeNonbondedUtil::c8;
// BigReal         ComputeNonbondedUtil::d0;
// fepb
Bool      ComputeNonbondedUtil::alchFepOn;
Bool      ComputeNonbondedUtil::alchThermIntOn;
Bool      ComputeNonbondedUtil::Fep_WCA_repuOn;
Bool      ComputeNonbondedUtil::Fep_WCA_dispOn;
Bool      ComputeNonbondedUtil::Fep_ElecOn;
Bool      ComputeNonbondedUtil::Fep_Wham; 
BigReal   ComputeNonbondedUtil::WCA_rcut1;
BigReal   ComputeNonbondedUtil::WCA_rcut2;
BigReal   ComputeNonbondedUtil::WCA_rcut3;
BigReal   ComputeNonbondedUtil::alchLambda2;
BigReal   ComputeNonbondedUtil::alchRepLambda;
BigReal   ComputeNonbondedUtil::alchDispLambda;
BigReal   ComputeNonbondedUtil::alchElecLambda;
BigReal   ComputeNonbondedUtil::alchVdwShiftCoeff;
Bool      ComputeNonbondedUtil::vdwForceSwitching;
Bool      ComputeNonbondedUtil::alchDecouple;
//fepe
Bool      ComputeNonbondedUtil::lesOn;
int       ComputeNonbondedUtil::lesFactor;
BigReal   ComputeNonbondedUtil::lesScaling;

BigReal*	ComputeNonbondedUtil::lambda_table = 0;

Bool            ComputeNonbondedUtil::pairInteractionOn;
Bool            ComputeNonbondedUtil::pairInteractionSelf;

Bool            ComputeNonbondedUtil::pressureProfileOn;
int             ComputeNonbondedUtil::pressureProfileSlabs;
int             ComputeNonbondedUtil::pressureProfileAtomTypes;
BigReal         ComputeNonbondedUtil::pressureProfileThickness;
BigReal         ComputeNonbondedUtil::pressureProfileMin;

Bool            ComputeNonbondedUtil::accelMDOn;

Bool            ComputeNonbondedUtil::drudeNbthole;

BigReal		ComputeNonbondedUtil::ewaldcof;
BigReal		ComputeNonbondedUtil::pi_ewaldcof;

int		ComputeNonbondedUtil::vdw_switch_mode;

// Ported by JLai -- JE - Go
Bool            ComputeNonbondedUtil::goGroPair;
Bool            ComputeNonbondedUtil::goForcesOn;
int             ComputeNonbondedUtil::goMethod; //6.3.11
// End of port -- JLai

void (*ComputeNonbondedUtil::calcPair)(nonbonded *);
void (*ComputeNonbondedUtil::calcPairEnergy)(nonbonded *);
void (*ComputeNonbondedUtil::calcSelf)(nonbonded *);
void (*ComputeNonbondedUtil::calcSelfEnergy)(nonbonded *);

void (*ComputeNonbondedUtil::calcFullPair)(nonbonded *);
void (*ComputeNonbondedUtil::calcFullPairEnergy)(nonbonded *);
void (*ComputeNonbondedUtil::calcFullSelf)(nonbonded *);
void (*ComputeNonbondedUtil::calcFullSelfEnergy)(nonbonded *);

void (*ComputeNonbondedUtil::calcMergePair)(nonbonded *);
void (*ComputeNonbondedUtil::calcMergePairEnergy)(nonbonded *);
void (*ComputeNonbondedUtil::calcMergeSelf)(nonbonded *);
void (*ComputeNonbondedUtil::calcMergeSelfEnergy)(nonbonded *);

void (*ComputeNonbondedUtil::calcSlowPair)(nonbonded *);
void (*ComputeNonbondedUtil::calcSlowPairEnergy)(nonbonded *);
void (*ComputeNonbondedUtil::calcSlowSelf)(nonbonded *);
void (*ComputeNonbondedUtil::calcSlowSelfEnergy)(nonbonded *);

// define splitting function
#define SPLIT_NONE	1
#define SPLIT_SHIFT	2
#define SPLIT_C1	3
#define SPLIT_XPLOR	4
#define SPLIT_C2	5
#define SPLIT_MARTINI	6

void ComputeNonbondedUtil::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_EXCLUSION_CHECKSUM) += data[exclChecksumIndex];
  reduction->item(REDUCTION_PAIRLIST_WARNINGS) += data[pairlistWarningIndex];
  reduction->item(REDUCTION_ELECT_ENERGY) += data[electEnergyIndex];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += data[fullElectEnergyIndex];
  reduction->item(REDUCTION_LJ_ENERGY) += data[vdwEnergyIndex];
  // Ported by JLai
  reduction->item(REDUCTION_GRO_LJ_ENERGY) += data[groLJEnergyIndex];
  reduction->item(REDUCTION_GRO_GAUSS_ENERGY) += data[groGaussEnergyIndex];
  reduction->item(REDUCTION_GO_NATIVE_ENERGY) += data[goNativeEnergyIndex];
  reduction->item(REDUCTION_GO_NONNATIVE_ENERGY) += data[goNonnativeEnergyIndex];
  // End of port -- JLai
//fepb
  reduction->item(REDUCTION_ELECT_ENERGY_F) += data[electEnergyIndex_s];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F) += data[fullElectEnergyIndex_s];
  reduction->item(REDUCTION_LJ_ENERGY_F) += data[vdwEnergyIndex_s];
  reduction->item(REDUCTION_LJ_ENERGY_F_LEFT) += data[vdwEnergyIndex_s_Left];

  reduction->item(REDUCTION_ELECT_ENERGY_TI_1) += data[electEnergyIndex_ti_1];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW_TI_1) += data[fullElectEnergyIndex_ti_1];
  reduction->item(REDUCTION_LJ_ENERGY_TI_1) += data[vdwEnergyIndex_ti_1];
  reduction->item(REDUCTION_ELECT_ENERGY_TI_2) += data[electEnergyIndex_ti_2];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW_TI_2) += data[fullElectEnergyIndex_ti_2];
  reduction->item(REDUCTION_LJ_ENERGY_TI_2) += data[vdwEnergyIndex_ti_2];
//fepe
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NBOND,data,virialIndex);
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_SLOW,data,fullElectVirialIndex);
  ADD_VECTOR(reduction,REDUCTION_PAIR_VDW_FORCE,data,pairVDWForceIndex);
  ADD_VECTOR(reduction,REDUCTION_PAIR_ELECT_FORCE,data,pairElectForceIndex);
  reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
}

void ComputeNonbondedUtil::submitPressureProfileData(BigReal *data,
  SubmitReduction *reduction)
{
  if (!reduction) return;
  int numAtomTypes = pressureProfileAtomTypes;
  // For ease of calculation we stored interactions between types
  // i and j in (ni+j).  For efficiency now we coalesce the
  // cross interactions so that just i<=j are stored.
  const int arraysize = 3*pressureProfileSlabs;
  size_t nelems = arraysize*(numAtomTypes*(numAtomTypes+1))/2;
  BigReal *arr = new BigReal[nelems];
  memset(arr, 0, nelems*sizeof(BigReal));

  int i, j;
  for (i=0; i<numAtomTypes; i++) {
    for (j=0; j<numAtomTypes; j++) {
      int ii=i;
      int jj=j;
      if (ii > jj) { int tmp=ii; ii=jj; jj=tmp; }
      const int reductionOffset = (ii*numAtomTypes - (ii*(ii+1))/2 + jj)*arraysize;
      for (int k=0; k<arraysize; k++) {
        arr[reductionOffset+k] += data[k];
      }
      data += arraysize;
    }
  }
  // copy into reduction
  reduction->add(nelems, arr);
  delete [] arr;
}
  
void ComputeNonbondedUtil::calc_error(nonbonded *) {
  NAMD_bug("Tried to call missing nonbonded compute routine.");
}
  
void ComputeNonbondedUtil::select(void)
{
  if ( CkMyRank() ) return;

  // These defaults die cleanly if nothing appropriate is assigned.
  ComputeNonbondedUtil::calcPair = calc_error;
  ComputeNonbondedUtil::calcPairEnergy = calc_error;
  ComputeNonbondedUtil::calcSelf = calc_error;
  ComputeNonbondedUtil::calcSelfEnergy = calc_error;
  ComputeNonbondedUtil::calcFullPair = calc_error;
  ComputeNonbondedUtil::calcFullPairEnergy = calc_error;
  ComputeNonbondedUtil::calcFullSelf = calc_error;
  ComputeNonbondedUtil::calcFullSelfEnergy = calc_error;
  ComputeNonbondedUtil::calcMergePair = calc_error;
  ComputeNonbondedUtil::calcMergePairEnergy = calc_error;
  ComputeNonbondedUtil::calcMergeSelf = calc_error;
  ComputeNonbondedUtil::calcMergeSelfEnergy = calc_error;
  ComputeNonbondedUtil::calcSlowPair = calc_error;
  ComputeNonbondedUtil::calcSlowPairEnergy = calc_error;
  ComputeNonbondedUtil::calcSlowSelf = calc_error;
  ComputeNonbondedUtil::calcSlowSelfEnergy = calc_error;

  SimParameters * simParams = Node::Object()->simParameters;
  Parameters * params = Node::Object()->parameters;

  table_ener = params->table_ener;
  rowsize = params->rowsize;
  columnsize = params->columnsize;

  commOnly = simParams->commOnly;
  fixedAtomsOn = ( simParams->fixedAtomsOn && ! simParams->fixedAtomsForces );

  qmForcesOn = simParams->qmForcesOn ;
  
  cutoff = simParams->cutoff;
  cutoff2 = cutoff*cutoff;
  cutoff2_f = cutoff2;

//fepb
  alchFepOn = simParams->alchFepOn;
  Fep_WCA_repuOn = simParams->alchFepWCARepuOn;
  Fep_WCA_dispOn = simParams->alchFepWCADispOn;
  Fep_ElecOn = simParams->alchFepElecOn;
  Fep_Wham = simParams->alchFepWhamOn;
  alchThermIntOn = simParams->alchThermIntOn;
  alchLambda2 = 0;
  lesOn = simParams->lesOn;
  lesScaling = lesFactor = 0;
  Bool tabulatedEnergies = simParams->tabulatedEnergies;
  alchVdwShiftCoeff = simParams->alchVdwShiftCoeff;
  vdwForceSwitching = simParams->vdwForceSwitching;
  WCA_rcut1 = simParams->alchFepWCArcut1;
  WCA_rcut2 = simParams->alchFepWCArcut2;
  WCA_rcut3 = simParams->alchFepWCArcut3;

  alchRepLambda = simParams->alchRepLambda;
  alchDispLambda = simParams->alchDispLambda;
  alchElecLambda = simParams->alchElecLambda;

  alchDecouple = simParams->alchDecouple;

  delete [] lambda_table;
  lambda_table = 0;

  pairInteractionOn = simParams->pairInteractionOn;
  pairInteractionSelf = simParams->pairInteractionSelf;
  pressureProfileOn = simParams->pressureProfileOn;

  // Ported by JLai -- Original JE - Go
  goGroPair = simParams->goGroPair;
  goForcesOn = simParams->goForcesOn;
  goMethod = simParams->goMethod; 
  // End of port

  accelMDOn = simParams->accelMDOn;

  drudeNbthole = simParams->drudeOn && (simParams->drudeNbtholeCut > 0.0);

  if ( drudeNbthole ) {
#ifdef NAMD_CUDA
    NAMD_die("drudeNbthole is not supported in CUDA version");
#endif
    if ( alchFepOn )
      NAMD_die("drudeNbthole is not supported with alchemical free-energy perturbation");
    if ( alchThermIntOn )
      NAMD_die("drudeNbthole is not supported with alchemical thermodynamic integration");
    if ( lesOn )
      NAMD_die("drudeNbthole is not supported with locally enhanced sampling");
    if ( pairInteractionOn )
      NAMD_die("drudeNbthole is not supported with pair interaction calculation");
    if ( pressureProfileOn )
      NAMD_die("drudeNbthole is not supported with pressure profile calculation");
  }

  if ( alchFepOn ) {
#ifdef NAMD_CUDA
    NAMD_die("Alchemical free-energy perturbation is not supported in CUDA version");
#endif
    alchLambda2 = simParams->alchLambda2;
    ComputeNonbondedUtil::calcPair = calc_pair_energy_fep;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_fep;
    ComputeNonbondedUtil::calcSelf = calc_self_energy_fep;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_fep;
    ComputeNonbondedUtil::calcFullPair = calc_pair_energy_fullelect_fep;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_fep;
    ComputeNonbondedUtil::calcFullSelf = calc_self_energy_fullelect_fep;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_fep;
    ComputeNonbondedUtil::calcMergePair = calc_pair_energy_merge_fullelect_fep;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_fep;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_energy_merge_fullelect_fep;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_fep;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_energy_slow_fullelect_fep;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect_fep;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_energy_slow_fullelect_fep;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect_fep;
  }  else if ( alchThermIntOn ) {
#ifdef NAMD_CUDA
    NAMD_die("Alchemical thermodynamic integration is not supported in CUDA version");
#endif
    alchLambda2 = simParams->alchLambda2;
    ComputeNonbondedUtil::calcPair = calc_pair_ti;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_ti;
    ComputeNonbondedUtil::calcSelf = calc_self_ti;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_ti;
    ComputeNonbondedUtil::calcFullPair = calc_pair_fullelect_ti;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_ti;
    ComputeNonbondedUtil::calcFullSelf = calc_self_fullelect_ti;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_ti;
    ComputeNonbondedUtil::calcMergePair = calc_pair_merge_fullelect_ti;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_ti;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_merge_fullelect_ti;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_ti;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_slow_fullelect_ti;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect_ti;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_slow_fullelect_ti;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect_ti;
  } else if ( lesOn ) {
#ifdef NAMD_CUDA
    NAMD_die("Locally enhanced sampling is not supported in CUDA version");
#endif
    lesFactor = simParams->lesFactor;
    lesScaling = 1.0 / (double)lesFactor;
    lambda_table = new BigReal[(lesFactor+1)*(lesFactor+1)];
    for ( int ip=0; ip<=lesFactor; ++ip ) {
      for ( int jp=0; jp<=lesFactor; ++jp ) {
        BigReal lambda_pair = 1.0;
        if (ip || jp ) {
          if (ip && jp && ip != jp) {
            lambda_pair = 0.0;
          } else {
            lambda_pair = lesScaling;
          }
        }
        lambda_table[(lesFactor+1)*ip+jp] = lambda_pair;
      }
    }
    ComputeNonbondedUtil::calcPair = calc_pair_les;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_les;
    ComputeNonbondedUtil::calcSelf = calc_self_les;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_les;
    ComputeNonbondedUtil::calcFullPair = calc_pair_fullelect_les;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_les;
    ComputeNonbondedUtil::calcFullSelf = calc_self_fullelect_les;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_les;
    ComputeNonbondedUtil::calcMergePair = calc_pair_merge_fullelect_les;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_les;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_merge_fullelect_les;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_les;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_slow_fullelect_les;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect_les;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_slow_fullelect_les;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect_les;
  } else if ( pressureProfileOn) {
#ifdef NAMD_CUDA
    NAMD_die("Pressure profile calculation is not supported in CUDA version");
#endif
    pressureProfileSlabs = simParams->pressureProfileSlabs;
    pressureProfileAtomTypes = simParams->pressureProfileAtomTypes;

    ComputeNonbondedUtil::calcPair = calc_pair_pprof;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_pprof;
    ComputeNonbondedUtil::calcSelf = calc_self_pprof;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_pprof;
    ComputeNonbondedUtil::calcFullPair = calc_pair_fullelect_pprof;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_pprof;
    ComputeNonbondedUtil::calcFullSelf = calc_self_fullelect_pprof;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_pprof;
    ComputeNonbondedUtil::calcMergePair = calc_pair_merge_fullelect_pprof;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_pprof;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_merge_fullelect_pprof;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_pprof;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_slow_fullelect_pprof;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect_pprof;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_slow_fullelect_pprof;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect_pprof;
  } else if ( pairInteractionOn ) {
#ifdef NAMD_CUDA
    NAMD_die("Pair interaction calculation is not supported in CUDA version");
#endif
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_int;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_int;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_int;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_int;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_int;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_int;
  } else if ( tabulatedEnergies ) {
#ifdef NAMD_CUDA
    NAMD_die("Tabulated energies is not supported in CUDA version");
#endif
    ComputeNonbondedUtil::calcPair = calc_pair_tabener;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_tabener;
    ComputeNonbondedUtil::calcSelf = calc_self_tabener;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_tabener;
    ComputeNonbondedUtil::calcFullPair = calc_pair_fullelect_tabener;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_tabener;
    ComputeNonbondedUtil::calcFullSelf = calc_self_fullelect_tabener;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_tabener;
    ComputeNonbondedUtil::calcMergePair = calc_pair_merge_fullelect_tabener;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_tabener;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_merge_fullelect_tabener;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_tabener;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_slow_fullelect_tabener;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect_tabener;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_slow_fullelect_tabener;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect_tabener;
  } else if ( goForcesOn ) {
#ifdef NAMD_CUDA
    NAMD_die("Go forces is not supported in CUDA version");
#endif
    ComputeNonbondedUtil::calcPair = calc_pair_go;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy_go;
    ComputeNonbondedUtil::calcSelf = calc_self_go;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy_go;
    ComputeNonbondedUtil::calcFullPair = calc_pair_fullelect_go;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect_go;
    ComputeNonbondedUtil::calcFullSelf = calc_self_fullelect_go;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect_go;
    ComputeNonbondedUtil::calcMergePair = calc_pair_merge_fullelect_go;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect_go;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_merge_fullelect_go;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect_go;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_slow_fullelect_go;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect_go;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_slow_fullelect_go;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect_go;
  } else {
    ComputeNonbondedUtil::calcPair = calc_pair;
    ComputeNonbondedUtil::calcPairEnergy = calc_pair_energy;
    ComputeNonbondedUtil::calcSelf = calc_self;
    ComputeNonbondedUtil::calcSelfEnergy = calc_self_energy;
    ComputeNonbondedUtil::calcFullPair = calc_pair_fullelect;
    ComputeNonbondedUtil::calcFullPairEnergy = calc_pair_energy_fullelect;
    ComputeNonbondedUtil::calcFullSelf = calc_self_fullelect;
    ComputeNonbondedUtil::calcFullSelfEnergy = calc_self_energy_fullelect;
    ComputeNonbondedUtil::calcMergePair = calc_pair_merge_fullelect;
    ComputeNonbondedUtil::calcMergePairEnergy = calc_pair_energy_merge_fullelect;
    ComputeNonbondedUtil::calcMergeSelf = calc_self_merge_fullelect;
    ComputeNonbondedUtil::calcMergeSelfEnergy = calc_self_energy_merge_fullelect;
    ComputeNonbondedUtil::calcSlowPair = calc_pair_slow_fullelect;
    ComputeNonbondedUtil::calcSlowPairEnergy = calc_pair_energy_slow_fullelect;
    ComputeNonbondedUtil::calcSlowSelf = calc_self_slow_fullelect;
    ComputeNonbondedUtil::calcSlowSelfEnergy = calc_self_energy_slow_fullelect;
  }

//fepe

  dielectric_1 = 1.0/simParams->dielectric;
  if ( ! ljTable ) ljTable = new LJTable;
  mol = Node::Object()->molecule;
  scaling = simParams->nonbondedScaling;
  if ( simParams->exclude == SCALED14 )
  {
    scale14 = simParams->scale14;
  }
  else
  {
    scale14 = 1.;
  }
  if ( simParams->switchingActive )
  {
    switchOn = simParams->switchingDist;
    switchOn_1 = 1.0/switchOn;
    // d0 = 1.0/(cutoff-switchOn);
    switchOn2 = switchOn*switchOn;
    c0 = 1.0/(cutoff2-switchOn2);

    if ( simParams->vdwForceSwitching ) {
      double switchOn3 = switchOn * switchOn2;
      double cutoff3 = cutoff * cutoff2;
      double switchOn6 = switchOn3 * switchOn3;
      double cutoff6 = cutoff3 * cutoff3;
      v_vdwa_f = v_vdwa = -1. / ( switchOn6 * cutoff6 );
      v_vdwb_f = v_vdwb = -1. / ( switchOn3 * cutoff3 );
      k_vdwa_f = k_vdwa = cutoff6 / ( cutoff6 - switchOn6 );
      k_vdwb_f = k_vdwb = cutoff3 / ( cutoff3 - switchOn3 );
      cutoff_3_f = cutoff_3 = 1. / cutoff3;
      cutoff_6_f = cutoff_6 = 1. / cutoff6;

    } else if ( simParams->martiniSwitching ) { // switching fxn for Martini RBCG

      BigReal p6 = 6;
      BigReal A6 = p6 * ((p6+1)*switchOn-(p6+4)*cutoff)/(pow(cutoff,p6+2)*pow(cutoff-switchOn,2));
      BigReal B6 = -p6 * ((p6+1)*switchOn-(p6+3)*cutoff)/(pow(cutoff,p6+2)*pow(cutoff-switchOn,3));        
      BigReal C6 = 1.0/pow(cutoff,p6)-A6/3.0*pow(cutoff-switchOn,3)-B6/4.0*pow(cutoff-switchOn,4);

      BigReal p12 = 12;
      BigReal A12 = p12 * ((p12+1)*switchOn-(p12+4)*cutoff)/(pow(cutoff,p12+2)*pow(cutoff-switchOn,2));
      BigReal B12 = -p12 * ((p12+1)*switchOn-(p12+3)*cutoff)/(pow(cutoff,p12+2)*pow(cutoff-switchOn,3));
      BigReal C12 = 1.0/pow(cutoff,p12)-A12/3.0*pow(cutoff-switchOn,3)-B12/4.0*pow(cutoff-switchOn,4);

      A6_f =  A6;  B6_f  = B6;  C6_f =  C6;
      A12_f = A12; B12_f = B12; C12_f = C12;
      switchOn_f = switchOn;

    }

  }
  else
  {
    switchOn = cutoff;
    switchOn_1 = 1.0/switchOn;
    // d0 = 0.;  // avoid division by zero
    switchOn2 = switchOn*switchOn;
    c0 = 0.;  // avoid division by zero
  }
  c1 = c0*c0*c0;
  c3 = 3.0 * (cutoff2 - switchOn2);
  c5 = 0;
  c6 = 0;
  c7 = 0;
  c8 = 0;

  const int PMEOn = simParams->PMEOn;
  const int MSMOn = simParams->MSMOn;
  const int MSMSplit = simParams->MSMSplit;

  if ( PMEOn ) {
    ewaldcof = simParams->PMEEwaldCoefficient;
    BigReal TwoBySqrtPi = 1.12837916709551;
    pi_ewaldcof = TwoBySqrtPi * ewaldcof;
  }

  int splitType = SPLIT_NONE;
  if ( simParams->switchingActive ) splitType = SPLIT_SHIFT;
  if ( simParams->martiniSwitching ) splitType = SPLIT_MARTINI;
  if ( simParams->fullDirectOn || simParams->FMAOn || PMEOn || MSMOn ||
      simParams->FMMOn ) {
    switch ( simParams->longSplitting ) {
      case C2:
      splitType = SPLIT_C2;
      break;

      case C1:
      splitType = SPLIT_C1;
      break;

      case XPLOR:
      NAMD_die("Sorry, XPLOR splitting not supported.");
      break;

      case SHARP:
      NAMD_die("Sorry, SHARP splitting not supported.");
      break;

      default:
      NAMD_die("Unknown splitting type found!");

    }
  }

  BigReal r2_tol = 0.1;
  
  r2_delta = 1.0;
  r2_delta_exp = 0;
  while ( r2_delta > r2_tol ) { r2_delta /= 2.0; r2_delta_exp += 1; }
  r2_delta_1 = 1.0 / r2_delta;

  if ( ! CkMyPe() ) {
    iout << iINFO << "NONBONDED TABLE R-SQUARED SPACING: " <<
				r2_delta << "\n" << endi;
  }

  BigReal r2_tmp = 1.0;
  int cutoff2_exp = 0;
  while ( (cutoff2 + r2_delta) > r2_tmp ) { r2_tmp *= 2.0; cutoff2_exp += 1; }

  int i;
  int n = (r2_delta_exp + cutoff2_exp) * 64 + 1;
  #if defined(NAMD_MIC)
    int n_16 = (n + 15) & (~15);
  #endif

  if ( ! CkMyPe() ) {
    iout << iINFO << "NONBONDED TABLE SIZE: " <<
				n << " POINTS\n" << endi;
  }

  if ( table_alloc ) delete [] table_alloc;
  #if defined(NAMD_MIC)
    table_alloc = new BigReal[61*n_16+16];
    BigReal *table_align = table_alloc;
    while ( ((long)table_align) % 128 ) ++table_align;
    mic_table_base_ptr = table_align;
    mic_table_n = n;
    mic_table_n_16 = n_16;
    table_noshort = table_align;
    table_short = table_align + 16*n_16;
    slow_table = table_align + 32*n_16;
    fast_table = table_align + 36*n_16;
    scor_table = table_align + 40*n_16;
    corr_table = table_align + 44*n_16;
    full_table = table_align + 48*n_16;
    vdwa_table = table_align + 52*n_16;
    vdwb_table = table_align + 56*n_16;
    r2_table = table_align + 60*n_16;
  #else
  table_alloc = new BigReal[61*n+16];
  BigReal *table_align = table_alloc;
  while ( ((long)table_align) % 128 ) ++table_align;
  table_noshort = table_align;
  table_short = table_align + 16*n;
  slow_table = table_align + 32*n;
  fast_table = table_align + 36*n;
  scor_table = table_align + 40*n;
  corr_table = table_align + 44*n;
  full_table = table_align + 48*n;
  vdwa_table = table_align + 52*n;
  vdwb_table = table_align + 56*n;
  r2_table = table_align + 60*n;
  #endif
  BigReal *fast_i = fast_table + 4;
  BigReal *scor_i = scor_table + 4;
  BigReal *slow_i = slow_table + 4;
  BigReal *vdwa_i = vdwa_table + 4;
  BigReal *vdwb_i = vdwb_table + 4;
  BigReal *r2_i = r2_table;  *(r2_i++) = r2_delta;
  BigReal r2_limit = simParams->limitDist * simParams->limitDist;
  if ( r2_limit < r2_delta ) r2_limit = r2_delta;
  int r2_delta_i = 0;  // entry for r2 == r2_delta

#ifdef NAMD_KNL
 if ( knl_table_alloc ) delete [] knl_table_alloc;
 knl_table_alloc = new float[10*KNL_TABLE_SIZE];
 knl_fast_ener_table = knl_table_alloc;
 knl_fast_grad_table = knl_table_alloc + KNL_TABLE_SIZE;
 knl_scor_ener_table = knl_table_alloc + 2*KNL_TABLE_SIZE;
 knl_scor_grad_table = knl_table_alloc + 3*KNL_TABLE_SIZE;
 knl_slow_ener_table = knl_table_alloc + 4*KNL_TABLE_SIZE;
 knl_slow_grad_table = knl_table_alloc + 5*KNL_TABLE_SIZE;
 knl_corr_ener_table = knl_table_alloc + 6*KNL_TABLE_SIZE;
 knl_corr_grad_table = knl_table_alloc + 7*KNL_TABLE_SIZE;
 knl_full_ener_table = knl_table_alloc + 8*KNL_TABLE_SIZE;
 knl_full_grad_table = knl_table_alloc + 9*KNL_TABLE_SIZE;
 knl_fast_ener_table[0] = 0.;
 knl_fast_grad_table[0] = 0.;
 knl_scor_ener_table[0] = 0.;
 knl_scor_grad_table[0] = 0.;
 knl_slow_ener_table[0] = 0.;
 knl_slow_grad_table[0] = 0.;
 knl_corr_ener_table[0] = 0.;
 knl_corr_grad_table[0] = 0.;
 knl_full_ener_table[0] = 0.;
 knl_full_grad_table[0] = 0.;
 for ( int knl_table = 0; knl_table < 2; ++knl_table ) {
  int nn = n;
  if ( knl_table ) {
    nn = KNL_TABLE_SIZE-1;
  }
  for ( i=1; i<nn; ++i ) {
#else
  // fill in the table, fix up i==0 (r2==0) below
  for ( i=1; i<n; ++i ) {
#endif

    const BigReal r2_base = r2_delta * ( 1 << (i/64) );
    const BigReal r2_del = r2_base / 64.0;
    BigReal r2 = r2_base - r2_delta + r2_del * (i%64);

    BigReal r = sqrt(r2);

#ifdef NAMD_KNL
    if ( knl_table ) {
      r = (double)(nn-1)/(double)(i);
      r2 = r*r;
    } else
#endif
    if ( r2 <= r2_limit ) r2_delta_i = i;

    const BigReal r_1 = 1.0/r;
    const BigReal r_2 = 1.0/r2;

    // fast_ is defined as (full_ - slow_)
    // corr_ and fast_ are both zero at the cutoff, full_ is not
    // all three are approx 1/r at short distances

    // for actual interpolation, we use fast_ for fast forces and
    // scor_ = slow_ + corr_ - full_ and slow_ for slow forces
    // since these last two are of small magnitude

    BigReal fast_energy, fast_gradient;
    BigReal scor_energy, scor_gradient;
    BigReal slow_energy, slow_gradient;

    // corr_ is PME direct sum, or similar correction term
    // corr_energy is multiplied by r until later
    // corr_gradient is multiplied by -r^2 until later
    BigReal corr_energy, corr_gradient;

    
    if ( PMEOn ) {
      BigReal tmp_a = r * ewaldcof;
      BigReal tmp_b = erfc(tmp_a);
      corr_energy = tmp_b;
      corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
    } else if ( MSMOn ) {
      BigReal a_1 = 1.0/cutoff;
      BigReal r_a = r * a_1;
      BigReal g, dg;
      SPOLY(&g, &dg, r_a, MSMSplit);
      corr_energy = 1 - r_a * g;
      corr_gradient = 1 + r_a*r_a * dg;
    } else {
      corr_energy = corr_gradient = 0;
    }

    switch(splitType) {
      case SPLIT_NONE:
        fast_energy = 1.0/r;
        fast_gradient = -1.0/r2;
        scor_energy = scor_gradient = 0;
        slow_energy = slow_gradient = 0;
	break;
      case SPLIT_SHIFT: {
	BigReal shiftVal = r2/cutoff2 - 1.0;
	shiftVal *= shiftVal;
	BigReal dShiftVal = 2.0 * (r2/cutoff2 - 1.0) * 2.0*r/cutoff2;
        fast_energy = shiftVal/r;
        fast_gradient = dShiftVal/r - shiftVal/r2;
        scor_energy = scor_gradient = 0;
        slow_energy = slow_gradient = 0;
        } 
	break;
      case SPLIT_MARTINI: { 
        // in Martini, the Coulomb switching distance is zero
        const BigReal COUL_SWITCH = 0.;
        // Gromacs shifting function
        const BigReal p1 = 1.;
        BigReal A1 = p1 * ((p1+1)*COUL_SWITCH-(p1+4)*cutoff)/(pow(cutoff,p1+2)*pow(cutoff-COUL_SWITCH,2));
        BigReal B1 = -p1 * ((p1+1)*COUL_SWITCH-(p1+3)*cutoff)/(pow(cutoff,p1+2)*pow(cutoff-COUL_SWITCH,3));
        BigReal X1 = 1.0/pow(cutoff,p1)-A1/3.0*pow(cutoff-COUL_SWITCH,3)-B1/4.0*pow(cutoff-COUL_SWITCH,4);
        BigReal r12 = (r-COUL_SWITCH)*(r-COUL_SWITCH);
        BigReal r13 = (r-COUL_SWITCH)*(r-COUL_SWITCH)*(r-COUL_SWITCH);
        BigReal shiftVal = -(A1/3.0)*r13 - (B1/4.0)*r12*r12 - X1;
        BigReal dShiftVal = -A1*r12 - B1*r13;
        fast_energy = (1/r) + shiftVal;
        fast_gradient = -1/(r2) + dShiftVal;
        scor_energy = scor_gradient = 0;
        slow_energy = slow_gradient = 0;
        } 
	break;
      case SPLIT_C1:
	// calculate actual energy and gradient
	slow_energy = 0.5/cutoff * (3.0 - (r2/cutoff2));
	slow_gradient = -1.0/cutoff2 * (r/cutoff);
	// calculate scor from slow and corr
	scor_energy = slow_energy + (corr_energy - 1.0)/r;
	scor_gradient = slow_gradient - (corr_gradient - 1.0)/r2;
	// calculate fast from slow
	fast_energy = 1.0/r - slow_energy;
	fast_gradient = -1.0/r2 - slow_gradient;
	break;
      case SPLIT_C2:
        //
        // Quintic splitting function contributed by
        // Bruce Berne, Ruhong Zhou, and Joe Morrone
        //
	// calculate actual energy and gradient
        slow_energy = r2/(cutoff*cutoff2) * (6.0 * (r2/cutoff2)
            - 15.0*(r/cutoff) + 10.0);
        slow_gradient = r/(cutoff*cutoff2) * (24.0 * (r2/cutoff2)
            - 45.0 *(r/cutoff) + 20.0);
	// calculate scor from slow and corr
        scor_energy = slow_energy + (corr_energy - 1.0)/r;
        scor_gradient = slow_gradient - (corr_gradient - 1.0)/r2;
	// calculate fast from slow
	fast_energy = 1.0/r - slow_energy;
	fast_gradient = -1.0/r2 - slow_gradient;
	break;
    }

    // foo_gradient is calculated as ( d foo_energy / d r )
    // and now divided by 2r to get ( d foo_energy / d r2 )

    fast_gradient *= 0.5 * r_1;
    scor_gradient *= 0.5 * r_1;
    slow_gradient *= 0.5 * r_1;

    // let modf be 1 if excluded, 1-scale14 if modified, 0 otherwise,
    // add scor_ - modf * slow_ to slow terms and
    // add fast_ - modf * fast_ to fast terms.

    BigReal vdwa_energy, vdwa_gradient;
    BigReal vdwb_energy, vdwb_gradient;

    const BigReal r_6 = r_2*r_2*r_2;
    const BigReal r_12 = r_6*r_6;

    // Lennard-Jones switching function
  if ( simParams->vdwForceSwitching ) {  // switch force
    vdw_switch_mode = VDW_SWITCH_MODE_FORCE;

    // from Steinbach & Brooks, JCC 15, pgs 667-683, 1994, eqns 10-13
    if ( r2 > switchOn2 ) {
      BigReal tmpa = r_6 - cutoff_6;
      vdwa_energy = k_vdwa * tmpa * tmpa;
      BigReal tmpb = r_1 * r_2 - cutoff_3;
      vdwb_energy = k_vdwb * tmpb * tmpb;
      vdwa_gradient = -6.0 * k_vdwa * tmpa * r_2 * r_6;
      vdwb_gradient = -3.0 * k_vdwb * tmpb * r_2 * r_2 * r_1;
    } else {
      vdwa_energy = r_12 + v_vdwa;
      vdwb_energy = r_6 + v_vdwb;
      vdwa_gradient = -6.0 * r_2 * r_12;
      vdwb_gradient = -3.0 * r_2 * r_6;
    }
  } else if ( simParams->martiniSwitching ) { // switching fxn for Martini RBCG
    vdw_switch_mode = VDW_SWITCH_MODE_MARTINI;

    BigReal r12 = (r-switchOn)*(r-switchOn);        BigReal r13 = (r-switchOn)*(r-switchOn)*(r-switchOn);

    BigReal p6 = 6;
    BigReal A6 = p6 * ((p6+1)*switchOn-(p6+4)*cutoff)/(pow(cutoff,p6+2)*pow(cutoff-switchOn,2));
    BigReal B6 = -p6 * ((p6+1)*switchOn-(p6+3)*cutoff)/(pow(cutoff,p6+2)*pow(cutoff-switchOn,3));        
    BigReal C6 = 1.0/pow(cutoff,p6)-A6/3.0*pow(cutoff-switchOn,3)-B6/4.0*pow(cutoff-switchOn,4);

    BigReal p12 = 12;
    BigReal A12 = p12 * ((p12+1)*switchOn-(p12+4)*cutoff)/(pow(cutoff,p12+2)*pow(cutoff-switchOn,2));
    BigReal B12 = -p12 * ((p12+1)*switchOn-(p12+3)*cutoff)/(pow(cutoff,p12+2)*pow(cutoff-switchOn,3));
    BigReal C12 = 1.0/pow(cutoff,p12)-A12/3.0*pow(cutoff-switchOn,3)-B12/4.0*pow(cutoff-switchOn,4);

    BigReal LJshifttempA = -(A12/3)*r13 - (B12/4)*r12*r12 - C12;
    BigReal LJshifttempB = -(A6/3)*r13 - (B6/4)*r12*r12 - C6;
    const BigReal shiftValA =         // used for Lennard-Jones
                        ( r2 > switchOn2 ? LJshifttempA : -C12);
    const BigReal shiftValB =         // used for Lennard-Jones
                        ( r2 > switchOn2 ? LJshifttempB : -C6);

    BigReal LJdshifttempA = -A12*r12 - B12*r13;
    BigReal LJdshifttempB = -A6*r12 - B6*r13;
    const BigReal dshiftValA =         // used for Lennard-Jones
                        ( r2 > switchOn2 ? LJdshifttempA*0.5*r_1 : 0 );
    const BigReal dshiftValB =         // used for Lennard-Jones
                        ( r2 > switchOn2 ? LJdshifttempB*0.5*r_1 : 0 );




    //have not addressed r > cutoff

    //  dshiftValA*= 0.5*r_1;
    //  dshiftValB*= 0.5*r_1;

    vdwa_energy = r_12 + shiftValA;
    vdwb_energy = r_6 + shiftValB;
   
    vdwa_gradient = -6/pow(r,14) + dshiftValA ;
    vdwb_gradient = -3/pow(r,8) + dshiftValB;

  } else {  // switch energy
    vdw_switch_mode = VDW_SWITCH_MODE_ENERGY;

    const BigReal c2 = cutoff2-r2;
    const BigReal c4 = c2*(c3-2.0*c2);
    const BigReal switchVal =         // used for Lennard-Jones
                        ( r2 > switchOn2 ? c2*c4*c1 : 1.0 );
    const BigReal dSwitchVal =        // d switchVal / d r2
                        ( r2 > switchOn2 ? 2*c1*(c2*c2-c4) : 0.0 );

    vdwa_energy = switchVal * r_12;
    vdwb_energy = switchVal * r_6;

    vdwa_gradient = ( dSwitchVal - 6.0 * switchVal * r_2 ) * r_12;
    vdwb_gradient = ( dSwitchVal - 3.0 * switchVal * r_2 ) * r_6;
  }


#ifdef NAMD_KNL
   if ( knl_table ) {
    knl_fast_ener_table[i] = -1.*fast_energy;
    knl_fast_grad_table[i] = -2.*fast_gradient;
    knl_scor_ener_table[i] = -1.*scor_energy;
    knl_scor_grad_table[i] = -2.*scor_gradient;
    knl_slow_ener_table[i] = -1.*slow_energy;
    knl_slow_grad_table[i] = -2.*slow_gradient;
    knl_corr_ener_table[i] = -1.*(fast_energy + scor_energy);
    knl_corr_grad_table[i] = -2.*(fast_gradient + scor_gradient);
    knl_full_ener_table[i] = -1.*(fast_energy + slow_energy);
    knl_full_grad_table[i] = -2.*(fast_gradient + slow_gradient);
    if ( i == nn-1 ) {
      knl_fast_ener_table[nn] = knl_fast_ener_table[i];
      knl_fast_grad_table[nn] = knl_fast_grad_table[i];
      knl_scor_ener_table[nn] = knl_scor_ener_table[i];
      knl_scor_grad_table[nn] = knl_scor_grad_table[i];
      knl_slow_ener_table[nn] = knl_slow_ener_table[i];
      knl_slow_grad_table[nn] = knl_slow_grad_table[i];
      knl_corr_ener_table[nn] = knl_corr_ener_table[i];
      knl_corr_grad_table[nn] = knl_corr_grad_table[i];
      knl_full_ener_table[nn] = knl_full_ener_table[i];
      knl_full_grad_table[nn] = knl_full_grad_table[i];
    }
   } else {
#endif
    *(fast_i++) = fast_energy;
    *(fast_i++) = fast_gradient;
    *(fast_i++) = 0;
    *(fast_i++) = 0;
    *(scor_i++) = scor_energy;
    *(scor_i++) = scor_gradient;
    *(scor_i++) = 0;
    *(scor_i++) = 0;
    *(slow_i++) = slow_energy;
    *(slow_i++) = slow_gradient;
    *(slow_i++) = 0;
    *(slow_i++) = 0;
    *(vdwa_i++) = vdwa_energy;
    *(vdwa_i++) = vdwa_gradient;
    *(vdwa_i++) = 0;
    *(vdwa_i++) = 0;
    *(vdwb_i++) = vdwb_energy;
    *(vdwb_i++) = vdwb_gradient;
    *(vdwb_i++) = 0;
    *(vdwb_i++) = 0;
    *(r2_i++) = r2 + r2_delta;
#ifdef NAMD_KNL
   }
#endif

  }
#ifdef NAMD_KNL
 } // knl_table loop
#endif

  if ( ! r2_delta_i ) {
    NAMD_bug("Failed to find table entry for r2 == r2_limit\n");
  }
  if ( r2_table[r2_delta_i] > r2_limit + r2_delta ) {
    NAMD_bug("Found bad table entry for r2 == r2_limit\n");
  }

  int j;
  const char *table_name = "XXXX";
  int smooth_short = 0;
  for ( j=0; j<5; ++j ) {
    BigReal *t0 = 0;
    switch (j) {
      case 0: 
        t0 = fast_table;
        table_name = "FAST";
        smooth_short = 1;
      break;
      case 1: 
        t0 = scor_table;
        table_name = "SCOR";
        smooth_short = 0;
      break;
      case 2: 
        t0 = slow_table;
        table_name = "SLOW";
        smooth_short = 0;
      break;
      case 3: 
        t0 = vdwa_table;
        table_name = "VDWA";
        smooth_short = 1;
      break;
      case 4: 
        t0 = vdwb_table;
        table_name = "VDWB";
        smooth_short = 1;
      break;
    }
    // patch up data for i=0
    t0[0] = t0[4] - t0[5] * ( r2_delta / 64.0 );  // energy
    t0[1] = t0[5];  // gradient
    t0[2] = 0;
    t0[3] = 0;
    if ( smooth_short ) {
      BigReal energy0 = t0[4*r2_delta_i];
      BigReal gradient0 = t0[4*r2_delta_i+1];
      BigReal r20 = r2_table[r2_delta_i];
      t0[0] = energy0 - gradient0 * (r20 - r2_table[0]);  // energy
      t0[1] = gradient0;  // gradient
    }
    BigReal *t;
    for ( i=0,t=t0; i<(n-1); ++i,t+=4 ) {
      BigReal x = ( r2_delta * ( 1 << (i/64) ) ) / 64.0;
      if ( r2_table[i+1] != r2_table[i] + x ) {
        NAMD_bug("Bad table delta calculation.\n");
      }
      if ( smooth_short && i+1 < r2_delta_i ) {
        BigReal energy0 = t0[4*r2_delta_i];
        BigReal gradient0 = t0[4*r2_delta_i+1];
        BigReal r20 = r2_table[r2_delta_i];
        t[4] = energy0 - gradient0 * (r20 - r2_table[i+1]);  // energy
        t[5] = gradient0;  // gradient
      }
      BigReal v1 = t[0];
      BigReal g1 = t[1];
      BigReal v2 = t[4];
      BigReal g2 = t[5];
      // explicit formulas for v1 + g1 x + c x^2 + d x^3
      BigReal c = ( 3.0 * (v2 - v1) - x * (2.0 * g1 + g2) ) / ( x * x );
      BigReal d = ( -2.0 * (v2 - v1) + x * (g1 + g2) ) / ( x * x * x );
      // since v2 - v1 is imprecise, we refine c and d numerically
      // important because we need accurate forces (more than energies!)
      for ( int k=0; k < 2; ++k ) {
        BigReal dv = (v1 - v2) + ( ( d * x + c ) * x + g1 ) * x;
        BigReal dg = (g1 - g2) + ( 3.0 * d * x + 2.0 * c ) * x;
        c -= ( 3.0 * dv - x * dg ) / ( x * x );
        d -= ( -2.0 * dv + x * dg ) / ( x * x * x );
      }
      // store in the array;
      t[2] = c;  t[3] = d;
    }

    if ( ! CkMyPe() ) {
    BigReal dvmax = 0;
    BigReal dgmax = 0;
    BigReal dvmax_r = 0;
    BigReal dgmax_r = 0;
    BigReal fdvmax = 0;
    BigReal fdgmax = 0;
    BigReal fdvmax_r = 0;
    BigReal fdgmax_r = 0;
    BigReal dgcdamax = 0;
    BigReal dgcdimax = 0;
    BigReal dgcaimax = 0;
    BigReal dgcdamax_r = 0;
    BigReal dgcdimax_r = 0;
    BigReal dgcaimax_r = 0;
    BigReal fdgcdamax = 0;
    BigReal fdgcdimax = 0;
    BigReal fdgcaimax = 0;
    BigReal fdgcdamax_r = 0;
    BigReal fdgcdimax_r = 0;
    BigReal fdgcaimax_r = 0;
    BigReal gcm = fabs(t0[1]);  // gradient magnitude running average
    for ( i=0,t=t0; i<(n-1); ++i,t+=4 ) {
      const BigReal r2_base = r2_delta * ( 1 << (i/64) );
      const BigReal r2_del = r2_base / 64.0;
      const BigReal r2 = r2_base - r2_delta + r2_del * (i%64);
      const BigReal r = sqrt(r2);
      if ( r > cutoff ) break;
      BigReal x = r2_del;
      BigReal dv = ( ( t[3] * x + t[2] ) * x + t[1] ) * x + t[0] - t[4];
      BigReal dg = ( 3.0 * t[3] * x + 2.0 * t[2] ) * x + t[1] - t[5];
      if ( t[4] != 0. && fabs(dv/t[4]) > fdvmax ) {
        fdvmax = fabs(dv/t[4]); fdvmax_r = r;
      }
      if ( fabs(dv) > dvmax ) {
        dvmax = fabs(dv); dvmax_r = r;
      }
      if ( t[5] != 0. && fabs(dg/t[5]) > fdgmax ) {
        fdgmax = fabs(dg/t[5]); fdgmax_r = r;
      }
      if ( fabs(dg) > dgmax ) {
        dgmax = fabs(dg); dgmax_r = r;
      }
      BigReal gcd = (t[4] - t[0]) / x;  // centered difference gradient
      BigReal gcd_prec = (fabs(t[0]) + fabs(t[4])) * 1.e-15 / x;  // roundoff
      gcm = 0.9 * gcm + 0.1 * fabs(t[5]);  // magnitude running average
      BigReal gca = 0.5  * (t[1] + t[5]);  // centered average gradient
      BigReal gci = ( 0.75 * t[3] * x + t[2] ) * x + t[1];  // interpolated
      BigReal rc = sqrt(r2 + 0.5 * x);
      BigReal dgcda = gcd - gca;
      if ( dgcda != 0. && fabs(dgcda) < gcd_prec ) {
        // CkPrintf("ERROR %g < PREC %g AT %g AVG VAL %g\n", dgcda, gcd_prec, rc, gca);
        dgcda = 0.;
      }
      BigReal dgcdi = gcd - gci;
      if ( dgcdi != 0. && fabs(dgcdi) < gcd_prec ) {
        // CkPrintf("ERROR %g < PREC %g AT %g INT VAL %g\n", dgcdi, gcd_prec, rc, gci);
        dgcdi = 0.;
      }
      BigReal dgcai = gca - gci;
      if ( t[1]*t[5] > 0. && gcm != 0. && fabs(dgcda/gcm) > fdgcdamax ) {
        fdgcdamax = fabs(dgcda/gcm); fdgcdamax_r = rc;
      }
      if ( fabs(dgcda) > fdgcdamax ) {
        dgcdamax = fabs(dgcda); dgcdamax_r = rc;
      }
      if ( t[1]*t[5] > 0. && gcm != 0. && fabs(dgcdi/gcm) > fdgcdimax ) {
        fdgcdimax = fabs(dgcdi/gcm); fdgcdimax_r = rc;
      }
      if ( fabs(dgcdi) > fdgcdimax ) {
        dgcdimax = fabs(dgcdi); dgcdimax_r = rc;
      }
      if ( t[1]*t[5] > 0. && gcm != 0. && fabs(dgcai/gcm) > fdgcaimax ) {
        fdgcaimax = fabs(dgcai/gcm); fdgcaimax_r = rc;
      }
      if ( fabs(dgcai) > fdgcaimax ) {
        dgcaimax = fabs(dgcai); dgcaimax_r = rc;
      }
#if 0
      CkPrintf("TABLE %s %g %g %g %g\n",table_name,rc,dgcda/gcm,dgcda,gci);
      if (dv != 0.) CkPrintf("TABLE %d ENERGY ERROR %g AT %g (%d)\n",j,dv,r,i);
      if (dg != 0.) CkPrintf("TABLE %d FORCE ERROR %g AT %g (%d)\n",j,dg,r,i);
#endif
    }
    if ( dvmax != 0.0 ) {
      iout << iINFO << "ABSOLUTE IMPRECISION IN " << table_name <<
        " TABLE ENERGY: " << dvmax << " AT " << dvmax_r << "\n" << endi;
    }
    if ( fdvmax != 0.0 ) {
      iout << iINFO << "RELATIVE IMPRECISION IN " << table_name <<
        " TABLE ENERGY: " << fdvmax << " AT " << fdvmax_r << "\n" << endi;
    }
    if ( dgmax != 0.0 ) {
      iout << iINFO << "ABSOLUTE IMPRECISION IN " << table_name <<
        " TABLE FORCE: " << dgmax << " AT " << dgmax_r << "\n" << endi;
    }
    if ( fdgmax != 0.0 ) {
      iout << iINFO << "RELATIVE IMPRECISION IN " << table_name <<
        " TABLE FORCE: " << fdgmax << " AT " << fdgmax_r << "\n" << endi;
    }
    if (fdgcdamax != 0.0 ) {
      iout << iINFO << "INCONSISTENCY IN " << table_name <<
        " TABLE ENERGY VS FORCE: " << fdgcdamax << " AT " << fdgcdamax_r << "\n" << endi;
      if ( fdgcdamax > 0.1 ) {
        iout << iERROR << "\n";
        iout << iERROR << "CALCULATED " << table_name <<
          " FORCE MAY NOT MATCH ENERGY! POSSIBLE BUG!\n";
        iout << iERROR << "\n";
      }
    }
    if (0 && fdgcdimax != 0.0 ) {
      iout << iINFO << "INCONSISTENCY IN " << table_name <<
        " TABLE ENERGY VS FORCE: " << fdgcdimax << " AT " << fdgcdimax_r << "\n" << endi;
    }
    if ( 0 && fdgcaimax != 0.0 ) {
      iout << iINFO << "INCONSISTENCY IN " << table_name <<
        " TABLE AVG VS INT FORCE: " << fdgcaimax << " AT " << fdgcaimax_r << "\n" << endi;
    }
    }

  }

  for ( i=0; i<4*n; ++i ) {
    corr_table[i] = fast_table[i] + scor_table[i];
    full_table[i] = fast_table[i] + slow_table[i];
  }

#if 0  
  for ( i=0; i<n; ++i ) {
   for ( int j=0; j<4; ++j ) {
    table_short[16*i+6-2*j] = table_noshort[16*i+6-2*j] = vdwa_table[4*i+j];
    table_short[16*i+7-2*j] = table_noshort[16*i+7-2*j] = vdwb_table[4*i+j];
    table_short[16*i+8+3-j] = fast_table[4*i+j];
    table_short[16*i+12+3-j] = scor_table[4*i+j];
    table_noshort[16*i+8+3-j] = corr_table[4*i+j];
    table_noshort[16*i+12+3-j] = full_table[4*i+j];
   }
  }
#endif 

  for ( i=0; i<n; ++i ) {
    table_short[16*i+0] = table_noshort[16*i+0] = -6.*vdwa_table[4*i+3];
    table_short[16*i+1] = table_noshort[16*i+1] = -4.*vdwa_table[4*i+2];
    table_short[16*i+2] = table_noshort[16*i+2] = -2.*vdwa_table[4*i+1];
    table_short[16*i+3] = table_noshort[16*i+3] = -1.*vdwa_table[4*i+0];
    
    table_short[16*i+4] = table_noshort[16*i+4] = -6.*vdwb_table[4*i+3];
    table_short[16*i+5] = table_noshort[16*i+5] = -4.*vdwb_table[4*i+2];
    table_short[16*i+6] = table_noshort[16*i+6] = -2.*vdwb_table[4*i+1];
    table_short[16*i+7] = table_noshort[16*i+7] = -1.*vdwb_table[4*i+0];
    
    table_short[16*i+8]  = -6.*fast_table[4*i+3];
    table_short[16*i+9]  = -4.*fast_table[4*i+2];
    table_short[16*i+10] = -2.*fast_table[4*i+1];
    table_short[16*i+11] = -1.*fast_table[4*i+0];

    table_noshort[16*i+8]  = -6.*corr_table[4*i+3];
    table_noshort[16*i+9]  = -4.*corr_table[4*i+2];
    table_noshort[16*i+10] = -2.*corr_table[4*i+1];
    table_noshort[16*i+11] = -1.*corr_table[4*i+0];

    table_short[16*i+12] = -6.*scor_table[4*i+3];
    table_short[16*i+13] = -4.*scor_table[4*i+2];
    table_short[16*i+14] = -2.*scor_table[4*i+1];
    table_short[16*i+15] = -1.*scor_table[4*i+0];

    table_noshort[16*i+12] = -6.*full_table[4*i+3];
    table_noshort[16*i+13] = -4.*full_table[4*i+2];
    table_noshort[16*i+14] = -2.*full_table[4*i+1];
    table_noshort[16*i+15] = -1.*full_table[4*i+0];
  }

#if 0
  char fname[100];
  sprintf(fname,"/tmp/namd.table.pe%d.dat",CkMyPe());
  FILE *f = fopen(fname,"w");
  for ( i=0; i<(n-1); ++i ) {
    const BigReal r2_base = r2_delta * ( 1 << (i/64) );
    const BigReal r2_del = r2_base / 64.0;
    const BigReal r2 = r2_base - r2_delta + r2_del * (i%64);
    BigReal *t;
    if ( r2 + r2_delta != r2_table[i] ) fprintf(f,"r2 error! ");
    fprintf(f,"%g",r2);
    t = fast_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = scor_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = slow_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = corr_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = full_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = vdwa_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = vdwb_table + 4*i;
    fprintf(f,"   %g %g %g %g", t[0], t[1], t[2], t[3]);
    fprintf(f,"\n");
  }
  fclose(f);
#endif

  //Flip slow table to match table_four_i
  for ( i=0; i<n; ++i ) {
    BigReal tmp0, tmp1, tmp2, tmp3;
    tmp0 = slow_table [i*4 + 0];
    tmp1 = slow_table [i*4 + 1];
    tmp2 = slow_table [i*4 + 2];
    tmp3 = slow_table [i*4 + 3];

    slow_table [i*4 + 0] = tmp3;
    slow_table [i*4 + 1] = tmp2;
    slow_table [i*4 + 2] = tmp1;
    slow_table [i*4 + 3] = tmp0;
  }

#ifdef NAMD_CUDA
  if (!simParams->useCUDA2) {
    send_build_cuda_force_table();
  }
#endif

  #ifdef NAMD_MIC
    send_build_mic_force_table();
  #endif
}

