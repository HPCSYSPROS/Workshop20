/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifndef COMPUTENONBONDEDUTIL_H
#define COMPUTENONBONDEDUTIL_H

#include "NamdTypes.h"
#include "ReductionMgr.h"
#include "Molecule.h"

class LJTable;
class Molecule;
class Random;

#ifdef NAMD_CUDA
void register_cuda_compute_self(ComputeID c, PatchID pid);
void register_cuda_compute_pair(ComputeID c, PatchID pid[], int t[]);
void unregister_cuda_compute(ComputeID c);
#endif
#ifdef NAMD_MIC
void register_mic_compute_self(ComputeID c, PatchID pid, int part, int numParts);
void register_mic_compute_pair(ComputeID c, PatchID pid[], int t[], int part, int numParts);
void unregister_mic_compute(ComputeID c);
#endif

typedef unsigned short plint;

class Pairlists {
  enum {initsize = 10};
  plint *data;
  int curpos;
  int size;
  Pairlists(const Pairlists&) { ; }
  Pairlists& operator=(const Pairlists&) { return *this; }
public:
  Pairlists() : size(initsize) { data = new plint[initsize]; }
  ~Pairlists() { delete [] data; }
  plint *newlist(int max_size) {  // get a new list w/ room for max_size
    int reqnewsize = curpos + max_size + 1;
    int newsize = size;
    while ( newsize < reqnewsize ) { newsize += newsize >> 1; }
    if ( newsize > size ) {
      plint *newdata = new plint[newsize];
      CmiMemcpy(newdata,data,curpos*sizeof(plint));
      delete [] data;
      data = newdata;
      size = newsize;
    }
    return &data[curpos+1];
  }

  // don't specify size if previous allocation should have extra space
  plint *newlist() {  // get a new list assuming already allocated
    return &data[curpos+1];
  }

  void newsize(int list_size) {  // set the size of the last list gotten
    data[curpos] = list_size;
    curpos += list_size + 1;
  }
  void reset() { curpos = 0; }  // go back to the beginning
  void nextlist(plint **list, int *list_size) {  // get next list and size
    *list = &data[curpos+1];
    curpos += ( *list_size = data[curpos] ) + 1;
  }
  int getSize() { return size; }

  void addIndex() {  // assume space for index already allocated
    curpos++;
  }
  void setIndexValue(plint i) {  // assume no newsize since addIndex
    data[curpos-1] = i;
  }

  plint getIndexValue() {
    return data[curpos++];
  }

};

#define NBWORKARRAYSINIT(ARRAYS) \
  ComputeNonbondedWorkArrays* const computeNonbondedWorkArrays = ARRAYS;

#ifdef __INTEL_COMPILER
#define NBWORKARRAY(TYPE,NAME,SIZE) \
  computeNonbondedWorkArrays->NAME.resize(SIZE); \
  TYPE * const NAME = computeNonbondedWorkArrays->NAME.begin(); \
  __assume_aligned(NAME,64);
#else
#define NBWORKARRAY(TYPE,NAME,SIZE) \
  computeNonbondedWorkArrays->NAME.resize(SIZE); \
  TYPE * const NAME = computeNonbondedWorkArrays->NAME.begin();
#endif

class ComputeNonbondedWorkArrays {
public:
  ResizeArray<int> pairlisti;
  ResizeArray<BigReal> r2list;
#ifdef NAMD_KNL
  ResizeArray<float> r2list_f;
  ResizeArray<float> xlist;
  ResizeArray<float> ylist;
  ResizeArray<float> zlist;
#endif

  // DMK - Atom Sort
  // NOTE : TODO : For pair nonbonded compute objects, these arrays are needed
  //   for the atom sorting code.  However, grouplist and fixglist will not be
  //   used if these arrays are so they could overlap in memory.  TODO: Find a
  //   way to allow them to use the same memory.  For now, it's not much memory
  //   so just keep them separate because of the type differences.  I haven't
  //   looked through all the details yet, but could probably just use the
  //   larger type and cast the pointer type as needed for the other one to
  //   ensure there is enough memory in both cases.
  #if NAMD_ComputeNonbonded_SortAtoms != 0
    ResizeArray<SortEntry> atomSort_0_sortValues__;
    ResizeArray<SortEntry> atomSort_1_sortValues__;
    ResizeArray<BigReal> p_0_sortValues;
  #endif

  ResizeArray<int> grouplist;
  ResizeArray<int> fixglist;
  ResizeArray<int> goodglist;
  ResizeArray<plint> pairlistx;
  ResizeArray<plint> pairlistm;

  ResizeArray<plint> pairlistnAlch;
  ResizeArray<plint> pairlistnA0;
  ResizeArray<plint> pairlistnA1;
  ResizeArray<plint> pairlistxA1;
  ResizeArray<plint> pairlistmA1;
  ResizeArray<plint> pairlistnA2;
  ResizeArray<plint> pairlistxA2;
  ResizeArray<plint> pairlistmA2;
  
  ResizeArray<int> pairlist;
  ResizeArray<int> pairlist2;
  ResizeArray<Force> f_0;
  ResizeArray<Force> fullf_0;
};

//struct sent to CalcGBIS
struct GBISParamStruct {
  int cid;
  int patchID[2];
  int sequence;
  int doSmoothing;
  int doGBIS;
  int numPatches;
  int gbisPhase;
  Real *intRad[2];
  GBReal *psiSum[2];
  Real *bornRad[2];
  GBReal *dEdaSum[2];
  Real *dHdrPrefix[2];
  BigReal epsilon_s;
  BigReal epsilon_p;
  BigReal rho_0;
  BigReal kappa;
  BigReal cutoff;
  BigReal a_cut;
  BigReal fsMax;
  BigReal delta;
  BigReal beta;
  BigReal gamma;
  BigReal alpha_max;
  BigReal gbSelfEnergy;
  BigReal gbInterEnergy;
  int doFullElectrostatics;
  int doEnergy;
  Pairlists *gbisStepPairlists[4];
  BigReal maxGroupRadius;
};

// function arguments
struct nonbonded {
  CompAtom* p[2];
#ifdef NAMD_KNL
  CompAtomFlt *pFlt[2];
#endif
  CompAtomExt *pExt[2];
  // BEGIN LA
  CompAtom* v[2];
  // END LA
  Force* ff[2];
  // for full electrostatics
  Force* fullf [2];

  int numAtoms[2];

  Vector offset;
  Vector offset_f;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int numWaterAtoms[2];
  #endif

  BigReal *reduction;
  BigReal *pressureProfileReduction;

  Parameters *parameters;
  SimParameters *simParameters;
  Random *random;

  ComputeNonbondedWorkArrays *workArrays;

  Pairlists *pairlists;
  int savePairlists;
  int usePairlists;
  BigReal plcutoff;
  BigReal groupplcutoff;

  // BEGIN LA
  int doLoweAndersen;
  // END LA
  
  int minPart;
  int maxPart;
  int numParts;
  int step;

  #if NAMD_ComputeNonbonded_SortAtoms != 0
    Vector projLineVec;
  #endif
};

class ComputeNonbondedUtil {

public:

  ComputeNonbondedUtil() {}
  ~ComputeNonbondedUtil() {}
  static void select(void);

  static void (*calcPair)(nonbonded *);
  static void (*calcPairEnergy)(nonbonded *);
  static void (*calcSelf)(nonbonded *);
  static void (*calcSelfEnergy)(nonbonded *);

  static void (*calcFullPair)(nonbonded *);
  static void (*calcFullPairEnergy)(nonbonded *);
  static void (*calcFullSelf)(nonbonded *);
  static void (*calcFullSelfEnergy)(nonbonded *);

  static void (*calcMergePair)(nonbonded *);
  static void (*calcMergePairEnergy)(nonbonded *);
  static void (*calcMergeSelf)(nonbonded *);
  static void (*calcMergeSelfEnergy)(nonbonded *);

  static void (*calcSlowPair)(nonbonded *);
  static void (*calcSlowPairEnergy)(nonbonded *);
  static void (*calcSlowSelf)(nonbonded *);
  static void (*calcSlowSelfEnergy)(nonbonded *);

  // Below modified by JLai -- Look at bottom 3rd line
  enum { exclChecksumIndex, pairlistWarningIndex,
	 electEnergyIndex, fullElectEnergyIndex, vdwEnergyIndex,
	 goNativeEnergyIndex, goNonnativeEnergyIndex, groLJEnergyIndex, groGaussEnergyIndex,
//sd-db
	 electEnergyIndex_s, fullElectEnergyIndex_s, vdwEnergyIndex_s, vdwEnergyIndex_s_Left,
	 electEnergyIndex_ti_1, fullElectEnergyIndex_ti_1, vdwEnergyIndex_ti_1,
	 electEnergyIndex_ti_2, fullElectEnergyIndex_ti_2, vdwEnergyIndex_ti_2,
//sd-de
	 TENSOR(virialIndex), TENSOR(fullElectVirialIndex),
         VECTOR(pairVDWForceIndex), VECTOR(pairElectForceIndex),
	 reductionDataSize };
  static void submitReductionData(BigReal*,SubmitReduction*);
  static void submitPressureProfileData(BigReal*,SubmitReduction*);

  static Bool commOnly;
  static Bool fixedAtomsOn;
  static Bool qmForcesOn ;
  static BigReal cutoff;
  static BigReal cutoff2;
  static float cutoff2_f;
  static BigReal dielectric_1;
  static const LJTable* ljTable;
  static const Molecule* mol;
  static BigReal r2_delta, r2_delta_1;
  static int rowsize;
  static int columnsize;
  static int r2_delta_exp;
  static BigReal *table_alloc;
  static BigReal *table_ener;
  static BigReal *table_short;
  static BigReal *table_noshort;
  static BigReal *fast_table;
  static BigReal *scor_table;
  static BigReal *slow_table;
  static BigReal *corr_table;
  static BigReal *full_table;
  static BigReal *vdwa_table;
  static BigReal *vdwb_table;
  static BigReal *r2_table;
  #if defined(NAMD_MIC)
    static BigReal *mic_table_base_ptr; // DMK - NOTE : Duplicate but, use so that nothing breaks if the ordering of sub-arrays changes
    static int mic_table_n;
    static int mic_table_n_16;
  #endif
  #ifdef NAMD_KNL
  #define KNL_TABLE_SIZE (4080+2)
  static float *knl_table_alloc;
  static float *knl_fast_ener_table;
  static float *knl_fast_grad_table;
  static float *knl_scor_ener_table;
  static float *knl_scor_grad_table;
  static float *knl_slow_ener_table;
  static float *knl_slow_grad_table;
  static float *knl_corr_ener_table;
  static float *knl_corr_grad_table;
  static float *knl_full_ener_table;
  static float *knl_full_grad_table;
  #endif
  static BigReal scaling;
  static BigReal scale14;
  static BigReal switchOn;
  static BigReal switchOn_1;
  static BigReal switchOn2;
  static BigReal v_vdwa;
  static BigReal v_vdwb;
  static BigReal k_vdwa;
  static BigReal k_vdwb;
  static BigReal cutoff_3;
  static BigReal cutoff_6;
  static float   v_vdwa_f;
  static float   v_vdwb_f;
  static float   k_vdwa_f;
  static float   k_vdwb_f;
  static float   cutoff_3_f;
  static float   cutoff_6_f;
  static float   switchOn_f;
  static float   A6_f;
  static float   B6_f;
  static float   C6_f;
  static float   A12_f;
  static float   B12_f;
  static float   C12_f;
  static BigReal c0;
  static BigReal c1;
  static BigReal c3;
  static BigReal c5;
  static BigReal c6;
  static BigReal c7;
  static BigReal c8;
  // static BigReal d0;
//sd-db
  static Bool alchFepOn;
  static Bool alchThermIntOn;
  static BigReal alchLambda2;
  static BigReal alchVdwShiftCoeff;
  static Bool vdwForceSwitching;
  static Bool Fep_WCA_repuOn;
  static Bool Fep_WCA_dispOn;
  static Bool Fep_ElecOn;
  static Bool Fep_Wham;
  static BigReal WCA_rcut1;
  static BigReal WCA_rcut2; 
  static BigReal WCA_rcut3; 
  static Bool alchDecouple;
  static BigReal alchRepLambda;
  static BigReal alchDispLambda;
  static BigReal alchElecLambda;
//sd-de
  static Bool lesOn;
  static int lesFactor;
  static BigReal lesScaling;

  static BigReal *lambda_table;

  static Bool pairInteractionOn;
  static Bool pairInteractionSelf;

  static Bool pressureProfileOn;
  static int pressureProfileSlabs;
  static int pressureProfileAtomTypes;
  static BigReal pressureProfileThickness;
  static BigReal pressureProfileMin;

  static Bool accelMDOn;

  static Bool drudeNbthole;

  // for particle mesh Ewald
  static BigReal ewaldcof;
  static BigReal pi_ewaldcof;

  // need macros for preprocessor
  #define VDW_SWITCH_MODE_ENERGY  0
  #define VDW_SWITCH_MODE_MARTINI 1
  #define VDW_SWITCH_MODE_FORCE   2
  static int vdw_switch_mode;

  // Ported by JLai -- JE - Go
  static Bool goGroPair;
  static Bool goForcesOn;
  static int goMethod; //6.3.11
  // End of port -- JL

  // for simplifying some common functions
  inline static BigReal square(	const BigReal &x,
 				const BigReal &y,
				const BigReal &z)
	{
	return(x*x+y*y+z*z);
	}

  static void calc_error(nonbonded *);

  static void calc_pair(nonbonded *);
  static void calc_pair_energy(nonbonded *);
  static void calc_pair_fullelect(nonbonded *);
  static void calc_pair_energy_fullelect(nonbonded *);
  static void calc_pair_merge_fullelect(nonbonded *);
  static void calc_pair_energy_merge_fullelect(nonbonded *);
  static void calc_pair_slow_fullelect(nonbonded *);
  static void calc_pair_energy_slow_fullelect(nonbonded *);

  static void calc_self(nonbonded *);
  static void calc_self_energy(nonbonded *);
  static void calc_self_fullelect(nonbonded *);
  static void calc_self_energy_fullelect(nonbonded *);
  static void calc_self_merge_fullelect(nonbonded *);
  static void calc_self_energy_merge_fullelect(nonbonded *);
  static void calc_self_slow_fullelect(nonbonded *);
  static void calc_self_energy_slow_fullelect(nonbonded *);

//alchemical fep calcualtion
  static void calc_pair_energy_fep(nonbonded *);
  static void calc_pair_energy_fullelect_fep (nonbonded *);
  static void calc_pair_energy_merge_fullelect_fep (nonbonded *);
  static void calc_pair_energy_slow_fullelect_fep (nonbonded *);
  static void calc_self_energy_fep (nonbonded *);
  static void calc_self_energy_fullelect_fep (nonbonded *);
  static void calc_self_energy_merge_fullelect_fep (nonbonded *);
  static void calc_self_energy_slow_fullelect_fep (nonbonded *);

  static void calc_pair_energy_ti(nonbonded *);
  static void calc_pair_ti(nonbonded *);
  static void calc_pair_energy_fullelect_ti (nonbonded *);
  static void calc_pair_fullelect_ti (nonbonded *);
  static void calc_pair_energy_merge_fullelect_ti (nonbonded *);
  static void calc_pair_merge_fullelect_ti (nonbonded *);
  static void calc_pair_energy_slow_fullelect_ti (nonbonded *);
  static void calc_pair_slow_fullelect_ti (nonbonded *);
  static void calc_self_energy_ti (nonbonded *);
  static void calc_self_ti (nonbonded *);
  static void calc_self_energy_fullelect_ti (nonbonded *);
  static void calc_self_fullelect_ti (nonbonded *);
  static void calc_self_energy_merge_fullelect_ti (nonbonded *);
  static void calc_self_merge_fullelect_ti (nonbonded *);
  static void calc_self_energy_slow_fullelect_ti (nonbonded *);
  static void calc_self_slow_fullelect_ti (nonbonded *);
  
//locally enhanced sampling calcualtion
  static void calc_pair_les(nonbonded *);
  static void calc_pair_energy_les(nonbonded *);
  static void calc_pair_fullelect_les (nonbonded *);
  static void calc_pair_energy_fullelect_les (nonbonded *);
  static void calc_pair_merge_fullelect_les (nonbonded *);
  static void calc_pair_energy_merge_fullelect_les (nonbonded *);
  static void calc_pair_slow_fullelect_les (nonbonded *);
  static void calc_pair_energy_slow_fullelect_les (nonbonded *);
  static void calc_self_les (nonbonded *);
  static void calc_self_energy_les (nonbonded *);
  static void calc_self_fullelect_les (nonbonded *);
  static void calc_self_energy_fullelect_les (nonbonded *);
  static void calc_self_merge_fullelect_les (nonbonded *);
  static void calc_self_energy_merge_fullelect_les (nonbonded *);
  static void calc_self_slow_fullelect_les (nonbonded *);
  static void calc_self_energy_slow_fullelect_les (nonbonded *);

//pair interaction calcualtion
  static void calc_pair_energy_int(nonbonded *);
  static void calc_pair_energy_fullelect_int(nonbonded *);
  static void calc_pair_energy_merge_fullelect_int(nonbonded *);
  static void calc_self_energy_int (nonbonded *);
  static void calc_self_energy_fullelect_int(nonbonded *);
  static void calc_self_energy_merge_fullelect_int(nonbonded *);

//pressure profile calcualtion
  static void calc_pair_pprof(nonbonded *);
  static void calc_pair_energy_pprof(nonbonded *);
  static void calc_pair_fullelect_pprof(nonbonded *);
  static void calc_pair_energy_fullelect_pprof(nonbonded *);
  static void calc_pair_merge_fullelect_pprof(nonbonded *);
  static void calc_pair_energy_merge_fullelect_pprof(nonbonded *);
  static void calc_pair_slow_fullelect_pprof(nonbonded *);
  static void calc_pair_energy_slow_fullelect_pprof(nonbonded *);
  static void calc_self_pprof(nonbonded *);
  static void calc_self_energy_pprof(nonbonded *);
  static void calc_self_fullelect_pprof (nonbonded *);
  static void calc_self_energy_fullelect_pprof (nonbonded *);
  static void calc_self_merge_fullelect_pprof (nonbonded *);
  static void calc_self_energy_merge_fullelect_pprof (nonbonded *);
  static void calc_self_slow_fullelect_pprof (nonbonded *);
  static void calc_self_energy_slow_fullelect_pprof (nonbonded *);
 
  //tabulated potential nonbonded calcs
  static void calc_pair_tabener(nonbonded *);
  static void calc_pair_energy_tabener(nonbonded *);
  static void calc_pair_fullelect_tabener(nonbonded *);
  static void calc_pair_energy_fullelect_tabener(nonbonded *);
  static void calc_pair_merge_fullelect_tabener(nonbonded *);
  static void calc_pair_energy_merge_fullelect_tabener(nonbonded *);
  static void calc_pair_slow_fullelect_tabener(nonbonded *);
  static void calc_pair_energy_slow_fullelect_tabener(nonbonded *);
  static void calc_self_tabener(nonbonded *);
  static void calc_self_energy_tabener(nonbonded *);
  static void calc_self_fullelect_tabener (nonbonded *);
  static void calc_self_energy_fullelect_tabener (nonbonded *);
  static void calc_self_merge_fullelect_tabener (nonbonded *);
  static void calc_self_energy_merge_fullelect_tabener (nonbonded *);
  static void calc_self_slow_fullelect_tabener (nonbonded *);
  static void calc_self_energy_slow_fullelect_tabener (nonbonded *);

  //Go forces nonbonded calcs
  static void calc_pair_go(nonbonded *);
  static void calc_pair_energy_go(nonbonded *);
  static void calc_pair_fullelect_go(nonbonded *);
  static void calc_pair_energy_fullelect_go(nonbonded *);
  static void calc_pair_merge_fullelect_go(nonbonded *);
  static void calc_pair_energy_merge_fullelect_go(nonbonded *);
  static void calc_pair_slow_fullelect_go(nonbonded *);
  static void calc_pair_energy_slow_fullelect_go(nonbonded *);
  static void calc_self_go(nonbonded *);
  static void calc_self_energy_go(nonbonded *);
  static void calc_self_fullelect_go(nonbonded *);
  static void calc_self_energy_fullelect_go(nonbonded *);
  static void calc_self_merge_fullelect_go(nonbonded *);
  static void calc_self_energy_merge_fullelect_go(nonbonded *);
  static void calc_self_slow_fullelect_go(nonbonded *);
  static void calc_self_energy_slow_fullelect_go(nonbonded *);

  void calcGBIS(nonbonded *params, GBISParamStruct *gbisParams);
};

#endif

