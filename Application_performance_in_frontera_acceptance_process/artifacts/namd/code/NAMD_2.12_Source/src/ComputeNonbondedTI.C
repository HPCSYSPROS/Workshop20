/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#include "ComputeNonbondedInl.h"

// 3 inline functions to handle explicit calculation of separation-shifted
// vdW for FEP and TI, and a shifted electrostatics potential for decoupling

/* ********************************** */
/* vdW energy, force and dU/dl for TI */
/* ********************************** */
inline void ti_vdw_force_energy_dUdl (BigReal A, BigReal B, BigReal r2, 
  BigReal myVdwShift, BigReal switchdist2, BigReal cutoff2, 
  BigReal myVdwLambda, BigReal alchVdwShiftCoeff, BigReal switchfactor, 
  Bool vdwForceSwitching, BigReal* alch_vdw_energy, BigReal* alch_vdw_force,
  BigReal* alch_vdw_dUdl) {
  //myVdwShift already multplied by relevant (1-vdwLambda)
  const BigReal r2_1 = 1./(r2 + myVdwShift);
  const BigReal r6_1 = r2_1*r2_1*r2_1;
    
  // switching function (this is correct whether switching is active or not)
  const BigReal switchmul = (r2 > switchdist2 ? switchfactor*(cutoff2 - r2) \
			     *(cutoff2 - r2) \
			     *(cutoff2 - 3.*switchdist2 + 2.*r2) : 1.);

  // separation-shifted vdW force and energy
  const BigReal U = A*r6_1*r6_1 - B*r6_1; // NB: unscaled! for shorthand only!
  *alch_vdw_energy = myVdwLambda*U;
  *alch_vdw_force = myVdwLambda*switchmul*(12.*U + 6.*B*r6_1)*r2_1;
  *alch_vdw_dUdl = U + myVdwLambda*alchVdwShiftCoeff*(6.*U + 3.*B*r6_1)*r2_1;

  // BKR - separation-shifted vdW force switching and potential shifting
  if(!vdwForceSwitching){ // add on chain rule for switch function
    const BigReal switchmul2 = (r2 > switchdist2 ?                      \
                                12.*switchfactor*(cutoff2 - r2)         \
                                *(r2 - switchdist2) : 0.);
    *alch_vdw_energy *= switchmul;
    *alch_vdw_force += myVdwLambda*switchmul2*U;
    *alch_vdw_dUdl *= switchmul;
  }else{ // add potential shifts and additional dU/dl terms
    const BigReal switchdist6 = switchdist2*switchdist2*switchdist2;
    const BigReal cutoff6 = cutoff2*cutoff2*cutoff2;
    const BigReal switchdist3 = switchdist2*sqrt(switchdist2);
    const BigReal cutoff3 = cutoff2*sqrt(cutoff2);
    if(r2 > switchdist2) {
      const BigReal k_vdwa = A*cutoff6 / (cutoff6 - switchdist6);
      const BigReal k_vdwb = -B*cutoff3 / (cutoff3 - switchdist3);
      const BigReal tmpa = r6_1 - (1./cutoff6);
      const BigReal r_1 = sqrt(r2_1);
      const BigReal tmpb = r2_1*r_1 - (1./cutoff3);
      const BigReal Uh = k_vdwa*tmpa*tmpa + k_vdwb*tmpb*tmpb; //"harmonic" U
      *alch_vdw_energy = myVdwLambda*Uh;
      *alch_vdw_dUdl = (Uh + alchVdwShiftCoeff*myVdwLambda      \
                        *(6.*k_vdwa*tmpa*r6_1*r2_1              \
                          + 3.*k_vdwb*tmpb*r2_1*r2_1*r_1));
    }else{
      const BigReal v_vdwa = -A / (switchdist6*cutoff6);
      const BigReal v_vdwb = B / (switchdist3*cutoff3);
      const BigReal dU = v_vdwa + v_vdwb; //deltaV2 from Steinbach & Brooks
      *alch_vdw_energy += myVdwLambda*dU;
      *alch_vdw_dUdl += dU;
    }
  }
}


/*************THERMODYNAMIC INTEGRATION*************/
#define TIFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef TIFLAG

