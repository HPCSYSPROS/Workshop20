/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDSPLITTING_H
#define COMPUTENONBONDEDSPLITTING_H

// Several special cases are defined:
//   NBPAIR, NBSELF, NBEXCL switch environment (mutually exclusive)
//   FULLELECT full electrostatics calculation?

#include "LJTable.h"
#include "Molecule.h"
#include "ComputeNonbondedUtil.h"

// ************************************************************
// Various switching functions (for inlining)
inline void ComputeNonbondedUtil::shifting (
	BigReal &shiftVal, BigReal &dShiftVal,
	const BigReal &r, const BigReal &r2,
	const BigReal &c5, const BigReal &c6)
{
  // Basic electrostatics shifting function for cutoff simulations
  shiftVal = 1.0 - r2*c5;
  dShiftVal = c6*shiftVal*r;
  shiftVal *= shiftVal;
}

inline void ComputeNonbondedUtil::xplorsplitting (
	BigReal &shiftVal, BigReal &dShiftVal,
	const BigReal &switchVal, const BigReal &dSwitchVal)
{
  // X-plor electrostatics splitting function for multiple timestepping
  // Same as X-plor VdW switching function so copy from above.
  shiftVal = switchVal;
  dShiftVal = dSwitchVal;
}

inline void ComputeNonbondedUtil::c1splitting (
	BigReal &shiftVal, BigReal &dShiftVal,
	const BigReal &r, const BigReal &d0, const BigReal &switchOn)
{
  // C1 electrostatics splitting function for multiple timestepping
  dShiftVal = 0;  // formula only correct for forces
  if (r > switchOn)
    {
    const BigReal d1 = d0*(r-switchOn);
    shiftVal = 1.0 + d1*d1*(2.0*d1-3.0);
    }
  else
    {
    shiftVal = 1;
    }
}

#endif

