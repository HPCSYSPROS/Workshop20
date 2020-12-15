/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDCUDAEXCL_INL
#define COMPUTENONBONDEDCUDAEXCL_INL

#include "ComputeNonbondedCUDAExcl.h"

inline ExclElem::ExclElem() { ; }

inline ExclElem::ExclElem(AtomID atom0, const TupleSignature *sig, const int *v){
    atomID[0] = atom0;
    atomID[1] = atom0 + sig->offset[0];
    modified = sig->tupleParamType;
}

inline ExclElem::ExclElem(const Exclusion *a, const int *v)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    modified = a->modified;
  }

inline ExclElem::ExclElem(AtomID atom0, AtomID atom1)
  {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }

inline int ExclElem::operator==(const ExclElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

inline int ExclElem::operator<(const ExclElem &a) const
  {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }

#endif

