/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEBOND_INL
#define COMPUTEBOND_INL

#include "ComputeBonds.h"

inline BondElem::BondElem() { ; }

inline BondElem::BondElem(AtomID atom0, const TupleSignature *sig, const BondValue *v){
    atomID[0] = atom0;
    atomID[1] = atom0 + sig->offset[0];
    value = &v[sig->tupleParamType];
}

inline BondElem::BondElem(const Bond *a, const BondValue *v)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    value = &v[a->bond_type];
  }

inline BondElem::BondElem(AtomID atom0, AtomID atom1)
  {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }

inline int BondElem::operator==(const BondElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

inline int BondElem::operator<(const BondElem &a) const
  {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }

#endif

