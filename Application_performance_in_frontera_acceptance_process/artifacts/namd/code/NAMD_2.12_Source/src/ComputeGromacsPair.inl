/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEGROMACSPAIR_INL
#define COMPUTEGROMACSPAIR_INL

#include "ComputeGromacsPair.h"

inline GromacsPairElem::GromacsPairElem() { ; }

inline GromacsPairElem::GromacsPairElem(AtomID atom0, const TupleSignature *sig, const GromacsPairValue *v){
    atomID[0] = atom0;
    atomID[1] = atom0 + sig->offset[0];
    value = &v[sig->tupleParamType];
}

inline GromacsPairElem::GromacsPairElem(const GromacsPair *a, const GromacsPairValue *v)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    value = &v[a->gromacsPair_type];
  }

inline GromacsPairElem::GromacsPairElem(AtomID atom0, AtomID atom1)
  {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }

inline int GromacsPairElem::operator==(const GromacsPairElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

inline int GromacsPairElem::operator<(const GromacsPairElem &a) const
  {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }

#endif

