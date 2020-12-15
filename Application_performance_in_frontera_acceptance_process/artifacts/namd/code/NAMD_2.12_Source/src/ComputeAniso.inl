/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEANISO_INL
#define COMPUTEANISO_INL

#include "ComputeAniso.h"

inline AnisoElem::AnisoElem() { ; }

inline AnisoElem::AnisoElem(AtomID atom0, const TupleSignature *sig, const AnisoValue *v){
    NAMD_die("Can't use Aniso with memory optimized version of NAMD.");
    // atomID[0] = atom0;
    // atomID[1] = atom0 + sig->offset[0];
    // atomID[2] = atom0 + sig->offset[1];
    // atomID[3] = atom0 + sig->offset[2];
    // value = &v[sig->tupleParamType];
}

inline AnisoElem::AnisoElem(const Aniso *a, const AnisoValue *v)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    value = a;  // expect v to be NULL
  }

inline AnisoElem::AnisoElem(AtomID atom0, AtomID atom1,
				  AtomID atom2, AtomID atom3)
  {
    // do not rearrange atom ordering of Aniso
    // the first atom is special
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
    atomID[3] = atom3;
  }

inline int AnisoElem::operator==(const AnisoElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2] && a.atomID[3] == atomID[3]);
  }

inline int AnisoElem::operator<(const AnisoElem &a) const
  {
    return  (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
            (atomID[2] < a.atomID[2] ||
            (atomID[2] == a.atomID[2] &&
             atomID[3] < a.atomID[3] 
	     ))))));
  }

#endif

