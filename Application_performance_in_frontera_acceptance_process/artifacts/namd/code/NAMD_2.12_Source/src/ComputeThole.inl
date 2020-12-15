/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTETHOLE_INL
#define COMPUTETHOLE_INL

#include "ComputeThole.h"

inline TholeElem::TholeElem() { ; }

inline TholeElem::TholeElem(AtomID atom0, const TupleSignature *sig, const TholeValue *v){
    NAMD_die("Can't use Thole with memory optimized version of NAMD.");
    // atomID[0] = atom0;
    // atomID[1] = atom0 + sig->offset[0];
    // atomID[2] = atom0 + sig->offset[1];
    // atomID[3] = atom0 + sig->offset[2];
    // value = &v[sig->tupleParamType];
}

inline TholeElem::TholeElem(const Thole *a, const TholeValue *v)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    value = a;  // expect v to be NULL
  }

inline TholeElem::TholeElem(AtomID atom0, AtomID atom1,
				  AtomID atom2, AtomID atom3)
  {
    // atoms arranged:  HEAVY DRUDE HEAVY DRUDE
    if (atom0 > atom2) {  // Swap heavy atoms so lowest is first!
      AtomID tmp = atom2; atom2 = atom0; atom0 = tmp; 
      tmp = atom1; atom1 = atom3; atom3 = tmp;
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
    atomID[3] = atom3;
  }

inline int TholeElem::operator==(const TholeElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2] && a.atomID[3] == atomID[3]);
  }

inline int TholeElem::operator<(const TholeElem &a) const
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

