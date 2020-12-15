/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "common.h"

/*
 * 12/27/2005: switch to Harasima contour, in which the interaction is split 
 * evenly between the pairs of atoms, rather than distributed among the slabs
 * between them.
 */

inline void pp_clamp(int &n, int nslabs) {
  int a = n < 0 ? nslabs : 0;
  int b = n >= nslabs ? nslabs : 0;
  n += a-b;
}


inline void pp_reduction(int nslabs, int n1, int n2, 
                        int atype1, int atype2, int numtypes,
                        BigReal vxx, BigReal vyy, BigReal vzz,
                        BigReal *reduction) {

  int slaboffset = atype1*numtypes + atype2;
  reduction += slaboffset * 3*nslabs;

  vxx *= 0.5;
  vyy *= 0.5;
  vzz *= 0.5;
  reduction[3*n1  ] += vxx;
  reduction[3*n1+1] += vyy;
  reduction[3*n1+2] += vzz;
  reduction[3*n2  ] += vxx;
  reduction[3*n2+1] += vyy;
  reduction[3*n2+2] += vzz;
}

