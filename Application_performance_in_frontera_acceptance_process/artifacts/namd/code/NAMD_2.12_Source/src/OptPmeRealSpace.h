/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef OPT_PME_REAL_SPACE_H__
#define OPT_PME_REAL_SPACE_H__

#include "PmeBase.h"
#include "Vector.h"

class OptPmeRealSpace {
  
public:
  OptPmeRealSpace(PmeGrid grid, int natoms);
  ~OptPmeRealSpace();
  
  void fill_charges(double **q_arr, PmeParticle p[], int zstart, int zlen); 
  void compute_forces(const double * const *q_arr, const PmeParticle p[], Vector f[], int zstart, int zlen, int start=0, int end=0);

  const int N;
  const PmeGrid myGrid;
  double *M, *dM;
};


#endif

