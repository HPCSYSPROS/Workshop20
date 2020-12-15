/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_K_SPACE_H__
#define PME_K_SPACE_H__

#include "PmeBase.h"
#include "Lattice.h"

#include "ComputeMoa.h"
#include "ComputeMoaMgr.decl.h" 

class PmeKSpace {
#ifdef OPENATOM_VERSION
  friend class ComputeMoaMgr;
#endif //OPENATOM_VERSION

public:
  PmeKSpace(PmeGrid grid, int K2_start, int K2_end, int K3_start, int K3_end);
#ifdef OPENATOM_VERSION
  PmeKSpace(PmeGrid grid, int K2_start, int K2_end, int K3_start, int K3_end, CProxy_ComputeMoaMgr moaProxy);
#endif //OPENATOM_VERSION
  ~PmeKSpace();

  double compute_energy(float q_arr[], const Lattice &lattice, double ewald, double virial[], int useCkLoop);
  double compute_energy_orthogonal_helper(float q_arr[], const Lattice &lattice, double ewald, double virial[]);
  void compute_energy_orthogonal_subset(float q_arr[], double *recips, double partialVirial[], double *partialEnergy, int k1from, int k1to);
  

private:
  // b-spline moduli
  double *bm1, *bm2, *bm3; 
  double *exp1, *exp2, *exp3;
  double i_pi_volume, piob;

  const PmeGrid myGrid;
  const int k2_start, k2_end, k3_start, k3_end;

  void init_exp(double *xp, int K, int k_start, int k_end, double recip);
};

#endif
