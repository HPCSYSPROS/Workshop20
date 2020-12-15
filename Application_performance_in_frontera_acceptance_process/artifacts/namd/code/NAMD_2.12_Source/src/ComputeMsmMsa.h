/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMSMMSA_H
#define COMPUTEMSMMSA_H

#ifdef CHARM_HAS_MSA

#include <vector>
#include "ComputeHomePatches.h"
#include "NamdTypes.h"
#include "ComputeMoa.h"  // needed for Int3 declaration
#include "msa/msa.h"

typedef MSA::MSA3D<float, DefaultEntry<float>,
        MSA_DEFAULT_ENTRIES_PER_PAGE> MsmMsaGrid;

#if 0
struct Int3 {
  int nx, ny, nz;
  Int3() : nx(0), ny(0), nz(0) { }
  Int3(int mx, int my, int mz) : nx(mx), ny(my), nz(mz) { }
  void pup(PUP::er &p) { p|nx, p|ny, p|nz; }
};
#endif

class SubmitReduction;
typedef Force MsmMsaForce;

class ComputeMsmMsa : public ComputeHomePatches {
public:
  ComputeMsmMsa(ComputeID c);
  virtual ~ComputeMsmMsa();
  void doWork();
  void saveResults(int n, const MsmMsaForce [], double self_energy);

private:
  double qscaling;  // charge scaling constant
  SubmitReduction *reduction;
};

struct MsmMsaData {
  int ispx, ispy, ispz;
  float hx_1, hy_1, hz_1;
  float a;

  float origin_x, origin_y, origin_z;

  int nlevels, maxlevels, toplevel;

  int approx;
  int split;

  std::vector<MsmMsaGrid> qh;
  std::vector<MsmMsaGrid> eh;
  std::vector<Int3> grid_len;  // grid points in each dimension for each level
  std::vector<Int3> grid_idstart;  // starting index for each level

  std::vector<float> scaling;  // scaling factor for each grid level
  Int3 gc_len;            // length of grid cutoff stencil in each dimension
  Int3 gc_idstart;        // starting index of grid cutoff stencil
  std::vector<float> gc;  // grid cutoff stencil

  Int3 gctop_len;
  Int3 gctop_idstart;
  std::vector<float> gctop;  // grid cutoff stencil for top level

  std::vector<int> num_clients_qh;  // number client chares for each qh grid
  std::vector<int> num_clients_eh;  // number client chares for each eh grid

  int num_anterpolation_chares;  // number of chares doing anterpolation
  int num_interpolation_chares;  // number of chares doing interpolation
  std::vector<int> num_restriction_chares;  // number restrictions per level
  std::vector<int> num_prolongation_chares; // number prolongations per level
  std::vector<int> num_gridcutoff_chares;   // number grid-cutoff-s per level
  std::vector<Int3> dim_gridcutoff_chares;  // grid cutoff chare dim per level
  std::vector<Int3> dim_gridtransfer_chares; // grid trans chare dim per level
  int num_total_restriction_chares;
  int num_total_prolongation_chares;
  int num_total_gridcutoff_chares;
  int num_energy_chares;  // number of energy summing chares

  Int3 num_points_per_chare;  // size of grid point sub-cubes

  double self_energy_const;

  void pup(PUP::er &p);  // for parameter marshalling
  void print();          // for debugging
};

#endif // CHARM_HAS_MSA

#endif // COMPUTEMSMMSA_H

