/* *
 * ***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
 * ***  The Board of Trustees of the University of Illinois.
 * ***  All rights reserved.
 * **
 * =====================================================================================
 *
 *       Filename:  ComputeMoa.h
 *
 *    Description:  
 *
 *        Version:  Stub File
 *        Created:  03/08/2012 03:55:27 PM
 *       Revision:  
 *       Compiler:  charm++ 
 *
 *         Author:  Christopher B Harrison, Ph.D.
 *         Email:   charris5@gmail.com, char1@illinois.edu, char@ks.uiuc.edu
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef COMPUTEOA_H
#define COMPUTEOA_H

#include <vector>
#include "ComputeHomePatches.h"
#include "NamdTypes.h"

#ifdef CHARM_HAS_MSA

#include "msa/msa.h"

typedef MSA::MSA1D<double, DefaultEntry<double>,
        MSA_DEFAULT_ENTRIES_PER_PAGE> Moa1Grid;

typedef MSA::MSA3D<double, DefaultEntry<double>,
        MSA_DEFAULT_ENTRIES_PER_PAGE> Moa3Grid;

#else

typedef int Moa1Grid;
typedef int Moa3Grid;

#endif // CHARM_HAS_MSA

struct Int2 {
  int nx, ny;
  Int2() : nx(0), ny(0) { }
  Int2(int mx, int my) : nx(mx), ny(my) { }
  void pup(PUP::er &p) { p|nx, p|ny; }
};

// /* Experimental */
// struct Vec3 {
//   std::vector<int> b1;  // b1 array 
//   std::vector<int> b2;  // b1 array 
//   std::vector<int> b3;  // b1 array 
//   void pup(PUP::er &p) { p|b1, p|b2, p|b3; }
// }


struct Int3 {
  int nx, ny, nz;
  Int3() : nx(0), ny(0), nz(0) { }
  Int3(int mx, int my, int mz) : nx(mx), ny(my), nz(mz) { }
  void pup(PUP::er &p) { p|nx, p|ny, p|nz; }
};

class SubmitReduction;

class ComputeMoa : public ComputeHomePatches {
public:
  ComputeMoa(ComputeID c);
  virtual ~ComputeMoa();
  void doWork();

private:
  SubmitReduction *reduction;
};

struct MoaData {

  int K1, K2, K3;
  int order;
  float orig_x, orig_y, orig_z;

  int num_clients;
  std::vector<int> num_clients_qh;  // number client chares for each qh grid
  std::vector<int> num_clients_bh;  // number client chares for each bh grid
  std::vector<int> num_clients_sh;  // number client chares for each bh grid


  std::vector<Moa3Grid> qh;      // charge grid
  std::vector<Moa3Grid> sh;      // s-value grid
//  std::vector<Moa1Grid> b1;      // b-value grid
//  std::vector<Moa1Grid> b2;      // b-value grid
//  std::vector<Moa1Grid> b3;      // b-value grid
  std::vector<Moa3Grid> bh;      // b-value grid

  std::vector<Int2> k1r, k2r, k3r;
  Int3 gbqs;
  std::vector<int> hasLQ;
  std::vector<int> hasLB;
  std::vector<int> hasLS;



  void pup(PUP::er &p);  // for parameter marshalling
  void print();          // for debugging
};

#endif

