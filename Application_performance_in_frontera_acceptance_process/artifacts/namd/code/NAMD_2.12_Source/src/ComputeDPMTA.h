/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEDPMTA_H
#define COMPUTEDPMTA_H

#include "ComputeHomePatches.h"

class SubmitReduction;

#ifdef DPMTA

extern "C"
  {
  #include "dpmta.h"
  }

class ComputeDPMTA : public ComputeHomePatches {
private:
  int *slavetids;	//  PID for slave processes
  int totalAtoms;	//  Total number of atoms being dealt with
  PmtaPartInfo *fmaResults;	//  Results from the PMTA code
  PmtaPartInfo *ljResults;	//  Results from the PMTA code
  Vector boxsize;	// FMA box size, set by get_FMA_cube()
  Vector boxcenter;	// FMA box center, set by get_FMA_cube()
  int usePBC;		// flag for PBC
  Vector initLattice;	// initial system lattice dimensions
  SubmitReduction *reduction;

  void get_FMA_cube(int resize);

public:
  ComputeDPMTA(ComputeID c);
  void initialize();
  virtual ~ComputeDPMTA();
  void doWork();
};

#endif
#endif

