/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDPAIR_H
#define COMPUTENONBONDEDPAIR_H

#include "ComputePatchPair.h"
#include "ComputeNonbondedUtil.h"

class ComputeNonbondedPair : public ComputePatchPair, private ComputeNonbondedUtil {

public:
  ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[],
	ComputeNonbondedWorkArrays* _workArrays,
	int minPartition = 0, int maxPartition = 1, int numPartitions = 1);
  ~ComputeNonbondedPair();
  nonbonded params;
  GBISParamStruct gbisParams;

protected :
  virtual void initialize();
  virtual int noWork();
  virtual void doForce(CompAtom* p[2], CompAtomExt* pExt[2], Results* r[2]);
  Box<Patch,CompAtom> *avgPositionBox[2];
  // BEGIN LA
  Box<Patch,CompAtom> *velocityBox[2];
  // END LA

  Real *intRad[2];
  Box<Patch,Real> *intRadBox[2];//write
  Box<Patch,GBReal> *psiSumBox[2];//write
  Box<Patch,Real> *bornRadBox[2];//read
  Box<Patch,GBReal> *dEdaSumBox[2];//write
  Box<Patch,Real> *dHdrPrefixBox[2];//read
  static const int numGBISPairlists = 4;
  Pairlists gbisStepPairlists[numGBISPairlists];//lasts a step

  BigReal reductionData[reductionDataSize];
  SubmitReduction *reduction;
  SubmitReduction *pressureProfileReduction;
  BigReal *pressureProfileData;

  ComputeNonbondedWorkArrays* const workArrays;

  Pairlists pairlists;
  int pairlistsValid;
  BigReal pairlistTolerance;

  int minPart, maxPart, numParts;

};

#endif

