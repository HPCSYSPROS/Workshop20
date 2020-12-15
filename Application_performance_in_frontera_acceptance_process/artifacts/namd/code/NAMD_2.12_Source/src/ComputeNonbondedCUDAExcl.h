/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDCUDAEXCL_H
#define COMPUTENONBONDEDCUDAEXCL_H

#include "common.h"
#include "NamdTypes.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"
#include "ComputeNonbondedUtil.h"


class TuplePatchElem;

class ExclElem : public ComputeNonbondedUtil {
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(ExclElem*, int, BigReal*, BigReal *);

    static void getMoleculePointers(Molecule*, int*, int32***, Exclusion**);
    static void getParameterPointers(Parameters*, const int**);
    static void getTupleInfo(ExclusionSignature* sig, int *count, TupleSignature** t) {
#ifdef NAMD_CUDA
	*count = sig->allExclCnt;
	*t = sig->allTuples;
#endif
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    int hash() const { return 0x7FFFFFFF & ( (atomID[0]<<16) + (atomID[1])); }

    // Internal data
    int modified;

  enum { vdwEnergyIndex, electEnergyIndex, fullElectEnergyIndex, TENSOR(virialIndex),
           TENSOR(slowVirialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_EXCLUSION_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline ExclElem();
  inline ExclElem(AtomID atom0, const TupleSignature *sig, const int *v);
  inline ExclElem(const Exclusion *a, const int *v);
  inline ExclElem(AtomID atom0, AtomID atom1);
  ~ExclElem() {};

  inline int operator==(const ExclElem &a) const;
  inline int operator<(const ExclElem &a) const;
};

class ComputeExcls : public ComputeHomeTuples<ExclElem,Exclusion,int>
{
public:

  ComputeExcls(ComputeID c, PatchIDList &p) : ComputeHomeTuples<ExclElem,Exclusion,int>(c,p) { ; }

};

class ComputeSelfExcls : public ComputeSelfTuples<ExclElem,Exclusion,int>
{
public:

  ComputeSelfExcls(ComputeID c, PatchID p) : ComputeSelfTuples<ExclElem,Exclusion,int>(c,p) { ; }

};

#include "ComputeNonbondedCUDAExcl.inl"

#endif

