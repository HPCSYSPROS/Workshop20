/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/


#ifndef COMPUTEGROMACSPAIR_H
#define COMPUTEGROMACSPAIR_H

#include "common.h"
#include "NamdTypes.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"

class TuplePatchElem;

//class Molecule;
//class GromacsPairValue;

class GromacsPairElem { 
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(GromacsPairElem*, int, BigReal*, BigReal *);
    //void doWork();

    static void getMoleculePointers(Molecule*, int*, int32***, GromacsPair**);
    //static void getParameterPointers(Parameters*, const GromacsPairValue**);
    static void getParameterPointers(Parameters*, const GromacsPairValue**);
    static void getTupleInfo(AtomSignature* sig, int *count, TupleSignature** t) {
	*count = sig->gromacsPairCnt;
	*t = sig->gromacsPairSigs;
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    int hash() const { return 0x7FFFFFFF & ( (atomID[0]<<16) + (atomID[1])); }

    // Internal data
    const GromacsPairValue *value;

  enum { gromacsPairEnergyIndex, TENSOR(virialIndex), reductionDataSize };

  enum { reductionChecksumLabel = REDUCTION_GRO_LJ_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);
  inline GromacsPairElem();
  inline GromacsPairElem(AtomID atom0, const TupleSignature *sig, const GromacsPairValue *v);
  inline GromacsPairElem(const GromacsPair *a, const GromacsPairValue *v);
  inline GromacsPairElem(AtomID atom0, AtomID atom1);
  ~GromacsPairElem() {};

  inline int operator==(const GromacsPairElem &a) const;
  inline int operator<(const GromacsPairElem &a) const;
};

class ComputeGromacsPair : public ComputeHomeTuples<GromacsPairElem,GromacsPair,GromacsPairValue>
{
public:

  ComputeGromacsPair(ComputeID c, PatchIDList &p) : ComputeHomeTuples<GromacsPairElem,GromacsPair,GromacsPairValue>(c,p) { ; }

};

class ComputeSelfGromacsPair : public ComputeSelfTuples<GromacsPairElem,GromacsPair,GromacsPairValue>
{
public:

  ComputeSelfGromacsPair(ComputeID c, PatchID p) : ComputeSelfTuples<GromacsPairElem,GromacsPair,GromacsPairValue>(c,p) { ; }

};

#include "ComputeGromacsPair.inl"

#endif
