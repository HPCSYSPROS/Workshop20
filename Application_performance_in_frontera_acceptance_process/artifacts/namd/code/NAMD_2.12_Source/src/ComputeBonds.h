/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEBOND_H
#define COMPUTEBOND_H

#include "common.h"
#include "NamdTypes.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"


class TuplePatchElem;

class BondElem {
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(BondElem*, int, BigReal*, BigReal *);

    static void getMoleculePointers(Molecule*, int*, int32***, Bond**);
    static void getParameterPointers(Parameters*, const BondValue**);
    static void getTupleInfo(AtomSignature* sig, int *count, TupleSignature** t) {
	*count = sig->bondCnt;
	*t = sig->bondSigs;
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    int hash() const { return 0x7FFFFFFF & ( (atomID[0]<<16) + (atomID[1])); }

    // Internal data
    const BondValue *value;

  enum { bondEnergyIndex, bondEnergyIndex_f, bondEnergyIndex_ti_1, 
         bondEnergyIndex_ti_2, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_BOND_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline BondElem();
  inline BondElem(AtomID atom0, const TupleSignature *sig, const BondValue *v);
  inline BondElem(const Bond *a, const BondValue *v);
  inline BondElem(AtomID atom0, AtomID atom1);
  ~BondElem() {};

  inline int operator==(const BondElem &a) const;
  inline int operator<(const BondElem &a) const;
};

class ComputeBonds : public ComputeHomeTuples<BondElem,Bond,BondValue>
{
public:

  ComputeBonds(ComputeID c, PatchIDList &p) : ComputeHomeTuples<BondElem,Bond,BondValue>(c,p) { ; }

};

class ComputeSelfBonds : public ComputeSelfTuples<BondElem,Bond,BondValue>
{
public:

  ComputeSelfBonds(ComputeID c, PatchID p) : ComputeSelfTuples<BondElem,Bond,BondValue>(c,p) { ; }

};

#include "ComputeBonds.inl"

#endif

