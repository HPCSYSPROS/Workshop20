/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTETHOLE_H
#define COMPUTETHOLE_H

#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"
#include "ReductionMgr.h"

class Molecule;
//class TholeValue;
typedef Thole TholeValue;

class TholeElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(TholeElem*, int, BigReal*, BigReal *);

    static void getMoleculePointers(Molecule*, int*, int32***, Thole**);
    static void getParameterPointers(Parameters*, const TholeValue**);
    static void getTupleInfo(AtomSignature* sig, int *count, TupleSignature** t) {
        NAMD_die("Can't use Thole with memory optimized version of NAMD.");
	// *count = sig->tholeCnt;
	// *t = sig->tholeSigs;
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    // Internal data
    const TholeValue *value;

  int hash() const { 
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }

  enum { tholeEnergyIndex, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_THOLE_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline TholeElem();
  inline TholeElem(AtomID atom0, const TupleSignature *sig, const TholeValue *v);
  inline TholeElem(const Thole *a, const TholeValue *v);
  inline TholeElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3);
  ~TholeElem() {};

  inline int operator==(const TholeElem &a) const;
  inline int operator<(const TholeElem &a) const;
};

class ComputeThole : public ComputeHomeTuples<TholeElem,Thole,TholeValue>
{
public:

  ComputeThole(ComputeID c, PatchIDList &p) : ComputeHomeTuples<TholeElem,Thole,TholeValue>(c,p) { ; }

};

class ComputeSelfThole : public ComputeSelfTuples<TholeElem,Thole,TholeValue>
{
public:

  ComputeSelfThole(ComputeID c, PatchID p) : ComputeSelfTuples<TholeElem,Thole,TholeValue>(c,p) { ; }

};

#include "ComputeThole.inl"

#endif

