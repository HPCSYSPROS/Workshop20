/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEANISO_H
#define COMPUTEANISO_H

#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"
#include "ReductionMgr.h"

class Molecule;
//class AnisoValue;
typedef Aniso AnisoValue;

class AnisoElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(AnisoElem*, int, BigReal*, BigReal *);

    static void getMoleculePointers(Molecule*, int*, int32***, Aniso**);
    static void getParameterPointers(Parameters*, const AnisoValue**);
    static void getTupleInfo(AtomSignature* sig, int *count, TupleSignature** t) {
        NAMD_die("Can't use Aniso with memory optimized version of NAMD.");
	// *count = sig->ansioCnt;
	// *t = sig->anisoSigs;
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    // Internal data
    const AnisoValue *value;

  int hash() const { 
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }

  enum { anisoEnergyIndex, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_ANISO_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline AnisoElem();
  inline AnisoElem(AtomID atom0, const TupleSignature *sig, const AnisoValue *v);
  inline AnisoElem(const Aniso *a, const AnisoValue *v);
  inline AnisoElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3);
  ~AnisoElem() {};

  inline int operator==(const AnisoElem &a) const;
  inline int operator<(const AnisoElem &a) const;
};

class ComputeAniso : public ComputeHomeTuples<AnisoElem,Aniso,AnisoValue>
{
public:

  ComputeAniso(ComputeID c, PatchIDList &p) : ComputeHomeTuples<AnisoElem,Aniso,AnisoValue>(c,p) { ; }

};

class ComputeSelfAniso : public ComputeSelfTuples<AnisoElem,Aniso,AnisoValue>
{
public:

  ComputeSelfAniso(ComputeID c, PatchID p) : ComputeSelfTuples<AnisoElem,Aniso,AnisoValue>(c,p) { ; }

};

#include "ComputeAniso.inl"

#endif

