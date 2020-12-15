/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEDIHEDRALS_H
#define COMPUTEDIHEDRALS_H

#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"
#include "ReductionMgr.h"

class Molecule;
class DihedralValue;

class DihedralElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(DihedralElem*, int, BigReal*, BigReal *);

    static void getMoleculePointers(Molecule*, int*, int32***, Dihedral**);
    static void getParameterPointers(Parameters*, const DihedralValue**);
    static void getTupleInfo(AtomSignature* sig, int *count, TupleSignature** t) {
	*count = sig->dihedralCnt;
	*t = sig->dihedralSigs;
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    // Internal data
    const DihedralValue *value;

  int hash() const { 
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }

  enum { dihedralEnergyIndex, dihedralEnergyIndex_f, dihedralEnergyIndex_ti_1, 
         dihedralEnergyIndex_ti_2,TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_DIHEDRAL_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline DihedralElem();
  inline DihedralElem(AtomID atom0, const TupleSignature *sig, const DihedralValue *v);
  inline DihedralElem(const Dihedral *a, const DihedralValue *v);
  inline DihedralElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3);
  ~DihedralElem() {};

  inline int operator==(const DihedralElem &a) const;
  inline int operator<(const DihedralElem &a) const;
};

class ComputeDihedrals : public ComputeHomeTuples<DihedralElem,Dihedral,DihedralValue>
{
public:

  ComputeDihedrals(ComputeID c, PatchIDList &p) : ComputeHomeTuples<DihedralElem,Dihedral,DihedralValue>(c,p) { ; }

};

class ComputeSelfDihedrals : public ComputeSelfTuples<DihedralElem,Dihedral,DihedralValue>
{
public:

  ComputeSelfDihedrals(ComputeID c, PatchID p) : ComputeSelfTuples<DihedralElem,Dihedral,DihedralValue>(c,p) { ; }

};

#include "ComputeDihedrals.inl"

#endif

