/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEIMPROPERS_H
#define COMPUTEIMPROPERS_H

#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"
#include "ReductionMgr.h"

class Molecule;
class ImproperValue;

class ImproperElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    static void computeForce(ImproperElem*, int, BigReal*, BigReal *);

    static void getMoleculePointers(Molecule*, int*, int32***, Improper**);
    static void getParameterPointers(Parameters*, const ImproperValue**);
    static void getTupleInfo(AtomSignature* sig, int *count, TupleSignature** t) {
	*count = sig->improperCnt;
	*t = sig->improperSigs;
    }

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    // Internal data
    const ImproperValue *value;

  int hash() const {
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }
  enum { improperEnergyIndex, improperEnergyIndex_f, improperEnergyIndex_ti_1, 
         improperEnergyIndex_ti_2, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_IMPROPER_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  ImproperElem() { ; }

  ImproperElem(AtomID atom0, const TupleSignature *sig, const ImproperValue *v){
      atomID[0] = atom0;
      atomID[1] = atom0 + sig->offset[0];
      atomID[2] = atom0 + sig->offset[1];
      atomID[3] = atom0 + sig->offset[2];
      value = &v[sig->tupleParamType];
  }

  ImproperElem(const Improper *a, const ImproperValue *v) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    value = &v[a->improper_type];
  }

  ImproperElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3) {
    if (atom0 > atom3) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom3; atom3 = atom0; atom0 = tmp; 
      tmp = atom1; atom1 = atom2; atom2 = tmp;
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
    atomID[3] = atom3;
  }
  ~ImproperElem() {};

  int operator==(const ImproperElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2] && a.atomID[3] == atomID[3]);
  }

  int operator<(const ImproperElem &a) const {
    return  (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
            (atomID[2] < a.atomID[2] ||
            (atomID[2] == a.atomID[2] &&
             atomID[3] < a.atomID[3] 
	     ))))));
  }
};

class ComputeImpropers : public ComputeHomeTuples<ImproperElem,Improper,ImproperValue>
{
public:

  ComputeImpropers(ComputeID c, PatchIDList &p) : ComputeHomeTuples<ImproperElem,Improper,ImproperValue>(c,p) { ; }

};

class ComputeSelfImpropers : public ComputeSelfTuples<ImproperElem,Improper,ImproperValue>
{
public:

  ComputeSelfImpropers(ComputeID c, PatchID p) : ComputeSelfTuples<ImproperElem,Improper,ImproperValue>(c,p) { ; }

};

#endif

