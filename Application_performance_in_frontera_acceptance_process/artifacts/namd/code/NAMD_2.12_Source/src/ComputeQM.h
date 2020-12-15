/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEQM_H
#define COMPUTEQM_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;
class QMForceMsg;
class QMCoordMsg;
class ComputeQMAtom;

struct patchDataStrc {
    Box<Patch, CompAtom >* posBoxP;
    CompAtom* compAtomP;
    HomePatch* homePatchP;
    patchDataStrc(Box<Patch, CompAtom >* pbP, CompAtom* caP, HomePatch* hpP) {
        posBoxP = pbP;
        compAtomP = caP;
        homePatchP = hpP;
    }
};

struct LSSSubsDat {
    
    int origID, newID;
    int newVdWType;
    Real newCharge;
    
    LSSSubsDat() {}
    LSSSubsDat(const LSSSubsDat &ref) {
        origID = ref.origID;
        newID = ref.newID;
        newVdWType = ref.newVdWType;
        newCharge = ref.newCharge;
    }
    LSSSubsDat(int newOrigID) {
        origID = newOrigID;
    }
    LSSSubsDat(int newOrigID, int ID, int vdw, Real charge) {
        origID = newOrigID;
        newID = ID;
        newVdWType = vdw;
        newCharge = charge;
    }
    
    LSSSubsDat& operator=(const LSSSubsDat& ref) {
        origID = ref.origID;
        newID = ref.newID;
        newVdWType = ref.newVdWType;
        newCharge = ref.newCharge;
        return *this ;
    }
    
    bool operator<(const LSSSubsDat& ref) {
        return (origID < ref.origID);
    }
    bool operator==(const LSSSubsDat& ref) {
        return (origID == ref.origID) ;
    }
} ;

struct meMMQMGrp {
    
    int mmIndx;
    Real qmGrp;
    
    meMMQMGrp() {}
    meMMQMGrp(const meMMQMGrp &ref) {
        mmIndx = ref.mmIndx;
        qmGrp = ref.qmGrp;
    }
    meMMQMGrp(int newmmIndx) {
        mmIndx = newmmIndx;
    }
    meMMQMGrp(int newmmIndx, Real newQMGrp) {
        mmIndx = newmmIndx;
        qmGrp = newQMGrp;
    }
    
    meMMQMGrp& operator=(const meMMQMGrp& ref) {
        mmIndx = ref.mmIndx;
        qmGrp = ref.qmGrp;
        return *this ;
    }
    
    bool operator<(const meMMQMGrp& ref) {
        return (mmIndx < ref.mmIndx);
    }
    bool operator==(const meMMQMGrp& ref) {
        return (mmIndx == ref.mmIndx) ;
    }
} ;

class ComputeQM : public ComputeHomePatches {
public:
  ComputeQM(ComputeID c);
  virtual ~ComputeQM();
  
  void initialize();
  
  void doWork();
  
  void saveResults(QMForceMsg *);
  void processFullQM(QMCoordMsg *) ;

 private:
  SubmitReduction *reduction;
  
  SimParameters* simParams;
  Molecule *molPtr;
  
  int numQMAtms;
  int numQMGrps;
  const Real *qmAtomGroup ;
  const Real *qmGrpIDArray;
  const Real *qmAtmChrg ;
  const int *qmAtmIndx ;
  
  Bool noPC;
  int meNumMMIndx ;
  int *meMMindx;
  Real *meQMGrp;
  SortedArray< meMMQMGrp > meQMBonds;
  
  Bool customPC;
  ResizeArray< SortedArray<int> > customPCLists ;
  
  BigReal cutoff;
  
  ExtForce *oldForces;
  
  std::vector<patchDataStrc> patchData;
  
  // This is only used in case we have a stride in point charge selection.
  SortedArray<int> pcIDSortList ;
  
};


SortedArray<LSSSubsDat> &lssSubs(ComputeQMMgr *mgr) ;

#endif

