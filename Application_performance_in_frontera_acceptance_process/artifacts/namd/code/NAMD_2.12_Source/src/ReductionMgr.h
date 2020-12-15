/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ProcessorPrivate.h"

#define VECTOR(A) A ## _X, A ## _Y, A ## _Z
#define TENSOR(A) A ## _XX, A ## _XY, A ## _XZ, \
                  A ## _YX, A ## _YY, A ## _YZ, \
                  A ## _ZX, A ## _ZY, A ## _ZZ

#define ADD_VECTOR(R,RL,D,DL) \
  R->item( RL ## _X ) += D[ DL ## _X ]; \
  R->item( RL ## _Y ) += D[ DL ## _Y ]; \
  R->item( RL ## _Z ) += D[ DL ## _Z ]

#define ADD_VECTOR_OBJECT(R,RL,D) \
  R->item( RL ## _X ) += D.x; \
  R->item( RL ## _Y ) += D.y; \
  R->item( RL ## _Z ) += D.z

#define ADD_TENSOR(R,RL,D,DL) \
  R->item( RL ## _XX) += D[ DL ## _XX ]; \
  R->item( RL ## _XY) += D[ DL ## _XY ]; \
  R->item( RL ## _XZ) += D[ DL ## _XZ ]; \
  R->item( RL ## _YX) += D[ DL ## _YX ]; \
  R->item( RL ## _YY) += D[ DL ## _YY ]; \
  R->item( RL ## _YZ) += D[ DL ## _YZ ]; \
  R->item( RL ## _ZX) += D[ DL ## _ZX ]; \
  R->item( RL ## _ZY) += D[ DL ## _ZY ]; \
  R->item( RL ## _ZZ) += D[ DL ## _ZZ ]

#define ADD_TENSOR_OBJECT(R,RL,D) \
  R->item( RL ## _XX) += D.xx; \
  R->item( RL ## _XY) += D.xy; \
  R->item( RL ## _XZ) += D.xz; \
  R->item( RL ## _YX) += D.yx; \
  R->item( RL ## _YY) += D.yy; \
  R->item( RL ## _YZ) += D.yz; \
  R->item( RL ## _ZX) += D.zx; \
  R->item( RL ## _ZY) += D.zy; \
  R->item( RL ## _ZZ) += D.zz

#define GET_VECTOR(O,R,A) \
  O.x = R->item( A ## _X ); \
  O.y = R->item( A ## _Y ); \
  O.z = R->item( A ## _Z )

#define GET_TENSOR(O,R,A) \
  O.xx = R->item( A ## _XX); \
  O.xy = R->item( A ## _XY); \
  O.xz = R->item( A ## _XZ); \
  O.yx = R->item( A ## _YX); \
  O.yy = R->item( A ## _YY); \
  O.yz = R->item( A ## _YZ); \
  O.zx = R->item( A ## _ZX); \
  O.zy = R->item( A ## _ZY); \
  O.zz = R->item( A ## _ZZ)

typedef enum
{
 // energy
  REDUCTION_ANGLE_ENERGY,
  REDUCTION_BOND_ENERGY,
  REDUCTION_BONDED_ENERGY_F,
  REDUCTION_BONDED_ENERGY_TI_1,
  REDUCTION_BONDED_ENERGY_TI_2,
  REDUCTION_DIHEDRAL_ENERGY,
  REDUCTION_ELECT_ENERGY,
  REDUCTION_ELECT_ENERGY_F,
  REDUCTION_ELECT_ENERGY_TI_1,
  REDUCTION_ELECT_ENERGY_TI_2,
  REDUCTION_ELECT_ENERGY_SLOW,
  REDUCTION_ELECT_ENERGY_SLOW_F,
  REDUCTION_ELECT_ENERGY_SLOW_TI_1,
  REDUCTION_ELECT_ENERGY_SLOW_TI_2,
  REDUCTION_ELECT_ENERGY_PME_TI_1,
  REDUCTION_ELECT_ENERGY_PME_TI_2,
  REDUCTION_IMPROPER_ENERGY,
  // REDUCTION_THOLE_ENERGY - Drude model "correction" to electrostatic energy
  // REDUCTION_ANISO_ENERGY - Drude model add into bond energy
  REDUCTION_CROSSTERM_ENERGY,
  REDUCTION_HALFSTEP_KINETIC_ENERGY,
  REDUCTION_CENTERED_KINETIC_ENERGY,
  REDUCTION_INT_HALFSTEP_KINETIC_ENERGY,
  REDUCTION_INT_CENTERED_KINETIC_ENERGY,
  REDUCTION_DRUDECOM_CENTERED_KINETIC_ENERGY,
  REDUCTION_DRUDEBOND_CENTERED_KINETIC_ENERGY,
  REDUCTION_LJ_ENERGY,
  REDUCTION_LJ_ENERGY_F,
  REDUCTION_LJ_ENERGY_F_LEFT,
  REDUCTION_LJ_ENERGY_TI_1,
  REDUCTION_LJ_ENERGY_TI_2,
  REDUCTION_BC_ENERGY,
  REDUCTION_MISC_ENERGY,
  REDUCTION_GRO_LJ_ENERGY,
  REDUCTION_GRO_GAUSS_ENERGY,
  REDUCTION_GO_NATIVE_ENERGY,
  REDUCTION_GO_NONNATIVE_ENERGY,
 // pressure
  TENSOR(REDUCTION_VIRIAL_NORMAL),
  TENSOR(REDUCTION_VIRIAL_NBOND),
  TENSOR(REDUCTION_VIRIAL_SLOW),
  TENSOR(REDUCTION_VIRIAL_AMD_DIHE),
#ifdef ALTVIRIAL
  TENSOR(REDUCTION_ALT_VIRIAL_NORMAL),
  TENSOR(REDUCTION_ALT_VIRIAL_NBOND),
  TENSOR(REDUCTION_ALT_VIRIAL_SLOW),
#endif
  TENSOR(REDUCTION_INT_VIRIAL_NORMAL),
  TENSOR(REDUCTION_INT_VIRIAL_NBOND),
  TENSOR(REDUCTION_INT_VIRIAL_SLOW),
  VECTOR(REDUCTION_EXT_FORCE_NORMAL),
  VECTOR(REDUCTION_EXT_FORCE_NBOND),
  VECTOR(REDUCTION_EXT_FORCE_SLOW),
 // momentum
  VECTOR(REDUCTION_MOMENTUM),
  TENSOR(REDUCTION_MOMENTUM_SQUARED),    // Multigrator
  VECTOR(REDUCTION_ANGULAR_MOMENTUM),
  VECTOR(REDUCTION_HALFSTEP_MOMENTUM),
  REDUCTION_MOMENTUM_MASS,
 // used for minimization
  REDUCTION_MIN_F_DOT_F,
  REDUCTION_MIN_F_DOT_V,
  REDUCTION_MIN_V_DOT_V,
  REDUCTION_MIN_HUGE_COUNT,
 // used for pair interaction calculations
  VECTOR(REDUCTION_PAIR_VDW_FORCE),
  VECTOR(REDUCTION_PAIR_ELECT_FORCE),
 // checksum
  REDUCTION_ATOM_CHECKSUM,
  REDUCTION_COMPUTE_CHECKSUM,
  REDUCTION_BOND_CHECKSUM,
  REDUCTION_ANGLE_CHECKSUM,
  REDUCTION_DIHEDRAL_CHECKSUM,
  REDUCTION_IMPROPER_CHECKSUM,
  REDUCTION_THOLE_CHECKSUM,  // Drude model
  REDUCTION_ANISO_CHECKSUM,  // Drude model
  REDUCTION_CROSSTERM_CHECKSUM,
  REDUCTION_GRO_LJ_CHECKSUM,
  REDUCTION_EXCLUSION_CHECKSUM,
#ifdef NAMD_CUDA
  REDUCTION_EXCLUSION_CHECKSUM_CUDA,
#endif
  REDUCTION_MARGIN_VIOLATIONS,
  REDUCTION_PAIRLIST_WARNINGS,
  REDUCTION_STRAY_CHARGE_ERRORS,
 // semaphore (must be last)
  REDUCTION_MAX_RESERVED
} ReductionTag;

typedef enum {
  MULTIGRATOR_REDUCTION_KINETIC_ENERGY,
  TENSOR(MULTIGRATOR_REDUCTION_MOMENTUM_SQUARED),
  // semaphore
  MULTIGRATOR_REDUCTION_MAX_RESERVED
} MultigratorReductionTag;

// Later this can be dynamic
enum {
  REDUCTIONS_BASIC,
  REDUCTIONS_MINIMIZER,
  REDUCTIONS_PPROF_BONDED,
  REDUCTIONS_PPROF_NONBONDED,
  REDUCTIONS_PPROF_INTERNAL,
  REDUCTIONS_PPROF_KINETIC,
  REDUCTIONS_AMD,   // for accelMD
  REDUCTIONS_USER1,
  REDUCTIONS_USER2,
  REDUCTIONS_MULTIGRATOR,
 // semaphore (must be last)
  REDUCTION_MAX_SET_ID
};

// Later this can be dynamic
#define REDUCTION_MAX_CHILDREN 4

class ReductionRegisterMsg;
class ReductionSubmitMsg;
class SubmitReduction;
class RequireReduction;

// Queue element which stores data for a particular sequence number
class ReductionSetData {
public:
  int sequenceNumber;
  int submitsRecorded;
  BigReal *data;
  ReductionSetData *next;
  ReductionSetData(int seqNum, int size) {
    sequenceNumber = seqNum;
    submitsRecorded = 0;
    data = new BigReal[size];
    for ( int i = 0; i < size; ++i ) { data[i] = 0; }
    next = 0;
  }
  ~ReductionSetData() {
    delete [] data;
  }
};

// Stores the submit queue for a particular set of reductions
class ReductionSet {
public:
  int reductionSetID;
  int nextSequenceNumber;
  int submitsRegistered;
  int dataSize;
  ReductionSetData *dataQueue;
  ReductionSetData* getData(int seqNum);
  ReductionSetData* removeData(int seqNum);  // removes from queue
  int requireRegistered;  // is a thread subscribed on this node?
  int threadIsWaiting;  // is there a thread waiting on this?
  int waitingForSequenceNumber;  // sequence number waited for
  CthThread waitingThread;
  ReductionSet(int setID, int size,int numChildren);
  ~ReductionSet();
  int *addToRemoteSequenceNumber;
};

// Top level class
// don't derive from CBase_ReductionMgr to avoid .decl.h file
class ReductionMgr : public Group
{
private:
  friend class SubmitReduction;
  friend class RequireReduction;

  ReductionSet * (reductionSets[REDUCTION_MAX_SET_ID]);

  int myParent;  // parent node or -1 if none
#if 0
  int firstChild, lastChild;  // firstChild <= children < lastChild
  int isMyChild(int nodeID) const {
    return ( nodeID >= firstChild && nodeID < lastChild );
  }
#endif
  int numChildren;
  int *children;
  int isMyChild(int nodeID) const {
    int i;
    for(i=0;i<numChildren;i++)
      if (children[i]==nodeID) return 1;
    return 0;
  }      
  int childIndex(int nodeID) const {
    int i;
    for(i=0;i<numChildren;i++)
      if (children[i]==nodeID) return i;
    return -1;
  }
  int isRoot(void) const { return ( myParent == -1 ); }

  ReductionSet* getSet(int setID, int size);
  void delSet(int setID);

  void mergeAndDeliver(ReductionSet *set, int seqNum);

  void submit(SubmitReduction*);
  void remove(SubmitReduction*);

  void require(RequireReduction*);
  void remove(RequireReduction*);

public:

  // Singleton Access method
  inline static ReductionMgr *Object(void) {
    return CkpvAccess(ReductionMgr_instance);
  }

  ReductionMgr();
  ~ReductionMgr();

  void buildSpanTree(const int pe, 
                     const int max_intranode_children,
                     const int max_internode_children,
                     int* parent, 
                     int* num_children, 
                     int** children);

  // client interface
  SubmitReduction* willSubmit(int setID, int size = -1);
  RequireReduction* willRequire(int setID, int size = -1);

  // message entry points
  void remoteRegister(ReductionRegisterMsg *msg);
  void remoteUnregister(ReductionRegisterMsg *msg);
  void remoteSubmit(ReductionSubmitMsg *msg);

};

// Client handle for submissions
class SubmitReduction {
private:
  friend class ReductionMgr;
  int reductionSetID;
  int sequenceNumber;
  ReductionMgr *master;
  BigReal *data;  // managed explicitly by master
public:
  inline BigReal& item(int i) {
    return data[i];
  }
  inline void max(int i, BigReal v) {
    if ( v > data[i] ) {
      data[i] = v;
    }
  }
  void add(int nitems, const BigReal *arr) {
    for (int i=0; i<nitems; i++) data[i] += arr[i];
  }
  void submit(void) {
    master->submit(this);
  }
  ~SubmitReduction(void) { master->remove(this); }
};

// Client handle for requires
class RequireReduction {
private:
  friend class ReductionMgr;
  int reductionSetID;
  int sequenceNumber;
  ReductionMgr *master;
  ReductionSetData *currentData;
  BigReal *data;
  RequireReduction(void) { currentData = 0; data = 0; }
public:
  BigReal item(int i) const { return data[i]; }
  void require(void) {
    master->require(this);
  }
  ~RequireReduction(void) { delete currentData; master->remove(this); }
};


#endif

