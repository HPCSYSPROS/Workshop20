/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

#include "common.h"
#include "Vector.h"
#ifndef __CUDACC__
#include "ResizeArray.h"
#endif

class Patch;
class Compute;

typedef Vector Position;
typedef Vector Velocity;

//#ifdef ARCH_POWERPC
//typedef AlignVector Force;
//#else
typedef Vector Force;
//#endif

typedef int AtomID;
typedef int AtomType;
typedef float Mass;
typedef float Charge;

typedef double Coordinate;

struct Transform
{
  signed char i,j,k;
  Transform(void) { i=0; j=0; k=0; }
};

/*
 * 1. "position" field in this structure is very important since it
 * needs to be sent to every patch after every timestep.
 * 2. Anything that is static (value is decided before computation)
 * or only changes after atom migration should be put into the CompAtomExt structure
 * 3. This data structure is 32-byte long which is particularly optimized for some machines
 * (including BG/L) for better cache and message performance. Therefore, changes
 * to this structure should be cautioned for the sake of performance.
 */

struct CompAtom {
  Position position;
  Charge charge;
  short vdwType;
  unsigned char partition;
  unsigned int nonbondedGroupSize : 3;
  unsigned int hydrogenGroupSize : 4;  // could be 3 if unsigned
  unsigned int isWater : 1;  // 0 = particle is not in water, 1 = is in water
};

#ifdef NAMD_KNL
struct CompAtomFlt {
  FloatVector position;
  int32 vdwType;
};
#endif

//CompAtomExt is now needed even in normal case
//for changing the packed msg type related to
//ProxyPatch into varsize msg type where
// two types of proxy msgs (originally, the msg 
// for the step where atoms migrate (ProxyAllMsg), 
// and  the msg for normal steps (ProxyDataMsg))
// are declared as the same class (ProxyDataMsg).
// Note that in normal case, the class is needed
// just for passing the compilation, but not involved
// in the actual force calculation.
// --Chao Mei

typedef int SigIndex;
typedef int AtomSigID;
typedef int ExclSigID;

struct CompAtomExt {
  #ifdef MEM_OPT_VERSION
  AtomSigID sigId;
  ExclSigID exclId;
  #endif
  #if defined(NAMD_CUDA) || defined(NAMD_MIC)
  int sortOrder;  // used to reorder atoms for CUDA
  #endif
  int id : 30;  // minimum for 100M atoms is 28 signed, 27 unsigned
  unsigned int atomFixed : 1;
  unsigned int groupFixed : 1;
};

struct FullAtom : CompAtom, CompAtomExt{
  Velocity velocity;
  Position fixedPosition;
  Mass mass;
  union{
      Real langevinParam;
#ifdef MEM_OPT_VERSION
      int hydVal;
#endif      
  };  
  int32 status;
  Transform transform;
  int migrationGroupSize;
  Real rigidBondLength;

#ifdef MEM_OPT_VERSION
  int outputRank;
#endif

#ifdef MEM_OPT_VERSION
  //a HACK to re-sort FullAtom list used in Parallel IO
  //When every home patch processor receives its atoms list for a patch,
  //the atoms inside this patch may not sorted according to hydList value
  //To save space, use anonymous union data structure to share the space
  //of "langevinParam" to store "hydList" from an InputAtom and then sort the 
  //atom list. The "langevinParam" value is not initialized until home 
  //patch creation -Chao Mei
  int operator < (const FullAtom &a) const {
      return hydVal < a.hydVal;
  }
#endif
};

//InputAtom is used to contain the info of the atoms
//loaded into input processors.
struct InputAtom: FullAtom{
	bool isValid;
	short isGP;
	short isMP;
	int hydList;
	int GPID;
	int MPID;
    	
	int operator < (const InputAtom &a) const{
		return hydList < a.hydList;
	}
};

struct CudaAtom {
  float x,y,z,q;
};

struct CudaForce {
  float x, y, z;
};

#ifndef __CUDACC__
typedef ResizeArray<CudaAtom> CudaAtomList;
typedef ResizeArray<CompAtom> CompAtomList;
typedef ResizeArray<CompAtomExt> CompAtomExtList;
#ifdef NAMD_KNL
typedef ResizeArray<CompAtomFlt> CompAtomFltList;
#endif
typedef ResizeArray<FullAtom> FullAtomList;
typedef ResizeArray<InputAtom> InputAtomList;
typedef ResizeArray<Position> PositionList;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArray<Force> ForceList;
typedef ResizeArray<Transform> TransformList;

typedef ResizeArray<AtomID> AtomIDList;
typedef ResizeArray<BigReal> BigRealList;
typedef ResizeArray<Real> RealList;
typedef float GBReal;
typedef ResizeArray<GBReal> GBRealList;
typedef ResizeArray<int> IntList;

typedef int PatchID;
typedef int ComputeID;
typedef int NodeID;

typedef ResizeArray<PatchID> PatchIDList;
typedef ResizeArray<Patch *> PatchList;

typedef ResizeArray<Compute *> ComputeList;

// See AtomMap
struct LocalID
{
  PatchID pid;
  int index;
};

typedef ResizeArray<NodeID> NodeIDList;

struct ExtForce {
  int replace;
  Force force;
  ExtForce() : replace(0) {;}
};


// DMK - Atom Sort
#if NAMD_ComputeNonbonded_SortAtoms != 0

  typedef struct __sort_entry {
    int index;  // Index of atom in CompAtom array
    BigReal sortValue;   // Distance of PAp from P0 (see calculation code)
  } SortEntry;

#endif

//This class represents a tree node of proxy spanning tree
//All pes in this array have the same "nodeID". In other words,
//all those pes are in the same physical node.
//This is a structure for adapting NAMD to multicore processors
struct proxyTreeNode{
    int nodeID;
    int *peIDs;
    int numPes;

    proxyTreeNode(){
        nodeID = -1;
        peIDs = NULL;
        numPes = 0;
    }
    proxyTreeNode(int nid, int numPes_, int *pes){
        nodeID = nid;
        numPes = numPes_;
        peIDs = new int[numPes];
        memcpy(peIDs, pes, sizeof(int)*numPes);
    }

    inline proxyTreeNode(const proxyTreeNode &n){
        nodeID = n.nodeID;
        numPes = n.numPes;
        if(numPes==0) {
            peIDs = NULL;
        }else{
            peIDs = new int[n.numPes];
            memcpy(peIDs, n.peIDs, sizeof(int)*numPes);
        }
    }
    inline proxyTreeNode &operator=(const proxyTreeNode &n){
        nodeID = n.nodeID;
        numPes = n.numPes;
        delete [] peIDs;
        if(numPes==0) {
            peIDs = NULL;
            return (*this);
        }
        peIDs = new int[n.numPes];
        memcpy(peIDs, n.peIDs, sizeof(int)*numPes);
        return (*this);
    }
    ~proxyTreeNode(){
        delete [] peIDs;
    }
};

typedef ResizeArray<proxyTreeNode> proxyTreeNodeList;
#endif // __CUDACC__

#endif /* NAMDTYPES_H */

