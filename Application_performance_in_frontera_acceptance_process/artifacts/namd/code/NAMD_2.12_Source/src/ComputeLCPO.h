/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Solvent Accessible Surface Area calculation by LCPO
   Linear Combination of Pairwise Overlaps
*/

#ifndef COMPUTELCPO_H
#define COMPUTELCPO_H

#include "Compute.h"
#include "PatchTypes.h"
#include "Box.h"
#include "OwnerBox.h"
#include "ComputeNonbondedUtil.h"
#include "NamdTypes.h"

struct LCPOAtom {
  float x, y, z, r;
  Vector *f; // where to write force 
};

class LCPONeighborList {
  enum {initsize = 10};
  int *nnfa; // number of neighbors for atom
  int curAtom;
  int maxAtoms;
  LCPOAtom *neighbors;
  int curNeighbor;
  int maxNeighbors;
  LCPONeighborList(const LCPONeighborList&) { ; }
  LCPONeighborList& operator=(const LCPONeighborList&) { return *this; }
public:
  LCPONeighborList() :
      maxNeighbors(initsize), maxAtoms(initsize),
      curNeighbor(0), curAtom(0) {
    neighbors = new LCPOAtom[initsize];
    nnfa = new int[initsize];
  }
  ~LCPONeighborList() {
    delete [] neighbors;
    delete [] nnfa;
  }
  LCPOAtom *newlist(int max_size) {  // get a new list w/ room for max_size
    //do we need to make room for more neighbors
    int reqNewSize = curNeighbor + max_size;
    int newSize = maxNeighbors;
    while ( newSize < reqNewSize ) { newSize += newSize >> 1; }
    if ( newSize > maxNeighbors ) {
      LCPOAtom *newNeighbors = new LCPOAtom[newSize];
      CmiMemcpy(newNeighbors,neighbors,curNeighbor*sizeof(LCPOAtom));
      delete [] neighbors;
      neighbors = newNeighbors;
      maxNeighbors = newSize;
    }
    //do we need to make room for more atoms
    if (curAtom == maxAtoms) {
      newSize = maxAtoms + (maxAtoms >> 1);
      int *newNnfa = new int[newSize];
      CmiMemcpy(newNnfa,nnfa,curAtom*sizeof(int));
      delete [] nnfa;
      nnfa = newNnfa;
      maxAtoms = newSize;
    }
    return &neighbors[curNeighbor];
  }
  // don't specify size if previous allocation should have extra space
  LCPOAtom *newlist() {  // get a new list assuming already allocated
    return &neighbors[curNeighbor];
  }
  void newsize(int list_size) {  // set the size of the last list gotten
    nnfa[curAtom] = list_size;
    curAtom++;
    curNeighbor += list_size;
  }
  void reset() {  // go back to the beginning
    curNeighbor = 0;
    curAtom = 0;
  }
  void nextlist(LCPOAtom **list, int *list_size) {  // get next list and size
    *list = &neighbors[curNeighbor];
    *list_size = nnfa[curAtom];
    curNeighbor += nnfa[curAtom];
    curAtom ++;
  }
  int getSize() { return maxNeighbors; }
};




class Patch;
class Node;
class PatchMap;

class ComputeLCPO: public Compute, private ComputeNonbondedUtil {

public:
  ComputeLCPO(ComputeID c, PatchID pid[], int t[],
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions, int numPatches);

  virtual ~ComputeLCPO();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();
  virtual int noWork();

protected :
  int numAtoms[8];
  int valid[8][8];
  //0 if patch is invalid due to only 1 patch in that dimension
  int invalidPatch[8];
  CompAtomExt *posExt[8];
  CompAtom *pos[8];
  Results *force[8];
  int *lcpoType[8];
  int step;

  virtual void doForce();
  Patch *patch[8];

  PatchID patchID[8];
  int trans[8];
  Box<Patch,CompAtom> *positionBox[8];
  Box<Patch,Results> *forceBox[8];
  Box<Patch,int> *lcpoTypeBox[8];

  ComputeNonbondedWorkArrays* const workArrays;
  int minPart, maxPart, numParts;

  SubmitReduction *reduction;

  private:
    BigReal bounds[3][2];
    int periodic[3];
    int oob[3];
    Vector offset[8];
    int minIg[8];
    int strideIg;//stride through partitions

    //index "i" is patch; index "j" is valid atoms in patch
    Pairlists inAtomsPl;
    Real surfTen;
    Real maxAtomRadius;
    Real cut2;
    LCPONeighborList lcpoNeighborList;

    static const Real lcpoParams[23][5];

    int isInBounds( Real x, Real y, Real z );
};

#endif
