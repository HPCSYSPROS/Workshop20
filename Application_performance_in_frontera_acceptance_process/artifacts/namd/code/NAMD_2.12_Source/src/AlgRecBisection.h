/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ALGROB_H
#define ALGROB_H

//#include "elements.h"
#include "PatchMap.h"
#include "Rebalancer.h"

class AlgRecBisection : public Rebalancer 
{
private: 

class Partition {
public:
  int refno;
  double load;				// total load in this set
  int origin[3];
  int corner[3];
  int  count;
  int node, mapped;
public:
  Partition(): refno(0), load(0.0), node(-1), mapped(0) {};
  inline int operator==(const Partition &p) const { return origin[0]==p.origin[0] && origin[1]==p.origin[1] && origin[2]==p.origin[2] && corner[0]==p.corner[0] && corner[1]==p.corner[1] && corner[2]==p.corner[2]; }
};

typedef struct {
  int id;
  int v[3];
  double load;
  int  refno;
  int  tv;
} ComputeLoad;


typedef struct {
  int v;
  int id;
} VecArray;

enum {XDIR=0, YDIR, ZDIR};

ComputeLoad *computeLoad;
VecArray  *(vArray[3]);
Partition *partitions;
Partition top_partition;
int npartition;
int currentp, refno;

void strategy();
void rec_divide(int, Partition&);
void setVal(int x, int y, int z);
int sort_partition(int x, int p, int r);
void qsort(int x, int p, int r);
void quicksort(int x);
void mapPartitionsToNodes();


public:
AlgRecBisection(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes);
};

#endif




