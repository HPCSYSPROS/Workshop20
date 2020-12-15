/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ALG7_H
#define ALG7_H

//#include "elements.h"
#include "Rebalancer.h"

class Alg7 : public Rebalancer 
{
private: 
void strategy();
void togrid(processorInfo* goodP[3][3][2], processorInfo* poorP[3][3][2],
                        processorInfo *p, computeInfo *c);

public:
Alg7(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes);
};

#endif




