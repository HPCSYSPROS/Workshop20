/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ALGNBOR_H
#define ALGNBOR_H

//#include "elements.h"
#include "Rebalancer.h"

class AlgNbor : public Rebalancer 
{
private: 
int  mype;
int  nNbors;
void strategy();


public:
AlgNbor(int pe, computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes, int nNbs);
};

#endif




