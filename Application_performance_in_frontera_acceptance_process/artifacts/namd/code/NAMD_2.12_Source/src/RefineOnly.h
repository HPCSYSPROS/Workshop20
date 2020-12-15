/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef _REFINEONLY_H_
#define _REFINEONLY_H_

#include "elements.h"
#include "Rebalancer.h"

class RefineOnly : public Rebalancer 
{
private: 
void strategy();


public:
RefineOnly(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes);
};

#endif
