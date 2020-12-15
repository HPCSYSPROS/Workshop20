/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEGRIDFORCE_H
#define COMPUTEGRIDFORCE_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"
#include "GridForceGrid.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "Molecule.h"

class ComputeGridForce : public ComputeHomePatch
{
protected:
    template <class T> void do_calc(T *grid, int gridnum, FullAtom *p, int numAtoms, Molecule *mol, Force *forces, BigReal &energy, Force &extForce, Tensor &extVirial);

public:
    ComputeGridForce(ComputeID c, PatchID pid); 	//  Constructor
    virtual ~ComputeGridForce();			//  Destructor
    
    void doForce(FullAtom* p, Results* r);
    
    SubmitReduction *reduction;
};

#endif
