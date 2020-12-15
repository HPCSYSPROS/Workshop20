/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTESPHERICALBC_H
#define COMPUTESPHERICALBC_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"

class ComputeSphericalBC : public ComputeHomePatch
{
private:
	BigReal r1;			//  Radius of first sphere
	BigReal r1_2;			//  Radius of first sphere squared
	BigReal k1;			//  First force constant
	BigReal r2;			//  Radius of second sphere (-1 if inactive)
	BigReal r2_2;			//  Raidus of second sphere squared
	BigReal k2;			//  Second force constant
	int exp1;			//  Exponent for first boundary condition
	int exp2;			//  Exponent for second boundary condition
	Bool twoForces;			//  Are there two potentials or just one
	BigReal energy;			//  Energy computed for the current timestep
	Vector center;			//  Center of spheres

public:
	ComputeSphericalBC(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeSphericalBC();			//  Destructor
	virtual void doForce(FullAtom* p, Results* r);
	SubmitReduction *reduction;

};

#endif







