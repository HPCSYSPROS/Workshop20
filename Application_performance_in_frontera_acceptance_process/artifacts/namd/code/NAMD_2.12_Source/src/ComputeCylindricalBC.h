/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTECYLINDRICALBC_H
#define COMPUTECYLINDRICALBC_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"

class ComputeCylindricalBC : public ComputeHomePatch
{
private:
	char axis;			//  'x', 'y', or 'z'
	BigReal r1;			//  Radius of first cylinder
	BigReal r1_2;			//  Radius of first cylinder squared
	BigReal l1;			//  Length of First cylinder
	BigReal l1_2;			//  Length of first cylinder, squared
	BigReal k1;			//  First force constant
	BigReal r2;			//  Radius of second cylinder (-1 if inactive)
	BigReal r2_2;			//  Raidus of second cylinder squared
	BigReal k2;			//  Second force constant
	BigReal l2;			//  Length of second cylinder
	BigReal l2_2;			//  Length of second cylinder, squared
	int exp1;			//  Exponent for first boundary condition
	int exp2;			//  Exponent for second boundary condition
	Bool twoForces;			//  Are there two potentials or just one
	Vector center;			//  Center of cylinder

public:
	ComputeCylindricalBC(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeCylindricalBC();			//  Destructor

	virtual void doForce(FullAtom* p, Results* r);
	SubmitReduction *reduction;

};

#endif

