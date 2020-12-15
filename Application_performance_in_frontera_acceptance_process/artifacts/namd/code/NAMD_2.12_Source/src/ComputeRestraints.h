/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTERESTRAINTS_H
#define COMPUTERESTRAINTS_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"

class ComputeRestraints : public ComputeHomePatch
{
private:
	int consExp;		//  Exponent for energy function from SimParameters
	//****** BEGIN selective restraints (X,Y,Z) changes 
	Bool consSelectOn;      // Selection of Cartesian components active?
	Bool consSelectX, consSelectY,
	     consSelectZ;       // which components are active?
	//****** END selective restraints (X,Y,Z) changes 
	//****** BEGIN moving constraints changes 
	Bool consMoveOn;        //  Are the moving constraints on?
        Vector moveVel;         // velocity of the constraint movement (A/timestep).
	//****** END moving constraints changes 
	//****** BEGIN rotating constraints changes 
	// rotating constraints. 
	// Ref. pos. of all the atoms that are constrained will rotate
	Bool consRotOn;         // Are the rotating constraints on?
	Vector rotAxis;         // Axis of rotation
        Vector rotPivot;        // Pivot point of rotation
        BigReal rotVel;         // Rotation velocity (deg/timestep);
	//****** END rotating constraints changes 

public:
	ComputeRestraints(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeRestraints();			//  Destructor
	virtual void doForce(FullAtom* p, Results* r);
	SubmitReduction *reduction;

};

#endif

