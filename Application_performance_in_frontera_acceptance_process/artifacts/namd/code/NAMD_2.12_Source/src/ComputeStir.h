/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000,2001 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* Barry Isralewitz, July 2001 */

#ifndef COMPUTESTIR_H
#define COMPUTESTIR_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"
#include "Vector.h"
#include "Tensor.h"


class ComputeStir : public ComputeHomePatch
{

 public:
	ComputeStir(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeStir();			//  Destructor
	
	virtual void doForce(FullAtom* p, Results* r);
	
	SubmitReduction *reduction;
 private:  

	Bool stirOn;
	BigReal omega; //rotation speed in radians/timestep
	BigReal startingTheta; //rotation offset in radians, normally zero
                         // even for restarted runs 
                         // since offset is automatically calculated 
                         // as omega * firsttimestep 
                         // useful to change speeds midrun
                         // or start torque mid-run, without losing
                         // track of timestep for handy bookeeping
	Tensor matMult (Tensor , Tensor );   //multiplies two matrices
	Tensor arbRotMat (BigReal theta); //make arbRotation matrix
                                          // needs pre/post translation
	Vector axisUnit;    //rotation axis direction
	Vector pivot;      //rotation axis pivot
	void printTensor (Tensor &);  //for debug printing
	void initVars ();     //initialize the matrices needed
	BigReal distanceToRay(Vector,Vector,Vector);// find dist. to Ray
	BigReal findHeight(Vector);      //find height along vector
	Vector projectionOnRay (Vector,Vector,Vector); //find projection on ray
	Vector placeThetaRadius (BigReal,BigReal,BigReal); //find point given theta
	                                            //and radius
	  
	BigReal findTheta (Vector refPos);
	Tensor leftMat, rightMat; // for arbitary rotations and related
                                  // coordinate operations

	  ////	BigReal* startTheta[10000];
	  ////[Node::Object()->configList->find("numStirredAtoms")];
	  ////array of starting thetas for
	  ////stirred atoms.  Ideally, do this more cleverly so so storage size varies
	  ////based on numStirredAtoms.
};

#endif







