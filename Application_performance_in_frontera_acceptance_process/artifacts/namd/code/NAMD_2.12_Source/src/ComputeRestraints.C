/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeRestraints.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "Patch.h"
#include "NamdOneTools.h"

/************************************************************************/
/*									*/
/*			FUNCTION ComputeRestraints			*/
/*									*/
/************************************************************************/

ComputeRestraints::ComputeRestraints(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{
	reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

	SimParameters *simParams = Node::Object()->simParameters;

	//  Get parameters from the SimParameters object
	consExp = simParams->constraintExp;
	
	//****** BEGIN selective restraints (X,Y,Z) changes 
	consSelectOn = simParams->selectConstraintsOn;
	if (consSelectOn) {
	  consSelectX = simParams->constrXOn;
	  consSelectY = simParams->constrYOn;
	  consSelectZ = simParams->constrZOn;
	}
	//****** END selective restraints (X,Y,Z) changes 
	
	//****** BEGIN moving constraints changes 
	consMoveOn = simParams->movingConstraintsOn;
	if (consMoveOn) {
	  moveVel = simParams->movingConsVel;
	}
	//****** END moving constraints changes 
	//****** BEGIN rotating constraints changes 
	consRotOn = simParams->rotConstraintsOn;
	if (consRotOn) {
	  rotVel = simParams->rotConsVel;
	  rotAxis = simParams->rotConsAxis;
	  rotPivot = simParams->rotConsPivot;
	}
	//****** END rotating constraints changes 

}
/*			END OF FUNCTION ComputeRestraints		*/

/************************************************************************/
/*									*/
/*			FUNCTION ~ComputeRestraints			*/
/*									*/
/*	This is the destructor for the ComputeRestraints force object.	*/
/*									*/
/************************************************************************/

ComputeRestraints::~ComputeRestraints()

{
	delete reduction;
}
/*			END OF FUNCTION ~ComputeRestraints		*/

/************************************************************************/
/*									*/
/*				FUNCTION force				*/
/*									*/
/************************************************************************/
void ComputeRestraints::doForce(FullAtom* p, Results* res)
{
	Molecule *molecule = Node::Object()->molecule;
	Real k;			//  Force constant
	SimParameters *simParams = Node::Object()->simParameters;
	BigReal scaling = simParams->constraintScaling;
	Vector refPos;		//  Reference position
	BigReal r, r2; 	//  r=distance between atom position and the
			//  reference position, r2 = r^2
	Vector Rij;	//  vector between current position and reference pos
	BigReal value;	//  Current calculated value

	// aliases to work with old code
	Force *f = res->f[Results::normal];
	BigReal energy = 0;
	BigReal m[9];
	Tensor virial;
	Vector netForce = 0;

	// BEGIN moving and rotating constraint changes ******

	// This version only allows one atom to be moved 
	// and only ALL ref positions to be rotated

	int currentTime = patch->flags.step;
	if (consRotOn) {
	  vec_rotation_matrix(rotVel * currentTime, rotAxis, m);
	}

	// END moving and rotating constraint changes ******

	  
	if (scaling != 0.) for (int localID=0; localID<numAtoms; ++localID)
	{
	  if (molecule->is_atom_constrained(p[localID].id))
	  {
#ifndef MEM_OPT_VERSION
	    molecule->get_cons_params(k, refPos, p[localID].id);
#else
	    k = 1.0;
	    refPos = p[localID].fixedPosition;
#endif

	    k *= scaling;

	    // BEGIN moving and rotating constraint changes ******
	    
	    if (consMoveOn) {
	      refPos += currentTime * moveVel;
	    }
	    else if(consRotOn) {
	      refPos = mat_multiply_vec(refPos - rotPivot, m) + rotPivot;
	    }

	    // END moving and rotating constraint changes *******

	    if (simParams->sphericalConstraintsOn) {
	      BigReal refRad = (refPos - simParams->sphericalConstrCenter).length();
	      Vector relPos = patch->lattice.delta(p[localID].position, simParams->sphericalConstrCenter);
	      refPos = simParams->sphericalConstrCenter + relPos * (refRad/relPos.length());
            }

	    Rij = patch->lattice.delta(refPos,p[localID].position);
	    Vector vpos = refPos - Rij;

	    //****** BEGIN selective restraints (X,Y,Z) changes 
	    if (consSelectOn) { // check which components we want to restrain:
	      if (!consSelectX) {Rij.x=0.0;}  // by setting the appropriate
	      if (!consSelectY) {Rij.y=0.0;}  // Cartesian component to zero
	      if (!consSelectZ) {Rij.z=0.0;}  // we will avoid a restoring force
	    }                                 // along that component, and the
	                                      // energy will also be correct
	    //****** END selective restraints (X,Y,Z) changes 

	    
	    //  Calculate the distance and the distance squared
	    r2 = Rij.length2();
	    r = sqrt(r2);

	    //  Only calculate the energy and the force if the distance is
	    //  non-zero.   Otherwise, you could end up dividing by 0, which
	    //  is bad
	    if (r>0.0)
	    {
	      value=k;
      
	      //  Loop through and multiple k by r consExp times.
	      //  i.e.  calculate kr^e
	      //  I know, I could use pow(), but I don't trust it.
	      for (int k=0; k<consExp; ++k)
	      {
		value *= r;
	      }

	      //  Add to the energy total
	      energy += value;
      
	      //  Now calculate the force, which is ekr^(e-1).  Also, divide
	      //  by another factor of r to normalize the vector before we
	      //  multiple it by the magnitude
	      value *= consExp;
	      value /= r2;
      
	      Rij *= value;
              //iout << iINFO << "restraining force" << Rij        << "\n"; 

	      f[localID] += Rij;
	      netForce += Rij;
	      virial += outer(Rij,vpos);
	    }
	  }
	}

	reduction->item(REDUCTION_BC_ENERGY) += energy;
	ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
	ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,netForce);
	reduction->submit();

}
/*			END OF FUNCTION force				*/

