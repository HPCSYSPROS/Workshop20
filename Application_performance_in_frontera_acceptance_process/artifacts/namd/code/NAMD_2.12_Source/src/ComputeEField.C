/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeEField.h"
#include "Node.h"
#include "SimParameters.h"
#include "HomePatch.h"


ComputeEField::ComputeEField(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{

	reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}
/*			END OF FUNCTION ComputeEField		*/


ComputeEField::~ComputeEField()

{
	delete reduction;
}
/*			END OF FUNCTION ~ComputeEField		*/


void ComputeEField::doForce(FullAtom* p, Results* r) {

  SimParameters *simParams = Node::Object()->simParameters;
  Vector eField = simParams->eField;
  // Calculate the angular frequency in 1/fs.
  BigReal omega = TWOPI * simParams->eFieldFreq / 1000.;
  BigReal phi = PI/180.* simParams->eFieldPhase;
  BigReal t = patch->flags.step * simParams->dt;
  Vector eField1 = cos(omega * t - phi) * eField;

  const int normalized = simParams->eFieldNormalized;
  if ( normalized ) {
    Lattice &l = homePatch->lattice;
    eField1 = Vector(l.a_r()*eField1, l.b_r()*eField1, l.c_r()*eField1);
  }

  Force *forces = r->f[Results::normal];
  BigReal energy = 0;
  Force extForce = 0.;
  Tensor extVirial;

  //  Loop through and check each atom
  for (int i=0; i<numAtoms; i++) {
    Force force = p[i].charge * eField1; 
    forces[i] += force;
    Position vpos = homePatch->lattice.reverse_transform(
		p[i].position, p[i].transform );
    energy -= force * (vpos - homePatch->lattice.origin());
    if ( ! normalized ) {
      extForce += force;
      extVirial += outer(force,vpos);
    }
  }

  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  if ( ! normalized ) {
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
  }
  reduction->submit();

}
/*			END OF FUNCTION force				*/
