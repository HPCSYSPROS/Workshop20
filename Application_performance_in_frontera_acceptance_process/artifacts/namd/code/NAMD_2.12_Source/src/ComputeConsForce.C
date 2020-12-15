
#include "InfoStream.h"
#include "SimParameters.h"
#include "ComputeConsForce.h"
#include "Molecule.h"
#include "Node.h"
#include "HomePatch.h"

//////////////////////////////////////////////////////////
// compute constant force
ComputeConsForce::ComputeConsForce(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeConsForce::~ComputeConsForce()
{
  delete reduction;
}

void ComputeConsForce::doForce(FullAtom* p, Results* r)
{ int localID,forceID;
  Molecule *molecule = Node::Object()->molecule;
  BigReal scaling = Node::Object()->simParameters->consForceScaling;
  int32 *index = molecule->consForceIndexes;  // Indexes into the force array
  Vector *cf = molecule->consForce;  // Force array
  Vector *forces = r->f[Results::normal];
  Force extForce = 0.;
  Tensor extVirial;

  for (localID=0; localID<numAtoms; ++localID) {
    // When the index is -1, it means there's no constant force on this atom
    if ((forceID=index[p[localID].id]) != -1) {
      Vector sf = scaling * cf[forceID];
      forces[localID] += sf;
      extForce += sf;
      Position vpos = homePatch->lattice.reverse_transform(
		p[localID].position, p[localID].transform );
      extVirial += outer(sf,vpos);
    }
  }

  ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
  ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
  reduction->submit();
}

//////////////////////////////////////////////////////////
// compute "constant" torque
ComputeConsTorque::ComputeConsTorque(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeConsTorque::~ComputeConsTorque()
{
  delete reduction;
}

void ComputeConsTorque::doForce(FullAtom* p, Results* r)
{ int localID,torqueID;
  Molecule *molecule = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;

  int32 *index = molecule->consTorqueIndexes;  // Indexes into the torque array
  Vector *forces = r->f[Results::normal];
  Force extForce = 0.;
  Tensor extVirial;
  const BigReal consTorqueGlobVal = simParams->consTorqueGlobVal;
  BigReal consTorqueVal;
  Vector consTorqueAxis, consTorquePivot;  
  Vector atomRadius;
  Vector torque;

  for (localID=0; localID<numAtoms; ++localID) {
    // When the index is -1, it means there's no constant torque on this atom
    if ((torqueID=index[p[localID].id]) != -1) {
      // compute the torqueing force and add it to the net force
      molecule->get_constorque_params(consTorqueVal, consTorqueAxis, consTorquePivot, p[localID].id);
      consTorqueAxis /= consTorqueAxis.length();
      atomRadius = p[localID].position - consTorquePivot;
      torque = cross(consTorqueAxis, atomRadius) * consTorqueVal * consTorqueGlobVal;
      forces[localID] += torque;
      extForce += torque;
      Position vpos = homePatch->lattice.reverse_transform(
		p[localID].position, p[localID].transform );
      extVirial += outer(torque,vpos);
    }
  }

  ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
  ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
  reduction->submit();
}
