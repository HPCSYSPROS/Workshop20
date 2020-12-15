/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeGromacsPair.h"
#include "Molecule.h"
#include "Parameters.h"
#include "LJTable.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"
#include <iostream>

// static initialization
int GromacsPairElem::pressureProfileSlabs = 0;
int GromacsPairElem::pressureProfileAtomTypes = 1;
BigReal GromacsPairElem::pressureProfileThickness = 0;
BigReal GromacsPairElem::pressureProfileMin = 0;


void GromacsPairElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, GromacsPair** structarray)
{
#ifdef MEM_OPT_VERSION
    NAMD_die("Should not be called in GromacsPairElem::getMoleculePointers in memory optimized version!");
#else
    *count = mol->numLJPair;
    *byatom = mol->gromacsPairByAtom;
    *structarray = mol->gromacsPair;
#endif
}

//void GromacsPairElem::getParameterPointers(Parameters *p, const GromacsPairValue **v) {
void GromacsPairElem::getParameterPointers(Parameters *p, const GromacsPairValue **v) {
    // JLai
    *v=p->gromacsPair_array;
}

void GromacsPairElem::computeForce(GromacsPairElem *tuples, int ntuple, BigReal *reduction, 
                            BigReal *pressureProfileData)
{
 const Lattice & lattice = tuples[0].p[0]->p->lattice;

 for ( int ituple=0; ituple<ntuple; ++ituple ) {
  const GromacsPairElem &tup = tuples[ituple];
  enum { size = 2 };
  const AtomID (&atomID)[size](tup.atomID);
  const int    (&localIndex)[size](tup.localIndex);
  TuplePatchElem * const(&p)[size](tup.p);
  const Real (&scale)(tup.scale);
  const GromacsPairValue * const(&value)(tup.value);

    DebugM(1, "::computeforce() localIndex = " << localIndex [0] << " "
	   << localIndex[1] << std::endl);

    // compute vectors between atoms and their distances
    const Vector r12 = lattice.delta(p[0]->x[localIndex[0]].position,
				     p[1]->x[localIndex[1]].position);

    if ( p[0]->patchID == p[1]->patchID && localIndex[0] == localIndex[1] ) {
	continue;
    }
    SimParameters *simParams = Node::Object()->simParameters;
    BigReal cutoff = simParams->cutoff;
    BigReal r12_len = r12.length2();
    if ( r12_len == 0 ) continue;
    //if ( r12_len > cutoff) {
	//continue;
    //}
    BigReal ri2 = 1.0/r12_len;
    BigReal ri = sqrt(ri2);
    BigReal ri6 = ri2*ri2*ri2;
    BigReal ri7 = ri*ri6;
    BigReal ri8 = ri2*ri2*ri2*ri2;
    BigReal ri12 = ri6*ri6;
    BigReal ri13 = ri12*ri;
    BigReal ri14 = ri12*ri2;

    BigReal energy = 0;
    BigReal diff = 0;
    BigReal pairC12 = value->pairC12;
    BigReal pairC6 = value->pairC6;

    // Add the energy for the 12-6 LJ interaction to the total energy
    energy = (pairC12*ri12) - (pairC6*ri6);
    // This is a dirty hack; currently the code is looping over N^2 instead of N(N-1)/2 
    // This is happening because the LJ list is 2x long as it has to be
    //energy *= 0.5; 

    // Determine the magnitude of the force
    diff = ((12*pairC12*ri14) - (6*pairC6*ri8));
    // This is a dirty hack; currently the code is looping over N^2 instead of N(N-1)/2 
    // This is happening because the LJ list is 2x long as it has to be
    //diff *= 0.5; 
    //std::cout << "Force: " << diff << " " << pairC12 << " " << pairC6 << " " << r12.length() << "\n";

    //Scale the force vector accordingly
    const Force f12 = diff * r12;
    //std::cout << "Atoms: " << localIndex[0] << " " << localIndex[1] << " " << f12.length() << " " << f12 << "\n";
    //std::cout << "Force2: " << f12 << "\n";

    // Now add the forces to each force vector
    p[0]->f[localIndex[0]] += f12;
    p[1]->f[localIndex[1]] -= f12;

    DebugM(3, "::computeForce() -- ending with delta energy " << energy << std::endl);
    reduction[gromacsPairEnergyIndex] += energy;
    reduction[virialIndex_XX] += f12.x * r12.x;
    reduction[virialIndex_XY] += f12.x * r12.y;
    reduction[virialIndex_XZ] += f12.x * r12.z;
    reduction[virialIndex_YX] += f12.y * r12.x;
    reduction[virialIndex_YY] += f12.y * r12.y;
    reduction[virialIndex_YZ] += f12.y * r12.z;
    reduction[virialIndex_ZX] += f12.z * r12.x;
    reduction[virialIndex_ZY] += f12.z * r12.y;
    reduction[virialIndex_ZZ] += f12.z * r12.z;

    if (pressureProfileData) {
	BigReal z1 = p[0]->x[localIndex[0]].position.z;
	BigReal z2 = p[1]->x[localIndex[1]].position.z;
	int n1 = (int)floor((z1-pressureProfileMin)/pressureProfileThickness);
	int n2 = (int)floor((z2-pressureProfileMin)/pressureProfileThickness);
	pp_clamp(n1, pressureProfileSlabs);
	pp_clamp(n2, pressureProfileSlabs);
	int p1 = p[0]->x[localIndex[0]].partition;
	int p2 = p[1]->x[localIndex[1]].partition;
	int pn = pressureProfileAtomTypes;
	pp_reduction(pressureProfileSlabs,
		     n1, n2, 
		     p1, p2, pn, 
		     f12.x * r12.x, f12.y * r12.y, f12.z * r12.z,
		     pressureProfileData);
    } 

 }
}

void GromacsPairElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
    // JLai
    reduction->item(REDUCTION_GRO_LJ_ENERGY) += data[gromacsPairEnergyIndex];
    ADD_TENSOR(reduction,REDUCTION_VIRIAL_NBOND,data,virialIndex);
    return;
}
