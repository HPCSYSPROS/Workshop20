/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeNonbondedCUDAExcl.h"
#include "Molecule.h"
#include "Parameters.h"
#include "LJTable.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"


// static initialization
int ExclElem::pressureProfileSlabs = 0;
int ExclElem::pressureProfileAtomTypes = 1;
BigReal ExclElem::pressureProfileThickness = 0;
BigReal ExclElem::pressureProfileMin = 0;

void ExclElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Exclusion** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in ExclElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numExclusions;
  *byatom = mol->exclusionsByAtom;
  *structarray = mol->exclusions;
#endif
}

void ExclElem::getParameterPointers(Parameters *p, const int **v) {
  *v = 0;
}

void ExclElem::computeForce(ExclElem *tuples, int ntuple, BigReal *reduction, 
                            BigReal *pressureProfileData)
{
  const Lattice & lattice = tuples[0].p[0]->p->lattice;
  const Flags &flags = tuples[0].p[0]->p->flags;
  if ( ! flags.doNonbonded ) return;
  const int doFull = flags.doFullElectrostatics;
  const int doEnergy = flags.doEnergy;

  for ( int ituple=0; ituple<ntuple; ++ituple ) {
    const ExclElem &tup = tuples[ituple];
    enum { size = 2 };
    const AtomID (&atomID)[size](tup.atomID);
    const int    (&localIndex)[size](tup.localIndex);
    TuplePatchElem * const(&p)[size](tup.p);
    const Real (&scale)(tup.scale);
    const int (&modified)(tup.modified);

    const CompAtom &p_i = p[0]->x[localIndex[0]];
    const CompAtom &p_j = p[1]->x[localIndex[1]];

    // compute vectors between atoms and their distances
    const Vector r12 = lattice.delta(p_i.position, p_j.position);
    BigReal r2 = r12.length2();

    if ( r2 > cutoff2 ) continue;

    if ( modified && r2 < 1.0 ) r2 = 1.0;  // match CUDA interpolation

    r2 += r2_delta;

    union { double f; int64 i; } r2i;
    r2i.f = r2;
    const int r2_delta_expc = 64 * (r2_delta_exp - 1023);
    int table_i = (r2i.i >> (32+14)) + r2_delta_expc;  // table_i >= 0

    const BigReal* const table_four_i = table_noshort + 16*table_i;

    BigReal diffa = r2 - r2_table[table_i];

    BigReal fast_a = 0., fast_b = 0., fast_c = 0., fast_d = 0.;
    BigReal slow_a, slow_b, slow_c, slow_d;

  if ( modified ) {

    const LJTable::TableEntry * lj_pars =
            ljTable->table_row(p_i.vdwType) + 2 * p_j.vdwType;

    // modified - normal = correction
    const BigReal A = scaling * ( (lj_pars+1)->A - lj_pars->A );
    const BigReal B = scaling * ( (lj_pars+1)->B - lj_pars->B );

    BigReal vdw_d = A * table_four_i[0] - B * table_four_i[4];
    BigReal vdw_c = A * table_four_i[1] - B * table_four_i[5];
    BigReal vdw_b = A * table_four_i[2] - B * table_four_i[6];
    BigReal vdw_a = A * table_four_i[3] - B * table_four_i[7];

    const BigReal kqq = (1.0 - scale14) *
            COULOMB * p_i.charge * p_j.charge * scaling * dielectric_1;

    fast_a =      kqq * fast_table[4*table_i+0];  // not used!
    fast_b = 2. * kqq * fast_table[4*table_i+1];
    fast_c = 4. * kqq * fast_table[4*table_i+2];
    fast_d = 6. * kqq * fast_table[4*table_i+3];

    if ( doFull ) {
      slow_a =      kqq * slow_table[4*table_i+3];  // not used!
      slow_b = 2. * kqq * slow_table[4*table_i+2];
      slow_c = 4. * kqq * slow_table[4*table_i+1];
      slow_d = 6. * kqq * slow_table[4*table_i+0];
    }

    if ( doEnergy ) {
      reduction[vdwEnergyIndex] -=
		( ( diffa * (1./6.)*vdw_d + 0.25*vdw_c ) * diffa
			+ 0.5*vdw_b ) * diffa + vdw_a;
      reduction[electEnergyIndex] -=
		( ( diffa * (1./6.)*fast_d + 0.25*fast_c ) * diffa
			+ 0.5*fast_b ) * diffa + fast_a;
      if ( doFull ) {
        reduction[fullElectEnergyIndex] -=
		( ( diffa * (1./6.)*slow_d + 0.25*slow_c ) * diffa
			+ 0.5*slow_b ) * diffa + slow_a;
      }
    }

    fast_a += vdw_a;
    fast_b += vdw_b;
    fast_c += vdw_c;
    fast_d += vdw_d;

  } else if ( doFull ) {  // full exclusion

    const BigReal kqq = 
            COULOMB * p_i.charge * p_j.charge * scaling * dielectric_1;

    slow_d = kqq * ( table_four_i[8]  - table_four_i[12] );
    slow_c = kqq * ( table_four_i[9]  - table_four_i[13] );
    slow_b = kqq * ( table_four_i[10] - table_four_i[14] );
    slow_a = kqq * ( table_four_i[11] - table_four_i[15] );  // not used!

    if ( doEnergy ) {
      reduction[fullElectEnergyIndex] -=
		( ( diffa * (1./6.)*slow_d + 0.25*slow_c ) * diffa
			+ 0.5*slow_b ) * diffa + slow_a;
    }
  }

  register BigReal fast_dir =
                  (diffa * fast_d + fast_c) * diffa + fast_b;

  const Force f12 = fast_dir * r12;

  //  Now add the forces to each force vector
  p[0]->r->f[Results::nbond][localIndex[0]] += f12;
  p[1]->r->f[Results::nbond][localIndex[1]] -= f12;

  // reduction[nonbondedEnergyIndex] += energy;
  reduction[virialIndex_XX] += f12.x * r12.x;
  reduction[virialIndex_XY] += f12.x * r12.y;
  reduction[virialIndex_XZ] += f12.x * r12.z;
  reduction[virialIndex_YX] += f12.y * r12.x;
  reduction[virialIndex_YY] += f12.y * r12.y;
  reduction[virialIndex_YZ] += f12.y * r12.z;
  reduction[virialIndex_ZX] += f12.z * r12.x;
  reduction[virialIndex_ZY] += f12.z * r12.y;
  reduction[virialIndex_ZZ] += f12.z * r12.z;

  if ( doFull ) {
    register BigReal slow_dir =
                  (diffa * slow_d + slow_c) * diffa + slow_b;

    const Force slow_f12 = slow_dir * r12;

    p[0]->r->f[Results::slow][localIndex[0]] += slow_f12;
    p[1]->r->f[Results::slow][localIndex[1]] -= slow_f12;

    // reduction[nonbondedEnergyIndex] += energy;
    reduction[slowVirialIndex_XX] += slow_f12.x * r12.x;
    reduction[slowVirialIndex_XY] += slow_f12.x * r12.y;
    reduction[slowVirialIndex_XZ] += slow_f12.x * r12.z;
    reduction[slowVirialIndex_YX] += slow_f12.y * r12.x;
    reduction[slowVirialIndex_YY] += slow_f12.y * r12.y;
    reduction[slowVirialIndex_YZ] += slow_f12.y * r12.z;
    reduction[slowVirialIndex_ZX] += slow_f12.z * r12.x;
    reduction[slowVirialIndex_ZY] += slow_f12.z * r12.y;
    reduction[slowVirialIndex_ZZ] += slow_f12.z * r12.z;
  }

  }
}

void ExclElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ELECT_ENERGY) += data[electEnergyIndex];
  reduction->item(REDUCTION_LJ_ENERGY) += data[vdwEnergyIndex];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += data[fullElectEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NBOND,data,virialIndex);
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_SLOW,data,slowVirialIndex);
}

