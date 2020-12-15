/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeGridForce.h"
#include "GridForceGrid.h"
#include "Node.h"

#define MIN_DEBUG_LEVEL 2
//#define DEBUGM
#include "Debug.h"

#include "GridForceGrid.inl"
#include "MGridforceParams.h"

//#define GF_FORCE_OUTPUT
//#define GF_FORCE_OUTPUT_FREQ 100
#define GF_OVERLAPCHECK_FREQ 1000


ComputeGridForce::ComputeGridForce(ComputeID c, PatchID pid)
    : ComputeHomePatch(c,pid)
{

    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}
/*			END OF FUNCTION ComputeGridForce	*/


ComputeGridForce::~ComputeGridForce()
{
    delete reduction;
}
/*			END OF FUNCTION ~ComputeGridForce	*/

template <class T> void ComputeGridForce::do_calc(T *grid, int gridnum, FullAtom *p, int numAtoms, Molecule *mol, Force *forces, BigReal &energy, Force &extForce, Tensor &extVirial)
{
    Real scale;			// Scaling factor
    Charge charge;		// Charge
    Vector dV;
    float V;
    
    Vector gfScale = grid->get_scale();

    //  Loop through and check each atom
    for (int i = 0; i < numAtoms; i++) {
	if (mol->is_atom_gridforced(p[i].id, gridnum))
	{
	    DebugM(1, "Atom " << p[i].id << " is gridforced\n" << endi);
	    
	    mol->get_gridfrc_params(scale, charge, p[i].id, gridnum);
	    
	    // Wrap coordinates using grid center
	    Position pos = grid->wrap_position(p[i].position, homePatch->lattice);
	    DebugM(1, "pos = " << pos << "\n" << endi);
	    
	    // Here's where the action happens
	    int err = grid->compute_VdV(pos, V, dV);
	    
	    if (err) {
		DebugM(2, "V = 0\n" << endi);
		DebugM(2, "dV = 0 0 0\n" << endi);
		continue;  // This means the current atom is outside the potential
	    }
	    
	    //Force force = scale * Tensor::diagonal(gfScale) * (-charge * dV);
	    Force force = -charge * scale * Vector(gfScale.x * dV.x, gfScale.y * dV.y, gfScale.z * dV.z);
	    
#ifdef DEBUGM
	    DebugM(2, "scale = " << scale << " gfScale = " << gfScale << " charge = " << charge << "\n" << endi);
	    
	    DebugM(2, "V = " << V << "\n" << endi);
	    DebugM(2, "dV = " << dV << "\n" << endi);
	    DebugM(2, "grid = " << gridnum << " force = " << force << " pos = " << pos << " V = " << V << " dV = " << dV << " step = " << homePatch->flags.step << " index = " << p[i].id << "\n" << endi);
	    
	    DebugM(1, "transform = " << (int)p[i].transform.i << " "
		   << (int)p[i].transform.j << " " << (int)p[i].transform.k << "\n" << endi);
	    
	    if (V != V) {
		iout << iWARN << "V is NaN!\natomid = " << p[i].id << " loc = " << p[i].position << " V = " << V << "\n" << endi;
	    }
#endif
	    
	    forces[i] += force;
	    extForce += force;
	    Position vpos = homePatch->lattice.reverse_transform(p[i].position, p[i].transform);
	    
	    //energy -= force * (vpos - homePatch->lattice.origin());
	    if (gfScale.x == gfScale.y && gfScale.x == gfScale.z)
	    {
		// only makes sense when scaling is isotropic
		energy += scale * gfScale.x * (charge * V);
		
		// add something when we're off the grid? I'm thinking no
	    }
	    extVirial += outer(force,vpos);
	}
    }
}

void ComputeGridForce::doForce(FullAtom* p, Results* r)
{
    SimParameters *simParams = Node::Object()->simParameters;
    Molecule *mol = Node::Object()->molecule;
    
    Force *forces = r->f[Results::normal];
    BigReal energy = 0;
    Force extForce = 0.;
    Tensor extVirial;
    
    int numAtoms = homePatch->getNumAtoms();

    if ( mol->numGridforceGrids < 1 ) NAMD_bug("No grids loaded in ComputeGridForce::doForce()");
    
    for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
	GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
	
	if (homePatch->flags.step % GF_OVERLAPCHECK_FREQ == 0) {
	    // only check on node 0 and every GF_OVERLAPCHECK_FREQ steps
	  if (simParams->langevinPistonOn || simParams->berendsenPressureOn) {
		// check for grid overlap if pressure control is on
		// not needed without pressure control, since the check is also performed on startup
      if (!grid->fits_lattice(homePatch->lattice)) {
        char errmsg[512];
        if (grid->get_checksize()) {
          sprintf(errmsg, "Warning: Periodic cell basis too small for Gridforce grid %d.  Set gridforcechecksize off in configuration file to ignore.\n", gridnum);
          NAMD_die(errmsg);      
        }
      }
	 }
	}
	
	Position center = grid->get_center();
	
	if (homePatch->flags.step % 100 == 1) {
	    DebugM(3, "center = " << center << "\n" << endi);
	    DebugM(3, "e = " << grid->get_e() << "\n" << endi);
	}
	
	if (grid->get_grid_type() == GridforceGrid::GridforceGridTypeFull) {
	    GridforceFullMainGrid *g = (GridforceFullMainGrid *)grid;
	    do_calc(g, gridnum, p, numAtoms, mol, forces, energy, extForce, extVirial);
	} else if (grid->get_grid_type() == GridforceGrid::GridforceGridTypeLite) {
	    GridforceLiteGrid *g = (GridforceLiteGrid *)grid;
	    do_calc(g, gridnum, p, numAtoms, mol, forces, energy, extForce, extVirial);
	}
    }
    reduction->item(REDUCTION_MISC_ENERGY) += energy;
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    reduction->submit();
}
/*			END OF FUNCTION force				*/
