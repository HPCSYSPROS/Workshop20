/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream>
#include <typeinfo>

#include "GridForceGrid.h"
#include "Vector.h"
#include "Node.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "InfoStream.h"
#include "common.h"
#include "ComputeGridForce.h"

#include "MGridforceParams.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

#include "GridForceGrid.inl"


/*****************/
/* GRIDFORCEGRID */
/*****************/

// Class methods

GridforceGrid * GridforceGrid::new_grid(int gridnum, char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams)
{
    GridforceGrid *grid = NULL;
    if (mgridParams->gridforceLite) {
	grid = new GridforceLiteGrid(gridnum);
    } else {
	grid = new GridforceFullMainGrid(gridnum);
    }
    
    grid->initialize(potfilename, simParams, mgridParams);
    
    return grid;
}

GridforceGrid::~GridforceGrid() { ; }

void GridforceGrid::pack_grid(GridforceGrid *grid, MOStream *msg)
{
    // Abstract interface for packing a grid into a message.  This
    // could easily be a non-static function as it was before, but is
    // static to make it similar to unpack_grid below, which does need
    // to be static since it creates a GridforceGrid object.
    msg->put(grid->get_grid_type());
    grid->pack(msg);
}

GridforceGrid * GridforceGrid::unpack_grid(int gridnum, MIStream *msg)
{
    // Abstract interface for unpacking a grid from a message.
    GridforceGrid *grid = NULL;
    int type;
    
    msg->get(type);
    
    switch (type) {
    case GridforceGridTypeFull:
	grid = new GridforceFullMainGrid(gridnum);
	break;
    case GridforceGridTypeLite:
	grid = new GridforceLiteGrid(gridnum);
	break;
    default:
	NAMD_bug("GridforceGrid::unpack_grid called with unknown grid type!");
    }
    
    grid->unpack(msg);
    
    return grid;
}

bool GridforceGrid::fits_lattice(const Lattice &lattice)
{
    // Idea: Go through each grid corner and wrap it to the grid center;
    // if the position moves, then the grid is too big and we return false
    DebugM(4, "Checking whether grid fits periodic cell\n");
    Position center = get_center();
    for (int i = 0; i < 8; i++) {
	Position pos = get_corner(i);
	Position pos_wrapped = wrap_position(pos, lattice);
	if ((pos - pos_wrapped).length() > 1.) {
	    DebugM(5, "(" << pos << ") != (" << pos_wrapped << ")\n" << endi);
	    return false;
	}
    }
    return true;
}

Position GridforceGrid::get_corner(int idx)
{
    // idx -> (x,y,z) (cell basis coordinates)
    // 0 -> (0,0,0)
    // 1 -> (1,0,0)
    // 2 -> (0,1,0)
    // 3 -> (1,1,0)
    // 4 -> (0,0,1)
    // 5 -> (1,0,1)
    // 6 -> (0,1,1)
    // 7 -> (1,1,1)
    Position pos;
    if (idx >= 8 || idx < 0) {
	// invalid index
	pos = Vector();	// default value of Vector() is strange enough to be a decent "undefined" value (-99999, -99999, -999999)
    } else if (corners[idx] != Vector()) {
	// use cached value if possible
	pos = corners[idx];
    } else {
	// must calculate
	Tensor e = get_e();
	pos = get_origin();
	if (idx & (1 << 0)) pos += e * Vector(get_k0()-1, 0, 0);
	if (idx & (1 << 1)) pos += e * Vector(0, get_k1()-1, 0);
	if (idx & (1 << 2)) pos += e * Vector(0, 0, get_k2()-1);
	corners[idx] = pos;	// cache for future use
	DebugM(4, "corner " << idx << " = " << pos << "\n" << endi);
    }
    return pos;
}


/*************************/
/* GRIDFORCEFULLBASEGRID */
/*************************/

GridforceFullBaseGrid::GridforceFullBaseGrid(void)
{
    cont[0] = cont[1] = cont[2] = FALSE;
    grid = NULL;
    numSubgrids = 0;
    subgrids = NULL;
}

GridforceFullBaseGrid::~GridforceFullBaseGrid()
{
    delete[] grid;
    for (int i = 0; i < numSubgrids; i++) {
	delete subgrids[i];
    }
    delete[] subgrids;
}


void GridforceFullBaseGrid::pack(MOStream *msg) const
{
    DebugM(2, "Packing message\n" << endi);
    
    msg->put(numSubgrids);
    msg->put(generation);
    
    msg->put(3*sizeof(int), (char*)k);
    msg->put(3*sizeof(int), (char*)k_nopad);
    msg->put(size);
    msg->put(size_nopad);
    msg->put(3*sizeof(long int), (char*)dk);
    msg->put(3*sizeof(long int), (char*)dk_nopad);
    msg->put(factor);
    
    msg->put(sizeof(Vector), (char*)&origin);
    msg->put(sizeof(Vector), (char*)&center);
    msg->put(sizeof(Tensor), (char*)&e);
    msg->put(sizeof(Tensor), (char*)&inv);
	     
//    msg->put(3*sizeof(float), (char*)pad_p);
//    msg->put(3*sizeof(float), (char*)pad_n);
    msg->put(3*sizeof(Bool), (char*)cont);
    msg->put(3*sizeof(float), (char*)offset);
    msg->put(3*sizeof(float), (char*)gap);
    msg->put(3*sizeof(float), (char*)gapinv);
    msg->put(sizeof(Vector), (char*)&scale);
    msg->put(sizeof(Bool), (char*)&checksize);
    
    DebugM(2, "Packing grid, size = " << size << "\n" << endi);
    
    msg->put(size*sizeof(float), (char*)grid);
    
    DebugM(2, "Packing subgrids\n" << endi);
    
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->pack(msg);
    }
}

void GridforceFullBaseGrid::unpack(MIStream *msg)
{
    DebugM(3, "Unpacking message\n" << endi);
//    iout << iINFO << CkMyPe() << " Unpacking message\n" << endi;

    delete[] grid;
    grid = NULL;
    for (int i = 0; i < numSubgrids; i++) {
	delete subgrids[i];
    }
    numSubgrids = 0;
    delete[] subgrids;
    subgrids = NULL;
    
    msg->get(numSubgrids);
    msg->get(generation);
    
    DebugM(3, "numSubgrids = " << numSubgrids << "\n");
    DebugM(3, "generation = " << generation << "\n" << endi);
    
    msg->get(3*sizeof(int), (char*)k);
    msg->get(3*sizeof(int), (char*)k_nopad);
    msg->get(size);
    msg->get(size_nopad);
    msg->get(3*sizeof(long int), (char*)dk);
    msg->get(3*sizeof(long int), (char*)dk_nopad);
    msg->get(factor);
    
    DebugM(3, "size = " << size << "\n" << endi);
    
    msg->get(sizeof(Vector), (char*)&origin);
    msg->get(sizeof(Vector), (char*)&center);
    msg->get(sizeof(Tensor), (char*)&e);
    msg->get(sizeof(Tensor), (char*)&inv);
	     
//    msg->get(3*sizeof(float), (char*)pad_p);
//    msg->get(3*sizeof(float), (char*)pad_n);
    msg->get(3*sizeof(Bool), (char*)cont);
    msg->get(3*sizeof(float), (char*)offset);
    msg->get(3*sizeof(float), (char*)gap);
    msg->get(3*sizeof(float), (char*)gapinv);
    msg->get(sizeof(Vector), (char*)&scale);
    msg->get(sizeof(Bool), (char*)&checksize);
    
    if (size) {
	DebugM(3, "allocating grid, size = " << size << "\n" << endi);
	grid = new float[size];
	msg->get(size*sizeof(float), (char*)grid);
    }
    
    if (numSubgrids) {
	DebugM(3, "Creating subgrids array, size " << numSubgrids << "\n" << endi);
	subgrids = new GridforceFullSubGrid *[numSubgrids];
	for (int i = 0; i < numSubgrids; i++) {
	    subgrids[i] = new GridforceFullSubGrid(this);
	    subgrids[i]->unpack(msg);
	}
    }
}


void GridforceFullBaseGrid::readHeader(SimParameters *simParams, MGridforceParams *mgridParams)
{
  char line[256];
  long int poten_offset;
  do {
      poten_offset = ftell(poten_fp);
      fgets(line, 256, poten_fp);	// Read comment lines
      DebugM(4, "Read line: " << line << endi);
  } while (line[0] == '#');
  fseek(poten_fp, poten_offset, SEEK_SET);
  
  // read grid dimensions
  fscanf(poten_fp, "object %*d class gridpositions counts %d %d %d\n",
         &k_nopad[0], &k_nopad[1], &k_nopad[2]);
  size_nopad = k_nopad[0] * k_nopad[1] * k_nopad[2];
  
  // Read origin
  fscanf(poten_fp, "origin %lf %lf %lf\n",
         &origin.x, &origin.y, &origin.z);
        
  // Read delta (unit vectors)
  // These are column vectors, so must fill gridfrcE tensor to reflect this
  fscanf(poten_fp, "delta %lf %lf %lf\n", &e.xx, &e.yx, &e.zx);
  fscanf(poten_fp, "delta %lf %lf %lf\n", &e.xy, &e.yy, &e.zy);
  fscanf(poten_fp, "delta %lf %lf %lf\n", &e.xz, &e.yz, &e.zz);
  
  center = origin + e * 0.5 
           * Position(k_nopad[0]-1, k_nopad[1]-1, k_nopad[2]-1);
  
  fscanf(poten_fp, "object %*d class gridconnections counts %*lf %*lf %*lf\n");
  fscanf(poten_fp, "object %*d class array type double rank 0 items %*d data follows\n");
    
  // Calculate inverse tensor
  BigReal det;
  det = e.xx*(e.yy*e.zz - e.yz*e.zy) - e.xy*(e.yx*e.zz - e.yz*e.zx) 
        + e.xz*(e.yx*e.zy - e.yy*e.zx);
  inv.xx =  (e.yy*e.zz - e.yz*e.zy)/det;
  inv.xy = -(e.xy*e.zz - e.xz*e.zy)/det;
  inv.xz =  (e.xy*e.yz - e.xz*e.yy)/det;
  inv.yx = -(e.yx*e.zz - e.yz*e.zx)/det;
  inv.yy =  (e.xx*e.zz - e.xz*e.zx)/det;
  inv.yz = -(e.xx*e.yz - e.xz*e.yx)/det;
  inv.zx =  (e.yx*e.zy - e.yy*e.zx)/det;
  inv.zy = -(e.xx*e.zy - e.xy*e.zx)/det;
  inv.zz =  (e.xx*e.yy - e.xy*e.yx)/det;
  
  DebugM(4, "origin = " << origin << "\n");
  DebugM(4, "e = " << e << "\n");
  DebugM(4, "inv = " << inv << "\n" << endi);
}

void GridforceFullBaseGrid::readSubgridHierarchy(FILE *poten, int &totalGrids)
{
    DebugM(4, "Beginning of readSubgridHierarchy, generation = " << generation << ", totalGrids = " << totalGrids << "\n" << endi);
    
    int elems, generation_in;
    
    subgrids = new GridforceFullSubGrid *[numSubgrids];
    
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i] = new GridforceFullSubGrid(this);
	elems = fscanf(poten_fp, "# namdnugrid subgrid %d generation %d min %d %d %d max %d %d %d subgrids count %d\n",
	       &subgrids[i]->subgridIdx, &generation_in,
	       &subgrids[i]->pmin[0], &subgrids[i]->pmin[1], &subgrids[i]->pmin[2],
	       &subgrids[i]->pmax[0], &subgrids[i]->pmax[1], &subgrids[i]->pmax[2],
	       &subgrids[i]->numSubgrids);
	if (elems < 9) {
	    char msg[256];
	    sprintf(msg, "Problem reading Gridforce potential file! (%d < 9)", elems);
	    NAMD_die(msg);
	}
	
	totalGrids++;
	
	if (subgrids[i]->subgridIdx != (totalGrids - 1)) {
	    char msg[256];
	    sprintf(msg, "Problem reading Gridforce potential file! (%d != %d)", subgrids[i]->subgridIdx, totalGrids - 1);
	    NAMD_die(msg);
	}
	if (subgrids[i]->generation != generation_in) {
	    char msg[256];
	    sprintf(msg, "Problem reading Gridforce potential file! (%d != %d)", subgrids[i]->generation, generation_in);
	    NAMD_die(msg);
	}
	
// 	DebugM(3, "setting maingrid\n");
// 	subgrids[i]->maingrid->subgrids_flat[subgrids[i]->subgridIdx] = subgrids[i];
// 	DebugM(3, "reading subgrid hierarchy\n");
	
	subgrids[i]->readSubgridHierarchy(poten, totalGrids);
    }
}


/*************************/
/* GRIDFORCEFULLMAINGRID */
/*************************/

GridforceFullMainGrid::GridforceFullMainGrid(int gridnum)
{
    mygridnum = gridnum;
    generation = 0;
    subgrids_flat = NULL;
    type = GridforceGridTypeFull;
}


GridforceFullMainGrid::~GridforceFullMainGrid()
{
    delete[] subgrids_flat;
}


void GridforceFullMainGrid::pack(MOStream *msg) const
{
    DebugM(4, "Packing maingrid\n" << endi);
    
//     msg->put(3*sizeof(float), (char*)pad_p);
//     msg->put(3*sizeof(float), (char*)pad_n);
    msg->put(totalGrids);
    msg->put(mygridnum);
    msg->put(129*sizeof(char), (char*)filename);
    
    DebugM(3, "calling GridforceFullBaseGrid::pack\n" << endi);
    
    GridforceFullBaseGrid::pack(msg);
}


void GridforceFullMainGrid::unpack(MIStream *msg)
{
    DebugM(4, "Unpacking maingrid\n" << endi);
    
//     msg->get(3*sizeof(float), (char*)pad_p);
//     msg->get(3*sizeof(float), (char*)pad_n);
    msg->get(totalGrids);
    msg->get(mygridnum);
    msg->get(129*sizeof(char), (char*)filename);
    
    GridforceFullBaseGrid::unpack(msg);
    
    DebugM(4, "size  = " << size << "\n");
    DebugM(4, "numSubgrids = " << numSubgrids << "\n");
    DebugM(4, "gapinv = " << gapinv[0] << " " << gapinv[2] << " " << gapinv[2] << " " << "\n");
    DebugM(4, "generation = " << generation << "\n" << endi);
    DebugM(4, "filename = " << filename << "\n" << endi);
    
    buildSubgridsFlat();
}


void GridforceFullMainGrid::buildSubgridsFlat(void)
{
    DebugM(4, "buildSubgridsFlat() called, totalGrids-1 = " << totalGrids-1 << "\n" << endi);
    delete[] subgrids_flat;
    subgrids_flat = new GridforceFullSubGrid *[totalGrids-1];
    for (int i = 0; i < numSubgrids; i++) {
	DebugM(3, "adding to subgridsFlat\n" << endi);
	subgrids[i]->addToSubgridsFlat();
	DebugM(3, "success!\n" << endi);
    }
    for (int i = 0; i < totalGrids-1; i++) {
	DebugM(4, "subgrids_flat[" << i << "]->numSubgrids = " << subgrids_flat[i]->numSubgrids << "\n" << endi);
    }
    for (int i = 0; i < numSubgrids; i++) {
	DebugM(4, "subgrids[" << i << "]->numSubgrids = " << subgrids[i]->numSubgrids << "\n" << endi);
    }
}


void GridforceFullMainGrid::initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams, int brd)
{
    if (brd >= 0) {
	border = brd;
    } else {
	border = default_border;
    }
    
    // FROM init1
    //FILE *poten = Fopen(potfilename, "r");
    poten_fp = Fopen(potfilename, "r");
    if (!poten_fp) {
	NAMD_die("Problem reading grid force potential file");
    }
    
    // save file name so that grid can be re-read via Tcl
    strcpy(filename, potfilename);
    
    // Read special comment fields and create subgrid objects
    totalGrids = 1;
    char line[256];
    Bool flag = FALSE;
    numSubgrids = 0;
    float version;
    long int poten_offset;
    do {
	poten_offset = ftell(poten_fp);
	fgets(line, 256, poten_fp);	// Read comment lines
	//flag = sscanf(line, "# maingrid subgrids count %d\n", &numSubgrids);
	flag = sscanf(line, "# namdnugrid version %f\n", &version);
    } while (line[0] == '#' && !flag);
    
    if (flag) {
	if (version != 1.0) {
	    NAMD_die("Unsupported version of non-uniform grid file format!");
	}
	fscanf(poten_fp, "# namdnugrid maingrid subgrids count %d\n", &numSubgrids);
	readSubgridHierarchy(poten_fp, totalGrids);
	buildSubgridsFlat();
    } else {
	fseek(poten_fp, poten_offset, SEEK_SET);
    }
    
    // Read header
    readHeader(simParams, mgridParams);
    
    factor = 1.0;
    if (mgridParams->gridforceVolts)
    {
	factor /= 0.0434;  // convert V -> kcal/mol*e
    }
    scale = mgridParams->gridforceScale;
    checksize = mgridParams->gridforceCheckSize;

    // Allocate storage for potential and read it
    float *grid_nopad = new float[size_nopad];
    
    float tmp2;
    for (long int count = 0; count < size_nopad; count++) {
	int err = fscanf(poten_fp, "%f", &tmp2);
	if (err == EOF || err == 0) {
	    NAMD_die("Grid force potential file incorrectly formatted");
	}
	grid_nopad[count] = tmp2 * factor;	// temporary, so just store flat
    }
    fscanf(poten_fp, "\n");
    
    // Shortcuts for accessing 1-D array with four indices
    dk_nopad[0] = k_nopad[1] * k_nopad[2];
    dk_nopad[1] = k_nopad[2];
    dk_nopad[2] = 1;
    
    Vector Kvec[3];
    Kvec[0] = e * Position(k_nopad[0]-1, 0, 0);
    Kvec[1] = e * Position(0, k_nopad[1]-1, 0);
    Kvec[2] = e * Position(0, 0, k_nopad[2]-1);
    Vector Avec[3];
    Avec[0] = simParams->lattice.a();
    Avec[1] = simParams->lattice.b();
    Avec[2] = simParams->lattice.c();
    
    // Decide whether we're wrapping
    for (int i0 = 0; i0 < 3; i0++) {
	if (mgridParams->gridforceCont[i0])
	{
	    Bool found = FALSE;
	    for (int i1 = 0; i1 < 3; i1++) {
		if (cross(Avec[i0].unit(), Kvec[i1].unit()).length() < 1e-4) {
		    found = TRUE;
		    cont[i1] = TRUE;
		    offset[i1] = mgridParams->gridforceVOffset[i0] * factor;
		    // want in grid-point units (normal = 1)
		    gap[i1] = (inv * (Avec[i0] - Kvec[i1])).length();	
		    gapinv[i1] = 1.0/gap[i1];
		    
		    if (gap[i1] < 0) {
			NAMD_die("Gridforce Grid overlap!");
		    }
		    
		    DebugM(4, "cont[" << i1 << "] = " << cont[i1] << "\n");
		    DebugM(4, "gap[" << i1 << "] = " << gap[i1] << "\n");
		    DebugM(4, "gapinv[" << i1 << "] = " << gapinv[i1] << "\n" << endi);
		}
	    }
	    
	    if (!found) {
		NAMD_die("No Gridforce unit vector found parallel to requested continuous grid direction!");
	    }
	} else {
	    // check for grid overlap in non-wrapping dimensions
	    // occurs below
	}
    }
    
    // Figure out size of true grid (padded on non-periodic sides)
    Vector delta = 0;
    for (int i = 0; i < 3; i++) {
	if (cont[i]) {
	    k[i] = k_nopad[i];
	} else {
	    k[i] = k_nopad[i] + 2*border;
	    delta[i] -= border;
	}
    }
    DebugM(4, "delta = " << e * delta << " (" << delta << ")\n" << endi);
    origin += e * delta;
    
    // Check for grid overlap
    if (!fits_lattice(simParams->lattice)) {
      char errmsg[512];
      if (checksize) {
        sprintf(errmsg, "Warning: Periodic cell basis too small for Gridforce grid %d.  Set gridforcechecksize off in configuration file to ignore.\n", mygridnum);
        NAMD_die(errmsg);      
      }
    }
    
    size = k[0] * k[1] * k[2];
    dk[0] = k[1] * k[2];
    dk[1] = k[2];
    dk[2] = 1;
    
    DebugM(3, "size = " << size << ", size_nopad = " << size_nopad << "\n" << endi);
    
    delete[] grid;
    grid = new float[size];
    
    n_sum[0] = n_sum[1] = n_sum[2] = 0;
    p_sum[0] = p_sum[1] = p_sum[2] = 0;
    for (int i0 = 0; i0 < k_nopad[0]; i0++) {
	for (int i1 = 0; i1 < k_nopad[1]; i1++) {
	    for (int i2 = 0; i2 < k_nopad[2]; i2++) {
		// Edges are special cases -- take force there to be
		// zero for smooth transition across potential
		// boundary
		
		long int ind_nopad = i0*dk_nopad[0] + i1*dk_nopad[1] + i2*dk_nopad[2];
		int j0 = (cont[0]) ? i0 : i0 + border;
		int j1 = (cont[1]) ? i1 : i1 + border;
		int j2 = (cont[2]) ? i2 : i2 + border;
		long int ind = j0*dk[0] + j1*dk[1] + j2*dk[2];
		
		if (i0 == 0)			n_sum[0] += grid_nopad[ind_nopad];
		else if (i0 == k_nopad[0]-1)	p_sum[0] += grid_nopad[ind_nopad];
		if (i1 == 0)			n_sum[1] += grid_nopad[ind_nopad];
		else if (i1 == k_nopad[1]-1)	p_sum[1] += grid_nopad[ind_nopad];
		if (i2 == 0)			n_sum[2] += grid_nopad[ind_nopad];
		else if (i2 == k_nopad[2]-1)	p_sum[2] += grid_nopad[ind_nopad];
		
		//grid[ind] = grid_nopad[ind_nopad];
		set_grid(j0, j1, j2, grid_nopad[ind_nopad]);
	    }
	}
    }
    
    const BigReal modThresh = 1.0;
    
    BigReal n_avg[3], p_avg[3];
    int i0;
    for (int i0 = 0; i0 < 3; i0++) {
	int i1 = (i0 + 1) % 3;
	int i2 = (i0 + 2) % 3;
	n_avg[i0] = n_sum[i0] / (k_nopad[i1] * k_nopad[i2]);
	p_avg[i0] = p_sum[i0] / (k_nopad[i1] * k_nopad[i2]);
	
	if (cont[i0] && fabs(offset[i0] - (p_avg[i0]-n_avg[i0])) > modThresh) 
	{
	    iout << iWARN << "GRID FORCE POTENTIAL DIFFERENCE IN K" << i0
		 << " DIRECTION IS " 
		 << offset[i0] - (p_avg[i0]-n_avg[i0]) 
		 << " KCAL/MOL*E\n" << endi;
	}
    }
    
    Bool twoPadVals = (cont[0] + cont[1] + cont[2] == 2);
    double padVal = 0.0;
    long int weight = 0;
    if (!twoPadVals) {
	// Determine pad value (must average)
	if (!cont[0]) {
	    padVal += p_sum[0] + n_sum[0];
	    weight += 2 * k_nopad[1] * k_nopad[2];
	}
	if (!cont[1]) {
	    padVal += p_sum[1] + n_sum[1];
	    weight += 2 * k_nopad[0] * k_nopad[2];
	}
	if (!cont[2]) {
	    padVal += p_sum[2] + n_sum[2];
	    weight += 2 * k_nopad[0] * k_nopad[1];
	}
	padVal /= weight;
    }
    
    for (int i = 0; i < 3; i++) {
	pad_n[i] = (cont[i]) ? 0.0 : (twoPadVals) ? n_avg[i] : padVal;
	pad_p[i] = (cont[i]) ? 0.0 : (twoPadVals) ? p_avg[i] : padVal;
	DebugM(4, "pad_n[" << i << "] = " << pad_n[i] << "\n");
	DebugM(4, "pad_p[" << i << "] = " << pad_p[i] << "\n" << endi);
    }
    
    if (cont[0] && cont[1] && cont[2]) {
	// Nothing to do
	return;
    }
    
    // Now fill in rest of new grid
    for (int i0 = 0; i0 < k[0]; i0++) {
	for (int i1 = 0; i1 < k[1]; i1++) {
	    for (int i2 = 0; i2 < k[2]; i2++) {
		if ( (cont[0] || (i0 >= border && i0 < k[0]-border)) 
		     && (cont[1] || (i1 >= border && i1 < k[1]-border)) 
		     && (cont[2] || i2 == border) )
		{
		    i2 += k_nopad[2]-1;
		    continue;
		}
		
		long int ind = i0*dk[0] + i1*dk[1] + i2*dk[2];

		Position pos = e * Position(i0, i1, i2);
		int var[3] = {i0, i1, i2};
		
		for (int dir = 0; dir < 3; dir++) {
		    if (cont[dir]) 
			continue;
  		    
		    if (var[dir] < border) {
			//grid[ind] = pad_n[dir];
			set_grid(i0, i1, i2, pad_n[dir]);
		    } else if (var[dir] >= k[dir]-border) {
			//grid[ind] = pad_p[dir];
			set_grid(i0, i1, i2, pad_p[dir]);
		    }
		}
		
// 		DebugM(2, "grid[" << ind << "; " << i0 << ", " << i1
// 		       << ", " << i2 << "] = " << get_grid(ind)
// 		       << "\n" << endi);
	    }
	}
    }
    
    for (int i0 = 0; i0 < k[0]; i0++) {
	for (int i1 = 0; i1 < k[1]; i1++) {
	    for (int i2 = 0; i2 < k[2]; i2++) {
		DebugM(1, "grid[" << i0 << ", " << i1 << ", " << i2 << "] = " << get_grid(i0,i1,i2) << "\n" << endi);
	    }
	}
    }
    
    // Clean up
    DebugM(3, "clean up\n" << endi);
    delete[] grid_nopad;
    
    // Call initialize for each subgrid
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->poten_fp = poten_fp;
	subgrids[i]->initialize(simParams, mgridParams);
    }
    
    // close file pointer
    fclose(poten_fp);
}


void GridforceFullMainGrid::reinitialize(SimParameters *simParams, MGridforceParams *mgridParams)
{
    DebugM(4, "reinitializing grid\n" << endi);
    initialize(filename, simParams, mgridParams);
}


long int GridforceFullMainGrid::get_all_gridvals(float **all_gridvals) const
{
    // Creates a flat array of all grid values, including subgrids,
    // and puts it in the value pointed to by the 'grids'
    // argument. Returns the resulting array size. Caller is
    // responsible for destroying the array via 'delete[]'
    
    DebugM(4, "get_all_gridvals called\n" << endi);
    
    long int sz = 0;
    sz += size;
    for (int i = 0; i < totalGrids-1; i++) {
	sz += subgrids_flat[i]->size;
    }
    DebugM(4, "size = " << sz << "\n" << endi);
    
    float *grid_vals = new float[sz];
    long int idx = 0;
    for (long int i = 0; i < size; i++) {
	grid_vals[idx++] = grid[i];
    }
    for (int j = 0; j < totalGrids-1; j++) {
	for (long int i = 0; i < subgrids_flat[j]->size; i++) {
	    grid_vals[idx++] = subgrids_flat[j]->grid[i];
	}
    }
    CmiAssert(idx == sz);
    
    *all_gridvals = grid_vals;
    
    DebugM(4, "get_all_gridvals finished\n" << endi);
    
    return sz;
}


void GridforceFullMainGrid::set_all_gridvals(float *all_gridvals, long int sz)
{
    DebugM(4, "set_all_gridvals called\n" << endi);
    
    long int sz_calc = 0;
    sz_calc += size;
    for (int i = 0; i < totalGrids-1; i++) {
	sz_calc += subgrids_flat[i]->size;
    }
    CmiAssert(sz == sz_calc);
    
    long int idx = 0;
    for (long int i = 0; i < size; i++) {
	DebugM(1, "all_gridvals[" << idx << "] = " << all_gridvals[idx] << "\n" << endi);
	grid[i] = all_gridvals[idx++];
    }
    for (int j = 0; j < totalGrids-1; j++) {
	for (long int i = 0; i < subgrids_flat[j]->size; i++) {
	    DebugM(1, "all_gridvals[" << idx << "] = " << all_gridvals[idx] << "\n" << endi);
	    subgrids_flat[j]->grid[i] = all_gridvals[idx++];
	}
    }
    CmiAssert(idx == sz);

    DebugM(4, "set_all_gridvals finished\n" << endi);
}


void GridforceFullMainGrid::compute_b(float *b, int *inds, Vector gapscale) const
{
    for (int i0 = 0; i0 < 8; i0++) {
	int inds2[3];
	int zero_derivs = FALSE;
	
	float voff = 0.0;
	int bit = 1;	// bit = 2^i1 in the below loop
	for (int i1 = 0; i1 < 3; i1++) {
	    inds2[i1] = (inds[i1] + ((i0 & bit) ? 1 : 0)) % k[i1];
	    
	    // Deal with voltage offsets
	    if (cont[i1] && inds[i1] == (k[i1]-1) && inds2[i1] == 0) {
		voff += offset[i1];
		DebugM(3, "offset[" << i1 << "] = " << offset[i1] << "\n" << endi);
	    }
	    
	    bit <<= 1;	// i.e. multiply by 2
	}
	
 	DebugM(1, "inds2 = " << inds2[0] << " " << inds2[1] << " " << inds2[2] << "\n" << endi);
	
	// NOTE: leaving everything in terms of unit cell coordinates for now,
	// eventually will multiply by inv tensor when applying the force
	
	// First set variables 'dk_{hi,lo}' (glob notation). The 'hi'
	// ('lo') variable in a given dimension is the number added (subtracted)
	// to go up (down) one grid point in that dimension; both are normally
	// just the corresponding 'dk[i]'. However, if we are sitting on a
	// boundary and we are using a continuous grid, then we want to map the
	// next point off the grid back around to the other side. e.g. point
	// (k[0], i1, k) maps to point (0, i1, k), which would be
	// accomplished by changing 'dk1_hi' to -(k[0]-1)*dk1.
	
	int d_hi[3] = {1, 1, 1};
	int d_lo[3] = {1, 1, 1};
	float voffs[3];
	float dscales[3] = {0.5, 0.5, 0.5};
	for (int i1 = 0; i1 < 3; i1++) {
	    if (inds2[i1] == 0) {
		if (cont[i1]) {
		    d_lo[i1] = -(k[i1]-1);
		    voffs[i1] = offset[i1];
		    dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
		}
		else zero_derivs = TRUE;
	    }
	    else if (inds2[i1] == k[i1]-1) {
		if (cont[i1]) {
		    d_hi[i1] = -(k[i1]-1);
		    voffs[i1] = offset[i1];
		    dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
		}
		else zero_derivs = TRUE;
	    }
	    else {
		voffs[i1] = 0.0;
	    }
	}
	
// 	DebugM(2, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << "\n" << endi);
// 	DebugM(2, "zero_derivs = " << zero_derivs << "\n" << endi);
// 	DebugM(2, "d_hi = " << d_hi[0] << " " << d_hi[1] << " " << d_hi[2] << "\n" << endi);
// 	DebugM(2, "d_lo = " << d_lo[0] << " " << d_lo[1] << " " << d_lo[2] << "\n" << endi);
 	DebugM(1, "dscales = " << dscales[0] << " " << dscales[1] << " " << dscales[2] << "\n" << endi);
 	DebugM(1, "voffs = " << voffs[0] << " " << voffs[1] << " " << voffs[2] << "\n" << endi);
	
	// V
	b[i0] = get_grid(inds2[0],inds2[1],inds2[2]) + voff;
	
	if (zero_derivs) {
	    DebugM(2, "zero_derivs\n" << endi);
	    b[8+i0] = 0.0;
	    b[16+i0] = 0.0;
	    b[24+i0] = 0.0;
	    b[32+i0] = 0.0;
	    b[40+i0] = 0.0;
	    b[48+i0] = 0.0;
	    b[56+i0] = 0.0;
	} else {
	    b[8+i0]  = dscales[0] * (get_grid_d(inds2[0]+d_hi[0],inds2[1],inds2[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1],inds2[2]) + voffs[0]);	//  dV/dx
	    b[16+i0] = dscales[1] * (get_grid_d(inds2[0],inds2[1]+d_hi[1],inds2[2]) - get_grid_d(inds2[0],inds2[1]-d_lo[1],inds2[2]) + voffs[1]);	//  dV/dy
	    b[24+i0] = dscales[2] * (get_grid_d(inds2[0],inds2[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0],inds2[1],inds2[2]-d_lo[2]) + voffs[2]);	//  dV/dz
	    b[32+i0] = dscales[0] * dscales[1]
		* (get_grid_d(inds2[0]+d_hi[0],inds2[1]+d_hi[1],inds2[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1]+d_hi[1],inds2[2]) -
		   get_grid_d(inds2[0]+d_hi[0],inds2[1]-d_lo[1],inds2[2]) + get_grid_d(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]));	//  d2V/dxdy
	    b[40+i0] = dscales[0] * dscales[2]
		* (get_grid_d(inds2[0]+d_hi[0],inds2[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1],inds2[2]+d_hi[2]) -
		   get_grid_d(inds2[0]+d_hi[0],inds2[1],inds2[2]-d_lo[2]) + get_grid_d(inds2[0]-d_lo[0],inds2[1],inds2[2]-d_lo[2]));	//  d2V/dxdz
	    b[48+i0] = dscales[1] * dscales[2]
		* (get_grid_d(inds2[0],inds2[1]+d_hi[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0],inds2[1]-d_lo[1],inds2[2]+d_hi[2]) -
		   get_grid_d(inds2[0],inds2[1]+d_hi[1],inds2[2]-d_lo[2]) + get_grid_d(inds2[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));	//  d2V/dydz
	
	    b[56+i0] = dscales[0] * dscales[1] * dscales[2]					// d3V/dxdydz
		* (get_grid_d(inds2[0]+d_hi[0],inds2[1]+d_hi[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]+d_hi[0],inds2[1]+d_hi[1],inds2[2]-d_lo[2]) -
		   get_grid_d(inds2[0]+d_hi[0],inds2[1]-d_lo[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1]+d_hi[1],inds2[2]+d_hi[2]) +
		   get_grid_d(inds2[0]+d_hi[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]) + get_grid_d(inds2[0]-d_lo[0],inds2[1]+d_hi[1],inds2[2]-d_lo[2]) +
		   get_grid_d(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));
	}
	
	DebugM(1, "V = " << b[i0] << "\n");
	
	DebugM(1, "dV/dx = " << b[8+i0] << "\n");
	DebugM(1, "dV/dy = " << b[16+i0] << "\n");
	DebugM(1, "dV/dz = " << b[24+i0] << "\n");
	
	DebugM(1, "d2V/dxdy = " << b[32+i0] << "\n");
	DebugM(1, "d2V/dxdz = " << b[40+i0] << "\n");
	DebugM(1, "d2V/dydz = " << b[48+i0] << "\n");
	
	DebugM(1, "d3V/dxdydz = " << b[56+i0] << "\n" << endi);
    }
}


/************************/
/* GRIDFORCEFULLSUBGRID */
/************************/

GridforceFullSubGrid::GridforceFullSubGrid(GridforceFullBaseGrid *parent_in) {
    parent = parent_in;
    generation = parent->generation + 1;
    GridforceFullBaseGrid *tmp = parent;
    while (tmp->generation > 0) {
	tmp = ((GridforceFullSubGrid *)tmp)->parent;
    }
    maingrid = (GridforceFullMainGrid *)tmp;
    DebugM(4, "generation = " << generation << "\n" << endi);
}


void GridforceFullSubGrid::initialize(SimParameters *simParams, MGridforceParams *mgridParams)
{
    int tmp;
    char line[256];
    long int poten_offset;
    
    // Skip 'attribute's
    DebugM(3, "Skipping 'attribute' keywords...\n" << endi);
    char str[256];
    do {
	poten_offset = ftell(poten_fp);
	fscanf(poten_fp, "%s", str);
	fgets(line, 256, poten_fp);
	DebugM(4, "Read line " << str << " " << line << endi);
    } while (strcmp(str, "attribute") == 0);
    fseek(poten_fp, poten_offset, SEEK_SET);
    
    // Skip 'field' object
    DebugM(3, "Skipping 'field' object\n" << endi);
    fscanf(poten_fp, "object");
    int n;
    n = fscanf(poten_fp, "\"%[^\"]\" class field\n", str);
    if (n == 0) {
	n = fscanf(poten_fp, "%d class field\n", &tmp);
    }
    
    if (n == 0) {
	NAMD_die("Error reading gridforce grid! Could not find field object!\n");
    }
    
    // Skip 'component's
    DebugM(3, "Skipping 'component' keywords\n" << endi);
    do {
	poten_offset = ftell(poten_fp);
	fscanf(poten_fp, "%s", str);
	fgets(line, 256, poten_fp);
    } while (strcmp(str, "component") == 0);
    fseek(poten_fp, poten_offset, SEEK_SET);
    
    // Read header
    readHeader(simParams, mgridParams);
    
    factor = 1.0;
    if (mgridParams->gridforceVolts)
    {
	factor /= 0.0434;  // convert V -> kcal/mol*e
    }
    scale = mgridParams->gridforceScale;
    checksize = mgridParams->gridforceCheckSize;
    
    for (int i = 0; i < 3; i++) {
	k[i] = k_nopad[i];	// subgrids aren't padded
    }
    
    // Make sure that each subgrid dimension is an integral
    // number of spanned supergrid cells. This is to ensure that no
    // supergrid nodes land in the middle of a subgrid, because in
    // this case forces could not be matched properly.
    for (int i = 0; i < 3; i++) {
	if ((k[i] - 1) % (pmax[i] - pmin[i] + 1) != 0) {
	    iout << (k[i] - 1) << " % " << (pmax[i] - pmin[i] + 1) << " != 0\n" << endi;
	    NAMD_die("Error reading gridforce grid! Subgrid dimensions must be an integral number spanned parent cells!");
	}
    }
    
    for (int i = 0; i < 3; i++) {
	if (parent->cont[i]) {
	    cont[i] = (pmin[i] == 0 && pmax[i] == parent->k[i]-2) ? TRUE : FALSE;
	    DebugM(3, "pmin[" << i << "] = " << pmin[i] << " pmax[" << i << "] = " << pmax[i] << " parent->k[" << i << "] = " << parent->k[i] << " cont[" << i << "] = " << cont[i] << "\n" << endi);
	} else {
	    cont[i] = false;
	    if (parent->generation == 0) {
		// Add to pmin, pmax since parent got extra gridpoint layer(s) (maybe)
		int brd = parent->get_border();
		pmin[i] += brd;
		pmax[i] += brd;
	    }
	}		
    }
    
    DebugM(4, "pmin = " << pmin[0] << " " << pmin[1] << " " << pmin[2] << "\n");
    DebugM(4, "pmax = " << pmax[0] << " " << pmax[1] << " " << pmax[2] << "\n" << endi);
    
    Vector origin2 = parent->origin + parent->e * Position(pmin[0], pmin[1], pmin[2]);
    Vector escale, invscale;
    for (int i = 0; i < 3; i++) {
	escale[i] = double(pmax[i] - pmin[i] + 1)/(k[i]-1);
	invscale[i] = 1.0/escale[i];
	if (cont[i]) { pmax[i]++; }
    }
    Tensor e2 = tensorMult(parent->e, Tensor::diagonal(escale));
    
    // Check that lattice parameters agree with min and max numbers
    // from subgrid hierarchy.
    double TOL2 = 1e-4;	// Totally arbitrary
    if (pow(origin2.x-origin.x, 2) > TOL2 ||
	pow(origin2.y-origin.y, 2) > TOL2 ||
	pow(origin2.z-origin.z, 2) > TOL2 ||
	pow(e2.xx-e.xx, 2) > TOL2 ||
	pow(e2.xy-e.xy, 2) > TOL2 ||
	pow(e2.xz-e.xz, 2) > TOL2 ||
	pow(e2.yx-e.yx, 2) > TOL2 ||
	pow(e2.yy-e.yy, 2) > TOL2 ||
	pow(e2.yz-e.yz, 2) > TOL2 ||
	pow(e2.zx-e.zx, 2) > TOL2 ||
	pow(e2.zy-e.zy, 2) > TOL2 ||
	pow(e2.zz-e.zz, 2) > TOL2)
    {
	NAMD_die("Error reading gridforce grid! Subgrid lattice does not match!");
    }
    
    // Overwrite what was read from the header
    origin = origin2;
    e = e2;
    
    inv = tensorMult(Tensor::diagonal(invscale), parent->inv);
    for (int i = 0; i < 3; i++) {
	gap[i] = escale[i] * parent->gap[i];
	gapinv[i] = invscale[i] * parent->gapinv[i];
	offset[i] = parent->offset[i];
    }
    center = origin + e * 0.5 * Position(k[0], k[1], k[2]);
    
    DebugM(4, "origin = " << origin << "\n");
    DebugM(4, "e = " << e << "\n");
    DebugM(4, "inv = " << inv << "\n");
    DebugM(4, "gap = " << gap[0] << " " << gap[2] << " " << gap[2] << " " << "\n");
    DebugM(4, "gapinv = " << gapinv[0] << " " << gapinv[2] << " " << gapinv[2] << " " << "\n");
    DebugM(4, "numSubgrids = " << numSubgrids << "\n");
    DebugM(4, "k = " << k[0] << " " << k[1] << " " << k[2] << "\n");
    DebugM(4, "escale = " << escale << "\n");
    DebugM(4, "invscale = " << invscale << "\n" << endi);
    
    /*** Set members ***/
    size = k[0] * k[1] * k[2];
    dk[0] = k[1] * k[2];
    dk[1] = k[2];
    dk[2] = 1;
    
    scale_dV = Tensor::diagonal(escale);
    scale_d2V = Tensor::diagonal(Vector(escale.x*escale.y, escale.x*escale.z, escale.y*escale.z));
    scale_d3V = escale.x * escale.y * escale.z;
    
    DebugM(4, "scale_dV = " << scale_dV << "\n");
    DebugM(4, "scale_d2V = " << scale_d2V << "\n");
    DebugM(4, "scale_d3V = " << scale_d3V << "\n" << endi);
    
    // Allocate storage for potential and read it
    float *grid_tmp = new float[size];
    
    float tmp2;
    DebugM(3, "size_nopad = " << size_nopad << "\n");
    for (long int count = 0; count < size_nopad; count++) {
// 	poten_offset = ftell(poten_fp);
// 	fscanf(poten_fp, "%s", str);
// 	fgets(line, 256, poten_fp);
// 	DebugM(4, "Read line " << str << " " << line << endi);
// 	fseek(poten_fp, poten_offset, SEEK_SET);
	
	int err = fscanf(poten_fp, "%f", &tmp2);
	if (err == EOF || err == 0) {
	    NAMD_die("Grid force potential file incorrectly formatted");
	}
	grid_tmp[count] = tmp2 * factor;
    }
    fscanf(poten_fp, "\n");
    
    // Set real grid
    DebugM(3, "allocating grid\n" << endi);
    delete[] grid;
    grid = new float[size];
    for (int i0 = 0; i0 < k_nopad[0]; i0++) {
	for (int i1 = 0; i1 < k_nopad[1]; i1++) {
	    for (int i2 = 0; i2 < k_nopad[2]; i2++) {
		long int ind = i0*dk[0] + i1*dk[1] + i2*dk[2];
		set_grid(i0, i1, i2, grid_tmp[ind]);
	    }
	}
    }
    
    for (int i0 = 0; i0 < k[0]; i0++) {
	for (int i1 = 0; i1 < k[1]; i1++) {
	    for (int i2 = 0; i2 < k[2]; i2++) {
		DebugM(1, "grid[" << i0 << ", " << i1 << ", " << i2 << "] = " << get_grid(i0,i1,i2) << "\n" << endi);
	    }
	}
    }
    
    // Clean up
    delete[] grid_tmp;
    
    // Call initialize for each subgrid
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->initialize(simParams, mgridParams);
    }
}


void GridforceFullSubGrid::pack(MOStream *msg) const
{
    DebugM(4, "Packing subgrid\n" << endi);
    
    msg->put(sizeof(Tensor), (char*)&scale_dV);
    msg->put(sizeof(Tensor), (char*)&scale_d2V);
    msg->put(sizeof(float), (char*)&scale_d3V);
    
    msg->put(3*sizeof(int), (char*)pmin);
    msg->put(3*sizeof(int), (char*)pmax);
    msg->put(subgridIdx);
    
    DebugM(3, "calling GridforceFullBaseGrid::pack\n" << endi);
    
    GridforceFullBaseGrid::pack(msg);
}


void GridforceFullSubGrid::unpack(MIStream *msg)
{
    DebugM(4, "Unpacking subgrid\n" << endi);
    
    msg->get(sizeof(Tensor), (char*)&scale_dV);
    msg->get(sizeof(Tensor), (char*)&scale_d2V);
    msg->get(sizeof(float), (char*)&scale_d3V);
    
    msg->get(3*sizeof(int), (char*)pmin);
    msg->get(3*sizeof(int), (char*)pmax);
    msg->get(subgridIdx);
    
    GridforceFullBaseGrid::unpack(msg);
    
    DebugM(4, "size  = " << size << "\n");
    DebugM(4, "numSubgrids = " << numSubgrids << "\n");
    DebugM(4, "gapinv = " << gapinv[0] << " " << gapinv[2] << " " << gapinv[2] << " " << "\n");
    DebugM(4, "generation = " << generation << "\n" << endi);
}


void GridforceFullSubGrid::addToSubgridsFlat(void)
{
    DebugM(4, "addToSubgridsFlat() called, subgridIdx = " << subgridIdx << ", maingrid->numSubgrids = " << maingrid->numSubgrids << "\n" << endi);
    maingrid->subgrids_flat[subgridIdx-1] = this;
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->addToSubgridsFlat();
    }
}


void GridforceFullSubGrid::compute_b(float *b, int *inds, Vector gapscale) const
{
    for (int i0 = 0; i0 < 8; i0++) {
	int inds2[3];
	
	float voff = 0.0;
	int bit = 1;	// bit = 2^i1 in the below loop
	for (int i1 = 0; i1 < 3; i1++) {
	    inds2[i1] = (inds[i1] + ((i0 & bit) ? 1 : 0)) % k[i1];
	    
	    // Deal with voltage offsets
	    if (cont[i1] && inds[i1] == (k[i1]-1) && inds2[i1] == 0) {
		voff += offset[i1];
		DebugM(3, "offset[" << i1 << "] = " << offset[i1] << "\n" << endi);
	    }
	    
	    bit <<= 1;	// i.e. multiply by 2
	}
	
 	DebugM(3, "inds2 = " << inds2[0] << " " << inds2[1] << " " << inds2[2] << "\n" << endi);
	
	int d_hi[3] = {1, 1, 1}; 
	int d_lo[3] = {1, 1, 1};
	float voffs[3];
	float dscales[3] = {0.5, 0.5, 0.5};
	for (int i1 = 0; i1 < 3; i1++) {
	    if (inds2[i1] == 0 && cont[i1]) {
		d_lo[i1] = -(k[i1]-1);
		voffs[i1] = offset[i1];
		dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
	    }
	    else if (inds2[i1] == k[i1]-1 && cont[i1]) {
		d_hi[i1] = -(k[i1]-1);
		voffs[i1] = offset[i1];
		dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
	    }
	    else {
		voffs[i1] = 0.0;
	    }
	}
	
	bool edge = false;
	for (int i1 = 0; i1 < 3; i1++) {
	    if (!cont[i1] && (inds2[i1] == 0 || inds2[i1] == k[i1]-1)) {
		edge = true;
	    }
	}
	
	if (inds2[2] == 0) {
// 	    DebugM(3, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << " d_hi = " << d_hi[0] << " " << d_hi[1] << " " << d_hi[2] << " d_lo = " << d_lo[0] << " " << d_lo[1] << " " << d_lo[2] << " dscales = " << dscales[0] << " " << dscales[1] << " " << dscales[2] << "\n" << endi);
	    DebugM(3, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << "\n" << endi);
	}
	
	if (edge) {
	    DebugM(2, "Edge!\n" << endi);
	    
	    // Must get derivatives from parent
	    Position pos = e * Vector(inds2[0], inds2[1], inds2[2]) + origin;	// Gridpoint position in realspace
	    Vector g = parent->inv * (pos - parent->origin);	// Gridpoint position in parent's gridspace
	    Vector dg;
	    int inds3[3];
	    
	    DebugM(2, "g = " << g << "\n" << endi);
	    
	    for (int i = 0; i < 3; i++) {
		inds3[i] = (int)floor(g[i]);
		dg[i] = g[i] - inds3[i];
	    }
    
	    float x[4], y[4], z[4];
	    x[0] = 1; y[0] = 1; z[0] = 1;
	    for (int j = 1; j < 4; j++) {
		x[j] = x[j-1] * dg.x;
		y[j] = y[j-1] * dg.y;
		z[j] = z[j-1] * dg.z;
		DebugM(1, "x[" << j << "] = " << x[j] << "\n");
		DebugM(1, "y[" << j << "] = " << y[j] << "\n");
		DebugM(1, "z[" << j << "] = " << z[j] << "\n" << endi);
	    }
	    
	    // Compute parent matrices
	    float b_parent[64];
	    parent->compute_b(b_parent, inds3, gapscale);
	    
	    float a_parent[64];
	    parent->compute_a(a_parent, b_parent);
	    
	    // Compute parent derivatives
	    float V = parent->compute_V(a_parent, x, y, z);
	    Vector dV = scale_dV * parent->compute_dV(a_parent, x, y, z);
	    Vector d2V = scale_d2V * parent->compute_d2V(a_parent, x, y, z);
	    float d3V = scale_d3V * parent->compute_d3V(a_parent, x, y, z);
	    
	    b[i0] = V;
	    b[8+i0] = dV[0];
	    b[16+i0] = dV[1];
	    b[24+i0] = dV[2];
	    b[32+i0] = d2V[0];
	    b[40+i0] = d2V[1];
	    b[48+i0] = d2V[2];
	    b[56+i0] = d3V;
	} else {
	    b[i0] = get_grid(inds2[0],inds2[1],inds2[2]) + voff;	// V
	    
	    b[8+i0]  = dscales[0] * (get_grid_d(inds2[0]+d_hi[0],inds2[1],inds2[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1],inds2[2]) + voffs[0]);	//  dV/dx
	    b[16+i0] = dscales[1] * (get_grid_d(inds2[0],inds2[1]+d_hi[1],inds2[2]) - get_grid_d(inds2[0],inds2[1]-d_lo[1],inds2[2]) + voffs[1]);	//  dV/dy
	    b[24+i0] = dscales[2] * (get_grid_d(inds2[0],inds2[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0],inds2[1],inds2[2]-d_lo[2]) + voffs[2]);	//  dV/dz
	    b[32+i0] = dscales[0] * dscales[1]
		* (get_grid_d(inds2[0]+d_hi[0],inds2[1]+d_hi[1],inds2[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1]+d_hi[1],inds2[2]) -
		   get_grid_d(inds2[0]+d_hi[0],inds2[1]-d_lo[1],inds2[2]) + get_grid_d(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]));	//  d2V/dxdy
	    b[40+i0] = dscales[0] * dscales[2]
		* (get_grid_d(inds2[0]+d_hi[0],inds2[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1],inds2[2]+d_hi[2]) -
		   get_grid_d(inds2[0]+d_hi[0],inds2[1],inds2[2]-d_lo[2]) + get_grid_d(inds2[0]-d_lo[0],inds2[1],inds2[2]-d_lo[2]));	//  d2V/dxdz
	    b[48+i0] = dscales[1] * dscales[2]
		* (get_grid_d(inds2[0],inds2[1]+d_hi[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0],inds2[1]-d_lo[1],inds2[2]+d_hi[2]) -
		   get_grid_d(inds2[0],inds2[1]+d_hi[1],inds2[2]-d_lo[2]) + get_grid_d(inds2[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));	//  d2V/dydz
	
	    b[56+i0] = dscales[0] * dscales[1] * dscales[2]					// d3V/dxdydz
		* (get_grid_d(inds2[0]+d_hi[0],inds2[1]+d_hi[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]+d_hi[0],inds2[1]+d_hi[1],inds2[2]-d_lo[2]) -
		   get_grid_d(inds2[0]+d_hi[0],inds2[1]-d_lo[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1]+d_hi[1],inds2[2]+d_hi[2]) +
		   get_grid_d(inds2[0]+d_hi[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]) + get_grid_d(inds2[0]-d_lo[0],inds2[1]+d_hi[1],inds2[2]-d_lo[2]) +
		   get_grid_d(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]+d_hi[2]) - get_grid_d(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));
	}
	
	if (inds2[0] == 1 && inds2[1] == 1 && inds2[2] == 0) {
	    DebugM(1, "Sub V = " << b[i0] << "\n");
	
	    DebugM(1, "Sub dV/dx = " << b[8+i0] << "\n");
	    DebugM(1, "Sub dV/dy = " << b[16+i0] << "\n");
	    DebugM(1, "Sub dV/dz = " << b[24+i0] << "\n");
	
	    DebugM(1, "Sub d2V/dxdy = " << b[32+i0] << "\n");
	    DebugM(1, "Sub d2V/dxdz = " << b[40+i0] << "\n");
	    DebugM(1, "Sub d2V/dydz = " << b[48+i0] << "\n");
	
	    DebugM(1, "Sub d3V/dxdydz = " << b[56+i0] << "\n" << endi);
	}
    }
}


/*********************/
/* GRIDFORCELITEGRID */
/*********************/

GridforceLiteGrid::GridforceLiteGrid(int gridnum)
{
    mygridnum = gridnum;
    grid = NULL;
    type = GridforceGridTypeLite;
}


GridforceLiteGrid::~GridforceLiteGrid()
{
    delete[] grid;
}


void GridforceLiteGrid::initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams)
{
    // cheat and use GridforceFullMainGrid to read the file
    GridforceFullMainGrid *tmp_grid = new GridforceFullMainGrid(mygridnum);
    tmp_grid->initialize(potfilename, simParams, mgridParams, 1);
    
    if (tmp_grid->get_total_grids() != 1) {
	NAMD_die("Cannot use gridforcelite option with multi-resolution grid!");
    }
    
    // save file name so that grid can be re-read via Tcl
    strcpy(filename, potfilename);
    
    // copy parameters
    k[0] = tmp_grid->get_k0();
    k[1] = tmp_grid->get_k1();
    k[2] = tmp_grid->get_k2();
    k[3] = 4;	// for V, dV/dx, dV/dy, dV/dz grids
    origin = tmp_grid->get_origin();
    center = tmp_grid->get_center();
    e = tmp_grid->get_e();
    inv = tmp_grid->get_inv();
    scale = tmp_grid->get_scale();
    
    // calculate rest of parameters
    size = k[0] * k[1] * k[2] * k[3];
    dk[0] = k[1] * k[2] * k[3];
    dk[1] = k[2] * k[3];
    dk[2] = k[3];
    dk[3] = 1;
    
    // copy the potential grid
    delete[] grid;
    grid = new float[size];
    for (int i0 = 0; i0 < k[0]; i0++) {
	for (int i1 = 0; i1 < k[1]; i1++) {
	    for (int i2 = 0; i2 < k[2]; i2++) {
		float V = tmp_grid->get_grid(i0, i1, i2);
		set_grid(i0, i1, i2, 0, V);
		DebugM(1, "V[" << i0 << "," << i1 << "," << i2 << "] = " << get_grid(i0, i1, i2, 0) << "(" << V << ")\n" << endi);
	    }
	}
    }
    
    delete tmp_grid;
    
    compute_derivative_grids();
}


void GridforceLiteGrid::compute_derivative_grids(void)
{
    // calculate derivative grids
    // separate loop so all V values have been set already
    for (int i0 = 0; i0 < k[0]; i0++) {
	for (int i1 = 0; i1 < k[1]; i1++) {
	    for (int i2 = 0; i2 < k[2]; i2++) {
		float dx, dy, dz;
		if (i0 == 0 || i0 == k[0]-1 || i1 == 0 || i1 == k[1]-1 || i2 == 0 || i2 == k[2]-1) {
		    // on edge, set ALL derivatives to zero (make up for lack of padding)
		    dx = 0;
		    dy = 0;
		    dz = 0;
		} else {
		    dx = 0.5 * (get_grid_d(i0+1,i1,i2,0) - get_grid_d(i0-1,i1,i2,0));
		    dy = 0.5 * (get_grid_d(i0,i1+1,i2,0) - get_grid_d(i0,i1-1,i2,0));
		    dz = 0.5 * (get_grid_d(i0,i1,i2+1,0) - get_grid_d(i0,i1,i2-1,0));
		}
		set_grid(i0, i1, i2, 1, dx);
		set_grid(i0, i1, i2, 2, dy);
		set_grid(i0, i1, i2, 3, dz);
		DebugM(1, "dx[" << i0 << "," << i1 << "," << i2 << "] = " << get_grid(i0, i1, i2, 1) << "(" << dx << ")\n" << endi);
		DebugM(1, "dy[" << i0 << "," << i1 << "," << i2 << "] = " << get_grid(i0, i1, i2, 2) << "(" << dy << ")\n" << endi);
		DebugM(1, "dz[" << i0 << "," << i1 << "," << i2 << "] = " << get_grid(i0, i1, i2, 3) << "(" << dz << ")\n" << endi);
	    }
	}
    }
}


void GridforceLiteGrid::reinitialize(SimParameters *simParams, MGridforceParams *mgridParams)
{
    initialize(filename, simParams, mgridParams);
}


void GridforceLiteGrid::pack(MOStream *msg) const
{
    msg->put(4*sizeof(int), (char*)k);
    msg->put(size);
    msg->put(4*sizeof(long int), (char*)dk);
    
    msg->put(sizeof(Vector), (char*)&origin);
    msg->put(sizeof(Vector), (char*)&center);
    msg->put(sizeof(Tensor), (char*)&e);
    msg->put(sizeof(Tensor), (char*)&inv);
    msg->put(sizeof(Vector), (char*)&scale);
    msg->put(sizeof(Bool), (char*)&checksize);
    
    msg->put(129*sizeof(char), (char*)filename);
    
    msg->put(size*sizeof(float), (char*)grid);
}


void GridforceLiteGrid::unpack(MIStream *msg)
{
    delete[] grid;
    grid = NULL;

    msg->get(4*sizeof(int), (char*)k);
    msg->get(size);
    msg->get(4*sizeof(long int), (char*)dk);
    
    msg->get(sizeof(Vector), (char*)&origin);
    msg->get(sizeof(Vector), (char*)&center);
    msg->get(sizeof(Tensor), (char*)&e);
    msg->get(sizeof(Tensor), (char*)&inv);
    msg->get(sizeof(Vector), (char*)&scale);
    msg->get(sizeof(Bool), (char*)&checksize);
    
    msg->get(129*sizeof(char), (char*)filename);
    
    if (size) {
	grid = new float[size];
	msg->get(size*sizeof(float), (char*)grid);
    }
}


long int GridforceLiteGrid::get_all_gridvals(float** all_gridvals) const
{
    // Creates a flat array of all grid values and puts it in the
    // value pointed to by the 'all_gridvals' argument. Returns the
    // resulting array size. Caller is responsible for destroying the
    // array via 'delete[]'
    
    DebugM(4, "GridforceLiteGrid::get_all_gridvals called\n" << endi);
    
    long int sz = size;
    DebugM(4, "size = " << sz << "\n" << endi);
    
    float *grid_vals = new float[sz];
    long int idx = 0;
    for (long int i = 0; i < size; i++) {
	grid_vals[idx++] = grid[i];
    }
    CmiAssert(idx == sz);
    
    *all_gridvals = grid_vals;
    
    DebugM(4, "GridforceLiteGrid::get_all_gridvals finished\n" << endi);
    
    return sz;
}


void GridforceLiteGrid::set_all_gridvals(float* all_gridvals, long int sz)
{
    DebugM(4, "GridforceLiteGrid::set_all_gridvals called\n" << endi);
    
    long int sz_calc = size;
    CmiAssert(sz == sz_calc);
    
    long int idx = 0;
    for (long int i = 0; i < size; i++) {
	grid[i] = all_gridvals[idx++];
    }
    CmiAssert(idx == sz);
    
    //compute_derivative_grids();	// not needed if we're sending all 4 grids
    
    DebugM(4, "GridforceLiteGrid::set_all_gridvals finished\n" << endi);
}
