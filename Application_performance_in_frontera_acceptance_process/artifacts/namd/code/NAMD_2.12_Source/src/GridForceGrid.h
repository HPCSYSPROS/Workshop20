/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GRIDFORCEGRID_H
#define GRIDFORCEGRID_H

#include <set>

#include "Vector.h"
#include "Tensor.h"
#include "SimParameters.h"
#include "NamdTypes.h"
#include "MStream.h"
#include "charm++.h"
//#include "ComputeGridForce.h"

#include "MGridforceParams.h"


class GridforceFullMainGrid;
class GridforceFullSubGrid;


// GridforceGrid is now an abstract class to act as an interface to both GridforceFullMainGrid's and GridforceLiteGrid's
class GridforceGrid {
public:
    static GridforceGrid * new_grid(int gridnum, char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams);
    virtual ~GridforceGrid();
    
    virtual void initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams) = 0;
    virtual void reinitialize(SimParameters *simParams, MGridforceParams *mgridParams) = 0;
    
    virtual Position get_center(void) const = 0;
    virtual Position get_origin(void) const = 0;
    virtual Tensor get_e (void) const = 0;
    virtual Tensor get_inv(void) const = 0;
    virtual Vector get_scale(void) const = 0;
    virtual Bool get_checksize(void) const = 0;
    virtual int get_k0(void) const = 0;
    virtual int get_k1(void) const = 0; 
    virtual int get_k2(void) const = 0;
    virtual int get_total_grids(void) const = 0;
    
    virtual long int get_all_gridvals(float** all_gridvals) const = 0;
    virtual void set_all_gridvals(float* all_gridvals, long int sz) = 0;
    virtual void set_scale(Vector s) = 0;
    
    Position wrap_position(const Position &pos, const Lattice &lattice);
    bool fits_lattice(const Lattice &lattice);
    
    inline int compute_VdV(Position pos, float &V, Vector &dV) const { return -1; }

    static void pack_grid(GridforceGrid *grid, MOStream *msg);
    static GridforceGrid * unpack_grid(int gridnum, MIStream *msg);
    
    typedef enum {
	GridforceGridTypeUndefined = 0,
	GridforceGridTypeFull,
	GridforceGridTypeLite
    } GridforceGridType;
    
    inline GridforceGridType get_grid_type(void) { return type; }

protected:    
    virtual void pack(MOStream *msg) const = 0;
    virtual void unpack(MIStream *msg) = 0;
    
    Position get_corner(int idx);
    
    GridforceGrid() { type = GridforceGridTypeUndefined; }
    GridforceGridType type;
    int mygridnum;
    
private:
    Vector corners[8];
};


class GridforceFullBaseGrid {
    friend class GridforceFullMainGrid;
    friend class GridforceFullSubGrid;
    
public:
    GridforceFullBaseGrid(void);
    virtual ~GridforceFullBaseGrid();
    
//    int request_box(Vector pos);
//    int get_box(Box *box, Vector pos) const;
  
    inline Position get_center(void) const { return center; }
    inline Position get_origin(void) const { return origin; }
    inline Tensor get_e (void) const { return e; }
    inline Tensor get_inv(void) const { return inv; }
    inline Vector get_scale(void) const { return scale; }
    inline Bool get_checksize(void) const { return checksize; }
    virtual int get_border(void) const = 0;
    
    inline float get_grid(int i0, int i1, int i2) const {
	return grid[grid_index(i0, i1, i2)];
    }
    inline double get_grid_d(int i0, int i1, int i2) const {
	return double(get_grid(i0, i1, i2));
    }
    inline void set_grid(int i0, int i1, int i2, float V) {
	grid[grid_index(i0, i1, i2)] = V;
    }
    
    inline void set_scale(Vector s) { scale = s; }
    
    int compute_VdV(Position pos, float &V, Vector &dV) const;
    
    inline int get_k0(void) const { return k[0]; }
    inline int get_k1(void) const { return k[1]; }
    inline int get_k2(void) const { return k[2]; }
    
protected:
    virtual void pack(MOStream *msg) const;
    virtual void unpack(MIStream *msg);
    
    struct GridIndices {
      int inds2;
      int dk_hi;
      int dk_lo;
      Bool zero_derivs;
    };
   
    // Utility functions
    void readHeader(SimParameters *simParams, MGridforceParams *mgridParams);
    
    inline long int grid_index(int i0, int i1, int i2) const {
	register int inds[3] = {i0, i1, i2};
#ifdef DEBUGM
	if (i0 < 0 || i0 >= k[0] || i1 < 0 || i1 >= k[1] || i2 < 0 || i2 >= k[2]) {
	    char buffer[256];
	    sprintf(buffer, "Bad grid index! (%d %d %d)", i0, i1, i2);
	    NAMD_bug(buffer);
	}
#endif
	return inds[0]*dk[0] + inds[1]*dk[1] + inds[2]*dk[2];
    }
    
    //virtual int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const = 0;
    int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_a(float *a, float *b) const;
    virtual void compute_b(float *b, int *inds, Vector gapscale)  const = 0;
    float compute_V(float *a, float *x, float *y, float *z) const;
    Vector compute_dV(float *a, float *x, float *y, float *z) const;
    Vector compute_d2V(float *a, float *x, float *y, float *z) const;
    float compute_d3V(float *a, float *x, float *y, float *z) const;
    
    void readSubgridHierarchy(FILE *poten, int &totalGrids);
    
    FILE *poten_fp;
    float *grid;	// Actual grid
    
    GridforceFullSubGrid **subgrids;
    int numSubgrids;
    int generation;	// Subgrid level (0 = main grid)
    
    // should move 'nopad' versions to maingrid only ... or not have them as ivars at all, why are they here?
    int k[3];		// Grid dimensions
    int k_nopad[3];	// Grid dimensions
    long int size;
    long int size_nopad;
    long int dk[3];
    long int dk_nopad[3];
    float factor;
    
    Position origin;	// Grid origin
    Position center;	// Center of grid (for wrapping)
    Tensor e;		// Grid unit vectors
    Tensor inv;		// Inverse of unit vectors
    
    double p_sum[3];     // Accumulators for sums
    double n_sum[3];
    double pad_p[3];	// Pad values (p = positive side, n = negative side) for each dimension
    double pad_n[3];
    Bool cont[3];	// Whether grid is continuous in each dimension
    float offset[3];	// Potential offset in each dimension
    float gap[3];	// Gap between images of grid in grid units for each dimension
    float gapinv[3];	// 1.0/gap

    Vector scale;
    Bool checksize;
};


class GridforceFullMainGrid : public GridforceGrid, public GridforceFullBaseGrid {
    friend class GridforceFullBaseGrid;
    friend class GridforceFullSubGrid;

public:
    explicit GridforceFullMainGrid(int gridnum);
    virtual ~GridforceFullMainGrid();
    
    void initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams, int border);
    inline void initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams) {
	initialize(potfilename, simParams, mgridParams, default_border);
    }
    void reinitialize(SimParameters *simParams, MGridforceParams *mgridParams);
    
    inline Position get_center(void) const { return GridforceFullBaseGrid::get_center(); };
    inline Position get_origin(void) const { return GridforceFullBaseGrid::get_origin(); };
    inline Tensor get_e (void) const { return GridforceFullBaseGrid::get_e(); };
    inline Tensor get_inv(void) const { return GridforceFullBaseGrid::get_inv(); };
    inline Vector get_scale(void) const { return GridforceFullBaseGrid::get_scale(); };
    inline Bool get_checksize(void) const { return GridforceFullBaseGrid::get_checksize(); };
    inline int get_k0(void) const { return GridforceFullBaseGrid::get_k0(); };
    inline int get_k1(void) const { return GridforceFullBaseGrid::get_k1(); };
    inline int get_k2(void) const { return GridforceFullBaseGrid::get_k2(); };
    inline int get_border(void) const { return border; }
    
    inline int compute_VdV(Position pos, float &V, Vector &dV) const { return GridforceFullBaseGrid::compute_VdV(pos, V, dV); };
    
    inline int get_total_grids(void) const { return totalGrids; }    
    inline void set_scale(Vector s) { scale = s; }

protected:
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    long int get_all_gridvals(float **all_gridvals) const;
    void set_all_gridvals(float *all_gridvals, long int sz);
    
    //int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_b(float *b, int *inds, Vector gapscale)  const;
    void buildSubgridsFlat(void);
    
    char filename[129];
    int totalGrids;
    GridforceFullSubGrid **subgrids_flat;
    //int mygridnum;
    
    static const int default_border = 1;
    int border;
};


class GridforceFullSubGrid : public GridforceFullBaseGrid {
    friend class GridforceFullBaseGrid;
    friend class GridforceFullMainGrid;

public:
    GridforceFullSubGrid(GridforceFullBaseGrid *parent_in);
    
    void initialize(SimParameters *simParams, MGridforceParams *mgridParams);
    
    inline int get_border(void) const { return 0; }
    
    inline Tensor tensorMult (const Tensor &t1, const Tensor &t2) {
	Tensor tmp;
	tmp.xx = t1.xx * t2.xx + t1.xy * t2.yx + t1.xz * t2.zx;
	tmp.xy = t1.xx * t2.xy + t1.xy * t2.yy + t1.xz * t2.zy;
	tmp.xz = t1.xx * t2.xz + t1.xy * t2.yz + t1.xz * t2.zz;
	tmp.yx = t1.yx * t2.xx + t1.yy * t2.yx + t1.yz * t2.zx;
	tmp.yy = t1.yx * t2.xy + t1.yy * t2.yy + t1.yz * t2.zy;
	tmp.yz = t1.yx * t2.xz + t1.yy * t2.yz + t1.yz * t2.zz;
	tmp.zx = t1.zx * t2.xx + t1.zy * t2.yx + t1.zz * t2.zx;
	tmp.zy = t1.zx * t2.xy + t1.zy * t2.yy + t1.zz * t2.zy;
	tmp.zz = t1.zx * t2.xz + t1.zy * t2.yz + t1.zz * t2.zz;
	return tmp;
    }
    
protected:
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    //int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_b(float *b, int *inds, Vector gapscale) const;
    void addToSubgridsFlat(void);
    
    // Utility numbers
    Tensor scale_dV;
    Tensor scale_d2V;
    float scale_d3V;
    
    GridforceFullBaseGrid *parent;
    int pmin[3], pmax[3];
    GridforceFullMainGrid *maingrid;
    int subgridIdx;
};


class GridforceLiteGrid : public GridforceGrid {
public:
    explicit GridforceLiteGrid(int gridnum);
    virtual ~GridforceLiteGrid();
    
    void initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams);
    void reinitialize(SimParameters *simParams, MGridforceParams *mgridParams);
    
    inline Position get_center(void) const { return center; }
    inline Position get_origin(void) const { return origin; }
    inline Tensor get_e (void) const { return e; }
    inline Tensor get_inv(void) const { return inv; }
    inline Vector get_scale(void) const { return scale; }
    inline Bool get_checksize(void) const { return checksize; }
    inline int get_k0(void) const { return k[0]; }
    inline int get_k1(void) const { return k[1]; }
    inline int get_k2(void) const { return k[2]; }
    inline int get_total_grids(void) const { return 1; }
    inline void set_scale(Vector s) { scale = s; }
    
    inline float get_grid(int i0, int i1, int i2, int i3) const {
	return grid[grid_index(i0, i1, i2, i3)];
    }
    inline double get_grid_d(int i0, int i1, int i2, int i3) const {
	return double(grid[grid_index(i0, i1, i2, i3)]);
    }
    inline void set_grid(int i0, int i1, int i2, int i3, float V) {
	grid[grid_index(i0, i1, i2, i3)] = V;
    }
    
    long int get_all_gridvals(float** all_gridvals) const;
    void set_all_gridvals(float* all_gridvals, long int sz);
    
    int compute_VdV(Position pos, float &V, Vector &dV) const;
    
protected:
    void compute_derivative_grids(void);
    void compute_wts(float *wts, const Vector &dg) const;
    int get_inds(Position pos, int *inds, Vector &dg) const;
    float linear_interpolate(int i0, int i1, int i2, int i3, const float *wts) const;
    
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    inline long int grid_index(int i0, int i1, int i2, int i3) const {
	// 'i3' is an index for the grid itself (0=V, 1=dV/dx, 2=dV/dy, 3=dV/dz)
	register int inds[4] = {i0, i1, i2, i3};
	return inds[0]*dk[0] + inds[1]*dk[1] + inds[2]*dk[2] + inds[3]*dk[3];
    }
    
    float *grid;
    
    int k[4];		// Grid dimensions ... 4th is always 4, for the different grid types
    long int size;
    long int dk[4];
    
    Position origin;	// Grid origin
    Position center;	// Center of grid (for wrapping)
    Tensor e;		// Grid unit vectors
    Tensor inv;		// Inverse of unit vectors
    
    Vector scale;
    Bool checksize;
    
    char filename[129];
};

#endif
