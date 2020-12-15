/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GRIDFORCEGRID_INL
#define GRIDFORCEGRID_INL

#include "GridForceGrid.h"

inline int GridforceFullBaseGrid::compute_VdV(Position pos, float &V, Vector &dV) const
{
    //SimParameters *simParams = Node::Object()->simParameters;
    int inds[3];
    Vector g, dg;
    Vector gapscale = Vector(1, 1, 1);
    
    int err = get_inds(pos, inds, dg, gapscale);
    if (err) {
	return -1;
    }
    
    DebugM(1, "gapscale = " << gapscale << "\n");
    DebugM(1, "dg = " << dg << "\n");
    DebugM(1, "ind + dg = " << inds[0]+dg[0] << " " << inds[1]+dg[1] << " " << inds[2]+dg[2] << "\n");
    DebugM(3, "compute_VdV: generation = " << generation << "\n" << endi);
    
    // Pass to subgrid if one exists here
    for (int i = 0; i < numSubgrids; i++) {
	if (((inds[0] >= subgrids[i]->pmin[0] && inds[0] <= subgrids[i]->pmax[0]) || subgrids[i]->cont[0]) &&
	    ((inds[1] >= subgrids[i]->pmin[1] && inds[1] <= subgrids[i]->pmax[1]) || subgrids[i]->cont[1]) &&
	    ((inds[2] >= subgrids[i]->pmin[2] && inds[2] <= subgrids[i]->pmax[2]) || subgrids[i]->cont[2]))
	{
	    return subgrids[i]->compute_VdV(pos, V, dV);
	}
    }
    
    // Compute b
    float b[64];	// Matrix of values at 8 box corners
    compute_b(b, inds, gapscale);
    for (int j = 0; j < 64; j++) DebugM(1, "b[" << j << "] = " << b[j] << "\n" << endi);
    
    // Compute a
    float a[64];
    compute_a(a, b);
    for (int j = 0; j < 64; j++) DebugM(1, "a[" << j << "] = " << a[j] << "\n" << endi);
	    
    // Calculate powers of x, y, z for later use
    // e.g. x[2] = x^2
    float x[4], y[4], z[4];
    x[0] = 1; y[0] = 1; z[0] = 1;
    for (int j = 1; j < 4; j++) {
	x[j] = x[j-1] * dg.x;
	y[j] = y[j-1] * dg.y;
	z[j] = z[j-1] * dg.z;
    }
    
    V = compute_V(a, x, y, z);
    dV = Tensor::diagonal(gapscale) * (compute_dV(a, x, y, z) * inv);
    
    return 0;
}


inline int GridforceLiteGrid::compute_VdV(Position pos, float &V, Vector &dV) const
{
    int inds[3];
    Vector g, dg;
    
    int err = get_inds(pos, inds, dg);
    if (err) {
	return -1;
    }
    
    float wts[4][8];
    float results[4];
    
    // compute_wts(wts, dg);
    // wts[0][0] = (1-dg.x) * (1-dg.y) * (1-dg.z);
    // wts[0][1] = (1-dg.x) * (1-dg.y) *   dg.z;
    // wts[0][2] = (1-dg.x) *   dg.y   * (1-dg.z);
    // wts[0][3] = (1-dg.x) *   dg.y   *   dg.z;
    // wts[0][4] =   dg.x   * (1-dg.y) * (1-dg.z);
    // wts[0][5] =   dg.x   * (1-dg.y) *   dg.z;
    // wts[0][6] =   dg.x   *   dg.y   * (1-dg.z);
    // wts[0][7] =   dg.x   *   dg.y   *   dg.z;

    int i = 1;
    wts[i][0] = -(1-dg.y) * (1-dg.z);
    wts[i][1] = -(1-dg.y) *   dg.z;
    wts[i][2] = -  dg.y   * (1-dg.z);
    wts[i][3] = -  dg.y   *   dg.z;
    for (int j=0; j<4; j++) wts[i][j+4] = -wts[i][j];

    i = 2;
    wts[i][0] = -(1-dg.x) * (1-dg.z);
    wts[i][1] = -(1-dg.x) *   dg.z;
    wts[i][2] = -wts[i][0];
    wts[i][3] = -wts[i][1];
    wts[i][4] =   - dg.x  * (1-dg.z);
    wts[i][5] =   - dg.x  *   dg.z;
    wts[i][6] = -wts[i][4];
    wts[i][7] = -wts[i][5];

    i = 3;
    wts[i][0] = - (1-dg.x) * (1-dg.y);
    wts[i][1] = -wts[i][0];
    wts[i][2] = - (1-dg.x) *   dg.y  ;
    wts[i][3] = -wts[i][2];
    wts[i][4] = - dg.x     * (1-dg.y);
    wts[i][5] = -wts[i][4];
    wts[i][6] = - dg.x     *   dg.y  ;
    wts[i][7] = -wts[i][6];

    i = 0;
    for (int j=0; j<4; j++) wts[i][j]   = (1-dg.x) * wts[i+1][j+4];
    for (int j=0; j<4; j++) wts[i][j+4] =   dg.x   * wts[i+1][j+4];    

    for (i = 0; i < 4; i++) {
	results[i] = linear_interpolate(inds[0], inds[1], inds[2], 0, wts[i]);
    }
    
    V = results[0];
    dV = Vector(results[1], results[2], results[3]) * inv;
    
    return 0;
}


inline int GridforceFullBaseGrid::get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const
{
    Vector p = pos - origin;
    Vector g;
    
    g = inv * p;
    
    for (int i = 0; i < 3; i++) {
	inds[i] = (int)floor(g[i]);
	dg[i] = g[i] - inds[i];
    }
    
    for (int i = 0; i < 3; i++) {
	if (inds[i] < 0 || inds[i] >= k[i]-1) {
	    if (cont[i]) inds[i] = k[i]-1;
	    else return -1;	// Outside potential and grid is not continuous
	}
	if (cont[i] && inds[i] == k[i]-1) {
	    // Correct for non-unit spacing between continuous grid images
	    gapscale[i] *= gapinv[i];
	    if (g[i] < 0.0) dg[i] = 1.0 + g[i]*gapinv[i]; // = (gap[i] + g[i]) * gapinv[i]
	    else dg[i] = (g[i] - inds[i]) * gapinv[i];
	}
    }
    
    return 0;
}


inline float GridforceFullBaseGrid::compute_V(float *a, float *x, float *y, float *z) const
{
    float V = 0.0;
    long int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		V += a[ind] * x[j] * y[k] * z[l];
		ind++;
	    }
	}
    }
    return V;
}


inline Vector GridforceFullBaseGrid::compute_dV(float *a, float *x, float *y, float *z) const
{
    Vector dV = 0;
    long int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		if (j > 0) dV.x += a[ind] * j * x[j-1] * y[k]   * z[l];		// dV/dx
		if (k > 0) dV.y += a[ind] * k * x[j]   * y[k-1] * z[l];		// dV/dy
		if (l > 0) dV.z += a[ind] * l * x[j]   * y[k]   * z[l-1];	// dV/dz
		ind++;
	    }
	}
    }
    return dV;
}


inline Vector GridforceFullBaseGrid::compute_d2V(float *a, float *x, float *y, float *z) const
{
    Vector d2V = 0;
    int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		if (j > 0 && k > 0) d2V.x += a[ind] * j * k * x[j-1] * y[k-1] * z[l];	// d2V/dxdy
		if (j > 0 && l > 0) d2V.y += a[ind] * j * l * x[j-1] * y[k]   * z[l-1];	// d2V/dxdz
		if (k > 0 && l > 0) d2V.z += a[ind] * k * l * x[j]   * y[k-1] * z[l-1];	// d2V/dydz
		ind++;
	    }
	}
    }
    return d2V;
}


inline float GridforceFullBaseGrid::compute_d3V(float *a, float *x, float *y, float *z) const
{
    float d3V = 0.0;
    long int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		if (j > 0 && k > 0 && l > 0) d3V += a[ind] * j * k * l * x[j-1] * y[k-1] * z[l-1];	// d3V/dxdydz
		ind++;
	    }
	}
    }
    return d3V;
}


inline void GridforceFullBaseGrid::compute_a(float *a, float *b) const
{
    // Static sparse 64x64 matrix times vector ... nicer looking way than this?
    a[0] = b[0];
    a[1] = b[8];
    a[2] = -3*b[0] + 3*b[1] - 2*b[8] - b[9];
    a[3] = 2*b[0] - 2*b[1] + b[8] + b[9];
    a[4] = b[16];
    a[5] = b[32];
    a[6] = -3*b[16] + 3*b[17] - 2*b[32] - b[33];
    a[7] = 2*b[16] - 2*b[17] + b[32] + b[33];
    a[8] = -3*b[0] + 3*b[2] - 2*b[16] - b[18];
    a[9] = -3*b[8] + 3*b[10] - 2*b[32] - b[34];
    a[10] = 9*b[0] - 9*b[1] - 9*b[2] + 9*b[3] + 6*b[8] + 3*b[9] - 6*b[10] - 3*b[11]
	+ 6*b[16] - 6*b[17] + 3*b[18] - 3*b[19] + 4*b[32] + 2*b[33] + 2*b[34] + b[35];
    a[11] = -6*b[0] + 6*b[1] + 6*b[2] - 6*b[3] - 3*b[8] - 3*b[9] + 3*b[10] + 3*b[11]
	- 4*b[16] + 4*b[17] - 2*b[18] + 2*b[19] - 2*b[32] - 2*b[33] - b[34] - b[35];
    a[12] = 2*b[0] - 2*b[2] + b[16] + b[18];
    a[13] = 2*b[8] - 2*b[10] + b[32] + b[34];
    a[14] = -6*b[0] + 6*b[1] + 6*b[2] - 6*b[3] - 4*b[8] - 2*b[9] + 4*b[10] + 2*b[11]
	- 3*b[16] + 3*b[17] - 3*b[18] + 3*b[19] - 2*b[32] - b[33] - 2*b[34] - b[35];
    a[15] = 4*b[0] - 4*b[1] - 4*b[2] + 4*b[3] + 2*b[8] + 2*b[9] - 2*b[10] - 2*b[11]
	+ 2*b[16] - 2*b[17] + 2*b[18] - 2*b[19] + b[32] + b[33] + b[34] + b[35];
    a[16] = b[24];
    a[17] = b[40];
    a[18] = -3*b[24] + 3*b[25] - 2*b[40] - b[41];
    a[19] = 2*b[24] - 2*b[25] + b[40] + b[41];
    a[20] = b[48];
    a[21] = b[56];
    a[22] = -3*b[48] + 3*b[49] - 2*b[56] - b[57];
    a[23] = 2*b[48] - 2*b[49] + b[56] + b[57];
    a[24] = -3*b[24] + 3*b[26] - 2*b[48] - b[50];
    a[25] = -3*b[40] + 3*b[42] - 2*b[56] - b[58];
    a[26] = 9*b[24] - 9*b[25] - 9*b[26] + 9*b[27] + 6*b[40] + 3*b[41] - 6*b[42] - 3*b[43]
	+ 6*b[48] - 6*b[49] + 3*b[50] - 3*b[51] + 4*b[56] + 2*b[57] + 2*b[58] + b[59];
    a[27] = -6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 3*b[40] - 3*b[41] + 3*b[42] + 3*b[43]
	- 4*b[48] + 4*b[49] - 2*b[50] + 2*b[51] - 2*b[56] - 2*b[57] - b[58] - b[59];
    a[28] = 2*b[24] - 2*b[26] + b[48] + b[50];
    a[29] = 2*b[40] - 2*b[42] + b[56] + b[58];
    a[30] = -6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 4*b[40] - 2*b[41] + 4*b[42] + 2*b[43]
	- 3*b[48] + 3*b[49] - 3*b[50] + 3*b[51] - 2*b[56] - b[57] - 2*b[58] - b[59];
    a[31] = 4*b[24] - 4*b[25] - 4*b[26] + 4*b[27] + 2*b[40] + 2*b[41] - 2*b[42] - 2*b[43]
	+ 2*b[48] - 2*b[49] + 2*b[50] - 2*b[51] + b[56] + b[57] + b[58] + b[59];
    a[32] = -3*b[0] + 3*b[4] - 2*b[24] - b[28];
    a[33] = -3*b[8] + 3*b[12] - 2*b[40] - b[44];
    a[34] = 9*b[0] - 9*b[1] - 9*b[4] + 9*b[5] + 6*b[8] + 3*b[9] - 6*b[12] - 3*b[13]
	+ 6*b[24] - 6*b[25] + 3*b[28] - 3*b[29] + 4*b[40] + 2*b[41] + 2*b[44] + b[45];
    a[35] = -6*b[0] + 6*b[1] + 6*b[4] - 6*b[5] - 3*b[8] - 3*b[9] + 3*b[12] + 3*b[13]
	- 4*b[24] + 4*b[25] - 2*b[28] + 2*b[29] - 2*b[40] - 2*b[41] - b[44] - b[45];
    a[36] = -3*b[16] + 3*b[20] - 2*b[48] - b[52];
    a[37] = -3*b[32] + 3*b[36] - 2*b[56] - b[60];
    a[38] = 9*b[16] - 9*b[17] - 9*b[20] + 9*b[21] + 6*b[32] + 3*b[33] - 6*b[36] - 3*b[37]
	+ 6*b[48] - 6*b[49] + 3*b[52] - 3*b[53] + 4*b[56] + 2*b[57] + 2*b[60] + b[61];
    a[39] = -6*b[16] + 6*b[17] + 6*b[20] - 6*b[21] - 3*b[32] - 3*b[33] + 3*b[36] + 3*b[37]
	- 4*b[48] + 4*b[49] - 2*b[52] + 2*b[53] - 2*b[56] - 2*b[57] - b[60] - b[61];
    a[40] = 9*b[0] - 9*b[2] - 9*b[4] + 9*b[6] + 6*b[16] + 3*b[18] - 6*b[20] - 3*b[22]
	+ 6*b[24] - 6*b[26] + 3*b[28] - 3*b[30] + 4*b[48] + 2*b[50] + 2*b[52] + b[54];
    a[41] = 9*b[8] - 9*b[10] - 9*b[12] + 9*b[14] + 6*b[32] + 3*b[34] - 6*b[36] - 3*b[38]
	+ 6*b[40] - 6*b[42] + 3*b[44] - 3*b[46] + 4*b[56] + 2*b[58] + 2*b[60] + b[62];
    a[42] = -27*b[0] + 27*b[1] + 27*b[2] - 27*b[3] + 27*b[4] - 27*b[5] - 27*b[6] + 27*b[7]
	- 18*b[8] - 9*b[9] + 18*b[10] + 9*b[11] + 18*b[12] + 9*b[13] - 18*b[14] - 9*b[15]
	- 18*b[16] + 18*b[17] - 9*b[18] + 9*b[19] + 18*b[20] - 18*b[21] + 9*b[22] - 9*b[23]
	- 18*b[24] + 18*b[25] + 18*b[26] - 18*b[27] - 9*b[28] + 9*b[29] + 9*b[30] - 9*b[31]
	- 12*b[32] - 6*b[33] - 6*b[34] - 3*b[35] + 12*b[36] + 6*b[37] + 6*b[38] + 3*b[39]
	- 12*b[40] - 6*b[41] + 12*b[42] + 6*b[43] - 6*b[44] - 3*b[45] + 6*b[46] + 3*b[47]
	- 12*b[48] + 12*b[49] - 6*b[50] + 6*b[51] - 6*b[52] + 6*b[53] - 3*b[54] + 3*b[55]
	- 8*b[56] - 4*b[57] - 4*b[58] - 2*b[59] - 4*b[60] - 2*b[61] - 2*b[62] - b[63];
    a[43] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
	+ 9*b[8] + 9*b[9] - 9*b[10] - 9*b[11] - 9*b[12] - 9*b[13] + 9*b[14] + 9*b[15]
	+ 12*b[16] - 12*b[17] + 6*b[18] - 6*b[19] - 12*b[20] + 12*b[21] - 6*b[22] + 6*b[23]
	+ 12*b[24] - 12*b[25] - 12*b[26] + 12*b[27] + 6*b[28] - 6*b[29] - 6*b[30] + 6*b[31]
	+ 6*b[32] + 6*b[33] + 3*b[34] + 3*b[35] - 6*b[36] - 6*b[37] - 3*b[38] - 3*b[39]
	+ 6*b[40] + 6*b[41] - 6*b[42] - 6*b[43] + 3*b[44] + 3*b[45] - 3*b[46] - 3*b[47]
	+ 8*b[48] - 8*b[49] + 4*b[50] - 4*b[51] + 4*b[52] - 4*b[53] + 2*b[54] - 2*b[55]
	+ 4*b[56] + 4*b[57] + 2*b[58] + 2*b[59] + 2*b[60] + 2*b[61] + b[62] + b[63];
    a[44] = -6*b[0] + 6*b[2] + 6*b[4] - 6*b[6] - 3*b[16] - 3*b[18] + 3*b[20] + 3*b[22]
	- 4*b[24] + 4*b[26] - 2*b[28] + 2*b[30] - 2*b[48] - 2*b[50] - b[52] - b[54];
    a[45] = -6*b[8] + 6*b[10] + 6*b[12] - 6*b[14] - 3*b[32] - 3*b[34] + 3*b[36] + 3*b[38]
	- 4*b[40] + 4*b[42] - 2*b[44] + 2*b[46] - 2*b[56] - 2*b[58] - b[60] - b[62];
    a[46] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
	+ 12*b[8] + 6*b[9] - 12*b[10] - 6*b[11] - 12*b[12] - 6*b[13] + 12*b[14] + 6*b[15]
	+ 9*b[16] - 9*b[17] + 9*b[18] - 9*b[19] - 9*b[20] + 9*b[21] - 9*b[22] + 9*b[23]
	+ 12*b[24] - 12*b[25] - 12*b[26] + 12*b[27] + 6*b[28] - 6*b[29] - 6*b[30] + 6*b[31]
	+ 6*b[32] + 3*b[33] + 6*b[34] + 3*b[35] - 6*b[36] - 3*b[37] - 6*b[38] - 3*b[39]
	+ 8*b[40] + 4*b[41] - 8*b[42] - 4*b[43] + 4*b[44] + 2*b[45] - 4*b[46] - 2*b[47]
	+ 6*b[48] - 6*b[49] + 6*b[50] - 6*b[51] + 3*b[52] - 3*b[53] + 3*b[54] - 3*b[55]
	+ 4*b[56] + 2*b[57] + 4*b[58] + 2*b[59] + 2*b[60] + b[61] + 2*b[62] + b[63];
    a[47] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
	- 6*b[8] - 6*b[9] + 6*b[10] + 6*b[11] + 6*b[12] + 6*b[13] - 6*b[14] - 6*b[15]
	- 6*b[16] + 6*b[17] - 6*b[18] + 6*b[19] + 6*b[20] - 6*b[21] + 6*b[22] - 6*b[23]
	- 8*b[24] + 8*b[25] + 8*b[26] - 8*b[27] - 4*b[28] + 4*b[29] + 4*b[30] - 4*b[31]
	- 3*b[32] - 3*b[33] - 3*b[34] - 3*b[35] + 3*b[36] + 3*b[37] + 3*b[38] + 3*b[39]
	- 4*b[40] - 4*b[41] + 4*b[42] + 4*b[43] - 2*b[44] - 2*b[45] + 2*b[46] + 2*b[47]
	- 4*b[48] + 4*b[49] - 4*b[50] + 4*b[51] - 2*b[52] + 2*b[53] - 2*b[54] + 2*b[55]
	- 2*b[56] - 2*b[57] - 2*b[58] - 2*b[59] - b[60] - b[61] - b[62] - b[63];
    a[48] = 2*b[0] - 2*b[4] + b[24] + b[28];
    a[49] = 2*b[8] - 2*b[12] + b[40] + b[44];
    a[50] = -6*b[0] + 6*b[1] + 6*b[4] - 6*b[5] - 4*b[8] - 2*b[9] + 4*b[12] + 2*b[13]
	- 3*b[24] + 3*b[25] - 3*b[28] + 3*b[29] - 2*b[40] - b[41] - 2*b[44] - b[45];
    a[51] = 4*b[0] - 4*b[1] - 4*b[4] + 4*b[5] + 2*b[8] + 2*b[9] - 2*b[12] - 2*b[13]
	+ 2*b[24] - 2*b[25] + 2*b[28] - 2*b[29] + b[40] + b[41] + b[44] + b[45];
    a[52] = 2*b[16] - 2*b[20] + b[48] + b[52];
    a[53] = 2*b[32] - 2*b[36] + b[56] + b[60];
    a[54] = -6*b[16] + 6*b[17] + 6*b[20] - 6*b[21] - 4*b[32] - 2*b[33] + 4*b[36] + 2*b[37]
	- 3*b[48] + 3*b[49] - 3*b[52] + 3*b[53] - 2*b[56] - b[57] - 2*b[60] - b[61];
    a[55] = 4*b[16] - 4*b[17] - 4*b[20] + 4*b[21] + 2*b[32] + 2*b[33] - 2*b[36] - 2*b[37]
	+ 2*b[48] - 2*b[49] + 2*b[52] - 2*b[53] + b[56] + b[57] + b[60] + b[61];
    a[56] = -6*b[0] + 6*b[2] + 6*b[4] - 6*b[6] - 4*b[16] - 2*b[18] + 4*b[20] + 2*b[22]
	- 3*b[24] + 3*b[26] - 3*b[28] + 3*b[30] - 2*b[48] - b[50] - 2*b[52] - b[54];
    a[57] = -6*b[8] + 6*b[10] + 6*b[12] - 6*b[14] - 4*b[32] - 2*b[34] + 4*b[36] + 2*b[38]
	- 3*b[40] + 3*b[42] - 3*b[44] + 3*b[46] - 2*b[56] - b[58] - 2*b[60] - b[62];
    a[58] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
	+ 12*b[8] + 6*b[9] - 12*b[10] - 6*b[11] - 12*b[12] - 6*b[13] + 12*b[14] + 6*b[15]
	+ 12*b[16] - 12*b[17] + 6*b[18] - 6*b[19] - 12*b[20] + 12*b[21] - 6*b[22] + 6*b[23]
	+ 9*b[24] - 9*b[25] - 9*b[26] + 9*b[27] + 9*b[28] - 9*b[29] - 9*b[30] + 9*b[31]
	+ 8*b[32] + 4*b[33] + 4*b[34] + 2*b[35] - 8*b[36] - 4*b[37] - 4*b[38] - 2*b[39]
	+ 6*b[40] + 3*b[41] - 6*b[42] - 3*b[43] + 6*b[44] + 3*b[45] - 6*b[46] - 3*b[47]
	+ 6*b[48] - 6*b[49] + 3*b[50] - 3*b[51] + 6*b[52] - 6*b[53] + 3*b[54] - 3*b[55]
	+ 4*b[56] + 2*b[57] + 2*b[58] + b[59] + 4*b[60] + 2*b[61] + 2*b[62] + b[63];
    a[59] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
	- 6*b[8] - 6*b[9] + 6*b[10] + 6*b[11] + 6*b[12] + 6*b[13] - 6*b[14] - 6*b[15]
	- 8*b[16] + 8*b[17] - 4*b[18] + 4*b[19] + 8*b[20] - 8*b[21] + 4*b[22] - 4*b[23]
	- 6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 6*b[28] + 6*b[29] + 6*b[30] - 6*b[31]
	- 4*b[32] - 4*b[33] - 2*b[34] - 2*b[35] + 4*b[36] + 4*b[37] + 2*b[38] + 2*b[39]
	- 3*b[40] - 3*b[41] + 3*b[42] + 3*b[43] - 3*b[44] - 3*b[45] + 3*b[46] + 3*b[47]
	- 4*b[48] + 4*b[49] - 2*b[50] + 2*b[51] - 4*b[52] + 4*b[53] - 2*b[54] + 2*b[55]
	- 2*b[56] - 2*b[57] - b[58] - b[59] - 2*b[60] - 2*b[61] - b[62] - b[63];
    a[60] = 4*b[0] - 4*b[2] - 4*b[4] + 4*b[6] + 2*b[16] + 2*b[18] - 2*b[20] - 2*b[22]
	+ 2*b[24] - 2*b[26] + 2*b[28] - 2*b[30] + b[48] + b[50] + b[52] + b[54];
    a[61] = 4*b[8] - 4*b[10] - 4*b[12] + 4*b[14] + 2*b[32] + 2*b[34] - 2*b[36] - 2*b[38]
	+ 2*b[40] - 2*b[42] + 2*b[44] - 2*b[46] + b[56] + b[58] + b[60] + b[62];
    a[62] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
	- 8*b[8] - 4*b[9] + 8*b[10] + 4*b[11] + 8*b[12] + 4*b[13] - 8*b[14] - 4*b[15]
	- 6*b[16] + 6*b[17] - 6*b[18] + 6*b[19] + 6*b[20] - 6*b[21] + 6*b[22] - 6*b[23]
	- 6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 6*b[28] + 6*b[29] + 6*b[30] - 6*b[31]
	- 4*b[32] - 2*b[33] - 4*b[34] - 2*b[35] + 4*b[36] + 2*b[37] + 4*b[38] + 2*b[39]
	- 4*b[40] - 2*b[41] + 4*b[42] + 2*b[43] - 4*b[44] - 2*b[45] + 4*b[46] + 2*b[47]
	- 3*b[48] + 3*b[49] - 3*b[50] + 3*b[51] - 3*b[52] + 3*b[53] - 3*b[54] + 3*b[55]
	- 2*b[56] - b[57] - 2*b[58] - b[59] - 2*b[60] - b[61] - 2*b[62] - b[63];
    a[63] = 8*b[0] - 8*b[1] - 8*b[2] + 8*b[3] - 8*b[4] + 8*b[5] + 8*b[6] - 8*b[7]
	+ 4*b[8] + 4*b[9] - 4*b[10] - 4*b[11] - 4*b[12] - 4*b[13] + 4*b[14] + 4*b[15]
	+ 4*b[16] - 4*b[17] + 4*b[18] - 4*b[19] - 4*b[20] + 4*b[21] - 4*b[22] + 4*b[23]
	+ 4*b[24] - 4*b[25] - 4*b[26] + 4*b[27] + 4*b[28] - 4*b[29] - 4*b[30] + 4*b[31]
	+ 2*b[32] + 2*b[33] + 2*b[34] + 2*b[35] - 2*b[36] - 2*b[37] - 2*b[38] - 2*b[39]
	+ 2*b[40] + 2*b[41] - 2*b[42] - 2*b[43] + 2*b[44] + 2*b[45] - 2*b[46] - 2*b[47]
	+ 2*b[48] - 2*b[49] + 2*b[50] - 2*b[51] + 2*b[52] - 2*b[53] + 2*b[54] - 2*b[55]
	+ b[56] + b[57] + b[58] + b[59] + b[60] + b[61] + b[62] + b[63];
}


inline int GridforceLiteGrid::get_inds(Position pos, int *inds, Vector &dg) const
{
    Vector p = pos - origin;
    Vector g;
    
    g = inv * p;
    
    for (int i = 0; i < 3; i++) {
	inds[i] = (int)floor(g[i]);
	dg[i] = g[i] - inds[i];
    }
    
    for (int i = 0; i < 3; i++) {
	if (inds[i] < 0 || inds[i] >= k[i]-1) {
	    return -1;	// Outside potential and grid is not continuous
	}
    }
    
    return 0;
}


inline void GridforceLiteGrid::compute_wts(float *wts, const Vector &dg) const
{
    wts[0] = (1-dg.x) * (1-dg.y) * (1-dg.z);
    wts[1] = (1-dg.x) * (1-dg.y) *   dg.z;
    wts[2] = (1-dg.x) *   dg.y   * (1-dg.z);
    wts[3] = (1-dg.x) *   dg.y   *   dg.z;
    wts[4] =   dg.x   * (1-dg.y) * (1-dg.z);
    wts[5] =   dg.x   * (1-dg.y) *   dg.z;
    wts[6] =   dg.x   *   dg.y   * (1-dg.z);
    wts[7] =   dg.x   *   dg.y   *   dg.z;
    DebugM(2, "dg = " << dg << "\n" << endi);
}


inline float GridforceLiteGrid::linear_interpolate(int i0, int i1, int i2, int i3, const float *wts) const
{
#ifdef DEBUGM
    float vals[8];
    vals[0] = get_grid(i0,   i1,   i2,   i3);
    vals[1] = get_grid(i0,   i1,   i2+1, i3);
    vals[2] = get_grid(i0,   i1+1, i2,   i3);
    vals[3] = get_grid(i0,   i1+1, i2+1, i3);
    vals[4] = get_grid(i0+1, i1,   i2,   i3);
    vals[5] = get_grid(i0+1, i1,   i2+1, i3);
    vals[6] = get_grid(i0+1, i1+1, i2,   i3);
    vals[7] = get_grid(i0+1, i1+1, i2+1, i3);
    
    switch (i3) {
    case 0:
	DebugM(2, "V\n" << endi);
	break;
    case 1:
	DebugM(2, "dV/dx\n" << endi);
	break;
    case 2:
	DebugM(2, "dV/dy\n" << endi);
	break;
    case 3:
	DebugM(2, "dV/dz\n" << endi);
	break;
    }
    
    for (int i = 0; i < 8; i++) {
	DebugM(2, "vals[" << i << "] = " << vals[i] << " wts[" << i << "] = " << wts[i] << "\n" << endi);
    }
#endif
    
    float result =
	wts[0] * get_grid(i0,   i1,   i2,   i3) +
	wts[1] * get_grid(i0,   i1,   i2+1, i3) +
	wts[2] * get_grid(i0,   i1+1, i2,   i3) +
	wts[3] * get_grid(i0,   i1+1, i2+1, i3) +
	wts[4] * get_grid(i0+1, i1,   i2,   i3) +
	wts[5] * get_grid(i0+1, i1,   i2+1, i3) +
	wts[6] * get_grid(i0+1, i1+1, i2,   i3) +
	wts[7] * get_grid(i0+1, i1+1, i2+1, i3);
    
    DebugM(2, "result = " << result << "\n" << endi);
    
    return result;
}


inline Position GridforceGrid::wrap_position(const Position &pos, const Lattice &lattice)
{
    // Wrap 'pos' about grid center, using periodic cell information in 'lattice'
    // Position pos_wrapped = pos;
    // Position center = get_center();
    // pos_wrapped += lattice.wrap_delta(pos);
    // pos_wrapped += lattice.delta(pos_wrapped, center) - (pos_wrapped - center);
    
    Position pos_wrapped = pos + lattice.wrap_delta(pos - get_center() + lattice.origin());
    
    return pos_wrapped;
}

#endif
