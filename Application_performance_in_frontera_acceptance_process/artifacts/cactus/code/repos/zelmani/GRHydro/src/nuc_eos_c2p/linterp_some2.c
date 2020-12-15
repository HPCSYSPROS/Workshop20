#include "nuc_eos_c.h"

void nuc_eosK_C_linterp_some2(double x, double y, double z,
			    double* f, double* ft, 
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt,struct nuc_eos_vars *struct_ptr,double* dlepsdlrho, double* dlepsdlt, double* dlPdlrho, double* dlPdlt
 			  ) 

{

//!
//!     purpose: interpolation of a function of three variables in an
//!              equidistant(!!!) table.
//!
//!     method:  8-point Lagrange linear interpolation formula
//!
//!     x        input vector of first  variable
//!     y        input vector of second variable
//!     z        input vector of third  variable
//!
//!     f        output vector of interpolated function values
//!
//!     ft       3d array of tabulated function values
//!     nx       x-dimension of table
//!     ny       y-dimension of table
//!     nz       z-dimension of table
//!     xt       vector of x-coordinates of table
//!     yt       vector of y-coordinates of table
//!     zt       vector of z-coordinates of table
//!
//!---------------------------------------------------------------------

  // helper variables
  double fh[11][nvars], delx, dely, delz, a[11];  // temp fix for opencl restriction on variable length arrays
  double dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi;
  int ix,iy,iz;
  
  // determine spacing parameters of equidistant (!!!) table
#if 0
  dx = (xt[nx-1] - xt[0]) / (1.0*(nx-1));
  dy = (yt[ny-1] - yt[0]) / (1.0*(ny-1));
  dz = (zt[nz-1] - zt[0]) / (1.0*(nz-1));
  
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
#endif
  
#if 1
  dx = struct_ptr->drho;
  dy = struct_ptr->dtemp;
  dz = struct_ptr->dye;
  dxi = struct_ptr->drhoi;
  dyi = struct_ptr->dtempi;
  dzi = struct_ptr->dyei;
#endif
  
  dxyi = dxi * dyi;
  dxzi = dxi * dzi;
  dyzi = dyi * dzi;

  dxyzi = dxi * dyi * dzi;

  // determine location in table

  ix = 1 + (int)( (x - xt[0] - 1.0e-10) * dxi );
  iy = 1 + (int)( (y - yt[0] - 1.0e-10) * dyi );
  iz = 1 + (int)( (z - zt[0] - 1.0e-10) * dzi );

  ix = MAX( 1, MIN( ix, nx-2 ) );
  iy = MAX( 1, MIN( iy, ny-2 ) );
  iz = MAX( 1, MIN( iz, nz-2 ) );

  // set up aux vars for interpolation

  delx = xt[ix] - x;
  dely = yt[iy] - y;
  delz = zt[iz] - z;

  int idx[11];

  idx[0] = NSUBTABLES*(ix + nx*(iy + ny*iz));
  idx[1] = NSUBTABLES*((ix-1) + nx*(iy + ny*iz));
  idx[2] = NSUBTABLES*(ix + nx*((iy-1) + ny*iz));
  idx[3] = NSUBTABLES*(ix + nx*(iy + ny*(iz-1)));
  idx[4] = NSUBTABLES*((ix-1) + nx*((iy-1) + ny*iz));
  idx[5] = NSUBTABLES*((ix-1) + nx*(iy + ny*(iz-1)));
  idx[6] = NSUBTABLES*(ix + nx*((iy-1) + ny*(iz-1)));
  idx[7] = NSUBTABLES*((ix-1) + nx*((iy-1) + ny*(iz-1)));
  idx[8] = NSUBTABLES*((ix+1) + nx*(iy + ny*iz));
  idx[9] = NSUBTABLES*(ix + nx*((iy+1) + ny*iz));
  idx[10] = NSUBTABLES*(ix + nx*(iy + ny*(iz+1)));

  int iv;
  
  for(iv=0;iv<nvars;iv++) {

    // set up aux vars for interpolation
    // assuming array ordering (iv, ix, iy, iz)

    fh[0][ivs[iv]] = ft[ivs[iv]+idx[0]];
    fh[1][ivs[iv]] = ft[ivs[iv]+idx[1]];
    fh[2][ivs[iv]] = ft[ivs[iv]+idx[2]];
    fh[3][ivs[iv]] = ft[ivs[iv]+idx[3]];
    fh[4][ivs[iv]] = ft[ivs[iv]+idx[4]];
    fh[5][ivs[iv]] = ft[ivs[iv]+idx[5]];
    fh[6][ivs[iv]] = ft[ivs[iv]+idx[6]];
    fh[7][ivs[iv]] = ft[ivs[iv]+idx[7]];
    fh[8][ivs[iv]] = ft[ivs[iv]+idx[8]];
    fh[9][ivs[iv]] = ft[ivs[iv]+idx[9]];
    fh[10][ivs[iv]] = ft[ivs[iv]+idx[10]];

    // set up coeffs of interpolation polynomical and
    // evaluate function values

    a[0] = fh[0][ivs[iv]];
    a[1] = dxi *   ( fh[1][ivs[iv]] - fh[0][ivs[iv]] );
    a[2] = dyi *   ( fh[2][ivs[iv]] - fh[0][ivs[iv]] );
    a[3] = dzi *   ( fh[3][ivs[iv]] - fh[0][ivs[iv]] );
    a[4] = dxyi *  ( fh[4][ivs[iv]] - fh[1][ivs[iv]] - fh[2][ivs[iv]] + fh[0][ivs[iv]] );
    a[5] = dxzi *  ( fh[5][ivs[iv]] - fh[1][ivs[iv]] - fh[3][ivs[iv]] + fh[0][ivs[iv]] );
    a[6] = dyzi *  ( fh[6][ivs[iv]] - fh[2][ivs[iv]] - fh[3][ivs[iv]] + fh[0][ivs[iv]] );
    a[7] = dxyzi * ( fh[7][ivs[iv]] - fh[0][ivs[iv]] + fh[1][ivs[iv]] + fh[2][ivs[iv]] + 
		     fh[3][ivs[iv]] - fh[4][ivs[iv]] - fh[5][ivs[iv]] - fh[6][ivs[iv]] );

//    a[8] = 0.5 * dxi *   ( fh[1][ivs[iv]] - fh[8][ivs[iv]] );
//    a[9] = 0.5 * dyi *   ( fh[2][ivs[iv]] - fh[9][ivs[iv]] );
    a[8] = 0.5 * dxi *   ( fh[1][ivs[iv]] - fh[8][ivs[iv]] );
    a[9] = 0.5 * dyi *   ( fh[2][ivs[iv]] - fh[9][ivs[iv]] );

    f[iv] = a[0] + a[1] * delx
      + a[2] * dely
      + a[3] * delz
      + a[4] * delx * dely
      + a[5] * delx * delz
      + a[6] * dely * delz
      + a[7] * delx * dely * delz;

    if (iv==0) {
    *dlPdlt = -a[9];
    *dlPdlrho = -a[8];
    }
    if (iv==1) {
    *dlepsdlt = -a[9];
    *dlepsdlrho = -a[8];
    }
  }

  return;
}
  
