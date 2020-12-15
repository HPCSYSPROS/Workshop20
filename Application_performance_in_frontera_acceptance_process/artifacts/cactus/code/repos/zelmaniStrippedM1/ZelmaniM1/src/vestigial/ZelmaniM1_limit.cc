//  Adapted from Charon: full-Boltzmann transport code for Cactus
//  Copyright (C) 2012, David Radice
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ZelmaniM1.hh"
//#include "utils.hh"

namespace {

CCTK_REAL UTILS_SIGN(CCTK_REAL a) {
  if (a<0.0) return -1.0;
  return 1.0;
}

CCTK_REAL minmod(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c) {
    return UTILS_SIGN(a)*static_cast<CCTK_REAL>(
            std::abs(UTILS_SIGN(a) + UTILS_SIGN(b) + UTILS_SIGN(c)) == 3)*
                std::min(std::abs(a), std::min(std::abs(b), std::abs(c)));
}

CCTK_REAL minmod2(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c) {
    CCTK_REAL const sa = UTILS_SIGN(a);
    CCTK_REAL const sb = UTILS_SIGN(b);
    CCTK_REAL const sc = UTILS_SIGN(c);

    CCTK_REAL const aa = std::abs(a);
    CCTK_REAL const ab = std::abs(b);
    CCTK_REAL const ac = std::abs(c);

    CCTK_REAL const s = sa*static_cast<CCTK_REAL>(std::abs(sa + sb + sc) == 3);
    CCTK_REAL const mabac = std::min(ab, ac);
    if(aa < 2.0*mabac) {
        return s*aa;
    }
    else {
        return s*mabac;
    }
}

CCTK_REAL none(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c) {
    return a;
}

CCTK_REAL step(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c) {
    return 0;
}

}

extern "C" void zm1_Limit(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(zm1_verbose) {
        CCTK_INFO("Charon_Limit");
    }

    // Select the slope limiter
    CCTK_REAL (*phi)(CCTK_REAL, CCTK_REAL, CCTK_REAL) = minmod2;

    //if(CCTK_Equals(zm1_DG_limiter, "step")) {
    //    phi = step;
    //}
    //else if(CCTK_Equals(zm1_DG_limiter, "minmod")) {
    //    phi = minmod;
    //}
    //else if(CCTK_Equals(zm1_DG_limiter, "none")) {
    //    phi = none;
    //}

    CCTK_REAL const idelta[3] = {
        cctk_levfac[0] / cctk_delta_space[0],
        cctk_levfac[1] / cctk_delta_space[1],
        cctk_levfac[2] / cctk_delta_space[2]
    };

    // Number of elements
    int const nelem[3] = {
        cctk_lsh[0] / 2,
        cctk_lsh[1] / 2,
        cctk_lsh[2] / 2
    };

    CCTK_REAL * avg = new CCTK_REAL[nelem[0]*nelem[1]*nelem[2]];
    // This is used to access avg
    int const stride[3] = {
        1,
        nelem[0],
        nelem[0]*nelem[1]
    };

#pragma omp parallel
    {
        for(int ig = 0; ig < ngroups*nspecies; ++ig) 
	  for (int ivar = 0; ivar < 4; ++ivar) {
	    
	    // Select the variable to limit
	    CCTK_REAL * u;
	    if ((*zm1_RKstep == 1) && do_m1_RK2) {
	      if (ivar == 0) u = enuh;
	      if (ivar == 1) u = fnuxh;
	      if (ivar == 2) u = fnuyh;
	      if (ivar == 3) u = fnuzh;
            } else {
	      if (ivar == 0) u = enu;
	      if (ivar == 1) u = fnux;
	      if (ivar == 2) u = fnuy;
	      if (ivar == 3) u = fnuz;
	    }

            // Step 1: compute element averages
#pragma omp for collapse(3)
            for(int ke = 0; ke < nelem[2]; ++ke)
            for(int je = 0; je < nelem[1]; ++je)
            for(int ie = 0; ie < nelem[0]; ++ie) {
                int const eidx = ke*stride[2] + je*stride[1] + ie*stride[0];
                avg[eidx] = 0;
                for(int kl = 0; kl < 2; ++kl)
                for(int jl = 0; jl < 2; ++jl)
                for(int il = 0; il < 2; ++il) {
                    int const i = il + ie*2;
                    int const j = jl + je*2;
                    int const k = kl + ke*2;

                    int const ijk = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, ig);
                    avg[eidx] += u[ijk];
                }
                avg[eidx] = 0.125*avg[eidx];
            }

            // Step 2: do the actual limiting
#pragma omp for collapse(3)
            for(int ke = 1; ke < nelem[2] - 1; ++ke)
            for(int je = 1; je < nelem[1] - 1; ++je)
            for(int ie = 1; ie < nelem[0] - 1; ++ie) {
                // Neighbouring elements
                int const eidx = ke*stride[2] + je*stride[1] + ie*stride[0];
                int const eidxm[3] = {
                    ke*stride[2] + je*stride[1] + (ie - 1)*stride[0],
                    ke*stride[2] + (je - 1)*stride[1] + ie*stride[0],
                    (ke - 1)*stride[2] + je*stride[1] + ie*stride[0],
                };
                int const eidxp[3] = {
                    ke*stride[2] + je*stride[1] + (ie + 1)*stride[0],
                    ke*stride[2] + (je + 1)*stride[1] + ie*stride[0],
                    (ke + 1)*stride[2] + je*stride[1] + ie*stride[0],
                };

                // First point of the element
                int const i = ie*2;
                int const j = je*2;
                int const k = ke*2;

                // Compute interior slopes
                CCTK_REAL islope[3];
                islope[0] = 0.25*idelta[0]*(
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j, k, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j+1, k, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j+1, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j, k+1, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j+1, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j+1, k+1, ig)]));
                islope[1] = 0.25*idelta[1]*(
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i, j+1, k, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j+1, k, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i, j+1, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j, k+1, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j+1, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j, k+1, ig)]));
                islope[2] = 0.25*idelta[2]*(
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i, j, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i, j+1, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i, j+1, k, ig)]) +
                        (u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j+1, k+1, ig)] -
                                u[CCTK_VECTGFINDEX3D(cctkGH, i+1, j+1, k, ig)]));

                // Compute exterior slopes
                CCTK_REAL eslope[3];
                for(int dir = 0; dir < 3; ++dir) {
                    // D_-, D, D_+
                    CCTK_REAL const Dm = 0.5*idelta[dir]*(avg[eidx] -
                            avg[eidxm[dir]]);
                    CCTK_REAL const Dp = 0.5*idelta[dir]*(avg[eidxp[dir]] -
                            avg[eidx]);

                    eslope[dir] = phi(islope[dir], Dm, Dp);
                
		}

                // Check if limiting is necessary
                bool limit = false;
		double ZM1_LIMIT_THRESHOLD = 0.e-10;
                for(int dir = 0; dir < 3; ++dir) {
                    // This makes the code more robust with respect to
                    // floating point errors
                    limit |= (std::abs(islope[dir]) >
                              std::abs(eslope[dir]) + ZM1_LIMIT_THRESHOLD);
                }

                // Reconstruction
                if(limit) {
                    int const ijk  = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    int const ijkp = CCTK_GFINDEX3D(cctkGH, i+1, j+1, k+1);
                    CCTK_REAL const xm = 0.5*(x[ijk] + x[ijkp]);
                    CCTK_REAL const ym = 0.5*(y[ijk] + y[ijkp]);
                    CCTK_REAL const zm = 0.5*(z[ijk] + z[ijkp]);

                    for(int kl = 0; kl < 2; ++kl)
                    for(int jl = 0; jl < 2; ++jl)
                    for(int il = 0; il < 2; ++il) {
                        int const ii = i + il;
                        int const jj = j + jl;
                        int const kk = k + kl;
                        int const iijjkk = CCTK_VECTGFINDEX3D(cctkGH, ii, jj, kk, ig);

                        u[iijjkk] = avg[eidx] + eslope[0]*(x[iijjkk] - xm) +
                            eslope[1]*(y[iijjkk] - ym) +
                            eslope[2]*(z[iijjkk] - zm);
                    }
                }
            }
        }
    }

    delete[] avg;
}
