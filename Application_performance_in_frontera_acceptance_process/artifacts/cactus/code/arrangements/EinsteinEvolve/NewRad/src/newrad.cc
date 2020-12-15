#include <cassert>
#include <cmath>

#include <cctk.h>

#define KRANC_C
#include <GenericFD.h>

using namespace std;



// Adapted from BSSN_MoL's files NewRad.F and newrad.h
// Erik Schnetter: This code was probably originally written by Miguel
// Alcubierre.
static
void newrad_kernel (cGH const* restrict const cctkGH,
                    int const* restrict const bmin,
                    int const* restrict const bmax,
                    int const* restrict const dir,
                    CCTK_REAL const* restrict const var,
                    CCTK_REAL      * restrict const rhs,
                    CCTK_REAL const* restrict const x,
                    CCTK_REAL const* restrict const y,
                    CCTK_REAL const* restrict const z,
                    CCTK_REAL const* restrict const r,
                    CCTK_REAL const& var0,
                    CCTK_REAL const& v0,
                    int const radpower)
{
  int const ni = cctkGH->cctk_lsh[0];
  int const nj = cctkGH->cctk_lsh[1];
  int const nk = cctkGH->cctk_lsh[2];
  
  int const ai = cctkGH->cctk_ash[0];
  int const aj = cctkGH->cctk_ash[1];
  
  int const si = dir[0];
  int const sj = dir[1];
  int const sk = dir[2];
  
  int const di = 1;
  int const dj = ai;
  int const dk = ai*aj;
  
  CCTK_REAL const dx = cctkGH->cctk_delta_space[0] / cctkGH->cctk_levfac[0];
  CCTK_REAL const dy = cctkGH->cctk_delta_space[1] / cctkGH->cctk_levfac[1];
  CCTK_REAL const dz = cctkGH->cctk_delta_space[2] / cctkGH->cctk_levfac[2];
  CCTK_REAL const idx = 1.0/dx;
  CCTK_REAL const idy = 1.0/dy;
  CCTK_REAL const idz = 1.0/dz;
  
  int imin[3], imax[3], idir[3];
  for (int d=0; d<3; ++d) {
    if (dir[d]<0) {
      // lower boundary
      assert (bmin[d] >= 0);
      assert (bmax[d] + 2 <= cctkGH->cctk_lsh[d]);
      imin[d] = bmax[d]-1;
      imax[d] = bmin[d]-1;
      idir[d] = -1;
    } else if (dir[d]>0) {
      // upper boundary
      assert (bmin[d] - 2 >= 0);
      assert (bmax[d] <= cctkGH->cctk_lsh[d]);
      imin[d] = bmin[d];
      imax[d] = bmax[d];
      idir[d] = +1;
    } else {
      // interior
      assert (bmin[d] - 1 >= 0);
      assert (bmax[d] + 1 <= cctkGH->cctk_lsh[d]);
      imin[d] = bmin[d];
      imax[d] = bmax[d];
      idir[d] = +1;
    }
  }
  
  // Warning: these loops are not parallel, since previously
  // calculated RHS are accessed for radpower>=0
  for (int k=imin[2]; k!=imax[2]; k+=idir[2]) {
    for (int j=imin[1]; j!=imax[1]; j+=idir[1]) {
      for (int i=imin[0]; i!=imax[0]; i+=idir[0]) {
        int const ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
        
        // Test looping directions
        if (i==0) assert (idir[0]<0);
        if (j==0) assert (idir[1]<0);
        if (k==0) assert (idir[2]<0);
        if (i==ni-1) assert (idir[0]>0);
        if (j==nj-1) assert (idir[1]>0);
        if (k==nk-1) assert (idir[2]>0);
        
        if (si==0) {
          assert (i-1>=0 and i+1<ni);
        } else {
          assert (i-2*si>=0 and i-2*si<ni);
        }
        if (sj==0) {
          assert (j-1>=0 and j+1<nj);
        } else {
          assert (j-2*sj>=0 and j-2*sj<nj);
        }
        if (sk==0) {
          assert (k-1>=0 and k+1<nk);
        } else {
          assert (k-2*sk>=0 and k-2*sk<nk);
        }
        
        {
          // The main part of the boundary condition assumes that we
          // have an outgoing radial wave with some speed v0:
          //
          //    var  =  var0 + u(r-v0*t)/r
          //
          // This implies the following differential equation:
          //
          //    d_t var  =  - v^i d_i var  -  v0 (var - var0) / r
          //
          // where  vi = v0 xi/r
          
          if (si==0) {
            assert (i-1>=0 and i+1<ni);
          } else {
            assert (i-2*si>=0 and i-2*si<ni);
          }
          if (sj==0) {
            assert (j-1>=0 and j+1<nj);
          } else {
            assert (j-2*sj>=0 and j-2*sj<nj);
          }
          if (sk==0) {
            assert (k-1>=0 and k+1<nk);
          } else {
            assert (k-2*sk>=0 and k-2*sk<nk);
          }
          
          // Find local wave speeds
          CCTK_REAL const rp = r[ind];
          CCTK_REAL const rpi = 1.0/rp;
          
          CCTK_REAL const vx = v0*x[ind]*rpi;
          CCTK_REAL const vy = v0*y[ind]*rpi;
          CCTK_REAL const vz = v0*z[ind]*rpi;
          
          // Find x derivative
          CCTK_REAL derivx;
          if (si==0) {
            derivx = 0.5*(var[ind+di]-var[ind-di])*idx;
          } else {
            derivx = si*0.5*(3*var[ind] - 4*var[ind-si*di]
                             + var[ind-2*si*di])*idx;
          }
          
          // Find y derivative
          CCTK_REAL derivy;
          if (sj==0) {
            derivy = 0.5*(var[ind+dj]-var[ind-dj])*idy;
          } else {
            derivy = sj*0.5*(3*var[ind] - 4*var[ind-sj*dj]
                             + var[ind-2*sj*dj])*idy;
          }
          
          // Find z derivative
          CCTK_REAL derivz;
          if (sk==0) {
            derivz = 0.5*(var[ind+dk]-var[ind-dk])*idz;
          } else {
            derivz = sk*0.5*(3*var[ind] - 4*var[ind-sk*dk]
                             + var[ind-2*sk*dk])*idz;
          }
          
          // Calculate source term
          rhs[ind] =
            - vx*derivx - vy*derivy - vz*derivz - v0*(var[ind] - var0)*rpi;
          
        }
        
        if (radpower >= 0) {
          // *****************************************
          // ***   EXTRAPOLATION OF MISSING PART   ***
          // *****************************************
          //
          // Here we try to extrapolate for the part of the boundary
          // that does not behave as a pure wave (i.e. Coulomb type
          // terms caused by infall of the coordinate lines).
          //
          // This we do by comparing the source term one grid point
          // away from the boundary (which we already have), to what
          // we would have obtained if we had used the boundary
          // condition there.  The difference gives us an idea of the
          // missing part and we extrapolate that to the boundary
          // assuming a power-law decay.
          
          int const ip = i-si;
          int const jp = j-sj;
          int const kp = k-sk;
          assert (ip>=0 and ip<ni);
          assert (jp>=0 and jp<nj);
          assert (kp>=0 and kp<nk);
          
          if (si==0) {
            assert (ip-1>=0 and ip+1<ni);
          } else {
            assert (ip-2*si>=0 and ip-2*si<ni);
          }
          if (sj==0) {
            assert (jp-1>=0 and jp+1<nj);
          } else {
            assert (jp-2*sj>=0 and jp-2*sj<nj);
          }
          if (sk==0) {
            assert (kp-1>=0 and kp+1<nk);
          } else {
            assert (kp-2*sk>=0 and kp-2*sk<nk);
          }
          
          // Find local wave speeds
          int const indp = CCTK_GFINDEX3D(cctkGH, ip,jp,kp);
          
          CCTK_REAL const rp = r[indp];
          CCTK_REAL const rpi = 1.0/rp;
    
          CCTK_REAL const vx = v0*x[indp]*rpi;
          CCTK_REAL const vy = v0*y[indp]*rpi;
          CCTK_REAL const vz = v0*z[indp]*rpi;
          
          // Find x derivative
          CCTK_REAL derivx;
          if (si==0) {
            derivx = 0.5*(var[indp+di]-var[indp-di])*idx;
          } else {
            derivx = si*0.5*(3*var[indp] - 4*var[indp-si*di]
                             + var[indp-2*si*di])*idx;
          }
          
          // Find y derivative
          CCTK_REAL derivy;
          if (sj==0) {
            derivy = 0.5*(var[indp+dj]-var[indp-dj])*idy;
          } else {
            derivy = sj*0.5*(3*var[indp] - 4*var[indp-sj*dj]
                             + var[indp-2*sj*dj])*idy;
          }
          
          // Find z derivative
          CCTK_REAL derivz;
          if (sk==0) {
            derivz = 0.5*(var[indp+dk]-var[indp-dk])*idz;
          } else {
            derivz = sk*0.5*(3*var[indp] - 4*var[indp-sk*dk]
                             + var[indp-2*sk*dk])*idz;
          }
          
          // Find difference in sources
          CCTK_REAL const aux =
            rhs[indp] +
            vx*derivx + vy*derivy + vz*derivz + v0*(var[indp] - var0)*rpi;
          
          // Extrapolate difference and add it to source in boundary
          rhs[ind] += aux*pow(rp/r[ind],radpower);
          
        } // if radpower>=0
        
      } // for i j k
    }
  }
}



// Adapted from Kranc's KrancNumericalTools/GenericFD's file
// GenericFD.c
static
void newrad_loop (cGH const* restrict const cctkGH,
                  CCTK_REAL const* restrict const var,
                  CCTK_REAL      * restrict const rhs,
                  CCTK_REAL const* restrict const x,
                  CCTK_REAL const* restrict const y,
                  CCTK_REAL const* restrict const z,
                  CCTK_REAL const* restrict const r,
                  CCTK_REAL const& var0,
                  CCTK_REAL const& v0,
                  int const radpower)
{
  int imin[3], imax[3], is_symbnd[6], is_physbnd[6], is_ipbnd[6];
  GenericFD_GetBoundaryInfo
    (cctkGH, cctkGH->cctk_ash, cctkGH->cctk_lsh, cctkGH->cctk_bbox,
     cctkGH->cctk_nghostzones, 
     imin, imax, is_symbnd, is_physbnd, is_ipbnd);
  
  // Loop over all faces:
  // Loop over faces first, then corners, and then edges, so that the
  // stencil only sees points that have already been treated.
  // ifec means: interior-face-edge-corner.
  for (int ifec=1; ifec<=3; ++ifec) {
    for (int dir2=-1; dir2<=+1; ++dir2) {
      for (int dir1=-1; dir1<=+1; ++dir1) {
        for (int dir0=-1; dir0<=+1; ++dir0) {
          int const dir[3] = { dir0, dir1, dir2 };
          
          int nnz = 0;
          for (int d=0; d<3; ++d) {
            if (dir[d]) ++nnz;
          }
          if (nnz == ifec) {
            
            // one of the faces is a boundary
            bool have_bnd = false;
            // at least one boundary face is a physical boundary
            bool any_physbnd = false;
            // all boundary faces are not inter-processor boundaries
            bool all_not_ipbnd = true;
            
            int bmin[3], bmax[3];
            for (int d=0; d<3; ++d) {
              switch (dir[d]) {
              case -1:
                bmin[d] = 0;
                bmax[d] = imin[d];
                have_bnd = true;
                any_physbnd = any_physbnd or is_physbnd[2*d+0];
                all_not_ipbnd = all_not_ipbnd and not is_ipbnd[2*d+0];
                break;
              case 0:
                bmin[d] = imin[d];
                bmax[d] = imax[d];
                break;
              case +1:
                bmin[d] = imax[d];
                bmax[d] = cctkGH->cctk_lsh[d];
                have_bnd = true;
                any_physbnd = any_physbnd or is_physbnd[2*d+1];
                all_not_ipbnd = all_not_ipbnd and not is_ipbnd[2*d+1];
                break;
              }
            }
            assert (have_bnd);  // must be true since nnz>0
            
            if (have_bnd and any_physbnd and all_not_ipbnd) {
              newrad_kernel (cctkGH, bmin, bmax, dir,
                             var, rhs, x,y,z,r, var0, v0, radpower);
            }
            
          }
        } // for dir0 dir1 dir2
      }
    }
  }
}



extern "C"
CCTK_INT NewRad_Apply1 (CCTK_POINTER_TO_CONST const cctkGH_,
                        CCTK_REAL const* restrict const var,
                        CCTK_REAL      * restrict const rhs,
                        CCTK_REAL const var0,
                        CCTK_REAL const v0,
                        CCTK_INT const radpower)
{
  cGH const* restrict const cctkGH = static_cast<cGH const*> (cctkGH_);
  if (not cctkGH) {
    CCTK_ERROR ("cctkGH is NULL");
  }
  
#if 0
  CCTK_REAL const* restrict const var = CCTK_VarDataPtr (cctkGH, 0, varname);
  if (not var) {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot access variable \"%s\"", varname);
  }
  CCTK_REAL      * restrict const rhs = CCTK_VarDataPtr (cctkGH, 0, rhsname);
  if (not rhs) {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot access RHS variable \"%s\"", rhsname);
  }
#endif
  
  if (not var) {
    CCTK_ERROR ("Pointer to variable is NULL");
  }
  if (not rhs) {
    CCTK_ERROR ("Pointer to RHS is NULL");
  }
  
  CCTK_REAL const* restrict const x =
    static_cast<CCTK_REAL const*> (CCTK_VarDataPtr (cctkGH, 0, "grid::x"));
  CCTK_REAL const* restrict const y =
    static_cast<CCTK_REAL const*> (CCTK_VarDataPtr (cctkGH, 0, "grid::y"));
  CCTK_REAL const* restrict const z =
    static_cast<CCTK_REAL const*> (CCTK_VarDataPtr (cctkGH, 0, "grid::z"));
  CCTK_REAL const* restrict const r =
    static_cast<CCTK_REAL const*> (CCTK_VarDataPtr (cctkGH, 0, "grid::r"));
  if (not x or not y or not z or not z) {
    CCTK_ERROR ("Cannot access coordinate variables x, y, z, and r");
  }
  
  newrad_loop (cctkGH, var, rhs, x,y,z,r, var0, v0, radpower);
  
  return 0;
}
