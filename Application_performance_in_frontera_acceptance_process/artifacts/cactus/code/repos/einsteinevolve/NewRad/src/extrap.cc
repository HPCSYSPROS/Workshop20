#include <cassert>
#include <cmath>

#include <cctk.h>

#define KRANC_C
#include <GenericFD.h>

using namespace std;



// Adapted from BSSN_MoL's files Init.F
// Erik Schnetter: This code was probably originally written by Miguel
// Alcubierre.
static
void extrap_kernel (cGH const* restrict const cctkGH,
                    int const* restrict const bmin,
                    int const* restrict const bmax,
                    int const* restrict const dir,
                    CCTK_REAL* restrict const var)
{
  int const ni = cctkGH->cctk_lsh[0];
  int const nj = cctkGH->cctk_lsh[1];
  int const nk = cctkGH->cctk_lsh[2];
  
  int const ai = cctkGH->cctk_ash[0];
  int const aj = cctkGH->cctk_ash[1];
  int const ak = cctkGH->cctk_ash[2];
  
  int const si = -dir[0];
  int const sj = -dir[1];
  int const sk = -dir[2];
  
  int const di = 1;
  int const dj = ai;
  int const dk = ai*aj;
  int const np = ai*aj*ak;
  
  int const dind = si*di + sj*dj + sk*dk;
  
  int imin[3], imax[3], idir[3];
  for (int d=0; d<3; ++d) {
    if (dir[d]<0) {
      // lower boundary
      imin[d] = bmax[d]-1;
      imax[d] = bmin[d]-1;
      idir[d] = -1;
    } else {
      // interior and upper boundary
      imin[d] = bmin[d];
      imax[d] = bmax[d];
      idir[d] = +1;
    }
  }
  
  // Note: These loops are not parallel, since each iteration may
  // access grid points set by previous iterations.
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
        
        // Apply boundary conditions to get e.g. Gammas on physical
        // boundaries.  Notice that we only want to apply boundary
        // conditions to the Gammas, all other variables have been
        // defined point-wise all the way to the boundaries from the
        // original ADM quantities (which SHOULD be correctly defined
        // all the way to the boundaries by the initial data
        // routines).  Here I use cubic extrapolation (though rational
        // extrapolation might be better).
        
        assert (ind       >=0 and ind       <np);
        assert (ind+4*dind>=0 and ind+4*dind<np);
        var[ind] = (4*var[ind+dind] - 6*var[ind+2*dind] + 4*var[ind+3*dind]
                    - var[ind+4*dind]);
        
      } // for i j k
    }
  }
}



// Adapted from Kranc's KrancNumericalTools/GenericFD's file
// GenericFD.c
void newrad_extrap_loop (cGH const* restrict const cctkGH,
			 CCTK_REAL* restrict const var)
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
              extrap_kernel (cctkGH, bmin, bmax, dir, var);
            }
            
          }
        } // for dir0 dir1 dir2
      }
    }
  }
}



extern "C"
CCTK_INT ExtrapolateGammas1 (CCTK_POINTER_TO_CONST const cctkGH_,
                             CCTK_REAL* restrict const var)
{
  cGH const* restrict const cctkGH = static_cast<cGH const*> (cctkGH_);
  if (not cctkGH) {
    CCTK_WARN (CCTK_WARN_ABORT,
               "cctkGH is NULL");
  }
  
#if 0
  CCTK_REAL* restrict const var = CCTK_VarDataPtr (cctkGH, 0, varname);
  if (not var) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot access variable \"%s\"", varname);
  }
#endif
  
  if (not var) {
    CCTK_WARN (CCTK_WARN_ABORT,
               "Pointer to variable is NULL");
  }
  
  newrad_extrap_loop (cctkGH, var);
  
  return 0;
}
