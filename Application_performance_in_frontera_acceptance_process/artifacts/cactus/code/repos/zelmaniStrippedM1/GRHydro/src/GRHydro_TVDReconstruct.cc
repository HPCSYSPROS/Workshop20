#include <cmath>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_TVDReconstruct.hh"
#include "GRHydro_Reconstruct_drv_impl.hh"


using namespace std;

/**
   TVD reconstruction operator.
   limiter_type: 0 = minmod
                 1 = mc2
                 2 = superbee
*/
template <tvd::limiter_type limiter>
CCTK_REAL GRHydro_TVDReconstruct1d<limiter>::
minmod_func(const CCTK_REAL dupw, const CCTK_REAL dloc)
{
  return 0.5 * (copysign(1.0,dupw) + copysign(1.0,dloc)) *
    fmin(fabs(dupw),fabs(dloc));
}

template <tvd::limiter_type limiter>
CCTK_REAL GRHydro_TVDReconstruct1d<limiter>::
min3(const CCTK_REAL a, const CCTK_REAL b, const CCTK_REAL c)
{
  return fmin(a,fmin(b,c));
}

template <tvd::limiter_type limiter>
template <int dir>
void GRHydro_TVDReconstruct1d<limiter>::
apply(const int nx, const CCTK_REAL* const restrict a,
      CCTK_REAL* const restrict aminus, CCTK_REAL* const restrict aplus,
      const cGH* const cctkGH, const int j, const int k)
{
  DECLARE_CCTK_PARAMETERS;

#define A(i_) (a[ijk[i_]])
#define Aplus(i_) (aplus[ijk[i_]])
#define Aminus(i_) (aminus[ijk[i_]])


  for (int i=GRHydro_stencil-1; i < nx-GRHydro_stencil + 1; ++i)
   {

      const int ijk[3] = {
	dir ==0 ? CCTK_GFINDEX3D(cctkGH, i-1, j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i-1, k) : CCTK_GFINDEX3D(cctkGH, j, k, i-1),
	dir ==0 ? CCTK_GFINDEX3D(cctkGH, i  , j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i  , k) : CCTK_GFINDEX3D(cctkGH, j, k, i  ),
	dir ==0 ? CCTK_GFINDEX3D(cctkGH, i+1, j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i+1, k) : CCTK_GFINDEX3D(cctkGH, j, k, i+1),
      };

      const CCTK_REAL dupw = A(1)-A(0);
      const CCTK_REAL dloc = A(2)-A(1);
      CCTK_REAL delta = 0.0;

      static_assert(limiter >= tvd::minmod && limiter < tvd::numlimiters,
                    "Unknown limiter");
      if(limiter == tvd::minmod) {
	// minmod
	delta = minmod_func(dupw,dloc);
      } else if (limiter == tvd::mc2) {
	// Van Leer MC slope limiter (mc2)
	if(dupw*dloc < 0.0)
	  delta = 0.0;
	else
	  delta = copysign(min3(2.0*fabs(dupw),2.0*fabs(dloc),
				0.5*(fabs(dupw)+fabs(dloc))),dupw+dloc);

      } else if (limiter == tvd::superbee) {
	// superbee
	if(dupw*dloc < 0.0)
	  delta = 0.0;
	else
	  delta = copysign(fmax(fmin(2.0*fabs(dupw),fabs(dloc)),
				fmin(2.0*fabs(dloc),fabs(dupw))),dupw+dloc);

      }
      Aminus(1) = A(1) - 0.5*delta;
      Aplus(1) = A(1) + 0.5*delta;
   }
}

// instantiate all copies we need, this way different operators can be compiled
// in parallel. This must match the select routine in GRHydro_Reconstruct.cc
template class GRHydro_TVDReconstruct1d<tvd::minmod>;
template class GRHydro_TVDReconstruct1d<tvd::mc2>;
template class GRHydro_TVDReconstruct1d<tvd::superbee>;

INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_TVDReconstruct1d<tvd::minmod>)
INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_TVDReconstruct1d<tvd::mc2>)
INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_TVDReconstruct1d<tvd::superbee>)

