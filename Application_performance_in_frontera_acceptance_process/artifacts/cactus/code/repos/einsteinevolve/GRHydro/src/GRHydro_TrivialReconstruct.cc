#include <cmath>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_TrivialReconstruct.hh"
#include "GRHydro_Reconstruct_drv_impl.hh"

using namespace std;


/**
   Trivial first-order reconstruction operator.
*/
template <int dir>
void GRHydro_TrivialReconstruct1d::
apply(const int nx, const CCTK_REAL* const restrict a,
      CCTK_REAL* const restrict aminus, CCTK_REAL* const restrict aplus,
      const cGH* const cctkGH, const int j, const int k)
{
   for (int i=0; i < nx; ++i)
   {
      const int ijk = dir ==0 ? CCTK_GFINDEX3D(cctkGH, i  , j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i  , k) : CCTK_GFINDEX3D(cctkGH, j, k, i  );

      // Trivial reconstruction!
      aminus[ijk] = a[ijk];
      aplus[ijk]  = a[ijk];
   }
}

// instantiate all copies we need, this way different operators can be compiled
// in parallel. This must match the select routine in GRHydro_Reconstruct.cc
INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_TrivialReconstruct1d)
