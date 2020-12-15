#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "bits.h"

extern "C" void MaskBase_SetMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    int const reflevel = GetRefinementLevel(cctkGH);
    CCTK_VInfo(CCTK_THORNSTRING, "Finalise the weight on level %d", reflevel);
  }

  unsigned const bits = BMSK(cctk_dim);
  CCTK_REAL const factor = 1.0 / bits;
#pragma omp parallel
  CCTK_LOOP3_ALL(MaskBase_SetMask, cctkGH, i, j, k) {
    int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
    weight[ind] = factor * BCNT(iweight[ind]);
    one[ind] = 1.0;
  }
  CCTK_ENDLOOP3_ALL(MaskBase_SetMask);
}
