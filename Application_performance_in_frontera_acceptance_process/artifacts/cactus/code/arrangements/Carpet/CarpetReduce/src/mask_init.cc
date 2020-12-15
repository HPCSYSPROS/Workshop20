#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "bits.h"

extern "C" void MaskBase_InitMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Initialise the weight to 1 everywhere */
  if (verbose) {
    int const reflevel = GetRefinementLevel(cctkGH);
    CCTK_VInfo(CCTK_THORNSTRING, "Initialising weight to 1 on level %d",
               reflevel);
  }

  unsigned const bits = BMSK(cctk_dim);
  unsigned const allbits = BMSK(bits) - 1;
#pragma omp parallel
  CCTK_LOOP3_ALL(MaskBase_InitMask, cctkGH, i, j, k) {
    int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
    iweight[ind] = allbits;
  }
  CCTK_ENDLOOP3_ALL(MaskBase_InitMask);
}
