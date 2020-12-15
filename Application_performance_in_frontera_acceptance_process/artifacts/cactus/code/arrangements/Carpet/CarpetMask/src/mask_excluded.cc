#include <cmath>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "mask_excluded.hh"

#include <CarpetReduce_bits.h>

namespace CarpetMask {

using namespace std;

/**
 * Set the weight in the excluded regions to zero.
 */

void CarpetExcludedSetup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT *restrict const iweight = static_cast<CCTK_INT *>(
      CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::iweight"));

  CCTK_REAL *restrict const excised_cells = static_cast<CCTK_REAL *>(
      CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::excised_cells"));

  if (not iweight or not excised_cells) {
    CCTK_WARN(CCTK_WARN_ABORT, "CarpetReduce is not active, or "
                               "CarpetReduce::iweight does not have storage");
  }

  // Volume of a grid cell on this level, in terms of coarse grid
  // cells
  CCTK_REAL const cell_volume =
      1.0 / (cctk_levfac[0] * cctk_levfac[1] * cctk_levfac[2]);

  unsigned const bits = BMSK(cctk_dim);
  CCTK_REAL const factor = 1.0 / bits;

  for (int n = 0; n < 10; ++n) {

    CCTK_REAL const r0 = excluded_radius[n];
    if (r0 >= 0.0) {

      CCTK_REAL const x0 = excluded_centre_x[n];
      CCTK_REAL const y0 = excluded_centre_y[n];
      CCTK_REAL const z0 = excluded_centre_z[n];

      CCTK_REAL const r2 = pow(r0, 2);

      bool const exterior = exclude_exterior[n];

      CCTK_REAL local_excised = 0.0;
#pragma omp parallel reduction(+ : local_excised)
      CCTK_LOOP3_ALL(CarpetExcludedSetup, cctkGH, i, j, k) {
        int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL const dx2 = pow(x[ind] - x0, 2);
        CCTK_REAL const dy2 = pow(y[ind] - y0, 2);
        CCTK_REAL const dz2 = pow(z[ind] - z0, 2);

        if (exterior ? dx2 + dy2 + dz2 >= r2 : dx2 + dy2 + dz2 <= r2) {
          // Tally up the weight we are removing
          local_excised += cell_volume * factor * BCNT(iweight[ind]);
          iweight[ind] = 0;
        }
      }
      CCTK_ENDLOOP3_ALL(CarpetExcludedSetup);
      *excised_cells += local_excised;

    } // if r>=0
  }   // for n
}

} // namespace CarpetMask
