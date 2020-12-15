#include <limits>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

namespace Carpet {

using namespace std;

void CarpetRefineTimeStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Find the smallest CFL factor for all refinement levels
  CCTK_REAL min_cfl = numeric_limits<CCTK_REAL>::max();
  for (int rl = 0; rl < reflevels; ++rl) {
    CCTK_REAL const level_cfl =
        timereffacts.AT(rl) / maxval(spacereffacts.AT(rl));
    min_cfl = min(min_cfl, level_cfl);
  }

  // Reduce the time step correspondingly
  cctkGH->cctk_delta_time *= min_cfl;
  CCTK_VInfo(CCTK_THORNSTRING, "Timestep reduced by a factor of %g to %g",
             double(1 / min_cfl),
             double(cctkGH->cctk_delta_time / cctkGH->cctk_timefac));
}

} // namespace Carpet
