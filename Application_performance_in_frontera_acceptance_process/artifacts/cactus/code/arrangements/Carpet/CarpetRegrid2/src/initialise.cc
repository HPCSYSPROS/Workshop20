#include <cassert>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "indexing.hh"

namespace CarpetRegrid2 {

extern "C" {
void CarpetRegrid2_Initialise(CCTK_ARGUMENTS);
}

void CarpetRegrid2_Initialise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Initialise meta-information
  *last_iteration = -1;
  *last_map = -1;

  // Initialise refinement information
  for (int n = 0; n < 10; ++n) {
    num_levels[n] = 0;
  }

  int lsh[2];
  getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);

#define INIT_CENTRE(N)                                                         \
  do {                                                                         \
    if (num_centres >= N) {                                                    \
      num_levels[N - 1] = num_levels_##N;                                      \
      active[N - 1] = active_##N;                                              \
      position_x[N - 1] = position_x_##N;                                      \
      position_y[N - 1] = position_y_##N;                                      \
      position_z[N - 1] = position_z_##N;                                      \
      for (int rl = 0; rl < 30; ++rl) {                                        \
        int const ind = index2(lsh, rl, N - 1);                                \
        radius[ind] = radius_##N[rl];                                          \
        radius_x[ind] = radius_x_##N[rl];                                      \
        radius_y[ind] = radius_y_##N[rl];                                      \
        radius_z[ind] = radius_z_##N[rl];                                      \
      }                                                                        \
      old_position_x[N - 1] = position_x[N - 1];                               \
      old_position_y[N - 1] = position_y[N - 1];                               \
      old_position_z[N - 1] = position_z[N - 1];                               \
    }                                                                          \
  } while (0)

  INIT_CENTRE(1);
  INIT_CENTRE(2);
  INIT_CENTRE(3);
  INIT_CENTRE(4);
  INIT_CENTRE(5);
  INIT_CENTRE(6);
  INIT_CENTRE(7);
  INIT_CENTRE(8);
  INIT_CENTRE(9);
  INIT_CENTRE(10);

#undef INIT_CENTRE

  if (verbose or veryverbose) {
    for (int n = 0; n < num_centres; ++n) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Initialising position of centre %d to [%g,%g,%g]", n,
                 static_cast<double>(position_x[n]),
                 static_cast<double>(position_y[n]),
                 static_cast<double>(position_z[n]));
    }
  }
}

} // namespace CarpetRegrid2
