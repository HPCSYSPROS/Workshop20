#include <cassert>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

namespace Carpet {

using namespace std;

/** Ensure that the parameters have legal values.
 *
 * Note that this checking happens only after most of Carpet has
 * already been set up.
 */
void CarpetParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if ((CCTK_ParameterQueryTimesSet("periodic", "Carpet") and periodic) or
      (CCTK_ParameterQueryTimesSet("periodic_x", "Carpet") and periodic_x) or
      (CCTK_ParameterQueryTimesSet("periodic_y", "Carpet") and periodic_y) or
      (CCTK_ParameterQueryTimesSet("periodic_z", "Carpet") and periodic_z)) {
    CCTK_PARAMWARN("Some of the parameters \"Carpet::periodic*\" have been set "
                   "to \"yes\".  These parameters are there for compatibility "
                   "reasons only and must not be used.");
  }

  if (adaptive_stepsize and max_refinement_levels > 1) {
    CCTK_PARAMWARN("Adaptive time step sizes do not work with mesh refinement "
                   "yet.  Please use only a single level, and set "
                   "max_refinement_levels=1.");
  }

  if (max_refinement_levels > 1) {

    // InitBase

    enum setup_method_t {
      init_error,
      init_some_levels,
      init_single_level,
      init_two_levels,
      init_all_levels,
    };
    setup_method_t setup_method = init_error;

    if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels")) {
      setup_method = init_some_levels;
    } else if (CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
      setup_method = init_single_level;
    } else if (CCTK_EQUALS(initial_data_setup_method, "init_two_levels")) {
      setup_method = init_two_levels;
    } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {
      setup_method = init_all_levels;
    } else {
      CCTK_PARAMWARN(
          "Unknown value of parameter InitBase::initial_data_setup_method");
    }

    // Carpet

    enum init_method_t {
      error_timelevel,
      each_timelevel,   // Initialise each time level
      fill_timelevels,  // Copy the current to all past timelevels
      three_timelevels, // Carpet's 3 timelevel scheme
      all_timelevels,   // Initial data thorn initialises all timelevels
    };
    init_method_t init_method = error_timelevel;

    if (init_each_timelevel) {
      if (init_fill_timelevels) {
        CCTK_PARAMWARN("Cannot set parameters init_each_timelevel and "
                       "init_fill_timelevels at the same time");
      }
      if (init_3_timelevels) {
        CCTK_PARAMWARN("Cannot set parameters init_each_timelevel and "
                       "init_3_timelevels at the same time");
      }
      init_method = each_timelevel;
    } else if (init_fill_timelevels) {
      if (init_3_timelevels) {
        CCTK_PARAMWARN("Cannot set parameters init_fill_timelevels and "
                       "init_3_timelevels at the same time");
      }
      init_method = fill_timelevels;
    } else if (init_3_timelevels) {
      init_method = three_timelevels;
    } else {
      init_method = all_timelevels;
    }

    switch (init_method) {
    case each_timelevel:
      if (setup_method != init_single_level) {
        CCTK_PARAMWARN(
            "When you set Carpet::init_each_timelevel=yes, then you must also "
            "use InitBase::initial_data_setup_method=\"init_single_level\"");
      }
      break;
    case fill_timelevels:
      // Everything is allowed
      break;
    case three_timelevels:
      // Everything is allowed
      break;
    case all_timelevels:
      if (setup_method != init_all_levels) {
        CCTK_PARAMWARN("When you set neither Carpet::init_each_timelevel=yes "
                       "nor Carpet::init_fill_timelevels=yes nor "
                       "Carpet::init_3_timelevels=yes, then you must also use "
                       "InitBase::initial_data_setup_method=\"init_all_"
                       "levels\"");
      }
      break;
    default:
      assert(0);
    }

  } // if max_refinement_levels > 0
}

} // namespace Carpet
