#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// IRIX wants this before <time.h>
#if HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#elif HAVE_TIME_H
#include <time.h>
#endif
#endif

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <defs.hh>
#include <dist.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

// Small number to avoid division by zero
CCTK_REAL const eps = 1.0e-15;

// Return the current wall time
static CCTK_REAL get_walltime() {
#ifdef HAVE_TIME_GETTIMEOFDAY
  // get the current time
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / CCTK_REAL(1.0e6);
#else
  return CCTK_REAL(0.0);
#endif
}

// Calculate the number of updates for the current level
static void current_level_updates(cGH const *const cctkGH,
                                  CCTK_REAL &local_grid_updates,
                                  CCTK_REAL &global_grid_updates,
                                  CCTK_REAL &local_interior_updates,
                                  CCTK_REAL &global_interior_updates) {
  DECLARE_CCTK_PARAMETERS;

  // Count the weighted number of grid points
  // (int is not good enough for this calculation)
  CCTK_REAL local_num_grid_points = 0;
  CCTK_REAL global_num_grid_points = 0;
  CCTK_REAL local_num_interior_points = 0;
  CCTK_REAL global_num_interior_points = 0;
  for (int m = 0; m < maps; ++m) {
    assert(reflevel >= 0);
    int const rl = reflevel;
    for (int c = 0; c < vhh.AT(m)->components(rl); ++c) {
      assert(mglevel >= 0);
      int const ml = mglevel;

      // Regions
      dh::light_dboxes const &light_box =
          vdd.AT(m)->light_boxes.AT(ml).AT(rl).AT(c);
      ibbox const &interior = light_box.interior; // with outer boundary
      ibbox const &owned = light_box.owned;       // without outer boundary

      // Count the grid points
      CCTK_REAL const interior_size = interior.size();
      CCTK_REAL const owned_size = owned.size();

      if (vhh.AT(m)->is_local(rl, c)) {
        local_num_grid_points += interior_size;
        local_num_interior_points += owned_size;
      }
      global_num_grid_points += interior_size;
      global_num_interior_points += owned_size;

    } // for c
  }   // for m

  // Take number of RHS evaluations per time step into account
  static int int_steps = -1;
  if (int_steps == -1) {
    if (num_integrator_substeps != -1) {
      // if the user parameter is set, use it
      int_steps = num_integrator_substeps;
    } else if (CCTK_IsFunctionAliased("MoLNumIntegratorSubsteps")) {
      // if there is an aliased function, use it
      int_steps = MoLNumIntegratorSubsteps();
    } else {
      // use a sensible default, even if it is wrong -- it is better
      // to have timing information which is wrong by a constant
      // factor than to abort the code
      int_steps = 1;
    }
  }
  local_grid_updates = local_num_grid_points * int_steps;
  global_grid_updates = global_num_grid_points * int_steps;
  local_interior_updates = local_num_interior_points * int_steps;
  global_interior_updates = global_num_interior_points * int_steps;
}

// Time at which the simulation started
CCTK_REAL startup_walltime; // in seconds

// Time at which the evolution started
bool in_evolution = false;
CCTK_REAL initial_walltime; // in seconds
CCTK_REAL initial_phystime;

// Last starting time
enum timing_state_t { state_computing, state_communicating, state_io };
timing_state_t timing_state = state_computing;
CCTK_REAL time_start;

// Last starting time for this level
int timing_level = -1;
CCTK_REAL time_level_start;

// Initialise the timing variables (to be called before basegrid)
void InitTimingStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  startup_walltime = get_walltime();

  *physical_time_per_hour = 0.0;
  *current_physical_time_per_hour = 0.0;

  *time_total = 0.0;
  *time_evolution = 0.0;
  *time_computing = 0.0;
  *time_communicating = 0.0;
  *time_io = 0.0;

  *evolution_steps_count = 0.0;

  *local_grid_points_per_second = 0.0;
  *total_grid_points_per_second = 0.0;
  *local_grid_point_updates_count = 0.0;
  *total_grid_point_updates_count = 0.0;
  *local_interior_points_per_second = 0.0;
  *total_interior_points_per_second = 0.0;
  *local_interior_point_updates_count = 0.0;
  *total_interior_point_updates_count = 0.0;

  *io_per_second = 0.0;
  *io_bytes_per_second = 0.0;
  *io_bytes_ascii_per_second = 0.0;
  *io_bytes_binary_per_second = 0.0;
  *io_count = 0.0;
  *io_bytes_count = 0.0;
  *io_bytes_ascii_count = 0.0;
  *io_bytes_binary_count = 0.0;

  *comm_per_second = 0.0;
  *comm_bytes_per_second = 0.0;
  *comm_count = 0.0;
  *comm_bytes_count = 0.0;

  *time_levels = 0.0;

  *grid_points_per_second = 0.0;
  *grid_point_updates_count = 0.0;

  for (int rl = 0; rl < max_refinement_levels; ++rl) {
    time_level[rl] = 0.0;
    time_level_count[rl] = 0.0;
  }
}

// Begin timing (to be called after initialisation, just before the
// main evolution begins)
void BeginTimingEvolution(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  in_evolution = true;
  initial_walltime = get_walltime();
  initial_phystime = cctkGH->cctk_time;
}

// Take a step on the current refinement and multigrid level (to be
// called when EVOL is scheduled)
void StepTimingEvolution(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  assert(in_evolution);
  assert(timing_state == state_computing);

  CCTK_REAL local_grid_updates, global_grid_updates;
  CCTK_REAL local_interior_updates, global_interior_updates;
  current_level_updates(cctkGH, local_grid_updates, global_grid_updates,
                        local_interior_updates, global_interior_updates);

  ++*evolution_steps_count;

  *local_grid_point_updates_count += local_grid_updates;
  *total_grid_point_updates_count += global_grid_updates;
  *local_interior_point_updates_count += local_interior_updates;
  *total_interior_point_updates_count += global_interior_updates;

  *grid_point_updates_count = *local_grid_point_updates_count;
  *interior_point_updates_count = *local_interior_point_updates_count;
}

// Count time spent on individual levels (to be called from Carpet's
// initialisation, evolution, and shutdown drivers)
void BeginTimingLevel(cGH const *const cctkGH) {
  assert(reflevel != -1);
  assert(timing_level == -1);
  timing_level = reflevel;
  time_level_start = get_walltime();
}

void EndTimingLevel(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  assert(reflevel != -1);
  assert(timing_level == reflevel);
  timing_level = -1;
  CCTK_REAL const time_level_end = get_walltime();

  time_level[reflevel] += time_level_end - time_level_start;
  ++time_level_count[reflevel];
  *time_levels += time_level_end - time_level_start;
}

// Count some I/O (to be called from the I/O routine)
void BeginTimingIO(cGH const *const cctkGH) {
  assert(timing_state == state_computing);
  timing_state = state_io;
  time_start = get_walltime();
}

void EndTimingIO(cGH const *const cctkGH, CCTK_REAL const files,
                 CCTK_REAL const bytes, bool const is_binary) {
  DECLARE_CCTK_ARGUMENTS;

  assert(timing_state == state_io);
  timing_state = state_computing;
  CCTK_REAL const time_end = get_walltime();

  *time_io += time_end - time_start;

  *io_count += files;
  *io_bytes_count += bytes;
  *(is_binary ? io_bytes_binary_count : io_bytes_ascii_count) += bytes;
}

// Count some communication (to be called from the communication routine)
void BeginTimingCommunication(cGH const *const cctkGH) {
  assert(timing_state == state_computing);
  timing_state = state_communicating;
  time_start = get_walltime();
}

void EndTimingCommunication(cGH const *const cctkGH, CCTK_REAL const messages,
                            CCTK_REAL const bytes) {
  DECLARE_CCTK_ARGUMENTS;

  assert(timing_state == state_communicating);
  timing_state = state_computing;
  CCTK_REAL const time_end = get_walltime();

  *time_communicating += time_end - time_start;

  *comm_count += messages;
  *comm_bytes_count += bytes;
}

static void UpdateTimes(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  assert(timing_state == state_computing);

  // Measure the elapsed time
  double const walltime = get_walltime();
  *time_total = walltime - startup_walltime;
  *time_evolution = in_evolution ? walltime - initial_walltime : 0.0;

  *time_computing = *time_total - (*time_communicating + *time_io);
}

static void UpdateUpdatesPerSecond(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  // Calculate updates per second
  *local_grid_points_per_second =
      *local_grid_point_updates_count / max(*time_computing, eps);
  *total_grid_points_per_second =
      *total_grid_point_updates_count / max(*time_computing, eps);
  *local_interior_points_per_second =
      *local_interior_point_updates_count / max(*time_computing, eps);
  *total_interior_points_per_second =
      *total_interior_point_updates_count / max(*time_computing, eps);

  *grid_points_per_second = *local_grid_points_per_second;
  *interior_points_per_second = *local_interior_points_per_second;
}

static void UpdateIOStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  *io_per_second = *io_count / max(*time_io, eps);
  *io_bytes_per_second = *io_bytes_count / max(*time_io, eps);
  *io_bytes_ascii_per_second = *io_bytes_ascii_count / max(*time_io, eps);
  *io_bytes_binary_per_second = *io_bytes_binary_count / max(*time_io, eps);
}

static void UpdateCommunicationStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  *comm_per_second = *comm_count / max(*time_communicating, eps);
  *comm_bytes_per_second = *comm_bytes_count / max(*time_communicating, eps);
}

static void UpdatePhysicalTimePerHour(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (not in_evolution) {
    *physical_time_per_hour = 0.0;
    *current_physical_time_per_hour = 0.0;
    return;
  }

  static int last_iteration = -1;
  static size_t num_samples = 0;
  // static CCTK_REAL last_physical_time;
  static CCTK_REAL last_time_evolution;
  assert(cctk_iteration > last_iteration); // expect progress

  // Calculate elapsed physical time
  CCTK_REAL const physical_time = cctkGH->cctk_time - initial_phystime;

  // Calculate average physical time per hour
  *physical_time_per_hour = 3600.0 * physical_time / max(*time_evolution, eps);

  // Calculate current physical time per hour as moving average
  if (last_iteration < 0) {
    // No past data are available
    *current_physical_time_per_hour = *physical_time_per_hour;
  } else if (num_samples < 3 or
             *time_evolution < 0.01 * timing_average_window_minutes * 60.0) {
    // Less than three previous samples are available, or less than
    // one percent of a window of past data are available
    *current_physical_time_per_hour = *physical_time_per_hour;
  } else {
    CCTK_REAL const window =
        min(*time_evolution, (CCTK_REAL)(timing_average_window_minutes * 60.0));
    CCTK_REAL const alpha =
        exp(-(*time_evolution - last_time_evolution) / window);
    *current_physical_time_per_hour = (1.0 - alpha) * *physical_time_per_hour +
                                      (alpha) * *current_physical_time_per_hour;
  }

  // Remember last iteration
  last_iteration = cctk_iteration;
  ++num_samples;
  // last_physical_time  = physical_time;
  last_time_evolution = *time_evolution;
}

static void UpdateMemoryStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  // Various metadata, excluding the small change
  *metadata = 0.0 +
              // memoryof(main_timer_tree) +
              // memoryof(mode_timer_tree) +
              memoryof(origin_space) + memoryof(delta_space) +
              // memoryof(domainspecs) +
              memoryof(vhh) + memoryof(vdd) + memoryof(level_regridding_epochs);
  // memoryof(groupdata) +
  // memoryof(arrdata);

  // TODO: Add TimerNode.root(?), TimerSet.timers

  // Grid structure
  *grid_structure = 0.0;
  // Storage for grid arrays
  *grid_arrays = 0.0;
  // Storage for grid functions
  *grid_functions = 0.0;

  for (size_t m = 0; m < vhh.size(); ++m) {
    *grid_structure += memoryof(*vhh.AT(m));
  }
  for (size_t m = 0; m < vdd.size(); ++m) {
    *grid_structure += memoryof(*vdd.AT(m));
  }
  *grid_structure += memoryof(*tt);

  const int myproc = CCTK_MyProc(0);
  for (size_t gi = 0; gi < arrdata.size(); ++gi) {
    if (CCTK_GroupTypeI(gi) == CCTK_GF) {
      for (size_t m = 0; m < arrdata.AT(gi).size(); ++m) {
        const gh &hh = *vhh.AT(m);
        for (size_t vi = 0; vi < arrdata.AT(gi).AT(m).data.size(); ++vi) {
          const ggf *f = arrdata.AT(gi).AT(m).data.AT(vi);
          if (f) {
            *grid_structure += memoryof(*f);
            for (int ml = 0; ml < hh.mglevels(); ++ml) {
              for (int rl = 0; rl < hh.reflevels(); ++rl) {
                const int tls = f->timelevels(ml, rl);
                for (int lc = 0; lc < hh.local_components(rl); ++lc) {
                  for (int tl = 0; tl < tls; ++tl) {
                    const gdata *d = f->data_pointer(tl, rl, lc, ml);
                    *grid_structure += memoryof(*d);
                    // TODO: add memoryof(d->_memory);
                    if (d->has_storage()) {
                      assert(d->proc() == myproc);
                      *grid_functions += d->size() * d->elementsize();
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else { // not CCTK_GF
      const int m = 0;
      *grid_structure += 0.0 + memoryof(*arrdata.AT(gi).AT(m).hh) +
                         memoryof(*arrdata.AT(gi).AT(m).dd) +
                         memoryof(*arrdata.AT(gi).AT(m).tt);
      for (size_t vi = 0; vi < arrdata.AT(gi).AT(m).data.size(); ++vi) {
        const ggf *f = arrdata.AT(gi).AT(m).data.AT(vi);
        if (f) {
          *grid_structure += memoryof(*f);
          const int ml = 0;
          const int rl = 0;
          const int tls = f->timelevels(ml, rl);
          const int lc = 0;
          for (int tl = 0; tl < tls; ++tl) {
            const gdata *d = f->data_pointer(tl, rl, lc, ml);
            *grid_arrays += memoryof(*d);
            // TODO: add memoryof(d->_memory);
            if (d->has_storage()) {
              assert(d->proc() == myproc);
              *grid_functions += d->size() * d->elementsize();
            }
          }
        }
      }
    }
  }
}

// Calculate timing statistics (to be called before output)
void UpdateTimingStats(cGH const *const cctkGH) {
  UpdateTimes(cctkGH);
  UpdateUpdatesPerSecond(cctkGH);
  UpdateIOStats(cctkGH);
  UpdateCommunicationStats(cctkGH);
  UpdatePhysicalTimePerHour(cctkGH);
  UpdateMemoryStats(cctkGH);
}

static void PrintTimes(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  CCTK_VInfo(CCTK_THORNSTRING, "Total run       time: %g h   (iteration %d)",
             double(*time_total / 3600.0), cctk_iteration);
  CCTK_VInfo(CCTK_THORNSTRING, "Total evolution time: %g h   (%d steps)",
             double(*time_evolution / 3600.0), int(*evolution_steps_count));
  CCTK_VInfo(CCTK_THORNSTRING,
             "(Comp, Comm, I/O) fractions = (%3.1f%%, %3.1f%%, %3.1f%%)",
             double(100.0 * *time_computing / max(*time_total, eps)),
             double(100.0 * *time_communicating / max(*time_total, eps)),
             double(100.0 * *time_io / max(*time_total, eps)));
}

static void PrintUpdatesPerSecond(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING,
             "This processor's grid point updates per second (local)    : %g",
             double(*local_grid_points_per_second));
  CCTK_VInfo(CCTK_THORNSTRING,
             "Overall grid point updates per second (total)             : %g",
             double(*total_grid_points_per_second));
  CCTK_VInfo(CCTK_THORNSTRING,
             "This processor's interior point updates per second (local): %g",
             double(*local_interior_points_per_second));
  CCTK_VInfo(CCTK_THORNSTRING,
             "Overall interior point updates per second (total)         : %g",
             double(*total_interior_points_per_second));

#if 0
    CCTK_REAL const updates_per_second_2 = ipow (updates_per_second, 2);
    
    struct {
      CCTK_REAL ups, ups2;
    } local, global;
    local.ups  = updates_per_second;
    local.ups2 = updates_per_second_2;
    MPI_Allreduce (& local, & global, 2,
                   dist::datatype (global.ups), MPI_SUM, dist::comm());
    
    int const count = dist::size();
    CCTK_REAL const avg = global.ups / count;
    CCTK_REAL const stddev = sqrt (fabs (global.ups2 - ipow (avg,2)) / count);
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "Local updates per second:   %g",
                double (updates_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "Global updates per second:  %g", double (global.ups));
    
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Average updates per second: %g", double (avg));
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Standard deviation:         %g", double (stddev));
    }
#endif
}

static void PrintIOStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  CCTK_VInfo(CCTK_THORNSTRING, "I/O operations per second:     %g",
             double(*io_per_second));
  CCTK_VInfo(CCTK_THORNSTRING, "I/O bytes per second:          %g",
             double(*io_bytes_per_second));
  CCTK_VInfo(CCTK_THORNSTRING, "I/O bytes per second (ASCII):  %g",
             double(*io_bytes_ascii_per_second));
  CCTK_VInfo(CCTK_THORNSTRING, "I/O bytes per second (binary): %g",
             double(*io_bytes_binary_per_second));
}

static void PrintCommunicationStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  CCTK_VInfo(CCTK_THORNSTRING, "Communication operations per second: %g",
             double(*comm_per_second));
  CCTK_VInfo(CCTK_THORNSTRING, "Communicated bytes per second:       %g",
             double(*comm_bytes_per_second));
}

static void PrintPhysicalTimePerHour(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  CCTK_VInfo(CCTK_THORNSTRING, "Current physical time per hour: %g",
             double(*current_physical_time_per_hour));
  CCTK_VInfo(CCTK_THORNSTRING, "Average physical time per hour: %g",
             double(*physical_time_per_hour));
}

#define PRINTMEM(x) (memoryof(x) / 1.0e+6) << " MB"

static void PrintMemoryStats(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;

  cout << eol << "Memory statistics:" << eol << "   Grid hierarchy:" << eol;
  for (int m = 0; m < Carpet::maps; ++m) {
    cout << "   gh[" << m << "]: " << PRINTMEM(*vhh.AT(m)) << eol << "   dh["
         << m << "]: " << PRINTMEM(*vdd.AT(m)) << eol;
  }
  cout << "   th: " << PRINTMEM(*tt) << eol;
#if 0
    for (int g = 0; g < (int)arrdata.size(); ++ g) {
      if (CCTK_GroupTypeI(g) != CCTK_GF) {
        char * const groupname = CCTK_GroupName(g);
        for (int m = 0; m < (int)arrdata.AT(g).size(); ++ m) {
          cout << "   Group " << groupname << ":" << eol
               << "   gh[" << m << "]: " << PRINTMEM(*arrdata.AT(g).AT(m).hh) << eol
               << "   dh[" << m << "]: " << PRINTMEM(*arrdata.AT(g).AT(m).dd) << eol
               << "   th[" << m << "]: " << PRINTMEM(*arrdata.AT(g).AT(m).tt) << eol;
          for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++ v) {
            char * const fullname = CCTK_FullName(CCTK_FirstVarIndexI(g)+v);
            cout << "   Variable " << fullname << ":" << eol
                 << "   ggf[" << m << "]: ";
            if (arrdata.AT(g).AT(m).data.AT(v)) {
              cout << PRINTMEM(*arrdata.AT(g).AT(m).data.AT(v));
            } else {
              cout << "<null>";
            }
            cout << eol;
            free (fullname);
          }
        }
        free (groupname);
      }
    }
#endif
  cout << endl;
}

// Calculate and print some timing statistics (to be called at any
// time)
void PrintTimingStats(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  if (print_timestats_every > 0 and
      cctkGH->cctk_iteration % print_timestats_every == 0) {
    PrintTimes(cctkGH);
    PrintUpdatesPerSecond(cctkGH);
    PrintIOStats(cctkGH);
    PrintCommunicationStats(cctkGH);
    PrintPhysicalTimePerHour(cctkGH);
    PrintMemoryStats(cctkGH);
  }
}

} // namespace Carpet
