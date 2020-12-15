#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Termination.h>
#include <util_String.h>

#include <Requirements.hh>

#include <CactusTimerSet.hh>
#include <Timer.hh>

#include <dist.hh>
#include <th.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

static bool do_terminate(const cGH *cctkGH);

static void AdvanceTime(cGH *cctkGH);
static void CallRegrid(cGH *cctkGH);
static void CallEvol(cGH *cctkGH);
static void CallRestrict(cGH *cctkGH);
static void CallAnalysis(cGH *cctkGH);

static void print_internal_data();

static void ScheduleTraverse(char const *where, char const *name, cGH *cctkGH);
static void OutputGH(char const *where, cGH *cctkGH);

int Evolve(tFleshConfig *const fc) {
  DECLARE_CCTK_PARAMETERS;

  Waypoint("Starting evolution loop");

  int const convlev = 0;
  cGH *cctkGH = fc->GH[convlev];

  // Main loop
  BeginTimingEvolution(cctkGH);
  static Timers::Timer timer("Evolve");
  timer.start();
  while (not do_terminate(cctkGH)) {

    AdvanceTime(cctkGH);
    {
      int const do_every = maxtimereflevelfact / timereffacts.AT(reflevels - 1);
      if ((cctkGH->cctk_iteration - 1) % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, 0) {
          BEGIN_REFLEVEL_LOOP(cctkGH) { CallRegrid(cctkGH); }
          END_REFLEVEL_LOOP;
        }
        LEAVE_GLOBAL_MODE;
      }
    }
    CallEvol(cctkGH);
    CallRestrict(cctkGH);
    CallAnalysis(cctkGH);
    print_internal_data();

    // Print timer values
    {
      Timers::Timer timer("PrintTimers");
      timer.start();
      int const do_every = maxtimereflevelfact / timereffacts.AT(reflevels - 1);
      if (output_timers_every > 0 and
          cctkGH->cctk_iteration % output_timers_every == 0 and
          cctkGH->cctk_iteration % do_every == 0) {
        Timers::CactusTimerSet::writeData(cctkGH, timer_file);
      }

      if (output_timer_tree_every > 0 and
          cctkGH->cctk_iteration % output_timer_tree_every == 0 and
          cctkGH->cctk_iteration % do_every == 0) {
        Timers::Timer::outputTree("Evolve");
      }
      timer.stop();
    }

    // Ensure that all levels have consistent times
    {
      Timers::Timer timer("CheckLevelTimes");
      timer.start();
      CCTK_REAL const eps =
          pow(numeric_limits<CCTK_REAL>::epsilon(), CCTK_REAL(0.75));
      assert(fabs(cctkGH->cctk_time - global_time) <= eps * global_time);
      for (int ml = 0; ml < mglevels; ++ml) {
        for (int rl = 0; rl < reflevels; ++rl) {
          int const do_every =
              ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
          if (cctkGH->cctk_iteration % do_every == 0) {
            // assert (fabs (leveltimes.AT(ml).AT(rl) - global_time) <=
            //         eps * global_time);
            assert(fabs(tt->get_time(ml, rl, 0) - global_time) <=
                   eps * global_time);
          }
        }
      }
      timer.stop();
    }

  } // end main loop
  timer.stop();

  Waypoint("Done with evolution loop");

  return 0;
}

bool do_terminate(const cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("DoTerminate");
  timer.start();

  bool term;

  // Do not test on non-active reflevels to save the call to
  // MPI_Allreduce below
  int const do_every = maxtimereflevelfact / timereffacts.AT(reflevels - 1);
  if (cctkGH->cctk_iteration % do_every != 0) {

    term = false;

  } else {

    if (terminate_next or CCTK_TerminationReached(cctkGH)) {

      // Terminate if someone or something said so
      term = true;

    } else {

      // Test the various conditions
      bool const term_iter =
          cctk_itlast >= 0 and cctkGH->cctk_iteration >= cctk_itlast;
      bool const term_time =
          (delta_time > 0.0
               ? (cctkGH->cctk_time >=
                  cctk_final_time - 1.0e-8 * cctkGH->cctk_delta_time)
               : (cctkGH->cctk_time <=
                  cctk_final_time + 1.0e-8 * cctkGH->cctk_delta_time));
      bool const term_runtime =
          max_runtime > 0.0 and CCTK_RunTime() >= 60.0 * max_runtime;

      if (CCTK_Equals(terminate, "never")) {
        term = false;
      } else if (CCTK_Equals(terminate, "iteration")) {
        term = term_iter;
      } else if (CCTK_Equals(terminate, "time")) {
        term = term_time;
        if (term) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Terminating due to cctk_final_time at t = %f",
                     cctkGH->cctk_time);
        }

      } else if (CCTK_Equals(terminate, "runtime")) {
        term = term_runtime;
      } else if (CCTK_Equals(terminate, "any")) {
        term = term_iter or term_time or term_runtime;
      } else if (CCTK_Equals(terminate, "all")) {
        term = term_iter and term_time and term_runtime;
      } else if (CCTK_Equals(terminate, "either")) {
        term = term_iter or term_time;
      } else if (CCTK_Equals(terminate, "both")) {
        term = term_iter and term_time;
      } else if (CCTK_Equals(terminate, "immediately")) {
        term = true;
      } else {
        CCTK_WARN(0, "Unsupported termination condition");
        abort(); // keep the compiler happy
      }
    }

    // Reduce termination condition
    int local, global;
    local = term;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LOR, dist::comm());
    term = global;
  }

  timer.stop();
  return term;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AdvanceTime(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("AdvanceTime");
  timer.start();

  Checkpoint("AdvanceTime");

  ++cctkGH->cctk_iteration;

  if (not adaptive_stepsize) {
    // Avoid accumulation of errors
    global_time = cctk_initial_time +
                  cctkGH->cctk_iteration * delta_time / maxtimereflevelfact;
    cctkGH->cctk_time = global_time;
  } else {
    // Take varying step sizes into account
    cctkGH->cctk_time += cctkGH->cctk_delta_time;
    delta_time = cctkGH->cctk_delta_time;
    global_time = cctkGH->cctk_time;
  }

  if ((cctkGH->cctk_iteration - 1) %
          (maxtimereflevelfact / timereffacts.AT(reflevels - 1)) ==
      0) {
    Waypoint("Evolving iteration %d at t=%g", cctkGH->cctk_iteration,
             (double)cctkGH->cctk_time);
  }

  timer.stop();
}

void CallRegrid(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CallRegrid";
  static Timers::Timer timer(where);
  timer.start();

  assert(is_level_mode());

  bool const old_do_global_mode = do_global_mode;
  bool const old_do_early_global_mode = do_early_global_mode;
  bool const old_do_late_global_mode = do_late_global_mode;
  bool const old_do_meta_mode = do_meta_mode;
  bool const old_do_early_meta_mode = do_early_meta_mode;
  bool const old_do_late_meta_mode = do_late_meta_mode;
  do_global_mode = true;
  do_early_global_mode = true;
  do_late_global_mode = true;
  do_meta_mode = true;
  do_early_meta_mode = true;
  do_late_meta_mode = true;

  // Preregrid
  Waypoint("Preregrid at iteration %d time %g%s%s", cctkGH->cctk_iteration,
           (double)cctkGH->cctk_time, (do_global_mode ? " (global)" : ""),
           (do_meta_mode ? " (meta)" : ""));
  ScheduleTraverse(where, "CCTK_PREREGRID", cctkGH);

  // Regrid
  Checkpoint("Regrid");
  int const oldreflevels = reflevels;
  bool const did_regrid = Regrid(cctkGH, false, true);
  bool const did_remove_level = reflevels < oldreflevels;
  assert(not did_remove_level or did_regrid);

  if (did_regrid) {
#ifdef REQUIREMENTS_HH
    Requirements::Regrid(reflevels);
#endif
    bool did_any_recompose = false;
    BEGIN_META_MODE(cctkGH) {

      bool have_done_global_mode = false;
      bool have_done_early_global_mode = false;
      bool have_done_late_global_mode = false;
      bool have_done_anything = false;

      for (int rl = 0; rl < reflevels; ++rl) {

        bool const did_recompose = Recompose(cctkGH, rl, true);
        did_any_recompose = did_any_recompose or did_recompose;
#ifdef REQUIREMENTS_HH
        Requirements::Recompose(cctkGH->cctk_iteration, rl,
                                not did_recompose
                                    ? Requirements::valid::everywhere
                                    : Requirements::valid::interior);
#endif

        // Carpet assumes that a regridding operation always changes
        // "level N and all finer levels" so we should call POSTREGRID
        // on all finer levels
        if (did_any_recompose or (did_remove_level and rl == reflevels - 1)) {
          BEGIN_MGLEVEL_LOOP(cctkGH) {
            ENTER_LEVEL_MODE(cctkGH, rl) {
              BeginTimingLevel(cctkGH);

              do_early_global_mode = not have_done_early_global_mode;
              do_late_global_mode = reflevel == reflevels - 1;
              do_early_meta_mode =
                  do_early_global_mode and mglevel == mglevels - 1;
              do_late_meta_mode = do_late_global_mode and mglevel == 0;
              do_global_mode = do_late_global_mode;
              do_meta_mode = do_late_meta_mode;
              assert(not(have_done_global_mode and do_global_mode));
              assert(not(have_done_early_global_mode and do_early_global_mode));
              assert(not(have_done_late_global_mode and do_late_global_mode));
              have_done_global_mode |= do_global_mode;
              have_done_early_global_mode |= do_early_global_mode;
              have_done_late_global_mode |= do_late_global_mode;
              have_done_anything = true;

              BEGIN_TIMELEVEL_LOOP(cctkGH) {

                Waypoint("Postregrid at iteration %d time %g timelevel %d%s%s",
                         cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                         timelevel, (do_global_mode ? " (global)" : ""),
                         (do_meta_mode ? " (meta)" : ""));

                // Postregrid
                ScheduleTraverse(where, "CCTK_POSTREGRID", cctkGH);
              }
              END_TIMELEVEL_LOOP;

              if (output_after_regridding) {
                // Output
                OutputGH(where, cctkGH);
              }

              EndTimingLevel(cctkGH);
            }
            LEAVE_LEVEL_MODE;
          }
          END_MGLEVEL_LOOP;

        } // if did_recompose

      } // for rl

      if (have_done_anything)
        assert(have_done_global_mode);
      if (have_done_anything)
        assert(have_done_early_global_mode);
      if (have_done_anything)
        assert(have_done_late_global_mode);
    }
    END_META_MODE;
#ifdef REQUIREMENTS_HH
    Requirements::RegridFree();
#endif
  } // if did_regrid

  RegridFree(cctkGH, true);

  do_global_mode = old_do_global_mode;
  do_early_global_mode = old_do_early_global_mode;
  do_late_global_mode = old_do_late_global_mode;
  do_meta_mode = old_do_meta_mode;
  do_early_meta_mode = old_do_early_meta_mode;
  do_late_meta_mode = old_do_late_meta_mode;

  timer.stop();
}

void CallEvol(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CallEvol";
  static Timers::Timer timer(where);
  timer.start();

  for (int ml = mglevels - 1; ml >= 0; --ml) {

    bool have_done_global_mode = false;
    bool have_done_anything = false;

    for (int rl = 0; rl < reflevels; ++rl) {
      int const do_every =
          ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
      if ((cctkGH->cctk_iteration - 1) % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, ml) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            BeginTimingLevel(cctkGH);

            do_early_global_mode = not have_done_global_mode;
            do_late_global_mode = reflevel == reflevels - 1;
            do_early_meta_mode =
                do_early_global_mode and mglevel == mglevels - 1;
            do_late_meta_mode = do_late_global_mode and mglevel == 0;
            do_global_mode = do_early_global_mode;
            do_meta_mode = do_early_meta_mode;
            assert(not(have_done_global_mode and do_global_mode));
            have_done_global_mode |= do_global_mode;
            have_done_anything = true;

            if (use_tapered_grids and reflevel > 0) {
              int const parent_do_every =
                  ipow(mgfact, mglevel) *
                  (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
              bool const parent_is_active =
                  (cctkGH->cctk_iteration - 1) % parent_do_every == 0;
              do_taper = not parent_is_active;
            }

            // Advance times
            CycleTimeLevels(cctkGH);
            if (not adaptive_stepsize) {
              cctkGH->cctk_time =
                  (global_time - delta_time / maxtimereflevelfact +
                   delta_time * mglevelfact / timereflevelfact);
            }
            tt->set_time(mglevel, reflevel, timelevel, cctkGH->cctk_time);

            Waypoint("Evolution I at iteration %d time %g%s%s%s",
                     cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                     (do_global_mode ? " (global)" : ""),
                     (do_meta_mode ? " (meta)" : ""),
                     (do_taper ? " (tapering)" : ""));

            // Checking
            CalculateChecksums(cctkGH, allbutcurrenttime);
            Poison(cctkGH, currenttimebutnotifonly);

            // Evolve
            ScheduleTraverse(where, "CCTK_PRESTEP", cctkGH);
            ScheduleTraverse(where, "CCTK_EVOL", cctkGH);

            // Checking
            PoisonCheck(cctkGH, currenttime);

            // Timing statistics
            StepTimingEvolution(cctkGH);

            do_taper = false;

            EndTimingLevel(cctkGH);
          }
          LEAVE_LEVEL_MODE;
        }
        LEAVE_GLOBAL_MODE;
      } // if do_every
    }   // for rl

    if (have_done_anything)
      assert(have_done_global_mode);

  } // for ml

  timer.stop();
}

void CallRestrict(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "Evolve::CallRestrict";
  static Timers::Timer timer("CallRestrict");
  timer.start();

  for (int ml = mglevels - 1; ml >= 0; --ml) {

    bool did_restrict = false;

    for (int rl = reflevels - 2; rl >= 0; --rl) {
      int const do_every =
          ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
      if (cctkGH->cctk_iteration % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, ml) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            BeginTimingLevel(cctkGH);

            Waypoint("Evolution/Restrict at iteration %d time %g",
                     cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

            Restrict(cctkGH);
            did_restrict = true;

            if (use_higher_order_restriction) {
              do_early_global_mode = false;
              do_late_global_mode = false;
              do_early_meta_mode = false;
              do_late_meta_mode = false;
              do_global_mode = false;
              do_meta_mode = false;

              if (use_tapered_grids and reflevel > 0) {
                int const parent_do_every =
                    ipow(mgfact, mglevel) *
                    (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
                bool const parent_is_active =
                    (cctkGH->cctk_iteration - 1) % parent_do_every == 0;
                do_taper = not parent_is_active;
              }

              // TODO: disable prolongation (but not
              // synchronization)
              Waypoint("Evolution/PostRestrict (intermediate) at iteration %d "
                       "time %g",
                       cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

              ScheduleTraverse(where, "CCTK_POSTRESTRICT", cctkGH);

              do_taper = false;
            }

            EndTimingLevel(cctkGH);
          }
          LEAVE_LEVEL_MODE;
        }
        LEAVE_GLOBAL_MODE;
      } // if do_every
    }   // for rl

    if (did_restrict) {

      bool have_done_global_mode = false;
      bool have_done_early_global_mode = false;
      bool have_done_late_global_mode = false;
      bool have_done_anything = false;

      for (int rl = 0; rl < reflevels; ++rl) {
        int const do_every =
            ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
        if (cctkGH->cctk_iteration % do_every == 0) {
          ENTER_GLOBAL_MODE(cctkGH, ml) {
            ENTER_LEVEL_MODE(cctkGH, rl) {
              BeginTimingLevel(cctkGH);

              do_early_global_mode = not have_done_early_global_mode;
              do_late_global_mode = reflevel == reflevels - 1;
              do_early_meta_mode =
                  do_early_global_mode and mglevel == mglevels - 1;
              do_late_meta_mode = do_late_global_mode and mglevel == 0;
              do_global_mode = do_late_global_mode;
              do_meta_mode = do_global_mode and do_late_meta_mode;
              assert(not(have_done_global_mode and do_global_mode));
              assert(not(have_done_early_global_mode and do_early_global_mode));
              assert(not(have_done_late_global_mode and do_late_global_mode));
              have_done_global_mode |= do_global_mode;
              have_done_early_global_mode |= do_early_global_mode;
              have_done_late_global_mode |= do_late_global_mode;
              have_done_anything = true;

              if (use_tapered_grids and reflevel > 0) {
                int const parent_do_every =
                    ipow(mgfact, mglevel) *
                    (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
                bool const parent_is_active =
                    (cctkGH->cctk_iteration - 1) % parent_do_every == 0;
                do_taper = not parent_is_active;
              }

              Waypoint("Evolution/PostRestrict at iteration %d time %g",
                       cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

              ScheduleTraverse(where, "CCTK_POSTRESTRICT", cctkGH);

              do_taper = false;

              EndTimingLevel(cctkGH);
            }
            LEAVE_LEVEL_MODE;
          }
          LEAVE_GLOBAL_MODE;
        } // if do_every
      }   // for rl

      if (have_done_anything)
        assert(have_done_global_mode);
      if (have_done_anything)
        assert(have_done_early_global_mode);
      if (have_done_anything)
        assert(have_done_late_global_mode);

    } // if did_restrict

  } // for ml

  timer.stop();
}

void CallAnalysis(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CallAnalysis";
  static Timers::Timer timer(where);
  timer.start();

  for (int ml = mglevels - 1; ml >= 0; --ml) {

    bool have_done_global_mode = false;
    bool have_done_early_global_mode = false;
    bool have_done_late_global_mode = false;
    bool have_done_anything = false;

    for (int rl = 0; rl < reflevels; ++rl) {
      int const do_every =
          ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
      if (cctkGH->cctk_iteration % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, ml) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            BeginTimingLevel(cctkGH);

            do_early_global_mode = not have_done_early_global_mode;
            do_late_global_mode = reflevel == reflevels - 1;
            do_early_meta_mode =
                do_early_global_mode and mglevel == mglevels - 1;
            do_late_meta_mode = do_late_global_mode and mglevel == 0;
            do_global_mode = do_late_global_mode;
            do_meta_mode = do_global_mode and do_late_meta_mode;
            assert(not(have_done_global_mode and do_global_mode));
            assert(not(have_done_early_global_mode and do_early_global_mode));
            assert(not(have_done_late_global_mode and do_late_global_mode));
            have_done_global_mode |= do_global_mode;
            have_done_early_global_mode |= do_early_global_mode;
            have_done_late_global_mode |= do_late_global_mode;
            have_done_anything = true;

            if (use_tapered_grids and reflevel > 0) {
              int const parent_do_every =
                  ipow(mgfact, mglevel) *
                  (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
              bool const parent_is_active =
                  (cctkGH->cctk_iteration - 1) % parent_do_every == 0;
              do_taper = not parent_is_active;
            }

            Waypoint("Evolution II at iteration %d time %g%s%s%s",
                     cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                     (do_global_mode ? " (global)" : ""),
                     (do_meta_mode ? " (meta)" : ""),
                     (do_taper ? " (tapering)" : ""));

#if 0
              if (reflevel < reflevels-1) {
                ScheduleTraverse (where, "CCTK_POSTRESTRICT", cctkGH);
              }
#endif

            // Poststep
            ScheduleTraverse(where, "CCTK_POSTSTEP", cctkGH);

            // Checking
            PoisonCheck(cctkGH, currenttime);
            CalculateChecksums(cctkGH, currenttime);

            // Checkpoint
            ScheduleTraverse(where, "CCTK_CHECKPOINT", cctkGH);

            // Analysis
            in_analysis_bin = true;
            ScheduleTraverse(where, "CCTK_ANALYSIS", cctkGH);
            in_analysis_bin = false;

            if (do_late_global_mode) {
              // Timing statistics
              UpdateTimingStats(cctkGH);
            }

            // Output
            OutputGH(where, cctkGH);

            // Checking
            CheckChecksums(cctkGH, alltimes);

            if (do_late_global_mode) {
              // Timing statistics
              PrintTimingStats(cctkGH);
            }

            do_taper = false;

            EndTimingLevel(cctkGH);
          }
          LEAVE_LEVEL_MODE;
        }
        LEAVE_GLOBAL_MODE;
      } // if do_every
    }   // for rl

    if (have_done_anything)
      assert(have_done_global_mode);
    if (have_done_anything)
      assert(have_done_early_global_mode);
    if (have_done_anything)
      assert(have_done_late_global_mode);

  } // for ml

  timer.stop();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void print_internal_data() {
  DECLARE_CCTK_PARAMETERS;

  if (output_internal_data) {
    CCTK_INFO("Internal data dump:");
    streamsize const oldprecision = cout.precision();
    cout.precision(17);
    cout << "   global_time: " << global_time << endl
         // << "   leveltimes: " << leveltimes << endl
         << "   delta_time: " << delta_time << endl;
    cout.precision(oldprecision);
  }
}

void ScheduleTraverse(char const *const where, char const *const name,
                      cGH *const cctkGH) {
  Timers::Timer timer(name);
  timer.start();
  ostringstream infobuf;
  infobuf << "Scheduling " << name;
  string const info = infobuf.str();
  Checkpoint(info.c_str());
  CCTK_ScheduleTraverse(name, cctkGH, CallFunction);
  timer.stop();
}

void OutputGH(char const *const where, cGH *const cctkGH) {
  static Timers::Timer timer("OutputGH");
  timer.start();
  CCTK_OutputGH(cctkGH);
  timer.stop();
}

} // namespace Carpet
