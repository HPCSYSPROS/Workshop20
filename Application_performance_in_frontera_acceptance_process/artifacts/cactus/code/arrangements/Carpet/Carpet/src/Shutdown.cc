#include <cstdio>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <CactusTimerSet.hh>
#include <Timer.hh>

#include <dist.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

int Shutdown(tFleshConfig *fc) {
  DECLARE_CCTK_PARAMETERS;

  Waypoint("Starting shutdown");

  const int convlev = 0;
  cGH *cctkGH = fc->GH[convlev];

  static Timers::Timer timer("Shutdown");
  timer.start();
  for (int rl = reflevels - 1; rl >= 0; --rl) {
    BEGIN_REVERSE_MGLEVEL_LOOP(cctkGH) {
      ENTER_LEVEL_MODE(cctkGH, rl) {
        BeginTimingLevel(cctkGH);

        do_early_global_mode = reflevel == reflevels - 1;
        do_late_global_mode = reflevel == 0;
        do_early_meta_mode = do_early_global_mode and mglevel == 0;
        do_late_meta_mode = do_late_global_mode and mglevel == mglevels - 1;
        do_global_mode = do_late_global_mode;
        do_meta_mode = do_late_meta_mode;

        Checkpoint("Shutdown at iteration %d time %g%s%s",
                   cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                   (do_global_mode ? " (global)" : ""),
                   (do_meta_mode ? " (meta)" : ""));

        // Terminate
        Checkpoint("Scheduling TERMINATE");
        CCTK_ScheduleTraverse("CCTK_TERMINATE", cctkGH, CallFunction);

        EndTimingLevel(cctkGH);
      }
      LEAVE_LEVEL_MODE;
    }
    END_REVERSE_MGLEVEL_LOOP;
  } // for rl

  // Stop all timers before shutdown, since timers may rely on data
  // structures which are destroyed during shutdown
  int const ierr = CCTK_TimerStop("CCTK total time");
  assert(not ierr);
  timer.stop();
  if (output_timers_every > 0) {
    Timers::CactusTimerSet::writeData(cctkGH, timer_file);
  }

  if (output_timer_tree_every > 0) {
    Timers::Timer::outputTree("Evolve");
  }

  if (output_xml_timer_tree) {
    Timers::Timer::outputTreeXML();
  }

  BEGIN_REVERSE_MGLEVEL_LOOP(cctkGH) {
    do_early_global_mode = true;
    do_late_global_mode = true;
    do_early_meta_mode = do_early_global_mode and mglevel == 0;
    do_late_meta_mode = do_late_global_mode and mglevel == mglevels - 1;
    do_global_mode = do_late_global_mode;
    do_meta_mode = do_late_meta_mode;

    // Shutdown
    Checkpoint("Scheduling SHUTDOWN");
    CCTK_ScheduleTraverse("CCTK_SHUTDOWN", cctkGH, CallFunction);
  }
  END_REVERSE_MGLEVEL_LOOP;

  // Free all memory, call all destructors
  for (size_t gi = 0; gi < arrdata.size(); ++gi) {
    for (size_t m = 0; m < arrdata.AT(gi).size(); ++m) {
      for (size_t vi = 0; vi < arrdata.AT(gi).AT(m).data.size(); ++vi) {
        ggf *&f = arrdata.AT(gi).AT(m).data.AT(vi);
        if (f) {
          delete f;
          f = 0;
        }
      }
    }
  }
  // for (int gi=0; gi<CCTK_NumGroups(); ++gi) {
  //   const int tls = 0;
  //   CCTK_GroupStorageDecrease(cctkGH, 1, &gi, &tls, 0);
  // }

  for (size_t gi = 0; gi < arrdata.size(); ++gi) {
    if (CCTK_GroupTypeI(gi) == CCTK_GF) {
      for (size_t m = 0; m < arrdata.AT(gi).size(); ++m) {
        arrdata.AT(gi).AT(m).tt = 0;
        arrdata.AT(gi).AT(m).dd = 0;
        arrdata.AT(gi).AT(m).hh = 0;
      }
    } else {
      for (size_t m = 0; m < arrdata.AT(gi).size(); ++m) {
        delete arrdata.AT(gi).AT(m).tt;
        arrdata.AT(gi).AT(m).tt = 0;
        delete arrdata.AT(gi).AT(m).dd;
        arrdata.AT(gi).AT(m).dd = 0;
        delete arrdata.AT(gi).AT(m).hh;
        arrdata.AT(gi).AT(m).hh = 0;
      }
    }
  }

  delete tt;
  tt = 0;
  for (size_t m = 0; m < vdd.size(); ++m) {
    delete vdd.AT(m);
    vdd.AT(m) = 0;
  }
  for (size_t m = 0; m < vhh.size(); ++m) {
    delete vhh.AT(m);
    vhh.AT(m) = 0;
  }

  // earlier checkpoint before finalising MPI
  Waypoint("Done with shutdown");

  return 0;
}

} // namespace Carpet
