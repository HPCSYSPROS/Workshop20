#include "cycleclock.h"

#ifdef HAVE_TICK_COUNTER

#include <cctk.h>

// Find a good wall clock timer
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

static double cycleclock_tick = -1.0; // uninitialised

#ifdef HAVE_SYS_TIME_H
static double get_sys_time() {
  timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec + 1.0e-6 * tp.tv_usec;
}
#endif

void measure_tick() {
  // Make a few warm-up measurements
  getticks();
  getticks();
  getticks();

#ifdef _OPENMP
  if (cycleclock_tick < 0.0) {
    // Use omp_get_sys_time to calibrate the timer
    CCTK_INFO("Measuring CycleClock tick via OpenMP...");
    ticks const rstart = getticks();
    double const wstart = omp_get_wtime();
    while (omp_get_wtime() < wstart + 0.1) {
      // do nothing, just wait
    }
    ticks const rend = getticks();
    double const wend = omp_get_wtime();
    cycleclock_tick = (wend - wstart) / elapsed(rend, rstart);
  }
#endif

#ifdef HAVE_CAPABILITY_MPI
  if (cycleclock_tick < 0.0) {
    // Use MPI_Wtime to calibrate the timer
    CCTK_INFO("Measuring CycleClock tick via MPI...");
    ticks const rstart = getticks();
    double const wstart = MPI_Wtime();
    while (MPI_Wtime() < wstart + 0.1) {
      // do nothing, just wait
    }
    ticks const rend = getticks();
    double const wend = MPI_Wtime();
    cycleclock_tick = (wend - wstart) / elapsed(rend, rstart);
  }
#endif

#ifdef HAVE_SYS_TIME_H
  if (cycleclock_tick < 0.0) {
    // Use gettimeofday to calibrate the timer
    CCTK_INFO("Measuring CycleClock tick via gettimeofday...");
    ticks const rstart = getticks();
    double const wstart = get_sys_time();
    while (get_sys_time() < wstart + 0.1) {
      // do nothing, just wait
    }
    ticks const rend = getticks();
    double const wend = get_sys_time();
    cycleclock_tick = (wend - wstart) / elapsed(rend, rstart);
  }
#endif

  if (cycleclock_tick < 0.0) {
    // Give up; just set the time scale to 1
    CCTK_INFO("Could not measure CycleClock tick");
    cycleclock_tick = 1.0;
  }

  CCTK_VInfo(CCTK_THORNSTRING,
             "Calibrated CycleClock: %g ns per clock tick (%g GHz)",
             1.0e9 * cycleclock_tick, 1.0e-9 / cycleclock_tick);
}

double seconds_per_tick() {
  // if (CCTK_BUILTIN_EXPECT(cycleclock_tick < 0.0, false)) {
  //   CCTK_WARN(CCTK_WARN_ALERT,
  //             "Thorn CycleClock has not been activated; the first
  //             measurements may be wrong");
  //   measure_tick();
  // }
  return cycleclock_tick;
}

#endif
