#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include <sys/time.h>
#include <unistd.h>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"

namespace CarpetLib {

using namespace std;

// Call a timer
static ticks call_timer() { return getticks(); }

// A global timer set
static TimerSet *timerSet = NULL;

// Add a timer
void TimerSet::add(Timer *const timer) { timers.push_back(timer); }

// Remove a timer
void TimerSet::remove(Timer *const timer) { timers.remove(timer); }

// Output all timer names
void TimerSet::outputNames(ostream &os) const {
  os << "Timer names:" << eol;
  int n = 0;
  for (list<Timer *>::const_iterator itimer = timers.begin();
       itimer != timers.end(); ++itimer) {
    os << "   [" << setw(4) << setfill('0') << n << "] " << (*itimer)->name()
       << eol;
    ++n;
  }
}

// Output all timer data
void TimerSet::outputData(ostream &os) const {
  for (list<Timer *>::const_iterator itimer = timers.begin();
       itimer != timers.end(); ++itimer) {
    os << *(*itimer);
  }
}

// Create a new timer with the given name
Timer::Timer(char const *const timername_) : timername(timername_) {
  assert(timername_);
  resetstats();
  if (not timerSet)
    timerSet = new TimerSet;
  timerSet->add(this);
}

// Destroy a timer
Timer::~Timer() {
  assert(timerSet);
  timerSet->remove(this);
}

// Reset the statistics
void Timer::resetstats() {
  wtime = 0.0;
  wtime2 = 0.0;
  wmin = 0.0;
  wmax = 0.0;

  bytes = 0.0;
  bytes2 = 0.0;
  bmin = 0.0;
  bmax = 0.0;

  count = 0.0;

  running = false;
}

// Add statistics of a timing operation
void Timer::addstat(double const t, double const b) {
  wtime += t;
  wtime2 += pow(t, 2);
  wmin = min(wmin, t);
  wmax = max(wmax, t);

  bytes += b;
  bytes2 += pow(b, 2);
  bmin = min(bmin, b);
  bmax = max(bmax, b);

  ++count;
}

// Start the timer
void Timer::start() {
  DECLARE_CCTK_PARAMETERS;
  assert(not running);
  running = true;
  if (use_ipm_timing_regions) {
    MPI_Pcontrol(+1, timername.c_str());
  }
  starttime = call_timer();
}

// Stop the timer
void Timer::stop(double const b) {
  DECLARE_CCTK_PARAMETERS;
  assert(running);
  running = false;
  ticks const endtime = call_timer();
  if (use_ipm_timing_regions) {
    MPI_Pcontrol(-1, timername.c_str());
  }
  addstat(elapsed(endtime, starttime), b);
}

// Reset the timer
void Timer::reset() { resetstats(); }

// Timer name
string Timer::name() const { return timername; }

// Output timer data
void Timer::outputData(ostream &os) const {
  double avg, stddev, bavg, bstddev;
  if (count == 0.0) {
    avg = 0.0;
    stddev = 0.0;
    bavg = 0.0;
    bstddev = 0.0;
  } else {
    avg = wtime / count;
    stddev = sqrt(max(0.0, wtime2 / count - pow(avg, 2)));
    bavg = bytes / count;
    bstddev = sqrt(max(0.0, bytes2 / count - pow(bavg, 2)));
  }

  os << timername << ":"
     << " cnt: " << count << "   time: sum: " << seconds_per_tick() * wtime
     << " avg: " << seconds_per_tick() * avg
     << " stddev: " << seconds_per_tick() * stddev
     << " min: " << seconds_per_tick() * wmin
     << " max: " << seconds_per_tick() * wmax << "   bytes: sum: " << bytes
     << " avg: " << bavg << " stddev: " << bstddev << " min: " << bmin
     << " max: " << bmax << eol;
}

// Fortran wrappers
extern "C" {

// In Fortran, a timer should be declared and used like this:
//
//    ! Save the timer handle, and initialise it to zero; this
//    ! ensures that the timer is created only once:
//    CCTK_POINTER, save :: timer = 0
//    call Timer_create (timer, "Name")
//
//    ! Start the timer:
//    call Timer_start (timer)
//
//    ! Stop the timer, and pass the number of bytes:
//    CCTK_REAL :: bytes
//    bytes = ...
//    call Timer_stop (timer, bytes)

void CCTK_FCALL CCTK_FNAME(Timer_create)(CCTK_POINTER *timer,
                                         ONE_FORTSTRING_ARG) {
  if (*timer != 0)
    return; // create the timer only once
  ONE_FORTSTRING_CREATE(timername);
  *timer = new Timer(timername);
  free(timername);
}

void CCTK_FCALL CCTK_FNAME(Timer_destroy)(CCTK_POINTER *timer) {
  if (*timer == 0)
    return; // delete the timer only if it has been created
  delete (Timer *)*timer;
  *timer = 0;
}

void CCTK_FCALL CCTK_FNAME(Timer_start)(CCTK_POINTER *timer) {
  assert(*timer != 0);
  ((Timer *)*timer)->start();
}

void CCTK_FCALL CCTK_FNAME(Timer_stop)(CCTK_POINTER *timer,
                                       CCTK_REAL const *b) {
  assert(*timer != 0);
  ((Timer *)*timer)->stop(*b);
}

void CCTK_FCALL CCTK_FNAME(Timer_reset)(CCTK_POINTER *timer) {
  assert(*timer != 0);
  ((Timer *)*timer)->reset();
}

} // extern "C"

extern "C" {
void CarpetLib_printtimestats(CCTK_ARGUMENTS);
}

void CarpetLib_printtimestats(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static int next_output = 0;
  if (print_timestats_every > 0 and cctk_iteration >= next_output) {
    next_output = cctk_iteration + print_timestats_every;

    ostringstream filenamebuf;
    filenamebuf << out_dir << "/" << timestat_file << "." << setw(4)
                << setfill('0') << dist::rank() << ".txt";
    string const filename = filenamebuf.str();

    ofstream file;
    static bool do_truncate = true;
    if (do_truncate) {
      if (not IO_TruncateOutputFiles(cctkGH)) {
        do_truncate = false;
      }
    }
    if (do_truncate) {
      do_truncate = false;
      file.open(filename.c_str(), ios::out | ios::trunc);
    } else {
      file.open(filename.c_str(), ios::out | ios::app);
    }

    static bool do_print_info = true;
    if (do_print_info) {
      do_print_info = false;
      if (CCTK_IsFunctionAliased("UniqueBuildID")) {
        char const *const build_id =
            static_cast<char const *>(UniqueBuildID(cctkGH));
        file << "Build ID: " << build_id << eol;
      }
      if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
        char const *const sim_id =
            static_cast<char const *>(UniqueSimulationID(cctkGH));
        file << "Simulation ID: " << sim_id << eol;
      }
      file << "Running with " << dist::size() << " processes and "
           << dist::total_num_threads() << " threads" << eol;
    } // if do_print_info

    if (not timerSet)
      timerSet = new TimerSet;
    file << "******************************************************************"
            "**************"
         << eol << "CarpetLib timing information at iteration "
         << cctkGH->cctk_iteration << " time " << cctkGH->cctk_time << ":"
         << eol << *timerSet;

    file.close();

  } // if print_timestats
}

} // namespace CarpetLib
