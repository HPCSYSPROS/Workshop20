#include <cassert>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <list>
#include <string>
#include <utility>
#include <vector>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_String.h>

#if HAVE_UNISTD_H
#include <fcntl.h>
#include <unistd.h>
#endif

#include <defs.hh>

#include "CactusTimer.hh"
#include "CactusTimerSet.hh"

namespace Timers {

using namespace std;

// Create a new Cactus timer with the given name
CactusTimer::CactusTimer(string timername) : running(false) {
  handle = CCTK_TimerCreate(timername.c_str());
  assert(handle >= 0);

  timerSet.add(this);
  msgCreate();
}

// Destroy a timer
CactusTimer::~CactusTimer() {
  timerSet.remove(this);
  check(not CCTK_TimerDestroyI(handle));
}

// Start the timer
void CactusTimer::start() {
  msgStart();
  running = true;
  CCTK_TimerStartI(handle);
}

// Stop the timer
void CactusTimer::stop() {
  CCTK_TimerStopI(handle);
  running = false;
  msgStop();
}

// Reset the timer
void CactusTimer::reset() { CCTK_TimerResetI(handle); }

// Timer name
string CactusTimer::name() const {
  const char *const timername = CCTK_TimerName(handle);
  assert(timername);
  return string(timername);
}

double CactusTimer::getTime() {
  DECLARE_CCTK_PARAMETERS;

  const bool was_running = running;
  if (was_running)
    stop();

  static cTimerData *timer = 0;
  if (not timer)
    timer = CCTK_TimerCreateData();
  assert(timer);

  CCTK_TimerI(handle, timer);
  const cTimerVal *tv = CCTK_GetClockValue(xml_clock, timer);
  double val;
  if (not tv) {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Clock \"%s\" not found for timer #%d \"%s\"", xml_clock, handle,
               CCTK_TimerName(handle));
    val = -1.0;
  } else {
    val = CCTK_TimerClockSeconds(tv);
  }
  msgRead(val);
  if (was_running)
    start();

  return val;
}

void CactusTimer::getGlobalTime(double &avg, double &max) {
  const cGH *const cctkGH = 0;

  int ierr;

  static int op_sum = -1;
  static int op_max = -1;
  if (op_sum < 0)
    op_sum = CCTK_ReductionArrayHandle("sum");
  if (op_max < 0)
    op_max = CCTK_ReductionArrayHandle("maximum");

  const double val = getTime();
  const CCTK_REAL val1 = val;

  CCTK_REAL sum1;
  ierr = CCTK_ReduceLocScalar(cctkGH, -1, op_sum, &val1, &sum1,
                              CCTK_VARIABLE_REAL);
  if (ierr) {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error in sum reduction");
  }
  avg = sum1 / CCTK_nProcs(cctkGH);

  CCTK_REAL max1;
  ierr = CCTK_ReduceLocScalar(cctkGH, -1, op_max, &val1, &max1,
                              CCTK_VARIABLE_REAL);
  if (ierr) {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error in maximum reduction");
  }
  max = max1;
}

vector<pair<string, string> > CactusTimer::getAllTimerNames() const {
  DECLARE_CCTK_PARAMETERS;

  static cTimerData *timer = NULL;
  if (not timer)
    timer = CCTK_TimerCreateData();
  assert(timer);

  CCTK_TimerI(handle, timer);

  vector<pair<string, string> > names(timer->n_vals);
  for (int i = 0; i < timer->n_vals; ++i) {
    names[i].first = timer->vals[i].heading;
    names[i].second = timer->vals[i].units;
  }

  return names;
}

vector<double> CactusTimer::getAllTimerValues() {
  DECLARE_CCTK_PARAMETERS;

  const bool was_running = running;
  if (was_running)
    stop();

  static cTimerData *timer = NULL;
  if (not timer)
    timer = CCTK_TimerCreateData();
  assert(timer);

  CCTK_TimerI(handle, timer);

  vector<double> vals(timer->n_vals);
  for (int i = 0; i < timer->n_vals; ++i) {
    switch (timer->vals[i].type) {
    case val_int:
      vals[i] = timer->vals[i].val.i;
      break;
    case val_long:
      vals[i] = timer->vals[i].val.l;
      break;
    case val_double:
      vals[i] = timer->vals[i].val.d;
      break;
    default:
      assert(0);
    }
  }

  if (was_running)
    start();

  return vals;
}

// Print timer data
void CactusTimer::printData() {
  const bool was_running = running;
  if (was_running)
    stop();

#if 0
    check (not CCTK_TimerPrintDataI (handle, -1)); // -1 means: all clocks
#endif

  static cTimerData *timer = 0;
  if (not timer)
    timer = CCTK_TimerCreateData();
  assert(timer);
  CCTK_TimerI(handle, timer);

  static bool firsttime = true;
  if (firsttime) {
    printf("# 1: timer name");
    for (int i = 0; i < timer->n_vals; ++i) {
      printf(" %d: %s [%s]", i + 2, timer->vals[i].heading,
             timer->vals[i].units);
    }
    printf("\n");
    firsttime = false;
  }

  printf("%s:", name().c_str());
  for (int i = 0; i < timer->n_vals; ++i) {
    switch (timer->vals[i].type) {
    case val_int:
      printf(" %d", timer->vals[i].val.i);
      break;
    case val_long:
      printf(" %ld", timer->vals[i].val.l);
      break;
    case val_double:
      printf(" %g", timer->vals[i].val.d);
      break;
    case val_none:
      break;
    default:
      assert(0);
    }
  }
  printf("\n");

  if (was_running)
    start();
}

// Output (debug) messages that a timer is starting or stopping
void CactusTimer::msgCreate() const {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Timer \"%s\" created", name().c_str());
  }
}

void CactusTimer::msgStart() const {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Timer \"%s\" starting", name().c_str());
  }
}

void CactusTimer::msgStop() const {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Timer \"%s\" stopping", name().c_str());
  }
}

void CactusTimer::msgRead(double val) const {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Timer \"%s\" read: %g", name().c_str(), val);
  }
}

ostream &CactusTimer::serialise(ostream &os) {
  os << scientific << setprecision(19) << getTime() << " " << name();
  return os;
}

} // namespace Carpet
