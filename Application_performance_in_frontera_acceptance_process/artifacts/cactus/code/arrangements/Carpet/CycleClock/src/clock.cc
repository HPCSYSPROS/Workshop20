#include "cycleclock.h"

#ifdef HAVE_TICK_COUNTER

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace CycleClock {

using namespace std;

class cycleclock_t {
  double sum;
  double sum2;
  double min;
  double max;
  double count;

  ticks last;

public:
  cycleclock_t() {
    reset();
    start();
  }

  ~cycleclock_t() {}

  void start() { last = getticks(); }

  void stop() {
    ticks const current = getticks();
    double const difference = elapsed(current, last);
    sum += difference;
    sum2 += pow(difference, 2.0);
    min = min == 0.0 ? difference : fmin(min, difference);
    max = fmax(min, difference);
    count += 1.0;
  }

  void reset() {
    sum = 0.0;
    sum2 = 0.0;
    min = 0.0; // numeric_limits<double>::max();
    max = 0.0;
    count = 0.0;
  }

  void get(cTimerVal *restrict const vals) const {
    double const tick = seconds_per_tick();

    // Sum time
    vals[0].type = val_double;
    vals[0].heading = "cycle";
    vals[0].units = "secs";
    vals[0].val.d = sum;
    vals[0].seconds = tick * vals[0].val.d;
    vals[0].resolution = tick;

    // Average
    vals[1].type = val_double;
    vals[1].heading = "cycle[avg]";
    vals[1].units = "secs";
    vals[1].val.d = count == 0.0 ? 0.0 : sum / count;
    vals[1].seconds = tick * vals[1].val.d;
    vals[1].resolution = tick;

    // Standard deviation
    vals[2].type = val_double;
    vals[2].heading = "cycle[sdv]";
    vals[2].units = "secs";
    vals[2].val.d =
        (count == 0.0 ? 0.0 : sqrt(fabs(sum2 * count - pow(sum, 2.0)) / count));
    vals[2].seconds = tick * vals[2].val.d;
    vals[2].resolution = tick;

    // Minimum
    vals[3].type = val_double;
    vals[3].heading = "cycle[min]";
    vals[3].units = "secs";
    vals[3].val.d = min;
    vals[3].seconds = tick * vals[3].val.d;
    vals[3].resolution = tick;

    // Maximum
    vals[4].type = val_double;
    vals[4].heading = "cycle[max]";
    vals[4].units = "secs";
    vals[4].val.d = max;
    vals[4].seconds = tick * vals[4].val.d;
    vals[4].resolution = tick;
  }

  void set(cTimerVal const *restrict const vals) {
    reset(); // punt
    sum = vals[0].val.d;
  }
};

void *cycleclock_create(int const timernum) { return new cycleclock_t; }

void cycleclock_destroy(int const timernum, void *const data) {
  if (!data)
    return;
  delete static_cast<cycleclock_t *>(data);
}

void cycleclock_start(int const timernum, void *const data) {
  static_cast<cycleclock_t *>(data)->start();
}

void cycleclock_stop(int const timernum, void *const data) {
  static_cast<cycleclock_t *>(data)->stop();
}

void cycleclock_reset(int const timernum, void *const data) {
  static_cast<cycleclock_t *>(data)->reset();
}

void cycleclock_get(int const timernum, void *const data,
                    cTimerVal *const vals) {
  static_cast<cycleclock_t const *>(data)->get(vals);
}

void cycleclock_set(int const timernum, void *const data,
                    cTimerVal *const vals) {
  static_cast<cycleclock_t *>(data)->set(vals);
}

void cycleclock_register() {
  cClockFuncs functions;
  functions.n_vals = 5;
  functions.create = cycleclock_create;
  functions.destroy = cycleclock_destroy;
  functions.start = cycleclock_start;
  functions.stop = cycleclock_stop;
  functions.reset = cycleclock_reset;
  functions.get = cycleclock_get;
  functions.set = cycleclock_set;
  CCTK_ClockRegister("cycle", &functions);
}

extern "C" int CycleClock_Setup() {
  DECLARE_CCTK_PARAMETERS;

  measure_tick();
  if (register_clock) {
    cycleclock_register();
  }
  return 0;
}

} // namespace CycleClock

#else // HAVE_TICK_COUNTER

namespace CycleClock {
extern "C" int CycleClock_Setup() { return 0; }
}
#endif
