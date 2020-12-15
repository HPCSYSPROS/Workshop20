#ifndef CYCLECLOCK_H
#define CYCLECLOCK_H

/* This defines:
 *    typedef XXX ticks;
 *    ticks getticks();
 *    double elapsed(ticks t1, ticks t0);
 * Use as:
 *    #include <cycleclock.h>
 *    ticks t0 = getticks();
 *    ...
 *    ticks t1 = getticks();
 *    double elapsed_ticks = elapsed(t1, t0);
 *    double elapsed_seconds = seconds_per_tick() * elapsed_ticks;
*/

#include <cctk.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "cycle.h"

#ifdef HAVE_TICK_COUNTER
double seconds_per_tick(void);
void measure_tick(void);
#else
#warning "tick counter not available"
#endif

#ifdef __cplusplus
}
#endif

#endif /* CYCLECLOCK_H */
