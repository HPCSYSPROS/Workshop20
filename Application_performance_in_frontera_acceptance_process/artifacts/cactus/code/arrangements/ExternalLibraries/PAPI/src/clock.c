#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_max_threads(void) { return 1; }
static int omp_get_thread_num(void) { return 0; }
#endif

#include <papi.h>

#include "stats.h"

static const char *const clock_name = "PAPI";

static int num_event_vals = -1;
static int num_clock_vals = -1;

static const int counter_nsec = 0;
static const int counter_cyc = 1;
static const int counter_count = 2;
static const int counter_begin = 3;

static int num_statistics_vals = -1;
static int num_cactus_clock_vals = -1;
static const char **counter_names = NULL; // [num_cactus_clock_vals]
static const char **counter_units = NULL; // [num_cactus_clock_vals]

typedef long long papi_counter_t;

typedef struct {
  // indexing: [num_clock_vals], laid out as [num_threads][num_event_vals]
  papi_counter_t *restrict accum;    // [num_clock_vals]
  papi_counter_t *restrict snapshot; // [num_clock_vals]
  papi_counter_t count;
  bool running;
} papi_clock_t;

static void read_papi_counters(papi_counter_t *restrict const counters) {
  int ierr;

/* for (int i=0; i<num_clock_vals; ++i) { */
/*   counters[i] = -1; */
/* } */
#pragma omp barrier
#pragma omp parallel private(ierr)
  {
    const int thread_num = omp_get_thread_num();

    papi_counter_t *restrict const my_counters =
        &counters[num_event_vals * thread_num];

    my_counters[counter_nsec] = PAPI_get_real_nsec();
    ierr = PAPI_read_ts(eventsets[thread_num], &my_counters[counter_begin],
                        &my_counters[counter_cyc]);
    // chkerr(ierr, "PAPI_read_ts");
    if (ierr < 0)
      my_counters[counter_cyc] = PAPI_get_real_cyc();
    ++my_counters[counter_count];
  }
  /* for (int i=0; i<num_clock_vals; ++i) { */
  /*   assert(counters[i] != -1); */
  /* } */
  /* for (int i=0; i<num_clock_vals; ++i) { */
  /*   assert(counters[i] >= 0); */
  /* } */
}

static void *clock_create(const int timernum) {
  papi_clock_t *restrict const data = malloc(sizeof *data);
  assert(data);
  data->accum = malloc(num_clock_vals * sizeof *data->accum);
  assert(data->accum);
  data->snapshot = malloc(num_clock_vals * sizeof *data->snapshot);
  assert(data->snapshot);
  for (int i = 0; i < num_clock_vals; ++i) {
    data->accum[i] = 0;
  }
  data->running = false;
  return data;
}

static void clock_destroy(const int timernum, void *restrict const data_) {
  papi_clock_t *restrict const data = data_;
  assert(data);
  free(data->accum);
  free(data->snapshot);
  free(data);
}

static void clock_start(const int timernum, void *restrict const data_) {
  papi_clock_t *restrict const data = data_;
  assert(!data->running);
  read_papi_counters(data->snapshot);
  data->running = true;
}

static void clock_stop(const int timernum, void *restrict const data_) {
  papi_clock_t *restrict const data = data_;
  // Apparently clocks can be created in a "running" state?
  if (!data->running)
    return;
  assert(data->running);
  papi_counter_t current[num_clock_vals];
  read_papi_counters(current);
  data->running = false;
  for (int i = 0; i < num_clock_vals; ++i) {
    data->accum[i] += current[i] - data->snapshot[i];
  }
  /* for (int i=0; i<num_clock_vals; ++i) { */
  /*   assert(data->accum[i] >= 0); */
  /* } */
}

static void clock_reset(const int timernum, void *restrict const data_) {
  papi_clock_t *restrict const data = data_;
  for (int i = 0; i < num_clock_vals; ++i) {
    data->accum[i] = 0;
  }
}

static void clock_get(const int timernum, void *restrict const data_,
                      cTimerVal *restrict const vals) {
  DECLARE_CCTK_PARAMETERS;

  const papi_clock_t *restrict const data = data_;

  int count = 0;
  for (int i = 0; i < num_event_vals; ++i) {

    double sum = 0.0, sum2 = 0.0, min = DBL_MAX, max = 0.0;
    for (int j = i; j < num_clock_vals; j += num_event_vals) {
      /* assert(data->accum[j] >= 0); */
      const double val = data->accum[j];
      sum += val;
      sum2 += val * val;
      min = fmin(min, val);
      max = fmax(max, val);
    }
    /* assert(sum >= 0.0); */

    if (!define_statistics_clocks) {

      vals[count].type = val_double;
      vals[count].heading = counter_names[count];
      vals[count].units = counter_units[count];
      vals[count].val.d = sum / 1.0e+9;
      vals[count].seconds = sum / 1.0e+9;
      vals[count].resolution = 1.0e-9;
      ++count;

    } else {

      const double avg = sum / num_threads;
      const double sdv = sqrt(fmax(0.0, sum2 / num_threads - avg * avg));

      vals[count].type = val_double;
      vals[count].heading = counter_names[count];
      vals[count].units = counter_units[count];
      vals[count].val.d = avg / 1.0e+9;
      vals[count].seconds = avg / 1.0e+9;
      vals[count].resolution = 1.0e-9;
      ++count;

      vals[count].type = val_double;
      vals[count].heading = counter_names[count];
      vals[count].units = counter_units[count];
      vals[count].val.d = sdv / 1.0e+9;
      vals[count].seconds = sdv / 1.0e+9;
      vals[count].resolution = 1.0e-9;
      ++count;

      vals[count].type = val_double;
      vals[count].heading = counter_names[count];
      vals[count].units = counter_units[count];
      vals[count].val.d = min / 1.0e+9;
      vals[count].seconds = min / 1.0e+9;
      vals[count].resolution = 1.0e-9;
      ++count;

      vals[count].type = val_double;
      vals[count].heading = counter_names[count];
      vals[count].units = counter_units[count];
      vals[count].val.d = max / 1.0e+9;
      vals[count].seconds = max / 1.0e+9;
      vals[count].resolution = 1.0e-9;
      ++count;
    }
  }
}

static void clock_set(const int timernum, void *restrict const data_,
                      cTimerVal *restrict const vals) {
  papi_clock_t *restrict const data = data_;

  for (int i = 0; i < num_event_vals; ++i) {
    assert(vals[i].type == val_double);
    const papi_counter_t val = llrint(vals[i].val.d * 1.0e+9);
    for (int j = i; j < num_clock_vals; j += num_event_vals) {
      data->accum[j] = val;
    }
  }
}

static void add_clock(const char **names, const char **units, const char *name,
                      const char *unit, int *clock, int nclocks) {
  DECLARE_CCTK_PARAMETERS;
  if (!define_statistics_clocks) {
    assert(*clock < nclocks);
    names[*clock] = strdup(name);
    units[*clock] = unit;
    ++*clock;
  } else {
    const char *const suffixes[] = {"_avg", "_sdv", "_min", "_max"};
    for (int s = 0; s < 4; ++s) {
      assert(*clock < nclocks);
      char buf[1000];
      snprintf(buf, sizeof buf, "%s%s", name, suffixes[s]);
      names[*clock] = strdup(buf);
      units[*clock] = unit;
      ++*clock;
    }
  }
}

void PAPI_register_clock(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("PAPI_register_clock");

  int ierr;

  num_event_vals = counter_begin + num_events;

  // PAPI clocks

  num_clock_vals = num_threads * num_event_vals;

  // Cactus clocks

  num_statistics_vals = define_statistics_clocks ? 4 : 1;
  num_cactus_clock_vals = num_event_vals * num_statistics_vals;

  counter_names = malloc(num_cactus_clock_vals * sizeof *counter_names);
  assert(counter_names);
  counter_units = malloc(num_cactus_clock_vals * sizeof *counter_units);
  assert(counter_units);

  int clock = 0;
  add_clock(counter_names, counter_units, "P_time", "sec", &clock,
            num_cactus_clock_vals);
  add_clock(counter_names, counter_units, "P_cycles", "Gcycles", &clock,
            num_cactus_clock_vals);
  add_clock(counter_names, counter_units, "P_count", "Gevents", &clock,
            num_cactus_clock_vals);
  for (int i = 0; i < num_events; ++i) {
    char eventname[PAPI_MAX_STR_LEN + 4];
    ierr = PAPI_event_code_to_name(events[i], eventname);
    chkerr(ierr, "PAPI_event_code_to_name");
    if (strlen(eventname) >= 4) {
      // Change prefix "PAPI" to "P" to save space
      for (char *p = eventname + 1; p[2]; ++p)
        p[0] = p[3];
    }
    add_clock(counter_names, counter_units, eventname, "Gevents", &clock,
              num_cactus_clock_vals);
  }
  assert(clock == num_cactus_clock_vals);

  cClockFuncs funcs;
  funcs.name = clock_name;
  funcs.n_vals = num_cactus_clock_vals;
  funcs.create = clock_create;
  funcs.destroy = clock_destroy;
  funcs.start = clock_start;
  funcs.stop = clock_stop;
  funcs.reset = clock_reset;
  funcs.get = clock_get;
  funcs.set = clock_set;
  funcs.seconds = NULL;
  ierr = CCTK_ClockRegister(clock_name, &funcs);
  assert(ierr >= 0);
}
