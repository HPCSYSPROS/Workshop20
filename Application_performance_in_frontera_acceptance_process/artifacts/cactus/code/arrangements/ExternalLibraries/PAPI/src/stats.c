#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <assert.h>
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

// For backward compatibility
#if PAPI_VERSION_MAJOR(PAPI_VERSION) <= 4

// Taken (and modified) from PAPI 5.1.0.1, file src/papi.c
static int PAPI_add_named_event(int EventSet, char *EventName) {
  int ret;
  int code;
  ret = PAPI_event_name_to_code(EventName, &code);
  if (ret != PAPI_OK)
    return ret;
  ret = PAPI_add_event(EventSet, code);
  return ret;
}
#endif

////////////////////////////////////////////////////////////////////////////////

void outinfo(const char *const function) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
#pragma omp master
    { printf("[%s]\n", function); }
  }
}

// Check for an error, and output an error message if there was one
void chkerr(const int ierr, const char *const function, ...) {
  if (ierr < 0) {
#pragma omp critical
    {
      va_list ap;
      va_start(ap, function);
      printf("ERROR %d [%s] in ", ierr, PAPI_strerror(ierr));
      vprintf(function, ap);
      printf("\n");
      va_end(ap);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

int num_threads;
int *restrict eventsets = NULL;
int num_events;
int *restrict events;

// The PAPI component of these eventsets
static int component = -1;

////////////////////////////////////////////////////////////////////////////////

// Remember timer values when PAPI counters were last reset
static long long *restrict last_nsec = NULL;
static long long *restrict last_cyc = NULL;

#if defined __linux__

// This is recommended for Linux, but e.g. doesn't work on OSX
#include <pthread.h>
static unsigned long thread_id(void) {
  const pthread_t tid = pthread_self();
  unsigned long ret;
  memcpy(&ret, &tid, sizeof ret);
  return ret;
}

#elif defined _OPENMP

// This works only if no threads are created/destroyed
// (if the number of threads doesn't change)
#include <omp.h>
static unsigned long thread_id(void) { return omp_get_thread_num(); }

#else

// Don't know how to obtain a thread id
static unsigned long thread_id(void) { return 0; }

#endif

void PAPI_init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  CCTK_INFO("Initialising PAPI");

  // OpenMP

  assert(per_thread_statistics);
  num_threads = omp_get_max_threads();

  // Initialise PAPI

  outinfo("PAPI_library_init");
  ierr = PAPI_library_init(PAPI_VER_CURRENT);
  chkerr(ierr, "PAPI_library_init");
  if (ierr != PAPI_VER_CURRENT) {
    CCTK_ERROR("PAPI library version mismatch");
  }

  if (per_thread_statistics) {
    outinfo("PAPI_thread_init");
// TODO: Investigate PAPI_INHERIT -- but this may work only with
// recent version of PAPI and/or Linux

// NOTE: Our implementation of thread_id probably cannot handle
// changes in the number of threads. That is, omp_set_num_threads
// should not be called after this point.

#if 0
    // Ensure all threads are created
    // (This is necessary for Intel's OpenMP.)
    volatile int dummy = 0;
#pragma omp parallel reduction(+ : dummy)
    {
      dummy += 1;
    }
#endif

    ierr = PAPI_thread_init(thread_id);
    chkerr(ierr, "PAPI_thread_init");
  }

  if (use_multiplexing) {
    outinfo("PAPI_multiplex_init");
    ierr = PAPI_multiplex_init();
    chkerr(ierr, "PAPI_multiplex_init");
  }

  // Get some basic information

  outinfo("PAPI_num_components");
  int num_components = ierr = PAPI_num_components();
  chkerr(ierr, "PAPI_num_components");
  if (ierr < 0)
    num_components = 0;
  if (verbose) {
    printf("There are %d PAPI components\n", num_components);
  }

  if (verbose) {
    for (int n = 0; n < num_components; ++n) {
      outinfo("PAPI_get_component_info");
      const PAPI_component_info_t *const info = PAPI_get_component_info(n);
      printf("PAPI component #%d:\n", n);
      printf("   name: %s\n", info->name);
#if PAPI_VERSION_MAJOR(PAPI_VERSION) >= 5
      printf("   short_name: %s\n", info->short_name);
      printf("   description: %s\n", info->description);
#endif
      printf("   num_cntrs: %d\n", info->num_cntrs);
      printf("   num_mpx_cntrs: %d\n", info->num_mpx_cntrs);
      printf("   num_preset_events: %d\n", info->num_preset_events);
      printf("   num_native_events: %d\n", info->num_native_events);
    }
  }

  component = 0;
  if (verbose) {
    printf("Using PAPI component %d\n", component);
  }

  outinfo("PAPI_num_counters");
  int num_counters = ierr = PAPI_num_counters();
  chkerr(ierr, "PAPI_num_counters");
  if (ierr < 0)
    num_counters = 0;
  if (verbose) {
    printf("There are %d PAPI counters\n", num_counters);
  }

  // Translate user event names to an eventset

  eventsets = malloc(num_threads * sizeof *eventsets);
  assert(eventsets);

#pragma omp parallel private(ierr)
  {
    const int thread_num = omp_get_thread_num();

    outinfo("PAPI_create_eventset");
    eventsets[thread_num] = PAPI_NULL;
    ierr = PAPI_create_eventset(&eventsets[thread_num]);
    chkerr(ierr, "PAPI_create_eventset");
    outinfo("PAPI_assign_eventset_component");
    ierr = PAPI_assign_eventset_component(eventsets[thread_num], component);
    chkerr(ierr, "PAPI_assign_eventset_component");
    if (use_multiplexing) {
      outinfo("PAPI_set_multiplex");
      ierr = PAPI_set_multiplex(eventsets[thread_num]);
      chkerr(ierr, "PAPI_set_multiplex");
    }

    // User-defined event sets
    const int num_eventsets = 5;
    const char *restrict const eventset_descriptions[] = {
        "flops", "ipc", "icache", "dcache", "memory",
    };
    assert(sizeof eventset_descriptions / sizeof *eventset_descriptions ==
           num_eventsets);
    const char *const eventset_eventnames[] = {
        events_flops, events_ipc, events_icache, events_dcache, events_memory,
    };
    assert(sizeof eventset_eventnames / sizeof *eventset_eventnames ==
           num_eventsets);

    for (int n = 0; n < num_eventsets; ++n) {
#pragma omp master
      if (verbose) {
        printf("Adding eventset %s:\n", eventset_descriptions[n]);
      }

      char *const eventnames = strdup(eventset_eventnames[n]);
      const char *const sep = ", ";
      char *lasts;
      for (char *event_name = strtok_r(eventnames, sep, &lasts); event_name;
           event_name = strtok_r(NULL, sep, &lasts)) {
#pragma omp master
        if (verbose) {
          printf("   event %s\n", event_name);
        }
        outinfo("PAPI_add_named_event");
        ierr = PAPI_add_named_event(eventsets[thread_num], event_name);
        chkerr(ierr, "PAPI_add_named_event[%s]", event_name);
      }
      free(eventnames);
    }
  }

  outinfo("PAPI_num_events");
  num_events = ierr = PAPI_num_events(eventsets[0]);
  chkerr(ierr, "PAPI_num_events");
  if (ierr < 0)
    num_events = 0;

  outinfo("PAPI_list_events");
  events = malloc(num_events * sizeof *events);
  ierr = PAPI_list_events(eventsets[0], events, &num_events);
  chkerr(ierr, "PAPI_list_events");

  last_nsec = malloc(num_threads * sizeof *last_nsec);
  assert(last_nsec);
  last_cyc = malloc(num_threads * sizeof *last_cyc);
  assert(last_cyc);

#pragma omp parallel private(ierr)
  {
    const int thread_num = omp_get_thread_num();

    outinfo("PAPI_start");
    ierr = PAPI_start(eventsets[thread_num]);
    chkerr(ierr, "PAPI_start");

    outinfo("PAPI_get_real_nsec");
    last_nsec[thread_num] = PAPI_get_real_nsec();
    outinfo("PAPI_get_real_cyc");
    last_cyc[thread_num] = PAPI_get_real_cyc();
  }
}

static void output_stats(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  CCTK_INFO("PAPI Statistics:");

  long long elapsed_nsec = 0;
  long long elapsed_cyc = 0;
  long long values[num_events];
  for (int i = 0; i < num_events; ++i) {
    values[i] = 0;
  }
#pragma omp parallel private(ierr)
  {
    const int thread_num = omp_get_thread_num();

    outinfo("PAPI_get_real_nsec");
    const long long my_nsec = PAPI_get_real_nsec();

    outinfo("PAPI_read_ts");
    long long my_cyc;
    long long my_values[num_events];
    ierr = PAPI_read_ts(eventsets[thread_num], my_values, &my_cyc);
    chkerr(ierr, "PAPI_read_ts");

    const long long my_elapsed_nsec = my_nsec - last_nsec[thread_num];
    const long long my_elapsed_cyc = my_cyc - last_cyc[thread_num];
    last_nsec[thread_num] = my_nsec;
    last_cyc[thread_num] = my_cyc;

#pragma omp critical(PAPI_output_stats_add)
    {
      elapsed_nsec += my_elapsed_nsec;
      elapsed_cyc += my_elapsed_cyc;
      for (int i = 0; i < num_events; ++i) {
        values[i] += my_values[i];
      }
    }
  }

  printf("   CPU time          %20.9f sec\n", elapsed_nsec / 1.0e+9);
  printf("   CPU cycles        %20.9f Gcyc\n", elapsed_cyc / 1.0e+9);
  for (int i = 0; i < num_events; ++i) {
    outinfo("PAPI_event_code_to_name");
    char eventname[PAPI_MAX_STR_LEN];
    ierr = PAPI_event_code_to_name(events[i], eventname);
    chkerr(ierr, "PAPI_event_code_to_name");
    printf("   %-15s   %20.9f cyc^(-1)\n", eventname,
           1.0 * values[i] / elapsed_cyc);
  }
}

void PAPI_output_stats_analysis(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_every == -1)
    return;
  if (cctk_iteration > 0) {
    if (out_every == 0)
      return;
    if (cctk_iteration % out_every > 0)
      return;
  }

  output_stats(CCTK_PASS_CTOC);
}

void PAPI_output_stats_terminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_every == -1)
    return;
  if (out_every > 0)
    return;

  output_stats(CCTK_PASS_CTOC);
}
