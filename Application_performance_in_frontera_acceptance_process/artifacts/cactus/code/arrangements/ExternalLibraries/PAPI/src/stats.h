#ifndef STATS_H
#define STATS_H

#include <cctk.h>
#include <cctk_Parameters.h>

#include <stdarg.h>
#include <stdio.h>

#define num_threads PAPI_Cactus_num_threads
#define eventsets PAPI_Cactus_eventsets
#define num_events PAPI_Cactus_num_events
#define events PAPI_Cactus_events

// Number of threads
extern int num_threads;
// Eventsets for each thread
extern int *restrict eventsets; // [num_threads]

// Number of events in these eventsets
extern int num_events;
// Events in these eventsets
extern int *restrict events; // [num_events]

// Output a a debug message
#define outinfo PAPI_Cactus_outinfo
void outinfo(const char *const function);

// Check for an error, and output an error message if there was one
#define chkerr PAPI_Cactus_chkerr
void chkerr(const int ierr, const char *const function, ...)
    CCTK_ATTRIBUTE_FORMAT(printf, 2, 3);

#endif // #ifndef STATS_H
