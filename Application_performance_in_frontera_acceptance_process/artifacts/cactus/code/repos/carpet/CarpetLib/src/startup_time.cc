#include <cctk.h>

#include <cstdio>
#include <cstdlib>
#include <string>

// IRIX wants this before <time.h>
#if HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#elif HAVE_TIME_H
#include <time.h>
#endif
#endif

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "startup_time.hh"

namespace CarpetLib {

using namespace std;

// Return the current wall time
static double get_walltime() {
#ifdef HAVE_TIME_GETTIMEOFDAY
  // get the current time
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1.0e6;
#else
  return 0.0;
#endif
}

void output_startup_time() {
  char *const cactus_starttime = getenv("CACTUS_STARTTIME");
  if (not cactus_starttime) {
    CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not determine Cactus startup time (environment variable "
               "CACTUS_STARTTIME is not set; it should be set to the output of "
               "\"date +%%s\")");
    return;
  }

  double starttime;
  int const iret = sscanf(cactus_starttime, "%lf", &starttime);
  if (iret != 1) {
    CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not determine Cactus startup time (environment variable "
               "CACTUS_STARTTIME has illegal value \"%s\"; it should instead "
               "be set to the output of \"date +%%s\", which is a single "
               "number)",
               cactus_starttime);
    return;
  }

  double const currenttime = get_walltime();
  double const startuptime = currenttime - starttime;

  CCTK_VInfo(CCTK_THORNSTRING, "Process startup time was %.3g seconds",
             startuptime);
}

} // namespace CarpetLib
