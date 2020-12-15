#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Termination.h"
#include "cctk_Timers.h"

void TerminationTrigger_StartTimer(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* only one processor needs to query the elapsed runtime */
  if (CCTK_MyProc(cctkGH) != 0) {
    return;
  }

  *watchminutes = output_remtime_every_minutes;

  CCTK_VInfo(CCTK_THORNSTRING,
             "Reminding you every %g minutes about remaining walltime",
             (double)output_remtime_every_minutes);
}

void TerminationTrigger_ResetMinutes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* only one processor needs to query the elapsed runtime */
  if (CCTK_MyProc(cctkGH) != 0) {
    return;
  }

  *watchminutes = output_remtime_every_minutes;
}

void TerminationTrigger_CheckWalltime(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL time;

  /* if the maximum wall time or the remaining wall time have not been
     set, then don't terminate */
  if (max_walltime == 0.0 || on_remaining_walltime == 0.0) {
    return;
  }

  /* only one processor needs to query the elapsed runtime */
  if (CCTK_MyProc(cctkGH) != 0) {
    return;
  }

  /* get walltime in seconds */
  time = CCTK_RunTime();

  if ((time / 60.0 > *watchminutes) && *watchminutes != 0) {
    *watchminutes += output_remtime_every_minutes;
    CCTK_INFO("***********************************************************");
    CCTK_VInfo(CCTK_THORNSTRING,
               "Remaining wallclock time for your job is %g minutes",
               (double)(max_walltime * 60.0 - time / 60.0));
    CCTK_INFO("***********************************************************");
  }

  if (time / 60.0 >= (max_walltime * 60.0 - on_remaining_walltime)) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Remaining wallclock time for your job is %g minutes.  "
               "Triggering termination...",
               (double)(max_walltime * 60.0 - time / 60.0));
    CCTK_TerminateNext(cctkGH);
  }
}
