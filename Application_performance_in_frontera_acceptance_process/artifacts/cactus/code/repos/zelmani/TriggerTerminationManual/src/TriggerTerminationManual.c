#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_String.h"

void TriggerTerminationManual_SayAllGood(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char filename[10000];
  FILE *outfile;

  /* only one processor needs to check the file */
  if (CCTK_MyProc (cctkGH) != 0)
  {
    return;
  }

  /* add IO::out_dir to filename */
  snprintf (filename, sizeof filename,
		 "%s/%s", out_dir, "allgood");
  outfile = fopen(filename, "w");
  if (outfile != NULL)
  {
    fprintf(outfile,"1\n");
    fclose(outfile);
  }
  else
  {
    CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not open file '%s' to write termination code: %s",
                filename, strerror(errno));
  }
}

void TriggerTerminationManual_StartTimer (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int TimerIndex;

  /* only one processor needs to query the elapsed runtime */
  if (CCTK_MyProc (cctkGH) != 0)
  {
    return;
  }

  /* Create timer */
  TimerIndex = CCTK_TimerCreate("WatchWalltime");

  /* Start timer */
  CCTK_TimerStart("WatchWalltime");

  *watchminutes = output_remtime_every_minutes;

  CCTK_VInfo (CCTK_THORNSTRING, "Started Timer");
  CCTK_VInfo (CCTK_THORNSTRING,
              "Reminding you every %d minutes about remaining walltime",
	      (int) output_remtime_every_minutes);
}



void TriggerTerminationManual_ResetMinutes (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *watchminutes = output_remtime_every_minutes;
}



void TriggerTermination_CheckWalltime (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  cTimerData *info;
  const cTimerVal *walltime;
  CCTK_REAL time;

  /* only one processor needs to query the elapsed runtime */
  if (CCTK_MyProc (cctkGH) != 0)
  {
    return;
  }

  info = CCTK_TimerCreateData();
  CCTK_Timer("WatchWalltime",info);

  /* stop timer */
  CCTK_TimerStop("WatchWalltime");
  /* get walltime */
  walltime = CCTK_GetClockValue("gettimeofday",info);
  time = CCTK_TimerClockSeconds(walltime);
  CCTK_TimerDestroyData(info);
  /* start timer again */
  CCTK_TimerStart("WatchWalltime");


  if ( (time/60.0 > *watchminutes) && *watchminutes != 0) {
    *watchminutes = (*watchminutes)+output_remtime_every_minutes;
    CCTK_INFO ("***********************************************************");
    CCTK_VInfo (CCTK_THORNSTRING,
                "Remaining wallclock time for your job is %1.2f minutes",
                (double) (max_walltime*60.0-time/60.0));
    CCTK_INFO ("***********************************************************");
  } 
  

  if (time/60.0e0 >= (max_walltime*60.0 - on_remaining_walltime)) {
    CCTK_VInfo (CCTK_THORNSTRING, 
                "Remaining wallclock time for your job is %1.2f minutes.  "
                "Triggering termination...",
                (double) (max_walltime*60.0-time/60.0));
    CCTK_TerminateNext (cctkGH);
  }
}
