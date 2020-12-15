 /*@@
   @file      InitialiseCactus.c
   @date      Fri Sep 18 14:04:02 1998
   @author    Tom Goodale
   @desc
              Responsible for doing any cactus specific initialisations
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <time.h>

#include "cctk_Flesh.h"
#include "cctk_Parameter.h"
#include "cctk_Schedule.h"
#include "cctk_WarnLevel.h"
#include "cctk_Misc.h"

#include "cctki_Banner.h"
#include "cctki_Bindings.h"
#include "cctki_Schedule.h"


static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_InitialiseCactus_c);

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
static int CCTKi_InitialiseScheduler (tFleshConfig *ConfigData);
static void CCTKi_RecoverParameters (tFleshConfig *config);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/
int CCTKBindings_RegisterThornFunctions (void);
int CCTKi_InitialiseSubsystemDefaults (void);
int CCTKi_BindingsParameterRecoveryInitialise (void);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/
static time_t startuptime;


 /*@@
   @routine    CCTKi_InitialiseCactus
   @date       Fri Sep 18 14:05:21 1998
   @author     Tom Goodale
   @desc
               Calls all CCTK-internal initialization routines.
   @enddesc
   @calls      CCTKi_InitialiseSubsystemDefaults
               CCTKi_ProcessEnvironment
               CCTKi_ProcessCommandLine
               CCTKi_CactusBanner
               CCTKi_InitialiseDataStructures
               CCTKi_ProcessParameterDatabase
               CCTKi_BindingsVariablesInitialise
               CCTKi_InitialiseScheduler
               CCTKi_CallStartupFunctions
               CCTKi_PrintBanners

   @var        argc
   @vdesc      The number of command line arguments
   @vtype      int *
   @vio        inout
   @endvar
   @var        argv
   @vdesc      The command line arguments
   @vtype      char ***
   @vio        inout
   @endvar
   @var        config
   @vdesc      Flesh configuration data
   @vtype      tFleshConfig *
   @vio        inout
   @endvar

   @returntype int
   @returndesc
               0  - success
   @endreturndesc
@@*/
int CCTKi_InitialiseCactus (int *argc, char ***argv, tFleshConfig *config)
{
  startuptime = time (NULL);

  CCTKi_InitialiseSubsystemDefaults ();

  CCTKi_ProcessEnvironment (argc, argv, config);

  CCTKi_ProcessCommandLine (argc, argv, config);

  CCTKi_CactusBanner ();

  CCTKi_InitialiseDataStructures (config);

  CCTKi_ProcessParameterDatabase (config);

  CCTKi_RecoverParameters (config);

  CCTKi_BindingsVariablesInitialise ();

  if (CCTKBindings_RegisterThornFunctions ())
  {
    CCTK_Warn (0, __LINE__, __FILE__, "Cactus",
               "Failed to register/provide aliased functions for active thorns");
  }

  CCTKi_InitialiseScheduler (config);

  CCTKi_CallStartupFunctions (config);

  CCTKi_PrintBanners ();

  return (0);
}


 /*@@
   @routine    CCTK_RunTime
   @date       Tue Oct 3 2000
   @author     Gabrielle Allen
   @desc
               Seconds since startup
   @enddesc

   @returntype int
   @returndesc
               The number of seconds since the run started.
   @endreturndesc
@@*/
int CCTK_RunTime (void)
{
  return ((int) (time (NULL) - startuptime));
}

void CCTK_FCALL CCTK_FNAME (CCTK_RunTime)
                           (int *runtime)
{
  *runtime = CCTK_RunTime ();
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
 /*@@
   @routine    CCTKi_InitialiseScheduler
   @date       Fri Sep 17 19:34:55 1999
   @author     Tom Goodale
   @desc
               Initialise all scheduled items
   @enddesc
   @calls      CCTKi_SetParameterSetMask
               CCTKi_BindingsParameterRecoveryInitialise
               CCTKi_BindingsScheduleInitialise
               CCTKi_DoScheduleSortAllGroups
               CCTK_SchedulePrint

   @var        config
   @vdesc      Flesh configuration data
   @vtype      tFleshConfig *
   @vio        out
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine CCTKi_DoScheduleSortAllGroups
   @endreturndesc
@@*/
static int CCTKi_InitialiseScheduler (tFleshConfig *config)
{
  int retval;
  const CCTK_INT *cctk_show_schedule;


  CCTKi_BindingsScheduleInitialise ();

  retval = CCTKi_DoScheduleSortAllGroups ();

  cctk_show_schedule = (const CCTK_INT *)
                       CCTK_ParameterGet ("cctk_show_schedule", "Cactus", NULL);

  if (*cctk_show_schedule)
  {
    puts ("-------------------------------------------------------------------"
          "-------------");
    CCTK_SchedulePrint (NULL);
    puts ("-------------------------------------------------------------------"
          "-------------");
  }

  return (retval);
}

 /*@@
   @routine    CCTKi_RecoverParameters
   @date       Sat May  3 18:49:25 PDT 2014
   @author     Roland Haas
   @desc
               Recover parameters from checkpoint
               This needs to be called right after
               CCTKi_ProcessParameterDatabase so that we always use the
               parameter values from the checkpoint and never the default
               values if the current parameter file does not set all
               parameters.
   @enddesc
   @calls      CCTKi_SetParameterSetMask
               CCTKi_BindingsParameterRecoveryInitialise
@@*/
static void CCTKi_RecoverParameters (tFleshConfig *config)
{
  extern void CCTKi_SetParameterSetMask (int mask);


  CCTKi_SetParameterSetMask (PARAMETER_RECOVERY_IN);

  config->recovered = CCTKi_BindingsParameterRecoveryInitialise ();
  if (config->recovered < 0)
  {
    CCTK_Warn (0, __LINE__, __FILE__, "Cactus", "Failed to recover parameters");
  }

  CCTKi_SetParameterSetMask (PARAMETER_RECOVERY_POST);
}
