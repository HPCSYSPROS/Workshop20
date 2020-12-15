 /*@@
   @file      CactusDefaultInitialise.c
   @date      Tue Sep 29 12:45:04 1998
   @author    Tom Goodale
   @desc
              Default Cactus initialisation routine.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>

#include "cctk_Misc.h"
#include "cctk_Flesh.h"
#include "cctk_Parameter.h"

#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "CactusMainDefaults.h"
#include "CactusCommFunctions.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_CactusDefaultInitialise_c);

/*# define DEBUG_CCTK 1 */

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
static void CactusInitialiseGH (const tFleshConfig *config, cGH *GH);


 /*@@
   @routine    CactusDefaultInitialise
   @date       Tue Sep 29 12:45:04 1998
   @author     Tom Goodale
   @desc
               Default initialisation routine.
   @enddesc
   @calls      CCTK_SetupGH
               CCTKi_AddGH
               CactusInitialiseGH

   @var        config
   @vdesc      flesh configuration structure
   @vtype      tFleshConfig *
   @vio        inout
   @endvar

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
int CactusDefaultInitialise (tFleshConfig *config)
{
  cGH *GH;
  int convergence_level;

#if 0
  CactusResetTimer (config->timer[INITIALISATION]);
  CactusResetTimer (config->timer[EVOLUTION]);
  CactusResetTimer (config->timer[ELLIPTIC]);

  CactusStartTimer (config->timer[INITIALISATION]);
#endif

  convergence_level = 0;
  while ((GH = CCTK_SetupGH (config, convergence_level)))
  {
    CCTKi_AddGH (config, convergence_level, GH);

    CactusInitialiseGH (config, GH);

    convergence_level++;
  };

#if 0
  CactusStopTimer (config->timer[INITIALISATION]);
#endif

  return (0);
}


 /*@@
   @routine    CactusInitialiseGH
   @date       Mon Feb  1 12:13:09 1999
   @author     Tom Goodale
   @desc
               Responsible for initialising a GH.
   @enddesc
   @calls      CCTKi_ScheduleGHInit
               CCTKi_InitGHExtensions
               CCTKi_FinaliseParamWarn
               CCTK_Traverse

   @var        config
   @vdesc      flesh configuration structure
   @vtype      tFleshConfig *
   @vio        inout
   @endvar
   @var        GH
   @vdesc      the GH to initialize
   @vtype      cGH *
   @vio        inout
   @endvar
@@*/
static void CactusInitialiseGH (const tFleshConfig *config, cGH *GH)
{
  const char *recovery_mode;


  recovery_mode = *(const char *const *)
                  CCTK_ParameterGet ("recovery_mode", "Cactus", NULL);

  /* Initialise time */
  GH->cctk_time = *(const CCTK_REAL *)
                  CCTK_ParameterGet ("cctk_initial_time", "Cactus", NULL);

  /* Initialise iteration number */
  GH->cctk_iteration = 0;

#ifdef DEBUG_CCTK
  CCTK_PRINTSEPARATOR
  printf ("In Cactus_Initialise\n--------------------\n");
  printf ("  Initializing GH->cctk_time = %f\n", GH->cctk_time);
  printf ("  Initializing GH->cctk_iteration = %u\n", GH->cctk_iteration);
  CCTK_PRINTSEPARATOR
#endif

  /* Do the schedule initialisation on this GH */
  CCTKi_ScheduleGHInit (GH);

  /* Initialise all the extensions. */
  CCTKi_InitGHExtensions (GH);

  CCTK_Traverse (GH, "CCTK_WRAGH");

  /* FIXME : PARAM_CHECK SHOULD BE BEFORE HERE */
  CCTK_Traverse (GH, "CCTK_PARAMCHECK");
  CCTKi_FinaliseParamWarn ();

  CCTK_Traverse (GH, "CCTK_BASEGRID");

  if (! (config->recovered && CCTK_Equals (recovery_mode, "strict")))
  {
    /* Traverse routines setting up initial data */
    CCTK_Traverse (GH, "CCTK_INITIAL");

    /* Traverse poststep initial routines which should only be done once */
    CCTK_Traverse (GH, "CCTK_POSTINITIAL");
    CCTK_Traverse (GH, "CCTK_POSTPOSTINITIAL");
    CCTK_Traverse (GH, "CCTK_POSTSTEP");
  }

  /* Traverse recovery and post-recovery routines */
  if (config->recovered)
  {
    CCTK_Traverse (GH, "CCTK_RECOVER_VARIABLES");
    CCTK_Traverse (GH, "CCTK_POST_RECOVER_VARIABLES");
  }

  /* Traverse ID checkpoint routines */
  CCTK_Traverse (GH, "CCTK_CPINITIAL");
}
