 /*@@
   @file      CactusDefaultEvolve.c
   @date      Thu Oct  8 17:28:46 1998
   @author    Tom Goodale
   @desc
              Default cactus evolution stuff.
   @enddesc
   @version   $Id$
 @@*/

/* #define DEBUG_CCTK 1 */

#include <stdio.h>
#include <stdlib.h>

#include "definethisthorn.h"

#include "cctk_Flesh.h"
#include "cctk_Parameters.h"
#include "cctk_Groups.h"
#include "cctk_WarnLevel.h"
#include "cctk_Termination.h"
#include "cctk_Main.h"
#include "cctk_Misc.h"
#include "cctk_IO.h"

#ifdef HAVE_TIME_GETTIMEOFDAY
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <unistd.h>
#endif

static const char *rcsid="$Header$";

CCTK_FILEVERSION(main_CactusDefaultEvolve_c);

/* Define some macros for convenience. */

#define ForallConvLevels(iteration, conv_level)                               \
        {                                                                     \
          int factor = 1;                                                     \
          for (conv_level = 0; conv_level < config->nGHs; conv_level++)       \
          {                                                                   \
            if (iteration % factor == 0)                                      \
            {

#define EndForallConvLevels                                                   \
            }                                                                 \
            factor *= 2;                                                      \
          }                                                                   \
        }

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int DoneMainLoop (const cGH *GH, CCTK_REAL cctk_time, int iteration);
static void StepGH (cGH *GH);

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

int CactusDefaultEvolve (tFleshConfig *config);

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    CactusDefaultEvolve
   @date       Thu Oct  8 17:30:15 1998
   @author     Tom Goodale
   @desc
               The default cactus evolution routine.
   @enddesc
   @calls      CCTK_Traverse
               CCTK_OutputGH
               DoneMainLoop
               StepGH
   @history
   @hdate      Fri May 12 2000
   @hauthor    Thomas Radke
   @hdesc      Moved evolution loop termination check into DoneMainLoop()
   @endhistory

   @var        config
   @vdesc      pointer to flesh configuration structure
   @vtype      tFleshConfig *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               always returns 0
   @endreturndesc
@@*/
int CactusDefaultEvolve (tFleshConfig *config)
{
  int var;
  unsigned int convergence_level, iteration;

  /*** Call OUTPUT for this GH (this routine    ***/
  /*** checks if output is necessary) and makes ***/
  /*** a Traverse with CCTK_ANALYSIS      ***/
  iteration = CCTK_MainLoopIndex ();
  ForallConvLevels (iteration, convergence_level)
  {
    CCTK_Traverse (config->GH[convergence_level], "CCTK_ANALYSIS");
    CCTK_OutputGH (config->GH[convergence_level]);
  }
  EndForallConvLevels;

  while (! DoneMainLoop (config->GH[0], config->GH[0]->cctk_time, iteration))
  {
    if (iteration == 0)
    {
      /* Can only use CactusDefaultEvolve with one timelevel */
      for (var = CCTK_NumVars () - 1; var >= 0; var--)
      {
        if (CCTK_DeclaredTimeLevelsVI (var) > 1)
        {
          CCTK_VWarn (0,__LINE__,__FILE__,"Cactus",
                      "Variable '%s' has multiple timelevels, default Cactus "
                      "evolve routine cannot rotate",
                      CCTK_VarName (var));
        }
      }
    }
    /* HERE ROTATE TIMELEVELS FOR ALL CONVERGENCE LEVELS */

    iteration = CCTK_SetMainLoopIndex (++iteration);

    /* Step each convergence level */
    ForallConvLevels (CCTK_MainLoopIndex (), convergence_level)
    {
      StepGH (config->GH[convergence_level]);
      /*
      CCTK_InfoOutput (config->GH[convergence_level], convergence_level);
      */
    }
    EndForallConvLevels;

    /* Dump out checkpoint data on all levels */
    ForallConvLevels (iteration, convergence_level)
    {
      CCTK_Traverse (config->GH[convergence_level], "CCTK_CHECKPOINT");
    }
    EndForallConvLevels;

    /*** Call OUTPUT for this GH (this routine    ***/
    /*** checks if output is necessary) and makes ***/
    /*** an Traverse with CCTK_ANALYSIS      ***/
    ForallConvLevels (iteration, convergence_level)
    {
        CCTK_Traverse (config->GH[convergence_level], "CCTK_ANALYSIS");
        CCTK_OutputGH (config->GH[convergence_level]);
    }
    EndForallConvLevels;

  } /*** END OF MAIN ITERATION LOOP ***/

  return (0);
}


 /*@@
   @routine    DoneMainLoop
   @date       Fri May 12 2000
   @author     Thomas Radke
   @desc
               Check the termination conditions for the evolution loop
   @enddesc
   @calls      CCTK_TerminationReached
   @history
   @hdate      Thu 6 Nov 2002
   @hauthor    Thomas Radke
   @hdesc      Added max_runtime condition test
   @endhistory

   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        simulation_time
   @vdesc      current physical simulation time
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        iteration
   @vdesc      current iteration
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               true (1) or false (0) for done/not done with main loop
   @endreturndesc
@@*/
static int DoneMainLoop (const cGH *GH, CCTK_REAL simulation_time,int iteration)
{
  int max_iteration_reached, max_simulation_time_reached, max_runtime_reached;
  int retval;
#ifdef HAVE_TIME_GETTIMEOFDAY
  struct timeval runtime;
  static struct timeval starttime = {0, 0};
#endif
  DECLARE_CCTK_PARAMETERS


#ifdef HAVE_TIME_GETTIMEOFDAY
  /* on the first time through, get the start time */
  if (starttime.tv_sec == 0 && starttime.tv_usec == 0)
  {
    gettimeofday (&starttime, NULL);
  }
#endif

  retval = terminate_next || CCTK_TerminationReached (GH);
  if (! retval && ! CCTK_Equals (terminate, "never"))
  {
    max_iteration_reached = iteration >= cctk_itlast;

    if (cctk_initial_time < cctk_final_time)
    {
      max_simulation_time_reached = simulation_time >= cctk_final_time;
    }
    else
    {
      max_simulation_time_reached = simulation_time <= cctk_final_time;
    }

    max_runtime_reached = 0;
#ifdef HAVE_TIME_GETTIMEOFDAY
    if (max_runtime > 0)
    {
      /* get the elapsed runtime in minutes and compare with max_runtime */
      gettimeofday (&runtime, NULL);
      runtime.tv_sec -= starttime.tv_sec;
      max_runtime_reached = ((CCTK_REAL) runtime.tv_sec / 60.0) >= max_runtime;
    }
#endif

    if (CCTK_Equals (terminate, "iteration"))
    {
      retval = max_iteration_reached;
    }
    else if (CCTK_Equals (terminate, "time"))
    {
      retval = max_simulation_time_reached;
    }
    else if (CCTK_Equals (terminate, "runtime"))
    {
      retval = max_runtime_reached;
    }
    else if (CCTK_Equals (terminate, "any"))
    {
      retval = max_iteration_reached || max_simulation_time_reached ||
               max_runtime_reached;
    }
    else if (CCTK_Equals (terminate, "all"))
    {
      retval = max_iteration_reached && max_simulation_time_reached &&
               max_runtime_reached;
    }
    /* the following two conditions are deprecated in BETA14 */
    else if (CCTK_Equals (terminate, "either"))
    {
      retval = max_iteration_reached || max_simulation_time_reached;
    }
    else /* if (CCTK_Equals (terminate, "both")) */
    {
      retval = max_iteration_reached && max_simulation_time_reached;
    }
  }

  return (retval);
}


 /*@@
   @routine    StepGH
   @date       Fri Aug 14 12:39:49 1998
   @author     Gerd Lanfermann
   @desc
               The full set of routines used to execute all schedule points
               in the main iteration loop. Makes calls to the individual
               routines for each schedule point.
   @enddesc
   @calls      CCTK_Traverse

   @var        GH
   @vdesc      pointer to CCTK grid hierachy
   @vtype      cGH *
   @vio        in
   @endvar
 @@*/
static void StepGH (cGH *GH)
{
  GH->cctk_time = GH->cctk_time + GH->cctk_delta_time;
  GH->cctk_iteration++;

  CCTK_Traverse (GH, "CCTK_PRESTEP");
  CCTK_Traverse (GH, "CCTK_EVOL");
  CCTK_Traverse (GH, "CCTK_POSTSTEP");
}

