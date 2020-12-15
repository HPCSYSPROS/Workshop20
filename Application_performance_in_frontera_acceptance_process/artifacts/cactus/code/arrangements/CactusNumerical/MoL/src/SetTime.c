 /*@@
   @file      SetTime.c
   @date      Mon May 20 09:45:45 2002
   @author    Ian Hawke
   @desc 
   Sets the time and dt depending on the ODE method and 
   position in the loop.
   @enddesc 
   @version   $Header$
 @@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_SetTime_c);

/* #define MOLDEBUG 1 */

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_SetTime(CCTK_ARGUMENTS);

void MoL_ResetTime(CCTK_ARGUMENTS);

void MoL_ResetDeltaTime(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/* RK45 Fehlberg coefficients */
static const CCTK_REAL alpha_array_F[6] = {
  0.0,
  1.0/4.0,
  3.0/8.0,
  12.0/13.0,
  1.0,
  1.0/2.0,
};

/* RK45 Cash-Karp coefficients */
static const CCTK_REAL alpha_array_CK[6] = {
  0.0,
  1.0/5.0,
  3.0/10.0,
  3.0/5.0,
  1.0,
  7.0/8.0,
};

/* RK65 coefficients */
static const CCTK_REAL alpha_array65[8] = {
  0.0,
  1.0/10.0,
  2.0/9.0,
  3.0/7.0,
  3.0/5.0,
  4.0/5.0,
  1.0,
  1.0
};

/* RK87 coefficients */
 static const CCTK_REAL alpha_array87[13] = {
  0.0,
  1.0/18.0,
  1.0/12.0,
  1.0/8.0,
  5.0/16.0,
  3.0/8.0,
  59.0/400.0,
  93.0/200.0,
  5490023248.0/9719169821.0,
  13.0/20.0,
  1201146811.0/1299019798.0,
  1.0,
  1.0
};

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    MoL_SetTime
   @date       Mon May 20 09:48:55 2002
   @author     Ian Hawke
   @desc 
   Sets the time and timestep before the MoL evolution loop starts.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_SetTime(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL beta;

  if (adaptive_stepsize && ! CCTK_EQUALS(verbose, "none"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Integrating from %g to %g with step size %g",
                (double)(cctkGH->cctk_time - cctkGH->cctk_delta_time),
                (double)cctkGH->cctk_time,
                (double)cctkGH->cctk_delta_time);
  }

  *Original_Time = cctkGH->cctk_time;
  *Original_Delta_Time = cctkGH->cctk_delta_time;
  cctkGH->cctk_time -= cctkGH->cctk_delta_time / cctkGH->cctk_timefac;

  if (CCTK_EQUALS(ODE_Method,"ICN"))
  {
    cctkGH->cctk_delta_time = 0.5*(*Original_Delta_Time);
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN-avg"))
  {
    cctkGH->cctk_delta_time = *Original_Delta_Time;
  }
  else if (CCTK_EQUALS(ODE_Method,"Generic"))
  {
    beta = RKBetaCoefficients[0];

    cctkGH->cctk_delta_time = beta*(*Original_Delta_Time);
  }
}

 /*@@
   @routine    MoL_ResetTime
   @date       Mon May 20 09:49:41 2002
   @author     Ian Hawke
   @desc 
   Sets the time during the MoL evolution loop.
   At the last time all methods should end up with the original
   values for time and timestep.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ResetTime(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL *previous_times;
  CCTK_INT alphaindex, i, j;
    
  previous_times = (CCTK_REAL*)malloc(MoL_Intermediate_Steps*sizeof(CCTK_REAL));
  if (!previous_times)
  {
    CCTK_ERROR("Failed to allocate memory for very small array!");
  }

  if (*MoL_Intermediate_Step == 0) 
  {
    cctkGH->cctk_time = (*Original_Time);
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN"))
  {
    cctkGH->cctk_time = (*Original_Time) - 
      0.5*(*Original_Delta_Time)/cctkGH->cctk_timefac;
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN-avg"))
  {
    cctkGH->cctk_time = (*Original_Time);
  }
  else if (CCTK_EQUALS(ODE_Method,"Generic"))
  {
    previous_times[0] = (*Original_Time) -
      (*Original_Delta_Time)/cctkGH->cctk_timefac;
    for (i = MoL_Intermediate_Steps - 1; i > *MoL_Intermediate_Step - 1; i--)
    {
      previous_times[MoL_Intermediate_Steps - i] = 
        RKBetaCoefficients[MoL_Intermediate_Steps - i - 1] * 
        (*Original_Delta_Time)/cctkGH->cctk_timefac;
      for (j = MoL_Intermediate_Steps; j > i; j--)
      {
        alphaindex = (MoL_Intermediate_Steps - i - 1) * 
          MoL_Intermediate_Steps + MoL_Intermediate_Steps - j;
        previous_times[MoL_Intermediate_Steps - i] += 
          RKAlphaCoefficients[alphaindex] *
          previous_times[MoL_Intermediate_Steps - j];
#ifdef MOLDEBUG
        printf("i %d j %d is %d index %d t %g dt %g alpha %g beta %g\n",
               i, j, MoL_Intermediate_Steps, alphaindex,
               previous_times[MoL_Intermediate_Steps - i], 
               (*Original_Delta_Time)/cctkGH->cctk_timefac,
               RKAlphaCoefficients[alphaindex], 
               RKBetaCoefficients[MoL_Intermediate_Steps - i - 1]);
#endif
      }
    }
#ifdef MOLDEBUG
    printf("MoL says the previous times are ");
    for (i = 0; i < MoL_Intermediate_Steps - *MoL_Intermediate_Step
           + 1; i++)
    {
      printf("%g ", previous_times[i]);
    }
    printf("\n");
#endif
    cctkGH->cctk_time = previous_times[MoL_Intermediate_Steps - 
                                      *MoL_Intermediate_Step];
  }
  else if (CCTK_EQUALS(ODE_Method,"Euler"))
  {
    if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_time = (*Original_Time);
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2"))
  {
    if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_time = (*Original_Time);
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2-central"))
  {
    if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_time = (*Original_Time) - 
        0.5*(*Original_Delta_Time)/cctkGH->cctk_timefac;
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK3"))
  {
    if (*MoL_Intermediate_Step == 2)
    {
      cctkGH->cctk_time = (*Original_Time);
    }
    else if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_time = (*Original_Time) - 
        0.5*(*Original_Delta_Time)/cctkGH->cctk_timefac;
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK4"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);

    CCTK_REAL dt = (*Original_Delta_Time)/cctkGH->cctk_timefac;
    switch (substep)
    {
      case 1:
      case 2:
        dt *= 0.5;
        break;
      default:
        dt = 0;
    }
    cctkGH->cctk_time = (*Original_Time) - dt;
  }
  else if (CCTK_EQUALS(ODE_Method,"RK45") || CCTK_EQUALS(ODE_Method,"RK45CK"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);
    const CCTK_REAL * alpha_array;
    if (CCTK_EQUALS(ODE_Method, "RK45"))
    {
      alpha_array = alpha_array_F;
    }
    else if (CCTK_EQUALS(ODE_Method, "RK45CK"))
    {
      alpha_array = alpha_array_CK;
    }
    else
    {
      CCTK_ERROR ("internal error");
      /* Avoid compiler warning */
      alpha_array = NULL;
    }
    cctkGH->cctk_time
      = ((* Original_Time)
         + ((alpha_array[substep] - 1)
            * (* Original_Delta_Time) / cctkGH->cctk_timefac));
  }
  else if (CCTK_EQUALS(ODE_Method,"RK65"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);
    cctkGH->cctk_time
      = ((* Original_Time)
         + ((alpha_array65[substep] - 1)
            * (* Original_Delta_Time) / cctkGH->cctk_timefac));
  }
  else if (CCTK_EQUALS(ODE_Method,"RK87"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);
    cctkGH->cctk_time
      = ((* Original_Time)
         + ((alpha_array87[substep] - 1)
            * (* Original_Delta_Time) / cctkGH->cctk_timefac));
  }
  else if (CCTK_EQUALS(ODE_Method,"AB"))
  {
    CCTK_ERROR ("internal error");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2-MR-2:1"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);

    CCTK_REAL dt = (*Original_Delta_Time)/cctkGH->cctk_timefac;
    switch (substep)
    {
      case 1:
      case 2:
        dt *= 0.5;
        break;
      default:
        dt = 0;
    }
    cctkGH->cctk_time = (*Original_Time) - dt;
  }
  else if (CCTK_EQUALS(ODE_Method,"RK4-MR-2:1"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);

    CCTK_REAL dt = (*Original_Delta_Time)/cctkGH->cctk_timefac;
    switch (substep)
    {
      case 1:
      case 2:
        dt *= 3.0/4.0;
        break;
      case 3:
      case 4:
      case 5:
        dt *= 1.0/2.0;
        break;
      case 6:
      case 7:
        dt *= 1.0/4.0;
        break;
      default:
        dt = 0;
    }
    cctkGH->cctk_time = (*Original_Time) - dt;
  }
  else if (CCTK_EQUALS(ODE_Method,"RK4-RK2"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);

    CCTK_REAL dt = (*Original_Delta_Time)/cctkGH->cctk_timefac;
    switch (substep)
    {
      case 1:
      case 2:
        dt *= 0.5;
        break;
      default:
        dt = 0;
    }
    cctkGH->cctk_time = (*Original_Time) - dt;
  }
#ifdef MOLDEBUG
  printf("MoL has once more reset t (%d): %f.\n", 
         *MoL_Intermediate_Step, cctkGH->cctk_time);
  fflush(stdout);
#endif

  free(previous_times);
  previous_times = NULL;

  if (adaptive_stepsize && CCTK_EQUALS(verbose, "extreme"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Evaluating RHS at %g",
                (double)cctkGH->cctk_time);
  }
}

 /*@@
   @routine    MoL_ResetDeltaTime
   @date       Mon May 20 09:49:41 2002
   @author     Ian Hawke
   @desc 
   Sets the timestep during the MoL evolution loop.
   At the last time all methods should end up with the original
   values for time and timestep.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ResetDeltaTime(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (*MoL_Intermediate_Step == 0) 
  {
    cctkGH->cctk_delta_time = (*Original_Delta_Time);
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN"))
  {
    if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_delta_time = (*Original_Delta_Time);
    }
    else
    {
      cctkGH->cctk_delta_time = 0.5*(*Original_Delta_Time);
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN-avg"))
  {
    cctkGH->cctk_delta_time = (*Original_Delta_Time);
  }
  else if (CCTK_EQUALS(ODE_Method,"Generic"))
  {
    cctkGH->cctk_delta_time = RKBetaCoefficients[MoL_Intermediate_Steps - 
                                                (*MoL_Intermediate_Step)] *
      (*Original_Delta_Time);
  }
  else if (CCTK_EQUALS(ODE_Method,"Euler"))
  {
    if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_delta_time = (*Original_Delta_Time);
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2"))
  {
    if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_delta_time = 0.5*(*Original_Delta_Time);
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK3"))
  {
    if (*MoL_Intermediate_Step == 2)
    {
      cctkGH->cctk_delta_time = 0.25*(*Original_Delta_Time);
    }
    else if (*MoL_Intermediate_Step == 1)
    {
      cctkGH->cctk_delta_time = 2.0/3.0*(*Original_Delta_Time);
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK45") || CCTK_EQUALS(ODE_Method,"RK45CK"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);
    const CCTK_REAL * alpha_array;
    if (CCTK_EQUALS(ODE_Method, "RK45"))
    {
      alpha_array = alpha_array_F;
    }
    else if (CCTK_EQUALS(ODE_Method, "RK45CK"))
    {
      alpha_array = alpha_array_CK;
    }
    else
    {
      CCTK_ERROR ("internal error"); 
      /* Avoid compiler warning */
      alpha_array = NULL;
    }
    cctkGH->cctk_delta_time
      = ((alpha_array[substep + 1] - alpha_array[substep])
         * (* Original_Delta_Time));
  }
  else if (CCTK_EQUALS(ODE_Method,"RK65"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);
    cctkGH->cctk_delta_time
      = ((alpha_array65[substep + 1] - alpha_array65[substep])
         * (* Original_Delta_Time));
  }
  else if (CCTK_EQUALS(ODE_Method,"RK87"))
  {
    const int substep = MoL_Intermediate_Steps - (* MoL_Intermediate_Step);
    cctkGH->cctk_delta_time
      = ((alpha_array87[substep + 1] - alpha_array87[substep])
         * (* Original_Delta_Time));
  }
#ifdef MOLDEBUG
  printf("MoL has once more reset dt (%d): %f.\n", 
         *MoL_Intermediate_Step, 
         cctkGH->cctk_delta_time/cctkGH->cctk_timefac);
  fflush(stdout);
#endif

  if (adaptive_stepsize && CCTK_EQUALS(verbose, "extreme"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Evaluating RHS with a time step of %g",
                (double)cctkGH->cctk_delta_time);
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
