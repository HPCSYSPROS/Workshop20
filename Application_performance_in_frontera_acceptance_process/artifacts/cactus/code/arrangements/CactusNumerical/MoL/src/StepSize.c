 /*@@
   @file      StepSize.c
   @date      Tue Sep 07 2004
   @author    Erik Schnetter
   @desc 
   Control the time step size.
   @enddesc 
   @version   $Header$
 @@*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_StepSize_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_StartLoop(CCTK_ARGUMENTS);

void MoL_InitAdaptiveError(CCTK_ARGUMENTS);
void MoL_FindAdaptiveError(CCTK_ARGUMENTS);
void MoL_ReduceAdaptiveError(CCTK_ARGUMENTS);

void MoL_SetEstimatedDt(CCTK_ARGUMENTS);

void MoL_FinishLoop(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    MoL_StartLoop
   @date       Tue Sep 07 2004
   @author     Erik Schnetter
   @desc 
   Start the step size control loop, so that at least one iteration is done.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void
MoL_StartLoop(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  *MoL_Stepsize_Bad = 1;

  if (adaptive_stepsize)
  {
    *EstimatedDt = cctkGH->cctk_delta_time;
  }
  
}

 /*@@
   @routine    MoL_InitAdaptiveError
   @date       Tue Sep 07 2004
   @author     Erik Schnetter
   @desc 
   Initialize error counters for adaptive stepsize control
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

static inline CCTK_REAL
square (CCTK_REAL const x)
{
  return x * x;
}

void MoL_InitAdaptiveError(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Initialise global error */
  *Error = 0;
  *Count = 0;
}

 /*@@
   @routine    MoL_FindAdaptiveError
   @date       Thu Jan 27 10:22:26 2005
   @author     Erik Schnetter
   @desc 
   Compute the error local to this component/patch/...
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FindAdaptiveError(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* TODO: exclude symmetry boundaries */
  assert (cctk_dim <= 3);
  int imin[3], imax[3];
  for (int d = 0; d < cctk_dim; d++)
  {
    imin[d] = cctk_bbox[2*d] ? 0 : cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - (cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d]);
  }

  /* Calculate absolute error */
  CCTK_REAL local_error = 0.0;
  for (int var = 0; var < MoLNumEvolvedVariables; var++)
  {

    CCTK_REAL const * restrict const
      UpdateVar = CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);
    CCTK_REAL const * restrict const
      RHSVar = CCTK_VarDataPtrI(cctkGH, 0, RHSVariableIndex[var]);
    CCTK_REAL const * restrict const
      ErrorVar
      = CCTK_VarDataPtrI(cctkGH, var,
                         CCTK_FirstVarIndex("MOL::ERRORESTIMATE"));

    CCTK_REAL const rhs_relative_error
      = maximum_relative_error * RHS_error_weight * (*Original_Delta_Time);

    assert (cctk_dim == 3);
#pragma omp parallel for reduction(+: local_error)
    for (int k = imin[2]; k < imax[2]; k++)
    {
      for (int j = imin[1]; j < imax[1]; j++)
      {
        for (int i = imin[0]; i < imax[0]; i++)
        {
          int const index = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL const scale
            = (square(maximum_absolute_error)
               + square(maximum_relative_error * UpdateVar[index])
               + square(rhs_relative_error * RHSVar[index]));
          local_error += square(ErrorVar[index]) / scale;
        }
      }
    }

  } /* for var */

  *Error += local_error;
  *Count
    += (MoLNumEvolvedVariables
        * (imax[0] - imin[0]) * (imax[1] - imin[1]) * (imax[2] - imin[2]));
}

 /*@@
   @routine    MoL_ReduceAdaptiveError
   @date       Thu Jan 27 10:23:14 2005
   @author     Erik Schnetter
   @desc 
   Find the global error estimate. 
   Change the timestep based on the error estimate.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ReduceAdaptiveError(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int redop;

  CCTK_REAL red_local[2], red_global[2];
  int p1, p2;
  int ierr;

  /* Get global result over all processors */
  redop = CCTK_ReductionHandle ("sum");
  assert (redop >= 0);

  red_local[0] = *Error;
  red_local[1] = *Count;
  ierr = CCTK_ReduceLocArrayToArray1D
    (cctkGH, -1, redop, red_local, red_global, 2, CCTK_VARIABLE_REAL);
  assert (ierr == 0);
  *Error = red_global[0];
  *Count = red_global[1];

  /* Calculate L2-norm */
  *Error = sqrt(*Error / *Count);
  if (! CCTK_EQUALS(verbose, "none"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Integration accuracy quotient is %g", (double)*Error);
  }
  if (! isfinite(*Error))
  {
    CCTK_VWarn (CCTK_WARN_ALERT,__LINE__,__FILE__,CCTK_THORNSTRING,
                "Integration accuracy quotient is %g, which is not a finite number -- reducing the step size", (double)*Error);
  }

  if ( CCTK_EQUALS(ODE_Method,"RK45") || CCTK_EQUALS(ODE_Method,"RK45CK") )
  {
    p1 = 5;
    p2 = 4;
  }
  else if ( CCTK_EQUALS(ODE_Method,"RK65") )
  {
    p1 = 6;
    p2 = 5;
  }
  else if ( CCTK_EQUALS(ODE_Method,"RK87") )
  {
    p1 = 8;
    p2 = 7;
  }
  else
  {
    CCTK_ERROR ("unsupported ODE_Method in stepsize control");
    /* Avoid compiler warnings */
    p1 = 0;
    p2 = 0;
  }

  /* Decide whether to accept this step */
  *MoL_Stepsize_Bad = ! isfinite(*Error) || *Error > 1;

  if (*MoL_Stepsize_Bad)
  {
    /* The error is too large; reject the time step and reduce the
       step size */
    cctkGH->cctk_time -= cctkGH->cctk_delta_time;
    if (! CCTK_EQUALS(verbose, "none"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "*** REJECTING TIME STEP ***");
    }

    if (isfinite(*Error))
    {
      cctkGH->cctk_delta_time
        = safety_factor * (*Original_Delta_Time) / pow(*Error, 1.0/p1);
    }
    else
    {
      cctkGH->cctk_delta_time = (*Original_Delta_Time) / maximum_decrease;
    }
    /* if (! CCTK_EQUALS(verbose, "none")) */
    /* { */
    /*   CCTK_VInfo (CCTK_THORNSTRING, "Setting time step to %g", (double)cctkGH->cctk_delta_time); */
    /* } */

    if (cctkGH->cctk_delta_time < (*Original_Delta_Time) / maximum_decrease)
    {
      /* No more than a factor of 10 decrease */
      cctkGH->cctk_delta_time = (*Original_Delta_Time) / maximum_decrease;
      if (! CCTK_EQUALS(verbose, "none"))
      {
        CCTK_VInfo (CCTK_THORNSTRING, "   Time step reduction too large; clamping time step to %g", (double)cctkGH->cctk_delta_time);
      }
    }
    if (cctkGH->cctk_delta_time == (CCTK_REAL)0.0)
      /* yes, we want to compare to zero exactly, to catch underflows */
      {
        CCTK_ERROR ("New step size would be zero -- aborting");
      }

    cctkGH->cctk_time += cctkGH->cctk_delta_time;
  }
  else
  {
    /* The error is acceptable; estimate the next step size */
    *EstimatedDt = 
      safety_factor * (*Original_Delta_Time) / pow(*Error, 1.0/p2);

    if (*EstimatedDt > (*Original_Delta_Time) * maximum_increase)
    {
      /* No more than a factor of 5 increase */
      *EstimatedDt = (*Original_Delta_Time) * maximum_increase;
      if (! CCTK_EQUALS(verbose, "none"))
      {
        CCTK_VInfo (CCTK_THORNSTRING, "   Time step increase too large; clamping time step to %g", (double)(*EstimatedDt));
      }
    }
  }
}

 /*@@
   @routine    MoL_SetEstimatedDt
   @date       Thu Jan 27 14:05:08 2005
   @author     Ian Hawke
   @desc 
   Actually set the timestep in PostStep. 
   This avoids problems when the timestep is changed in the middle of the
   evolution loop.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_SetEstimatedDt(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  cctkGH->cctk_delta_time = *EstimatedDt;

  if (! CCTK_EQUALS(verbose, "none"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Setting time step to %g", 
                (double)cctkGH->cctk_delta_time);
  }

}

 /*@@
   @routine    MoL_FinishLoop
   @date       Thu Jan 27 10:24:19 2005
   @author     Erik Schnetter
   @desc 
   Loop control if adaptive timestepping is not used.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FinishLoop(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Keep time step size unchanged */
  *MoL_Stepsize_Bad = 0;
  
  // Set these flags to ONE outside of MoL-loop!
  *MoL_SlowPostStep = 1;
  *MoL_SlowStep = 1;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
