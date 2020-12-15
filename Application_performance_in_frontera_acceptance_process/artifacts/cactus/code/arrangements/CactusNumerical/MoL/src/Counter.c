 /*@@
   @file      Counter.c
   @date      Mon May 20 09:52:33 2002
   @author    Ian Hawke
   @desc 
   Routines setting and altering the loop counter.
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/* #include "carpet.h" */

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_Counter_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_SetCounter(CCTK_ARGUMENTS);

void MoL_DecrementCounter(CCTK_ARGUMENTS);

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
   @routine    MoL_SetCounter
   @date       Mon May 20 09:54:39 2002
   @author     Ian Hawke
   @desc 
   Initially set the counter to the number of intermediate steps.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_SetCounter(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *MoL_Intermediate_Step = MoL_Intermediate_Steps;

  // Always indicate execution of slow steps.
  // If multirate methods are used, they are set to zero
  // for certain substeps.
  // In any case, outside the MoL loop, these scalars are guaranteed to be 1
  *MoL_SlowStep = 1;
  *MoL_SlowPostStep = 1;
  
/* #ifdef HAVE_CARPET */
  if ((*MoL_Intermediate_Step))
  {
    /* Disable prolongating during the iterations */
/*     CarpetEnableProlongating (0); */
    if (disable_prolongation)
    {
      if (CCTK_IsFunctionAliased("EnableProlongating"))
      {
        EnableProlongating(0);
      }
      else
      {
        CCTK_WARN(CCTK_WARN_DEBUG, "Cannot disable prolongation as function"
                  " \"EnableProlongating\" is not provided by any thorn!");
      }
    }
  }
/* #endif */
}

 /*@@
   @routine    MoL_DecrementCounter
   @date       Mon May 20 09:55:13 2002
   @author     Ian Hawke
   @desc 
   During the loop decrement the counter
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_DecrementCounter(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  (*MoL_Intermediate_Step) --;


  if (CCTK_EQUALS(ODE_Method, "RK2-MR-2:1"))
  {
     switch (MoL_Intermediate_Steps - *MoL_Intermediate_Step)
     {
        case 2:
        case 4:
           if (*MoL_SlowStep)
              *MoL_SlowPostStep = 1;
           else
              *MoL_SlowPostStep = 0;
           // don't do slow step
           *MoL_SlowStep = 0;
           break;
        default:
           if (!(*MoL_SlowStep))
              *MoL_SlowPostStep = 0;
           else
              *MoL_SlowPostStep = 1;
           // do a slow step!
           *MoL_SlowStep = 1;
     }
  }

  if (CCTK_EQUALS(ODE_Method, "RK4-MR-2:1"))
  {
     switch (MoL_Intermediate_Steps - *MoL_Intermediate_Step)
     {
        case 2:
        case 4:
        case 7:
        case 9:
           // don't do slow step
           *MoL_SlowStep = 0;
           break;
        default:
           // do a slow step!
           *MoL_SlowStep = 1;
     }
  }

  if (CCTK_EQUALS(ODE_Method, "RK4-RK2"))
  {
     switch (MoL_Intermediate_Steps - *MoL_Intermediate_Step)
     {
        case 0:
           // do a slow step!
           *MoL_SlowStep = 1;
           // don't sync!
           *MoL_SlowPostStep = 0;
           break;
        case 1:
        case 2:
           // don't do slow step
           *MoL_SlowStep = 0;
           *MoL_SlowPostStep = 0;
           break;
        case 3:
           // do a slow step!
           *MoL_SlowStep = 1;
           // sync!
           *MoL_SlowPostStep = 1;
     }
  }

/* #ifdef HAVE_CARPET */
  if (! (*MoL_Intermediate_Step))
  {
    /* Re-enable prolongating before the final PostStep */
/*     CarpetEnableProlongating (1); */
    if (disable_prolongation)
    {
      if (CCTK_IsFunctionAliased("EnableProlongating"))
      {
        EnableProlongating(1);
      }
      else
      {
        CCTK_WARN(CCTK_WARN_DEBUG, "Cannot enable prolongation as function"
                  " \"EnableProlongating\" is not provided by any thorn!");
      }
    }
  }
/* #endif */
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
