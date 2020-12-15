 /*@@
   @file      IndexArrays.c
   @date      Mon Jun  3 13:15:30 2002
   @author    Ian Hawke
   @desc 
   Routines for dealing with the index arrays in 
   @seefile ExternalVariables.h
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_String.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_IndexArrays_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_SetupIndexArrays(CCTK_ARGUMENTS);

void MoL_FreeIndexArrays(CCTK_ARGUMENTS);

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
   @routine    MoL_SetupIndexArrays
   @date       Mon Jun  3 13:24:05 2002
   @author     Ian Hawke
   @desc 
   Allocates sufficient space for the index arrays.
   These arrays are defined in the external file
   @seefile ExternalVariables.h
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_SetupIndexArrays(CCTK_ARGUMENTS) 
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char *infoline;

  /* Initialize the time variables */

  *Original_Time = cctkGH->cctk_time;
  *Original_Delta_Time = cctkGH->cctk_delta_time;

  /* 
     Presumably a user will not try to register more variables than the total
     number of Cactus variables
   */

  MoLMaxNumRegisteredVariables = CCTK_NumVars();

  /* 
     We only want to set up the index arrays once.
     With mesh refinement this routine could be scheduled
     multiple times, leading to multiple copies of the
     index arrays used at different times!!!
  */

  if (EvolvedVariableIndex)
  {
    return;
  }


  if (MoLMaxNumRegisteredVariables)
  {
    EvolvedVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
					      sizeof(CCTK_INT));
    if (!EvolvedVariableIndex)
      {
	CCTK_ERROR("Failed to allocate the evolved variable index array");
      }

    RHSVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
					  sizeof(CCTK_INT));
    if (!RHSVariableIndex)
      {
	CCTK_ERROR("Failed to allocate the RHS variable index array");
      }
  }

  if (MoLMaxNumRegisteredVariables)
  {
    EvolvedVariableIndexSlow = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
					      sizeof(CCTK_INT));
    if (!EvolvedVariableIndexSlow)
      {
	CCTK_ERROR("Failed to allocate the slow evolved variable index array");
      }

    RHSVariableIndexSlow = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
					  sizeof(CCTK_INT));
    if (!RHSVariableIndexSlow)
      {
	CCTK_ERROR("Failed to allocate the slow RHS variable index array");
      }
  }

  if (MoLMaxNumRegisteredVariables)
  {
    ConstrainedVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
						  sizeof(CCTK_INT));
    if (!ConstrainedVariableIndex)
      {
	CCTK_ERROR("Failed to allocate the constrained variable index array");
      }
  }

  if (MoLMaxNumRegisteredVariables)
  {
    SandRVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
					    sizeof(CCTK_INT));
    if (!SandRVariableIndex)
    {
      CCTK_ERROR("Failed to allocate the save and restore "
		"variable index array");
    }
  }

  if (EvolvedArrayVariableIndex)
  {
    return;
  }

  if (MoLMaxNumRegisteredVariables)
  {
    EvolvedArrayVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
						   sizeof(CCTK_INT));
    if (!EvolvedArrayVariableIndex)
    {
      CCTK_ERROR("Failed to allocate the evolved array variable index array");
    }
  
    RHSArrayVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
					       sizeof(CCTK_INT));
    if (!RHSArrayVariableIndex)
    {
      CCTK_ERROR("Failed to allocate the RHS array variable index array");
    }
  }

  if (MoLMaxNumRegisteredVariables)
  {
    ConstrainedArrayVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
						       sizeof(CCTK_INT));
    if (!ConstrainedArrayVariableIndex)
    {
      CCTK_ERROR("Failed to allocate the constrained array "
                 "variable index array");
    }
  }

  if (MoLMaxNumRegisteredVariables)
  {
    SandRArrayVariableIndex = (CCTK_INT *)malloc(MoLMaxNumRegisteredVariables * 
						 sizeof(CCTK_INT));
    if (!SandRArrayVariableIndex)
    {
      CCTK_ERROR("Failed to allocate the save and restore "
                 "array variable index array");
    }
  }

  infoline = (char *)malloc((strlen(Generic_Method_Descriptor)+100)*
                            sizeof(char));
  if (!infoline) 
  {
    CCTK_ERROR("Failed to malloc 100 characters!");
  }
  if (CCTK_EQUALS(ODE_Method,"Generic")) 
  {
    if (CCTK_EQUALS(Generic_Type,"ICN"))
    {
      sprintf(infoline,"Generic Iterative Crank Nicholson with %i iterations", 
              (int)MoL_Intermediate_Steps);
    }
    else if (CCTK_EQUALS(Generic_Type,"RK"))
    {
      sprintf(infoline, "Generic Runge-Kutta %i",(int)MoL_Intermediate_Steps);
    }
    else if (CCTK_EQUALS(Generic_Type,"Table"))
    {
      sprintf(infoline, "Generic method, options:\n %s\n",
              Generic_Method_Descriptor);
    }
    else if (CCTK_EQUALS(Generic_Type,"Classic RK3")) 
    {
      sprintf(infoline, "Classic Runge-Kutta 3");
    }
    else
    {
      CCTK_ERROR("Generic_Type not recognized!");
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"Euler")) 
  {
    sprintf(infoline, "Euler");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2")) 
  {
    sprintf(infoline, "Runge-Kutta 2");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2-central")) 
  {
    sprintf(infoline, "Runge-Kutta 2 (central)");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK3")) 
  {
    sprintf(infoline, "Runge-Kutta 3");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK4")) 
  {
    sprintf(infoline, "Runge-Kutta 4");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK45")) 
  {
    sprintf(infoline, "Runge-Kutta 45 (Fehlberg)");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK45CK")) 
  {
    sprintf(infoline, "Runge-Kutta 45 (Cash-Karp)");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK65")) 
  {
    sprintf(infoline, "Runge-Kutta 65");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK87")) 
  {
    sprintf(infoline, "Runge-Kutta 87");
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN"))
  {
    sprintf(infoline, "Iterative Crank Nicholson with %i iterations", 
            (int)MoL_Intermediate_Steps);
  }
  else if (CCTK_EQUALS(ODE_Method,"ICN-avg"))
  {
    sprintf(infoline,
            "Averaging iterative Crank Nicholson with %i iterations", 
            (int)MoL_Intermediate_Steps);
  }
  else if (CCTK_EQUALS(ODE_Method,"AB"))
  {
    if (CCTK_EQUALS(AB_Type,"1"))
    {
      sprintf(infoline, "Adams-Bashforth of order 1");
    }
    else if (CCTK_EQUALS(AB_Type,"2"))
    {
      sprintf(infoline, "Adams-Bashforth of order 2");
    }
    else if (CCTK_EQUALS(AB_Type,"3"))
    {
      sprintf(infoline, "Adams-Bashforth of order 3");
    }
    else if (CCTK_EQUALS(AB_Type,"4"))
    {
      sprintf(infoline, "Adams-Bashforth of order 4");
    }
    else if (CCTK_EQUALS(AB_Type,"5"))
    {
      sprintf(infoline, "Adams-Bashforth of order 5");
    }
    else
    {
      CCTK_ERROR("AB_Type not recognized!");
    }
  }
  else if (CCTK_EQUALS(ODE_Method,"RK2-MR-2:1")) 
  {
    sprintf(infoline, "Multi-rate 2:1 Runge-Kutta 2");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK4-MR-2:1")) 
  {
    sprintf(infoline, "Multi-rate 2:1 Runge-Kutta 4");
  }
  else if (CCTK_EQUALS(ODE_Method,"RK4-RK2")) 
  {
    sprintf(infoline, "Multi-rate 2:1 Runge-Kutta 4 and Runge-Kutta 2");
  }
  else 
  {
    CCTK_ERROR("ODE_Method not recognized!");
  }
  
  CCTK_VInfo(CCTK_THORNSTRING, "Using %s as the time integrator.", infoline);
  
  free(infoline);
  infoline = NULL;

  // These scalars must be 1 oustide of the MoL loop.
  // They will only be zero for certain substeps when multirate methods are used.
  // Otherwise, they are guaranteed to always be ONE.
  *MoL_SlowPostStep = 1;
  *MoL_SlowStep = 1;

  return;

}

 /*@@
   @routine    MoL_FreeIndexArrays
   @date       Mon Jun  3 13:26:15 2002
   @author     Ian Hawke
   @desc 
   Frees the external index arrays. 
   These arrays are defined in the external file
   @seefile ExternalVariables.h
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FreeIndexArrays(CCTK_ARGUMENTS)
{
  
  if (EvolvedVariableIndex)
  {
    free(EvolvedVariableIndex);
    EvolvedVariableIndex = NULL;
  }

  if (RHSVariableIndex)
  {
    free(RHSVariableIndex);  
    RHSVariableIndex = NULL;
  }

  if (EvolvedVariableIndexSlow)
  {
    free(EvolvedVariableIndexSlow);
    EvolvedVariableIndexSlow = NULL;
  }

  if (RHSVariableIndexSlow)
  {
    free(RHSVariableIndexSlow);  
    RHSVariableIndexSlow = NULL;
  }

  if (ConstrainedVariableIndex)
  {
    free(ConstrainedVariableIndex);
    ConstrainedVariableIndex = NULL;
  }

  if (SandRVariableIndex)
  {
    free(SandRVariableIndex);
    SandRVariableIndex = NULL;
  }
  
  if (EvolvedArrayVariableIndex)
  {
    free(EvolvedArrayVariableIndex);
    EvolvedArrayVariableIndex = NULL;
  }

  if (RHSArrayVariableIndex)
  {
    free(RHSArrayVariableIndex);  
    RHSArrayVariableIndex = NULL;
  }

  if (ConstrainedArrayVariableIndex)
  {
    free(ConstrainedArrayVariableIndex);
    ConstrainedArrayVariableIndex = NULL;
  }

  if (SandRArrayVariableIndex)
  {
    free(SandRArrayVariableIndex);
    SandRArrayVariableIndex = NULL;
  }
  
  if (ArrayScratchSizes)
  {
    free(ArrayScratchSizes);
    ArrayScratchSizes = NULL;
  }
  
  if (ArrayScratchSpace)
  {
    free(ArrayScratchSpace);
    ArrayScratchSpace = NULL;
  }

  return;
  
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
