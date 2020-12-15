 /*@@
   @file      SandR.c
   @date      Sun May 26 03:35:58 2002
   @author    Ian Hawke
   @desc 
   Restores the Save and Restore variables to their original positions.
   @enddesc 
   @version   $Header$
 @@*/

#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"
#include "Operators.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_SandR_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RestoreSandR(CCTK_ARGUMENTS);

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
   @routine    MoL_RestoreSandR
   @date       Sun May 26 03:39:02 2002
   @author     Ian Hawke
   @desc 
   Save and Restore variables are those that the physics thorn may 
   need to know to calculate the RHS, but which may be evolved by
   something other than MoL. In order to get the timelevels correct,
   the previous timelevel is copied to the current before the MoL step.
   As we do not know whether the variable will be evolved before or
   after MoL we must save the data that was in the current timelevel, 
   and then restore it at the end of the MoL timestep. This routine
   restores the variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RestoreSandR(CCTK_ARGUMENTS) 
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int rl = 0;
  if (CCTK_IsFunctionAliased("GetRefinementLevel")) {
    rl = GetRefinementLevel(cctkGH);
  }
  const int tl = 0;

  const int scratchvarindex = CCTK_FirstVarIndex("MOL::SANDRSCRATCHSPACE");
  if (scratchvarindex < 0)
  {
    CCTK_ERROR("Internal error");
  }

  for (CCTK_INT var = 0; var < MoLNumSandRVariables; var++) 
  {
    const int       nsrc = 1;
    const CCTK_INT  srcs[1] = {scratchvarindex};
    const CCTK_INT  tls[1] = {var};
    const CCTK_REAL facts[1] = {1.0};
    
    MoL_LinearCombination(cctkGH, SandRVariableIndex[var], rl, tl, 0.0,
                          srcs, tls, facts, nsrc);
  }

  return;
  
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
