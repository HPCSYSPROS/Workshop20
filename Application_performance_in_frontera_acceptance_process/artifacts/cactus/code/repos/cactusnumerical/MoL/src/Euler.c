 /*@@
   @file      Euler.c
   @date      2012-12-05
   @author    Erik Schnetter
   @desc 
   An explicit Euler time integrator. This method is unstable in most
   cases, but is sometimes useful for debugging.
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"
#include "Operators.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_Euler_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_EulerAdd(CCTK_ARGUMENTS);

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
   @routine    MoL_EulerAdd
   @date       2012-12-05
   @author     Erik Schnetter
   @desc 
   Performs first order Euler time integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_EulerAdd(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int totalsize = 1;
  for (int arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }
  
  CCTK_INT rl = 0;
  if (CCTK_IsFunctionAliased("GetRefinementLevel")) {
    rl = GetRefinementLevel(cctkGH);
  }
  CCTK_INT tl = 0;
  
  switch (*MoL_Intermediate_Step)
  {
  case 1:
    {
      for (int var = 0; var < MoLNumEvolvedVariables; var++)
      {
        int const nsrcs = 2;
        CCTK_INT const srcs[] =
          {EvolvedVariableIndex[var], RHSVariableIndex[var]};
        CCTK_INT const tls[] = {1, 0};
        CCTK_REAL const facts[] = {1.0, CCTK_DELTA_TIME};
        MoL_LinearCombination(cctkGH,
                              EvolvedVariableIndex[var], rl, tl, 0.0,
                              srcs, tls, facts, nsrcs);
      }
      
      for (int var = 0; var < MoLNumEvolvedArrayVariables; var++)
      {
        int const nsrcs = 2;
        CCTK_INT const srcs[] =
          {EvolvedArrayVariableIndex[var], RHSArrayVariableIndex[var]};
        CCTK_INT const tls[] = {1, 0};
        CCTK_REAL const facts[] = {1.0, CCTK_DELTA_TIME};
        MoL_LinearCombination(cctkGH,
                              EvolvedArrayVariableIndex[var], rl, tl, 0.0,
                              srcs, tls, facts, nsrcs);
      }
      
      break;
    }
  default:
    {
      CCTK_ERROR("Euler expects MoL_Intermediate_Step to be in [1,1]. "
                 "This should be caught at ParamCheck - bug Ian!");
    }
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
