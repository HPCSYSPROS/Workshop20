 /*@@
   @file      RK3.c
   @date      Tue Jul 22 00:38:47 2003
   @author    Ian Hawke
   @desc 
   A specialized third order Runge-Kutta time integrator. This is
   the integrator that Shu refers to as the optimal TVD third 
   order method (see reference in documentation).
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"
#include "Operators.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_RK3_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK3Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK3Add
   @date       Tue Jul 22 00:39:55 2003
   @author     Ian Hawke
   @desc 
   Performs third order Runge-Kutta time integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

 @@*/

void MoL_RK3Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT arraydim;
  
  CCTK_INT var;
  CCTK_INT totalsize;
  
#ifdef MOLDEBUG
  printf("Inside RK3.\nStep %d.\nRefinement %d.\nTimestep %g.\n"
         "Spacestep %g.\nTime %g\n",
         MoL_Intermediate_Steps - *MoL_Intermediate_Step + 1,
         *cctk_levfac,
         CCTK_DELTA_TIME,
         CCTK_DELTA_SPACE(0),
         cctk_time);
#endif  

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
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

    case 3:
      {
        for (var = 0; var < MoLNumEvolvedVariables; var++)
        {
          int const nsrcs = 1;
          CCTK_INT const srcs[] = {RHSVariableIndex[var]};
          CCTK_INT const tls[] = {0};
          CCTK_REAL const facts[] = {CCTK_DELTA_TIME};
          MoL_LinearCombination(cctkGH,
                                EvolvedVariableIndex[var], rl, tl, 1.0,
                                srcs, tls, facts, nsrcs);
        }

        for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
        {
          int const nsrcs = 1;
          CCTK_INT const srcs[] = {RHSArrayVariableIndex[var]};
          CCTK_INT const tls[] = {0};
          CCTK_REAL const facts[] = {CCTK_DELTA_TIME};
          MoL_LinearCombination(cctkGH,
                                EvolvedArrayVariableIndex[var], rl, tl, 1.0,
                                srcs, tls, facts, nsrcs);
        }

        break;
      }
  
    case 2:
      {
        for (var = 0; var < MoLNumEvolvedVariables; var++)
        {
          int const nsrcs = 2;
          CCTK_INT const srcs[] =
            {EvolvedVariableIndex[var], RHSVariableIndex[var]};
          CCTK_INT const tls[] = {1, 0};
          CCTK_REAL const facts[] = {0.75, CCTK_DELTA_TIME};
          MoL_LinearCombination(cctkGH,
                                EvolvedVariableIndex[var], rl, tl, 0.25,
                                srcs, tls, facts, nsrcs);
        }

        for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
        {
          int const nsrcs = 2;
          CCTK_INT const srcs[] =
            {EvolvedArrayVariableIndex[var], RHSArrayVariableIndex[var]};
          CCTK_INT const tls[] = {1, 0};
          CCTK_REAL const facts[] = {0.75, CCTK_DELTA_TIME};
          MoL_LinearCombination(cctkGH,
                                EvolvedArrayVariableIndex[var], rl, tl, 0.25,
                                srcs, tls, facts, nsrcs);
        }
        break;
      }
    case 1:
      {
        for (var = 0; var < MoLNumEvolvedVariables; var++)
        {
          int const nsrcs = 2;
          CCTK_INT const srcs[] =
            {EvolvedVariableIndex[var], RHSVariableIndex[var]};
          CCTK_INT const tls[] = {1, 0};
          CCTK_REAL const facts[] = {1.0/3.0, CCTK_DELTA_TIME};
          MoL_LinearCombination(cctkGH,
                                EvolvedVariableIndex[var], rl, tl, 2.0/3.0,
                                srcs, tls, facts, nsrcs);
        }

        for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
        {
          int const nsrcs = 2;
          CCTK_INT const srcs[] =
            {EvolvedArrayVariableIndex[var], RHSArrayVariableIndex[var]};
          CCTK_INT const tls[] = {1, 0};
          CCTK_REAL const facts[] = {1.0/3.0, CCTK_DELTA_TIME};
          MoL_LinearCombination(cctkGH,
                                EvolvedArrayVariableIndex[var], rl, tl, 2.0/3.0,
                                srcs, tls, facts, nsrcs);
        }
        break;
      }
    default:
      {
        CCTK_ERROR("RK3 expects MoL_Intermediate_Step to be "
                   "in [1,3]. This should be caught at ParamCheck - bug Ian!");
        break;
      }
      
  }

  return;

}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
