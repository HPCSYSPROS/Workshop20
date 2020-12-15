 /*@@
   @file      RK4.c
   @date      Fri July 14, 2006
   @author    Yosef Zlochower
   @desc 
   A routine to perform RK4 evolution. Mostly copied from
   genericRK.c
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"
#include "Operators.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_RK4_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK4Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK4Add
   @date       
   @author     
   @desc 
   Performs a single step of a RK4 type time
   integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK4Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;

  static CCTK_INT scratchspace_firstindex = -99;
  CCTK_INT index, var;
  CCTK_INT totalsize;
  CCTK_REAL alpha, beta;
  CCTK_REAL * restrict UpdateVar;
  CCTK_REAL * restrict OldVar;
  CCTK_REAL const * restrict RHSVar;
  CCTK_REAL * restrict ScratchVar;

  CCTK_INT arrayscratchlocation;

  /* Keep a running total of alpha as we perform the substeps, so that
     we know the "real" alpha (including round-off errors) when we
     calculate the final result. */
  CCTK_REAL const time_rhs = 1.0;
  CCTK_REAL const old_time = 0.0;
  static CCTK_REAL time, scratch_time;
  
  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }

  if (scratchspace_firstindex == -99)
  {
    scratchspace_firstindex = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
  }

  switch (MoL_Intermediate_Steps - (*MoL_Intermediate_Step))
  {
    case 0:
      alpha = 1.0 / 3.0;
      beta  = 1.0 / 2.0;
      break;
    case 1:
      alpha = 2.0 / 3.0;
      beta  = 1.0 / 2.0;
      break;
    case 2:
      alpha = 1.0 / 3.0;
      beta  = 1.0;
      break;
    case 3:
      alpha = 1.0;
      beta  = 1.0 / 6.0;
      break;
  }

  if (MoL_Intermediate_Steps == (*MoL_Intermediate_Step))
  {
    time = 0.0;
  }
  
  /* Real GFs */

  CCTK_INT rl = 0;
  if (CCTK_IsFunctionAliased("GetRefinementLevel")) {
    rl = GetRefinementLevel(cctkGH);
  }
  CCTK_INT tl = 0;
  
  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    {
      int const nsrcs = 2;
      CCTK_INT const srcs[] =
        {EvolvedVariableIndex[var], RHSVariableIndex[var]};
      CCTK_INT const tls[] = {1, 0};
      CCTK_REAL const facts[] =
        {1.0, (*Original_Delta_Time) / cctkGH->cctk_timefac * beta};
      MoL_LinearCombination(cctkGH,
                            EvolvedVariableIndex[var], rl, tl, 0.0,
                            srcs, tls, facts, nsrcs);
      time = facts[0] * old_time + facts[1] * time_rhs;
    }

    /* scratch storage */
    if ((*MoL_Intermediate_Step) == MoL_Intermediate_Steps)
    {
      MoL_LinearCombination(cctkGH,
                            scratchspace_firstindex, rl, var, 0.0,
                            NULL, NULL, NULL, 0);
      scratch_time = 0.0;
    }

    if ((*MoL_Intermediate_Step)>1)
    {
      int const nsrcs = 1;
      CCTK_INT const srcs[] = {EvolvedVariableIndex[var]};
      CCTK_INT const tls[] = {0};
      CCTK_REAL const facts[] = {alpha};
      MoL_LinearCombination(cctkGH,
                            scratchspace_firstindex, rl, var, 1.0,
                            srcs, tls, facts, nsrcs);
      scratch_time += facts[0] * time;
    }
    else
    {
      int const nsrcs = 2;
      CCTK_INT const srcs[] =
        {scratchspace_firstindex, EvolvedVariableIndex[var]};
      CCTK_INT const tls[] = {var, 1};
      CCTK_REAL const facts[] = {1.0, -4.0/3.0};
      MoL_LinearCombination(cctkGH,
                            EvolvedVariableIndex[var], rl, tl, 1.0,
                            srcs, tls, facts, nsrcs);
      time += facts[0] * scratch_time + facts[1] * old_time;
    }
  }

  /* Real arrays */

  arrayscratchlocation = 0;

  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedArrayVariableIndex[var]);
    OldVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 1, 
                                              EvolvedArrayVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSArrayVariableIndex[var]);
    
    groupindex = CCTK_GroupIndexFromVarI(EvolvedArrayVariableIndex[var]);
    ierr = CCTK_GroupDynamicData(cctkGH, groupindex,
                                 &arraydata);
    if (ierr)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, 
                 "The driver does not return group information "
                 "for group '%s'.", 
                 CCTK_GroupName(groupindex));
    }
    arraytotalsize = 1;
    for (arraydim = 0; arraydim < arraydata.dim; arraydim++)
    {
      arraytotalsize *= arraydata.ash[arraydim];
    }

    ScratchVar = &ArrayScratchSpace[arrayscratchlocation];

#pragma omp parallel for
    for (index = 0; index < arraytotalsize; index++)
    {
      UpdateVar[index] = OldVar[index] +
        (*Original_Delta_Time) / cctkGH->cctk_timefac * beta * RHSVar[index];
    }

    if ((*MoL_Intermediate_Step) == MoL_Intermediate_Steps)
    {
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        ScratchVar[index] = 0;
      }
    }
    
    if ((*MoL_Intermediate_Step)>1)
    {
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        ScratchVar[index] += alpha * UpdateVar[index];
      }
    }
    else
    {
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        UpdateVar[index] += ScratchVar[index] - 4.0 / 3.0 * OldVar[index];
      }
    }
    arrayscratchlocation += arraytotalsize;
  }

  return;
}
