 /*@@
   @file      InitialCopy.c
   @date      Sun May 26 04:43:06 2002
   @author    Ian Hawke
   @desc 
   Performs the initial copy from the previous timelevel to the
   current. This is required because the driver has rotated the 
   timelevels, but the physics thorns are expecting data in the
   current.
   @enddesc 
   @version   $Header$
 @@*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"
#include "Operators.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_InitialCopy_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_InitialCopy(CCTK_ARGUMENTS);

void MoL_InitRHS(CCTK_ARGUMENTS);

void MoL_FillAllLevels(CCTK_ARGUMENTS);

void MoL_ReportNumberVariables(CCTK_ARGUMENTS);

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
   @routine    MoL_InitialCopy
   @date       ???
   @author     Ian Hawke
   @desc 
   Copy the previous time level to the current time level.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_InitialCopy(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;
  CCTK_INT totalarrayscratchsize;

  CCTK_INT var;
  CCTK_INT totalsize;

  CCTK_REAL       * restrict CurrentVar;
  CCTK_REAL const * restrict PreviousVar;
  CCTK_INT StorageOn;

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }

  int rl = 0;
  if (CCTK_IsFunctionAliased("GetRefinementLevel")) {
    rl = GetRefinementLevel(cctkGH);
  }
  const int tl = 0;

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    const int       nsrc = 1;
    const CCTK_INT  srcs[1] = {EvolvedVariableIndex[var]};
    const CCTK_INT  tls[1] = {1};
    const CCTK_REAL facts[1] = {1.0};
    
    StorageOn = CCTK_QueryGroupStorageI
      (cctkGH, CCTK_GroupIndexFromVarI(EvolvedVariableIndex[var]));
    
    if (StorageOn < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for index %i", 
                 (int)EvolvedVariableIndex[var]);
      CCTK_ERROR("The index passed does not correspond to a GF.");
    }
    else if (StorageOn == 0) {
#ifdef MOLDEBUG
      printf("Aargh! Vars %d var %d index %d name %s\n",
             MoLNumEvolvedVariables, var, EvolvedVariableIndex[var],  
                 CCTK_VarName(EvolvedVariableIndex[var]));
#endif
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for GF %s", 
                 CCTK_VarName(EvolvedVariableIndex[var]));
      CCTK_ERROR("The grid function does not have storage assigned.");
    }
    
    MoL_LinearCombination(cctkGH, EvolvedVariableIndex[var], rl, tl, 0.0,
                          srcs, tls, facts, nsrc);
  }

  for (var = 0; var < MoLNumEvolvedVariablesSlow; var++)
  {
    const CCTK_INT   nsrc = 1;
    const CCTK_INT   srcs[1] = {EvolvedVariableIndexSlow[var]};
    const CCTK_INT   tls[1] = {1};
    const CCTK_REAL facts[1] = {1.0};
    
    StorageOn = CCTK_QueryGroupStorageI
      (cctkGH,
       CCTK_GroupIndexFromVarI(EvolvedVariableIndexSlow[var]));
    
    if (StorageOn < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for index %i", 
                 (int)EvolvedVariableIndexSlow[var]);
      CCTK_ERROR("The index passed does not correspond to a GF.");
    }
    else if (StorageOn == 0) {
#ifdef MOLDEBUG
      printf("Aargh! Vars %d var %d index %d name %s\n",
             MoLNumEvolvedVariablesSlow, var, EvolvedVariableIndexSlow[var],  
                 CCTK_VarName(EvolvedVariableIndexSlow[var]));
#endif
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for GF %s", 
                 CCTK_VarName(EvolvedVariableIndexSlow[var]));
      CCTK_ERROR("The grid function does not have storage assigned.");
    }
    
    MoL_LinearCombination(cctkGH, EvolvedVariableIndexSlow[var], rl, tl, 0.0,
                          srcs, tls, facts, nsrc);
  }
    
  /* Set up the array sizes */

  if (MoLNumEvolvedArrayVariables)
  {
    if (!ArrayScratchSizes)
    {
      ArrayScratchSizes = (CCTK_INT *)malloc(MoLNumEvolvedArrayVariables * sizeof(CCTK_INT));
      if (!ArrayScratchSizes)
      {
        CCTK_ERROR("Failed to allocate the array scratch sizes array.");
      }
      for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
      {
        ArrayScratchSizes[var] = -1;
      }
    }
  }
  
  totalarrayscratchsize = 0;

  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    PreviousVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 1,
                                               EvolvedArrayVariableIndex[var]);
    CurrentVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                              EvolvedArrayVariableIndex[var]);
    
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

    ArrayScratchSizes[var] = arraytotalsize;
    totalarrayscratchsize += arraytotalsize;

    if (arraytotalsize)
    {  
      if (PreviousVar && CurrentVar)
      {      
        memcpy(CurrentVar, PreviousVar, arraytotalsize * sizeof(CCTK_REAL));
      }
      else
      {
        printf("The pointers are %p (prev) and %p (curr)\n.",
               PreviousVar, CurrentVar);
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Null pointer for variable %s", 
                   CCTK_VarName(EvolvedArrayVariableIndex[var]));
      }
    }
    
  }

  if (totalarrayscratchsize > CurrentArrayScratchSize)
  {
    if (ArrayScratchSpace)
    {
      free(ArrayScratchSpace);
      ArrayScratchSpace = NULL;
    }
    ArrayScratchSpace = 
      (CCTK_REAL*)malloc(totalarrayscratchsize * 
                         MoL_Num_Scratch_Levels * 
                         sizeof(CCTK_REAL));
    for (var = 0; var < totalarrayscratchsize * MoL_Num_Scratch_Levels; var++)
    {
      ArrayScratchSpace[var] = 0.0;
    }
    CurrentArrayScratchSize = totalarrayscratchsize;
  }
  

  /* 
     Now the Save and Restore variables. Shift the data in the 
     current level to the scratch space, then do the copy 
  */

  const int scratchvarindex = CCTK_FirstVarIndex("MOL::SANDRSCRATCHSPACE");
  if (scratchvarindex < 0)
  {
    CCTK_ERROR("Internal error");
  }

  for (var = 0; var < MoLNumSandRVariables; var++)
  {
    const int       nsrc_save = 1;
    const CCTK_INT  srcs_save[1] = {SandRVariableIndex[var]};
    const CCTK_INT  tls_save[1] = {0};
    const CCTK_REAL facts_save[1] = {1.0};
    const int       nsrc = 1;
    const CCTK_INT  srcs[1] = {SandRVariableIndex[var]};
    const CCTK_INT  tls[1] = {1};
    const CCTK_REAL facts[1] = {1.0};
    
    StorageOn = CCTK_QueryGroupStorageI
      (cctkGH,
       CCTK_GroupIndexFromVarI(SandRVariableIndex[var]));
    
    if (StorageOn < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for index %i", 
                 (int)SandRVariableIndex[var]);
      CCTK_ERROR("The index passed does not correspond to a GF.");
    }
    else if (StorageOn == 0) {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for GF %s", 
                 CCTK_VarName(SandRVariableIndex[var]));
      CCTK_ERROR("The grid function does not have storage assigned.");
    }

    MoL_LinearCombination(cctkGH, scratchvarindex, rl, var, 0.0,
                          srcs_save, tls_save, facts_save, nsrc_save);
    MoL_LinearCombination(cctkGH, SandRVariableIndex[var], rl, tl, 0.0,
                          srcs, tls, facts, nsrc);
  }  

  /*
    Now do the constrained variables.
  */

  for (var = 0; var < MoLNumConstrainedVariables; var++)
  {
    const int       nsrc = 1;
    const CCTK_INT  srcs[1] = {ConstrainedVariableIndex[var]};
    const CCTK_INT  tls[1] = {1};
    const CCTK_REAL facts[1] = {1.0};
    
    StorageOn = CCTK_QueryGroupStorageI
      (cctkGH,
       CCTK_GroupIndexFromVarI(ConstrainedVariableIndex[var]));
    
    if (StorageOn < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for index %i", 
                 (int)ConstrainedVariableIndex[var]);
      CCTK_ERROR("The index passed does not correspond to a GF.");
    }
    else if (StorageOn == 0) {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for GF %s", 
                 CCTK_VarName(ConstrainedVariableIndex[var]));
      CCTK_ERROR("The grid function does not have storage assigned.");
    }

    /* 
       Check that there is more than one timelevel.
       If not, copying is unnecessary.
    */

    StorageOn = CCTK_ActiveTimeLevelsVI(cctkGH,
                                        ConstrainedVariableIndex[var]);
    
    if (StorageOn > 1)
    {
      MoL_LinearCombination(cctkGH, ConstrainedVariableIndex[var], rl, tl, 0.0,
                            srcs, tls, facts, nsrc);
    }
    
  }

  return;
}

 /*@@
   @routine    MoL_InitRHS
   @date       Tue Mar 23 2004
   @author     Erik Schnetter
   @desc 
   Initialise all RHS variables with zero.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_InitRHS(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;

  CCTK_INT var;
  CCTK_INT index;
/*   CCTK_INT i,j,k; */
  CCTK_INT totalsize;

  CCTK_REAL * restrict RHSVar;
  CCTK_INT StorageOn;

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
  
  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    StorageOn = CCTK_QueryGroupStorageI(cctkGH,
                                       CCTK_GroupIndexFromVarI(RHSVariableIndex[var]));
    
    if (StorageOn < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for index %i", 
                 (int)RHSVariableIndex[var]);
      CCTK_ERROR("The index passed does not correspond to a GF.");
    }
    else if (StorageOn == 0) {
#ifdef MOLDEBUG
      printf("Aargh! Vars %d var %d index %d name %s\n",
             MoLNumEvolvedVariables, var, RHSVariableIndex[var],  
                 CCTK_VarName(RHSVariableIndex[var]));
#endif
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for GF %s", 
                 CCTK_VarName(RHSVariableIndex[var]));
      CCTK_ERROR("The grid function does not have storage assigned.");
    }
    
    MoL_LinearCombination(cctkGH, RHSVariableIndex[var], rl, tl, 0.0,
                          NULL, NULL, NULL, 0);
  }
  
  for (var = 0; var < MoLNumEvolvedVariablesSlow; var++)
  {
    StorageOn = CCTK_QueryGroupStorageI(cctkGH,
                                       CCTK_GroupIndexFromVarI(RHSVariableIndexSlow[var]));
    
    if (StorageOn < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for index %i", 
                 (int)RHSVariableIndexSlow[var]);
      CCTK_ERROR("The index passed does not correspond to a GF.");
    }
    else if (StorageOn == 0) {
#ifdef MOLDEBUG
      printf("Aargh! Vars %d var %d index %d name %s\n",
             MoLNumEvolvedVariablesSlow, var, RHSVariableIndexSlow[var],  
                 CCTK_VarName(RHSVariableIndexSlow[var]));
#endif
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for GF %s", 
                 CCTK_VarName(RHSVariableIndexSlow[var]));
      CCTK_ERROR("The grid function does not have storage assigned.");
    }
    
    MoL_LinearCombination(cctkGH, RHSVariableIndexSlow[var], rl, tl, 0.0,
                          NULL, NULL, NULL, 0);
  }
    
  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    RHSVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                          RHSArrayVariableIndex[var]);
    
    groupindex = CCTK_GroupIndexFromVarI(RHSArrayVariableIndex[var]);
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

    if (arraytotalsize)
    {
      if (RHSVar)
      {
#pragma omp parallel for
        for (index = 0; index < arraytotalsize; index++)
        {
          RHSVar[index] = 0;
        }
      }
      else
      {
        printf("The pointer is %p (rhs)\n.",
               RHSVar);
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Null pointer for variable %s", 
                   CCTK_VarName(RHSArrayVariableIndex[var]));
      }
    }
    
  }

  return;
}

 /*@@
   @routine    MoL_FillAllLevels
   @date       Fri Apr 25 16:11:18 2003
   @author     Ian Hawke
   @desc 
   This routine is a bit of a hack, and I'm still not convinced 
   it is really necessary. It fills the previous timelevels by 
   copying the data from the current timelevels, which should have
   been set up during the CCTK_INITIAL timebin.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FillAllLevels(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT var, level;
  CCTK_INT totalsize, arraydim;

  CCTK_REAL const * restrict CurrentVar;
  CCTK_REAL       * restrict PreviousVar;

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    CurrentVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndex[var]);
    for (level = 1; 
         level < CCTK_ActiveTimeLevelsVI(cctkGH,
                                         EvolvedVariableIndex[var]); 
         level++)
    {
      PreviousVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, level, 
                                                 EvolvedVariableIndex[var]);
      if (PreviousVar)
      {   
        memcpy(PreviousVar, CurrentVar, totalsize * sizeof(CCTK_REAL));
      }
    }
  }
  

  for (var = 0; var < MoLNumEvolvedVariablesSlow; var++)
  {
    CurrentVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndexSlow[var]);
    for (level = 1; 
         level < CCTK_ActiveTimeLevelsVI(cctkGH,
                                         EvolvedVariableIndexSlow[var]); 
         level++)
    {
      PreviousVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, level, 
                                                 EvolvedVariableIndexSlow[var]);
      if (PreviousVar)
      {   
        memcpy(PreviousVar, CurrentVar, totalsize * sizeof(CCTK_REAL));
      }
    }
  }
  

  for (var = 0; var < MoLNumConstrainedVariables; var++)
  {
    CurrentVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                              ConstrainedVariableIndex[var]);
    for (level = 1; level < CCTK_ActiveTimeLevelsVI(cctkGH,
                                                    ConstrainedVariableIndex[var]); level++)
    {
      PreviousVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, level, 
                                                 ConstrainedVariableIndex[var]);
      if (PreviousVar)
      {   
        memcpy(PreviousVar, CurrentVar, totalsize * sizeof(CCTK_REAL));
      }
    }
  }
  

  for (var = 0; var < MoLNumSandRVariables; var++)
  {
    CurrentVar = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, 0, 
                                              SandRVariableIndex[var]);
    for (level = 1; level < CCTK_ActiveTimeLevelsVI(cctkGH,
                                                    SandRVariableIndex[var]); 
         level++)
    {
      PreviousVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, level, 
                                                 SandRVariableIndex[var]);
      if (PreviousVar)
      {   
        memcpy(PreviousVar, CurrentVar, totalsize * sizeof(CCTK_REAL));
      }
    }
  }

  return;
}

 /*@@
   @routine    MoL_ReportNumberVariables
   @date       Thu Jul 17 18:35:54 2003
   @author     Ian Hawke
   @desc 
   Prints some useful information about the number of
   registered functions.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ReportNumberVariables(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT var;

  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of evolved variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,(int)MoLNumEvolvedVariables);

  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of slow evolved variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,(int)MoLNumEvolvedVariablesSlow);

  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of constrained variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,(int)MoLNumConstrainedVariables);
  
  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of SandR variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,(int)MoLNumSandRVariables);
  

  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of evolved array variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,
             (int)MoLNumEvolvedArrayVariables);
  
  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of constrained array variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,
             (int)MoLNumConstrainedArrayVariables);
  
  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum number of SandR array variables is %d. "
             "%d are registered.",
             (int)MoLMaxNumRegisteredVariables,
             (int)MoLNumSandRArrayVariables);
  

  CCTK_VInfo(CCTK_THORNSTRING,
             "The maximum size of any array variables is %d.",
             (int)MoL_Max_Evolved_Array_Size);

  if (CCTK_Equals(verbose, "register"))
  {

    if (MoLNumEvolvedVariables)
    {
      CCTK_INFO("The evolved variables are:");
      for (var = 0; var < MoLNumEvolvedVariables; var++)
      {
        CCTK_VInfo(CCTK_THORNSTRING,"   %d   :   %s",
                   (int)var, CCTK_VarName(EvolvedVariableIndex[var]));
      }
    }

    if (MoLNumConstrainedVariables)
    {    
      CCTK_INFO("The constrained variables are:");
      for (var = 0; var < MoLNumConstrainedVariables; var++)
      {
        CCTK_VInfo(CCTK_THORNSTRING,"   %d   :   %s",
                   (int)var, CCTK_VarName(ConstrainedVariableIndex[var]));
      }
    }
    
    if (MoLNumSandRVariables)
    {
      CCTK_INFO("The save and restore variables are:");
      for (var = 0; var < MoLNumSandRVariables; var++)
      {
        CCTK_VInfo(CCTK_THORNSTRING,"   %d   :   %s",
                   (int)var, CCTK_VarName(SandRVariableIndex[var]));
      }
    }
    
    if (MoLNumEvolvedArrayVariables)
    {
      CCTK_INFO("The evolved array variables are:");
      for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
      {
        CCTK_VInfo(CCTK_THORNSTRING,"   %d   :   %s",
                   (int)var, CCTK_VarName(EvolvedArrayVariableIndex[var]));
      }
    }

    if (MoLNumConstrainedArrayVariables)
    {    
      CCTK_INFO("The constrained array variables are:");
      for (var = 0; var < MoLNumConstrainedArrayVariables; var++)
      {
        CCTK_VInfo(CCTK_THORNSTRING,"   %d   :   %s",
                   (int)var, CCTK_VarName(ConstrainedArrayVariableIndex[var]));
      }
    }
    
    if (MoLNumSandRArrayVariables)
    {
      CCTK_INFO("The save and restore array variables are:");
      for (var = 0; var < MoLNumSandRArrayVariables; var++)
      {
        CCTK_VInfo(CCTK_THORNSTRING,"   %d   :   %s",
                   (int)var, CCTK_VarName(SandRArrayVariableIndex[var]));
      }
    }
    
  }

  return;
}
