 /*@@
   @file      GenericRK.c
   @date      Sun May 26 03:47:15 2002
   @author    Ian Hawke
   @desc 
   This routine performs a generic Runge-Kutta type integration
   given the set of coefficients defined in the RKAlphaCoefficients
   and RKBetaCoefficients arrays. See the article by Shu referenced
   in the documentation for more details.
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>

#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_GenericRK_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

CCTK_INT AlphaIndex(CCTK_INT Step_Number, CCTK_INT Scratch_Level);

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_GenericRKAdd(CCTK_ARGUMENTS);

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
   @routine    MoL_GenericRKAdd
   @date       Sun May 26 03:50:44 2002
   @author     Ian Hawke
   @desc 
   Performs a single step of a generic Runge-Kutta type time
   integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_GenericRKAdd(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  CCTK_INT arraytotalsize, arraydim;

  static CCTK_INT scratchspace_firstindex = -99;
  CCTK_INT index, var, scratchstep, alphaindex, scratchindex;
  CCTK_INT totalsize;
  CCTK_REAL alpha, beta;
  CCTK_REAL * restrict UpdateVar;
  CCTK_REAL const * restrict RHSVar;
  CCTK_REAL * restrict ScratchVar;

  CCTK_INT arrayscratchlocation;

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }

  if (scratchspace_firstindex == -99)
  {
    scratchspace_firstindex = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
  }

  beta = RKBetaCoefficients[MoL_Intermediate_Steps - 
                           (*MoL_Intermediate_Step)];

  /* Real GFs */

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSVariableIndex[var]);
/* #define MOLDEBUG 1 */
#ifdef MOLDEBUG
    printf("In generic RK. Variable %d (%s). RHS %d (%s). beta %g.\n",
           EvolvedVariableIndex[var],
           CCTK_VarName(EvolvedVariableIndex[var]),
           RHSVariableIndex[var],
           CCTK_VarName(RHSVariableIndex[var]),
           beta);
#endif

#pragma omp parallel for
    for (index = 0; index < totalsize; index++)
    {
      UpdateVar[index] = (*Original_Delta_Time) / cctkGH->cctk_timefac * beta * RHSVar[index];
#ifdef MOLDEBUG
      if (CCTK_EQUALS(verbose,"extreme"))
      {
        printf("Variable: %d. Index: %d. dt: %f. beta %f. RHS: %f. q: %f.\n",
               var, index, (*Original_Delta_Time) / cctkGH->cctk_timefac, beta, RHSVar[index], 
               UpdateVar[index]);
      }
#endif
    }
    
    for (scratchstep = 0; 
         scratchstep < MoL_Intermediate_Steps - (*MoL_Intermediate_Step) + 1;
         scratchstep++)
    {

      alphaindex = AlphaIndex(*MoL_Intermediate_Step, scratchstep);
      scratchindex = scratchstep - 1;

      alpha = RKAlphaCoefficients[alphaindex];
#ifdef MOLDEBUG
      printf("In generic RK. Variable %d (%s). RHS %d (%s). step %d. alpha %g.\n",
             EvolvedVariableIndex[var],
             CCTK_VarName(EvolvedVariableIndex[var]),
             RHSVariableIndex[var],
             CCTK_VarName(RHSVariableIndex[var]),
             scratchstep,
             alpha);
#endif

      if (scratchstep) 
      {
        /*
          The following would work if all drivers considered
          a vector group to have contiguous storage.
        */
        /*
          ScratchVar = &ScratchSpace[(var * MoL_Num_Scratch_Levels + 
                                   scratchindex) * totalsize];
        */
        ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                      scratchspace_firstindex
                                      + scratchindex);
#ifdef MOLDEBUG
        if (CCTK_EQUALS(verbose,"extreme"))
        {
          printf("Reading from scratch space var %d, initial index %d\n", 
                 var, scratchindex);
        }
#endif
      }
      else
      {
        ScratchVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 1, 
                                                  EvolvedVariableIndex[var]);
      }
      
      if ( (alpha > MoL_Tiny)||(alpha < -MoL_Tiny) )
      {
#pragma omp parallel for
        for (index = 0; index < totalsize; index++)
        {
          UpdateVar[index] += alpha * ScratchVar[index];
#ifdef MOLDEBUG
          if (CCTK_EQUALS(verbose,"extreme"))
          {
            printf("Variable: %d. Index: %d. step: %d. "
                   "alpha: %f. Scratch: %f. q: %f.\n",
                   var, index, (*MoL_Intermediate_Step), alpha, 
                   ScratchVar[index], UpdateVar[index]);
          }
#endif
        }
      }
      
    }
    
  }

  if (*MoL_Intermediate_Step > 1)
  {
    for (var = 0; var < MoLNumEvolvedVariables; var++)
    {
      UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                EvolvedVariableIndex[var]);
      /*
        The following would work if all drivers considered
        a vector group to have contiguous storage.
      */
      /*
        ScratchVar = &ScratchSpace[(var * MoL_Num_Scratch_Levels + 
                                 MoL_Intermediate_Steps -
                                 (*MoL_Intermediate_Step)) * totalsize];
      */
      ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                    scratchspace_firstindex
                                    + (MoL_Intermediate_Steps - (*MoL_Intermediate_Step)));
#ifdef MOLDEBUG
      printf("Writing to scratch space, var %d, index %d\n", 
             var, MoL_Intermediate_Steps - (*MoL_Intermediate_Step));
#endif
#pragma omp parallel for
      for (index = 0; index < totalsize; index++)
      {
        ScratchVar[index] = UpdateVar[index];
#ifdef MOLDEBUG
        if (CCTK_EQUALS(verbose,"extreme"))
        {
          printf("Variable: %d. Index: %d. step: %d. Scratch: %f.\n",
                 var, index, (*MoL_Intermediate_Step), ScratchVar[index]);
          fflush(stdout);
        }
#endif        
      }
    }
  }
  
  /* Real arrays */

/* #define MOLDEBUGARRAYS 1 */

  arrayscratchlocation = 0;

#ifdef MOLDEBUGARRAYS
  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, 
             "Array sizes are %d %d %d\n", MoL_Max_Evolved_Array_Size, arraytotalsize, singlearraysize);
#endif  

  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedArrayVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSArrayVariableIndex[var]);

    arraytotalsize = ArrayScratchSizes[var];

#pragma omp parallel for
    for (index = 0; index < arraytotalsize; index++)
    {
      UpdateVar[index] = (*Original_Delta_Time) / cctkGH->cctk_timefac * beta * RHSVar[index];
    }
    
    for (scratchstep = 0; 
         scratchstep < MoL_Intermediate_Steps - (*MoL_Intermediate_Step) + 1;
         scratchstep++)
    {

      alphaindex = AlphaIndex(*MoL_Intermediate_Step, scratchstep);
      scratchindex = scratchstep - 1;

      alpha = RKAlphaCoefficients[alphaindex];

      if (scratchstep) 
      {
        ScratchVar = &ArrayScratchSpace[scratchindex*CurrentArrayScratchSize + 
                                        arrayscratchlocation];
      }
      else
      {
        ScratchVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 1, 
                                                  EvolvedArrayVariableIndex[var]);
      }
      
      if ( (alpha > MoL_Tiny)||(alpha < -MoL_Tiny) )
      {
#pragma omp parallel for
        for (index = 0; index < arraytotalsize; index++)
        {
          UpdateVar[index] += alpha * ScratchVar[index];
#ifdef MOLDEBUGARRAYS
          if (CCTK_EQUALS(verbose,"extreme"))
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, 
                       "Variable: %d. Index: %d. step: %d. "
                       "alpha: %f. Scratch: %f. q: %f.\n",
                       var, index, (*MoL_Intermediate_Step), 
                       alpha, ScratchVar[index], UpdateVar[index]);
          }
#endif
        }
      }
      
    }
   
    arrayscratchlocation += arraytotalsize;
 
  }

  arrayscratchlocation = 0;

  if (*MoL_Intermediate_Step > 1)
  {
    for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
    {
      UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                EvolvedArrayVariableIndex[var]);
      ScratchVar = &ArrayScratchSpace[(MoL_Intermediate_Steps -
                                       (*MoL_Intermediate_Step)) * 
                                      CurrentArrayScratchSize +
/*                                      singlearraysize + */
/*                                      (MoL_Max_Evolved_Array_Size+1) +  */
                                     arrayscratchlocation];

      arraytotalsize = ArrayScratchSizes[var];


#ifdef MOLDEBUGARRAYS
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, 
                 "Writing to scratch space, initial address %ld, index %d \n", 
                 ScratchVar, (MoL_Intermediate_Steps -
                              (*MoL_Intermediate_Step)) * 
                 CurrentArrayScratchSize +
                 arrayscratchlocation);
#endif
#pragma omp parallel for
      for (index = 0; index < arraytotalsize; index++)
      {
        ScratchVar[index] = UpdateVar[index];
#ifdef MOLDEBUGARRAYS
        if (CCTK_EQUALS(verbose,"extreme"))
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, 
                     "Variable: %d. Index: %d. step: %d. Scratch: %f.\n",
                     var, index, (*MoL_Intermediate_Step), ScratchVar[index]);
        }
#endif        
      }
      arrayscratchlocation += arraytotalsize;
    }
  }
        
  return;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

CCTK_INT AlphaIndex(CCTK_INT Step_Number, CCTK_INT Scratch_Level)
{
  DECLARE_CCTK_PARAMETERS;
  
  return (MoL_Intermediate_Steps - Step_Number) * MoL_Intermediate_Steps + 
    Scratch_Level;
}
