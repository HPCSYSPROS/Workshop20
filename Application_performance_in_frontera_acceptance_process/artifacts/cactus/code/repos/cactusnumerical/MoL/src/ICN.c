 /*@@
   @file      ICN.c
   @date      Sun May 26 04:29:07 2002
   @author    Ian Hawke
   @desc 
   This implements the more efficient Iterative Crank Nicholson integrator.
   This follows the implementation of ICN in all AEI codes and is 
   equivalent to (but hopefully more efficient than) the generic ICN
   integrator also implemented by MoL.
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_ICN_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_ICNAdd(CCTK_ARGUMENTS);
void MoL_ICNAverage(CCTK_ARGUMENTS);

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
   @routine    MoL_ICNAdd
   @date       Sun May 26 04:17:23 2002
   @author     Ian Hawke
   @desc 
   Performs Iterative Crank Nicholson time integration. The number of
   steps is arbitrary.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ICNAdd(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;

  CCTK_INT index, var;
  CCTK_INT totalsize;
  CCTK_REAL *OldVar;
  CCTK_REAL *UpdateVar;
  CCTK_REAL *RHSVar;
  
#ifdef MOLDEBUG
  printf("Inside ICN.\nProcessor %d.\nStep %d.\nRefinement %d.\n"
         "Timestep %g.\nSpacestep %g.\nTime %g\n",
         CCTK_MyProc(cctkGH),
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

#ifdef MOLDEBUG
  printf("MoL: the ICN routine says dt = %f.\n", CCTK_DELTA_TIME);
#endif
  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    OldVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 1,
                                          EvolvedVariableIndex[var]);
    UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                             EvolvedVariableIndex[var]);
    RHSVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0, 
                                          RHSVariableIndex[var]);
    
#pragma omp parallel for
    for (index = 0; index < totalsize; index++)
    {
      UpdateVar[index] = OldVar[index] + CCTK_DELTA_TIME * RHSVar[index];
    }
  }
  
  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    OldVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 1,
                                             EvolvedArrayVariableIndex[var]);
    UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                             EvolvedArrayVariableIndex[var]);
    RHSVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0, 
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

/*     CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,  */
/*                  "This proc array total size is %d.",  */
/*                  arraytotalsize); */

#pragma omp parallel for
    for (index = 0; index < arraytotalsize; index++)
    {
      UpdateVar[index] = OldVar[index] + CCTK_DELTA_TIME * RHSVar[index];
    }
  }

  return;

}

 /*@@
   @routine    MoL_ICNAverage
   @date       Fri Jul 18 14:02:00 2003
   @author     Ian Hawke, Erik Schnetter
   @desc 
   Averages between the current and the previous time level.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_ICNAverage(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  cGroupDynamicData arraydata;
  CCTK_INT groupindex, ierr;
  CCTK_INT arraytotalsize, arraydim;

  CCTK_INT index, var;
  CCTK_INT totalsize;
  CCTK_REAL *OldVar;
  CCTK_REAL *UpdateVar;
  
  CCTK_REAL theta;

  theta = ICN_avg_theta;
  if (ICN_avg_swapped && (*MoL_Intermediate_Step%2))
  {
    theta = 1.0 - theta;
  }

#ifdef MOLDEBUG
  printf("Inside ICN.\nProcessor %d.\nStep %d.\nRefinement %d.\n"
         "Timestep %g.\nSpacestep %g.\nTime %g Theta %g\n",
         CCTK_MyProc(cctkGH),
         MoL_Intermediate_Steps - *MoL_Intermediate_Step + 1,
         *cctk_levfac,
         CCTK_DELTA_TIME,
         CCTK_DELTA_SPACE(0),
         cctk_time,
         theta);
#endif  

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }

#ifdef MOLDEBUG
  printf("MoL: the ICN routine says dt = %f.\n", CCTK_DELTA_TIME);
#endif
  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    OldVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 1,
                                          EvolvedVariableIndex[var]);
    UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
                                             EvolvedVariableIndex[var]);
    
#pragma omp parallel for
    for (index = 0; index < totalsize; index++)
    {
/*       UpdateVar[index] = 0.5 * (UpdateVar[index] + OldVar[index]); */
      UpdateVar[index] = (1.0 - theta) * UpdateVar[index] + 
                                theta  * OldVar[index];
    }
  }
  
  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    OldVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 1,
                                             EvolvedArrayVariableIndex[var]);
    UpdateVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, 0,
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

#pragma omp parallel for
    for (index = 0; index < arraytotalsize; index++)
    {
/*       UpdateVar[index] = 0.5 * (UpdateVar[index] + OldVar[index]); */
      UpdateVar[index] = (1.0 - theta) * UpdateVar[index] + 
                                theta  * OldVar[index];
    }
  }

  return;

}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
