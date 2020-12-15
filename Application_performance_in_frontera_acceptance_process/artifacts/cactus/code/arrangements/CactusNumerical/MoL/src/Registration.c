 /*@@
   @file      Registration.c
   @date      Thu May 30 11:21:44 2002
   @author    Ian Hawke
   @desc 
   The external functions called (via function aliasing) by physics
   thorns to tell MoL that they want these GFs to be treated as a
   given type.
   @enddesc 
   @version   $Header$
 @@*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <stdlib.h>
#ifdef MOLDEBUG
#include <stdio.h>
#endif

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_Registration_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
static char *VarListToString(const CCTK_INT varindices[], const int nvars);

/* support old versions of the flesh */
#ifndef HAVE_CCTK_DECLARED_TIMELEVELS
#define CCTK_DeclaredTimeLevelsVI(vi) CCTK_MaxTimeLevelsVI(vi)
#endif

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_SetScheduleStatus(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

CCTK_INT MoL_RegisterEvolved(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);

CCTK_INT MoL_RegisterEvolvedSlow(CCTK_INT EvolvedIndex, CCTK_INT RHSIndexSlow);

CCTK_INT MoL_RegisterConstrained(CCTK_INT ConstrainedIndex);

CCTK_INT MoL_RegisterSaveAndRestore(CCTK_INT SandRIndex);

CCTK_INT MoL_RegisterEvolvedGroup(CCTK_INT EvolvedGroupIndex, 
                                  CCTK_INT RHSGroupIndex);

CCTK_INT MoL_RegisterConstrainedGroup(CCTK_INT ConstrainedGroupIndex);

CCTK_INT MoL_RegisterSaveAndRestoreGroup(CCTK_INT SandRGroupIndex);


CCTK_INT MoL_RegisterEvolvedReal(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);

CCTK_INT MoL_RegisterEvolvedRealSlow(CCTK_INT EvolvedIndex, CCTK_INT RHSIndexSlow);

CCTK_INT MoL_RegisterConstrainedReal(CCTK_INT ConstrainedIndex);

CCTK_INT MoL_RegisterSaveAndRestoreReal(CCTK_INT SandRIndex);

CCTK_INT MoL_RegisterEvolvedRealGroup(CCTK_INT EvolvedGroupIndex, 
                                  CCTK_INT RHSGroupIndex);

CCTK_INT MoL_RegisterEvolvedRealGroupSlow(CCTK_INT EvolvedGroupIndex, 
                                          CCTK_INT RHSGroupIndexSlow);


CCTK_INT MoL_RegisterConstrainedRealGroup(CCTK_INT ConstrainedGroupIndex);

CCTK_INT MoL_RegisterSaveAndRestoreRealGroup(CCTK_INT SandRGroupIndex);


CCTK_INT MoL_RegisterEvolvedArray(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);

CCTK_INT MoL_RegisterConstrainedArray(CCTK_INT ConstrainedIndex);

CCTK_INT MoL_RegisterSaveAndRestoreArray(CCTK_INT SandRIndex);

CCTK_INT MoL_RegisterEvolvedArrayGroup(CCTK_INT EvolvedGroupIndex, 
                                       CCTK_INT RHSGroupIndex);

CCTK_INT MoL_RegisterConstrainedArrayGroup(CCTK_INT ConstrainedGroupIndex);

CCTK_INT MoL_RegisterSaveAndRestoreArrayGroup(CCTK_INT SandRGroupIndex);


CCTK_INT MoL_QueryEvolvedRHS(CCTK_INT EvolvedIndex);
 
CCTK_INT MoL_NumIntegratorSubsteps(void);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    MoL_SetScheduleStatus
   @date       Mon Jan 31 16:01:51 2005
   @author     Ian Hawke
   @desc 
   Set the ScheduleStatus flag. This is to ensure that we can catch
   calls to the registration routines that happen too early.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


void MoL_SetScheduleStatus(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;

  ScheduleStatus = 1;
}


 /*@@
   @routine    MoL_RegisterEvolved
   @date       Thu May 30 11:36:59 2002
   @author     Ian Hawke
   @desc 
   Given the index of the GF to be evolved and the RHS GF, it stores 
   the indexes for later use together with various error checking.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterEvolved(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex)
{

  CCTK_INT retval, GroupType;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register variable '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_VarName(EvolvedIndex));
    retval++;
  }

  GroupType = CCTK_GroupTypeFromVarI(EvolvedIndex);
  if (GroupType < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved index %i is not a valid variable index.",
               (int)EvolvedIndex);
    retval++;
  }

  if (!retval)
  {
    switch (CCTK_VarTypeI(EvolvedIndex))
    {
      case CCTK_VARIABLE_REAL:
        {
          switch (GroupType)
          {
            case CCTK_GF:
              {
                retval +=
                  MoL_RegisterEvolvedReal(EvolvedIndex, 
                                          RHSIndex);
                break;
              }
            case CCTK_ARRAY:
            case CCTK_SCALAR:
              {
                retval +=
                  MoL_RegisterEvolvedArray(EvolvedIndex, 
                                           RHSIndex);
                break;
              }
            default:
              {
                CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                           "The variable '%s' is not a GF or an array", 
                           CCTK_VarName(EvolvedIndex));
                retval++;
                break;
              }
          }
          break;
        }
      default:
        {
          CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                     "The variable '%s' is not REAL.", 
                     CCTK_VarName(EvolvedIndex));
          retval++;
          break;
        }
    }  
  }
  
  return retval;

}


 /*@@
   @routine    MoL_RegisterEvolvedSlow
   @date       Thu May 30 11:36:59 2002
   @author     Ian Hawke
   @desc 
   Given the index of the GF to be evolved and the RHS GF, it stores 
   the indexes for later use together with various error checking.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterEvolvedSlow(CCTK_INT EvolvedIndex, CCTK_INT RHSIndexSlow)
{

  CCTK_INT retval, GroupType;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register variable '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_VarName(EvolvedIndex));
    retval++;
  }

  GroupType = CCTK_GroupTypeFromVarI(EvolvedIndex);
  if (GroupType < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved index %i is not a valid variable index.",
               (int)EvolvedIndex);
    retval++;
  }

  if (!retval)
  {
    switch (CCTK_VarTypeI(EvolvedIndex))
    {
      case CCTK_VARIABLE_REAL:
        {
          switch (GroupType)
          {
            case CCTK_GF:
              {
                retval +=
                  MoL_RegisterEvolvedRealSlow(EvolvedIndex,
                                              RHSIndexSlow);
                break;
              }
            case CCTK_ARRAY:
            case CCTK_SCALAR:
              {
                CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                           "The variable '%s' is not supported for multirate methods", 
                           CCTK_VarName(EvolvedIndex));
                retval++;
                break;
              }
            default:
              {
                CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                           "The variable '%s' is not supported for multirate methods", 
                           CCTK_VarName(EvolvedIndex));
                retval++;
                break;
              }
          }
          break;
        }
      default:
        {
          CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                     "The variable '%s' is not REAL.", 
                     CCTK_VarName(EvolvedIndex));
          retval++;
          break;
        }
    }  
  }
  
  return retval;

}



 /*@@
   @routine    MoL_RegisterConstrained
   @date       Thu May 30 12:35:58 2002
   @author     Ian Hawke
   @desc 
   Given the index of the GF, register it as a constrained variable.
   If there's only one timelevel then ignore it as there will be no
   rotation and so MoL doesn't have to do anything.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterConstrained(CCTK_INT ConstrainedIndex)
{
  
  CCTK_INT retval, GroupType;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register variable '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_VarName(ConstrainedIndex));
    retval++;
  }

  GroupType = CCTK_GroupTypeFromVarI(ConstrainedIndex);
  if (GroupType < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Constrained index %i is not a valid variable index.",
               (int)ConstrainedIndex);
    retval++;
    return retval;
  }
  
  switch (CCTK_VarTypeI(ConstrainedIndex))
  {
    case CCTK_VARIABLE_REAL:
      {
        switch (GroupType)
        {
          case CCTK_GF:
            {
              retval +=
                MoL_RegisterConstrainedReal(ConstrainedIndex);
              break;
            }
          case CCTK_ARRAY:
          case CCTK_SCALAR:
            {
              retval +=
                MoL_RegisterConstrainedArray(ConstrainedIndex);
              break;
            }
          default:
            {
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The variable '%s' is not a GF or an array", 
                         CCTK_VarName(ConstrainedIndex));
              retval++;
              break;
            }
        }
        break;
      }
    default:
      {
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "The variable '%s' is not REAL.", 
                   CCTK_VarName(ConstrainedIndex));
        retval++;
        break;
      }
  }  
  
  return retval;

}

 /*@@
   @routine    MoL_RegisterSaveAndRestore
   @date       Thu May 30 12:37:40 2002
   @author     Ian Hawke
   @desc 
   Given a GF index store it for later use as a save and restore type.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterSaveAndRestore(CCTK_INT SandRIndex)
{
  
  CCTK_INT retval, GroupType;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register variable '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_VarName(SandRIndex));
    retval++;
  }

  GroupType = CCTK_GroupTypeFromVarI(SandRIndex);
  if (GroupType < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "SaveAndRestore index %i is not a valid variable index.",
               (int)SandRIndex);
    retval++;
    return retval;
  }
  
  switch (CCTK_VarTypeI(SandRIndex))
  {
    case CCTK_VARIABLE_REAL:
      {
        switch (GroupType)
        {
          case CCTK_GF:
            {
              retval +=
                MoL_RegisterSaveAndRestoreReal(SandRIndex);
              break;
            }
          case CCTK_ARRAY:
          case CCTK_SCALAR:
            {
              retval +=
                MoL_RegisterSaveAndRestoreArray(SandRIndex);
              break;
            }
          default:
            {
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The variable '%s' is not a GF or an array", 
                         CCTK_VarName(SandRIndex));
              retval++;
              break;
            }
        }
        break;
      }
    default:
      {
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "The variable '%s' is not REAL.", 
                   CCTK_VarName(SandRIndex));
        retval++;
        break;
      }
  }  
  
  return retval;
  
}

CCTK_INT MoL_RegisterEvolvedGroup(CCTK_INT EvolvedGroupIndex, 
                                  CCTK_INT RHSGroupIndex)
{

  CCTK_INT retval, GroupFirstVar;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register group '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_GroupName(EvolvedGroupIndex));
    retval++;
  }

  GroupFirstVar = CCTK_FirstVarIndexI(EvolvedGroupIndex);
  if (GroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved group index %i is not a valid group index.",
               (int)EvolvedGroupIndex);
    retval++;
  }
  
  switch (CCTK_VarTypeI(CCTK_FirstVarIndexI(EvolvedGroupIndex)))
  {
    case CCTK_VARIABLE_REAL:
      {
        switch (CCTK_GroupTypeI(EvolvedGroupIndex))
        {
          case CCTK_GF:
            {
              retval +=
                MoL_RegisterEvolvedRealGroup(EvolvedGroupIndex, 
                                             RHSGroupIndex);
              break;
            }
          case CCTK_ARRAY:
          case CCTK_SCALAR:
            {
              retval +=
                MoL_RegisterEvolvedArrayGroup(EvolvedGroupIndex, 
                                              RHSGroupIndex);
              break;
            }
          default:
            {
              char *groupname = CCTK_GroupName(EvolvedGroupIndex);
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The group '%s' is not a GF or an array", 
                         groupname);
              free(groupname);
              retval++;
              break;
            }
        }
        break;
      }
    default:
      {
        char *groupname = CCTK_GroupName(EvolvedGroupIndex);
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "The group '%s' is not REAL.", 
                   groupname);
        free(groupname);
        retval++;
        break;
      }
  }  
  
  return retval;

}


CCTK_INT MoL_RegisterEvolvedGroupSlow(CCTK_INT EvolvedGroupIndex, 
                                      CCTK_INT RHSGroupIndex2)
{

  CCTK_INT retval, GroupFirstVar;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register group '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_GroupName(EvolvedGroupIndex));
    retval++;
  }

  GroupFirstVar = CCTK_FirstVarIndexI(EvolvedGroupIndex);
  if (GroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved group index %i is not a valid group index.",
               (int)EvolvedGroupIndex);
    retval++;
  }
  
  switch (CCTK_VarTypeI(CCTK_FirstVarIndexI(EvolvedGroupIndex)))
  {
    case CCTK_VARIABLE_REAL:
      {
        switch (CCTK_GroupTypeI(EvolvedGroupIndex))
        {
          case CCTK_GF:
            {
              retval +=
                MoL_RegisterEvolvedRealGroupSlow(EvolvedGroupIndex, 
                                                 RHSGroupIndex2);
              break;
            }
          case CCTK_ARRAY:
          case CCTK_SCALAR:
            {
              char *groupname = CCTK_GroupName(EvolvedGroupIndex);
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The group '%s' is not supported by multirate RK", 
                         groupname);
              free(groupname);
              retval++;
              break;
            }
          default:
            {
              char *groupname = CCTK_GroupName(EvolvedGroupIndex);
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The group '%s' is not supported by multirate RK", 
                         groupname);
              free(groupname);
              retval++;
              break;
            }
        }
        break;
      }
    default:
      {
        char *groupname = CCTK_GroupName(EvolvedGroupIndex);
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "The group '%s' is not REAL.", 
                   groupname);
        free(groupname);
        retval++;
        break;
      }
  }  
  
  return retval;

}



CCTK_INT MoL_RegisterConstrainedGroup(CCTK_INT ConstrainedGroupIndex)
{

  CCTK_INT retval, GroupFirstVar;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register group '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_VarName(ConstrainedGroupIndex));
    retval++;
  }

  GroupFirstVar = CCTK_FirstVarIndexI(ConstrainedGroupIndex);
  if (GroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Constrained group index %i is not a valid group index.",
               (int)ConstrainedGroupIndex);
    retval++;
  }
  
  switch (CCTK_VarTypeI(CCTK_FirstVarIndexI(ConstrainedGroupIndex)))
  {
    case CCTK_VARIABLE_REAL:
      {
        switch (CCTK_GroupTypeI(ConstrainedGroupIndex))
        {
          case CCTK_GF:
            {
              retval +=
                MoL_RegisterConstrainedRealGroup(ConstrainedGroupIndex);
              break;
            }
          case CCTK_ARRAY:
            {
              retval +=
                MoL_RegisterConstrainedArrayGroup(ConstrainedGroupIndex);
              break;
            }
          default:
            {
              char *groupname = CCTK_GroupName(ConstrainedGroupIndex);
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The group '%s' is not a GF or an array", 
                         groupname);
              free(groupname);
              retval++;
              break;
            }
        }
        break;
      }
    default:
      {
        char *groupname = CCTK_GroupName(ConstrainedGroupIndex);
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "The group '%s' is not REAL.", 
                   groupname);
        free(groupname);
        retval++;
        break;
      }
  }  
  
  return retval;

}

CCTK_INT MoL_RegisterSaveAndRestoreGroup(CCTK_INT SandRGroupIndex)
{

  CCTK_INT retval, GroupFirstVar;
  
  retval = 0;

  if (!ScheduleStatus)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MoL registration routine called too early!\n"
               "Trying to register group '%s',"
               "Please ensure that all calls to MoL registration routines "
               "occur within the \"MoL_Register\" timebin.",
               CCTK_VarName(SandRGroupIndex));
    retval++;
  }

  GroupFirstVar = CCTK_FirstVarIndexI(SandRGroupIndex);
  if (GroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "SaveAndRestore group index %i is not a valid group index.",
               (int)SandRGroupIndex);
    retval++;
  }
  
  switch (CCTK_VarTypeI(CCTK_FirstVarIndexI(SandRGroupIndex)))
  {
    case CCTK_VARIABLE_REAL:
      {
        switch (CCTK_GroupTypeI(SandRGroupIndex))
        {
          case CCTK_GF:
            {
              retval +=
                MoL_RegisterSaveAndRestoreRealGroup(SandRGroupIndex);
              break;
            }
          case CCTK_ARRAY:
            {
              retval +=
                MoL_RegisterSaveAndRestoreArrayGroup(SandRGroupIndex);
              break;
            }
          default:
            {
              char *groupname = CCTK_GroupName(SandRGroupIndex);
              CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "The group '%s' is not a GF or an array", 
                         groupname);
              free(groupname);
              retval++;
              break;
            }
        }
        break;
      }
    default:
      {
        char *groupname = CCTK_GroupName(SandRGroupIndex);
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "The group '%s' is not REAL.", 
                   groupname);
        free(groupname);
        retval++;
        break;
      }
  }  
  
  return retval;

}

/*
  Temporary no-op functions for testing.
*/

 /*@@
   @routine    MoL_RegisterEvolved
   @date       Thu May 30 11:36:59 2002
   @author     Ian Hawke
   @desc 
   Given the index of the GF to be evolved and the RHS GF, it stores 
   the indexes for later use together with various error checking.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterEvolvedReal(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_INT /*  ierr, */ index, varused, numtimelevs1, numtimelevs2;
  
#ifdef MOLDEBUG
  printf("Arrived in MoLRegisterEvolved \n");
  printf("The indexes are %d and %d.\n",EvolvedIndex, RHSIndex);
  printf("These correspond to variables %s and %s.\n",
         CCTK_VarName(EvolvedIndex),CCTK_VarName(RHSIndex));
  printf("The pointer to EvolvedVariableIndex: %p\n",
         EvolvedVariableIndex);
#endif

  if (!(CCTK_GroupTypeFromVarI(EvolvedIndex)==CCTK_GF))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not a GF and so should "
               "not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(EvolvedIndex)); 
  }

  if (!(CCTK_GroupTypeFromVarI(RHSIndex)==CCTK_GF))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The rhs variable %s is not a GF and so should "
               "not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(RHSIndex)); 
  }    

  if (!(CCTK_VarTypeI(EvolvedIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not of type CCTK_REAL and so "
               "should not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(EvolvedIndex)); 
  }

  if (!(CCTK_VarTypeI(RHSIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The rhs variable %s is not of type CCTK_REAL and "
               "so should not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(RHSIndex)); 
  }

  numtimelevs1 = CCTK_DeclaredTimeLevelsVI(EvolvedIndex);
  numtimelevs2 = CCTK_DeclaredTimeLevelsVI(RHSIndex);
  
  if ( (numtimelevs1 < 0) || (numtimelevs2 < 0) ) 
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for variable index %i", (int)EvolvedIndex); 
    CCTK_ERROR("The index passed does not correspond to a GF.");
  }

  if (numtimelevs1 < 2)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for variable index %i name %s", 
               (int)EvolvedIndex, CCTK_VarName(EvolvedIndex));
    CCTK_ERROR("The GF passed only has one timelevel. "
              "It must have at least two.");
  }

  varused = 0;

  for (index = 0; (index < MoLNumEvolvedVariables)&&(!varused); index++)
  {
    varused = (EvolvedIndex == EvolvedVariableIndex[index]);
#ifdef MOLDEBUG
    printf("Registering %d. Checking index %d which is %d\n",EvolvedIndex,
           index,EvolvedVariableIndex[index]);
#endif
  }

  if (varused)
  {

    CCTK_VWarn(2,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The GF %s has already been registered "
               "as an evolved variable with RHS GF %s. "
               "The attempt to register with RHS GF %s will be ignored",
               CCTK_VarName(EvolvedIndex),
               CCTK_VarName(RHSVariableIndex[index-1]),
               CCTK_VarName(RHSIndex));

  }
  else
  {

    if (MoLNumEvolvedVariables+1 > MoLMaxNumRegisteredVariables)
    {
      const char *evolvedvarname = CCTK_VarName(EvolvedIndex);
      char *registeredvars = 
        VarListToString(EvolvedVariableIndex, MoLNumEvolvedVariables);
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "You tried to register more than %d evolved variables. Since "
                 "this should be the total number of Cactus variables, this "
                 "looks like a bug. You tried to register '%s'. The "
                 "registered variables are: %s.",
                 (int)MoLMaxNumRegisteredVariables,
                 evolvedvarname, registeredvars);
      free(registeredvars); /* NOTREACHED */
    }
  
    EvolvedVariableIndex[MoLNumEvolvedVariables] = EvolvedIndex;
    RHSVariableIndex[MoLNumEvolvedVariables] = RHSIndex;
  
    MoLNumEvolvedVariables++;

#ifdef MOLDEBUG
    printf("The max number is now %d. Just added %d (%s).\n",
           MoLNumEvolvedVariables, EvolvedIndex,
           CCTK_VarName(EvolvedIndex));
#endif

  }

  varused = -1;

  for (index = 0; (index < MoLNumConstrainedVariables)&&(!(varused+1)); 
       index++)
  {
    if (EvolvedIndex == ConstrainedVariableIndex[index])
    {
      varused = index;
    }
    
  }
  
  if ((varused+1))
  {
    for (index = varused; index < MoLNumConstrainedVariables-1; index++)
    {
      ConstrainedVariableIndex[index] = ConstrainedVariableIndex[index+1];
    }
    MoLNumConstrainedVariables--;
  }
    
  varused = -1;

  for (index = 0; (index < MoLNumSandRVariables)&&(!(varused+1)); index++)
  {
    if (EvolvedIndex == SandRVariableIndex[index])
    {
      varused = index;
    }
    
#ifdef MOLDEBUG
    printf("Checking SandR var %d. Index %d (evolvedindex %d).\n",
           index, SandRVariableIndex[index], EvolvedIndex);    
#endif

  }
  
  if ((varused+1))
  {
    for (index = varused; index < MoLNumSandRVariables-1; index++)
    {
      SandRVariableIndex[index] = SandRVariableIndex[index+1];

#ifdef MOLDEBUG
    printf("The registered evolved variable was SandR."
           " Now index %d is %d (%s).\n",
           index, SandRVariableIndex[index], 
           CCTK_VarName(SandRVariableIndex[index]));    
#endif

    }
    MoLNumSandRVariables--;
  }
  
  return 0;

}


 /*@@
   @routine    MoL_RegisterEvolvedSlow
   @date       Thu May 30 11:36:59 2002
   @author     Ian Hawke
   @desc 
   Given the index of the GF to be evolved and the RHS GF, it stores 
   the indexes for later use together with various error checking.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterEvolvedRealSlow(CCTK_INT EvolvedIndex, CCTK_INT RHSIndexSlow)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_INT /*  ierr, */ index, varused, numtimelevs1, numtimelevs2;
  
#ifdef MOLDEBUG
  printf("Arrived in MoLRegisterEvolvedSlow \n");
  printf("The indexes are %d and %d.\n",EvolvedIndex, RHSIndexSlow);
  printf("These correspond to variables %s and %s.\n",
         CCTK_VarName(EvolvedIndex),CCTK_VarName(RHSIndexSlow));
  printf("The pointer to EvolvedVariableIndex: %p\n",
         EvolvedVariableIndex);
#endif

  if (!(CCTK_GroupTypeFromVarI(EvolvedIndex)==CCTK_GF))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not a GF and so should "
               "not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(EvolvedIndex)); 
  }

  if (!(CCTK_GroupTypeFromVarI(RHSIndexSlow)==CCTK_GF))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The rhs variable %s is not a GF and so should "
               "not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(RHSIndexSlow)); 
  }    

  if (!(CCTK_VarTypeI(EvolvedIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not of type CCTK_REAL and so "
               "should not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(EvolvedIndex)); 
  }

  if (!(CCTK_VarTypeI(RHSIndexSlow)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The rhs variable %s is not of type CCTK_REAL and "
               "so should not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(RHSIndexSlow));
  }

  numtimelevs1 = CCTK_DeclaredTimeLevelsVI(EvolvedIndex);
  numtimelevs2 = CCTK_DeclaredTimeLevelsVI(RHSIndexSlow);
  
  if ( (numtimelevs1 < 0) || (numtimelevs2 < 0) ) 
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for variable index %i", (int)EvolvedIndex); 
    CCTK_ERROR("The index passed does not correspond to a GF.");
  }

  if (numtimelevs1 < 2)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for variable index %i name %s", 
               (int)EvolvedIndex, CCTK_VarName(EvolvedIndex));
    CCTK_ERROR("The GF passed only has one timelevel. "
              "It must have at least two.");
  }

  varused = 0;

  for (index = 0; (index < MoLNumEvolvedVariablesSlow)&&(!varused); index++)
  {
    varused = (EvolvedIndex == EvolvedVariableIndexSlow[index]);
#ifdef MOLDEBUG
    printf("Registering %d. Checking index %d which is %d\n",EvolvedIndex,
           index,EvolvedVariableIndexSlow[index]);
#endif
  }

  if (varused)
  {

    CCTK_VWarn(2,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The GF %s has already been registered "
               "as an evolved variable with RHS GF %s. "
               "The attempt to register with RHS GF %s will be ignored",
               CCTK_VarName(EvolvedIndex),
               CCTK_VarName(RHSVariableIndexSlow[index-1]),
               CCTK_VarName(RHSIndexSlow));

  }
  else
  {

    if (MoLNumEvolvedVariablesSlow+1 > MoLMaxNumRegisteredVariables)
    {
      const char *evolvedvarname = CCTK_VarName(EvolvedIndex);
      char *registeredvars = 
        VarListToString(EvolvedVariableIndex, MoLNumEvolvedVariablesSlow);
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "You tried to register more than %d slow evolved variables. "
                 "Since this should be the total number of Cactus variables, "
                 "this looks like a bug. You tried to register '%s'. The "
                 "registered variables are: %s.",
                 (int)MoLMaxNumRegisteredVariables,
                 evolvedvarname, registeredvars);
      free(registeredvars); /* NOTREACHED */
    }

    EvolvedVariableIndexSlow[MoLNumEvolvedVariablesSlow] = EvolvedIndex;
    RHSVariableIndexSlow[MoLNumEvolvedVariablesSlow] = RHSIndexSlow;
  
    MoLNumEvolvedVariablesSlow++;

#ifdef MOLDEBUG
    printf("The max number is now %d. Just added %d (%s).\n",
           MoLNumEvolvedVariablesSlow, EvolvedIndex,
           CCTK_VarName(EvolvedIndex));
#endif

  }

  varused = -1;

  for (index = 0; (index < MoLNumConstrainedVariables)&&(!(varused+1)); 
       index++)
  {
    if (EvolvedIndex == ConstrainedVariableIndex[index])
    {
      varused = index;
    }
    
  }
  
  if ((varused+1))
  {
    for (index = varused; index < MoLNumConstrainedVariables-1; index++)
    {
      ConstrainedVariableIndex[index] = ConstrainedVariableIndex[index+1];
    }
    MoLNumConstrainedVariables--;
  }
    
  varused = -1;

  for (index = 0; (index < MoLNumSandRVariables)&&(!(varused+1)); index++)
  {
    if (EvolvedIndex == SandRVariableIndex[index])
    {
      varused = index;
    }
    
#ifdef MOLDEBUG
    printf("Checking SandR var %d. Index %d (evolvedindex %d).\n",
           index, SandRVariableIndex[index], EvolvedIndex);    
#endif

  }
  
  if ((varused+1))
  {
    for (index = varused; index < MoLNumSandRVariables-1; index++)
    {
      SandRVariableIndex[index] = SandRVariableIndex[index+1];

#ifdef MOLDEBUG
    printf("The registered evolved variable was SandR."
           " Now index %d is %d (%s).\n",
           index, SandRVariableIndex[index], 
           CCTK_VarName(SandRVariableIndex[index]));    
#endif

    }
    MoLNumSandRVariables--;
  }
  
  return 0;

}




 /*@@
   @routine    MoL_RegisterConstrained
   @date       Thu May 30 12:35:58 2002
   @author     Ian Hawke
   @desc 
   Given the index of the GF, register it as a constrained variable.
   If there's only one timelevel then ignore it as there will be no
   rotation and so MoL doesn't have to do anything.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterConstrainedReal(CCTK_INT ConstrainedIndex)
{

  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT numtimelevs, varused, evolved, index;

  if (!(CCTK_GroupTypeFromVarI(ConstrainedIndex)==CCTK_GF))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"The variable %s is "
               "not a GF and so should not be registered with "
               "MoLRegisterConstrained.", 
               CCTK_VarName(ConstrainedIndex)); 
  }

  if (!(CCTK_VarTypeI(ConstrainedIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"The variable %s is "
               "not of type CCTK_REAL and so should not be "
               "registered with MoLRegisterConstrained.", 
               CCTK_VarName(ConstrainedIndex)); 
  }
  
  numtimelevs = CCTK_DeclaredTimeLevelsVI(ConstrainedIndex);

  if (numtimelevs < 1) {

    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for "
               "constrained variable index %i", (int)ConstrainedIndex);
    CCTK_ERROR("The index passed does not correspond to a GF.");

  }
  else if (numtimelevs > 1) {
    
    varused = 0;
    
    for (evolved = 0; (evolved < MoLNumEvolvedVariables)&&(!varused); 
         evolved++)
    {
      varused = (EvolvedVariableIndex[evolved] == ConstrainedIndex);
    }
    
    for (evolved = 0; (evolved < MoLNumConstrainedVariables)&&(!varused); 
         evolved++)
    {
      varused = (ConstrainedVariableIndex[evolved] == ConstrainedIndex);
    }
    
    if (!varused)
    {

      if (MoLNumConstrainedVariables+1 > MoLMaxNumRegisteredVariables)
      {
        const char *constrainedvarname = CCTK_VarName(ConstrainedIndex);
        char *registeredvars = 
          VarListToString(ConstrainedVariableIndex, MoLNumConstrainedVariables);
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "You tried to register more than %d constrained variables. Since "
                   "this should be the total number of Cactus variables, this "
                   "looks like a bug. You tried to register '%s'. The "
                   "registered variables are: %s.",
                   (int)MoLMaxNumRegisteredVariables, constrainedvarname,
                   registeredvars);
        free(registeredvars); /* NOTREACHED */
      }

      ConstrainedVariableIndex[MoLNumConstrainedVariables] = ConstrainedIndex;
      MoLNumConstrainedVariables++;

    }
    
    varused = -1;
    
    for (evolved = 0; (evolved < MoLNumSandRVariables)&&(!(varused+1)); 
         evolved++)
    {
      if (SandRVariableIndex[evolved] == ConstrainedIndex)
      {
        varused = evolved;
      }
      
    }
    
    if ((varused+1))
    {
      for (index = varused; index < MoLNumSandRVariables-1; index++)
      {
        SandRVariableIndex[index] = SandRVariableIndex[index+1];
      }
      MoLNumSandRVariables--;
    }
    
  }
  else
  {
    
    CCTK_VInfo(CCTK_THORNSTRING,
               "MoL will not treat variable %s as a constrained "
               "variable as it has only one timelevel. This should "
               "not cause problems with the evolution.", 
               CCTK_VarName(ConstrainedIndex));

  }
  
  return 0;
  
}

 /*@@
   @routine    MoL_RegisterSaveAndRestore
   @date       Thu May 30 12:37:40 2002
   @author     Ian Hawke
   @desc 
   Given a GF index store it for later use as a save and restore type.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_RegisterSaveAndRestoreReal(CCTK_INT SandRIndex)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_INT numtimelevs, varused, evolved;

  if (!(CCTK_GroupTypeFromVarI(SandRIndex)==CCTK_GF))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"The variable %s is not "
               "a GF and so should not be registered with "
               "MoLRegisterSaveAndRestore.", 
               CCTK_VarName(SandRIndex)); 
  }

  if (!(CCTK_VarTypeI(SandRIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"The variable %s is "
               "not of type CCTK_REAL and so should not be "
               "registered with MoLRegisterSaveAndRestore.", 
               CCTK_VarName(SandRIndex)); 
  }
  
  numtimelevs = CCTK_DeclaredTimeLevelsVI(SandRIndex);

  if (numtimelevs < 1) {

    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Warning for save "
               "and restore variable index %i", (int)SandRIndex);
    CCTK_ERROR("The index passed does not correspond to a GF.");

  }
  else if (numtimelevs > 1) {
    
    varused = 0;
    
    for (evolved = 0; (evolved < MoLNumEvolvedVariables)&&(!varused); 
         evolved++)
    {
      varused = (EvolvedVariableIndex[evolved] == SandRIndex);
    }
    
    for (evolved = 0; (evolved < MoLNumConstrainedVariables)&&(!varused); 
         evolved++)
    {
      varused = (ConstrainedVariableIndex[evolved] == SandRIndex);
    }
    
    for (evolved = 0; (evolved < MoLNumSandRVariables)&&(!varused); evolved++)
    {
      varused = (SandRVariableIndex[evolved] == SandRIndex);
    }

    if (!varused)
    {

      if (MoLNumSandRVariables+1 > MoLMaxNumRegisteredVariables)
      {
        const char *saveandrestorevarname = CCTK_VarName(SandRIndex);
        char *registeredvars = 
          VarListToString(SandRVariableIndex, MoLNumSandRVariables);
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "You tried to register more than %d safe and restore variables. Since "
                   "this should be the total number of Cactus variables, this "
                   "looks like a bug. You tried to register '%s'. The "
                   "registered variables are: %s.",
                   (int)MoLMaxNumRegisteredVariables, saveandrestorevarname,
                   registeredvars);
        free(registeredvars); /* NOTREACHED */
      }

      SandRVariableIndex[MoLNumSandRVariables] = SandRIndex;
      MoLNumSandRVariables++;

    }
        
  }
  else
  {
    
    CCTK_VInfo(CCTK_THORNSTRING,
               "MoL will not treat variable %s as a save and "
               "restore variable as it has only one timelevel. "
               "This should not cause problems with the evolution.", 
               CCTK_VarName(SandRIndex));

  }
  
  return 0;
  
}

CCTK_INT MoL_RegisterEvolvedRealGroup(CCTK_INT EvolvedGroupIndex, 
                                      CCTK_INT RHSGroupIndex)
{

  CCTK_INT EvolvedGroupFirstVar, RHSGroupFirstVar, GroupNumVars, retval;
  CCTK_INT EvolvedVar, RHSVar;
  
  EvolvedGroupFirstVar = CCTK_FirstVarIndexI(EvolvedGroupIndex);
  if (EvolvedGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved group index %i is not a valid group index.",
               (int)EvolvedGroupIndex);
  }
  
  RHSGroupFirstVar = CCTK_FirstVarIndexI(RHSGroupIndex);
  if (RHSGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "RHS group index %d is not a valid group index.",
               (int)RHSGroupIndex);
  }

  GroupNumVars = CCTK_NumVarsInGroupI(EvolvedGroupIndex);
  if (CCTK_NumVarsInGroupI(RHSGroupIndex) != GroupNumVars)
  {
    char *groupname = CCTK_GroupName(EvolvedGroupIndex);
    char *rhsgroupname = CCTK_GroupName(RHSGroupIndex);
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There are a different number of variables in "
               "evolved group %d ('%s') and RHS group %d ('%s').", 
               (int)EvolvedGroupIndex, groupname, (int)RHSGroupIndex,
               rhsgroupname);
    free(rhsgroupname);
    free(groupname);
  }
  
  retval = 0;
  
  for (EvolvedVar = EvolvedGroupFirstVar, RHSVar = RHSGroupFirstVar;
       EvolvedVar < EvolvedGroupFirstVar + GroupNumVars;
       EvolvedVar++, RHSVar++)
  {
    retval += MoL_RegisterEvolvedReal(EvolvedVar, RHSVar);
  }
  
  return retval;
}


CCTK_INT MoL_RegisterEvolvedRealGroupSlow(CCTK_INT EvolvedGroupIndex, 
                                          CCTK_INT RHSGroupIndexSlow)
{

  CCTK_INT EvolvedGroupFirstVar, RHSGroupFirstVar, GroupNumVars, retval;
  CCTK_INT EvolvedVar, RHSVar;
  
  EvolvedGroupFirstVar = CCTK_FirstVarIndexI(EvolvedGroupIndex);
  if (EvolvedGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved group index %i is not a valid group index.",
               (int)EvolvedGroupIndex);
  }
  
  RHSGroupFirstVar = CCTK_FirstVarIndexI(RHSGroupIndexSlow);
  if (RHSGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "RHS group index %d is not a valid group index.",
               (int)RHSGroupIndexSlow);
  }

  GroupNumVars = CCTK_NumVarsInGroupI(EvolvedGroupIndex);
  if (CCTK_NumVarsInGroupI(RHSGroupIndexSlow) != GroupNumVars)
  {
    char *groupname = CCTK_GroupName(EvolvedGroupIndex);
    char * rhsgroupname = CCTK_GroupName(RHSGroupIndexSlow);
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There are a different number of variables in "
               "evolved group %d ('%s') and RHS group %d ('%s').", 
               (int)EvolvedGroupIndex,
               groupname, (int)RHSGroupIndexSlow, rhsgroupname);
    free(rhsgroupname);
    free(groupname);
  }
  
  retval = 0;
  
  for (EvolvedVar = EvolvedGroupFirstVar, RHSVar = RHSGroupFirstVar;
       EvolvedVar < EvolvedGroupFirstVar + GroupNumVars;
       EvolvedVar++, RHSVar++)
  {
    retval += MoL_RegisterEvolvedRealSlow(EvolvedVar, RHSVar);
  }
  
  return retval;
}


CCTK_INT MoL_RegisterConstrainedRealGroup(CCTK_INT ConstrainedGroupIndex)
{
  
  CCTK_INT ConstrainedGroupFirstVar, GroupNumVars, retval;
  CCTK_INT ConstrainedVar;
  
  ConstrainedGroupFirstVar = CCTK_FirstVarIndexI(ConstrainedGroupIndex);
  if (ConstrainedGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Constrained group index %i is not a valid group index.",
               (int)ConstrainedGroupIndex);
  }
  
  GroupNumVars = CCTK_NumVarsInGroupI(ConstrainedGroupIndex);
  
  retval = 0;
  
  for (ConstrainedVar = ConstrainedGroupFirstVar;
       ConstrainedVar < ConstrainedGroupFirstVar + GroupNumVars;
       ConstrainedVar++)
  {
    retval += MoL_RegisterConstrainedReal(ConstrainedVar);
  }
  
  return retval;
}

CCTK_INT MoL_RegisterSaveAndRestoreRealGroup(CCTK_INT SandRGroupIndex)
{
  
  CCTK_INT SandRGroupFirstVar, GroupNumVars, retval;
  CCTK_INT SandRVar;
  
  SandRGroupFirstVar = CCTK_FirstVarIndexI(SandRGroupIndex);
  if (SandRGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Save and Restore group index %i is not a valid group index.",
               (int)SandRGroupIndex);
  }
  
  GroupNumVars = CCTK_NumVarsInGroupI(SandRGroupIndex);
  
  retval = 0;
  
  for (SandRVar = SandRGroupFirstVar;
       SandRVar < SandRGroupFirstVar + GroupNumVars;
       SandRVar++)
  {
    retval += MoL_RegisterSaveAndRestoreReal(SandRVar);
  }
  
  return retval;
}

CCTK_INT MoL_RegisterEvolvedArray(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_INT /*  ierr,  */ index, varused, numtimelevs1, numtimelevs2;
  
#ifdef MOLDEBUG
  printf("Arrived in MoLRegisterEvolvedArray \n");
  printf("The indexes are %d and %d.\n",EvolvedIndex, RHSIndex);
  printf("These correspond to variables %s and %s.\n",
         CCTK_VarName(EvolvedIndex),CCTK_VarName(RHSIndex));
  printf("The pointer to EvolvedVariableIndex: %p\n",
         EvolvedArrayVariableIndex);
#endif

  if (!(CCTK_GroupTypeFromVarI(EvolvedIndex)==CCTK_ARRAY))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not an array and so should "
               "not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(EvolvedIndex)); 
  }

  if (!(CCTK_GroupTypeFromVarI(RHSIndex)==CCTK_ARRAY))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The rhs variable %s is not an array and so "
               "should not be registered with MoLRegisterEvolved.", 
               CCTK_VarName(RHSIndex)); 
  }    

  if (!(CCTK_VarTypeI(EvolvedIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not of type CCTK_REAL "
               "and so should not be registered with "
               "MoLRegisterEvolvedArray.", 
               CCTK_VarName(EvolvedIndex)); 
  }

  if (!(CCTK_VarTypeI(RHSIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The rhs variable %s is not of type CCTK_REAL "
               "and so should not be registered with "
               "MoLRegisterEvolvedArray.", 
               CCTK_VarName(RHSIndex)); 
  }

  numtimelevs1 = CCTK_DeclaredTimeLevelsVI(EvolvedIndex);
  numtimelevs2 = CCTK_DeclaredTimeLevelsVI(RHSIndex);
  
  if ( (numtimelevs1 < 0) || (numtimelevs2 < 0) ) 
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for variable index %i name %s", 
               (int)EvolvedIndex, CCTK_VarName(EvolvedIndex)); 
    CCTK_ERROR("The index passed does not correspond to an array.");
  }

  if (numtimelevs1 < 2)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for variable index %i name %s", 
               (int)EvolvedIndex, CCTK_VarName(EvolvedIndex));
    CCTK_ERROR("The array passed only has one timelevel. "
              "It must have at least two.");
  }

  varused = 0;

  for (index = 0; (index < MoLNumEvolvedArrayVariables)&&(!varused); index++)
  {
    varused = (EvolvedIndex == EvolvedArrayVariableIndex[index]);
#ifdef MOLDEBUG
    printf("Registering %d. Checking index %d which is %d\n",EvolvedIndex,
           index,EvolvedArrayVariableIndex[index]);
#endif
  }

  if (varused)
  {

    CCTK_VWarn(2,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The array %s has already been registered "
               "as a array evolved variable with RHS array %s. "
               "The attempt to register with RHS array %s will be ignored",
               CCTK_VarName(EvolvedIndex),
               CCTK_VarName(RHSArrayVariableIndex[index-1]),
               CCTK_VarName(RHSIndex));

  }
  else
  {
  
    if (MoLNumEvolvedArrayVariables+1 > MoLMaxNumRegisteredVariables)
    {
      const char *evolvedvarname = CCTK_VarName(EvolvedIndex);
      char *registeredvars = 
        VarListToString(EvolvedArrayVariableIndex, MoLNumEvolvedArrayVariables);
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "You tried to register more than %d evolved array variables. "
                 "Since this should be the total number of Cactus variables, "
                 "this looks like a bug. You tried to register '%s'. The "
                 "registered variables are: %s.",
                 (int)MoLMaxNumRegisteredVariables,
                 evolvedvarname, registeredvars);
      free(registeredvars); /* NOTREACHED */
    }

    EvolvedArrayVariableIndex[MoLNumEvolvedArrayVariables] = EvolvedIndex;
    RHSArrayVariableIndex[MoLNumEvolvedArrayVariables] = RHSIndex;
  
    MoLNumEvolvedArrayVariables++;
    
  }

  varused = 0;

  for (index = 0; (index <
                   MoLNumConstrainedArrayVariables)&&(!varused);
       index++) 
  {
    varused = (EvolvedIndex == ConstrainedArrayVariableIndex[index]);
  }
  
  if (varused)
  {
    for (index = varused; index < MoLNumConstrainedArrayVariables-1; index++)
    {
      ConstrainedArrayVariableIndex[index] =
        ConstrainedArrayVariableIndex[index+1]; 
    }
    MoLNumConstrainedArrayVariables--;
  }
    
  varused = 0;

  for (index = 0; (index < MoLNumSandRArrayVariables)&&(!varused); index++)
  {
    varused = (EvolvedIndex == SandRArrayVariableIndex[index]);
  }
  
  if (varused)
  {
    for (index = varused; index < MoLNumSandRArrayVariables-1; index++)
    {
      SandRArrayVariableIndex[index] = SandRArrayVariableIndex[index+1];
    }
    MoLNumSandRArrayVariables--;
  }

  return 0;

}

CCTK_INT MoL_RegisterConstrainedArray(CCTK_INT ConstrainedIndex)
{

  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT numtimelevs, varused, evolved, index;
    
  if (!(CCTK_GroupTypeFromVarI(ConstrainedIndex)==CCTK_ARRAY))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not a ARRAY and so "
               "should not be registered with MoLRegisterConstrained.", 
               CCTK_VarName(ConstrainedIndex)); 
  }

  if (!(CCTK_VarTypeI(ConstrainedIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not of type CCTK_REAL "
               "and so should not be registered with "
               "MoLRegisterConstrainedArray.", 
               CCTK_VarName(ConstrainedIndex)); 
  }

  numtimelevs = CCTK_DeclaredTimeLevelsVI(ConstrainedIndex);

  if (numtimelevs < 1) {

    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for constrained array variable index %i name %s", 
               (int)ConstrainedIndex, CCTK_VarName(ConstrainedIndex));
    CCTK_ERROR("The index passed does not correspond to a array.");

  }
  else if (numtimelevs > 1) {
    
    varused = 0;
    
    for (evolved = 0; (evolved < MoLNumEvolvedArrayVariables)
           &&(!varused); evolved++)
    {
      varused = (EvolvedArrayVariableIndex[evolved] == ConstrainedIndex);
    }
    
    for (evolved = 0; (evolved <
                       MoLNumConstrainedArrayVariables)&&(!varused);
         evolved++) 
    {
      varused = (ConstrainedArrayVariableIndex[evolved] == ConstrainedIndex);
    }
    
    if (!varused)
    {

      if (MoLNumConstrainedArrayVariables+1 > MoLMaxNumRegisteredVariables)
      {
        const char *constrainedvarname = CCTK_VarName(ConstrainedIndex);
        char *registeredvars = 
          VarListToString(ConstrainedArrayVariableIndex, MoLNumConstrainedArrayVariables);
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "You tried to register more than %d constrained array "
                   "variables. Since " "this should be the total number of "
                   "Cactus variables, this looks like a bug. You tried to "
                   "register '%s'. The " "registered variables are: %s.",
                   (int)MoLMaxNumRegisteredVariables, constrainedvarname,
                   registeredvars);
        free(registeredvars); /* NOTREACHED */
      }

      ConstrainedArrayVariableIndex[MoLNumConstrainedArrayVariables] = 
        ConstrainedIndex;
      MoLNumConstrainedArrayVariables++;

    }
    
    varused = 0;
    
    for (evolved = 0; (evolved <
                       MoLNumSandRArrayVariables)&&(!varused);
         evolved++) 
    {
      varused = (SandRArrayVariableIndex[evolved] == ConstrainedIndex);
    }
    
    if (varused)
    {
      for (index = evolved; index < MoLNumSandRArrayVariables-1; index++)
      {
        SandRArrayVariableIndex[index] = SandRArrayVariableIndex[index+1];
      }
      MoLNumSandRArrayVariables--;
    }
    
  }
  else
  {
    
    CCTK_VInfo(CCTK_THORNSTRING,
               "MoL will not treat variable %s as a constrained "
               "array variable as it has only one timelevel. "
               "This should not cause problems with the evolution.", 
               CCTK_VarName(ConstrainedIndex));

  }  

  return 0;

}

CCTK_INT MoL_RegisterSaveAndRestoreArray(CCTK_INT SandRIndex)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_INT numtimelevs, varused, evolved;
  
  if (!(CCTK_GroupTypeFromVarI(SandRIndex)==CCTK_ARRAY))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not an array and so should not "
               "be registered with MoLRegisterSaveAndRestore.", 
               CCTK_VarName(SandRIndex)); 
  }

  if (!(CCTK_VarTypeI(SandRIndex)==CCTK_VARIABLE_REAL))
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "The variable %s is not of type CCTK_REAL and "
               "so should not be registered with "
               "MoLRegisterSaveAndRestoreArray.", 
               CCTK_VarName(SandRIndex)); 
  }

  numtimelevs = CCTK_DeclaredTimeLevelsVI(SandRIndex);

  if (numtimelevs < 1) {

    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Warning for save and restore array variable index %i", 
               (int)SandRIndex);
    CCTK_ERROR("The index passed does not correspond to a GF.");

  }
  else if (numtimelevs > 1) {
    
    varused = 0;
    
    for (evolved = 0; (evolved < MoLNumEvolvedArrayVariables)
           &&(!varused); evolved++)
    {
      varused = (EvolvedArrayVariableIndex[evolved] == SandRIndex);
    }
    
    for (evolved = 0; (evolved <
                       MoLNumConstrainedArrayVariables)&&(!varused);
         evolved++) 
    {
      varused = (ConstrainedArrayVariableIndex[evolved] == SandRIndex);
    }
    
    for (evolved = 0; (evolved <
                       MoLNumSandRArrayVariables)&&(!varused);
         evolved++) 
    {
      varused = (SandRArrayVariableIndex[evolved] == SandRIndex);
    }

    if (!varused)
    {

      if (MoLNumSandRArrayVariables+1 > MoLMaxNumRegisteredVariables)
      {
        const char *saveandrestorevarname = CCTK_VarName(SandRIndex);
        char *registeredvars = 
          VarListToString(SandRVariableIndex, MoLNumSandRArrayVariables);
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "You tried to register more than %d save and restore array "
                   "variables. Since this should be the total number of "
                   "Cactus variables, this looks like a bug. You tried to "
                   "register '%s'. The registered variables are: %s.",
                   (int)MoLMaxNumRegisteredVariables, saveandrestorevarname,
                   registeredvars);
        free(registeredvars); /* NOTREACHED */
      }

      SandRArrayVariableIndex[MoLNumSandRArrayVariables] = SandRIndex;
      MoLNumSandRArrayVariables++;

    }
        
  }
  else
  {
    
    CCTK_VInfo(CCTK_THORNSTRING,
               "MoL will not treat variable %s as a save "
               "and restore array variable as it has only one "
               "timelevel. This should not cause problems "
               "with the evolution.", 
               CCTK_VarName(SandRIndex));

  }  

  return 0;

}

CCTK_INT MoL_RegisterEvolvedArrayGroup(CCTK_INT EvolvedGroupIndex, 
                                       CCTK_INT RHSGroupIndex)
{

  CCTK_INT EvolvedGroupFirstVar, RHSGroupFirstVar, GroupNumVars, retval;
  CCTK_INT EvolvedVar, RHSVar;
  
  EvolvedGroupFirstVar = CCTK_FirstVarIndexI(EvolvedGroupIndex);
  if (EvolvedGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved array group index %i is not a valid group index.",
               (int)EvolvedGroupIndex);
  }
  
  RHSGroupFirstVar = CCTK_FirstVarIndexI(RHSGroupIndex);
  if (RHSGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "RHS array group index %d is not a valid group index.",
               (int)RHSGroupIndex);
  }

  GroupNumVars = CCTK_NumVarsInGroupI(EvolvedGroupIndex);
  if (CCTK_NumVarsInGroupI(RHSGroupIndex) != GroupNumVars)
  {
    char *groupname = CCTK_GroupName(EvolvedGroupIndex);
    char *rhsgroupname = CCTK_GroupName(RHSGroupIndex);
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There are a different number of variables in "
               "evolved array group %d ('%s') and RHS group %d ('%s').", 
               (int)EvolvedGroupIndex, groupname, (int)RHSGroupIndex,
               rhsgroupname);
    free(rhsgroupname);
    free(groupname);
  }
  
  retval = 0;
  
  for (EvolvedVar = EvolvedGroupFirstVar, RHSVar = RHSGroupFirstVar;
       EvolvedVar < EvolvedGroupFirstVar + GroupNumVars;
       EvolvedVar++, RHSVar++)
  {
    retval += MoL_RegisterEvolvedArray(EvolvedVar, RHSVar);
  }

  return retval;

}

CCTK_INT MoL_RegisterConstrainedArrayGroup(CCTK_INT ConstrainedGroupIndex)
{
  
  CCTK_INT ConstrainedGroupFirstVar, GroupNumVars, retval;
  CCTK_INT ConstrainedVar;
  
  ConstrainedGroupFirstVar = CCTK_FirstVarIndexI(ConstrainedGroupIndex);
  if (ConstrainedGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Constrained array group index %i is not a valid group index.",
               (int)ConstrainedGroupIndex);
  }
  
  GroupNumVars = CCTK_NumVarsInGroupI(ConstrainedGroupIndex);
  
  retval = 0;
  
  for (ConstrainedVar = ConstrainedGroupFirstVar;
       ConstrainedVar < ConstrainedGroupFirstVar + GroupNumVars;
       ConstrainedVar++)
  {
    retval += MoL_RegisterConstrainedArray(ConstrainedVar);
  }
  
  return retval;

}

CCTK_INT MoL_RegisterSaveAndRestoreArrayGroup(CCTK_INT SandRGroupIndex)
{
  
  CCTK_INT SandRGroupFirstVar, GroupNumVars, retval;
  CCTK_INT SandRVar;
  
  SandRGroupFirstVar = CCTK_FirstVarIndexI(SandRGroupIndex);
  if (SandRGroupFirstVar < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Save and Restore array group index %i is "
               "not a valid group index.",
               (int)SandRGroupIndex);
  }
  
  GroupNumVars = CCTK_NumVarsInGroupI(SandRGroupIndex);
  
  retval = 0;
  
  for (SandRVar = SandRGroupFirstVar;
       SandRVar < SandRGroupFirstVar + GroupNumVars;
       SandRVar++)
  {
    retval += MoL_RegisterSaveAndRestoreArray(SandRVar);
  }
  
  return retval;

}

 /*@@
   @routine    MoL_QueryEvolvedRHS
   @date       Thu Dec 25 2004
   @author     Erik Schnetter
   @desc 
   Given the index of a registered evolved GF,
   returns the corresponding RHS GF.
   Checks both regular and slow evolved variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_QueryEvolvedRHS(CCTK_INT EvolvedIndex)
{
  CCTK_INT index;
  
  if (EvolvedIndex < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Evolved variables index %i is not a valid variable index.",
               (int)EvolvedIndex);
  }
  
  for (index = 0; index < MoLNumEvolvedVariables; index++)
  {
    if (EvolvedVariableIndex[index] == EvolvedIndex)
    {
      return RHSVariableIndex[index];
    }
  }
  
  // If we haven't found any registered variable that matches,
  // we also try to check whether this has been registered with the slow multirate sector!
  
  for (index = 0; index < MoLNumEvolvedVariablesSlow; index++)
  {
    if (EvolvedVariableIndexSlow[index] == EvolvedIndex)
    {
      return RHSVariableIndexSlow[index];
    }
  }
  
  return -1;

}



 /*@@
   @routine    MoL_NumIntegratorSubsteps
   @date       2007-05-23
   @author     Erik Schnetter
   @desc 
   Return how many time integrator substeps there are.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT MoL_NumIntegratorSubsteps (void)
{
  DECLARE_CCTK_PARAMETERS;

  return MoL_Intermediate_Steps;
}



/*
  Old function names. Just calls the new version. 
  Included for compatibility
  Now removed.
*/

/*
CCTK_INT MoL_RegisterVar(CCTK_INT molvarindex,CCTK_INT molrhsvarindex)
{
  return MoL_RegisterEvolved(molvarindex, molrhsvarindex);
}

CCTK_INT MoL_RegisterPrimitive(CCTK_INT primitiveindex)
{
  return MoL_RegisterConstrained(primitiveindex);
}

CCTK_INT MoL_RegisterDepends(CCTK_INT dependsindex)
{
  return MoL_RegisterSaveAndRestore(dependsindex);
}

CCTK_INT MoL_RegisterVarGroup(CCTK_INT groupindex,CCTK_INT rhsgroupindex)
{
  return MoL_RegisterEvolvedGroup(groupindex, rhsgroupindex);
}

CCTK_INT MoL_RegisterPrimitiveGroup(CCTK_INT groupindex)
{
  return MoL_RegisterConstrainedGroup(groupindex);
}

CCTK_INT MoL_RegisterDependsGroup(CCTK_INT groupindex)
{
  return MoL_RegisterConstrainedGroup(groupindex);
}

CCTK_INT MoL_ChangeVarToEvolved(CCTK_INT varindex, CCTK_INT rhsindex)
{
  return MoL_ChangeToEvolved(varindex, rhsindex);
}

CCTK_INT MoL_ChangeVarToDependent(CCTK_INT dependsindex)
{
  return MoL_ChangeToSaveAndRestore(dependsindex);
}

CCTK_INT MoL_ChangeVarToPrimitive(CCTK_INT primitiveindex)
{
  return MoL_ChangeToConstrained(primitiveindex);
}

CCTK_INT MoL_ChangeVarToNone(CCTK_INT removeindex)
{
  return MoL_ChangeToNone(removeindex);
}
*/

/* 
  Fortran wrappers for the above functions. 
  Should be replaced by using function aliasing eventually.
  Now removed.
*/

/*
void CCTK_FCALL CCTK_FNAME(MoL_RegisterEvolved)(int *ierr, 
                                            CCTK_INT *EvolvedIndex,
                                            CCTK_INT *RHSIndex)
{
  *ierr = MoL_RegisterEvolved(*EvolvedIndex, *RHSIndex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterConstrained)(int *ierr, 
                                            CCTK_INT *EvolvedIndex)
{
  *ierr = MoL_RegisterConstrained(*EvolvedIndex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterSaveAndRestore)(int *ierr, 
                                            CCTK_INT *EvolvedIndex)
{
  *ierr = MoL_RegisterSaveAndRestore(*EvolvedIndex);
  return;
}
*/

/* And for MoL compatibility... */
/* Now removed */

/*
void CCTK_FCALL CCTK_FNAME(MoL_RegisterVar)(int *ierr, 
                                            CCTK_INT *molvarindex,
                                            CCTK_INT *molrhsvarindex)
{
  *ierr = MoL_RegisterVar(*molvarindex, *molrhsvarindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterDepends)(int *ierr, 
                                                CCTK_INT *moldependsindex)
{
  *ierr = MoL_RegisterDepends(*moldependsindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterPrimitive)(int *ierr, 
                                                CCTK_INT *molprimitiveindex)
{
  *ierr = MoL_RegisterPrimitive(*molprimitiveindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterVarGroup)(int *ierr, 
                                                 CCTK_INT *groupindex,
                                                 CCTK_INT *rhsgroupindex)
{
  *ierr = MoL_RegisterVarGroup(*groupindex, *rhsgroupindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterPrimitiveGroup)(int *ierr, 
                                                       CCTK_INT *groupindex)
{
  *ierr = MoL_RegisterPrimitiveGroup(*groupindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_RegisterDependsGroup)(int *ierr, 
                                                       CCTK_INT *groupindex)
{
  *ierr = MoL_RegisterDependsGroup(*groupindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_ChangeVarToEvolved)(int *ierr, 
                                                   CCTK_INT *varindex,
                                                   CCTK_INT *rhsindex)
{
  *ierr = MoL_ChangeVarToEvolved(*varindex, *rhsindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_ChangeVarToDependent)(int *ierr, 
                                                     CCTK_INT *dependsindex)
{
  *ierr = MoL_ChangeVarToDependent(*dependsindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_ChangeVarToPrimitive)(int *ierr, 
                                                     CCTK_INT *primitiveindex)
{
  *ierr = MoL_ChangeVarToPrimitive(*primitiveindex);
  return;
}

void CCTK_FCALL CCTK_FNAME(MoL_ChangeVarToNone)(int *ierr, 
                                                CCTK_INT *removeindex)
{
  *ierr = MoL_ChangeVarToNone(*removeindex);
  return;
}
*/

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

static char *VarListToString(const CCTK_INT varindices[], const int nvars)
{
  size_t off = 0, len = 1;
  char *vars = NULL;
  /* count how much room we need for the full output string */
  for(int index = 0 ; index < nvars ; index++)
  {
    char *varname = CCTK_FullName(varindices[index]);
    len += strlen(varname) + 1; /* ' ' or '\0' as separator */
    free(varname);
  }
  vars = malloc(len);
  assert(vars != NULL);
  /* construct list in buffer */
  for(int index = 0 ; index < nvars ; index++)
  {
    char *varname = CCTK_FullName(varindices[index]);
    const size_t written =
      snprintf(vars+off, len-off, "%s%s", off>0?" ":"", varname);
    assert(written < len-off);
    off = strlen(vars);
    free(varname);
  }
  vars[off] = '\0'; /* close string in case we had no variables at all */
  return vars;
}
