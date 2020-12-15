 /*@@
   @file      ScratchSpace.c
   @date      Tue Jul 10 11:26:48 PDT 2012
   @author    Roland Haas
   @desc 
   Allocates and deallocates memory for scratch levels.
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

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_ScratchSpace_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static void MoL_AllocateScratchLevelsForVar(const cGH *cctkGH, const char *scratchgroupname, int timelevels);

static void MoL_FreeScratchLevelsForVar(const cGH *cctkGH, const char *scratchgroupname);

/* support old versions of the flesh */
#ifndef HAVE_CCTK_DECLARED_TIMELEVELS
#define CCTK_DeclaredTimeLevelsGI(gi) CCTK_MaxTimeLevelsGI(gi)
#endif

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

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
   @routine    MoL_AllocateScratchSpace
   @date       Tue Jul 10 11:29:00 PDT 2012
   @author     Roland Haas
   @desc 
   Allocates scratch space for a all registered variables.
   @enddesc 
   @calls      MoL_AllocateScratchLevelsForVar
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_AllocateScratchSpace(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  const int need_error_estimate = CCTK_Equals(ODE_Method,"RK45") || 
    CCTK_Equals(ODE_Method,"RK45CK") || CCTK_Equals(ODE_Method,"RK65") || CCTK_Equals(ODE_Method,"RK87");

  /* we abuse the Cactus timelevels (since they are something we can change at
   * run time) to obtain mesh refined scratch space */

  /* CCTK_REAL GF */
  MoL_AllocateScratchLevelsForVar(cctkGH, "MoL::ScratchSpace", MoLNumEvolvedVariables);
  MoL_AllocateScratchLevelsForVar(cctkGH, "MoL::ScratchSpaceSlow", MoLNumEvolvedVariablesSlow);
  MoL_AllocateScratchLevelsForVar(cctkGH, "MoL::SandRScratchSpace", MoLNumSandRVariables);
  if (need_error_estimate) /* adaptive stepping is only supported for CCTK_REAL grid functions */
    MoL_AllocateScratchLevelsForVar(cctkGH, "MoL::ErrorEstimate", MoLNumEvolvedVariables);

  /* CCTK_REAL GA */
  /*MoL_AllocateScratchLevelsForVar(cctkGH, "MoL::ArrayScratchSpace", MoLNumEvolvedArrayVariables);*/
  /*MoL_AllocateScratchLevelsForVar(cctkGH, "MoL::ArraySandRScratchSpace", MoLNumSandRArrayVariables);*/

}

 /*@@
   @routine    MoL_FreeScratchSpace
   @date       Tue Jul 10 11:29:00 PDT 2012
   @author     Roland Haas
   @desc 
   Free all scratch space.
   @enddesc 
   @calls      MoL_FreeScratchLevelsForVar
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_FreeScratchSpace(CCTK_ARGUMENTS)
{
  /* CCTK_REAL GF */
  MoL_FreeScratchLevelsForVar(cctkGH, "MoL::ScratchSpace");
  MoL_FreeScratchLevelsForVar(cctkGH, "MoL::ScratchSpaceSlow");
  MoL_FreeScratchLevelsForVar(cctkGH, "MoL::SandRScratchSpace");
  MoL_FreeScratchLevelsForVar(cctkGH, "MoL::ErrorEstimate");

  /* CCTK_REAL GA */
  /*MoL_FreeScratchLevelsForVar(cctkGH, "MoL::ArrayScratchSpace");*/
  /*MoL_FreeScratchLevelsForVar(cctkGH, "MoL::ArraySandRScratchSpace");*/
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    MoL_AllocateScratchLevelsForVar
   @date       Tue Jul 10 11:29:00 PDT 2012
   @author     Roland Haas
   @desc 
   Allocates scratch space for a single variable by manipulating the number of
   timelevels present.
   @enddesc 
   @calls      CCTK_GroupStorageIncrease
   @calledby   
   @history 
 
   @endhistory 

@@*/

static void MoL_AllocateScratchLevelsForVar(const cGH *cctkGH, const char *scratchgroupname, int timelevels)
{
  int scratchgroup, ierr, activetimelevels;
  cGroup group;

  scratchgroup = CCTK_GroupIndex(scratchgroupname);
  if (scratchgroup < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to obtain group index for group %s, ierr = %d.", 
               scratchgroupname, scratchgroup);
  }
  /* ODE methods without scratch levels lead to SCRATCHSPACE having 0 variables */
  CCTK_GroupData(scratchgroup, &group);
  if (group.numvars > 0)
  {
    int maxtimelevels = CCTK_DeclaredTimeLevelsGI(scratchgroup);
    if (timelevels > maxtimelevels)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Too many (%d) scratch variables required. Only %d are supported. You will have to change interface.ccl",
                 timelevels, (int)maxtimelevels);
    }

    ierr = CCTK_GroupStorageIncrease(cctkGH, 1, &scratchgroup, &timelevels, NULL);
    ierr += CCTK_GroupStorageDecrease(cctkGH, 1, &scratchgroup, &timelevels, NULL);
    activetimelevels = CCTK_ActiveTimeLevelsGI(cctkGH, scratchgroup);
    if (ierr < 0 || timelevels != activetimelevels)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Failed to adjust space for %s, ierr = %d, %d variables of %d active.", 
                 scratchgroupname, ierr, activetimelevels, timelevels);
    }
  }
}

 /*@@
   @routine    MoL_FreeScratchLevelsForVar
   @date       Tue Jul 10 11:29:00 PDT 2012
   @author     Roland Haas
   @desc 
   Deallocates scratch space for a single variable by manipulating the number
   of timelevels present.
   @enddesc 
   @calls      CCTK_GroupStorageDecrease
   @calledby   
   @history 
 
   @endhistory 

@@*/

static void MoL_FreeScratchLevelsForVar(const cGH *cctkGH, const char *scratchgroupname)
{
  int scratchgroup, ierr, activetimelevels;
  int timelevels = 0;

  scratchgroup = CCTK_GroupIndex(scratchgroupname);
  if (scratchgroup < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to obtain group index for group %s, ierr = %d.", 
               scratchgroupname, scratchgroup);
  }
  ierr = CCTK_GroupStorageDecrease(cctkGH, 1, &scratchgroup, &timelevels, NULL);
  activetimelevels = CCTK_ActiveTimeLevelsGI(cctkGH, scratchgroup);
  if (ierr < 0 || timelevels != activetimelevels)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to deallocate space for %s, ierr = %d, %d variables of %d still active.", 
               scratchgroupname, ierr, activetimelevels, timelevels);
  }
}

