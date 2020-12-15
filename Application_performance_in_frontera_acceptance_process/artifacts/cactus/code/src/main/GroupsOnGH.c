 /*@@
   @file      GroupsOnGH.c
   @date      Tues April 6
   @author    Gabrielle Allen
   @desc
              GH specific Routines to deal with groups.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk_Comm.h"
#include "cctk_Constants.h"
#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_Groups.h"
#include "cctk_GroupsOnGH.h"
#include "cctk_Misc.h"
#include "cctk_WarnLevel.h"

#include "cctki_GroupsOnGH.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_GroupsOnGH_c);


void CCTK_FCALL CCTK_FNAME(CCTK_VarDataPtr)
     (CCTK_POINTER *res,
      const cGH **cctkGH,
      const int *timelevel,
      ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(CCTK_VarDataPtrI)
     (CCTK_POINTER *res,
      const cGH **cctkGH,
      const int *timelevel,
      const int *vindex);

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupashGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ash,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupashGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ash,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupashVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ash,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupashVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ash,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            const int *groupindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            const int *varindex);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsVI)
                           (int *num, const cGH **cctkGH, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsVN)
                           (int *num, const cGH **cctkGH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsGI)
                           (int *num, const cGH **cctkGH, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsGN)
                           (int *num, const cGH **cctkGH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevels)
                           (int *num, const cGH **cctkGH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsVI)
                           (int *num, const cGH **cctkGH, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsVN)
                           (int *num, const cGH **cctkGH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsGI)
                           (int *num, const cGH **cctkGH, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsGN)
                           (int *num, const cGH **cctkGH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevels)
                           (int *num, const cGH **cctkGH, ONE_FORTSTRING_ARG);



 /*@@
   @routine    CCTK_VarDataPtr
   @date       Tue 6th April 1999
   @author     Gabrielle Allen
   @desc
   Passes back a variable data pointer, given a full name and timelevel
   @enddesc

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar

   @var        varname
   @vdesc      Full name of the grid variable
   @vtype      const char *
   @vio        in
   @vcomment   Format <implementation>::<variable>
   @endvar

   @var        timelevel
   @vdesc      Index of timelevel on which data is required
   @vtype      int
   @vio        in
   @endvar

   @returntype void *
   @returndesc Pointer to the required data, should be cast to required type
   @endreturndesc
@@*/

void *CCTK_VarDataPtr(const cGH *GH, int timelevel, const char *varname)
{
  int vindex;
  void *retval;


  retval = NULL;
  vindex = CCTK_VarIndex(varname);
  if (vindex >= 0)
  {
    if (timelevel >= 0 && timelevel < CCTK_MaxActiveTimeLevelsVI (GH, vindex))
    {
      retval = GH->data[vindex][timelevel];
    }
    else
    {
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "Invalid timelevel %d for variable '%s' in CCTK_VarDataPtr",
                 timelevel, varname);
    }
  }
  else
  {
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "Invalid variable name '%s' in CCTK_VarDataPtr", varname);
  }

#ifdef DEBUG_GROUPS
  printf("In CCTK_VarDataPtr\n----------------------------\n");
  printf("  Data pointer for %s (%d) is %x\n",varname,vindex,retval);
#endif

  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_VarDataPtr)
     (CCTK_POINTER *res,
      const cGH **cctkGH,
      const int *timelevel,
      ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname);
  *res = CCTK_VarDataPtr (*cctkGH, *timelevel, varname);
  free (varname);
}

 /*@@
   @routine    CCTK_VarDataPtrI
   @date       Tue 6th April 1999
   @author     Gabrielle Allen
   @desc
   Passes back a variable data pointer, given a variable index and timelevel
   @enddesc
   @history
      A check for a valid index (i>=0) included by Gerd Lanfermann
   @endhistory

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar

   @var        vindex
   @vdesc      Index of grid variable
   @vtype      int
   @vio        in
   @endvar

   @var        timelevel
   @vdesc      Index of timelevel on which data is required
   @vtype      int
   @vio        in
   @endvar

   @returntype void *
   @returndesc Pointer to the required data, should be cast to required type
   @endreturndesc
@@*/
void *CCTK_VarDataPtrI(const cGH *GH, int timelevel, int vindex)
{
  int numtimelevels;
  void *retval;


  retval = NULL;
  numtimelevels = CCTK_MaxActiveTimeLevelsVI (GH, vindex);
  if (numtimelevels > 0)
  {
    if (timelevel >= 0 && timelevel < numtimelevels)
    {
      retval = GH->data[vindex][timelevel];
    }
    else
    {
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_VarDataPtrI: Invalid timelevel %d for variable '%s'",
                 timelevel, CCTK_VarName (vindex));
    }
  }
  else
  {
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_VarDataPtrI: Invalid index %d given", vindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_VarDataPtrI)
     (CCTK_POINTER *res,
      const cGH **cctkGH,
      const int *timelevel,
      const int *vindex)
{
  *res = CCTK_VarDataPtrI (*cctkGH, *timelevel, *vindex);
}

 /*@@
   @routine    CCTKi_VarDataPtrI
   @date       2012-10-26
   @author     Erik Schnetter
   @desc
   Passes back a variable data pointer, given a variable index and timelevel.
   This function does not output any warnings.
   @enddesc

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar

   @var        vindex
   @vdesc      Index of grid variable
   @vtype      int
   @vio        in
   @endvar

   @var        timelevel
   @vdesc      Index of timelevel on which data is required
   @vtype      int
   @vio        in
   @endvar

   @returntype void *
   @returndesc Pointer to the required data, should be cast to required type
   @endreturndesc
@@*/
void *CCTKi_VarDataPtrI(const cGH *GH, int timelevel, int vindex)
{
  int numvars = CCTK_NumVars();
  if (vindex < 0 || vindex >= numvars)
  {
    return NULL;
  }

  int numtimelevels = CCTK_MaxActiveTimeLevelsVI (GH, vindex);
  if (timelevel < 0 || timelevel >= numtimelevels)
  {
    return NULL;
  }
  return GH->data[vindex][timelevel];
}

 /*@@
   @routine    CCTK_VarDataPtrB
   @date       Tue 6th April 1999
   @author     Gabrielle Allen
   @desc
   Passes back a variable data pointer, given either a variable index
   or a full name and timelevel
   @enddesc

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar

   @var        vari
   @vdesc      Index of grid variable
   @vtype      int
   @vio        in
   @vcomment   Assumed to be in correct range
   @endvar

   @var        varn
   @vdesc      Full name of the grid variable
   @vtype      char *
   @vio        in
   @vcomment   Format <implementation>::<variable>
   @endvar

   @var        timelevel
   @vdesc      Index of timelevel on which data is required
   @vtype      int
   @vio        in
   @endvar

   @returntype void *
   @returndesc Pointer to the required data, should be cast to required type
   @endreturndesc
@@*/
void *CCTK_VarDataPtrB(const cGH *GH, int timelevel, int vindex, char *varname)
{
  void *result;


  result = varname ? CCTK_VarDataPtr(GH, timelevel, varname)
                   : CCTK_VarDataPtrI(GH, timelevel, vindex);

  return (result);
}

 /*@@
   @routine    CCTK_EnableGroupCommI
   @date       Sat Feb 13 17:06:30 1999
   @author     Tom Goodale
   @desc
   Enables communication for a group based upon its name.
   @enddesc
@@*/
int CCTK_EnableGroupCommI(const cGH *GH, int group)
{
  int retcode;
  char *group_name;

  retcode = 0;
  group_name = CCTK_GroupName(group);
  if(group_name)
  {
    retcode = CCTK_EnableGroupComm(GH, group_name);
    free(group_name);
  }

  return retcode;
}

 /*@@
   @routine    CCTK_EnableGroupStorageI
   @date       Sat Feb 13 17:06:30 1999
   @author     Tom Goodale
   @desc
   Enables storage for a group based upon its name.
   @enddesc
@@*/
int CCTK_EnableGroupStorageI(const cGH *GH, int group)
{
  int retcode;
  char *group_name;

  retcode = 0;
  group_name = CCTK_GroupName(group);
  if(group_name)
  {
#ifdef DEBUG_GROUPS
    printf("Turning on storage in %s for group %s (%d)\n",__FILE__,group_name,group);
#endif
    retcode = CCTK_EnableGroupStorage(GH, group_name);

    free(group_name);
  }

  return retcode;
}

 /*@@
   @routine    CCTK_DisableGroupCommI
   @date       Sat Feb 13 17:06:30 1999
   @author     Tom Goodale
   @desc
   Routine to switch communication off for a group based upon its index
   @enddesc
@@*/
int CCTK_DisableGroupCommI(const cGH *GH, int group)
{
  int retcode;
  char *group_name;

  retcode = 0;
  group_name = CCTK_GroupName(group);
  if(group_name)
  {
#ifdef DEBUG_GROUPS
    printf("Turning off comm in %s for group %s (%d)\n",__FILE__,group_name,group);
#endif
    retcode = CCTK_DisableGroupComm(GH, group_name);

    free(group_name);
  }

  return retcode;
}

 /*@@
   @routine    CCTK_DisableGroupStorageI
   @date       Sat Feb 13 17:06:30 1999
   @author     Tom Goodale
   @desc
   Routine to switch storage off for a group based upon its index
   @enddesc
@@*/
int CCTK_DisableGroupStorageI(const cGH *GH, int group)
{
  int retcode;
  char *group_name;

  retcode = 0;
  group_name = CCTK_GroupName(group);
  if(group_name)
  {
    retcode = CCTK_DisableGroupStorage(GH, group_name);

    free(group_name);
  }

  return retcode;
}


const int *CCTK_ArrayGroupSizeI(const cGH *GH, int dir, int groupi)
{
  return CCTK_ArrayGroupSizeB(GH,dir,groupi,NULL);
}

const int *CCTK_ArrayGroupSize(const cGH *GH, int dir, const char *groupn)
{
  return CCTK_ArrayGroupSizeB(GH,dir,-1,groupn);
}


/********************************************************************
 ********************    ActiveTimeLevels    ************************
 ********************************************************************/

 /*@@
   @routine    CCTK_ActiveTimeLevels
   @date       Thu July 17 2003
   @author     Gabrielle Allen
   @desc
   Return number of active timelevels for a group
   @enddesc
@@*/
int CCTK_ActiveTimeLevels(const cGH *GH, const char *groupname)
{
  int gindex;
  int timelevels;
  int increase=0;

  gindex = CCTK_GroupIndex(groupname);
  if (gindex < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_ActiveTimeLevels: invalid group name '%s'", groupname);
    return (-1);
  }

  timelevels = CCTK_GroupStorageIncrease (GH, 1, &gindex, &increase, NULL);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevels)
                           (int *timelevels,
                            const cGH **cctkGH,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *timelevels = CCTK_ActiveTimeLevels (*cctkGH, groupname);
  free (groupname);
}

 /*@@
   @routine    CCTK_ActiveTimeLevelsGN
   @date       Thu July 17 2003
   @author     Gabrielle Allen
   @desc
   Return number of active timelevels for a group
   @enddesc
@@*/
int CCTK_ActiveTimeLevelsGN(const cGH *GH, const char *groupname)
{
  int timelevels;

  timelevels = CCTK_ActiveTimeLevels(GH, groupname);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsGN)
                           (int *timelevels,
                            const cGH **cctkGH,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *timelevels = CCTK_ActiveTimeLevels (*cctkGH, groupname);
  free (groupname);
}

 /*@@
   @routine    CCTK_ActiveTimeLevelsGI
   @date       Thu July 17 2003
   @author     Gabrielle Allen
   @desc
   Return number of active timelevels for a group
   @enddesc
@@*/
int CCTK_ActiveTimeLevelsGI(const cGH *GH, int gindex)
{
  int increase=0;
  int timelevels;

  if (gindex < 0 || gindex >= CCTK_NumGroups ())
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_ActiveTimeLevelsGI: invalid group index %d given",
                gindex);
    return (-1);
  }

  timelevels = CCTK_GroupStorageIncrease (GH, 1, &gindex, &increase, NULL);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsGI)
                           (int *timelevels,
                            const cGH **cctkGH,
                            const int *gindex)
{
  *timelevels = CCTK_ActiveTimeLevelsGI (*cctkGH, *gindex);
}


 /*@@
   @routine    CCTK_ActiveTimeLevelsVN
   @date       Thu July 17 2003
   @author     Gabrielle Allen
   @desc
   Return number of active timelevels for a group
   @enddesc
@@*/
int CCTK_ActiveTimeLevelsVN(const cGH *GH, const char *varname)
{
  int timelevels;
  int gindex;

  gindex = CCTK_GroupIndexFromVar(varname);
  if (gindex < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_ActiveTimeLevelsVN: invalid variable name '%s'", varname);
    return (-1);
  }

  timelevels = CCTK_ActiveTimeLevelsGI(GH, gindex);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsVN)
                           (int *timelevels,
                            const cGH **cctkGH,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *timelevels = CCTK_ActiveTimeLevelsVN (*cctkGH, varname);
  free (varname);
}


 /*@@
   @routine    CCTK_ActiveTimeLevelsVI
   @date       Thu July 17 2003
   @author     Gabrielle Allen
   @desc
   Return number of active timelevels for a group
   @enddesc
@@*/
int CCTK_ActiveTimeLevelsVI(const cGH *GH, int vindex)
{
  int timelevels;
  int gindex;

  gindex = CCTK_GroupIndexFromVarI(vindex);
  if (gindex < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_ActiveTimeLevelsGI: invalid variable index %d", vindex);
    return (-1);
  }

  timelevels = CCTK_ActiveTimeLevelsGI(GH, gindex);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_ActiveTimeLevelsVI)
                           (int *timelevels,
                            const cGH **cctkGH,
                            const int *vindex)
{
  *timelevels = CCTK_ActiveTimeLevelsVI (*cctkGH, *vindex);
}


/********************************************************************
 ********************    MaxActiveTimeLevels    *********************
 ********************************************************************/


 /*@@
   @routine    CCTK_MaxActiveTimeLevels
   @date       Mon Mar 10 2014
   @author     Roland Haas
   @desc
   Return maximum number of active timelevels for a group
   @enddesc
@@*/
int CCTK_MaxActiveTimeLevels(const cGH *GH, const char *groupname)
{
  int gindex;
  int status;

  gindex = CCTK_GroupIndex(groupname);
  if (gindex < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_MaxActiveTimeLevels: invalid group name '%s'", groupname);
    return (-1);
  }

  CCTK_QueryMaxTimeLevels(GH, 1, &gindex, &status);

  return status;
}
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevels)
                           (int *timelevels,
                            const cGH **cctkGH,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *timelevels = CCTK_MaxActiveTimeLevels (*cctkGH, groupname);
  free (groupname);
}

 /*@@
   @routine    CCTK_MaxActiveTimeLevelsGN
   @date       Mon Mar 10 2014
   @author     Roland Haas
   @desc
   Return maximum number of active timelevels for a group
   @enddesc
@@*/
int CCTK_MaxActiveTimeLevelsGN(const cGH *GH, const char *groupname)
{
  int timelevels;

  timelevels = CCTK_MaxActiveTimeLevels(GH, groupname);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsGN)
                           (int *timelevels,
                            const cGH **cctkGH,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *timelevels = CCTK_MaxActiveTimeLevels (*cctkGH, groupname);
  free (groupname);
}

 /*@@
   @routine    CCTK_MaxActiveTimeLevelsGI
   @date       Mon Mar 10 2014
   @author     Roland Haas
   @desc
   Return maximum number of active timelevels for a group
   @enddesc
@@*/
int CCTK_MaxActiveTimeLevelsGI(const cGH *GH, int gindex)
{
  int status;

  if (gindex < 0 || gindex >= CCTK_NumGroups ())
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_MaxActiveTimeLevelsGI: invalid group index %d given",
                gindex);
    return (-1);
  }

  CCTK_QueryMaxTimeLevels(GH, 1, &gindex, &status);

  return status;
}
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsGI)
                           (int *timelevels,
                            const cGH **cctkGH,
                            const int *gindex)
{
  *timelevels = CCTK_MaxActiveTimeLevelsGI (*cctkGH, *gindex);
}


 /*@@
   @routine    CCTK_MaxActiveTimeLevelsVN
   @date       Mon Mar 10 2014
   @author     Roland Haas
   @desc
   Return maximum number of active timelevels for a group
   @enddesc
@@*/
int CCTK_MaxActiveTimeLevelsVN(const cGH *GH, const char *varname)
{
  int timelevels;
  int gindex;

  gindex = CCTK_GroupIndexFromVar(varname);
  if (gindex < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_MaxActiveTimeLevelsVN: invalid variable name '%s'", varname);
    return (-1);
  }

  timelevels = CCTK_MaxActiveTimeLevelsGI(GH, gindex);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsVN)
                           (int *timelevels,
                            const cGH **cctkGH,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *timelevels = CCTK_MaxActiveTimeLevelsVN (*cctkGH, varname);
  free (varname);
}


 /*@@
   @routine    CCTK_MaxActiveTimeLevelsVI
   @date       Mon Mar 10 2014
   @author     Roland Haas
   @desc
   Return maximum number of active timelevels for a group
   @enddesc
@@*/
int CCTK_MaxActiveTimeLevelsVI(const cGH *GH, int vindex)
{
  int timelevels;
  int gindex;

  gindex = CCTK_GroupIndexFromVarI(vindex);
  if (gindex < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_MaxActiveTimeLevelsVI: invalid variable index %d", vindex);
    return (-1);
  }

  timelevels = CCTK_MaxActiveTimeLevelsGI(GH, gindex);

  return timelevels;
}
void CCTK_FCALL CCTK_FNAME (CCTK_MaxActiveTimeLevelsVI)
                           (int *timelevels,
                            const cGH **cctkGH,
                            const int *vindex)
{
  *timelevels = CCTK_MaxActiveTimeLevelsVI (*cctkGH, *vindex);
}


/********************************************************************
 ********************    QueryGroupStorage   ************************
 ********************************************************************/

int CCTK_QueryGroupStorageI(const cGH *GH, int groupi)
{
  return CCTK_QueryGroupStorageB(GH,groupi,NULL);
}

int CCTK_QueryGroupStorage(const cGH *GH, const char *groupn)
{
  return CCTK_QueryGroupStorageB(GH, -1, groupn);
}



/********************************************************************
 ********************    Group Lower Bound    ***********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GrouplbndGI
   @date       Mon June 19 June 2000
   @author     Gabrielle
   @desc
   Returns the lower bound for a variable group
   @enddesc
@@*/
int CCTK_GrouplbndGI(const cGH *cctkGH,
                     int dim,
                     int *lbnd,
                     int groupindex)
{
  int retval = 0;
  int ierr;
  int usedim = dim;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.lbnd))
  {
    if (data.dim != dim)
    {
      retval = -1;
      usedim = (data.dim < dim) ? data.dim : dim;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GrouplbndGI: Incorrect dimension %d supplied, "
                 "group '%s' has dimension %d, copying %d integers",
                 dim,groupname,data.dim,usedim);
      free (groupname);
    }
    memcpy(lbnd,(const int *)data.lbnd,usedim*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplbndGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            const int *groupindex)
{
  *ierr = CCTK_GrouplbndGI (*cctkGH, *dim, lbnd, *groupindex);
}



int CCTK_GrouplbndGN(const cGH *cctkGH,
                     int dim,
                     int *lbnd,
                     const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GrouplbndGI(cctkGH,dim,lbnd,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplbndGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *ierr = CCTK_GrouplbndGN (*cctkGH, *dim, lbnd, groupname);
  free (groupname);
}



int CCTK_GrouplbndVI(const cGH *cctkGH,
                     int dim,
                     int *lbnd,
                     int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GrouplbndGI(cctkGH,dim,lbnd,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplbndVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            const int *varindex)
{
  *ierr = CCTK_GrouplbndVI (*cctkGH, *dim, lbnd, *varindex);
}



int CCTK_GrouplbndVN(const cGH *cctkGH,
                     int dim,
                     int *lbnd,
                     const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GrouplbndGI(cctkGH,dim,lbnd,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplbndVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplbndVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lbnd,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GrouplbndVN (*cctkGH, *dim, lbnd, varname);
  free (varname);
}



/********************************************************************
 ********************    Group Upper Bound    ***********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GroupubndGI
   @date       Mon June 19 June 2000
   @author     Gabrielle
   @desc
   Returns the upper bound for a variable group
   @enddesc
@@*/
int CCTK_GroupubndGI(const cGH *cctkGH,
                     int dim,
                     int *ubnd,
                     int groupindex)
{
  int retval = 0;
  int ierr;
  int usedim = dim;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.ubnd))
  {
    if (data.dim != dim)
    {
      retval = -1;
      usedim = (data.dim < dim) ? data.dim : dim;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GroupubndGI: Incorrect dimension %d supplied, "
                 "group '%s' has dimension %d, copying %d integers",
                 dim,groupname,data.dim,usedim);
      free (groupname);
    }
    memcpy(ubnd,(const int *)data.ubnd,usedim*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupubndGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            const int *groupindex)
{
  *ierr = CCTK_GroupubndGI (*cctkGH, *dim, ubnd, *groupindex);
}



int CCTK_GroupubndGN(const cGH *cctkGH,
                     int dim,
                     int *ubnd,
                     const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GroupubndGI(cctkGH,dim,ubnd,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupubndGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *ierr = CCTK_GroupubndGN (*cctkGH, *dim, ubnd, groupname);
  free (groupname);
}



int CCTK_GroupubndVI(const cGH *cctkGH,
                     int dim,
                     int *ubnd,
                     int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GroupubndGI(cctkGH,dim,ubnd,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupubndVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            const int *varindex)
{
  *ierr = CCTK_GroupubndVI (*cctkGH, *dim, ubnd, *varindex);
}



int CCTK_GroupubndVN(const cGH *cctkGH,
                     int dim,
                     int *ubnd,
                     const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GroupubndGI(cctkGH,dim,ubnd,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupubndVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupubndVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *ubnd,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GroupubndVN (*cctkGH, *dim, ubnd, varname);
  free (varname);
}



/********************************************************************
 ********************    Group Local Shape    ***********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GrouplshGI
   @date       Mon June 19 June 2000
   @author     Gabrielle
   @desc
   Returns the local shape for a variable group
   @enddesc
@@*/
int CCTK_GrouplshGI(const cGH *cctkGH,
                    int dim,
                    int *lsh,
                    int groupindex)
{
  int retval = 0;
  int ierr;
  int usedim = dim;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.lsh))
  {
    if (data.dim != dim)
    {
      retval = -1;
      usedim = (data.dim < dim) ? data.dim : dim;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GrouplshGI: Incorrect dimension %d supplied, "
                 "group '%s' has dimension %d, copying %d integers",
                 dim,groupname,data.dim,usedim);
      free (groupname);
    }
    memcpy(lsh,(const int *)data.lsh,usedim*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplshGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            const int *groupindex)
{
  *ierr = CCTK_GrouplshGI (*cctkGH, *dim, lsh, *groupindex);
}



int CCTK_GrouplshGN(const cGH *cctkGH,
                     int dim,
                     int *lsh,
                     const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GrouplshGI(cctkGH,dim,lsh,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplshGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *ierr = CCTK_GrouplshGN (*cctkGH, *dim, lsh, groupname);
  free (groupname);
}



int CCTK_GrouplshVI(const cGH *cctkGH,
                     int dim,
                     int *lsh,
                     int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GrouplshGI(cctkGH,dim,lsh,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplshVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            const int *varindex)
{
  *ierr = CCTK_GrouplshVI (*cctkGH, *dim, lsh, *varindex);
}



int CCTK_GrouplshVN(const cGH *cctkGH,
                     int dim,
                     int *lsh,
                     const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GrouplshGI(cctkGH,dim,lsh,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GrouplshVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GrouplshVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *lsh,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GrouplshVN (*cctkGH, *dim, lsh, varname);
  free (varname);
}



/********************************************************************
 ********************    Group local allocated shape    *************
 ********************************************************************/

 /*@@
   @routine    CCTK_GroupashGI
   @date       2011-03-16
   @author     Erik Schnetter
   @desc
   Returns the ash for a variable group
   @enddesc
@@*/
int CCTK_GroupashGI(const cGH *cctkGH,
                    int size,
                    int *ash,
                    int groupindex)
{
  int retval = 0;
  int ierr;
  int usesize = size;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.ash))
  {
    if (data.dim != size)
    {
      retval = -1;
      usesize = (data.dim < size) ? data.dim : size;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GroupashGI: Incorrect size %d supplied, "
                 "group '%s' has dimension %d, copying %d integers",
                 size,groupname,data.dim,usesize);
      free (groupname);
    }
    memcpy(ash,(const int *)data.ash,usesize*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupashGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupashGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *ash,
                            const int *groupindex)
{
  *ierr = CCTK_GroupashGI (*cctkGH, *size, ash, *groupindex);
}



int CCTK_GroupashGN(const cGH *cctkGH,
                    int size,
                    int *ash,
                    const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GroupashGI(cctkGH,size,ash,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupashGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupashGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *ash,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *ierr = CCTK_GroupashGN (*cctkGH, *size, ash, groupname);
  free (groupname);
}



int CCTK_GroupashVI(const cGH *cctkGH,
                    int size,
                    int *ash,
                    int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GroupashGI(cctkGH,size,ash,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupashVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupashVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *ash,
                            const int *varindex)
{
  *ierr = CCTK_GroupashVI (*cctkGH, *size, ash, *varindex);
}



int CCTK_GroupashVN(const cGH *cctkGH,
                    int size,
                    int *ash,
                    const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GroupashGI(cctkGH,size,ash,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupashVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupashVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *ash,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GroupashVN (*cctkGH, *size, ash, varname);
  free (varname);
}



/********************************************************************
 *******************    Group Global Shape    ***********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GroupgshGI
   @date       Mon June 19 June 2000
   @author     Gabrielle
   @desc
   Returns the global shape for a variable group
   @enddesc
@@*/
int CCTK_GroupgshGI(const cGH *cctkGH,
                    int dim,
                    int *gsh,
                    int groupindex)
{
  int retval = 0;
  int ierr;
  int usedim = dim;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.gsh))
  {
    if (data.dim != dim)
    {
      retval = -1;
      usedim = (data.dim < dim) ? data.dim : dim;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GroupgshGI: Incorrect dimension %d supplied, "
                 "group '%s' has dimension %d, copying %d integers",
                 dim,groupname,data.dim,usedim);
      free (groupname);
    }
    memcpy(gsh,(const int *)data.gsh,usedim*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupgshGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            const int *groupindex)
{
  *ierr = CCTK_GroupgshGI (*cctkGH, *dim, gsh, *groupindex);
}



int CCTK_GroupgshGN(const cGH *cctkGH,
                    int dim,
                    int *gsh,
                    const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GroupgshGI(cctkGH,dim,gsh,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupgshGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *ierr = CCTK_GroupgshGN (*cctkGH, *dim, gsh, groupname);
  free (groupname);
}



int CCTK_GroupgshVI(const cGH *cctkGH,
                    int dim,
                    int *gsh,
                    int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GroupgshGI(cctkGH,dim,gsh,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupgshVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            const int *varindex)
{
  *ierr = CCTK_GroupgshVI (*cctkGH, *dim, gsh, *varindex);
}



int CCTK_GroupgshVN(const cGH *cctkGH,
                    int dim,
                    int *gsh,
                    const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GroupgshGI(cctkGH,dim,gsh,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupgshVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupgshVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *gsh,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GroupgshVN (*cctkGH, *dim, gsh, varname);
  free (varname);
}



/********************************************************************
 ********************    Group nGhostzones    ***********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GroupnghostzonesGI
   @date       Mon June 19 June 2000
   @author     Gabrielle
   @desc
   Returns the ghostzone size for a variable group
   @enddesc
@@*/
int CCTK_GroupnghostzonesGI(const cGH *cctkGH,
                            int dim,
                            int *nghostzones,
                            int groupindex)
{
  int retval = 0;
  int ierr;
  int usedim = dim;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.nghostzones))
  {
    if (data.dim != dim)
    {
      retval = -1;
      usedim = (data.dim < dim) ? data.dim : dim;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GroupnghostzonesGI: Incorrect dimension %d supplied, "
                 "group '%s' has dimension %d, copying %d integers",
                 dim,groupname,data.dim,usedim);
      free (groupname);
    }
    memcpy(nghostzones,(const int *)data.nghostzones,usedim*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupnghostzonesGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            const int *groupindex)
{
  *ierr = CCTK_GroupnghostzonesGI (*cctkGH, *dim, nghostzones, *groupindex);
}



int CCTK_GroupnghostzonesGN(const cGH *cctkGH,
                            int dim,
                            int *nghostzones,
                            const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GroupnghostzonesGI(cctkGH,dim,nghostzones,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupnghostzonesGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *ierr = CCTK_GroupnghostzonesGN (*cctkGH, *dim, nghostzones, groupname);
  free (groupname);
}



int CCTK_GroupnghostzonesVI(const cGH *cctkGH,
                            int dim,
                            int *nghostzones,
                            int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GroupnghostzonesGI(cctkGH,dim,nghostzones,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupnghostzonesVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            const int *varindex)
{
  *ierr = CCTK_GroupnghostzonesVI (*cctkGH, *dim, nghostzones, *varindex);
}



int CCTK_GroupnghostzonesVN(const cGH *cctkGH,
                            int dim,
                            int *nghostzones,
                            const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GroupnghostzonesGI(cctkGH,dim,nghostzones,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupnghostzonesVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupnghostzonesVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *dim,
                            int *nghostzones,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GroupnghostzonesVN (*cctkGH, *dim, nghostzones, varname);
  free (varname);
}



/********************************************************************
 ********************        Group bbox       ***********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GroupbboxGI
   @date       Mon June 19 June 2000
   @author     Gabrielle
   @desc
   Returns the bbox array for a variable group
   @enddesc
@@*/
int CCTK_GroupbboxGI(const cGH *cctkGH,
                     int size,
                     int *bbox,
                     int groupindex)
{
  int retval = 0;
  int ierr;
  int usesize = size;  /* Actual number of integers copied */
  char *groupname;
  cGroupDynamicData data;

  ierr = CCTK_GroupDynamicData(cctkGH,groupindex,&data);

  if (ierr == 0 && (data.dim == 0 || data.bbox))
  {
    if (2*data.dim != size)
    {
      retval = -1;
      usesize = (2*data.dim < size) ? 2*data.dim : size;
      groupname = CCTK_GroupName (groupindex);
      CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
                 "CCTK_GroupbboxGI: Incorrect size %d supplied, "
                 "group %s has dimension %d, copying %d integers",
                 size,groupname,data.dim,usesize);
      free (groupname);
    }
    memcpy(bbox,(const int *)data.bbox,usesize*sizeof(int));
  }
  else
  {
    retval = -2;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupbboxGI: Data not available from driver thorn");
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxGI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            const int *groupindex)
{
  *ierr = CCTK_GroupbboxGI (*cctkGH, *size, bbox, *groupindex);
}



int CCTK_GroupbboxGN(const cGH *cctkGH,
                     int size,
                     int *bbox,
                     const char *groupname)
{
  int retval;
  int gindex = CCTK_GroupIndex(groupname);

  if (gindex > -1)
  {
    retval = CCTK_GroupbboxGI(cctkGH,size,bbox,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupbboxGN: Group index not found for %s",groupname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxGN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(groupname)
  *ierr = CCTK_GroupbboxGN(*cctkGH, *size, bbox, groupname);
  free(groupname);
}



int CCTK_GroupbboxVI(const cGH *cctkGH,
                     int size,
                     int *bbox,
                     int varindex)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVarI(varindex);

  if (gindex > -1)
  {
    retval = CCTK_GroupbboxGI(cctkGH,size,bbox,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupbboxVI: Group index not found for variable index %d",
               varindex);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxVI)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            const int *varindex)
{
  *ierr = CCTK_GroupbboxVI (*cctkGH, *size, bbox, *varindex);
}



int CCTK_GroupbboxVN(const cGH *cctkGH,
                     int size,
                     int *bbox,
                     const char *varname)
{
  int retval;
  int gindex = CCTK_GroupIndexFromVar(varname);

  if (gindex > -1)
  {
    retval = CCTK_GroupbboxGI(cctkGH,size,bbox,gindex);
  }
  else
  {
    retval = -4;
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_GroupbboxVN: Group index not found for %s",varname);
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupbboxVN)
                           (int *ierr,
                            const cGH **cctkGH,
                            const int *size,
                            int *bbox,
                            ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname)
  *ierr = CCTK_GroupbboxVN (*cctkGH, *size, bbox, varname);
  free (varname);
}
