 /*@@
   @file      CactusSync.c
   @date      Thu Sep  18 14:27:18 1999
   @author    Gerd Lanfermann
   @desc
              A collection of SyncGroup routines:
              Sync a group by the GROUP INDEX, by the group's VARIABLE NAME,
              by the group's VARIABLE INDEX.

              It ends up calling CCTK_SyncGroup(GH,groupname),
              which is overloaded (currently by PUGH).
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>

#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_Comm.h"
#include "cctk_Groups.h"
#include "cctk_Sync.h"

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(comm_CactusSync_c);


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/* prototypes for external C routines are declared in header cctk_Groups.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupI)
                           (int *ierror, const cGH **GH, const int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupWithVar)
                           (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupWithVarI)
                           (int *ierror, const cGH **GH, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupsI)
                           (int *ierror, const cGH **GH, const int *num_groups,
                            const int *groups);


 /*@@
   @routine    CCTK_SyncGroupI
   @date       Thu Sep  18 14:27:18 1999
   @author     Gerd Lanfermann
   @desc
               Synchronizes a group given by its index.
   @enddesc
   @calls      CCTK_SyncGroupsI

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        group
   @vdesc      group index
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0 for success, or<BR>
               -negative return code of @seeroutine CCTK_SyncGroupsI
   @endreturndesc
@@*/
int CCTK_SyncGroupI (const cGH *GH, int group)
{
  int retval;


  retval = CCTK_SyncGroupsI (GH, 1, &group);

  return (retval == 1 ? 0 : retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupI)
                           (int *ierror, const cGH **GH, const int *group)
{
  CCTK_SyncGroupI (*GH, *group);
  *ierror = 0;
}


 /*@@
   @routine    CCTK_SyncGroupWithVar
   @date       Thu Sep  18 14:27:18 1999
   @author     Gerd Lanfermann
   @desc
               Synchronizes a group which contains the given variable
               (given by its name).
   @enddesc
   @calls      CCTK_GroupIndexFromVarI
               CCTK_VarIndex
               CCTK_SyncGroupI

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        varname
   @vdesc      full variable name
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_SyncGroupI
   @endreturndesc
@@*/
int CCTK_SyncGroupWithVar (const cGH *GH, const char *varname)
{
  int idx;


  idx = CCTK_VarIndex (varname);
  idx = CCTK_GroupIndexFromVarI (idx);
  return (CCTK_SyncGroupI (GH, idx));
}

void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupWithVar)
                           (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (varname);
  *ierror = CCTK_SyncGroupWithVar (*GH, varname);
  free (varname);
}


 /*@@
   @routine    CCTK_SyncGroupWithVarI
   @date       Thu Sep  18 14:27:18 1999
   @author     Gerd Lanfermann
   @desc
               Synchronizes a group which contains the given variable
               (given by its index).
   @enddesc
   @calls      CCTK_GroupIndexFromVarI
               CCTK_SyncGroupI

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        var
   @vdesc      variable index
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_SyncGroupI
   @endreturndesc
@@*/
int CCTK_SyncGroupWithVarI (const cGH *GH, int var)
{
  int idx;


  idx = CCTK_GroupIndexFromVarI (var);
  return (CCTK_SyncGroupI (GH, idx));
}

void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupWithVarI)
                           (int *ierror, const cGH **GH, const int *var)
{
  *ierror = CCTK_SyncGroupWithVarI (*GH, *var);
}


 /*@@
   @routine    CCTK_SyncGroupsI
   @date       Thu Jan 27 18:00:15 2000
   @author     Tom Goodale
   @desc
               Synchronises a list of groups given by their group indices.
   @enddesc
   @calls      CCTK_SyncGroupsByDirI

   @var        GH
   @vdesc      Pointer to Grid Hierachy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        n_groups
   @vdesc      number of groups to synchronize
   @vtype      int
   @vio        in
   @endvar
   @var        groups
   @vdesc      list of group indices
   @vtype      const int *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the total number of groups synchronized, or
               -negative return code of @seeroutine CCTK_SyncGroupsByDirI
   @endreturndesc
@@*/
int CCTK_SyncGroupsI (const cGH *GH,
                      int n_groups,
                      const int *groups)
{
  int retval;


  /* passing NULL as the last argument means: synchronise all directions */
  retval = CCTK_SyncGroupsByDirI (GH, n_groups, groups, NULL);

  return (retval);
}


void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupsI)
     (int *ierror, const cGH **GH, const int *num_groups, const int *groups)
{
  *ierror = CCTK_SyncGroupsI (*GH, *num_groups, groups);
}
