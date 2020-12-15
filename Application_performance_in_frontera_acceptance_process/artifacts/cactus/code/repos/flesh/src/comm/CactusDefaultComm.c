 /*@@
   @file      CactusDefaultComm.c
   @date      Tue Sep 29 15:06:22 1998
   @author    Tom Goodale
   @desc
              Default communication routines.
   @enddesc
   @version   $Id$
 @@*/


#include "cctk_Capabilities.h"
#include "cctk_Flesh.h"
#include "cctk_Groups.h"
#include "cctk_Constants.h"
#include "cctk_Comm.h"
#include "cctk_Sync.h"
#include "cctk_WarnLevel.h"

#include "cctki_GHExtensions.h"

#include "cctk_ParamCheck.h"

#include "CactusMainDefaults.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_CAPABILITY_MPI
#  include <mpi.h>
#endif

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(comm_CactusDefaultComm_c);

/********************************************************************
 *********************     Global Data   *****************************
 ********************************************************************/

/* FIXME:  This should be in a header somewhere */
#ifdef HAVE_CAPABILITY_MPI
extern char cctki_MPI_Active;
#endif

/********************************************************************
 *********************     Local Definitions   **********************
 ********************************************************************/

#ifdef HAVE_CAPABILITY_MPI
#define CACTUS_MPI_ERROR(xf)                                                  \
          {                                                                   \
            int errcode;                                                      \
                                                                              \
                                                                              \
            if ((errcode = (xf)) != MPI_SUCCESS)                              \
            {                                                                 \
              char mpi_error_string[MPI_MAX_ERROR_STRING+1];                  \
              int resultlen;                                                  \
                                                                              \
                                                                              \
              MPI_Error_string(errcode, mpi_error_string, &resultlen);        \
              fprintf(stderr, "MPI Call %s returned error code %d (%s)\n",    \
                              #xf, errcode, mpi_error_string);                \
              fprintf(stderr, "At line %d of file %s\n",                      \
                              __LINE__, __FILE__);                            \
            }                                                                 \
          }
#endif

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
cGH *CactusDefaultSetupGH(tFleshConfig *config, int convergence_level);
int CactusDefaultMyProc(const cGH *GH);
int CactusDefaultnProcs(const cGH *GH);
int CactusDefaultExit(cGH *GH, int retval);
int CactusDefaultAbort(cGH *GH, int retval);
int CactusDefaultBarrier(const cGH *GH);

int CactusDefaultSyncGroup (const cGH *GH, const char *groupname);
int CactusDefaultSyncGroupsByDirI(const cGH *GH, int num_groups,
                                  const int *groups,
                                  const int *directions);
int CactusDefaultEnableGroupStorage(const cGH *GH, const char *group);
int CactusDefaultDisableGroupStorage(const cGH *GH, const char *group);
int CactusDefaultGroupStorageIncrease(const cGH *GH, int n_groups,
                                      const int *groups, const int *timelevels,
                                      int *status);
int CactusDefaultGroupStorageDecrease(const cGH *GH, int n_groups,
                                      const int *groups, const int *timelevels,
                                      int *status);
int CactusDefaultInterpGridArrays (const cGH *GH, int N_dims,
                                   int local_interp_handle,
                                   int param_table_handle,
                                   int coord_system_handle, int N_points,
                                   int interp_coords_type,
                                   const void *const interp_coords[],
                                   int N_input_arrays,
                                   const CCTK_INT input_array_indices[],
                                   int N_output_arrays,
                                   const CCTK_INT output_array_types[],
                                   void *const output_arrays[]);


 /*@@
   @routine    CactusDefaultSetupGH
   @date       Tue Sep 29 15:06:22 1998
   @author     Tom Goodale
   @desc
               Default cactus SetupGH routine.
   @enddesc
   @calls      CCTK_MaxGFDim
               CCTK_NumVars
               CCTK_DeclaredTimeLevelsVI
               CCTK_NumGroups
               CCTKi_SetupGHExtensions

   @var        config
   @vdesc      Pointer to flesh configuration environment
   @vtype      tFleshConfig *
   @vio        in
   @endvar
   @var        convergence_level
   @vdesc      convergence leve of new cGH
   @vtype      int
   @vio        in
   @endvar

   @returntype cGH *
   @returndesc
               the pointer to the new cGH structure,
               or NULL if memory allocation failed
   @endreturndesc
@@*/
cGH *CactusDefaultSetupGH(tFleshConfig *config, int convergence_level)
{
  cGH *thisGH;
  int n_groups, n_variables;
  int variable, ntimelevels, cctk_dim;


  /* Put this in for the moment until parameter stuff is done. */
  if (convergence_level > 0)
  {
    return (NULL);
  }

  /* Initialise this since it is used later and in exceptional
   * circumstances might not be initialsed beforehand.
   */
  variable = -1;

  /* Create a new Grid Hierarchy */
  thisGH = malloc(sizeof(cGH));
  if (thisGH)
  {
    thisGH->cctk_dim = CCTK_MaxGFDim();

    /* Need this to be at least one otherwise the memory allocation will fail */
    cctk_dim = thisGH->cctk_dim;
    if (thisGH->cctk_dim == 0)
    {
      cctk_dim = 1;
    }
    thisGH->cctk_iteration    = 0;
    thisGH->cctk_gsh          = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_lsh          = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_lbnd         = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_ubnd         = malloc(cctk_dim*sizeof(int));

    thisGH->cctk_ash          = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_to           = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_from         = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_bbox         = malloc(2*cctk_dim*sizeof(int));
    thisGH->cctk_nghostzones  = malloc(2*cctk_dim*sizeof(int));
    thisGH->cctk_levfac       = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_levoff       = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_levoffdenom  = malloc(cctk_dim*sizeof(int));
    thisGH->cctk_delta_space  = malloc(cctk_dim*sizeof(CCTK_REAL));
    /* FIXME : Next line goes when coords are done properly */
    thisGH->cctk_origin_space = malloc(cctk_dim*sizeof(CCTK_REAL));

    thisGH->cctk_delta_time = 1;
    thisGH->cctk_timefac = 1;
    thisGH->cctk_convlevel = 0;

    n_variables = CCTK_NumVars();

    /* Allocate memory for the variable data pointers.
     * Note we want at least one to prevent memory allocation from failing!
     */
    thisGH->data = malloc((n_variables ? n_variables:1)*sizeof(void **));
    if (thisGH->data)
    {
      for(variable = 0; variable < n_variables; variable++)
      {
        ntimelevels = CCTK_DeclaredTimeLevelsVI(variable);

        thisGH->data[variable] = calloc(ntimelevels, sizeof(void *));
        if (thisGH->data[variable] == NULL)
        {
          break;
        }
      }
    }

    thisGH->extensions = NULL;

    /* Allocate memory for the group data pointers.
     * Note we want at least one to prevent memory allocation from failing !
     */
    n_groups = CCTK_NumGroups();
    thisGH->GroupData = malloc((n_groups ? n_groups:1)*sizeof(cGHGroupData));
  }

  if (thisGH &&
     thisGH->cctk_gsh &&
     thisGH->cctk_lsh &&
     thisGH->cctk_lbnd &&
     thisGH->cctk_ubnd &&
     thisGH->cctk_ash &&
     thisGH->cctk_from &&
     thisGH->cctk_to &&
     thisGH->cctk_bbox &&
     thisGH->cctk_nghostzones &&
     thisGH->cctk_levfac &&
     thisGH->cctk_levoff &&
     thisGH->cctk_levoffdenom &&
     thisGH->cctk_delta_space &&
     thisGH->cctk_origin_space &&
     thisGH->data &&
     variable == n_variables &&
     thisGH->GroupData)
  {
    /* Traverse list of GH setup routines. */
    CCTKi_SetupGHExtensions(config, convergence_level, thisGH);
  }
  else
  {
    /* FIXME: should free potentially allocated memory for this GH */
    thisGH = NULL;
  }

  return (thisGH);
}


 /*@@
   @routine    CactusDefaultMyProc
   @date       Tue Jan 23 1999
   @author     Gabrielle Allen
   @desc
               Default Cactus MyProc routine.
   @enddesc
   @calls      CCTK_ParamChecking
               MPI_Comm_rank

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        unused
   @endvar

   @returntype int
   @returndesc
               the processor number of the caller
   @endreturndesc
@@*/
int CactusDefaultMyProc (const cGH *GH)
{
  int myproc;


  /* avoid compiler warning about unused parameter */
  (void) (GH + 0);

  myproc = 0;
#ifdef HAVE_CAPABILITY_MPI
  if (! CCTK_ParamChecking() && cctki_MPI_Active)
  {
    CACTUS_MPI_ERROR (MPI_Comm_rank (MPI_COMM_WORLD, &myproc));
  }
#endif

  return (myproc);
}


 /*@@
   @routine    CactusDefaultnProcs
   @date       Tue Jan 23 1999
   @author     Gabrielle Allen
   @desc
               Default Cactus nProcs routine.
   @enddesc
   @calls      CCTK_ParamCheckNProcs
               MPI_Comm_size

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        unused
   @endvar

   @returntype int
   @returndesc
               the total number of processors
   @endreturndesc
@@*/
int CactusDefaultnProcs (const cGH *GH)
{
  int nprocs;


  /* avoid compiler warning about unused parameter */
  (void) (GH + 0);

  if (CCTK_ParamChecking ())
  {
    nprocs = CCTK_ParamCheckNProcs ();
  }
  else
  {
    nprocs = 1;
#ifdef HAVE_CAPABILITY_MPI
    if (cctki_MPI_Active)
    {
      CACTUS_MPI_ERROR (MPI_Comm_size (MPI_COMM_WORLD, &nprocs));
    }
#endif
  }

  return (nprocs);
}


 /*@@
   @routine    CactusDefaultExit
   @date       Tue Apr 18 15:21:15 2000
   @author     Gerd Lanfermann
   @desc
               The default for when people call CCTK_Exit.
   @enddesc
   @calls      MPI_Finalize
               exit

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        unused
   @endvar
   @var        retval
   @vdesc      return code to exit with
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               This function should never return.
               But if it does it will return 0.
   @endreturndesc
@@*/
int CactusDefaultExit (cGH *GH, int retval)
{
  /* avoid compiler warning about unused parameter */
  (void) (GH + 0);

#ifdef HAVE_CAPABILITY_MPI
  if (cctki_MPI_Active)
  {
    CACTUS_MPI_ERROR (MPI_Finalize ());
  }
#endif
  exit (retval);
}


 /*@@
   @routine    CactusDefaultAbort
   @date       Saturday July 15 2000
   @author     Gabrielle Allen
   @desc
               The default for when people call CCTK_Abort.
   @enddesc
   @calls      MPI_Abort
               exit

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        unused
   @endvar
   @var        retval
   @vdesc      return code to abort with
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               This function should never return.
               But if it does it will return the return code to abort with.
   @endreturndesc
@@*/
int CactusDefaultAbort (cGH *GH, int retval)
{
  /* avoid compiler warning about unused parameter */
  (void) (GH + 0);

#ifdef HAVE_CAPABILITY_MPI
  if (cctki_MPI_Active)
  {
    /* flush stdout/stderr and then wait a few seconds before calling
       MPI_Abort()
       This cures a problem where processor 0 is slightly behind, eg. with
       activating thorns, and would print an error message (to stdout!!)
       about a thorn missing only a few milliseconds later. If other processors
       call CCTK_Abort() before that those messages wouldn't be seen. */
    fflush (stdout);
    fflush (stderr);

#ifdef HAVE_UNISTD_H
    sleep (5);
#endif

    CACTUS_MPI_ERROR (MPI_Abort (MPI_COMM_WORLD, retval));
  }
#else
  /* FIXME */
  /*abort();*/
#endif
  exit (retval);
}


 /*@@
   @routine    CactusDefaultBarrier
   @date       Tue Apr 18 15:21:42 2000
   @author     Tom Goodale
   @desc
               The default for when people call CCTK_Barrier
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        unused
   @endvar

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
int CactusDefaultBarrier (const cGH *GH)
{
  /* avoid compiler warning about unused parameter */
  (void) (GH + 0);

  return (0);
}


 /*@@
   @routine    CactusDefaultSyncGroup
   @date       Mon 3 July 2006
   @author     Thomas Radke
   @desc
               Default routine for synchronising a single group.

               If this function has not been overloaded by a driver,
               it will simply call CCTK_SyncGroupsByDirI() (which then must be
               overloaded by a driver).
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        groupname
   @vdesc      full name of the group to be synchronised
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success,
               or (negative) return code of @seeroutine CCTK_SyncGroupsByDirI
   @endreturndesc
 @@*/
int CactusDefaultSyncGroup (const cGH *GH, const char *groupname)
{
  int group, retval;


  group = CCTK_GroupIndex (groupname);
  retval = CCTK_SyncGroupsByDirI (GH, 1, &group, NULL);

  return (retval == 1 ? 0 : retval);
}


 /*@@
   @routine    CactusDefaultSyncGroupsByDirI
   @date       Mon 3 July 2006
   @author     Thomas Radke
   @desc
               Default groups synchronisation routine.

               If this function has not been overloaded by a driver,
               it will loop over all groups to be synchronised and
               for each of them call CCTK_SyncGroupI() which itself calls
               the routine CCTK_SyncGroup() (which then must be overloaded
               by a driver).
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_groups
   @vdesc      number of groups to be synchronised
   @vtype      int
   @vio        in
   @endvar
   @var        groups
   @vdesc      list of indices of groups to be synchronised
   @vtype      int *
   @vio        in
   @endvar
   @var        directions
   @vdesc      (optional) array of dimensions which should be synchronised
   @vtype      int *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the number of groups which were synchronised
   @endreturndesc
 @@*/
int CactusDefaultSyncGroupsByDirI (const cGH *GH, int num_groups,
                                   const int *groups,
                                   const int *directions)
{
  int group, retval = 0;
  char *groupname;
  static int user_has_been_notified = 0;


  /* individual directions aren't supported in the CCTK_SyncGroup* interface */
  if (directions != NULL)
  {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, "Cactus",
                "Synchronisation of individual directions isn't supported "
                "with no driver overloading CCTK_SyncGroupsByDirI().");
  }

  /* if CCTK_SyncGroup() hasn't been overloaded then this is a no-op */
  if (CCTK_SyncGroup != CactusDefaultSyncGroup)
  {
    /* on the first time through, warn the user if a driver overloaded
       the (deprecated) routine CCTK_SyncGroup() but not the newer one
       CCTK_SyncGroupsByDirI() */
    if (! user_has_been_notified)
    {
      CCTK_VWarn (CCTK_WARN_COMPLAIN, __LINE__, __FILE__, "Cactus",
                  "Overloading CCTK_SyncGroup() is deprecated. "
                  "Please have your driver thorn updated to overload "
                  "CCTK_SyncGroupsByDirI() instead !");
    }
    user_has_been_notified = 1;

    /* synchronise all groups one by one */
    for (group = 0; group < num_groups; group++)
    {
      groupname = CCTK_GroupName (groups[group]);
      if (CCTK_SyncGroup (GH, groupname) == 0)
      {
        retval++;
      }
      free (groupname);
    }
  }

  return (retval);
}


 /*@@
   @routine    CactusDefaultEnableGroupStorage
   @date       Wed Apr  3 17:01:22 2002
   @author     Tom Goodale
   @desc
               Default enable group storage routine.

               The enable group storage routine should allocate memory
               for a group and return the previous status of that memory.

               This default checks for the presence of the newer
               GroupStorageIncrease function, and if that is not available
               it flags an error. If it is available it makes a call to it,
               passing -1 as the timelevel argument, which is supposed to mean
               enable all timelevels, i.e. preserving this obsolete behaviour.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        inout
   @vcomment
               A driver should replace the appropriate GV pointers on this
               structure when they change the storage state of a GV.
   @endvar
   @var        groupname
   @vdesc      name of the group to allocate storage for
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               true(0) or false(1) if the group previously had storage or not,
               -1 if group increase storage routine wasn't overloaded
   @endreturndesc
 @@*/
int CactusDefaultEnableGroupStorage(const cGH *GH, const char *groupname)
{
  int group, timelevel, retval;


  /* Has the increase group storage routine been overloaded ? */
  if (CCTK_GroupStorageIncrease != CactusDefaultGroupStorageIncrease)
  {
    group = CCTK_GroupIndex(groupname);
    timelevel = -1;
    CCTK_GroupStorageIncrease(GH, 1, &group, &timelevel, &retval);
  }
  else
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "No driver thorn activated to provide storage for variables");
    retval = -1;
  }

  return retval;
}

 /*@@
   @routine    CactusDefaultDisableGroupStorage
   @date       Wed Apr  3 17:01:22 2002
   @author     Tom Goodale
   @desc
               Default disable group storage routine.

               The disable group storage routine should deallocate memory
               for a group and return the previous status of that memory.

               This default checks for the presence of the newer
               GroupStorageDecrease function, and if that is not available
               it flags an error. If it is available it makes a call to it,
               passing -1 as the timelevel argument, which is supposed to mean
               disable all timelevels, i.e. preserving this obsolete behaviour.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        inout
   @vcomment
               A driver should replace the appropriate GV pointers on this
               structure when they change the storage state of a GV.
   @endvar
   @var        groupname
   @vdesc      name of group to deallocate storage for
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               true(0) or false(1) if the group previously had storage or not,
               -1 if group decrease storage routine wasn't overloaded
   @endreturndesc
 @@*/
int CactusDefaultDisableGroupStorage(const cGH *GH, const char *groupname)
{
  int group, timelevel, retval;


  /* Has the decrease group storage routine been overloaded ? */
  if (CCTK_GroupStorageDecrease != CactusDefaultGroupStorageDecrease)
  {
    group = CCTK_GroupIndex(groupname);
    timelevel = -1;
    CCTK_GroupStorageDecrease(GH, 1, &group, &timelevel, &retval);
  }
  else
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "No driver thorn activated to provide storage for variables");
    retval = -1;
  }

  return retval;
}

 /*@@
   @routine    CactusDefaultGroupStorageIncrease
   @date       Wed Apr  3 17:01:22 2002
   @author     Tom Goodale
   @desc
               Default increase group storage routine.

               The increase group storage routine should increase the allocated
               memory to the specified number of timelevels of each listed
               group, returning the previous number of timelevels enabled for
               that group in the status array, if that is not NULL.
               It should never decrease the number of timelevels enabled,
               i.e. if it is asked to enable less timelevels than are already
               enabled it should not change the storage for that group.

               This default checks for the presence of the older
               EnableGroupStorage function, and if that is not available it
               flags an error. If it is available it makes a call to it,
               and puts its return value in the status flag for the group.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        inout
   @vcomment
               A driver should replace the appropriate GV pointers on this
               structure when they change the storage state of a GV.
   @endvar
   @var        n_groups
   @vdesc      number of groups in group array
   @vtype      int
   @vio        in
   @endvar
   @var        groups
   @vdesc      list of group indices to allocate storage for
   @vtype      const int *
   @vio        in
   @endvar
   @var        timelevels
   @vdesc      number of timelevels to allocate storage for for each group
   @vtype      const int *
   @vio        in
   @endvar
   @var        status
   @vdesc      optional return array which, if not NULL, will, on return,
               contain the number of timelevels which were previously allocated
               storage for each group
   @vtype      const int *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               The total number of timelevels with storage enabled for all
               groups queried or modified.
   @endreturndesc
 @@*/
int CactusDefaultGroupStorageIncrease (const cGH *GH, int n_groups,
                                       const int *groups, const int *timelevels,
                                       int *status)
{
  int i, value, retval;
  char *gname;


  /* Has the normal group storage been overloaded ? */
  if (CCTK_EnableGroupStorage != CactusDefaultEnableGroupStorage)
  {
    for(i = retval = 0; i < n_groups; i++)
    {
      if (groups[i] >= 0)
      {
        /* Since the old enable and disable group storage just returned true or
         * false and did all timelevels, only enable storage if timelevels is
         * not 0
         */
        if (CCTK_QueryGroupStorageI(GH, groups[i]))
        {
          value = CCTK_DeclaredTimeLevelsVI(groups[i]);
        }
        else
        {
          value = 0;
        }
        if (timelevels[i] != 0)
        {
          gname = CCTK_GroupName(groups[i]);
          CCTK_EnableGroupStorage(GH, gname);
          free (gname);
        }
        retval += value;
        if (status)
        {
          status[i] = value;
        }
      }
    }
  }
  else
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "No driver thorn activated to provide storage for variables");
    retval = -1;
  }

  return retval;
}

 /*@@
   @routine    CactusDefaultGroupStorageDecrease
   @date       Wed Apr  3 17:01:22 2002
   @author     Tom Goodale
   @desc
               Default decrease group storage routine.

               The decrease group storage routine should decrease the memory
               allocated to the specified number of timelevels for each listed
               group, returning the previous number of timelevels enabled for
               that group in the status array, if that is not NULL.
               It should never increase the number of timelevels enabled,
               i.e. if it is asked to reduce to more timelevels than are enabled
               it should not change the storage for that group.

               This default checks for the presence of the older
               DisableGroupStorage function, and if that is not available it
               flags an error. If it is available it makes a call to it,
               and puts its return value in the status flag for the group.

               A driver should replace the appropriate GV pointers on the cGH
               structure when they change the storage state of a GV.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        inout
   @endvar
   @var        n_groups
   @vdesc      number of groups in group array
   @vtype      int
   @vio        in
   @endvar
   @var        groups
   @vdesc      list of group indices to reduce storage for
   @vtype      const int *
   @vio        in
   @endvar
   @var        timelevels
   @vdesc      number of timelevels to reduce storage for for each group
   @vtype      const int *
   @vio        in
   @endvar
   @var        status
   @vdesc      optional return array which, if not NULL, will, on return,
               contain the number of timelevels which were previously allocated
               storage for each group
   @vtype      const int *
   @vio        out
   @endvar

   @returntype int
   @returndesc
               The total number of timelevels with storage enabled
               for all groups queried or modified.
   @endreturndesc
 @@*/
int CactusDefaultGroupStorageDecrease (const cGH *GH, int n_groups,
                                       const int *groups, const int *timelevels,
                                       int *status)
{
  int i, value, retval;
  char *gname;


  /* Has the normal group storage been overloaded ? */
  if (CCTK_DisableGroupStorage != CactusDefaultDisableGroupStorage)
  {
    for(i = retval = 0; i < n_groups; i++)
    {
      /* Bogus entries in group array are marked with -1.*/
      if (groups[i] >= 0)
      {
        /* Since the old enable and disable group storage just returned true or
         * false and did all timelevels, only disable storage if timelevels is
         * 0
         */
        if (CCTK_QueryGroupStorageI(GH, groups[i]))
        {
          value = CCTK_DeclaredTimeLevelsVI(groups[i]);
        }
        else
        {
          value = 0;
        }
        if (timelevels[i] == 0)
        {
          gname = CCTK_GroupName(groups[i]);
          CCTK_DisableGroupStorage(GH, gname);
          free (gname);
        }
        retval += value;
        if (status)
        {
          status[i] = value;
        }
      }
    }
  }
  else
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "No driver thorn activated to provide storage for variables");
    retval = -1;
  }

  return retval;
}

 /*@@
   @routine    CactusDefaultQueryMaxTimeLevels
   @date       Mon Mar 10 19:16:22 PDT 2014
   @author     Roland Haas
   @desc
               Default routine to query the size of cctkGH->data.

               Using GroupStorageIncrease any number of timelevels can be
               created, this routine returns the total number created.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        inout
   @endvar
   @var        n_groups
   @vdesc      number of groups in group array
   @vtype      int
   @vio        in
   @endvar
   @var        groups
   @vdesc      list of group indices to reduce storage for
   @vtype      const int *
   @vio        in
   @endvar
   @var        status
   @vdesc      on return,
               contain the number of timelevels for which storage has ever been
               allocated for each group
   @vtype      const int *
   @vio        out
   @endvar

   @returntype int
   @returndesc
               Negative uppon errors.
   @endreturndesc
 @@*/
int CactusDefaultQueryMaxTimeLevels (const cGH *GH, int n_groups,
                                     const int *groups, int *status)
{
  assert (status || n_groups == 0);
  assert (groups || n_groups == 0);
  assert (GH);

  /* A driver that supports arbitrary time levels needs to overload the
   * routine, otherwise this should be fine.
   */
  for(int g = 0; g < n_groups; g++)
    status[g] = CCTK_DeclaredTimeLevelsGI(groups[g]);

  (void)(GH);
  return 0;
}


/*@@
   @routine    CactusDefaultInterpGridArrays
   @date       Mon 16 Dec 2002
   @author     Thomas Radke
   @desc
               Default grid array interpolation routine.

               This routine must be overloaded by a driver thorn
               otherwise it will print an appropriate warning and stop the code.
   @enddesc

   @returntype int
   @returndesc
               -1 in all cases
   @endreturndesc
@@*/
int CactusDefaultInterpGridArrays (const cGH *GH, int N_dims,
                                   int local_interp_handle,
                                   int param_table_handle,
                                   int coord_system_handle,
                                   int N_points, int interp_coords_type,
                                   const void *const interp_coords[],
                                   int N_input_arrays,
                                   const CCTK_INT input_array_indices[],
                                   int N_output_arrays,
                                   const CCTK_INT output_array_types[],
                                   void *const output_arrays[])
{
  /* avoid warnings about unused parameters */
  (void) (GH + 0);
  (void) (N_dims + 0);
  (void) (local_interp_handle + 0);
  (void) (param_table_handle + 0);
  (void) (coord_system_handle + 0);
  (void) (N_points + 0);
  (void) (interp_coords_type + 0);
  (void) (interp_coords + 0);
  (void) (N_input_arrays + 0);
  (void) (input_array_indices + 0);
  (void) (N_output_arrays + 0);
  (void) (output_array_types + 0);
  (void) (output_arrays + 0);

  CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
              "No driver thorn activated to provide an interpolation routine "
              "for grid arrays");

  return (-1);
}
