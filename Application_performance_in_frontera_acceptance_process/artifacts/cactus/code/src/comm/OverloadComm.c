 /*@@
   @file      Overload.c
   @date      Wed Feb  3 23:27:18 1999
   @author    Tom Goodale
   @desc
              Contains routines to overload the communication functions.
              Uses the overload macros to make sure of consistency and
              to save typing !
   @enddesc
   @version   $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_WarnLevel.h"
#include "CactusRegister.h"
#include "OverloadMacros.h"

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(comm_OverloadComm_c);


/* Define the prototypes for the dummy functions. */
#define OVERLOADABLE(name) OVERLOADABLE_DUMMYPROTOTYPE(name)

  /* Deal seperately with the SetupGH routine */
#define CCTKi_DummyAbort    CactusDefaultAbort
#define CCTKi_DummyBarrier  CactusDefaultBarrier
#define CCTKi_DummyExit     CactusDefaultExit
#define CCTKi_DummyMyProc   CactusDefaultMyProc
#define CCTKi_DummynProcs   CactusDefaultnProcs
#define CCTKi_DummySetupGH  CactusDefaultSetupGH

  /* Need to do some stuff to make sure old routines still work. */
#define CCTKi_DummySyncGroup             CactusDefaultSyncGroup
#define CCTKi_DummySyncGroupsByDirI      CactusDefaultSyncGroupsByDirI
#define CCTKi_DummyEnableGroupStorage    CactusDefaultEnableGroupStorage
#define CCTKi_DummyDisableGroupStorage   CactusDefaultDisableGroupStorage
#define CCTKi_DummyGroupStorageIncrease  CactusDefaultGroupStorageIncrease
#define CCTKi_DummyGroupStorageDecrease  CactusDefaultGroupStorageDecrease
#define CCTKi_DummyQueryMaxTimeLevels    CactusDefaultQueryMaxTimeLevels

#define CCTKi_DummyInterpGridArrays  CactusDefaultInterpGridArrays

#include "CommOverloadables.h"

  /* Reset the #define to prevent complications. */
#undef CCTKi_DummySetupGH  
#undef CCTKi_DummyMyProc  
#undef CCTKi_DummynProcs
#undef CCTKi_DummyBarrier
#undef CCTKi_DummyExit
#undef CCTKi_DummyAbort

#undef CCTKi_DummySyncGroup
#undef CCTKi_DummySyncGroupsByDirI
#undef CCTKi_DummyEnableGroupStorage
#undef CCTKi_DummyDisableGroupStorage
#undef CCTKi_DummyGroupStorageIncrease
#undef CCTKi_DummyGroupStorageDecrease
#undef CCTKi_DummyQueryMaxTimeLevels

#undef CCTKi_DummyInterpGridArrays

#undef OVERLOADABLE

/* Create the overloadable function variables and the 
 * functions allowing the variables to be set.
 */
#define OVERLOADABLE(name) OVERLOADABLE_FUNCTION(name)

#include "CommOverloadables.h"

#undef OVERLOADABLE


 /*@@
   @routine    CCTKi_SetupCommFunctions
   @date       Thu Feb  4 08:21:26 1999
   @author     Tom Goodale
   @desc 
   Set any comm function which hasn't been overloaded to the default.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int CCTKi_SetupCommFunctions(void)
{

#define OVERLOADABLE(name) OVERLOADABLE_INITIALISE(name)

  /* Deal seperately with the SetupGH routine */
#define CCTKi_DummyAbort   CactusDefaultAbort
#define CCTKi_DummyBarrier CactusDefaultBarrier
#define CCTKi_DummyExit    CactusDefaultExit
#define CCTKi_DummyMyProc  CactusDefaultMyProc
#define CCTKi_DummynProcs  CactusDefaultnProcs
#define CCTKi_DummySetupGH CactusDefaultSetupGH

#define CCTKi_DummySyncGroup             CactusDefaultSyncGroup
#define CCTKi_DummySyncGroupsByDirI      CactusDefaultSyncGroupsByDirI
#define CCTKi_DummyEnableGroupStorage    CactusDefaultEnableGroupStorage
#define CCTKi_DummyDisableGroupStorage   CactusDefaultDisableGroupStorage
#define CCTKi_DummyGroupStorageIncrease  CactusDefaultGroupStorageIncrease
#define CCTKi_DummyGroupStorageDecrease  CactusDefaultGroupStorageDecrease
#define CCTKi_DummyQueryMaxTimeLevels    CactusDefaultQueryMaxTimeLevels

#define CCTKi_DummyInterpGridArrays  CactusDefaultInterpGridArrays

#include "CommOverloadables.h"

  /* Reset the #define to prevent complications. */
#undef CCTKi_DummyAbort
#undef CCTKi_DummyBarrier
#undef CCTKi_DummyExit
#undef CCTKi_DummyMyProc  
#undef CCTKi_DummynProcs
#undef CCTKi_DummySetupGH  

#undef CCTKi_DummySyncGroup
#undef CCTKi_DummySyncGroupsByDirI
#undef CCTKi_DummyEnableGroupStorage
#undef CCTKi_DummyDisableGroupStorage
#undef CCTKi_DummyGroupStorageIncrease
#undef CCTKi_DummyGroupStorageDecrease
#undef CCTKi_DummyQueryMaxTimeLevels

#undef CCTKi_DummyInterpGridArrays

#undef OVERLOADABLE

  return 0;
}

/* Create the dummy function prototypes. */
#define OVERLOADABLE(name) OVERLOADABLE_DUMMYPROTOTYPE(name)

#include "CommOverloadables.h"

#undef OVERLOADABLE


/* Create the dummy functions. */
#define OVERLOADABLE(name) OVERLOADABLE_DUMMY(name)

#include "CommOverloadables.h"

#undef OVERLOADABLE



/* Fortran bindings prototypes for the comm functions */
int CCTK_FCALL CCTK_FNAME (CCTK_nProcs) (const cGH **GH);
int CCTK_FCALL CCTK_FNAME (CCTK_MyProc) (const cGH **GH);
void CCTK_FCALL CCTK_FNAME (CCTK_Barrier) (int *ierror, const cGH **GH);
void CCTK_FCALL CCTK_FNAME (CCTK_Exit) (int *ierror, cGH **GH, const int *retval);
void CCTK_FCALL CCTK_FNAME (CCTK_Abort) (int *ierror, cGH **GH, const int *retval);
void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroup) (int *ierror, cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupsByDirI) (int *ierror, cGH **GH, const int *num_groups, const int *groups, const int *directions);
void CCTK_FCALL CCTK_FNAME (CCTK_EnableGroupComm) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_DisableGroupComm) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_EnableGroupStorage) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_QueryMaxTimeLevels) (int *ierror, const cGH **GH, const int *num_groups, const int *groups, int *status);
void CCTK_FCALL CCTK_FNAME (CCTK_DisableGroupStorage) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_QueryGroupStorage) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG);


/* Fortran bindings definitions for the comm functions */
int CCTK_FCALL CCTK_FNAME (CCTK_nProcs) (const cGH **GH)
{
  return (CCTK_nProcs (*GH));
}
 
int CCTK_FCALL CCTK_FNAME (CCTK_MyProc) (const cGH **GH)
{
  return (CCTK_MyProc (*GH));
}
 
void CCTK_FCALL CCTK_FNAME (CCTK_Barrier) (int *ierror, const cGH **GH)
{
  *ierror = CCTK_Barrier (*GH);
}

void CCTK_FCALL CCTK_FNAME (CCTK_Exit) (int *ierror, cGH **GH, const int *retval)
{
  *ierror = CCTK_Exit (*GH, *retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_Abort) (int *ierror, cGH **GH, const int *retval)
{
  *ierror = CCTK_Abort (*GH, *retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroup) (int *ierror, cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group_name)
  *ierror = CCTK_SyncGroup (*GH, group_name);
  free (group_name); 
}

void CCTK_FCALL CCTK_FNAME (CCTK_SyncGroupsByDirI) (int *ierror, cGH **GH,
                                                    const int *num_groups,
                                                    const int *groups,
                                                    const int *directions)
{
  *ierror = CCTK_SyncGroupsByDirI (*GH, *num_groups, groups, directions);
}

void CCTK_FCALL CCTK_FNAME (CCTK_EnableGroupComm) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group_name)
  *ierror = CCTK_EnableGroupComm (*GH, group_name); 
  free (group_name);
}

void CCTK_FCALL CCTK_FNAME (CCTK_DisableGroupComm) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group_name)
  *ierror = CCTK_DisableGroupComm (*GH, group_name); 
  free (group_name);
}

void CCTK_FCALL CCTK_FNAME (CCTK_EnableGroupStorage) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group_name)
  *ierror = CCTK_EnableGroupStorage (*GH, group_name);
  free (group_name);
}

void CCTK_FCALL CCTK_FNAME (CCTK_DisableGroupStorage) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group_name)
  *ierror = CCTK_DisableGroupStorage (*GH, group_name);
  free (group_name);
}

void CCTK_FCALL CCTK_FNAME (CCTK_QueryMaxTimeLevels) (int *ierror, const cGH **GH,
                                                      const int *num_groups,
                                                      const int *groups,
                                                      int *status)
{
  *ierror = CCTK_QueryMaxTimeLevels (*GH, *num_groups, groups, status);
}

void CCTK_FCALL CCTK_FNAME (CCTK_QueryGroupStorage) (int *ierror, const cGH **GH, ONE_FORTSTRING_ARG)
{
  extern int CCTK_QueryGroupStorage (const cGH *, const char *);
  ONE_FORTSTRING_CREATE (group_name)
  *ierror = CCTK_QueryGroupStorage (*GH, group_name);
  free (group_name);
}
