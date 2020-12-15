/*@@
   @header    cctk.h
   @date      Tue Jan 26 17:29:34 1999
   @author    Tom Goodale
   @desc
              Main include file for the CCTK.
              All thorns should include this.
   @enddesc
   @version   $Header$
 @@*/

#ifndef _CCTK_CORE_H_
#define _CCTK_CORE_H_ 1

/* Grab the main configuration info. */
#include "cctk_Config.h"

/* Include the constants */
#include "cctk_Constants.h"

/* Include definitions provided by capabilities */
#include "cctk_Capabilities.h"

/* Define some stuff */

#ifdef FCODE

#include "cctk_Faces.h"
#include "cctk_Interp.h"
#include "cctk_Loop.h"
#include "cctk_WarnLevel.h"

#define CCTK_PRINTSEPARATOR\
  print '("--------------------------------------------------------------------------------")'

#define _CCTK_FARGUMENTS \
cctk_dim,cctk_gsh,cctk_lsh,cctk_lbnd,cctk_ubnd,cctk_ash,cctk_from,cctk_to,\
cctk_bbox,cctk_delta_time,cctk_time,cctk_delta_space,cctk_origin_space,\
cctk_levfac,cctk_levoff,cctk_levoffdenom,cctk_timefac,cctk_convlevel,\
cctk_convfac,cctk_nghostzones,cctk_iteration,cctkGH,\
cctk_ash1,cctk_ash2,cctk_ash3

#define _DECLARE_CCTK_ARGUMENTS _DECLARE_CCTK_FARGUMENTS
#define _DECLARE_CCTK_FARGUMENTS &&\
        CCTK_DECLARE(INTEGER,cctk_dim,)&&\
        CCTK_DECLARE(INTEGER,cctk_gsh,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_lsh,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_lbnd,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_ubnd,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_ash,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_from,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_to,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_bbox,(2*cctk_dim))&&\
        CCTK_DECLARE(CCTK_REAL,cctk_delta_time,)&&\
        CCTK_DECLARE(CCTK_REAL,cctk_time,)&&\
        CCTK_DECLARE(CCTK_REAL,cctk_delta_space,(cctk_dim))&&\
        CCTK_DECLARE(CCTK_REAL,cctk_origin_space,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_levfac,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_levoff,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_levoffdenom,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_timefac,)&&\
        CCTK_DECLARE(INTEGER,cctk_convlevel,)&&\
        CCTK_DECLARE(INTEGER,cctk_convfac,)&&\
        CCTK_DECLARE(INTEGER,cctk_nghostzones,(cctk_dim))&&\
        CCTK_DECLARE(INTEGER,cctk_iteration,)&&\
        CCTK_DECLARE(CCTK_POINTER,cctkGH,)&&\
        CCTK_DECLARE(INTEGER,cctk_ash1,)&&\
        CCTK_DECLARE(INTEGER,cctk_ash2,)&&\
        CCTK_DECLARE(INTEGER,cctk_ash3,)&&

#define CCTK_WARN(a,b) CCTK_Warn(a,__LINE__,__FORTRANFILE__,CCTK_THORNSTRING,b)
#define CCTK_ERROR(b) CCTK_Error(__LINE__,__FORTRANFILE__,CCTK_THORNSTRING,b)

#define CCTK_CoordRegisterSystem(a,b,c) CCTKi_CoordRegisterSystem(a,b,CCTK_THORNSTRING,c)

/* John Shalf says that the operator .ne. needs to be enclosed by
   spaces, because ANSI C preprocessors otherwise interpret the
   character sequence ".0" as preprocessor token */
#define CCTK_EQUALS(a,b) (CCTK_Equals(a,b) .ne. 0)

#define CCTK_PASS_FTOF CCTK_FARGUMENTS

#define CCTK_ORIGIN_SPACE(x) (cctk_origin_space(x)+cctk_delta_space(x)/cctk_levfac(x)*cctk_levoff(x)/cctk_levoffdenom(x))
#define CCTK_DELTA_SPACE(x) (cctk_delta_space(x)/cctk_levfac(x))
#define CCTK_DELTA_TIME (cctk_delta_time/cctk_timefac)

#ifdef F90CODE

#define _DECLARE_CCTK_FUNCTIONS                   \
  external     CCTK_PointerTo                   &&\
  CCTK_POINTER CCTK_PointerTo                   &&\
  interface                                     &&\
     integer function CCTK_Equals (arg1, arg2)  &&\
       implicit none                            &&\
       CCTK_POINTER_TO_CONST arg1               &&\
       character(*) arg2                        &&\
     end function CCTK_Equals                   &&\
     integer function CCTK_MyProc (cctkGH)      &&\
       implicit none                            &&\
       CCTK_POINTER_TO_CONST cctkGH             &&\
     end function CCTK_MyProc                   &&\
     integer function CCTK_nProcs (cctkGH)      &&\
       implicit none                            &&\
       CCTK_POINTER_TO_CONST cctkGH             &&\
     end function CCTK_nProcs                   &&\
     integer function CCTK_IsThornActive (name) &&\
       implicit none                            &&\
       character(*) name                        &&\
     end function CCTK_IsThornActive            &&\
     CCTK_POINTER function CCTK_NullPointer ()  &&\
       implicit none                            &&\
     end function CCTK_NullPointer              &&\
  end interface                                 &&

#else /* ! F90CODE */

#define _DECLARE_CCTK_FUNCTIONS \
  integer      CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornActive &&\
  external     CCTK_Equals, CCTK_MyProc, CCTK_nProcs, CCTK_IsThornActive &&\
  CCTK_POINTER CCTK_PointerTo, CCTK_NullPointer &&\
  external     CCTK_PointerTo, CCTK_NullPointer &&

#endif /* ! F90CODE */

#endif /*FCODE*/

#ifdef CCODE

/* get the definition of ptrdiff_t */
#include <stddef.h>

#include "cGH.h"

#include "cctk_ActiveThorns.h"
#include "cctk_Banner.h"
#include "cctk_Coord.h"
#include "cctk_Comm.h"
#include "cctk_CommandLine.h"
#include "cctk_Complex.h"
#include "cctk_DebugDefines.h"
#include "cctk_Faces.h"
#include "cctk_File.h"
#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_Functions.h"
#include "cctk_GHExtensions.h"
#include "cctk_Groups.h"
#include "cctk_GroupsOnGH.h"
#include "cctk_Interp.h"
#include "cctk_IO.h"
#include "cctk_IOMethods.h"
#include "cctk_Loop.h"
#include "cctk_Main.h"
#include "cctk_Malloc.h"
#include "cctk_Math.h"
#include "cctk_Misc.h"
#include "cctk_Parameter.h"
#include "cctk_Reduction.h"
#include "cctk_Sync.h"
#include "cctk_Timers.h"
#include "cctk_Termination.h"
#include "cctk_WarnLevel.h"



/*
 * Routines to compute the linear index of a grid funtion element from
 * its i,j,k indices
 */

#ifdef CCTK_DEBUG

/*
 * For CCTK_DEBUG, call the external C routines defined in
 * DebugDefines.c
 */

#  define CCTK_GFINDEX0D(cctkGH) CCTK_GFIndex0D(cctkGH)
#  define CCTK_GFINDEX1D(cctkGH,i) CCTK_GFIndex1D(cctkGH,i)
#  define CCTK_GFINDEX2D(cctkGH,i,j) CCTK_GFIndex2D(cctkGH,i,j)
#  define CCTK_GFINDEX3D(cctkGH,i,j,k) CCTK_GFIndex3D(cctkGH,i,j,k)
#  define CCTK_GFINDEX4D(cctkGH,i,j,k,l) CCTK_GFIndex4D(cctkGH,i,j,k,l)

#  define CCTK_VECTGFINDEX0D(cctkGH,n) CCTK_VectGFIndex0D(cctkGH,n)
#  define CCTK_VECTGFINDEX1D(cctkGH,i,n) CCTK_VectGFIndex1D(cctkGH,i,n)
#  define CCTK_VECTGFINDEX2D(cctkGH,i,j,n) CCTK_VectGFIndex2D(cctkGH,i,j,n)
#  define CCTK_VECTGFINDEX3D(cctkGH,i,j,k,n) CCTK_VectGFIndex3D(cctkGH,i,j,k,n)
#  define CCTK_VECTGFINDEX4D(cctkGH,i,j,k,l,n) CCTK_VectGFIndex4D(cctkGH,i,j,k,l,n)

#else

/*
 * Without CCTK_DEBUG, define the calculations directly
 */

static inline int CCTK_GFINDEX0D (const cGH *restrict cctkGH)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_GFINDEX0D (const cGH *restrict cctkGH)
{
  return 0;
}

static inline int CCTK_GFINDEX1D (const cGH *restrict cctkGH,
                                  int i)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_GFINDEX1D (const cGH *restrict cctkGH,
                                  int i)
{
  return i;
}

static inline int CCTK_GFINDEX2D (const cGH *restrict cctkGH,
                                  int i, int j)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_GFINDEX2D (const cGH *restrict cctkGH,
                                  int i, int j)
{
  return i + cctkGH->cctk_ash[0] * j;
}

static inline int CCTK_GFINDEX3D (const cGH *restrict cctkGH,
                                  int i, int j, int k)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_GFINDEX3D (const cGH *restrict cctkGH,
                                  int i, int j, int k)
{
  return (i + cctkGH->cctk_ash[0] *
          (j + cctkGH->cctk_ash[1] * k));
}

static inline int CCTK_GFINDEX4D (const cGH *restrict cctkGH,
                                  int i, int j, int k, int l)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_GFINDEX4D (const cGH *restrict cctkGH,
                                  int i, int j, int k, int l)
{
  return (i + cctkGH->cctk_ash[0] *
          (j + cctkGH->cctk_ash[1] *
           (k + cctkGH->cctk_ash[2] * l)));
}

static inline int CCTK_VECTGFINDEX0D (const cGH *restrict cctkGH,
                                      int n)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_VECTGFINDEX0D (const cGH *restrict cctkGH,
                                      int n)
{
  return n;
}

static inline int CCTK_VECTGFINDEX1D (const cGH *restrict cctkGH,
                                      int i, int n)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_VECTGFINDEX1D (const cGH *restrict cctkGH,
                                      int i, int n)
{
  return i + cctkGH->cctk_ash[0] * n;
}

static inline int CCTK_VECTGFINDEX2D (const cGH *restrict cctkGH,
                                      int i, int j, int n)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_VECTGFINDEX2D (const cGH *restrict cctkGH,
                                      int i, int j, int n)
{
  return (i + cctkGH->cctk_ash[0] *
          (j + cctkGH->cctk_ash[1] * n));
}

static inline int CCTK_VECTGFINDEX3D (const cGH *restrict cctkGH,
                                      int i, int j, int k, int n)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_VECTGFINDEX3D (const cGH *restrict cctkGH,
                                      int i, int j, int k, int n)
{
  return (i + cctkGH->cctk_ash[0] *
          (j + cctkGH->cctk_ash[1] *
           (k + cctkGH->cctk_ash[2] * n)));
}

static inline int CCTK_VECTGFINDEX4D (const cGH *restrict cctkGH,
                                      int i, int j, int k, int l, int n)
  CCTK_ATTRIBUTE_NONNULL(1);
static inline int CCTK_VECTGFINDEX4D (const cGH *restrict cctkGH,
                                      int i, int j, int k, int l, int n)
{
  return (i + cctkGH->cctk_ash[0] *
          (j + cctkGH->cctk_ash[1] *
           (k + cctkGH->cctk_ash[2] *
            (l + cctkGH->cctk_ash[3] * n))));
}

#endif




#define CCTK_PRINTSEPARATOR \
  printf("--------------------------------------------------------------------------------\n");

#define _DECLARE_CCTK_ARGUMENTS _DECLARE_CCTK_CARGUMENTS
#define _DECLARE_CCTK_CARGUMENTS \
        CCTK_DECLARE_INIT(ptrdiff_t,cctki_dummy_int,0);\
        CCTK_DECLARE_INIT(int const,cctk_dim,cctkGH->cctk_dim);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_gsh,cctkGH->cctk_gsh);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_lsh,cctkGH->cctk_lsh);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_lbnd,cctkGH->cctk_lbnd);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_ubnd,cctkGH->cctk_ubnd);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_ash,cctkGH->cctk_ash);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_from,cctkGH->cctk_from);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_to,cctkGH->cctk_to);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_bbox,cctkGH->cctk_bbox);\
        CCTK_DECLARE_INIT(CCTK_REAL const,cctk_delta_time,cctkGH->cctk_delta_time);\
        CCTK_DECLARE_INIT(CCTK_REAL const,cctk_time,cctkGH->cctk_time);\
        CCTK_DECLARE_INIT(CCTK_REAL const *restrict const,cctk_delta_space,cctkGH->cctk_delta_space);\
        CCTK_DECLARE_INIT(CCTK_REAL const *restrict const,cctk_origin_space,cctkGH->cctk_origin_space);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_levfac,cctkGH->cctk_levfac);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_levoff,cctkGH->cctk_levoff);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_levoffdenom,cctkGH->cctk_levoffdenom);\
        CCTK_DECLARE_INIT(int const,cctk_timefac,cctkGH->cctk_timefac);\
        CCTK_DECLARE_INIT(int const,cctk_convlevel,cctkGH->cctk_convlevel);\
        CCTK_DECLARE_INIT(int const,cctk_convfac,cctkGH->cctk_convfac);\
        CCTK_DECLARE_INIT(int const *restrict const,cctk_nghostzones,cctkGH->cctk_nghostzones);\
        CCTK_DECLARE_INIT(int const,cctk_iteration,cctkGH->cctk_iteration);\

#define _INITIALISE_CCTK_C2F
#define _DECLARE_CCTK_C2F
#define _PASS_CCTK_C2F(xGH) &((xGH)->cctk_dim),\
                            (xGH)->cctk_gsh,(xGH)->cctk_lsh,\
                            (xGH)->cctk_lbnd,(xGH)->cctk_ubnd,\
                            (xGH)->cctk_ash,\
                            (xGH)->cctk_from,(xGH)->cctk_to,\
                            (xGH)->cctk_bbox,\
                            &((xGH)->cctk_delta_time), &((xGH)->cctk_time),\
                            (xGH)->cctk_delta_space, (xGH)->cctk_origin_space,\
                            (xGH)->cctk_levfac,\
                            (xGH)->cctk_levoff,\
                            (xGH)->cctk_levoffdenom,\
                            &((xGH)->cctk_timefac),\
                            &((xGH)->cctk_convlevel),\
                            &((xGH)->cctk_convfac),\
                            (xGH)->cctk_nghostzones,\
                            &((xGH)->cctk_iteration),\
                            &(xGH),\
                            &(xGH)->cctk_ash[0],\
                            &(xGH)->cctk_ash[1],\
                            &(xGH)->cctk_ash[2]
#define _CCTK_C2F_PROTO     int const *,\
                            int const *, int const *,\
                            int const *, int const *,\
                            int const *,\
                            int const *, int const *,\
                            int const *,\
                            CCTK_REAL const *, CCTK_REAL const *,\
                            CCTK_REAL const *, CCTK_REAL const *,\
                            int const *,\
                            int const *,\
                            int const *,\
                            int const *,\
                            int const *,\
                            int const *,\
                            int const *,\
                            int const *,\
                            cGH const *const *,\
                            int const *,\
                            int const *,\
                            int const *

#define CCTK_EQUALS(a,b) (CCTK_Equals((a),(b)))

#define CCTK_PASS_CTOC cctkGH

#define CCTK_ORIGIN_SPACE(x) (cctk_origin_space[x]+cctk_delta_space[x]/cctk_levfac[x]*cctk_levoff[x]/cctk_levoffdenom[x])
#define CCTK_DELTA_SPACE(x) (cctk_delta_space[x]/cctk_levfac[x])
#define CCTK_DELTA_TIME (cctk_delta_time/cctk_timefac)

#define CCTK_WARN(a,b) CCTK_Warn(a,__LINE__,__FILE__,CCTK_THORNSTRING,b)
#define CCTK_ERROR(b) CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,b)

#define CCTK_MALLOC(s) CCTKi_Malloc(s,__LINE__,__FILE__)
#define CCTK_FREE(p) CCTKi_Free(p)

#endif /*CCODE*/

#define CCTK_INFO(a) CCTK_Info(CCTK_THORNSTRING,(a))
#define CCTK_PARAMWARN(a) CCTK_ParamWarn(CCTK_THORNSTRING,(a))

#endif
