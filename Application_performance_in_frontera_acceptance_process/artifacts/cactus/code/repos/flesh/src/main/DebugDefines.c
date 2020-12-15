/*@@
   @file      DebugDefines.c
   @date      Tue 2 Jul 2001
   @author    Thomas Radke
   @desc
              Routines to provide some debugging support for the Cactus code.
   @enddesc
   @version   $Id$
 @@*/

#include "cctk_Config.h"
#include "cctk_Flesh.h"
#include "cctk_DebugDefines.h"
#include "cctk_WarnLevel.h"
#include "definethisthorn.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_DebugDefines_c);


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTK_GFIndex?D
   @date       Tue 2 Jul 2001
   @author     Thomas Radke
   @desc
               Compute the linear index of a grid function element
               from its spatial indices
   @enddesc

   @var        GH
   @vdesc      pointer to CCTK grid hierarchy extension
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        i, j, k, l
   @vdesc      spatial indices
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the linear index for the given spatial indices
   @endreturndesc
@@*/
int CCTK_GFIndex0D (const cGH *GH)
{
  return (0);
}

int CCTK_GFIndex1D (const cGH *GH, int i)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d cctk_lsh=[%d]",
                 i, GH->cctk_lsh[0]);
  }
#endif
  return (i);
}

int CCTK_GFIndex2D (const cGH *GH, int i, int j)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0] ||
      j < 0 || j >= GH->cctk_lsh[1])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d j=%d cctk_lsh=[%d,%d]",
                 i, j, GH->cctk_lsh[0], GH->cctk_lsh[1]);
  }
#endif
  return (i + GH->cctk_ash[0]*j);
}

int CCTK_GFIndex3D (const cGH *GH, int i, int j, int k)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0] ||
      j < 0 || j >= GH->cctk_lsh[1] ||
      k < 0 || k >= GH->cctk_lsh[2])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d j=%d k=%d cctk_lsh=[%d,%d,%d]",
                 i, j, k, GH->cctk_lsh[0], GH->cctk_lsh[1], GH->cctk_lsh[2]);
  }
#endif
  return (i + GH->cctk_ash[0]*(j + GH->cctk_ash[1]*k));
}

int CCTK_GFIndex4D (const cGH *GH, int i, int j, int k, int l)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0] ||
      j < 0 || j >= GH->cctk_lsh[1] ||
      k < 0 || k >= GH->cctk_lsh[2] ||
      l < 0 || l >= GH->cctk_lsh[3])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d j=%d k=%d l=%d cctk_lsh=[%d,%d,%d,%d]",
                 i, j, k, l, GH->cctk_lsh[0], GH->cctk_lsh[1], GH->cctk_lsh[2], GH->cctk_lsh[3]);
  }
#endif
  return (i + GH->cctk_ash[0]*(j + GH->cctk_ash[1]*(k + GH->cctk_ash[2] * l)));
}

int CCTK_VectGFIndex0D (const cGH *GH, int n)
{
#ifdef CCTK_DEBUG
  if (n < 0)
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Vector index out of bounds.  n=%d",
                 n);
  }
#endif
  return (n);
}

int CCTK_VectGFIndex1D (const cGH *GH, int i, int n)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d cctk_lsh=[%d]",
                 i, GH->cctk_lsh[0]);
  }
  if (n < 0)
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Vector index out of bounds.  n=%d",
                 n);
  }
#endif
  return (i + GH->cctk_ash[0]*n);
}

int CCTK_VectGFIndex2D (const cGH *GH, int i, int j, int n)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0] ||
      j < 0 || j >= GH->cctk_lsh[1])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d j=%d cctk_lsh=[%d,%d]",
                 i, j, GH->cctk_lsh[0], GH->cctk_lsh[1]);
  }
  if (n < 0)
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Vector index out of bounds.  n=%d",
                 n);
  }
#endif
  return (i + GH->cctk_ash[0]*(j + GH->cctk_ash[1]*n));
}

int CCTK_VectGFIndex3D (const cGH *GH, int i, int j, int k, int n)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0] ||
      j < 0 || j >= GH->cctk_lsh[1] ||
      k < 0 || k >= GH->cctk_lsh[2])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d j=%d k=%d cctk_lsh=[%d,%d,%d]",
                 i, j, k, GH->cctk_lsh[0], GH->cctk_lsh[1], GH->cctk_lsh[2]);
  }
  if (n < 0)
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Vector index out of bounds.  n=%d",
                 n);
  }
#endif
  return (i + GH->cctk_ash[0]*(j + GH->cctk_ash[1]*(k + GH->cctk_ash[2]*n)));
}

int CCTK_VectGFIndex4D (const cGH *GH, int i, int j, int k, int l, int n)
{
#ifdef CCTK_DEBUG
  if (i < 0 || i >= GH->cctk_lsh[0] ||
      j < 0 || j >= GH->cctk_lsh[1] ||
      k < 0 || k >= GH->cctk_lsh[2] ||
      l < 0 || l >= GH->cctk_lsh[3])
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid function index out of bounds.  i=%d j=%d k=%d l=%d cctk_lsh=[%d,%d,%d,%d]",
                 i, j, k, l, GH->cctk_lsh[0], GH->cctk_lsh[1], GH->cctk_lsh[2], GH->cctk_lsh[3]);
  }
  if (n < 0)
  {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Vector index out of bounds.  n=%d",
                 n);
  }
#endif
  return (i + GH->cctk_ash[0]*(j + GH->cctk_ash[1]*(k + GH->cctk_ash[2]*
                                                    (l + GH->cctk_ash[3]*n))));
}
