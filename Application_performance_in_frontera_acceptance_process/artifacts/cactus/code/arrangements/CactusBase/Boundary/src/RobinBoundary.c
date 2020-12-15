/*@@
  @file      RobinBoundary.c
  @date      July 6th 2000
  @author    Miguel Alcubierre, Gabrielle Allen, Gerd Lanfermann
  @desc
             Routines for Robin boundary conditions
  @enddesc
  @history
  @hdate     Tue 10 Apr 2001
  @hauthor   Thomas Radke
  @hdesc     BC routines generalized for applying to arbitrary CCTK data types
  @endhistory
  @version   $Id$
@@*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"

#include "Boundary.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_Boundary_RobinBoundary_c);

static int ApplyBndRobin(const cGH *GH, const CCTK_INT *stencil, CCTK_REAL finf,
                         int npow, int first_var, int num_vars);
static int OldApplyBndRobin(const cGH *GH, const int *stencil, CCTK_REAL finf,
                            int npow, int first_var, int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/*@@
   @routine    BndRobin
   @date       14 Feb 2003
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               this boundary condition
   @enddesc
   @calls      ApplyBndRobin
   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables passed in through var_indices[]
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        var_indices
   @vdesc      array of variable indicies to which to apply this boundary
               condition
   @vtype      CCTK_INT *
   @vio        in
   @endvar
   @var        faces
   @vdesc      array of set of faces to which to apply the bc
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        widths
   @vdesc      array of boundary widths for each variable
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        table_handles
   @vdesc      array of table handles which hold extra arguments
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @returntype CCTK_INT
   @returndesc
               return code of @seeroutine ApplyBndRobin
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/

CCTK_INT BndRobin(const cGH *GH, CCTK_INT num_vars, CCTK_INT *vars,
                  CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *tables) {
  int i, j, k, gi, err, gdim, max_gdim, retval;

  /* variables to pass to ApplyBndRobin */
  CCTK_INT *width_alldirs; /* width of boundary in all directions */
  CCTK_REAL finf;          /* value of function at infinity */
  CCTK_INT npow;           /* decay rate */

#ifdef DEBUG
  printf(
      "BndRobin(): got passed GH=%p, num_vars=%d, vars[0]=%d, tables[0]=%d\n",
      (const void *)GH, num_vars, vars[0], tables[0]);
#endif

  retval = 0;
  width_alldirs = NULL;
  max_gdim = 0;

  /* loop through variables, j at a time */
  for (i = 0; i < num_vars; i += j) {
    /* find other adjacent vars which are selected for identical bcs */
    j = 1;
    /* Since GFs are allowed to have different staggering, the best we
       can do is find variables of the same group which are selected
       for identical bcs.  If all GFs had the same staggering then we
       could groups many GFs together. */
    gi = CCTK_GroupIndexFromVarI(vars[i]);
    while (i + j < num_vars && vars[i + j] == vars[i] + j &&
           CCTK_GroupIndexFromVarI(vars[i + j]) == gi &&
           tables[i + j] == tables[i] && faces[i + j] == faces[i] &&
           widths[i + j] == widths[i]) {
      ++j;
    }

    /* Check to see if faces specification is valid */
    if (faces[i] != CCTK_ALL_FACES) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Faces specification %d for Robin boundary conditions on "
                 "%s is not implemented yet.  "
                 "Applying Robin bcs to all (external) faces.",
                 (int)faces[i], CCTK_VarName(vars[i]));
    }

    /* Set up default arguments for ApplyBndRobin */
    finf = 0;
    npow = 1;

    /* Look on table for possible non-default arguments
     * (If any of these table look-ups fail, the value will be unchanged
     * from its default value)
     */
    /* Asymptotic value of function at infinity */
    err = Util_TableGetReal(tables[i], &finf, "FINF");
    if (err == UTIL_ERROR_BAD_HANDLE) {
      CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Invalid table handle passed for Robin boundary "
                 "conditions for %s.  Using all default values.",
                 CCTK_VarName(vars[i]));
    } else {
      /* Decay power */
      Util_TableGetInt(tables[i], &npow, "DECAY_POWER");
    }

    /* Determine boundary width on all faces */
    /* allocate memory for buffer */
    gdim = CCTK_GroupDimI(gi);
    if (gdim > max_gdim) {
      width_alldirs =
          (CCTK_INT *)realloc(width_alldirs, 2 * gdim * sizeof(CCTK_INT));
      max_gdim = gdim;
    }

    /* fill it with values, either from table or the boundary_width
       parameter */
    if (widths[i] < 0) {
      err = Util_TableGetIntArray(tables[i], gdim, width_alldirs,
                                  "BOUNDARY_WIDTH");
      if (err < 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Error %d when reading boundary width array from table "
                   "for %s",
                   err, CCTK_VarName(vars[i]));
        return -21;
      } else if (err != 2 * gdim) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Boundary width array for %s has %d elements, but %d "
                   "expected",
                   CCTK_VarName(vars[i]), err, 2 * gdim);
        return -22;
      }
    } else {
      for (k = 0; k < 2 * gdim; ++k) {
        width_alldirs[k] = widths[i];
      }
    }

    /* Apply the boundary condition */
    if ((retval = ApplyBndRobin(GH, width_alldirs, finf, npow, vars[i], j)) <
        0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndRobin() returned %d", retval);
    }
  }
#ifdef DEBUG
  printf("BndRobin(): returning %d\n", retval);
#endif
  free(width_alldirs);

  return retval;
}

/* prototypes for external C routines are declared in header Boundary.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL CCTK_FNAME(BndRobinGI)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       const int *gi);
void CCTK_FCALL CCTK_FNAME(BndRobinGN)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndRobinVI)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       const int *vi);
void CCTK_FCALL CCTK_FNAME(BndRobinVN)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       ONE_FORTSTRING_ARG);

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndRobinGI
   @date       Tue Jul 18 18:08:28 2000
   @author     Gerd Lanfermann
   @desc
               Apply Robin boundary conditions by group index
   @enddesc
   @calls      ApplyBndRobin

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width array
   @vtype      int [ dimension of group ]
   @vio        in
   @endvar
   @var        finf
   @vdesc      value of f at infimum
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        npow
   @vdesc      power of decay rate
   @vtype      int
   @vio        in
   @endvar
   @var        gi
   @vdesc      index of group to apply boundary conditions to
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndRobin
               -1 if invalid group index was given
   @endreturndesc
@@*/
int BndRobinGI(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               int gi) {
  int first_vi, retval;

  first_vi = CCTK_FirstVarIndexI(gi);
  if (first_vi >= 0) {
    retval = OldApplyBndRobin(GH, stencil, finf, npow, first_vi,
                              CCTK_NumVarsInGroupI(gi));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group index %d in BndFlatGI", gi);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndRobinGI)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       const int *gi) {
  *ierr = BndRobinGI(*GH, stencil, *finf, *npow, *gi);
}

/*@@
   @routine    BndRobinGN
   @date       Tue Jul 18 18:08:28 2000
   @author     Gerd Lanfermann
   @desc
               Apply Robin boundary conditions by group name
   @enddesc
   @calls      BndRobinGI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width array
   @vtype      int [ dimension of group ]
   @vio        in
   @endvar
   @var        finf
   @vdesc      value of f at infimum
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        npow
   @vdesc      power of decay rate
   @vtype      int
   @vio        in
   @endvar
   @var        gname
   @vdesc      name of group to apply boundary conditions to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndRobinGI <BR>
               -1 if invalid group index was given
   @endreturndesc
@@*/
int BndRobinGN(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               const char *gname) {
  int gi, retval;

  gi = CCTK_GroupIndex(gname);
  if (gi >= 0) {
    retval = BndRobinGI(GH, stencil, finf, npow, gi);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group name '%s' in BndRobinGN", gname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndRobinGN)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(gname)
  *ierr = BndRobinGN(*GH, stencil, *finf, *npow, gname);
  free(gname);
}

/*@@
   @routine    BndRobinVI
   @date       Tue Jul 18 18:08:28 2000
   @author     Gerd Lanfermann
   @desc
               Apply Robin boundary conditions by variable index
   @enddesc
   @calls      ApplyBndRobin

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width array
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        finf
   @vdesc      value of f at infimum
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        npow
   @vdesc      power of decay rate
   @vtype      int
   @vio        in
   @endvar
   @var        vi
   @vdesc      index of variable to apply boundary conditions to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndRobin <BR>
               -1 if invalid group index was given
   @endreturndesc
@@*/
int BndRobinVI(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               int vi) {
  int retval;

  if (vi >= 0 && vi < CCTK_NumVars()) {
    retval = OldApplyBndRobin(GH, stencil, finf, npow, vi, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "BndRobinVI: Invalid variable index %d", vi);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndRobinVI)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       const int *vi) {
  *ierr = BndRobinVI(*GH, stencil, *finf, *npow, *vi);
}

/*@@
   @routine    BndRobinVN
   @date       Tue Jul 18 18:08:28 2000
   @author     Gerd Lanfermann
   @desc
               Apply Robin boundary conditions by variable name
   @enddesc
   @calls      BndRobinVI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width array
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        finf
   @vdesc      value of f at infimum
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        npow
   @vdesc      power of decay rate
   @vtype      int
   @vio        in
   @endvar
   @var        vname
   @vdesc      name of variable to apply boundary conditions to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndRobinVI <BR>
               -1 if invalid group index was given
   @endreturndesc
@@*/
int BndRobinVN(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               const char *vname) {
  int vi, retval;

  vi = CCTK_VarIndex(vname);
  if (vi >= 0) {
    retval = BndRobinVI(GH, stencil, finf, npow, vi);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable name '%s' in BndRobinVN", vname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndRobinVN)(int *ierr, const cGH **GH,
                                       const int *stencil,
                                       const CCTK_REAL *finf, const int *npow,
                                       ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(vname)
  *ierr = BndRobinVN(*GH, stencil, *finf, *npow, vname);
  free(vname);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/* maximum dimension we can deal with */
#define MAXDIM 3

/* macro to compute x*x */
#define SQR(x) ((x) * (x))

/*@@
   @routine    SET_LINEAR_INDICES
   @date       Thu 7 June 2001
   @author     Thomas Radke
   @desc
               Macro to set the linear indices for the source and destination
               element of the current grid variable
   @enddesc

   @var        i
   @vdesc      x index to use
   @vtype      int
   @vio        in
   @endvar
@@*/
#define SET_LINEAR_INDICES(i)                                                  \
  {                                                                            \
    dst = CCTK_GFINDEX3D(GH, i, j, k);                                         \
    src = CCTK_GFINDEX3D(GH, (i) + dx, j + dy, k + dz);                        \
    distance = dist[abs(dx) + 2 * abs(dy) + 4 * abs(dz)];                      \
  }

/*@@
   @routine    ROBIN_BOUNDARY_TYPED_3D
   @date       Thu 7 June 2001
   @author     Thomas Radke
   @desc
               Macro to apply Robin boundary conditions to a 3D variable
               of given datatype
   @enddesc

   @var        cctk_type
   @vdesc      CCTK datatype of the variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define ROBIN_BOUNDARY_TYPED_3D(cctk_type)                                     \
  {                                                                            \
    cctk_type *data;                                                           \
    double u_src, u_dst, aux;                                                  \
                                                                               \
    /* avoid the else branch with the expensive sqrt() operation if possible   \
     */                                                                        \
    if (abs(dx) + abs(dy) + abs(dz) == 1) {                                    \
      u_dst = fabs((double)(dx ? x[dst] : (dy ? y[dst] : z[dst])));            \
      u_src = fabs((double)(dx ? x[src] : (dy ? y[src] : z[src])));            \
    } else {                                                                   \
      u_dst = sqrt(SQR(dx * x[dst]) + SQR(dy * y[dst]) + SQR(dz * z[dst]));    \
      u_src = sqrt(SQR(dx * x[src]) + SQR(dy * y[src]) + SQR(dz * z[src]));    \
    }                                                                          \
                                                                               \
    aux = decay * distance * (u_src + u_dst) / SQR(r[src] + r[dst]);           \
                                                                               \
    data = (cctk_type *)GH->data[var][0];                                      \
    data[dst] =                                                                \
        (cctk_type)((2 * aux * finf + data[src] * (1 - aux)) / (1 + aux));     \
  }

/*@@
   @routine    ROBIN_BOUNDARY
   @date       Thu 7 June 2001
   @author     Thomas Radke
   @desc
               Macro to apply Robin boundary conditions to a variable
               of a given datatype in all directions
               Currently it is limited up to 3D variables only.
   @enddesc
   @calls      SET_LINEAR_INDICES
               ROBIN_BOUNDARY_TYPED_3D

   @var        cctk_type
   @vdesc      CCTK datatype of the variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define ROBIN_BOUNDARY(cctk_type)                                              \
  {                                                                            \
    int i, j, k;                                                               \
    int dx, dy, dz;                                                            \
    int src, dst;                                                              \
    double distance;                                                           \
                                                                               \
    /* check the dimensionality */                                             \
    if (gdim != 3) {                                                           \
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,                      \
                 "ApplyBndRobin: variable dimension of %d not supported",      \
                 gdim);                                                        \
      return (-5);                                                             \
    }                                                                          \
                                                                               \
    if (in_widths[0] == 2 || in_widths[1] == 2 || in_widths[2] == 2) {         \
      /* outermost loop over almost all z points */                            \
      for (k = 1; k < GH->cctk_lsh[2] - 1; k++) {                              \
        dz = 0;                                                                \
        if (k == 1 && doBC[4]) {                                               \
          dz = +1;                                                             \
        } else if (k == GH->cctk_lsh[2] - 2 && doBC[5]) {                      \
          dz = -1;                                                             \
        }                                                                      \
                                                                               \
        /* middle loop over all y points */                                    \
        for (j = 1; j < GH->cctk_lsh[1] - 1; j++) {                            \
          dy = 0;                                                              \
          if (j == 1 && doBC[2]) {                                             \
            dy = +1;                                                           \
          } else if (j == GH->cctk_lsh[1] - 2 && doBC[3]) {                    \
            dy = -1;                                                           \
          }                                                                    \
                                                                               \
          /* lower x */                                                        \
          dx = 0;                                                              \
          if (doBC[0]) {                                                       \
            dx = +1;                                                           \
          }                                                                    \
          if (dx || dy || dz) {                                                \
            SET_LINEAR_INDICES(1);                                             \
            ROBIN_BOUNDARY_TYPED_3D(cctk_type);                                \
          }                                                                    \
                                                                               \
          /* lower/upper y and/or z */                                         \
          if (dy || dz) {                                                      \
            dx = 0;                                                            \
            SET_LINEAR_INDICES(2);                                             \
            for (i = 2; i < GH->cctk_lsh[0] - 2; i++, src++, dst++) {          \
              ROBIN_BOUNDARY_TYPED_3D(cctk_type);                              \
            }                                                                  \
          }                                                                    \
                                                                               \
          /* upper x */                                                        \
          dx = 0;                                                              \
          if (doBC[1]) {                                                       \
            dx = -1;                                                           \
          }                                                                    \
          if (dx || dy || dz) {                                                \
            SET_LINEAR_INDICES(GH->cctk_lsh[0] - 2);                           \
            ROBIN_BOUNDARY_TYPED_3D(cctk_type);                                \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
                                                                               \
    /* outermost loop over all z points */                                     \
    for (k = 0; k < GH->cctk_lsh[2]; k++) {                                    \
      dz = 0;                                                                  \
      if (k == 0 && doBC[4]) {                                                 \
        dz = +1;                                                               \
      } else if (k == GH->cctk_lsh[2] - 1 && doBC[5]) {                        \
        dz = -1;                                                               \
      }                                                                        \
                                                                               \
      /* middle loop over all y points */                                      \
      for (j = 0; j < GH->cctk_lsh[1]; j++) {                                  \
        dy = 0;                                                                \
        if (j == 0 && doBC[2]) {                                               \
          dy = +1;                                                             \
        } else if (j == GH->cctk_lsh[1] - 1 && doBC[3]) {                      \
          dy = -1;                                                             \
        }                                                                      \
                                                                               \
        /* lower x */                                                          \
        dx = 0;                                                                \
        if (doBC[0]) {                                                         \
          dx = +1;                                                             \
        }                                                                      \
        if (dx || dy || dz) {                                                  \
          SET_LINEAR_INDICES(0);                                               \
          ROBIN_BOUNDARY_TYPED_3D(cctk_type);                                  \
        }                                                                      \
                                                                               \
        /* lower/upper y and/or z */                                           \
        if (dy || dz) {                                                        \
          dx = 0;                                                              \
          SET_LINEAR_INDICES(1);                                               \
          for (i = 1; i < GH->cctk_lsh[0] - 1; i++, src++, dst++) {            \
            ROBIN_BOUNDARY_TYPED_3D(cctk_type);                                \
          }                                                                    \
        }                                                                      \
                                                                               \
        /* upper x */                                                          \
        dx = 0;                                                                \
        if (doBC[1]) {                                                         \
          dx = -1;                                                             \
        }                                                                      \
        if (dx || dy || dz) {                                                  \
          SET_LINEAR_INDICES(GH->cctk_lsh[0] - 1);                             \
          ROBIN_BOUNDARY_TYPED_3D(cctk_type);                                  \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    ApplyBndRobin
   @date       Tue Jul 18 18:08:28 2000
   @author     Gerd Lanfermann
   @desc
               Apply Robin boundary conditions to a group of grid functions
               given by their indices
               This routine is called by the various BndRobinXXX wrappers.

               Although it is currently limited to handle 3D variables only
               it can easily be extended for higher dimensions
               by adapting the appropriate macros.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        in_widths
   @vdesc      boundary width array
   @vtype      CCTK_INT [ dimension of variable ]
   @vio        in
   @endvar
   @var        finf
   @vdesc      value of f at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        npow
   @vdesc      power of decay rate
   @vtype      int
   @vio        in
   @endvar
   @var        first_var
   @vdesc      index of first variable to apply boundary conditions to
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_VarTypeI
               CCTK_GroupDimFromVarI
               ROBIN_BOUNDARY
   @history
   @hdate      Tue 10 Apr 2001
   @hauthor    Thomas Radke
   @hdesc      Merged separate routines for 1D, 2D, and 3D
               into a single generic routine
   @endhistory

   @returntype int
   @returndesc
                0 for success
               -1 if variable dimension is not supported
               -2 if NULL pointer passed as boundary width array
               -3 if stencil width is other than 1
               -4 if variable type is not supported
               -5 if variable dimension is other than 3D
               -6 if no coordinate information is available
   @endreturndesc
@@*/
static int ApplyBndRobin(const cGH *GH, const CCTK_INT *in_widths,
                         CCTK_REAL finf, int npow, int first_var,
                         int num_vars) {
  int var, vtype, dim, gdim;
  int doBC[2 * MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];
  CCTK_INT ierr;
  char coord_system_name[20];
  double decay;
  const CCTK_REAL *x, *y, *z, *r;
  double dist[8];

  /* get the number of dimensions and the variables' type */
  gdim = CCTK_GroupDimI(CCTK_GroupIndexFromVarI(first_var));
  vtype = CCTK_VarTypeI(first_var);

  /* make sure we can deal with this number of dimensions */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndRobin: Variable dimension of %d not supported", gdim);
    return (-1);
  }

  /* check the boundary width */
  if (!in_widths) {
    CCTK_WARN(1, "ApplyBndRobin: NULL pointer passed for boundary width "
                 "array");
    return (-2);
  }

  for (dim = 0; dim < 2 * gdim; dim++) {
    if (in_widths[dim] != 1 && in_widths[dim] != 2) {
      CCTK_WARN(1, "ApplyBndRobin: Stencil width must be 1 or 2 "
                   "for Robin boundary conditions");
      return (-3);
    }
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths(GH, first_var, gdim, in_widths, "Robin");

  /* Robin boundaries need the underlying grid coordinates */
  sprintf(coord_system_name, "cart%dd", gdim);
  if (CCTK_CoordSystemHandle(coord_system_name) < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndRobin: Couldn't get coordinates from '%s'",
               coord_system_name);
    return (-6);
  }
  x = GH->data[CCTK_CoordIndex(-1, "x", coord_system_name)][0];
  y = GH->data[CCTK_CoordIndex(-1, "y", coord_system_name)][0];
  z = GH->data[CCTK_CoordIndex(-1, "z", coord_system_name)][0];

  sprintf(coord_system_name, "spher%dd", gdim);
  if (CCTK_CoordSystemHandle(coord_system_name) < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndRobin: Couldn't get coordinates from '%s'",
               coord_system_name);
    return (-6);
  }
  r = GH->data[CCTK_CoordIndex(-1, "r", coord_system_name)][0];

  /* see if we have a physical boundary */
  symtable = SymmetryTableHandleForGrid(GH);
  if (symtable < 0)
    CCTK_WARN(0, "internal error");
  ierr = Util_TableGetIntArray(symtable, 2 * gdim, symbnd, "symmetry_handle");
  if (ierr != 2 * gdim)
    CCTK_WARN(0, "internal error");
  for (dim = 0; dim < 2 * gdim; dim++) {
    is_physical[dim] = symbnd[dim] < 0;
  }

  /* get the decay rate as a double */
  decay = (double)npow;

  /* precompute the distance to all 8 neighbors in a 3D grid */
  dist[0] = 0; /* not used */
  dist[1] = GH->cctk_delta_space[0] / GH->cctk_levfac[0];
  dist[2] = GH->cctk_delta_space[1] / GH->cctk_levfac[1];
  dist[3] = sqrt(SQR(dist[1]) + SQR(dist[2]));
  dist[4] = GH->cctk_delta_space[2] / GH->cctk_levfac[2];
  dist[5] = sqrt(SQR(dist[1]) + SQR(dist[4]));
  dist[6] = sqrt(SQR(dist[2]) + SQR(dist[4]));
  dist[7] = sqrt(SQR(dist[1]) + SQR(dist[2]) + SQR(dist[4]));

  /* now loop over all variables */
  for (var = first_var; var < first_var + num_vars; var++) {
    /* Apply condition if:
       + boundary is a physical boundary
       + boundary is an outer boundary
       + have enough grid points
    */
    for (dim = 0; dim < 2 * gdim; dim++) {
      doBC[dim] = is_physical[dim];
    }
    for (dim = 0; dim < gdim; dim++) {
      doBC[dim * 2] &= GH->cctk_lsh[dim] > 1 && GH->cctk_bbox[dim * 2];
      doBC[dim * 2 + 1] &= GH->cctk_lsh[dim] > 1 && GH->cctk_bbox[dim * 2 + 1];
    }

    switch (vtype) {
    case CCTK_VARIABLE_REAL:
      ROBIN_BOUNDARY(CCTK_REAL);
      break;

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      ROBIN_BOUNDARY(CCTK_REAL4);
      break;
#endif

#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      ROBIN_BOUNDARY(CCTK_REAL8);
      break;
#endif

#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      ROBIN_BOUNDARY(CCTK_REAL16);
      break;
#endif

    default:
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndRobin: Unsupported variable type %d for "
                 "variable '%s'",
                 CCTK_VarTypeI(var), CCTK_VarName(var));
      return (-4);
    }
  }

  return (0);
}

/*@@
   @routine    OldApplyBndRobin
   @date       5 May 2003
   @author     David Rideout
   @desc
               The new boundary API expects a 2d-element array for the
               boundary_widths (d=dimension of grid variable), while
               the old API expects a d-element array.  This function
               converts the old array to the new format.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        in_widths
   @vdesc      boundary width array
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        finf
   @vdesc      value of f at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        npow
   @vdesc      power of decay rate
   @vtype      int
   @vio        in
   @endvar
   @var        first_var
   @vdesc      index of first variable to apply boundary conditions to
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupIndexFromVarI
               ApplyBndRobin
   @returntype int
   @returndesc
               returncode from @seeroutine ApplyBndRobin
   @endreturndesc
@@*/
static int OldApplyBndRobin(const cGH *GH, const int *in_widths, CCTK_REAL finf,
                            int npow, int first_var, int num_vars) {
  int i, dim, retval;
  CCTK_INT *boundary_widths;
  static int warned;

  /* Convert stencil_alldirs to new format */
  dim = CCTK_GroupDimFromVarI(first_var);
  boundary_widths = malloc(2 * dim * sizeof(CCTK_INT));
  for (i = 0; i < 2 * dim; ++i) {
    boundary_widths[i] = in_widths[i / 2];
  }

  /* Bug people for using the old interface */
  if (!warned) {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Copied older d-element array of boundary widths into the "
               "newer 2d-element format.  Please use the new boundary "
               "interface to avoid this.");
    warned = 1;
  }

  /* Call ApplyBnd... with new boundary width array */
  retval = ApplyBndRobin(GH, boundary_widths, finf, npow, first_var, num_vars);

  free(boundary_widths);
  return retval;
}
