/*@@
  @file      StaticBoundary.c
  @date      Sat Mar 16 15:09:00 2001
  @author    Gabrielle Allen
  @desc
             Routines for applying static-boundary conditions
  @enddesc
  @version   $Id$
@@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include "cctk_FortranString.h"

#include "Boundary.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_Boundary_StaticBoundary_c);

static int ApplyBndStatic(const cGH *GH, CCTK_INT stencil_dir,
                          const CCTK_INT *stencil_alldirs, int dir,
                          int first_var, int num_vars);
static int OldApplyBndStatic(const cGH *GH, int stencil_dir,
                             const int *stencil_alldirs, int dir, int first_var,
                             int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/*@@
   @routine    BndStatic
   @date       14 Feb 2003
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               this boundary condition
   @enddesc
   @calls      ApplyBndStatic
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
               return code of @seeroutine ApplyBndStatic
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/
CCTK_INT BndStatic(const cGH *GH, CCTK_INT num_vars, CCTK_INT *vars,
                   CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *tables) {
  int i, j, k, gi, err, gdim, max_gdim, retval;

  /* variables to pass to ApplyBndStatic */
  CCTK_INT *width_alldirs; /* width of boundary in all directions */
  int dir;                 /* direction in which to apply bc */

#ifdef DEBUG
  printf(
      "BndStatic(): got passed GH=%p, num_vars=%d, vars[0]=%d, tables[0]=%d\n",
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
                 "Faces specification %d for Static boundary conditions on "
                 "%s is not implemented yet.  "
                 "Applying Static bcs to all (external) faces.",
                 (int)faces[i], CCTK_FullName(vars[i]));
    }
    dir = 0;

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
                   err, CCTK_FullName(vars[i]));
        return -21;
      } else if (err != 2 * gdim) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Boundary width array for %s has %d elements, but %d "
                   "expected",
                   CCTK_FullName(vars[i]), err, 2 * gdim);
        return -22;
      }
    } else {
      for (k = 0; k < 2 * gdim; ++k) {
        width_alldirs[k] = widths[i];
      }
    }

    /* Apply the boundary condition */
    if ((retval = ApplyBndStatic(GH, 0, width_alldirs, dir, vars[i], j)) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndStatic() returned %d", retval);
    }
  }
#ifdef DEBUG
  printf("BndStatic(): returning %d\n", retval);
#endif
  free(width_alldirs);

  return retval;
}

/* prototypes for external C routines are declared in header Boundary.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL CCTK_FNAME(BndStaticDirVI)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, const int *vi);
void CCTK_FCALL CCTK_FNAME(BndStaticVI)(int *ierr, const cGH **GH,
                                        const int *stencil, const int *vi);
void CCTK_FCALL CCTK_FNAME(BndStaticDirGI)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, const int *gi);
void CCTK_FCALL CCTK_FNAME(BndStaticGI)(int *ierr, const cGH **GH,
                                        const int *stencil, const int *gi);
void CCTK_FCALL CCTK_FNAME(BndStaticDirVN)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndStaticVN)(int *ierr, const cGH **GH,
                                        const int *stencil, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndStaticDirGN)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndStaticGN)(int *ierr, const cGH **GH,
                                        const int *stencil, ONE_FORTSTRING_ARG);

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndStaticDirVI
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply static boundary routines by var index in given direction
   @enddesc
   @calls      ApplyBndStatic

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply boundary
   @vtype      int
   @vio        in
   @endvar
   @var        vi
   @vdesc      index of variable to apply static boundaries
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndStatic <BR>
               -1 if invalid variable indices are given
   @endreturndesc
@@*/
int BndStaticDirVI(const cGH *GH, int stencil_size, int dir, int vi) {
  int retval, num_vars;

  num_vars = CCTK_NumVars();
  if (vi >= 0 && vi < num_vars) {
    retval = ApplyBndStatic(GH, stencil_size, NULL, dir, vi, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable index %d in BndStaticDirVI", vi);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticDirVI)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, const int *vi) {
  *ierr = BndStaticDirVI(*GH, *stencil_size, *dir, *vi);
}

/*@@
   @routine    BndStaticVI
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply static boundary routines by var index
   @enddesc
   @calls      ApplyBndStatic

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        vi
   @vdesc      index of variable to apply static boundaries to
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndStatic <BR>
               -1 if invalid variable indices are given
   @endreturndesc
@@*/
int BndStaticVI(const cGH *GH, const int *stencil, int vi) {
  int retval, num_vars;

  num_vars = CCTK_NumVars();
  if (vi >= 0 && vi < num_vars) {
    retval = OldApplyBndStatic(GH, -1, stencil, 0, vi, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable index %d in BndStaticVI", vi);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticVI)(int *ierr, const cGH **GH,
                                        const int *stencil, const int *vi) {
  *ierr = BndStaticVI(*GH, stencil, *vi);
}

/* ====================================================== */

/*@@
  @routine    BndStaticDirGI
  @date       Sat Jan 20 2001
  @author     Gabrielle Allen
  @desc
              Apply static boundaries by group index in given direction
  @enddesc
  @calls      ApplyBndStatic

  @var        GH
  @vdesc      Pointer to CCTK grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        stencil_size
  @vdesc      stencil size in this direction
  @vtype      int
  @vio        in
  @endvar
  @var        dir
  @vdesc      direction to apply static boundaries
  @vtype      int
  @vio        in
  @endvar
  @var        gi
  @vdesc      index of group to apply static boundaries to
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              return code of @seeroutine ApplyBndStatic <BR>
              -1 if invalid group indices are given
  @endreturndesc
@@*/
int BndStaticDirGI(const cGH *GH, int stencil_size, int dir, int gi) {
  int first_vi, retval;

  first_vi = CCTK_FirstVarIndexI(gi);
  if (first_vi >= 0) {
    retval = ApplyBndStatic(GH, stencil_size, NULL, dir, first_vi,
                            CCTK_NumVarsInGroupI(gi));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group indices %d in BndStaticDirGI", gi);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticDirGI)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, const int *gi) {
  *ierr = BndStaticDirGI(*GH, *stencil_size, *dir, *gi);
}

/*@@
  @routine    BndStaticGI
  @date       Thu Mar  2 11:07:11 2000
  @author     Gerd Lanfermann
  @desc
              Apply static boundaries by group index
  @enddesc
  @calls      ApplyBndStatic

  @var        GH
  @vdesc      Pointer to CCTK grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        stencil
  @vdesc      stencil width
  @vtype      int [ dimension of group ]
  @vio        in
  @endvar
  @var        gi
  @vdesc      index of group to apply static boundaries to
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              return code of @seeroutine ApplyBndStatic <BR>
              -1 if invalid group indices are given
  @endreturndesc
@@*/
int BndStaticGI(const cGH *GH, const int *stencil, int gi) {
  int first_vi, retval;

  first_vi = CCTK_FirstVarIndexI(gi);
  if (first_vi >= 0) {
    retval = OldApplyBndStatic(GH, -1, stencil, 0, first_vi,
                               CCTK_NumVarsInGroupI(gi));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group index %d in BndStaticGI", gi);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticGI)(int *ierr, const cGH **GH,
                                        const int *stencil, const int *gi) {
  *ierr = BndStaticGI(*GH, stencil, *gi);
}

/* ======================================================= */

/*@@
   @routine    BndStaticDirGN
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply static boundary routines by group name in given direction
   @enddesc
   @calls      BndStaticDirGI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply static boundaries
   @vtype      int
   @vio        in
   @endvar
   @var        gname
   @vdesc      name of group to apply static boundaries to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndStaticDirGI <BR>
               -1 if invalid group names are given
   @endreturndesc
@@*/
int BndStaticDirGN(const cGH *GH, int stencil_size, int dir,
                   const char *gname) {
  int gi, num_groups, retval;

  gi = CCTK_GroupIndex(gname);
  num_groups = CCTK_NumGroups();

  if (gi >= 0 && gi < num_groups) {
    retval = BndStaticDirGI(GH, stencil_size, dir, gi);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group name '%s' in BndStaticDirGN", gname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticDirGN)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(gname)
  *ierr = BndStaticDirGN(*GH, *stencil_size, *dir, gname);
  free(gname);
}

/*@@
   @routine    BndStaticGN
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply static boundary routines by group name
   @enddesc
   @calls      BndStaticGI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of group ]
   @vio        in
   @endvar
   @var        gname
   @vdesc      name of group to apply static boundaries to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndStaticGI <BR>
               -1 if invalid group names are given
   @endreturndesc
@@*/
int BndStaticGN(const cGH *GH, const int *stencil, const char *gname) {
  int gi, num_groups, retval;

  gi = CCTK_GroupIndex(gname);
  num_groups = CCTK_NumGroups();

  if (gi >= 0 && gi < num_groups) {
    retval = BndStaticGI(GH, stencil, gi);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group name '%s' in BndStaticGN", gname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticGN)(int *ierr, const cGH **GH,
                                        const int *stencil,
                                        ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(gname)
  *ierr = BndStaticGN(*GH, stencil, gname);
  free(gname);
}

/* ======================================================= */

/*@@
   @routine    BndStaticDirVN
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply static boundary routines by group name in given direction
   @enddesc
   @calls      BndStaticDirVI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply static boundaries
   @vtype      int
   @vio        in
   @endvar
   @var        vname
   @vdesc      name of variable to apply static boundaries to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndStaticDirVI <BR>
               -1 if invalid variable names are given
   @endreturndesc
@@*/
int BndStaticDirVN(const cGH *GH, int stencil_size, int dir,
                   const char *vname) {
  int vi, num_vars, retval;

  vi = CCTK_VarIndex(vname);
  num_vars = CCTK_NumVars();

  if (vi >= 0 && vi < num_vars) {
    retval = BndStaticDirVI(GH, stencil_size, dir, vi);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable name '%s'  in BndStaticDirVN", vname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticDirVN)(int *ierr, const cGH **GH,
                                           const int *stencil_size,
                                           const int *dir, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(vname)
  *ierr = BndStaticDirVN(*GH, *stencil_size, *dir, vname);
  free(vname);
}

/*@@
   @routine    BndStaticVN
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply static boundary routines by variable name
   @enddesc
   @calls      BndStaticVI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        vname
   @vdesc      name of variable to apply static boundaries to
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndStaticVI <BR>
               -1 if invalid variable names are given
   @endreturndesc
@@*/
int BndStaticVN(const cGH *GH, const int *stencil, const char *vname) {
  int vi, num_vars, retval;

  vi = CCTK_VarIndex(vname);
  num_vars = CCTK_NumVars();

  if (vi >= 0 && vi < num_vars) {
    retval = BndStaticVI(GH, stencil, vi);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable name '%s' in BndStaticVN", vname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndStaticVN)(int *ierr, const cGH **GH,
                                        const int *stencil,
                                        ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(vname)
  *ierr = BndStaticVN(*GH, stencil, vname);
  free(vname);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/* maximum dimension we can deal with */
#define MAXDIM 3

/* macro to compute the linear index of a 3D point */
#define INDEX_3D(ash, i, j, k) ((i) + (ash)[0] * ((j) + (ash)[1] * (k)))

/*@@
   @routine    STATIC_BOUNDARY
   @date       Sat 20 Jan 2001
   @author     Thomas Radke
   @desc
               Macro to apply static boundary conditions to a variable
               Currently it is limited up to 3D variables only.
   @enddesc

   @var        doBC
   @vdesc      flag telling whether to apply boundary conditions or not
   @vtype      int
   @vio        in
   @endvar
   @var        iend, jend, kend
   @vdesc      upper ranges for the loopers
   @vtype      int
   @vio        in
   @endvar
   @var        ii, jj, kk
   @vdesc      indices of the current grid point
   @vtype      int
   @vio        in
   @endvar
@@*/
#define STATIC_BOUNDARY(doBC, iend, jend, kend, ii, jj, kk)                    \
  {                                                                            \
    if (doBC) {                                                                \
      for (k = 0; k < kend; k++) {                                             \
        for (j = 0; j < jend; j++) {                                           \
          for (i = 0; i < iend; i++) {                                         \
            int _index;                                                        \
                                                                               \
            _index = INDEX_3D(ash, ii, jj, kk) * vtypesize;                    \
            memcpy((char *)GH->data[var][timelvl_to] + _index,                 \
                   (char *)GH->data[var][timelvl_from] + _index, vtypesize);   \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    ApplyBndStatic
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply static boundary conditions to a group of grid functions
               given by their indices
               This routine is called by the various BndStaticXXX wrappers.

               Although it is currently limited to handle 1D, 2D, or 3D
               variables only it can easily be extended for higher dimensions
               by adapting the appropriate macros.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        width_dir
   @vdesc      boundary width in direction dir
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        in_widths
   @vdesc      boundary widths for all directions
   @vtype      CCTK_INT [ dimension of variable(s) ]
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction for static boundaries (0 for copying all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        first_var
   @vdesc      index of first variable to apply static boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupIndexFromVarI
               CCTK_GroupDimI
               CCTK_VarTypeI
               CCTK_GroupStaggerDirArrayGI
               STATIC_BOUNDARY
   @history
   @hdate      Sat 20 Jan 2001
   @hauthor    Thomas Radke
   @hdesc      Merged separate routines for 1D, 2D, and 3D
               into a single generic routine
   @endhistory

   @returntype int
   @returndesc
                0 for success
               -1 if dimension is not supported
               -2 if direction parameter is invalid
               -3 if stencil width array parameter is NULL
               -4 if there is only one timelevel
   @endreturndesc
@@*/
static int ApplyBndStatic(const cGH *GH, CCTK_INT width_dir,
                          const CCTK_INT *in_widths, int dir, int first_var,
                          int num_vars) {
  int ierr;
  int i, j, k;
  int timelvl_to, timelvl_from;
  int gindex, gdim;
  int var, vtypesize;
  int doBC[2 * MAXDIM], ash[MAXDIM], lsh[MAXDIM];
  CCTK_INT widths[2 * MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];

  /* Only apply boundary condition if more than one timelevel */
  if (CCTK_DeclaredTimeLevelsVI(first_var) <= 1) {
    return (-4);
  }

  /* get the group index of the target variable */
  gindex = CCTK_GroupIndexFromVarI(first_var);

  /* get the number of dimensions and the size of the variable's type */
  gdim = CCTK_GroupDimI(gindex);
  vtypesize = CCTK_VarTypeSize(CCTK_VarTypeI(first_var));

  /* make sure we can deal with this number of dimensions */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Variable dimension of %d not supported", gdim);
    return (-1);
  }

  /* check the direction parameter */
  if (abs(dir) > gdim) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndStatic: direction %d greater than dimension %d", dir,
               gdim);
    return (-2);
  }

  /* set up boundary width array */
  if (dir) {
    widths[2 * (abs(dir) - 1)] = width_dir;
    widths[2 * (abs(dir) - 1) + 1] = width_dir;
  } else if (in_widths) {
    memcpy(widths, in_widths, 2 * gdim * sizeof *widths);
  } else {
    CCTK_WARN(1, "ApplyBndStatic: NULL pointer passed for boundary width "
                 "array");
    return (-3);
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths(GH, first_var, gdim, widths, "Static");

  /* initialize arrays for variables with less dimensions than MAXDIM
     so that we can use the INDEX_3D macro later on */
  for (i = gdim; i < MAXDIM; i++) {
    ash[i] = 1;
    lsh[i] = 1;
  }

  /* get the current timelevel */
  timelvl_to = 0;
  timelvl_from = 1;

  /* see if we have a physical boundary */
  symtable = SymmetryTableHandleForGrid(GH);
  if (symtable < 0)
    CCTK_WARN(0, "internal error");
  ierr = Util_TableGetIntArray(symtable, 2 * gdim, symbnd, "symmetry_handle");
  if (ierr != 2 * gdim)
    CCTK_WARN(0, "internal error");
  for (i = 0; i < 2 * gdim; i++) {
    is_physical[i] = symbnd[i] < 0;
  }

  /* now loop over all variables */
  for (var = first_var; var < first_var + num_vars; var++) {
    if (CCTK_ActiveTimeLevelsVI(GH, var) < 2) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Static Boundary condition needs at least two timelevels "
                 "active, but %s only has %d.",
                 CCTK_FullName(var), CCTK_ActiveTimeLevelsVI(GH, var));
    }
    /* Apply condition if:
       + boundary is an outer boundary
       + have enough grid points
    */
    for (i = 0; i < 2 * gdim; i++) {
      doBC[i] = is_physical[i];
    }
    for (i = 0; i < gdim; i++) {
      ash[i] = GH->cctk_ash[i];
      lsh[i] = GH->cctk_lsh[i];
      doBC[i * 2] &= GH->cctk_lsh[i] > widths[i * 2] && GH->cctk_bbox[i * 2];
      doBC[i * 2 + 1] &=
          GH->cctk_lsh[i] > widths[i * 2 + 1] && GH->cctk_bbox[i * 2 + 1];
      if (dir != 0) {
        doBC[i * 2] &= (dir < 0 && (i + 1 == abs(dir)));
        doBC[i * 2 + 1] &= (dir > 0 && (i + 1 == abs(dir)));
      }
    }

    /* now copy the boundaries face by face */
    if (gdim > 0) {
      /* lower x */
      STATIC_BOUNDARY(doBC[0], widths[0], lsh[1], lsh[2], i, j, k);
      /* upper x */
      STATIC_BOUNDARY(doBC[1], widths[1], lsh[1], lsh[2], lsh[0] - i - 1, j, k);
    }
    if (gdim > 1) {
      /* lower y */
      STATIC_BOUNDARY(doBC[2], lsh[0], widths[2], lsh[2], i, j, k);
      /* upper y */
      STATIC_BOUNDARY(doBC[3], lsh[0], widths[3], lsh[2], i, lsh[1] - j - 1, k);
    }
    if (gdim > 2) {
      /* lower z */
      STATIC_BOUNDARY(doBC[4], lsh[0], lsh[1], widths[4], i, j, k);
      /* upper z */
      STATIC_BOUNDARY(doBC[5], lsh[0], lsh[1], widths[5], i, j, lsh[2] - k - 1);
    }
  }

  return (0);
}

/*@@
   @routine    OldApplyBndStatic
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
   @var        stencil_dir
   @vdesc      boundary width in direction dir
   @vtype      int
   @vio        in
   @endvar
   @var        stencil_alldirs
   @vdesc      boundary widths for all directions
   @vtype      int [ dimension of variable(s) ]
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction for static boundaries (0 for copying all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        first_var
   @vdesc      index of first variable to apply static boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupIndexFromVarI
               ApplyBndScalar
   @returntype int
   @returndesc
               returncode from @seeroutine ApplyBndScalar
   @endreturndesc
@@*/

static int OldApplyBndStatic(const cGH *GH, int stencil_dir,
                             const int *stencil_alldirs, int dir, int first_var,
                             int num_vars) {
  int retval, dim, i;
  CCTK_INT *boundary_widths;
  static int warned;

  /* Convert stencil_alldirs to new format */
  dim = CCTK_GroupDimFromVarI(first_var);
  boundary_widths = malloc(2 * dim * sizeof(CCTK_INT));
  for (i = 0; i < 2 * dim; ++i) {
    boundary_widths[i] = stencil_alldirs[i / 2];
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
  retval = ApplyBndStatic(GH, stencil_dir, boundary_widths, dir, first_var,
                          num_vars);

  free(boundary_widths);
  return retval;
}
