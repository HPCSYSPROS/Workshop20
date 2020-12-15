/*@@
  @file      CopyBoundary.c
  @date      Mon Mar 15 15:09:00 1999
  @author    Gerd Lanfermann, Gabrielle Allen
  @desc
             Routines for applying copying-boundary conditions
  @enddesc
  @history
  @hdate     Sun 25 Feb 2001
  @hauthor   Thomas Radke
  @hdesc     BC routines generalized for applying to arbitrary CCTK data types
  @endhistory
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
CCTK_FILEVERSION(CactusBase_Boundary_CopyBoundary_c);

static int ApplyBndCopy(const cGH *GH, CCTK_INT stencil_dir,
                        const CCTK_INT *stencil_alldirs, int dir,
                        int first_var_to, int first_var_from, int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/*@@
   @routine    BndCopy
   @date       13 Feb 2003
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               the Copy boundary condition
   @enddesc
   @calls      ApplyBndCopy
               CCTK_GroupDimFromVarI
               Util_TableGetIntArray
               Util_TableQueryValueInfo
               CCTK_VWarn
               Util_TableGetString
               CCTK_VarIndex
               Util_TableGetInt

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
               return code of @seeroutine ApplyBndCopy
               -11 invalid table handle
               -12 no "COPY_FROM" key in table
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/
CCTK_INT BndCopy(const cGH *GH, CCTK_INT num_vars, CCTK_INT *vars,
                 CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *tables) {
  int i, j, k, gi, gdim, max_gdim, err, retval;
  CCTK_INT value_type, value_size;
  char *copy_from_name;

  /* variables to pass to ApplyBndCopy */
  CCTK_INT *width_alldirs; /* width of boundary on each face */
  int dir;                 /* direction in which to apply bc */
  CCTK_INT
      copy_from; /* variable (index) from which to copy the boundary data */

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
                 "Faces specification %d for Copy boundary conditions on "
                 "%s is not implemented yet.  "
                 "Applying Copy bcs to all (external) faces.",
                 (int)faces[i], CCTK_VarName(vars[i]));
    }
    dir = 0; /* apply bc to all faces */

    /* Look on table for copy-from variable */
    err = Util_TableQueryValueInfo(tables[i], &value_type, &value_size,
                                   "COPY_FROM");
    if (err == UTIL_ERROR_BAD_HANDLE) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Invalid table handle passed for Copy boundary "
                 "conditions for %s.  Name or index of variable to copy from "
                 "must be provided via key \"COPY_FROM\".  Aborting.",
                 CCTK_VarName(vars[i]));
      return -11;
    } else if (err == 1) {
      if (value_type == CCTK_VARIABLE_STRING) {
        copy_from_name = malloc(value_size * sizeof(char));
        Util_TableGetString(tables[i], value_size, copy_from_name, "COPY_FROM");
        copy_from = CCTK_VarIndex(copy_from_name);
        free(copy_from_name);
      } else if (value_type == CCTK_VARIABLE_INT) {
        Util_TableGetInt(tables[i], &copy_from, "COPY_FROM");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid data type for key \"COPY_FROM\" "
                   "Please use CCTK_STRING for the variable name, "
                   "or CCTK_INT for the variable index.");
      }
    } else {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "No key \"COPY_FROM\" provided in table.  Please enter the "
                 "name or index of variable to copy from into the table "
                 "under this key.  Aborting.");
      return -12;
    }

    /* Determine boundary width on all faces */
    /* (re-)allocate memory for buffer */
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
    if (!retval &&
        (retval = ApplyBndCopy(GH, 0, width_alldirs, dir, vars[i], copy_from,
                               j)) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndCopy() returned %d", retval);
    }
  }
  free(width_alldirs);

  return retval;
}

/* prototypes for external C routines are declared in header Boundary.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL CCTK_FNAME(BndCopyDirVI)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, const int *vi_to,
                                         const int *vi_from);
void CCTK_FCALL CCTK_FNAME(BndCopyVI)(int *ierr, const cGH **GH,
                                      const int *stencil, const int *vi_to,
                                      const int *vi_from);
void CCTK_FCALL CCTK_FNAME(BndCopyDirGI)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, const int *gi_to,
                                         const int *gi_from);
void CCTK_FCALL CCTK_FNAME(BndCopyGI)(int *ierr, const cGH **GH,
                                      const int *stencil, const int *gi_to,
                                      const int *gi_from);
void CCTK_FCALL CCTK_FNAME(BndCopyDirVN)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, TWO_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndCopyVN)(int *ierr, const cGH **GH,
                                      const int *stencil, TWO_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndCopyDirGN)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, TWO_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(BndCopyGN)(int *ierr, const cGH **GH,
                                      const int *stencil, TWO_FORTSTRING_ARG);

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndCopyDirVI
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply copy boundary routines by var index in given direction
   @enddesc
   @calls      ApplyBndCopy

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
   @vdesc      direction to copy
   @vtype      int
   @vio        in
   @endvar
   @var        vi_to
   @vdesc      index of variable to copy boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        vi_from
   @vdesc      index of variable to copy boundaries from
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndCopy <BR>
               -1 if invalid variable indices are given
   @endreturndesc
@@*/
int BndCopyDirVI(const cGH *GH, int stencil_size, int dir, int vi_to,
                 int vi_from) {
  int retval, num_vars;

  num_vars = CCTK_NumVars();
  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = ApplyBndCopy(GH, stencil_size, NULL, dir, vi_to, vi_from, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable indices %d and/or %d in BndCopyDirVI", vi_to,
               vi_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyDirVI)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, const int *vi_to,
                                         const int *vi_from) {
  *ierr = BndCopyDirVI(*GH, *stencil_size, *dir, *vi_to, *vi_from);
}

/*@@
   @routine    BndCopyVI
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply copy boundary routines by var index
   @enddesc
   @calls      ApplyBndCopy

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
   @var        vi_to
   @vdesc      index of variable to copy boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        vi_from
   @vdesc      index of variable to copy boundaries from
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndCopy <BR>
               -1 if invalid variable indices are given
   @endreturndesc
@@*/
int BndCopyVI(const cGH *GH, const int *stencil, int vi_to, int vi_from) {
  int retval, num_vars, dim, i;
  CCTK_INT *boundary_widths; /* boundary widths as expected by ApplyBndCopy */

  /* Set up boundary_widths array */
  dim = CCTK_GroupDimFromVarI(vi_to);
  boundary_widths = malloc(2 * dim * sizeof(CCTK_INT));
  for (i = 0; i < 2 * dim; i += 2) {
    boundary_widths[i] = stencil[i / 2];
  }

  num_vars = CCTK_NumVars();
  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = ApplyBndCopy(GH, -1, boundary_widths, 0, vi_to, vi_from, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable indices %d and/or %d in BndCopyVI", vi_to,
               vi_from);
    retval = -1;
  }

  free(boundary_widths);

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyVI)(int *ierr, const cGH **GH,
                                      const int *stencil, const int *vi_to,
                                      const int *vi_from) {
  *ierr = BndCopyVI(*GH, stencil, *vi_to, *vi_from);
}

/* ====================================================== */

/*@@
  @routine    BndCopyDirGI
  @date       Sat Jan 20 2001
  @author     Gabrielle Allen
  @desc
              Apply copy boundaries by group index in given direction
  @enddesc
  @calls      ApplyBndCopy

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
  @vdesc      direction to copy boundaries
  @vtype      int
  @vio        in
  @endvar
  @var        gi_to
  @vdesc      index of group to copy boundaries to
  @vtype      int
  @vio        in
  @endvar
  @var        gi_from
  @vdesc      index of group to copy boundaries from
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              return code of @seeroutine ApplyBndCopy <BR>
              -1 if invalid group indices are given
  @endreturndesc
@@*/
int BndCopyDirGI(const cGH *GH, int stencil_size, int dir, int gi_to,
                 int gi_from) {
  int first_vi_to, first_vi_from, retval;

  first_vi_to = CCTK_FirstVarIndexI(gi_to);
  first_vi_from = CCTK_FirstVarIndexI(gi_from);
  if (first_vi_to >= 0 && first_vi_from >= 0) {
    retval = ApplyBndCopy(GH, stencil_size, NULL, dir, first_vi_to,
                          first_vi_from, CCTK_NumVarsInGroupI(gi_to));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group indices %d and/or %d in BndCopyDirGI", gi_to,
               gi_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyDirGI)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, const int *gi_to,
                                         const int *gi_from) {
  *ierr = BndCopyDirGI(*GH, *stencil_size, *dir, *gi_to, *gi_from);
}

/*@@
  @routine    BndCopyGI
  @date       Thu Mar  2 11:07:11 2000
  @author     Gerd Lanfermann
  @desc
              Apply copy boundaries by group index
  @enddesc
  @calls      ApplyBndCopy

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
  @var        gi_to
  @vdesc      index of group to copy boundaries to
  @vtype      int
  @vio        in
  @endvar
  @var        gi_from
  @vdesc      index of group to copy boundaries from
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              return code of @seeroutine ApplyBndCopy <BR>
              -1 if invalid group indices are given
  @endreturndesc
@@*/
int BndCopyGI(const cGH *GH, const int *stencil, int gi_to, int gi_from) {
  int first_vi_to, first_vi_from, retval;
  int i, dim;
  CCTK_INT *boundary_widths;

  /* Set up boundary_widths array */
  dim = CCTK_GroupDimI(gi_to);
  boundary_widths = malloc(2 * dim * sizeof(CCTK_INT));
  for (i = 0; i < 2 * dim; i += 2) {
    boundary_widths[i] = stencil[i / 2];
  }

  first_vi_to = CCTK_FirstVarIndexI(gi_to);
  first_vi_from = CCTK_FirstVarIndexI(gi_from);
  if (first_vi_to >= 0 && first_vi_from >= 0) {
    retval = ApplyBndCopy(GH, -1, boundary_widths, 0, first_vi_to,
                          first_vi_from, CCTK_NumVarsInGroupI(gi_to));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group indices %d and/or %d in BndCopyGI", gi_to,
               gi_from);
    retval = -1;
  }

  free(boundary_widths);

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyGI)(int *ierr, const cGH **GH,
                                      const int *stencil, const int *gi_to,
                                      const int *gi_from) {
  *ierr = BndCopyGI(*GH, stencil, *gi_to, *gi_from);
}

/* ======================================================= */

/*@@
   @routine    BndCopyDirGN
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply copy boundary routines by group name in given direction
   @enddesc
   @calls      BndCopyDirGI

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
   @vdesc      direction to copy boundaries
   @vtype      int
   @vio        in
   @endvar
   @var        gname_to
   @vdesc      name of group to copy boundaries to
   @vtype      const char *
   @vio        in
   @endvar
   @var        gname_from
   @vdesc      name of group to copy boundaries from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndCopyDirGI <BR>
               -1 if invalid group names are given
   @endreturndesc
@@*/
int BndCopyDirGN(const cGH *GH, int stencil_size, int dir, const char *gname_to,
                 const char *gname_from) {
  int gi_to, gi_from, num_groups, retval;

  gi_to = CCTK_GroupIndex(gname_to);
  gi_from = CCTK_GroupIndex(gname_from);
  num_groups = CCTK_NumGroups();

  if (gi_to >= 0 && gi_to < num_groups && gi_from >= 0 &&
      gi_from < num_groups) {
    retval = BndCopyDirGI(GH, stencil_size, dir, gi_to, gi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group names '%s' and/or '%s' in BndCopyDirGN", gname_to,
               gname_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyDirGN)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, TWO_FORTSTRING_ARG) {
  TWO_FORTSTRINGS_CREATE(gname_to, gname_from)
  *ierr = BndCopyDirGN(*GH, *stencil_size, *dir, gname_to, gname_from);
  free(gname_to);
  free(gname_from);
}

/*@@
   @routine    BndCopyGN
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply copy boundary routines by group name
   @enddesc
   @calls      BndCopyGI

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
   @var        gname_to
   @vdesc      name of group to copy boundaries to
   @vtype      const char *
   @vio        in
   @endvar
   @var        gname_from
   @vdesc      name of group to copy boundaries from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndCopyGI <BR>
               -1 if invalid group names are given
   @endreturndesc
@@*/
int BndCopyGN(const cGH *GH, const int *stencil, const char *gname_to,
              const char *gname_from) {
  int gi_to, gi_from, num_groups, retval;
  int i, dim, *boundary_widths;

  gi_to = CCTK_GroupIndex(gname_to);
  gi_from = CCTK_GroupIndex(gname_from);
  num_groups = CCTK_NumGroups();

  /* Set up boundary_widths array */
  dim = CCTK_GroupDimI(gi_to);
  boundary_widths = malloc(2 * dim * sizeof(int));
  for (i = 0; i < 2 * dim; i += 2) {
    boundary_widths[i] = stencil[i / 2];
  }

  if (gi_to >= 0 && gi_to < num_groups && gi_from >= 0 &&
      gi_from < num_groups) {
    retval = BndCopyGI(GH, boundary_widths, gi_to, gi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group names '%s' and/or '%s' in BndCopyGN", gname_to,
               gname_from);
    retval = -1;
  }

  free(boundary_widths);

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyGN)(int *ierr, const cGH **GH,
                                      const int *stencil, TWO_FORTSTRING_ARG) {
  TWO_FORTSTRING_CREATE(gname_to, gname_from)
  *ierr = BndCopyGN(*GH, stencil, gname_to, gname_from);
  free(gname_to);
  free(gname_from);
}

/* ======================================================= */

/*@@
   @routine    BndCopyDirVN
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply copy boundary routines by group name in given direction
   @enddesc
   @calls      BndCopyDirVI

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
   @vdesc      direction to copy boundaries
   @vtype      int
   @vio        in
   @endvar
   @var        vname_to
   @vdesc      name of variable to copy boundaries to
   @vtype      const char *
   @vio        in
   @endvar
   @var        vname_from
   @vdesc      name of variable to copy boundaries from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndCopyDirVI <BR>
               -1 if invalid variable names are given
   @endreturndesc
@@*/
int BndCopyDirVN(const cGH *GH, int stencil_size, int dir, const char *vname_to,
                 const char *vname_from) {
  int vi_to, vi_from, num_vars, retval;

  vi_to = CCTK_VarIndex(vname_to);
  vi_from = CCTK_VarIndex(vname_from);
  num_vars = CCTK_NumVars();

  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = BndCopyDirVI(GH, stencil_size, dir, vi_to, vi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable names '%s' and/or '%s' in BndCopyDirVN",
               vname_to, vname_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyDirVN)(int *ierr, const cGH **GH,
                                         const int *stencil_size,
                                         const int *dir, TWO_FORTSTRING_ARG) {
  TWO_FORTSTRINGS_CREATE(vname_to, vname_from)
  *ierr = BndCopyDirVN(*GH, *stencil_size, *dir, vname_to, vname_from);
  free(vname_to);
  free(vname_from);
}

/*@@
   @routine    BndCopyVN
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply copy boundary routines by variable name
   @enddesc
   @calls      BndCopyVI

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
   @var        vname_to
   @vdesc      name of variable to copy boundaries to
   @vtype      const char *
   @vio        in
   @endvar
   @var        vname_from
   @vdesc      name of variable to copy boundaries from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndCopyVI <BR>
               -1 if invalid variable names are given
   @endreturndesc
@@*/
int BndCopyVN(const cGH *GH, const int *stencil, const char *vname_to,
              const char *vname_from) {
  int vi_to, vi_from, num_vars, retval;
  int i, dim;
  int *boundary_widths;

  vi_to = CCTK_VarIndex(vname_to);
  vi_from = CCTK_VarIndex(vname_from);
  num_vars = CCTK_NumVars();

  /* Set up boundary_widths array */
  dim = CCTK_GroupDimFromVarI(vi_to);
  boundary_widths = malloc(2 * dim * sizeof(int));
  for (i = 0; i < 2 * dim; i += 2) {
    boundary_widths[i] = stencil[i / 2];
  }

  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = BndCopyVI(GH, boundary_widths, vi_to, vi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable names '%s' and/or '%s' in BndCopyVN", vname_to,
               vname_from);
    retval = -1;
  }

  free(boundary_widths);

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(BndCopyVN)(int *ierr, const cGH **GH,
                                      const int *stencil,
                                      TWO_FORTSTRINGS_ARGS) {
  TWO_FORTSTRINGS_CREATE(vname_to, vname_from)
  *ierr = BndCopyVN(*GH, stencil, vname_to, vname_from);
  free(vname_to);
  free(vname_from);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/* maximum dimension we can deal with */
#define MAXDIM 3

/* macro to compute the linear index of a 3D point */
#define INDEX_3D(ash, i, j, k) ((i) + (ash)[0] * ((j) + (ash)[1] * (k)))

/*@@
   @routine    COPY_BOUNDARY
   @date       Sat 20 Jan 2001
   @author     Thomas Radke
   @desc
               Macro to apply copy boundary conditions to a variable
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
#define COPY_BOUNDARY(doBC, iend, jend, kend, ii, jj, kk)                      \
  {                                                                            \
    if (doBC) {                                                                \
      for (k = 0; k < kend; k++) {                                             \
        for (j = 0; j < jend; j++) {                                           \
          for (i = 0; i < iend; i++) {                                         \
            int _index;                                                        \
                                                                               \
            _index = INDEX_3D(ash, ii, jj, kk) * vtypesize;                    \
            memcpy((char *)GH->data[var_to][timelvl_to] + _index,              \
                   (char *)GH->data[var_from][timelvl_from] + _index,          \
                   vtypesize);                                                 \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    ApplyBndCopy
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply copy boundary conditions to a group of grid functions
               given by their indices
               This routine is called by the various BndCopyXXX wrappers.

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
   @vdesc      direction to copy boundaries (0 for copying all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        first_var_to
   @vdesc      index of first variable to copy boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        first_var_from
   @vdesc      index of first variable to copy boundaries from
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
               COPY_BOUNDARY
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
               -3 if boundary width array parameter is NULL
   @endreturndesc
@@*/
static int ApplyBndCopy(const cGH *GH, CCTK_INT width_dir,
                        const CCTK_INT *in_widths, int dir, int first_var_to,
                        int first_var_from, int num_vars) {
  int i, j, k;
  int timelvl_to, timelvl_from;
  int gindex, gdim;
  int var_to, var_from, vtypesize;
  int doBC[2 * MAXDIM], ash[MAXDIM], lsh[MAXDIM];
  CCTK_INT widths[2 * MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];
  CCTK_INT ierr;

  /* get the group index of the target variable */
  gindex = CCTK_GroupIndexFromVarI(first_var_to);

  /* get the number of dimensions and the size of the variable's type */
  gdim = CCTK_GroupDimI(gindex);
  vtypesize = CCTK_VarTypeSize(CCTK_VarTypeI(first_var_to));

  /* make sure we can deal with this number of dimensions */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Variable dimension of %d not supported", gdim);
    return (-1);
  }

  /* check the direction parameter */
  if (abs(dir) > gdim) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndCopy: direction %d greater than dimension %d", dir,
               gdim);
    return (-2);
  }

  /* set up stencil width array */
  if (dir) {
    widths[2 * (abs(dir) - 1)] = width_dir;
    widths[2 * (abs(dir) - 1) + 1] = width_dir;
  } else if (in_widths) {
    memcpy(widths, in_widths, 2 * gdim * sizeof *widths);
  } else {
    CCTK_WARN(1, "ApplyBndCopy: NULL pointer passed for boundary width "
                 "array");
    return (-3);
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths(GH, first_var_to, gdim, widths, "Copy");

  /* initialize arrays for variables with less dimensions than MAXDIM
     so that we can use the INDEX_3D macro later on */
  for (i = gdim; i < MAXDIM; i++) {
    ash[i] = 1;
    lsh[i] = 1;
  }

  /* get the current timelevel */
  timelvl_to = 0;
  timelvl_from = 0;

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
  for (var_to = first_var_to, var_from = first_var_from;
       var_to < first_var_to + num_vars; var_to++, var_from++) {
    /* Apply condition if:
       + boundary is a physical boundary
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
      COPY_BOUNDARY(doBC[0], widths[0], lsh[1], lsh[2], i, j, k);
      /* upper x */
      COPY_BOUNDARY(doBC[1], widths[1], lsh[1], lsh[2], lsh[0] - i - 1, j, k);
    }
    if (gdim > 1) {
      /* lower y */
      COPY_BOUNDARY(doBC[2], lsh[0], widths[2], lsh[2], i, j, k);
      /* upper y */
      COPY_BOUNDARY(doBC[3], lsh[0], widths[3], lsh[2], i, lsh[1] - j - 1, k);
    }
    if (gdim > 2) {
      /* lower z */
      COPY_BOUNDARY(doBC[4], lsh[0], lsh[1], widths[4], i, j, k);
      /* upper z */
      COPY_BOUNDARY(doBC[5], lsh[0], lsh[1], widths[5], i, j, lsh[2] - k - 1);
    }
  }

  return (0);
}
