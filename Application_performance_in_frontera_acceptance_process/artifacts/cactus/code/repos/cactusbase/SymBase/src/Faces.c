/*@@
  @file      Faces.c
  @author    Erik Schnetter
  @date      2004/03/07 09:48:53
  @desc
             Register a symmetry condition for a face
  @version   $Header$
  @enddesc
@@*/

#include <stdlib.h>

#include "cctk.h"

#include "util_Table.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Faces_c);

/*@@
  @routine    SymBase_SymmetryRegisterFaces
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Register a symmetry for certain faces
  @enddesc
  @var        sym_table
  @vtype      CCTK_INT
  @vdesc      Table which describes the grid or grid array
  @vio        in
  @endvar
  @var        group_dim
  @vtype      CCTK_INT
  @vdesc      Dimension of the grid or grid array
  @vio        in
  @endvar
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry handle
  @vio        in
  @endvar
  @var        which_faces
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      The set of faces
  @vio        in
  @endvar
  @var        symmetry_zone_width
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      Symmetry boundary width
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              -1 if sym_table has an illegal value
              -9 if group_dim has an illegal value
              -2 if sym_handle has an illegal value
              -3 if which_faces has an illegal value
              -4 if symmetry_zone_width has an illegal value
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegisterFaces(CCTK_INT const sym_table,
                              CCTK_INT const group_dim,
                              CCTK_INT const sym_handle,
                              CCTK_INT const *const which_faces,
                              CCTK_INT const *const new_symmetry_zone_width) {
  CCTK_INT symmetry_handle[100];
  CCTK_INT symmetry_zone_width[100];
  int face;
  int ierr;

  /* Check arguments */
  if (sym_table < 0) {
    return -1; /* illegal argument */
  }
  if (group_dim < 0) {
    return -9; /* illegal argument */
  }
  if (sym_handle < 0 || sym_handle >= SymBase_num_symmetries) {
    return -2; /* illegal argument */
  }
  if (!which_faces) {
    return -3; /* illegal argument */
  }
  if (!new_symmetry_zone_width) {
    return -4; /* illegal argument */
  }

  /* Get table entries */
  if (2 * group_dim > 100) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetIntArray(sym_table, 2 * group_dim, symmetry_handle,
                               "symmetry_handle");
  if (ierr != 2 * group_dim) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetIntArray(sym_table, 2 * group_dim, symmetry_zone_width,
                               "symmetry_zone_width");
  if (ierr != 2 * group_dim) {
    CCTK_WARN(0, "internal error");
  }

  /* Update table entries */
  for (face = 0; face < 2 * group_dim; ++face) {
    if (which_faces[face]) {
      if (symmetry_handle[face] != -1) {
        return -5; /* The face is already taken */
      }
      symmetry_handle[face] = sym_handle;
      if (new_symmetry_zone_width[face] > 100) {
        CCTK_WARN(0, "Tried to register a symmetry face with a symmetry zone "
                     "width >100 -- this looks like an error.  This can be "
                     "caused by using a symmetry thorn (or the old built-in "
                     "CartGrid3D symmetries) which does not use CoordBase with "
                     "new versions of the Carpet driver.  Please use a "
                     "separate symmetry thorn (e.g. ReflectionSymmetry) "
                     "instead.");
      }
      symmetry_zone_width[face] = new_symmetry_zone_width[face];
    }
  }

  /* Set table entries */
  ierr = Util_TableSetIntArray(sym_table, 2 * group_dim, symmetry_handle,
                               "symmetry_handle");
  if (ierr != 1) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableSetIntArray(sym_table, 2 * group_dim, symmetry_zone_width,
                               "symmetry_zone_width");
  if (ierr != 1) {
    CCTK_WARN(0, "internal error");
  }

  return 0;
}

/*@@
  @routine    SymBase_SymmetryRegisterGrid
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Register a symmetry for certain faces of the grid hierarchy
  @enddesc
  @var        cctkGH
  @vtype      cGH *
  @vdesc      The grid hierarchy
  @vio        in
  @endvar
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry handle
  @vio        in
  @endvar
  @var        which_faces
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      The set of faces
  @vio        in
  @endvar
  @var        symmetry_zone_width
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      Symmetry boundary width
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              Error codes of SymBase_SymmetryRegisterFaces
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegisterGrid(CCTK_POINTER const cctkGH_,
                             CCTK_INT const sym_handle,
                             CCTK_INT const *const which_faces,
                             CCTK_INT const *const new_symmetry_zone_width) {
  cGH const *const cctkGH = cctkGH_;
  struct SymBase const *symdata;

  if (!cctkGH) {
    CCTK_WARN(0, "internal error");
  }

  symdata = CCTK_GHExtension(cctkGH, "SymBase");
  if (!symdata) {
    CCTK_WARN(0, "internal error");
  }

  return SymBase_SymmetryRegisterFaces(symdata->sym_table, cctkGH->cctk_dim,
                                       sym_handle, which_faces,
                                       new_symmetry_zone_width);
}

/*@@
  @routine    SymBase_SymmetryRegisterGI
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Register a symmetry for certain faces of grid arrays
  @enddesc
  @var        cctkGH
  @vtype      cGH *
  @vdesc      The grid hierarchy
  @vio        in
  @endvar
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry handle
  @vio        in
  @endvar
  @var        which_faces
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      The set of faces
  @vio        in
  @endvar
  @var        symmetry_zone_width
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      Symmetry boundary width
  @vio        in
  @endvar
  @var        group_index
  @vtype      CCTK_INT
  @vdesc      Grid array group
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              -6 if group_index has an illegal value
              -7 if the group has an illegal type
              -8 internal error
              Error codes of SymBase_SymmetryRegisterFaces
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegisterGI(CCTK_POINTER const cctkGH_,
                           CCTK_INT const sym_handle,
                           CCTK_INT const *const which_faces,
                           CCTK_INT const *const new_symmetry_zone_width,
                           CCTK_INT const group_index) {
  cGH const *const cctkGH = cctkGH_;
  struct SymBase const *symdata;

  if (!cctkGH) {
    CCTK_WARN(0, "internal error");
  }

  symdata = CCTK_GHExtension(cctkGH, "SymBase");
  if (!symdata) {
    CCTK_WARN(0, "internal error");
  }

  if (group_index < 0 || group_index >= CCTK_NumGroups()) {
    return -6; /* illegal argument */
  }

  switch (CCTK_GroupTypeI(group_index)) {
  case CCTK_GF:
    return -7; /* illegal group type */
  case CCTK_SCALAR:
  case CCTK_ARRAY:
    return SymBase_SymmetryRegisterFaces(
        symdata->array_sym_tables[group_index], CCTK_GroupDimI(group_index),
        sym_handle, which_faces, new_symmetry_zone_width);
  default:
    CCTK_WARN(0, "internal error");
  }

  return -8; /* internal error */
}

/*@@
  @routine    SymBase_SymmetryRegisterGI
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Register a symmetry for certain faces of grid arrays
  @enddesc
  @var        cctkGH
  @vtype      cGH *
  @vdesc      The grid hierarchy
  @vio        in
  @endvar
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry handle
  @vio        in
  @endvar
  @var        which_faces
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      The set of faces
  @vio        in
  @endvar
  @var        symmetry_zone_width
  @vtype      CCTK_INT[2*group_dim]
  @vdesc      Symmetry boundary width
  @vio        in
  @endvar
  @var        group_name
  @vtype      CCTK_STRING
  @vdesc      Grid array group
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              Error codes of CCTK_GroupIndex
              Error codes of SymBase_SymmetryRegisterGI
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegisterGN(CCTK_POINTER const cctkGH_,
                           CCTK_INT const sym_handle,
                           CCTK_INT const *const which_faces,
                           CCTK_INT const *const new_symmetry_zone_width,
                           CCTK_STRING const group_name) {
  int group_index;

  group_index = CCTK_GroupIndex(group_name);
  if (group_index < 0) {
    return group_index; /* illegal argument */
  }

  return SymBase_SymmetryRegisterGI(cctkGH_, sym_handle, which_faces,
                                    new_symmetry_zone_width, group_index);
}

/*@@
  @routine    SymBase_SymmetryRegisterInterpolatorFaces
  @author     Erik Schnetter
  @date       2004-05-25
  @desc
              Register a symmetry interpolator for certain faces
  @enddesc
  @var        sym_table
  @vtype      CCTK_INT
  @vdesc      Table which describes the grid or grid array
  @vio        in
  @endvar
  @var        group_dim
  @vtype      CCTK_INT
  @vdesc      Dimension of the grid or grid array
  @vio        in
  @endvar
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry handle
  @vio        in
  @endvar
  @var        symmetry_interpolate
  @vtype      CCTK_FPOINTER
  @vdesc      Routine that applies symmetries to the interpolation points
              and then calls back
              (may be NULL)
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              -1 if sym_table has an illegal value
              -9 if group_dim has an illegal value
              -2 if sym_handle has an illegal value
              -10 symmetry_interpolate is NULL
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegisterInterpolatorFaces(
    CCTK_INT const sym_table, CCTK_INT const group_dim,
    CCTK_INT const sym_handle,
    SymmetryInterpolateFPointer const new_symmetry_interpolate) {
  CCTK_INT symmetry_handle[100];
  CCTK_FPOINTER symmetry_interpolate[100];
  int face;
  int ierr;

  /* Check arguments */
  if (sym_table < 0) {
    return -1; /* illegal argument */
  }
  if (group_dim < 0) {
    return -9; /* illegal argument */
  }
  if (sym_handle < 0 || sym_handle >= SymBase_num_symmetries) {
    return -2; /* illegal argument */
  }
  if (!new_symmetry_interpolate) {
    return -10;
  }

  /* Get table entries */
  if (2 * group_dim > 100) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetIntArray(sym_table, 2 * group_dim, symmetry_handle,
                               "symmetry_handle");
  if (ierr != 2 * group_dim) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetFPointerArray(
      sym_table, 2 * group_dim, symmetry_interpolate, "symmetry_interpolate");
  if (ierr != 2 * group_dim) {
    CCTK_WARN(0, "internal error");
  }

  /* Update table entries */
  for (face = 0; face < 2 * group_dim; ++face) {
    if (symmetry_handle[face] == sym_handle) {
      symmetry_interpolate[face] = (CCTK_FPOINTER)new_symmetry_interpolate;
    }
  }

  /* Set table entries */
  ierr = Util_TableSetFPointerArray(
      sym_table, 2 * group_dim, symmetry_interpolate, "symmetry_interpolate");
  if (ierr != 1) {
    CCTK_WARN(0, "internal error");
  }

  return 0;
}

/*@@
  @routine    SymBase_SymmetryRegisterGridInterpolator
  @author     Erik Schnetter
  @date       2004-05-25
  @desc
              Register a symmetry interpolator
              for certain faces of the grid hierarchy
  @enddesc
  @var        cctkGH
  @vtype      cGH *
  @vdesc      The grid hierarchy
  @vio        in
  @endvar
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry handle
  @vio        in
  @endvar
  @var        symmetry_interpolate
  @vtype      CCTK_FPOINTER
  @vdesc      Routine that applies symmetries to the interpolation points
              and then calls back
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              Error codes of SymBase_SymmetryRegisterInterpolatorFaces
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegisterGridInterpolator(
    CCTK_POINTER const cctkGH_, CCTK_INT const sym_handle,
    SymmetryInterpolateFPointer const new_symmetry_interpolate) {
  cGH const *const cctkGH = cctkGH_;
  struct SymBase const *symdata;

  if (!cctkGH) {
    CCTK_WARN(0, "internal error");
  }

  symdata = CCTK_GHExtension(cctkGH, "SymBase");
  if (!symdata) {
    CCTK_WARN(0, "internal error");
  }

  return SymBase_SymmetryRegisterInterpolatorFaces(symdata->sym_table,
                                                   cctkGH->cctk_dim, sym_handle,
                                                   new_symmetry_interpolate);
}
