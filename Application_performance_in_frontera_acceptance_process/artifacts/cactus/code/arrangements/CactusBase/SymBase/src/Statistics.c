/*@@
  @file      Statistics.c
  @author    Erik Schnetter
  @date      2004-05-25
  @desc
             Provide information about the registered symmetry conditions
  @version   $Header$
  @enddesc
@@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_Table.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Statistics_c);

/*@@
  @routine    SymBase_Statistics
  @author     Erik Schnetter
  @date       2004-05-25
  @desc
              If verbose, print the symmetry boundary conditions for each
              face.
              Warn if symmetries do not have interpolation operators
              associated with them.
  @enddesc
@@*/

void SymBase_Statistics(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT sym_table;
  CCTK_INT symmetry_handle[100];
  CCTK_FPOINTER symmetry_interpolate[100];
  char const *symmetry_name;
  int face;
  int ierr;

  /* Get table */
  sym_table = SymmetryTableHandleForGrid(cctkGH);
  if (sym_table < 0)
    CCTK_WARN(0, "internal error");

  /* Get table entries */
  if (2 * cctkGH->cctk_dim > 100) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetIntArray(sym_table, 2 * cctkGH->cctk_dim, symmetry_handle,
                               "symmetry_handle");
  if (ierr != 2 * cctkGH->cctk_dim) {
    CCTK_WARN(0, "internal error");
  }
  ierr =
      Util_TableGetFPointerArray(sym_table, 2 * cctkGH->cctk_dim,
                                 symmetry_interpolate, "symmetry_interpolate");
  if (ierr != 2 * cctkGH->cctk_dim) {
    CCTK_WARN(0, "internal error");
  }

  /* Print information about the registered symmetry faces, if desired */

  if (verbose) {

    /* Loop over all faces */
    for (face = 0; face < 2 * cctkGH->cctk_dim; ++face) {
      if (symmetry_handle[face] >= 0) {
        symmetry_name = SymBase_SymmetryNameOfHandle(symmetry_handle[face]);
        if (!symmetry_name) {
          CCTK_WARN(0, "internal error");
        }
        if (face >= 6) {
          CCTK_WARN(0, "only 3 dimensions are supported so far");
        }
        CCTK_VInfo(CCTK_THORNSTRING, "Symmetry on %s %c-face: %s",
                   face % 2 == 0 ? "lower" : "upper", "xyz"[face / 2],
                   symmetry_name);
      }
    }
  }

  /* Warn about symmetries without interpolators */

  /* Loop over all registered symmetries */
  for (face = 0; face < 2 * cctkGH->cctk_dim; ++face) {
    if (symmetry_handle[face] >= 0) {
      if (symmetry_interpolate[face] == NULL) {
        symmetry_name = SymBase_SymmetryNameOfHandle(symmetry_handle[face]);
        if (!symmetry_name) {
          CCTK_WARN(0, "internal error");
        }
        if (face >= 6) {
          CCTK_WARN(0, "only 3 dimensions are supported so far");
        }
        CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The symmetry \"%s\" on the %s %c-face has not registered a "
                   "symmetry interpolator",
                   symmetry_name, face % 2 == 0 ? "lower" : "upper",
                   "xyz"[face / 2]);
      }
    }
  }
}
