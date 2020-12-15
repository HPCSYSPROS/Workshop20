/*@@
  @file      Check.c
  @author    Erik Schnetter
  @date      2004-06-19
  @desc
             Check whether the driver set up the grid consistently.
  @version   $Header$
  @enddesc
@@*/

#include "util_Table.h"

#include "cctk.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Check_c);

/*@@
  @routine    SymBase_Check
  @author     Erik Schnetter
  @date       2004-06-19
  @desc
              Check whether the driver set up the grid consistently.
  @enddesc
  @var        CCTK_ARGUMENTS
  @endvar
@@*/

void SymBase_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  struct SymBase const *symdata;
  CCTK_INT sym_table;
  CCTK_INT symmetry_handle[100];
  CCTK_INT symmetry_zone_width[100];
  int ierr;
  int d;

  /* Get table entries */
  symdata = CCTK_GHExtension(cctkGH, "SymBase");
  if (!symdata) {
    CCTK_WARN(0, "internal error");
  }

  sym_table = symdata->sym_table;

  if (2 * cctk_dim > 100) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetIntArray(sym_table, 2 * cctk_dim, symmetry_handle,
                               "symmetry_handle");
  if (ierr != 2 * cctk_dim) {
    CCTK_WARN(0, "internal error");
  }
  ierr = Util_TableGetIntArray(sym_table, 2 * cctk_dim, symmetry_zone_width,
                               "symmetry_zone_width");
  if (ierr != 2 * cctk_dim) {
    CCTK_WARN(0, "internal error");
  }

  ierr = 0;

  for (d = 0; d < cctk_dim; ++d) {

    /* Check lower boundary */
    if (symmetry_handle[2 * d] >= 0) {
      if (cctk_bbox[2 * d]) {
        if (cctk_lsh[d] - cctk_nghostzones[d] < symmetry_zone_width[2 * d]) {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "This processor's domain is too small for the boundary "
                     "condition at the lower boundary in direction %d",
                     d);
          ierr = 1;
        }
      }
    }

    /* Check upper boundary */
    if (symmetry_handle[2 * d + 1] >= 0) {
      if (cctk_bbox[2 * d + 1]) {
        if (cctk_lsh[d] - cctk_nghostzones[d] <
            symmetry_zone_width[2 * d + 1]) {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "This processor's domain is too small for the boundary "
                     "condition at the upper boundary in direction %d",
                     d);
          ierr = 1;
        }
      }
    }
  }

  if (ierr) {
    CCTK_WARN(0, "The grid setup is inconsistent with the boundary sizes.  One "
                 "or more processors' domains are too small.  The grid setup "
                 "is decided by the driver.  Try to make the driver lay out "
                 "the grids differently.");
  }
}
