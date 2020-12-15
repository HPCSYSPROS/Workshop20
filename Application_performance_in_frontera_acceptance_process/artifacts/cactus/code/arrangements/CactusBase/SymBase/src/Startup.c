/*@@
  @file      Startup.c
  @author    Erik Schnetter
  @date      2004-03-06
  @desc
             Thorn startup
  @version   $Id$
  @enddesc
@@*/

#include <stdlib.h>

#include "cctk.h"

#include "util_Table.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Startup_c);

/*@@
  @routine    SymBase_Startup
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Register a GH extension and initialise the global variables
  @enddesc
  @returntype int
  @returndesc
              status code
              0 for success
  @endreturndesc
@@*/

int SymBase_Startup(void) {
  int handle;
  int ierr;

  handle = CCTK_RegisterGHExtension("SymBase");
  if (handle < 0) {
    CCTK_WARN(
        0, "Internal Error: Failed to create GH extension for symmetry data");
  }

  ierr = CCTK_RegisterGHExtensionSetupGH(handle, SymBase_Setup);
  if (ierr != 1) { /* strange error code convention */
    CCTK_WARN(0,
              "Internal Error: Could not register GH Extension setup routine");
  }

  SymBase_num_symmetries = 0;
  SymBase_symmetry_names = NULL;

  return 0; /* no error */
}

/*@@
  @routine    SymBase_Setup
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Initialise the GH extension
  @enddesc
  @var        config
  @vtype      tFleshConfig *
  @vdesc      Internal flesh information, unused
  @vio        in
  @endvar
  @var        convlev
  @vtype      int
  @vdesc      Convergence level, unused
  @vio        in
  @endvar
  @var        cctkGH
  @vtype      cGH *
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
@@*/

void *SymBase_Setup(tFleshConfig *const config, int const convlev,
                    cGH *const cctkGH) {
  struct SymBase *symdata;
  CCTK_INT symmetry_handle[100];
  CCTK_INT symmetry_zone_width[100];
  CCTK_FPOINTER symmetry_interpolate[100];
  int group;
  int face;
  int ierr;
  const void *dummy;

  /* avoid compiler warnings about unused parameters */
  dummy = &config;
  dummy = &convlev;
  dummy = &cctkGH;
  dummy = &dummy;

  /* Create GH extension */
  symdata = malloc(sizeof *symdata);
  if (!symdata) {
    CCTK_WARN(0, "Failed to allocate memory for symmetry data");
  }

  /* Initialise handle and zone width arrays */
  for (face = 0; face < 100; ++face) {
    symmetry_handle[face] = -1;
    symmetry_zone_width[face] = 0;
    symmetry_interpolate[face] = NULL;
  }

  /* Create grid symmetry table */
  symdata->sym_table = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);

  if (symdata->sym_table < 0) {
    CCTK_WARN(0, "Internal Error: Failed to create symmetry table.");
  }

  if (2 * cctkGH->cctk_dim > 100) {
    CCTK_WARN(0, "Internal Error: Can currently only handle a hundred faces.");
  }
  ierr = Util_TableSetIntArray(symdata->sym_table, 2 * cctkGH->cctk_dim,
                               symmetry_handle, "symmetry_handle");
  if (ierr) {
    CCTK_WARN(0, "Internal Error: Failed to create symmetry_handle array");
  }
  ierr = Util_TableSetIntArray(symdata->sym_table, 2 * cctkGH->cctk_dim,
                               symmetry_zone_width, "symmetry_zone_width");
  if (ierr) {
    CCTK_WARN(0, "Internal Error: Failed to create symmetry_zone_width array");
  }
  ierr =
      Util_TableSetFPointerArray(symdata->sym_table, 2 * cctkGH->cctk_dim,
                                 symmetry_interpolate, "symmetry_interpolate");
  if (ierr) {
    CCTK_WARN(0, "Internal Error: Failed to create symmetry_fold_points array");
  }

  /* Create grid array symmetry tables */
  symdata->array_sym_tables =
      malloc(CCTK_NumGroups() * sizeof *symdata->array_sym_tables);
  if (!symdata->array_sym_tables) {
    CCTK_WARN(
        0,
        "Internal Error: Failed to create array for symmetry tables for GAs");
  }

  for (group = 0; group < CCTK_NumGroups(); ++group) {
    switch (CCTK_GroupTypeI(group)) {
    case CCTK_GF:

      /* No table */
      symdata->array_sym_tables[group] = -1;
      break;

    case CCTK_SCALAR:
    /* FALL_THROUGH */

    case CCTK_ARRAY:

      /* Create array symmetry table */
      symdata->array_sym_tables[group] =
          Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
      if (symdata->array_sym_tables[group] < 0) {
        CCTK_WARN(0, "Internal Error: Failed to create symmetry table for GA");
      }

      if (2 * CCTK_GroupDimI(group) > 100) {
        CCTK_WARN(0, "Internal Error: Can only handle 100 faces");
      }
      ierr = Util_TableSetIntArray(symdata->array_sym_tables[group],
                                   2 * CCTK_GroupDimI(group), symmetry_handle,
                                   "symmetry_handle");
      if (ierr) {
        CCTK_WARN(0, "Internal Error: Failed to create symmetry_handle array");
      }
      ierr = Util_TableSetIntArray(symdata->array_sym_tables[group],
                                   2 * CCTK_GroupDimI(group),
                                   symmetry_zone_width, "symmetry_zone_width");
      if (ierr) {
        CCTK_WARN(0,
                  "Internal Error: Failed to create symmetry_zone_width array");
      }
      ierr = Util_TableSetFPointerArray(
          symdata->array_sym_tables[group], 2 * CCTK_GroupDimI(group),
          symmetry_interpolate, "symmetry_interpolate");
      if (ierr) {
        CCTK_WARN(
            0, "Internal Error: Failed to create symmetry_fold_points array");
      }
      break;

    default:
      CCTK_WARN(0, "Internal Error: Unknown GV type");
    }
  }

  return symdata;
}
