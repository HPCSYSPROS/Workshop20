/*@@
  @file      Handles.c
  @author    Erik Schnetter
  @date      2004/03/07 09:48:53
  @desc
             Register a symmetry, and query the symmetry name/handle mapping
  @version   $Id$
  @enddesc
@@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Handles_c);

/* Number of registered symmetries */
int SymBase_num_symmetries;

/* The names of these symmetries */
const char **SymBase_symmetry_names;

/*@@
  @routine    SymBase_SymmetryRegister
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Register a symmetry
  @enddesc
  @var        sym_name
  @vtype      CCTK_STRING
  @vdesc      Name of the symmetry
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              >=0 symmetry handle
              -1 if sym_name has an illegal value
              -2 if a symmetry with the same name has already been registered
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryRegister(CCTK_STRING const sym_name) {
  int n;

  if (!sym_name) {
    return -1; /* illegal argument */
  }

  for (n = 0; n < SymBase_num_symmetries; ++n) {
    if (strcmp(sym_name, SymBase_symmetry_names[n]) == 0) {
      return -2; /* symmetry exists already */
    }
  }

  SymBase_symmetry_names =
      realloc(SymBase_symmetry_names,
              (SymBase_num_symmetries + 1) * sizeof *SymBase_symmetry_names);
  if (!SymBase_symmetry_names) {
    CCTK_WARN(0, "internal error");
  }

  SymBase_symmetry_names[SymBase_num_symmetries] = strdup(sym_name);
  if (!SymBase_symmetry_names[SymBase_num_symmetries]) {
    CCTK_WARN(0, "internal error");
  }

  ++SymBase_num_symmetries;

  return SymBase_num_symmetries - 1;
}

/*@@
  @routine    SymBase_SymmetryHandleOfName
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Map a symmetry name to a symmetry handle
  @enddesc
  @var        sym_name
  @vtype      CCTK_STRING
  @vdesc      Name of the symmetry
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              >=0 symmetry handle
              -1 if sym_name has an illegal value
              -2 if no symmetry with that name has been registered
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryHandleOfName(CCTK_STRING const sym_name) {
  int n;

  if (!sym_name) {
    return -1; /* illegal argument */
  }

  for (n = 0; n < SymBase_num_symmetries; ++n) {
    if (strcmp(sym_name, SymBase_symmetry_names[n]) == 0) {
      return n; /* found */
    }
  }

  return -2; /* not found */
}

/*@@
  @routine    SymBase_SymmetryNameOfHandle
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Map a symmetry handle to a symmetry name
  @enddesc
  @var        sym_handle
  @vtype      CCTK_INT
  @vdesc      Symmetry name
  @vio        in
  @endvar
  @returntype CCTK_POINTER_TO_CONST
  @returndesc
              char const * containing the symmetry name
              NULL if no symmetry with that handle has been registered
  @endreturndesc
@@*/

CCTK_POINTER_TO_CONST
SymBase_SymmetryNameOfHandle(CCTK_INT const sym_handle) {
  if (sym_handle < 0 || sym_handle >= SymBase_num_symmetries) {
    return NULL; /* illegal argument */
  }

  return SymBase_symmetry_names[sym_handle];
}
