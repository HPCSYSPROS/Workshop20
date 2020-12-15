/*@@
  @file      Table.c
  @author    Erik Schnetter
  @date      2004/03/07 09:48:53
  @desc
             Get the symmetry table handle for a grid or grid array
  @version   $Header$
  @enddesc
@@*/

#include "cctk.h"
#include "util_Table.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Table_c);

/*@@
  @routine    SymBase_SymmetryTableHandleForGrid
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Return the symmetry table handle for the grid hierarchy
  @enddesc
  @var        cctkGH
  @vtype      CCTK_POINTER_TO_CONST
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              >=0 symmetry table handle
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryTableHandleForGrid(CCTK_POINTER_TO_CONST const cctkGH_) {
  cGH const *const cctkGH = cctkGH_;
  struct SymBase const *symdata;

  if (!cctkGH) {
    CCTK_WARN(0, "internal error: cctkGH does not exist");
  }

  symdata = CCTK_GHExtension(cctkGH, "SymBase");
  if (!symdata) {
    CCTK_WARN(0, "internal error: SymBase GH extension does not exist");
  }

  return symdata->sym_table;
}

/*@@
  @routine    SymBase_SymmetryTableHandleForGI
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Return the symmetry table handle for a grid array
  @enddesc
  @var        cctkGH
  @vtype      CCTK_POINTER_TO_CONST
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
  @var        group_index
  @vtype      CCTK_INT
  @vdesc      Group
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              >=0 symmetry table handle
              -6 if group_index has an illegal value
              -7 if the group type has an illegal value
              -8 if there was an internal error
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryTableHandleForGI(CCTK_POINTER_TO_CONST const cctkGH_,
                                 CCTK_INT const group_index) {
  cGH const *const cctkGH = cctkGH_;
  struct SymBase const *symdata;

  if (!cctkGH) {
    CCTK_WARN(0, "internal error: cctkGH does not exist");
  }

  symdata = CCTK_GHExtension(cctkGH, "SymBase");
  if (!symdata) {
    CCTK_WARN(0, "internal error: SymBase GH extension does not exist");
  }

  if (group_index < 0 || group_index >= CCTK_NumGroups()) {
    return -6; /* illegal argument */
  }

  switch (CCTK_GroupTypeI(group_index)) {
  case CCTK_GF:
    return -7; /* illegal group type */
  case CCTK_SCALAR:
  case CCTK_ARRAY:
    return symdata->array_sym_tables[group_index];
  default:
    CCTK_WARN(0, "internal error");
  }

  return -8; /* internal error */
}

/*@@
  @routine    SymBase_SymmetryTableHandleForGN
  @author     Erik Schnetter
  @date       2004-03-06
  @desc
              Return the symmetry table handle for a grid array
  @enddesc
  @var        cctkGH
  @vtype      CCTK_POINTER_TO_CONST
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
  @var        group_name
  @vtype      CCTK_STRING
  @vdesc      Group
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              >=0 symmetry table handle
              Error codes of CCTK_GroupIndex
              Error codes of SymBase_SymmetryTableHandleForGI
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryTableHandleForGN(CCTK_POINTER_TO_CONST const cctkGH_,
                                 CCTK_STRING const group_name) {
  int group_index;

  group_index = CCTK_GroupIndex(group_name);
  if (group_index < 0) {
    return group_index; /* illegal argument */
  }

  return SymBase_SymmetryTableHandleForGI(cctkGH_, group_index);
}

/*@@
  @routine    SymBase_GetSymmetryBoundaries
  @author     Erik Schnetter
  @date       2006-05-16
  @desc
              Determine which boundaries are symmetry boundaries
  @enddesc
  @var        cctkGH
  @vtype      CCTK_POINTER_TO_CONST
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
  @var        size
  @vtype      CCTK_INT
  @vdesc      Array size
  @vio        in
  @endvar
  @var        symbnd
  @vtype      CCTK_INT [size]
  @vdesc      0 for outer boundary, 1 for symmetry boundary
  @vio        out
  @endvar
  @returntype CCTK_INT
  @returndesc
              0 Success
              -10 Wrong array size; must be 2*cctk_dim
              -11 Symmetry table handle does not exist
              -12 Wrong number of symmetry table entries (internal error)
  @endreturndesc
@@*/

CCTK_INT
SymBase_GetSymmetryBoundaries(CCTK_POINTER_TO_CONST const cctkGH_,
                              CCTK_INT const size,
                              CCTK_INT *restrict const symbnd) {
  cGH const *const cctkGH = cctkGH_;
  CCTK_INT symtable;
  CCTK_INT symmetry_handle[6];
  CCTK_INT iret;
  int d;

  /* Check the arguments */
  if (size != 2 * cctkGH->cctk_dim) {
    return -10;
  }

  /* Get the symmetry table */
  symtable = SymmetryTableHandleForGrid(cctkGH);
  if (symtable < 0) {
    return -11;
  }

  /* Get the symmetry handles for each face */
  iret = Util_TableGetIntArray(symtable, 2 * cctkGH->cctk_dim, symmetry_handle,
                               "symmetry_handle");
  if (iret != 2 * cctkGH->cctk_dim) {
    return -12;
  }

  /* A face has a symmetry boundary if there is a valid (non-negative)
     handle registered for that face */
  for (d = 0; d < 2 * cctkGH->cctk_dim; ++d) {
    symbnd[d] = symmetry_handle[d] >= 0;
  }

  /* Success */
  return 0;
}
