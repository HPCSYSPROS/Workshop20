/*@@
  @file      Interpolation.c
  @author    Erik Schnetter
  @date      2004-05-11
  @desc
             Apply symmetry conditions during interpolation
  @version   $Header$
  @enddesc
@@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "SymBase.h"

/* the rcs ID and its dummy function to use it */
static const char *const rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_SymBase_Interpolation_c);

/*@@
  @routine    SymBase_SymmetryInterpolate
  @author     Erik Schnetter
  @date       2004-05-11
  @desc
              Adjust the coordinates of the interpolation points,
              call the driver's interpolation function,
              and then adjust the interpolated tensor components.
              All arguments are equivalent to the CCTK_InterpGridArrays
              function call.
  @enddesc
  @var        cctkGH
  @vtype      CCTK_POINTER_TO_CONST
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
  @var        N_dims
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        local_interp_handle
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        param_table_handle
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        coord_system_handle
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        N_interp_points
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        interp_coords_type
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        interp_coords
  @vtype      CCTK_POINTER_TO_CONST[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        N_input_arrays
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        input_array_indices
  @vtype      CCTK_INT[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        N_output_arrays
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        output_array_types
  @vtype      CCTK_INT[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        output_arrays
  @vtype      CCTK_POINTER[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              -91 GH is NULL
              -92 N_dims is not equal to cctkGH->cctk_dim
              -93 Input array has illegal group type
              -94 Mixing grid functions and grid arrays
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryInterpolate(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[]) {
  cGH const *restrict const cctkGH = cctkGH_;
  CCTK_INT sym_table;
  CCTK_FPOINTER symmetry_interpolate[100];
  CCTK_INT faces;
  int has_grid_arrays, has_only_grid_arrays;
  int n;
  int f;
  int ierr;

  /* Check arguments */
  if (!cctkGH) {
    CCTK_WARN(1, "cctkGH is NULL");
    return -91;
  }

  /* If this interpolates grid arrays, then do not apply symmetries --
     there are no symmetries registered for grid arrays (yet).  */
  has_grid_arrays = 0;
  has_only_grid_arrays = 1;
  for (n = 0; n < N_input_arrays; ++n) {
    if (input_array_indices[n] != -1) {
      const int grouptype = CCTK_GroupTypeFromVarI(input_array_indices[n]);
      switch (grouptype) {
      case CCTK_GF:
        has_only_grid_arrays = 0;
        break;
      case CCTK_SCALAR:
      case CCTK_ARRAY:
        has_grid_arrays = 1;
        break;
      default: {
        char *groupname = CCTK_GroupName(input_array_indices[n]);
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not determine group type of input arrays for group "
                   "'%s' (%ld), received %d",
                   groupname, (long)input_array_indices[n], grouptype);
        free(groupname);
      }
        return -93;
      }
    }
  }
  if (has_grid_arrays && !has_only_grid_arrays) {
    CCTK_WARN(1, "The input arrays contain both grid function and grid arrays; "
                 "this is not possible");
    return -94;
  }
  if (has_grid_arrays) {
    /* Call the real interpolator */
    ierr = CCTK_IsFunctionAliased("DriverInterpolate");
    if (!ierr) {
      CCTK_WARN(
          0,
          "The aliased function \"DriverInterpolate\" has not been provided");
    }
    return DriverInterpolate(
        cctkGH, N_dims, local_interp_handle, param_table_handle,
        coord_system_handle, N_interp_points, interp_coords_type, interp_coords,
        N_input_arrays, input_array_indices, N_output_arrays,
        output_array_types, output_arrays);
  }

  /* Check arguments */
  if (N_dims != cctkGH->cctk_dim) {
    CCTK_WARN(1, "The number of dimensions is not equal to the GH's number of "
                 "dimensions");
    return -92;
  }

  /* Get table */
  sym_table = SymmetryTableHandleForGrid(cctkGH);
  if (sym_table < 0)
    CCTK_WARN(0, "internal error");

  /* Get table entries */
  if (2 * cctkGH->cctk_dim > 100) {
    CCTK_WARN(0, "internal error");
  }
  ierr =
      Util_TableGetFPointerArray(sym_table, 2 * cctkGH->cctk_dim,
                                 symmetry_interpolate, "symmetry_interpolate");
  if (ierr != 2 * cctkGH->cctk_dim) {
    CCTK_WARN(0, "internal error");
  }

  /* Find all faces that need symmetries applied */
  faces = 0;
  for (f = 0; f < 2 * cctkGH->cctk_dim; ++f) {
    if (symmetry_interpolate[f]) {
      faces |= (1 << f);
    }
  }

  /* Forward the call */
  return SymBase_SymmetryInterpolateFaces(
      cctkGH, N_dims, local_interp_handle, param_table_handle,
      coord_system_handle, N_interp_points, interp_coords_type, interp_coords,
      N_input_arrays, input_array_indices, N_output_arrays, output_array_types,
      output_arrays, faces);
}

/*@@
  @routine    SymBase_SymmetryInterpolateFaces
  @author     Erik Schnetter
  @date       2004-05-11
  @desc
              Adjust the coordinates of the interpolation points,
              call the driver's interpolation function,
              and then adjust the interpolated tensor components,
              but only for the selected faces.
              When no faces are selected, forward the call to the
              driver's interpolator throug the aliased function
              "DriverInterpolate".
              All arguments except "faces" are equivalent to the
              CCTK_InterpGridArrays function call.
  @enddesc
  @var        cctkGH
  @vtype      CCTK_POINTER_TO_CONST
  @vdesc      Grid hierarchy
  @vio        in
  @endvar
  @var        N_dims
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        local_interp_handle
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        param_table_handle
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        coord_system_handle
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        N_interp_points
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        interp_coords_type
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        interp_coords
  @vtype      CCTK_POINTER_TO_CONST[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        N_input_arrays
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        input_array_indices
  @vtype      CCTK_INT[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        N_output_arrays
  @vtype      CCTK_INT
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        output_array_types
  @vtype      CCTK_INT[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        output_arrays
  @vtype      CCTK_POINTER[]
  @vdesc      same as for the CCTK_InterpGridArrays call
  @vio        in
  @endvar
  @var        faces
  @vtype      CCTK_INT
  @vdesc      Bit mask selecting which faces' symmetries should be treated.
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              status code
              0 for success
              -91 GH is NULL
              -92 N_dims is not equal to cctkGH->cctk_dim
              -93 faces is NULL
              or the return value from DriverInterpolate
  @endreturndesc
@@*/

CCTK_INT
SymBase_SymmetryInterpolateFaces(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[],
    CCTK_INT const faces) {
  cGH const *restrict const cctkGH = cctkGH_;
  CCTK_INT sym_table;
  CCTK_FPOINTER symmetry_interpolate[100];
  int face;
  int ierr;

  /* Check arguments */

  if (!cctkGH) {
    CCTK_WARN(1, "cctkGH is NULL");
    return -91;
  }

  if (N_dims != cctkGH->cctk_dim) {
    CCTK_WARN(1, "The number of dimensions is not equal to the GH's number of "
                 "dimensions");
    return -92;
  }

  /* Get table */
  sym_table = SymmetryTableHandleForGrid(cctkGH);
  if (sym_table < 0)
    CCTK_WARN(0, "internal error");

  /* Get table entries */
  if (2 * cctkGH->cctk_dim > 100) {
    CCTK_WARN(0, "internal error");
  }
  ierr =
      Util_TableGetFPointerArray(sym_table, 2 * cctkGH->cctk_dim,
                                 symmetry_interpolate, "symmetry_interpolate");
  if (ierr != 2 * cctkGH->cctk_dim) {
    CCTK_WARN(0, "internal error");
  }

  /* Find a face that needs a symmetry applied */
  for (face = 0; face < 2 * N_dims; ++face) {
    if (faces & (1 << face)) {
      break;
    }
  }

  if (face == 2 * N_dims) {
    /* All faces are done */

    /* Call the real interpolator */
    ierr = CCTK_IsFunctionAliased("DriverInterpolate");
    if (!ierr) {
      CCTK_WARN(
          0,
          "The aliased function \"DriverInterpolate\" has not been provided");
    }
    ierr = DriverInterpolate(
        cctkGH, N_dims, local_interp_handle, param_table_handle,
        coord_system_handle, N_interp_points, interp_coords_type, interp_coords,
        N_input_arrays, input_array_indices, N_output_arrays,
        output_array_types, output_arrays);

  } else {
    /* Treat the face */

    /* Recursive call to a symmetry condition */
    if (!symmetry_interpolate[face]) {
      CCTK_WARN(0, "internal error");
    }
    ierr = ((SymmetryInterpolateFPointer)symmetry_interpolate[face])(
        cctkGH, N_dims, local_interp_handle, param_table_handle,
        coord_system_handle, N_interp_points, interp_coords_type, interp_coords,
        N_input_arrays, input_array_indices, N_output_arrays,
        output_array_types, output_arrays, faces);
  }

  return ierr;
}
