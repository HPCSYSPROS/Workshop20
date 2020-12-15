/*@@
  @file      Operator.c
  @date      Tue Apr 15 18:22:45 1997
  @author    Paul Walker
  @desc
             Definition of interpolation operators for regular uniform grids.
  @enddesc

  @history
  @date      Sun Jul 04 1999
  @author    Thomas Radke
  @hdesc     conversion to Cactus 4.0 (copied from pughGetPoints.c)
  @date      Wed 31 Jan 2001
  @author    Thomas Radke
  @hdesc     translation of fortran interpolators into C
  @date      22 Jan 2002
  @author    Jonathan Thornburg
  @hdesc     Move all local-interpolation code from LocalInterp to here
  @endhistory

  @version   $Id: Operator.c 180 2004-05-17 12:28:55Z goodale $
  @@*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_ErrorCodes.h"
#include "util_Table.h"
#include "Interpolate.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_LocalInterp_Operator_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int LocalInterp_InterpLocalUniform (int num_dims,
                                    int table,
                                    /***** coordinate system *****/
                                    const CCTK_REAL coord_origin[],
                                    const CCTK_REAL coord_delta[],
                                    /***** interpolation points *****/
                                    int num_interp_points,
                                    int interp_coords_type_code,
                                    const void *const interp_coords[],
                                    /***** input arrays *****/
                                    int num_input_arrays,
                                    const CCTK_INT input_array_dims[],
                                    const CCTK_INT input_array_type_codes[],
                                    const void *const input_arrays[],
                                    /***** output arrays *****/
                                    int num_output_arrays,
                                    const CCTK_INT output_array_type_codes[],
                                    void *const output_arrays[])
{
  int iterator, retval;
  char key[128];
  CCTK_INT order, type, nelems;


  /* check for invalid arguments */
  if (num_dims < 0 || num_interp_points < 0 ||
      num_input_arrays < 0 || num_output_arrays < 0)
  {
    return (UTIL_ERROR_BAD_INPUT);
  }

  /* this interpolation operator computes one output array per input array */
  if (num_input_arrays != num_output_arrays)
  {
    CCTK_WARN (1, "Number of input arrays must match number of output arrays");
    return (UTIL_ERROR_BAD_INPUT);
  }

  if (interp_coords_type_code != CCTK_VARIABLE_REAL)
  {
    CCTK_WARN (1, "Interpolation coordinates must be of type CCTK_REAL");
    return (UTIL_ERROR_BAD_INPUT);
  }

  /* check if there's anything to do at all */
  if (num_dims == 0 || num_input_arrays == 0)
  {
    return (0);
  }


  /* get the interpolation order from the user-supplied parameter table */
  order = 1;
  if (table >= 0)
  {
    /* loop through all table options */
    for (iterator = Util_TableItCreate (table);
         Util_TableItQueryIsNonNull (iterator) > 0 &&
         Util_TableItQueryKeyValueInfo (iterator, sizeof (key), key, &type,
                                        &nelems) > 0;
         Util_TableItAdvance (iterator))
    {
      if (CCTK_Equals (key, "order"))
      {
        if (type == CCTK_VARIABLE_INT && nelems == 1)
        {
          Util_TableGetInt (table, &order, "order");
        }
        else
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Invalid value for option 'order' in interpolation "
                      "parameter options table "
                      "(must be CCTK_INT scalar value)");
        }
      }
      else if (CCTK_Equals (key, "N_boundary_points_to_omit") ||
               CCTK_Equals (key, "boundary_off_centering_tolerance") ||
               CCTK_Equals (key, "boundary_extrapolation_tolerance") ||
               CCTK_Equals (key, "per_point_status") ||
               CCTK_Equals (key, "local_interpolator_status"))
      {
        /* warn about unsupported options */
        CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Option with key '%s' in interpolation parameter options "
                    "table is not not supported (will be ignored)", key);
      }
      else if (CCTK_Equals (key, "input_array_time_levels"))
      {
        /* silently ignore options which are meant for the global
           interpolator only */
      }
      else
      {
        /* warn about other options */
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Found option with unrecognized key '%s' in interpolation "
                    "parameter options table (will be ignored)", key);
      }
    }
    Util_TableItDestroy (iterator);
  }

  /* call the interpolator function */
  retval = LocalInterp_Interpolate (order, num_interp_points, num_dims,
                                    num_output_arrays, input_array_dims,
                                    (const CCTK_REAL *const *) interp_coords,
                                    coord_origin, coord_delta,
                                    input_array_type_codes, input_arrays,
                                    output_array_type_codes, output_arrays);

  return (retval);
}
