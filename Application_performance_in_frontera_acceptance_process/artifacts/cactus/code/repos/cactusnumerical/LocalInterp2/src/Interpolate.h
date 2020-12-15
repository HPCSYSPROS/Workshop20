/*@@
  @file      Interpolate.h
  @date      22 Jan 2002
  @author    Jonathan Thornburg
  @desc      function prototypes
  @enddesc
  @version   $Header$
  @history
  @date      22 Jan 2002
  @author    Jonathan Thornburg
  @hdesc     created by editing down LocalInterp::pughInterpGH.h
             (originally dated 4 July 1999, by Thomas Radke)
  @endhistory
  @@*/

#ifndef _LOCALINTERP_INTERPOLATE_H_
#define _LOCALINTERP_INTERPOLATE_H_  1

#ifdef __cplusplus
extern "C"
{
#endif

/* prototypes of interpolation operator to be registered */
int LocalInterp2_InterpLocalUniform (int num_dims,
                                     int param_table_handle,
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
                                     void *const output_arrays[]);

/* prototypes of routines used internally */
int LocalInterp2_Interpolate (int order,
                             int num_points,
                             int num_dims,
                             int num_arrays,
                             const CCTK_INT dims[ /* num_dims */ ],
                             const CCTK_REAL *const coords[],
                             const CCTK_REAL origin[ /* num_dims */ ],
                             const CCTK_REAL delta[ /* num_dims */ ],
                             const CCTK_INT in_types[ /* num_arrays */ ],
                             const void *const in_arrays[ /* num_arrays */ ],
                             const CCTK_INT out_types[ /* num_arrays */ ],
                             void *const out_arrays[ /* num_arrays */ ]);

#ifdef __cplusplus
}
#endif

#endif  /* _LOCALINTERP_INTERPOLATE_H_ */
