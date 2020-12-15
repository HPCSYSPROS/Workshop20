 /*@@
   @header  cctk_Interp.h
   @date    July 07 1999
   @author  Thomas Radke
   @desc
            Header file for using interpolation operators
   @enddesc

   @history
   @date    July 07 1999
   @author  Thomas Radke
   @hdesc   Just copied from cctk_Reduction.h

   @date    Thu Feb 21 14:36:02 CET 2002
   @author  Jonathan Thornburg <jthorn@aei.mpg.de>
   @hdesc   add more comments, add new stuff for new interpolator API
   @endhistory

   @version $Header$
 @@*/


#ifndef _CCTK_INTERP_H_
#define _CCTK_INTERP_H_

#ifdef CCODE

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * typedefs for interpolation operator routines
 */
typedef int (*cInterpOpLocalUniform) (int N_dims,
                                      int param_table_handle,
                                      /***** coordinate system *****/
                                      const CCTK_REAL coord_origin[],
                                      const CCTK_REAL coord_delta[],
                                      /***** interpolation points *****/
                                      int N_interp_points,
                                      int interp_coords_type_code,
                                      const void *const interp_coords[],
                                      /***** input arrays *****/
                                      int N_input_arrays,
                                      const CCTK_INT input_array_dims[],
                                      const CCTK_INT input_array_type_codes[],
                                      const void *const input_arrays[],
                                      /***** output arrays *****/
                                      int N_output_arrays,
                                      const CCTK_INT output_array_type_codes[],
                                      void *const output_arrays[]);

/*
 * prototypes for user-visible interpolation-registration API
 */
int CCTK_InterpHandle (const char *name);

int CCTK_InterpRegisterOpLocalUniform (cInterpOpLocalUniform operator_ptr,
                                       const char *operator_name,
                                       const char *thorn_name);

const char *CCTK_InterpOperatorImplementation (int handle);
const char *CCTK_InterpOperator (int handle);
int CCTK_NumInterpOperators (void);


/*
 * prototypes for user-visible interpolation API
 */
int CCTK_InterpLocalUniform (int N_dims,
                             int operator_handle,
                             int param_table_handle,
                             /***** coordinate system *****/
                             const CCTK_REAL coord_origin[],
                             const CCTK_REAL coord_delta[],
                             /***** interpolation points *****/
                             int N_interp_points,
                             int interp_coords_type_code,
                             const void *const interp_coords[],
                             /***** input arrays *****/
                             int N_input_arrays,
                             const CCTK_INT input_array_dims[],
                             const CCTK_INT input_array_type_codes[],
                             const void *const input_arrays[],
                             /***** output arrays *****/
                             int N_output_arrays,
                             const CCTK_INT output_array_type_codes[],
                             void *const output_arrays[]);

#ifdef __cplusplus
}
#endif

#endif /* ifdef CCODE */

/*
 * error codes for CCTK_InterpLocalUniform()
 */

/* the grid is too small for the selected interpolation molecule */
#define CCTK_ERROR_INTERP_GRID_TOO_SMALL (-1000)
/* ... old code for backwards compatability */
#define CCTK_ERROR_INTERP_GRID_TOO_TINY  CCTK_ERROR_INTERP_GRID_TOO_SMALL

/*
 * the (multiprocessor) grid's ghostzone size is too small for the selected
 * interpolation molecule (or this processor's chunk of the grid is too small)
 */
#define CCTK_ERROR_INTERP_GHOST_SIZE_TOO_SMALL (-1001)

/*
 * an interpolation point is outside (or too close to an edge of)
 * the input grid
 */
#define CCTK_ERROR_INTERP_POINT_OUTSIDE (-1002)
/* ... old code for backwards compatability */
#define CCTK_ERROR_INTERP_POINT_X_RANGE CCTK_ERROR_INTERP_POINT_OUTSIDE

/* an interpolation point is in (or too close to) an excised region */
#define CCTK_ERROR_INTERP_POINT_EXCISED (-1003)

/* an interpolation coordinate (or some other intermediate value in the */
/* interpolation process) was an IEEE NaN or other non-finite number */
#define CCTK_ERROR_INTERP_COORD_NAN	(-1004)

/* the grid spacing was specified as zero along at least one axis */
#define CCTK_ERROR_INTERP_DELTA_X_ZERO	(-1005)

#endif  /* _INTERP_H_ */
