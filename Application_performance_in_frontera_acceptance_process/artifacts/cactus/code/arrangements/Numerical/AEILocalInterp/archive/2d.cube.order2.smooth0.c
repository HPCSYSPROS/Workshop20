 /*@@
   @file      2d.cube.order2.smooth0.c
   @date      23 Oct 2001
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc
	Generalized interpolation for 2d, hypercube-shaped molecules,
	order=2, smoothing=0.  For details, see the header comments
	for "InterpLocalArrays.c" in this directory.
   @enddesc

   @version   $Id$
 @@*/

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "util_ErrorCodes.h"
#include "cctk.h"
#include "InterpLocalArrays.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(
   CactusPUGH_LocalInterp_GeneralizedPolynomial_2d_cube_order2_smooth0_c
		)

/******************************************************************************/

/*@@
   @routine    LocalInterp_ILA_2d_cube_ord2_sm0
   @date       23 Oct 2001
   @author     Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc
	This function does generalized interpolation of one or more
	2d arrays to arbitrary points.  For details, see the header
	comments for InterpLocalArrays() (in "InterpLocalArrays.c"
	in this same directory).

	This function's arguments are all a subset of those of
	InterpLocalArrays() ; the only difference is that this function
	takes all its arguments explicitly, whereas  InputLocalArrays()
	takes some of them indirectly via a key/value parameter table.

  @returntype   int
  @returndesc	This function's return result is the same as that of
		InterpLocalArrays():
  		 0 ==> successful interpolation
                -1 ==> in case of any errors
  @endreturndesc

  @@*/
int LocalInterp_ILA_2d_cube_ord2_sm0
	(int param_table_handle,
	 const CCTK_REAL coord_system_origin[],
	 const CCTK_REAL grid_spacing[],
	 int N_interp_points,
	 const CCTK_INT interp_coord_type_codes[],
	 const void *const interp_coords[],
	 int N_input_arrays,
	 const CCTK_INT input_array_offsets[],
	 const CCTK_INT input_array_strides[],
	 const CCTK_INT input_array_min_subscripts[],
	 const CCTK_INT input_array_max_subscripts[],
	 const CCTK_INT input_array_type_codes[],
	 const void *const input_arrays[],
	 int N_output_arrays,
	 const CCTK_INT output_array_type_codes[],
	 void *const output_arrays[],
	 const CCTK_INT operand_indices[], const CCTK_INT opcodes[])
{
/*
 * Implementation notes:
 * 
 * The basic outline of this function is as follows:
 *
 *	compute "which derivatives are wanted" flags
 *		for each interpolation point
 *		{
 *		declare all the coefficients
 *		declare all the data-values variables
 *		***fetch*** interpolation point coordinates
 *		compute coefficients for all derivatives which are wanted
 *			for each output array
 *			{
 *			int part;
 *				for (part = 0 ; part <= 1 ; ++part)
 *				{
 *				if (this output array is computed
 *				    using a different input array
 *				    than the previous one || part != 0)
 *				   then ***fetch*** the input array values
 *					            into local variables
 *				  {
 *				fp result;
 *				switch	(opcode)
 *					{
 *				case 0:
 *					result = evaluate the interpolant
 *					break;
 *				case 1:
 *					result = evaluate the interpolant
 *					break;
 *				case ...
 *					}
 *				***store*** result in output array
 *				bool complex_flag = is datatype complex?
 *				if (! complex_flag)
 *				   then break;
 *				  }
 *				}
 *			}
 *		}
 *
 * Here "***fetch***" and "***store***" are all actually switches on
 * the various array datatypes.  For complex datatypes they offset the
 * 1D array position by  part  to handle real/imaginary components of
 * the complex values.
 *
 * At present we do all floating-point computations in type "fp"
 * (typically a typedef for CCTK_REAL), so arrays of higher precision
 * than this will incur extra rounding errors.  In practice these should
 * be negligible compared to the "truncation" interpolation errors.
 */

/*
 * Naming conventions:
 * input, output = 0-origin indices each selecting an input/output array
 * point = 0-origin index selecting an interpolation point
 */

/*
 * these are compile-time constants here; InterpLocalArrays() decoded
 * them and called us (as opposed to another function) based in part
 * on these values
 */
#define N_DIMS		2
#define MOLECULE_SIZE	3

/* layout of axes in N_dims-element arrays */
#define X_AXIS	0
#define Y_AXIS	1

/* input array size, strides, and subscripting computation */
const int stride_i = input_array_strides[X_AXIS];
const int stride_j = input_array_strides[Y_AXIS];
#define SUB2(i,j)	(i*stride_i + j*stride_j)

/* macros used by machine-generated interpolation coefficient expressions */
#define RATIONAL(num,den)	(num/den)

/*
 * compute flags specifying which derivatives are wanted
 */
bool want_I = false;
bool want_dx = false, want_dy = false;
  {
int output;
	for (output = 0 ; output < N_output_arrays ; ++output)
	{
	switch	(opcodes[output])
		{
	case 0:		want_I = true;		break;
	case 1:		want_dx = true;		break;
	case 2:		want_dy = true;		break;
	default:
		CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
			   "Generalized interpolation opcode %d not supported",
			   opcodes[output]);			/*NOTREACHED*/
		return UTIL_ERROR_BAD_INPUT;		/*** ERROR RETURN ***/
		}
	}
  }

/*
 * interpolate at each point
 */
  {
int point;
	for (point = 0 ; point < N_interp_points ; ++point)
	{
	/* declare all the interpolation coefficients */
	#include "coeffs/2d.cube.order2.smooth0.I.dcl.c"
	#include "coeffs/2d.cube.order2.smooth0.dx.dcl.c"
	#include "coeffs/2d.cube.order2.smooth0.dy.dcl.c"

	/* declare all the data-values variables */
	#include "coeffs/2d.cube.size3.data-var.dcl.c"

	/*
	 * ***fetch*** interpolation point coordinates
	 */
	fp interp_coords_fp[N_DIMS];
	int axis;
		for (axis = 0 ; axis < N_DIMS ; ++axis)
		{
		/* pointer to array of interp coords for this axis */
		const void *const interp_coords_ptr = interp_coords[axis];

		switch	(interp_coord_type_codes[axis])
			{
		case CCTK_VARIABLE_REAL:
			  {
			const CCTK_REAL *const interp_coords_ptr_real
				= (const CCTK_REAL *) interp_coords_ptr;
			interp_coords_fp[axis] = interp_coords_ptr_real[point];
			break;
			  }
#ifdef HAVE_CCTK_REAL4
		case CCTK_VARIABLE_REAL4:
			  {
			const CCTK_REAL4 *const interp_coords_ptr_real4
				= (const CCTK_REAL4 *) interp_coords_ptr;
			interp_coords_fp[axis] = interp_coords_ptr_real4[point];
			break;
			  }
#endif
#ifdef HAVE_CCTK_REAL8
		case CCTK_VARIABLE_REAL8:
			  {
			const CCTK_REAL8 *const interp_coords_ptr_real8
				= (const CCTK_REAL8 *) interp_coords_ptr;
			interp_coords_fp[axis] = interp_coords_ptr_real8[point];
			break;
			  }
#endif
		default:
			CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
				   "interp-coords datatype %d not supported",
				   interp_coord_type_codes[axis]);
								/*NOTREACHED*/
			return UTIL_ERROR_BAD_INPUT;	/*** ERROR RETURN ***/
			}
		}

	/*
	 * locate interpolation molecules with respect to the grid,
	 * i.e. compute interp_rel_(x,y)
	 */
	  {
	fp interp_rel_x, interp_rel_y;	/* (x,y) coordinates of interpolation */
					/* point relative to molecule center */
					/* (in units of the grid spacing) */
	const int center_i
		= LocalInterp_molecule_posn(coord_system_origin[X_AXIS],
					   grid_spacing[X_AXIS],
					   input_array_min_subscripts[X_AXIS],
					   input_array_max_subscripts[X_AXIS],
					   MOLECULE_SIZE,
					   interp_coords_fp[X_AXIS],
					   &interp_rel_x,
					   (int *) NULL, (int *) NULL);
	const int center_j
		= LocalInterp_molecule_posn(coord_system_origin[Y_AXIS],
					   grid_spacing[Y_AXIS],
					   input_array_min_subscripts[Y_AXIS],
					   input_array_max_subscripts[Y_AXIS],
					   MOLECULE_SIZE,
					   interp_coords_fp[Y_AXIS],
					   &interp_rel_y,
					   (int *) NULL, (int *) NULL);
	const int center_sub = SUB2(center_i, center_j);

	/*
	 * compute the coefficients for whichever derivatives are wanted
	 * using machine-generated coefficient expressions
	 * ... these expressions are polynomials in (x,y)
	 *     ==> we need these names for the relative coordinates
	 *	   (copying to fresh local variables will likely also
	 *	    give better register allocation, since [xy]_rel
	 *	    had their addresses taken and so probably won't be
	 *	    register-allocated)
	 */
	const fp x = interp_rel_x;
	const fp y = interp_rel_y;
	if (want_I)
	   then {
		#include "coeffs/2d.cube.order2.smooth0.I.coeff.c"
		}
	if (want_dx)
	   then {
		#include "coeffs/2d.cube.order2.smooth0.dx.coeff.c"
		}
	if (want_dy)
	   then {
		#include "coeffs/2d.cube.order2.smooth0.dy.coeff.c"
		}

	/*
	 * compute each output array at this point
	 */
	  {
	int output;
	const void *input_array_ptr = NULL;
		for (output = 0 ; output < N_output_arrays ; ++output)
		{
		const int input = operand_indices[output];
		const int input_offset = input_array_offsets[input];

		/*
		 * for each real/imag part of complex data values
		 * ... for real we'll break out of this loop at the bottom
		 *     after only a single iteration;
		 *     for complex we'll do both iterations
		 */
		int part;
			for (part = 0 ; part <= 1 ; ++part)
			{
			if ( (input_arrays[input] != input_array_ptr)
			     || (part != 0) )
			   then {
				/*
				 * ***fetch*** the input array values
				 *             into local variables
				 */
				input_array_ptr = input_arrays[input];
				switch	(input_array_type_codes[input])
					{
case CCTK_VARIABLE_REAL:
	  {
	const CCTK_REAL *const input_array_ptr_real
		= (const CCTK_REAL *) input_array_ptr;
	#undef DATA
	#define DATA(i,j)	\
		input_array_ptr_real[input_offset + center_sub + SUB2(i,j)]
	#include "coeffs/2d.cube.size3.data-var.assign.c"
	break;
	  }
#ifdef HAVE_CCTK_REAL4
case CCTK_VARIABLE_REAL4:
	  {
	const CCTK_REAL4 *const input_array_ptr_real4
		= (const CCTK_REAL4 *) input_array_ptr;
	#undef DATA
	#define DATA(i,j)	\
		input_array_ptr_real4[input_offset + center_sub + SUB2(i,j)]
	#include "coeffs/2d.cube.size3.data-var.assign.c"
	break;
	  }
#endif
#ifdef HAVE_CCTK_REAL8
case CCTK_VARIABLE_REAL8:
	  {
	const CCTK_REAL8 *const input_array_ptr_real8
		= (const CCTK_REAL8 *) input_array_ptr;
	#undef DATA
	#define DATA(i,j)	\
		input_array_ptr_real8[input_offset + center_sub + SUB2(i,j)]
	#include "coeffs/2d.cube.size3.data-var.assign.c"
	break;
	  }
#endif
case CCTK_VARIABLE_COMPLEX:
	  {
	const CCTK_COMPLEX *const input_array_ptr_complex
		= (const CCTK_COMPLEX *) input_array_ptr;
	#undef DATA
	#define DATA(i,j)	\
   ( (const CCTK_REAL *)	\
     & input_array_ptr_complex[input_offset + center_sub + SUB2(i,j)] ) [part]
	#include "coeffs/2d.cube.size3.data-var.assign.c"
	break;
	  }
#ifdef HAVE_CCTK_COMPLEX8
case CCTK_VARIABLE_COMPLEX8:
	  {
	const CCTK_COMPLEX8 *const input_array_ptr_complex8
		= (const CCTK_COMPLEX8 *) input_array_ptr;
	#undef DATA
	#define DATA(i,j)	\
   ( (const CCTK_REAL4 *)	\
     & input_array_ptr_complex8[input_offset + center_sub + SUB2(i,j)] ) [part]
	#include "coeffs/2d.cube.size3.data-var.assign.c"
	break;
	  }
#endif
#ifdef HAVE_CCTK_COMPLEX16
case CCTK_VARIABLE_COMPLEX16:
	  {
	const CCTK_COMPLEX16 *const input_array_ptr_complex16
		= (const CCTK_COMPLEX16 *) input_array_ptr;
	#undef DATA
	#define DATA(i,j)	\
   ( (const CCTK_REAL8 *)	\
     & input_array_ptr_complex16[input_offset + center_sub + SUB2(i,j)] ) [part]
	#include "coeffs/2d.cube.size3.data-var.assign.c"
	break;
	  }
#endif
default:
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "input datatype %d not supported",
		   input_array_type_codes[input]);		/*NOTREACHED*/
	return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
					}
				}

			/*
			 * evaluate the interpolant
			 * as a linear combination of the data variables
			 */
			  {
			fp result;
			switch	(opcodes[output])
				{
			case 0:
				result =
	#include "coeffs/2d.cube.order2.smooth0.I.eval.c"
				break;
			case 1:
				result =
	#include "coeffs/2d.cube.order2.smooth0.dx.eval.c"
				break;
			case 2:
				result =
	#include "coeffs/2d.cube.order2.smooth0.dy.eval.c"
				break;
			default:
				CCTK_VWarn(1, __LINE__, __FILE__,
					   CCTK_THORNSTRING,
					   "opcode %d not supported",
					   opcodes[output]);	/*NOTREACHED*/
				return UTIL_ERROR_BAD_INPUT;
							/*** ERROR RETURN ***/
				}

			/*
			 * ***store*** the result in the output array
			 */
			switch	(output_array_type_codes[output])
				{
case CCTK_VARIABLE_REAL:
	  {
	CCTK_REAL *const output_array_ptr_real
		= (CCTK_REAL *) output_arrays[output];
	output_array_ptr_real[point] = result;
	break;
	  }
#ifdef HAVE_CCTK_REAL4
case CCTK_VARIABLE_REAL4:
	  {
	CCTK_REAL4 *const output_array_ptr_real4
		= (CCTK_REAL4 *) output_arrays[output];
	output_array_ptr_real4[point] = result;
	break;
	  }
#endif
#ifdef HAVE_CCTK_REAL8
case CCTK_VARIABLE_REAL8:
	  {
	CCTK_REAL8 *const output_array_ptr_real8
		= (CCTK_REAL8 *) output_arrays[output];
	output_array_ptr_real8[point] = result;
	break;
	  }
#endif
#ifdef HAVE_CCTK_REAL16
case CCTK_VARIABLE_REAL16:
	  {
	CCTK_REAL16 *const output_array_ptr_real16
		= (CCTK_REAL16 *) output_arrays[output];
	output_array_ptr_real16[point] = result;
	break;
	  }
#endif
case CCTK_VARIABLE_COMPLEX:
	  {
	CCTK_COMPLEX *const output_array_ptr_complex
		= (CCTK_COMPLEX *) output_arrays[output];
	((CCTK_REAL *)  & output_array_ptr_complex[point]) [part]
		= result;
	break;
	  }
#ifdef HAVE_CCTK_COMPLEX8
case CCTK_VARIABLE_COMPLEX8:
	  {
	CCTK_COMPLEX8 *const output_array_ptr_complex8
		= (CCTK_COMPLEX8 *) output_arrays[output];
	((CCTK_REAL4 *)  & output_array_ptr_complex8[point]) [part]
		= result;
	break;
	  }
#endif
#ifdef HAVE_CCTK_COMPLEX16
case CCTK_VARIABLE_COMPLEX16:
	  {
	CCTK_COMPLEX16 *const output_array_ptr_complex16
		= (CCTK_COMPLEX16 *) output_arrays[output];
	((CCTK_REAL8 *)  & output_array_ptr_complex16[point]) [part]
		= result;
	break;
	  }
#endif
#ifdef HAVE_CCTK_COMPLEX32
case CCTK_VARIABLE_COMPLEX32:
	  {
	CCTK_COMPLEX32 *const output_array_ptr_complex32
		= (CCTK_COMPLEX32 *) output_arrays[output];
	((CCTK_REAL16 *)  & output_array_ptr_complex32[point]) [part]
		= result;
	break;
	  }
#endif
default:
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "output datatype %d not supported",
		   output_array_type_codes[output]);		/*NOTREACHED*/
	return UTIL_ERROR_BAD_INPUT;			/*** ERROR RETURN ***/
				}

			/* decode datatype: is it real or complex? */
			  {
			bool complex_flag;
			switch	(output_array_type_codes[output])
				{
			case CCTK_VARIABLE_REAL:
			case CCTK_VARIABLE_REAL4:
			case CCTK_VARIABLE_REAL8:
			case CCTK_VARIABLE_REAL16:
				complex_flag = false;
				break;
			case CCTK_VARIABLE_COMPLEX:
			case CCTK_VARIABLE_COMPLEX8:
			case CCTK_VARIABLE_COMPLEX16:
			case CCTK_VARIABLE_COMPLEX32:
				complex_flag = true;
				break;
			default:
				CCTK_VWarn(1, __LINE__, __FILE__,
					   CCTK_THORNSTRING,
					   "output datatype %d not supported",
					   output_array_type_codes[output]);
								/*NOTREACHED*/
				return UTIL_ERROR_BAD_INPUT;
							/*** ERROR RETURN ***/
				}
			
			/* skip part=1 (imaginary part) for real datatypes */
			if (! complex_flag)
			   then break;
			  }
			  }
			/* end of  for (part = ...)  loop */
			}
		/* end of  for (output = ...)  loop */
		}
	  }
	  }

	/* end of  for (point = ...)  loop */
	}

return 0;						/*** NORMAL RETURN ***/
  }
}
