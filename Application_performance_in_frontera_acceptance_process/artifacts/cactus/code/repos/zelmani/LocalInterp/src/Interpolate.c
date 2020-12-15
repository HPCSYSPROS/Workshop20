/*@@
  @file      Interpolate.c
  @date      Wed 17 Jan 2001
  @author    Thomas Radke
  @desc
             Interpolation of arrays to arbitrary points

             This interpolator is based on the Cactus 3.x Fortran version
             written by Paul Walker.  It also contains some nice optimization
             features from Erik Schnetter.  Jonathan Thornburg added some
             additional comments in October 2001.
  @enddesc

  @history
  @date      Wed 17 Jan 2001
  @author    Thomas Radke
  @hdesc     Translation from Fortran to C
  @date      Thu 18 Oct 2001
  @author    Jonathan Thornburg
  @hdesc     Add lots of comments, LOCALINTERP_VERBOSE_DEBUG debugging code
  @date      22 Jan 2002
  @author    Jonathan Thornburg
  @hdesc     Move all local-interpolation code from LocalInterp to here
  @endhistory

  @version   $$
  @@*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "util_ErrorCodes.h"		/* defines error codes */
#include "cctk.h"
#include "cctk_Interp.h"		/* defines error codes */
#include "cctk_Parameters.h"
#include "Interpolate.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_LocalInterp_Interpolate_c)

#ifdef LOCALINTERP_VERBOSE_DEBUG
  /* if this is >= 0, we print verbose debugging information at this point */
  int LocalInterp_verbose_debug_n = -1;
#endif

/* the highest order of interpolation we support so far */
#define MAXORDER  3

/* the highest dimension for variables we can deal with (so far) */
#define MAXDIM    3

/* we return this if we are successful */
#define RETURN_SUCCESS	0

/******************************************************************************/

/*@@
  @routine INTERPOLATE (macro)
  @date    18 Oct 2001
  @author  code by ???, these comments by Jonathan Thornburg
  @desc
       This macro does the interpolation of in_array[] to compute
       a single value  out_array[n]  (actually  out_array[n]subpart ;
       see the comments below for details).  The data to be interpolated
       must be real numbers of some type, i.e. if the arrays are
       complex this macro must be called separately for the real
       and imaginary parts.
  @enddesc

  @var    cctk_type
  @vdesc  C type of input and output array elements (might be complex)
  @endvar

  @var    cctk_subtype
  @vdesc  C type of actual numbers being interpolated (must be real)
  @endvar

  @var    subpart
  @vdesc  string to be suffixed to input/output array element to get
          to get real number, i.e. empty string if cctk_type is real,
          .Re or .Im as appropriate if cctk_type is complex
  @endvar

  @var    in_array
  @vdesc  A pointer to array to be interpolated (strictly speaking, to
          the array's [0][0]...[0] element); this is typically passed
          as a  void *  pointer so we typecast it as necessary.
  @endvar

  @var    out_array
  @vdesc  A 1-dimensional array where the interpolation result should be
          stored; this is typically passed as a  void *  pointer so we
          typecast it as necessary.
  @endvar

  @var    order
  @vdesc  The order of the interpolation (1=linear, 2=quadratic, 3=cubic, ...)
  @endvar

  @var    point
  @vdesc  [MAXDIM] array of integers giving the integer grid coordinates
          of the closest grid point to the interpolation point; the
          interpolation stencil/molecule is centered at this point.
  @endvar

  @var    dims
  @vdesc  [MAXDIM] array of integers giving the dimensions of  in_array .
  @endvar

  @var    n
  @vdesc  Position in  out_array  where we should store the interpolation
          result.
  @endvar

  @var    coeff
  @vdesc  [MAXDIM][MAX_ORDER+1] array of (floating-point) interpolation
          coefficients; detailed semantics are that coeff[axis][m] is the
          coefficient of y[m] when the 1-dimensional Lagrange interpolation
          polynomial passing through the  order+1  points
             {(0,y[0]), (1,y[1]), ..., (order,y[order])}
          is evaluated at the position x=offset[axis].
  @endvar
  @@*/
/*
 * The basic idea here is that conceptually we first interpolate the
 * (say) 3D gridfn in the x direction at each y and z grid point,
 * then interpolate that 2D plane of values in the y direction at
 * each z grid point, and finally interpolate that 1D line of values
 * in the z direction.  The implementation actually interleaves the
 * different directions' interpolations so that only 3 scalar temporaries
 * are needed.
 */
#define INTERPOLATE(cctk_type, in_array, out_array,                           \
                    order, point, dims, n, coeff)                             \
    {                                                                         \
      int ii, jj, kk;                                                         \
      const cctk_type *fi;                                                    \
      cctk_type interp_result, fj, fk;                                        \
                                                                              \
                                                                              \
      interp_result = 0;                                                      \
                                                                              \
      /* NOTE-MAXDIM: support >3D arrays by adding more loops */              \
      for (kk = 0; kk <= order; kk++)                                         \
      {                                                                       \
        fk = 0;                                                               \
        for (jj = 0; jj <= order; jj++)                                       \
        {                                                                     \
          /* NOTE-MAXDIM: for >3D arrays adapt the index calculation here */  \
          fi = (const cctk_type *) in_array +                                 \
               point[0] + dims[0]*(point[1]+jj + dims[1]*(point[2]+kk));      \
                                                                              \
          fj = 0;                                                             \
          for (ii = 0; ii <= order; ii++)                                     \
          {                                                                   \
            fj += fi[ii] * coeff[0][ii];                                      \
          }                                                                   \
          /* at this point we have just computed */                           \
          /* fj = in_array[*][jj][kk] interpolated to x=offset[0] */          \
                                                                              \
          fk += fj * coeff[1][jj];                                            \
        }                                                                     \
        /* at this point we have just computed */                             \
        /* fk = fj[*][kk] interpolated to y=offset[1] */                      \
        /*       = in_array[*][*][kk] interpolated to */                      \
        /*                            x=offset[0], y=offset[1] */             \
                                                                              \
        interp_result += fk * coeff[2][kk];                                   \
      }                                                                       \
      /* at this point we have just computed */                               \
      /* interp_result = fk[*] interpolated to z=offset[2] */                 \
      /*               = in_array[*][*][*] interpolated to */                 \
      /*                 x=offset[0], y=offset[1], z=offset[2] */             \
                                                                              \
      /* assign the result */                                                 \
      ((cctk_type *) out_array)[n] = interp_result;                           \
    }                                                      /* end of macro */

/******************************************************************************/

/*@@
  @routine    LocalInterp_Interpolate
  @date       Wed 17 Jan 2001
  @author     Thomas Radke
  @desc
              This routine interpolates a set of input arrays
              to a set of output arrays (one-to-one) at arbitrary points
              which are given by their coordinates and the underlying
              regular, uniform grid.

              Current limitations of this implementation are:
              - arrays up to three (MAXDIM) dimensions only can be handled
              - interpolation orders up to three (MAXORDER) only are supported
              - coordinates must be given as CCTK_REAL types
              - input and output array types must be the same
                (no type casting of interpolation results supported)

              Despite of these limitations, the code was programmed almost
              generically in that it can easily be extended to support
              higher-dimensional arrays or more interpolation orders.
              Places where the code would need to be changed to do this,
              are marked with NOTE-MAXDIM and/or NOTE-MAXORDER comments
              as appropriate.
  @enddesc

  @var        num_points
  @vdesc      number of points to interpolate at
  @vtype      int
  @vio        in
  @endvar
  @var        num_dims
  @vdesc      dimensionality of the input arrays
  @vtype      int
  @vio        in
  @endvar
  @var        num_arrays
  @vdesc      number of input/output arrays
  @vtype      int
  @vio        in
  @endvar
  @var        dims
  @vdesc      dimensions of the input arrays
  @vtype      CCTK_INT[ num_dims ]
  @vio        in
  @endvar
  @var        coord
  @vdesc      list of coordinates to interpolate at
  @vtype      CCTK_REAL coord[ num_dims ][ num_points ]
  @vio        in
  @endvar
  @var        origin
  @vdesc      origin of the underlying grid
  @vtype      CCTK_REAL origin[ num_dims ]
  @vio        in
  @endvar
  @var        delta
  @vdesc      deltas of the underlying grid
  @vtype      CCTK_REAL delta[ num_dims ]
  @vio        in
  @endvar
  @var        in_types
  @vdesc      CCTK variable types of input arrays
  @vtype      CCTK_INT in_types[ num_arrays ]
  @vio        in
  @endvar
  @var        in_arrays
  @vdesc      list of input arrays
  @vtype      void *in_arrays[ num_arrays ]
  @vio        in
  @endvar
  @var        out_types
  @vdesc      CCTK variable types of output arrays
  @vtype      CCTK_INT out_types[ num_arrays ]
  @vio        in
  @endvar
  @var        out_arrays
  @vdesc      list of output arrays
  @vtype      void *out_arrays[ num_arrays ]
  @vio        out
  @endvar

  @returntype int
  @returndesc
              0  - successful interpolation
              negative in case of any errors (eg. the negative total number
              of out-of-bounds interpolation points)
  @endreturndesc
  @@*/
int LocalInterp_Interpolate (int order,
                            int num_points,
                            int num_dims,
                            int num_arrays,
                            const CCTK_INT dims[],
                            const CCTK_REAL *const coord[],
                            const CCTK_REAL origin[],
                            const CCTK_REAL delta[],
                            const CCTK_INT in_types[],
                            const void *const in_arrays[],
                            const CCTK_INT out_types[],
                            void *const out_arrays[])
{
  int retval;
  CCTK_REAL delta_inv[MAXDIM];


  retval = RETURN_SUCCESS;

  /*
   * verify parameters and check against our restrictions
   */
  if (num_dims < 1)
  {
    CCTK_WARN (1, "Number of dimensions must be positive");
    return UTIL_ERROR_BAD_INPUT;
  }

  if (num_dims > MAXDIM)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Interpolation of %d-dimensional arrays not implemented",
                num_dims);
    return UTIL_ERROR_BAD_INPUT;
  }

  if (order < 0)
  {
    CCTK_WARN (1, "Interpolation order must be non-negative");
    return UTIL_ERROR_BAD_INPUT;
  }

  if (order > MAXORDER)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Interpolation order %d not implemented", order);
    return UTIL_ERROR_BAD_INPUT;
  }

  /* check that the stencil/molecule isn't bigger than the grid */
  for (int i = 0; i < num_dims; i++)
  {
    if (order+1 > dims[i])      /* stencil/molecule size = order+1 */
    {
      return CCTK_ERROR_INTERP_GRID_TOO_SMALL;
    }
  }

  if (num_points < 0)
  {
    CCTK_WARN (1, "Negative number of points given");
    return UTIL_ERROR_BAD_INPUT;
  }


  /* also immediately return if there's nothing to do */
  if (num_points == 0)
  {
    return RETURN_SUCCESS;
  }

  /* avoid divisions by delta later on */
  for (int i = 0; i < num_dims; i++)
  {
    delta_inv[i] = 1.0 / delta[i];
  }


#pragma omp parallel
  {
  int shift, myretval;
  int out_of_bounds;
  CCTK_INT max_dims[MAXDIM], point[MAXDIM];
  CCTK_REAL below[MAXDIM];
  CCTK_REAL offset[MAXDIM];
  CCTK_REAL coeff[MAXDIM][MAXORDER + 1];
  // HACK: this assumes that there is a point at the origin
  CCTK_INT iorigin[MAXDIM];
  for(int i = 0 ; i < MAXDIM; i++)
    iorigin[i] = floor ((0.-origin[i]) * delta_inv[i]);
  assert(order == 1); // otherwise the "below" point moves

  /* duplicate the dims[] vector into one with MAXDIM-1 elements
     (with the remaining elements zeroed out)
     so that we can use nested loops over MAXDIM dimensions later on */
  memset (max_dims, 0, sizeof (max_dims));
  memcpy (max_dims, dims, (num_dims - 1) * sizeof (*max_dims));

  /* zero out the coefficients and set the elements with index 'order' to one
     so that we can use nested loops over MAXDIM dimensions later on */
  memset (coeff, 0, sizeof (coeff));
  for (int i = num_dims; i < MAXDIM; i++)
  {
    coeff[i][0] = 1;
  }

  /* zero out the iterator */
  memset (point, 0, sizeof (point));

  myretval = RETURN_SUCCESS;

  /* loop over all points to interpolate at */
  #pragma omp for
  for (int n = 0; n < num_points; n++)
  {
    /* reset the out-of-bounds flag */
    out_of_bounds = 0;

    /* loop over all dimensions */
    for (int i = 0; i < num_dims; i++)
    {
      /* closest grid point for stencil/molecule */
      point[i] = iorigin[i] + floor (coord[i][n] * delta_inv[i]
                        - 0.5 * (order - 1));

      /* test bounds */
      out_of_bounds |= point[i] < 0 || point[i]+order >= dims[i];

      /* if beyond lower bound shift the grid point to the right */
      shift = point[i];
      if (shift < 0)
      {
        point[i] -= shift;
      }

      /* if beyond upper bound shift the grid point to the left */
      shift = point[i] + order - (dims[i] - 1);
      if (shift > 0)
      {
        point[i] -= shift;
      }

      /* physical coordinate of that grid point */
      below[i] = (point[i] - iorigin[i]) * delta[i];

      /* offset from that grid point, in fractions of grid points */
      offset[i] = (coord[i][n] - below[i]) * delta_inv[i];
    }

#ifdef LOCALINTERP_VERBOSE_DEBUG
if (n == LocalInterp_verbose_debug_n)
#pragma omp critical
        {
        int ii;
        printf("out_of_bounds = %d\n", out_of_bounds);
                for (ii = 0 ; ii < num_dims ; ++ii)
                {
                printf("offset[%d] = %g\n", ii, (double) offset[ii]);
                }
        }
#endif /* LOCALINTERP_VERBOSE_DEBUG */

    /* check bounds */
    if (out_of_bounds)
    {
#pragma omp critical
      if (num_dims == 1)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Interpolation stencil/molecule out of bounds\n"
                    "     (the interpolation point is either outside the grid, "
                    "or inside but too close to a grid boundary)\n"
                    "     for interpolation order %d\n"
                    "     at interpolation point %d with x = (%f)\n"
                    "     and grid min/max = (%f/%f)\n"
                    "(this may be caused by a global interpolation with\n"
                    " driver::ghost_size too small)\n"
                    ,
                    order, n,
                    (double) coord[0][n],
                    (double) origin[0], (double) (origin[0] + (dims[0]-1)*delta[0]));
      }
      else if (num_dims == 2)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Interpolation stencil/molecule out of bounds\n"
                    "     (the interpolation point is either outside the grid, "
                    "or inside but too close to a grid boundary)\n"
                    "     for interpolation order %d\n"
                    "     at interpolation point %d with xy = (%f, %f)\n"
                    "     and grid min/max = (%f/%f, %f/%f)\n"
                    "(this may be caused by a global interpolation with\n"
                    " driver::ghost_size too small)\n"
                    ,
                    order, n,
                    (double) coord[0][n], (double) coord[1][n],
                    (double) origin[0], (double) (origin[0] + (dims[0]-1)*delta[0]),
                    (double) origin[1], (double) (origin[1] + (dims[1]-1)*delta[1]));
      }
      else if (num_dims == 3)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Interpolation stencil/molecule out of bounds\n"
                    "     (the interpolation point is either outside the grid, "
                    "or inside but too close to a grid boundary)\n"
                    "     for interpolation order %d\n"
                    "     at interpolation point %d with xyz = (%f, %f, %f)\n"
                    "     and grid min/max = (%f/%f, %f/%f, %f/%f)\n"
                    "(this may be caused by a global interpolation with\n"
                    " driver::ghost_size too small)\n"
                    ,
                    order, n,
                    (double) coord[0][n], (double) coord[1][n], (double) coord[2][n],
                    (double) origin[0], (double) (origin[0] + (dims[0]-1)*delta[0]),
                    (double) origin[1], (double) (origin[1] + (dims[1]-1)*delta[1]),
                    (double) origin[2], (double) (origin[2] + (dims[2]-1)*delta[2]));
      }
      else
      {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Internal error: %d dimensions aren't supported", num_dims);
      }

      myretval = CCTK_ERROR_INTERP_POINT_OUTSIDE;
      continue;				/* with next interpolation point */
    }

    /*
     * *** compute the interpolation coefficients according to the order ***
     *
     * (Thanks to Erik for formulating these so nicely.)
     *
     * These formulas are "just" the coefficients of the classical
     * Lagrange interpolation polynomials along each dimension.
     * For example, in 1 dimension the unique quadratic passing
     * through the 3 points {(x0,y0), (x1,y1), (x2,y2)} is:
     *    ( x-x1)( x-x2)        ( x-x0)( x-x2)        ( x-x0)( x-x1)
     *    -------------- y0  +  -------------- y1  +  -------------- y2
     *    (x0-x1)(x0-x2)        (x1-x0)(x1-x2)        (x2-x0)(x2-x1)
     * (It's easy to see this: each of the terms is yi if x=xi, or
     * zero if x=any other xj.)  To get the formulas below, just negate
     * each (x-x) factor, and substitute the values xi=i.
     */
    /*
     * NOTE-MAXORDER: support higher interpolation orders by adding the
     *                appropriate coefficients in another else branch
     */
    switch(order)
    {
      case 0:
        /* zeroeth order (copy) 1D interpolation */
        for (int i = 0; i < num_dims; i++)
        {
          coeff[i][0] = 1;
        }
        break;
      case 1:
        /* first order (linear) 1D interpolation */
        for (int i = 0; i < num_dims; i++)
        {
          CCTK_REAL x1 = (point[i]+1-iorigin[i]) * delta[i];
          CCTK_REAL x0 = (point[i]+0-iorigin[i]) * delta[i];
          coeff[i][0] = (x1 - coord[i][n]) * delta_inv[i]; //1 - offset[i];
          coeff[i][1] = (coord[i][n] - x0) * delta_inv[i]; //    offset[i];
        }
        break;
      case 2:
        /* second order (quadratic) 1D interpolation */
        for (int i = 0; i < num_dims; i++)
        {
          coeff[i][0] = (1-offset[i]) * (2-offset[i]) / (  2  *   1 );
          coeff[i][1] = ( -offset[i]) * (2-offset[i]) / (  1  * (-1));
          coeff[i][2] = ( -offset[i]) * (1-offset[i]) / ((-1) * (-2));
        }
        break;
      case 3:
        /* third order (cubic) 1D interpolation */
        for (int i = 0; i < num_dims; i++)
        {
          coeff[i][0] = (1-offset[i]) * (2-offset[i]) * (3-offset[i]) /
                        (  3  *   2  *   1 );
          coeff[i][1] = ( -offset[i]) * (2-offset[i]) * (3-offset[i]) /
                        (  2  *   1  * (-1));
          coeff[i][2] = ( -offset[i]) * (1-offset[i]) * (3-offset[i]) /
                        (  1  * (-1) * (-2));
          coeff[i][3] = ( -offset[i]) * (1-offset[i]) * (2-offset[i]) /
                        ((-1) * (-2) * (-3));
        }
        break;
      default:
#pragma omp critical
        {
        CCTK_VError (__LINE__,__FILE__,CCTK_THORNSTRING,
                     "Implementation error. Unexpected interpolation order %d",
                     order);
        }
        break;
    }

#ifdef LOCALINTERP_VERBOSE_DEBUG
if (n == LocalInterp_verbose_debug_n)
#pragma omp critical
        {
        int ii,mm;
                for (ii = 0 ; ii < num_dims ; ++ii)
                {
                        for (mm = 0 ; mm <= order ; ++mm)
                        {
                        printf("coeff[%d][%d] = %g\n",
                               ii, mm, (double) coeff[ii][mm]);
                        }
                }
        }
#endif /* LOCALINTERP_VERBOSE_DEBUG */

    /* now loop over all arrays to interpolate at the current point */
    for (int a = 0; a < num_arrays; a++)
    {
      /* check for valid input and output array type */
      if (in_types[a] < 0 || out_types[a] < 0)
      {
#pragma omp critical
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Datatype for input and/or output array with index %d "
                    "is invalid", a);
        myretval = UTIL_ERROR_BAD_INPUT;
        continue;
      }

      /* sorry, for now input and output arrays must be of same type */
      if (in_types[a] != out_types[a])
      {
#pragma omp critical
        CCTK_WARN (1, "Type casting of interpolation results not implemented");
        myretval = UTIL_ERROR_BAD_INPUT;
        continue;
      }

      /* skip this array if it's a query call only */
      if (! out_arrays[a])
      {
        continue;
      }

      /* now do the interpolation according to the array type
         we support all kinds of CCTK_REAL* and CCTK_COMPLEX* types here */
      if (in_types[a] == CCTK_VARIABLE_REAL)
      {
        INTERPOLATE (CCTK_REAL, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
      else if (in_types[a] == CCTK_VARIABLE_COMPLEX)
      {
        INTERPOLATE (CCTK_COMPLEX, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
#ifdef CCTK_REAL4
      else if (in_types[a] == CCTK_VARIABLE_REAL4)
      {
        INTERPOLATE (CCTK_REAL4, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
      else if (in_types[a] == CCTK_VARIABLE_COMPLEX8)
      {
        INTERPOLATE (CCTK_COMPLEX8, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
#endif
#ifdef CCTK_REAL8
      else if (in_types[a] == CCTK_VARIABLE_REAL8)
      {
        INTERPOLATE (CCTK_REAL8, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
      else if (in_types[a] == CCTK_VARIABLE_COMPLEX16)
      {
        INTERPOLATE (CCTK_COMPLEX16, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
#endif
#ifdef CCTK_REAL16
      else if (in_types[a] == CCTK_VARIABLE_REAL16)
      {
        INTERPOLATE (CCTK_REAL16, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
      else if (in_types[a] == CCTK_VARIABLE_COMPLEX32)
      {
        INTERPOLATE (CCTK_COMPLEX32, in_arrays[a],
                     out_arrays[a], order, point, max_dims, n, coeff);
      }
#endif
      else
      {
#pragma omp critical
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Unsupported variable type %d", (int) in_types[a]);
      }
    } /* end of loop over all arrays */

  } /* end of loop over all points to interpolate at */
#pragma omp critical
  if (myretval != RETURN_SUCCESS && retval == RETURN_SUCCESS)
  {
    retval = myretval; /* return one of the error conditions that occured */
  }
  } /* end of parallel section */

  /* we're done */
  return (retval);
}
