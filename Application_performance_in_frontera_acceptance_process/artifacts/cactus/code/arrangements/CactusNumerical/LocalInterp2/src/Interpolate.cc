#include "cctk.h"
#include "cctk_Interp.h"		/* defines error codes */
#include "util_ErrorCodes.h"		/* defines error codes */
#include "cctk_Parameters.h"
#include "Interpolate.h"

#include "LagrangeInterp.hh"

/* the highest order of interpolation we support so far */
#define MAXORDER  3
/* the highest dimension for variables we can deal with (so far) */
#define MAXDIM    3

/* we return this if we are successful */
#define RETURN_SUCCESS	0

namespace {

#define INTERP_ID(order, num_dims)                                            \
  ((order)*MAXDIM + (num_dims) - 1)

template<int interp_id>
class interp_t {
  public:
    // interp_id = order * MAXDIM + num_dims
    enum { order    = interp_id / MAXDIM };
    enum { num_dims = interp_id % MAXDIM + 1 };
};

template<int order, int num_dims>
int __do_interpolate(
    int num_points,
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
  bool error_point_outside = false;
  bool error_bad_input     = false;
#pragma omp parallel for reduction(||: error_point_outside, error_bad_input)
  for(int n = 0; n < num_points; ++n) {
    CCTK_REAL point[num_dims];
    for(int d = 0; d < num_dims; ++d) {
      point[d] = coord[d][n];
    }

    LagrangeInterpND<order, num_dims> interp(origin, delta, dims, point);

    if(interp.out_of_bounds) {
#pragma omp critical
      {
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
      } // end of critical section
      error_point_outside = true;
      continue;
    } // end of if(interp.out_of_bounds)

    for(int a = 0; a < num_arrays; ++a) {
      /* check for valid input and output array type */
      if (in_types[a] < 0 || out_types[a] < 0)
      {
#pragma omp critical
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Datatype for input and/or output array with index %d "
                    "is invalid", a);
        error_bad_input = true;
        continue;
      }

      /* sorry, for now input and output arrays must be of same type */
      if (in_types[a] != out_types[a])
      {
#pragma omp critical
        CCTK_WARN (1, "Type casting of interpolation results not implemented");
        error_bad_input = true;
        continue;
      }

      /* skip this array if it's a query call only */
      if (! out_arrays[a])
      {
        continue;
      }

      switch(in_types[a]) {
#define TYPECASE(CCTK_TYPE_NAME, CCTK_TYPE)                                   \
        case CCTK_TYPE_NAME:                                                  \
          reinterpret_cast<CCTK_TYPE * const>(out_arrays[a])[n] = interp.     \
            LagrangeInterpND<order, num_dims>::template eval<CCTK_TYPE>(      \
              reinterpret_cast<CCTK_TYPE const * const>(in_arrays[a]));       \
          break
        TYPECASE(CCTK_VARIABLE_REAL, CCTK_REAL);
        TYPECASE(CCTK_VARIABLE_COMPLEX, CCTK_COMPLEX);
#ifdef HAVE_CCTK_REAL4
        TYPECASE(CCTK_VARIABLE_REAL4, CCTK_REAL4);
        TYPECASE(CCTK_VARIABLE_COMPLEX8, CCTK_COMPLEX8);
#endif
#ifdef HAVE_CCTK_REAL8
        TYPECASE(CCTK_VARIABLE_REAL8, CCTK_REAL8);
        TYPECASE(CCTK_VARIABLE_COMPLEX16, CCTK_COMPLEX16);
#endif
#ifdef HAVE_CCTK_REAL16
        TYPECASE(CCTK_VARIABLE_REAL16, CCTK_REAL16);
        TYPECASE(CCTK_VARIABLE_COMPLEX32, CCTK_COMPLEX32);
#endif
#undef TYPECASE
        default:
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Unsupported variable type %d", (int) in_types[a]);
      }
    }
  } // end of parallel for loop
  if(error_bad_input) {
    return UTIL_ERROR_BAD_INPUT;
  }
  else if(error_point_outside) {
    return CCTK_ERROR_INTERP_POINT_OUTSIDE;
  }
  return RETURN_SUCCESS;
}

int interpolate(
    int order,
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

  int id = INTERP_ID(order, num_dims);
  switch(id)
  {
#define TYPECASE(ID)                                                          \
    case ID:                                                                  \
      return                                                                  \
        __do_interpolate<interp_t<ID>::order, interp_t<ID>::num_dims>(        \
            num_points, num_arrays, dims, coord, origin, delta, in_types,     \
            in_arrays, out_types, out_arrays);                                \
      break
    TYPECASE(0);
    TYPECASE(1);
    TYPECASE(2);
    TYPECASE(3);
    TYPECASE(4);
    TYPECASE(5);
    TYPECASE(6);
    TYPECASE(7);
    TYPECASE(8);
    TYPECASE(9);
    TYPECASE(10);
    TYPECASE(11); // MAXDIM * (MAXORDER + 1) - 1
#undef TYPECASE
    default:
      CCTK_ERROR("This is a bug in Interpolate.cc");
  }

  // If we ever get here, this is a bug in the code
  return -1;
}

} // anonymous namespace

/*@@
  @routine    LocalInterp2_Interpolate
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
extern "C"
int LocalInterp2_Interpolate (int order,
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
                              void *const out_arrays[]) {
  return interpolate(order, num_points, num_dims, num_arrays, dims, coord,
      origin, delta, in_types, in_arrays, out_types, out_arrays);
}
