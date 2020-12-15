 /*@@
   @file      Interp.c
   @date      July 07 1999
   @author    Thomas Radke
   @desc
              Registration and invocation routines for interpolation operators.
   @enddesc

   @history
   @date      July 07 1999
   @author    Thomas Radke
   @hdesc     Just copied from Reduction.c

   @date      Wed Feb 20 17:52:55 CET 2002
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @hdesc     * revise registration routines to pull out common code
                into macro and new static function
              * permute arguments to the CCTKi_* registration functions
              * add CCTK_InterpRegisterOpLocalUniform() to register
                new-API interpolators
   @date      Tue May 11 13:01:51 CEST 2004
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @hdesc     change  CCTK_InterpLocal()  and  CCTK_InterpGV()
              to give level 1 warnings that these APIs are obsolescent
              and will be phased out soon (and to point to their
              replacement APIs ).
   @date      Sat 19 June 2004
   @author    Thomas Radke
   @hdesc     finally removed old interpolator API code
   @endhistory

   @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cGH.h"
#include "cctk_FortranString.h"
#include "StoreHandledData.h"
#include "cctk_Interp.h"
#include "cctk_Comm.h"
#include "cctk_WarnLevel.h"
#include "cctk_ActiveThorns.h"
#include "util_ErrorCodes.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(comm_Interp_c);


/******************************************************************************
 *************************    External Routines   *****************************
 ******************************************************************************/
/* prototypes for external C routines are declared in header cctk_Interp.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL CCTK_FNAME (CCTK_InterpHandle)
                           (int *handle,
                            ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_InterpGridArrays)
                           (int *ierror,
                            const cGH **GH,
                            const int *N_dims,
                            const int *local_interp_handle,
                            const int *param_table_handle,
                            const int *coord_system_handle,
                            const int *N_interp_points,
                              const int *interp_coords_type,
                              const void *const interp_coords[],
                            const int *N_input_arrays,
                              const CCTK_INT input_array_indices[],
                            const int *N_output_arrays,
                              const CCTK_INT output_array_types[],
                              void *const output_arrays[]);
void CCTK_FCALL CCTK_FNAME (CCTK_InterpLocalUniform)
                           (int *ierror,
                            const int *N_dims,
                            const int *operator_handle,
                            const int *param_table_handle,
                            /***** coordinate system *****/
                            const CCTK_REAL coord_origin[],
                            const CCTK_REAL coord_delta[],
                            /***** interpolation points *****/
                            const int *N_interp_points,
                            const int *interp_coords_type_code,
                            const void *const interp_coords[],
                            /***** input arrays *****/
                            const int *N_input_arrays,
                            const CCTK_INT input_array_dims[],
                            const CCTK_INT input_array_type_codes[],
                            const void *const input_arrays[],
                            /***** output arrays *****/
                            const int *N_output_arrays,
                            const CCTK_INT output_array_type_codes[],
                            void *const output_arrays[]);

/******************************************************************************
 *************************    Internal Data Structures ************************
 ******************************************************************************/

/* structure holding the routines for a registered interpolation operator */
typedef struct
{
  const char *thorn_name;
  const char *operator_name;
  cInterpOpLocalUniform interp_op_local_uniform;
} interp_op_t;

/* typedef for a generic operator function pointer */
typedef int (*operator_fn_t) (void);

/******************************************************************************
 *************************    Static Variables   ******************************
 ******************************************************************************/

/* static data: interpolation operator database and counter for registered
                operators */
static cHandledData *interp_operators = NULL;
static int num_interp_operators = 0;

/******************************************************************************
 ****************** Prototypes for Functions Local to this File ***************
 ******************************************************************************/

static int InterpRegisterOpLocal (operator_fn_t operator_fn,
                                  const char *operator_name,
                                  const char *thorn_name,
                                  const char *fn_name);

/******************************************************************************
 ************************* Registration Functions *****************************
 ******************************************************************************/

 /*@@
   @routine    CCTK_NumInterpOperators
   @date       Mon Oct 22 2001
   @author     Gabrielle Allen
   @desc
               The number of interp operators registered
   @enddesc
   @returntype int
   @returndesc
               number of interpolation operators
   @endreturndesc
@@*/
int CCTK_NumInterpOperators (void)
{
  return (num_interp_operators);
}

/******************************************************************************/

 /*@@
   @routine    CCTK_InterpOperatorImplementation
   @date       Mon Oct 22 2001
   @author     Gabrielle Allen
   @desc
               Provide the name of the implementation
               which provides an interpolation operator
   @enddesc
   @returntype const char *
   @returndesc
               implementation name
   @endreturndesc
@@*/
const char *CCTK_InterpOperatorImplementation (int handle)
{
  const char *name;
  interp_op_t *operator;


  operator = Util_GetHandledData (interp_operators, handle);
  name = operator ? CCTK_ThornImplementation (operator->thorn_name) : NULL;

  return (name);
}

/******************************************************************************/

/*@@
  @routine    CCTK_InterpRegisterOpLocalUniform
  @date       Thu Feb 21 16:02:05 CET 2002
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc       Registers a user-specified function as a local uniform
              interpolation operator.
  @enddesc

  @var        operator_fn
  @vdesc      function pointer pointing to the interpolation operator
  @vtype      cInterpOpLocalUniform
  @vio        in
  @endvar

  @var        operator_name
  @vdesc      character-string name identifying the interpolation operator
  @vtype      const char *
  @vio        in
  @endvar

  @var        thorn_name
  @vdesc      character-string name identifying which thorn provides
              the operator being registered
  @vtype      const char *
  @vio        in
  @endvar

  @returntype int
  @returndesc
              the handle for the newly registered operator, or<p>
              -1 NULL pointer was passed as interpolation operator routine<p>
              -2 failed to allocate memory<p>
              -3 interpolation operator by given name already exists
  @endreturndesc

  @history
  @date       Mon 12 Feb 2001
  @author     Thomas Radke
  @hdesc      Original version

  @date       Thu Feb 21 16:03:25 CET 2002
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @hdesc      * move common logic in all interpolator-registration
                functions into new  InterpRegisterOpLocal()  function
              * convert remaining boilerplate which differs from one
                registration function to another, into this macro

  @date       Sat 19 June 2004
  @author     Thomas Radke
  @hdesc      Worked CCTK_INTERP_REGISTER_FN_BODY macro into
              InterpRegisterOpLocal()
  @endhistory
 @@*/
int CCTK_InterpRegisterOpLocalUniform (cInterpOpLocalUniform operator_fn,
                                       const char *operator_name,
                                       const char *thorn_name)
{
  int handle;


  handle = InterpRegisterOpLocal ((operator_fn_t) operator_fn, operator_name,
                                  thorn_name,
                                  "CCTK_InterpRegisterOpLocalUniform");

  return (handle);
}

/******************************************************************************/

/*@@
  @routine    InterpRegisterOpLocal
  @date       Thu Feb 21 14:41:35 CET 2002
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc       This is an internal worker routine used as part of the
              process of registering an interpolation operator.  It
              gets the  interpolator handle and the  interp_op_t  structure
              in which the operator information will be stored.

              If some interpolation operator is already registered
              under the specified operator name, we use the existing
               interp_op_t  structure .  Otherwise, we allocate a new
              interp_op_t  structure and set it up (thorn, implementation,
              and operator names assigned from this function's arguments,
              all operator pointers set to NULL), then use it.
  @enddesc

  @var        operator_fn
  @vdesc      function pointer pointing to the interpolation operator
  @vtype      operator_fn_t
  @vio        in
  @endvar

  @var        operator_name
  @vdesc      character-string name identifying the interpolation operator
  @vtype      const char *
  @vio        in
  @endvar

  @var        thorn_name
  @vdesc      character-string name identifying which thorn provides
              the operator being registered
  @vtype      const char *
  @vio        in
  @endvar

  @var        fn_name
  @vdesc      character-string name identifying the calling registration
              routine (used to choose the corresponding function pointer
              within the  interp_op_t  structure)
  @vtype      const char *
  @vio        in
  @endvar

  @returntype int
              the handle for the newly registered operator, or<p>
              -1 NULL pointer was passed as interpolation operator routine<p>
              -2 failed to allocate memory<p>
              -3 interpolation operator by given name already exists
  @endreturndesc
 @@*/
static int InterpRegisterOpLocal (operator_fn_t operator_fn,
                                  const char *operator_name,
                                  const char *thorn_name,
                                  const char *fn_name)
{
  int handle;
  interp_op_t *operator;


  if (! operator_fn)
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "%s (called from thorn '%s'):\n"
                "   NULL function pointer passed for interpolation operator "
                "'%s' !", fn_name, thorn_name, operator_name);
    return (-1);
  }

  /* has some operator already been registered under this operator name? */
  handle = Util_GetHandle (interp_operators, operator_name, (void **)&operator);
  if (handle < 0)
  {
    /* no ==> set up a new  interp_op_t  structure for the registration */
    /*        be sure to nullify all function pointers */
    operator = calloc (1, sizeof (*operator));
    if (! operator)
    {
      return (-2);
    }

    operator->thorn_name    = thorn_name;
    operator->operator_name = operator_name;

    handle = Util_NewHandle (&interp_operators, operator_name, operator);
    num_interp_operators++;
  }

  /* register the operator only if it wasn't already */
  if (! strcmp (fn_name, "CCTK_InterpRegisterOpLocalUniform") &&
      ! operator->interp_op_local_uniform)
  {
    operator->interp_op_local_uniform = (cInterpOpLocalUniform) operator_fn;
  }
#if 0
  else if (! strcmp (fn_name, "CCTK_InterpRegisterOpLocalNonUniform") &&
           ! operator->interp_op_local_nonuniform)
  {
    operator->interp_op_local_nonuniform = (cInterpOpLocalNonUniform) operator_fn;
  }
#endif
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "%s (called from thorn '%s'):\n"
                "   Ignoring attempt to register operator '%s' because it has "
                "already been registered by thorn '%s' !",
                fn_name, thorn_name, operator_name, operator->thorn_name);
    handle = -3;
  }

  return (handle);
}

/******************************************************************************
 ************* User Functions to Get Interpolator Handle/Name/etc *************
 ******************************************************************************/

 /*@@
   @routine    CCTK_InterpHandle
   @date       July 07 1999
   @author     Thomas Radke
   @desc
               Returns the handle of a given interpolation operator
   @enddesc
   @var        name
   @vdesc      String containing name of interpolation operator
   @vtype      const char *
   @vio        in
   @vcomment
   @endvar

   @returntype int
   @returndesc
               the handle for the operator registered by this name
               or negative otherwise
   @endreturndesc
@@*/
int CCTK_InterpHandle (const char *name)
{
  int handle;


  handle = Util_GetHandle (interp_operators, name, NULL);
  if (handle < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "No handle found for interpolation operator '%s'", name);
  }

  return (handle);
}

/******************************************************************************/

void CCTK_FCALL CCTK_FNAME (CCTK_InterpHandle)
                           (int *handle, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (name)
  *handle = CCTK_InterpHandle (name);
  free (name);
}

/******************************************************************************/

 /*@@
   @routine    CCTK_InterpOperator
   @date       December 27 2001
   @author     Gabrielle Allen
   @desc
               Returns the name of a interpolation operator
   @enddesc
   @var        handle
   @vdesc      Handle for interpolation operator
   @vtype      int
   @vio        in
   @vcomment
   @endvar

   @returntype const char *
   @returndesc
   The name of the interpolation operator, or NULL if the handle is invalid
   @endreturndesc
@@*/
const char *CCTK_InterpOperator (int handle)
{
  interp_op_t *operator;


  operator = Util_GetHandledData (interp_operators, handle);
  if (! operator)
  {
    CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                "CCTK_InterpHandle: Handle %d invalid", handle);
  }

  return (operator ? operator->operator_name : NULL);
}


/******************************************************************************
 ****************** User Functions to Do Interpolation ************************
 ******************************************************************************/

 /*@@
   @routine    CCTK_InterpGridArrays
   @date       Mon 16 Dec 2002
   @author     Thomas Radke
   @desc
               The general CCTK interpolation routine for grid variables

               Here only the fortran wrapper is defined which calls the
               the C routine CCTK_InterpGridArrays(). This is an overloadable
               routine defined in src/comm/OverloadComm.c.
   @enddesc
   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        N_dims
   @vdesc      (reference to) number of dimensions for the interpolation
   @vtype      const int *
   @vio        in
   @endvar
   @endvar
   @var        local_interp_handle
   @vdesc      (reference to) the handle which specifies the local interpolator
               to use
   @vtype      const int *
   @vio        in
   @endvar
   @var        param_table_handle
   @vdesc      (reference to) the parameter table handle for passing optional
               parameters to the interpolator routine
   @vtype      const int *
   @vio        in
   @endvar
   @var        coord_system_handle
   @vdesc      (reference to) the handle for the underlying coordinate system
   @vtype      const int *
   @vio        in
   @endvar
   @var        N_interp_points
   @vdesc      (reference to) the number of points to interpolate at
   @vtype      const int *
   @vio        in
   @endvar
   @var        interp_coords_type
   @vdesc      (reference to) the CCTK datatype of the coordinate arrays as
               passed via <interp_coords> (common datatype for all arrays)
   @vtype      const int *
   @vio        in
   @endvar
   @var        interp_coords
   @vdesc      list of <N_dims> arrays with coordinate for <N_interp_points>
               points to interpolate at
   @vtype      const void *const []
   @vio        in
   @endvar
   @var        N_input_arrays
   @vdesc      (reference to) the number of input arrays
   @vtype      const int *
   @vio        in
   @endvar
   @var        input_array_indices
   @vdesc      list of <N_input_arrays> grid variables (given by their indices)
               to interpolate
   @vtype      const CCTK_INT []
   @vio        in
   @endvar
   @var        N_output_arrays
   @vdesc      (reference to) the number of output arrays
   @vtype      const int *
   @vio        in
   @endvar
   @var        out_array_types
   @vdesc      list of <N_output_arrays> requested CCTK datatypes for the
               output arrays
   @vtype      const CCTK_INT []
   @vio        in
   @endvar
   @var        output_arrays
   @vdesc      list of <N_output_arrays> output arrays (given by their pointers)
               which receive the interpolation results
   @vtype      void *const []
   @vio        out
   @endvar

   @returntype int
   @returndesc
               return code from routine which overloades CCTK_InterpGridArrays()
               The return code is passed back in <ierror>.
   @endreturndesc
@@*/
void CCTK_FCALL CCTK_FNAME (CCTK_InterpGridArrays)
                           (int *ierror,
                            const cGH **GH,
                            const int *N_dims,
                            const int *local_interp_handle,
                            const int *param_table_handle,
                            const int *coord_system_handle,
                            const int *N_interp_points,
                              const int *interp_coords_type,
                              const void *const interp_coords[],
                            const int *N_input_arrays,
                              const CCTK_INT input_array_indices[],
                            const int *N_output_arrays,
                              const CCTK_INT output_array_types[],
                              void *const output_arrays[])
{
  *ierror = CCTK_InterpGridArrays (*GH, *N_dims, *local_interp_handle,
                                   *param_table_handle, *coord_system_handle,
                                   *N_interp_points, *interp_coords_type,
                                   interp_coords,
                                   *N_input_arrays,input_array_indices,
                                   *N_output_arrays, output_array_types,
                                   output_arrays);
}


/******************************************************************************/

/*@@
  @routine    CCTK_InterpLocalUniform
  @date       22 Oct 2001
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc
        This API interpolates a set of data arrays defined on a uniform
        N-dimensional tensor-product grid, to a set of interpolation points.
        A key-value table is used to pass options to the interpolation
        operator.
  @enddesc

  ***** misc arguments *****

  @var          N_dims
  @vdesc        dimensionality of the interpolation
  @vtype        int N_dims                      (must be >= 1)
  @endvar

  @var          operator_handle
  @vdesc        handle to the interpolation operator
  @vtype        int operator_handle             (must be >= 0)
  @endvar

  @var          param_table_handle
  @vdesc        handle to a key-value table giving additonal parameters
                for the interpolation
  @vtype        int param_table_handle          (must be >= 0)
  @endvar

  ***** arguments specifying the coordinate system *****

  @var          coord_origin
  @vdesc        (pointer to) array[N_dims] of values giving the
                x,y,z,... coordinates which correspond to the
                integer input-array subscripts (0,0,0,...,0)
                (note there is no implication here that such a
                grid point need actually exist; the arrays just
                give the coordinates it would have if it did exist)
  @vtype        const CCTK_REAL coord_origin[N_dims]
  @endvar

  @var          coord_delta
  @vdesc        (pointer to) array[N_dims] of values giving the
                coordinate spacing of the grid
  @vtype        const CCTK_REAL coord_delta[N_dims]
  @endvar

  ***** arguments specifying the interpolation points *****

  @var          N_interp_points
  @vdesc number of interpolation points
  @vtype        int N_interp_points             (must be >= 0)
  @endvar

  @var          interp_coords_type_code
  @vdesc        one of the CCTK_VARIABLE_* codes giving the data
                type of the arrays pointed to by  interp_coords[]
  @vtype        int
  @endvar

  @var          interp_coords
  @vdesc        (pointer to) array[N_dims] of pointers
                to arrays[N_interp_points] giving
                x,y,z,... coordinates of interpolation points
  @vtype        const void *const interp_coords[N_dims]
  @endvar

  ***** arguments specifying the inputs (the data to be interpolated) *****

  @var          N_input_arrays
  @vdesc        number of arrays input to the interpolation
  @vtype        int N_input_arrays              (must be >= 0)
  @endvar

  @var          input_array_dims
  @vdesc        dimensions of the input arrays: unless overridden by
                entries in the parameter table, all input arrays are
                taken to have these dimensions, with [0] the most contiguous
                axis and [N_dims-1] the least contiguous axis, and
                array subscripts in the range 0 <= subscript < dims[axis]
  @vtype        const int input_array_dims[N_dims]
  @endvar

  @var          input_array_type_codes
  @vdesc        (pointer to) array of CCTK_VARIABLE_* codes
                giving the data types of the input arrays
  @vtype        const int input_array_type_codes[N_input_arrays]
  @endvar

  @var          input_arrays
  @vdesc        (pointer to) array[N_input_arrays] of pointers to input arrays
  @vtype        const void *const input_arrays[N_input_arrays]
  @endvar

  ***** arguments specifying the outputs (the interpolation results) *****

  @var          N_output_arrays
  @vdesc        number of arrays output from the interpolation
  @vtype        int N_output_arrays             (must be >= 0)
  @endvar

  @var          output_array_type_codes
  @vdesc        (pointer to) array of CCTK_VARIABLE_* codes
                giving the data types of the output arrays
  @vtype        const int output_array_type_codes[N_output_arrays]
  @endvar

  @var          output_arrays
  @vdesc        (pointer to) array[N_output_arrays] of pointers to output arrays
  @vtype        void *const output_arrays[N_output_arrays]
  @endvar

  ***** return result *****

  @returntype   int
  @returndesc    0 ==> successful, otherwise
                 UTIL_ERROR_BAD_HANDLE ==> bad operator_handle
  @endreturndesc
  @@*/
int CCTK_InterpLocalUniform(int N_dims,
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
                            void *const output_arrays[])
{
  int retval;
  const interp_op_t *operator;
  cInterpOpLocalUniform operator_fn;


  operator = Util_GetHandledData (interp_operators, operator_handle);
  operator_fn = operator ? operator->interp_op_local_uniform : NULL;

  if (operator_fn)
  {
    retval = operator_fn (N_dims, param_table_handle,
                          /*** coordinate system ***/
                          coord_origin, coord_delta,
                          /*** interpolation points ***/
                          N_interp_points, interp_coords_type_code,
                          interp_coords,
                          /***** input arrays *****/
                          N_input_arrays, input_array_dims,
                          input_array_type_codes, input_arrays,
                          /***** output arrays *****/
                          N_output_arrays, output_array_type_codes,
                          output_arrays);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_InterpLocalUniform(): no interpolation operator is "
                "registered under handle %d (did you activate LocalInterp "
                "or some other thorn providing this interpolation operator ?)",
                operator_handle);
    retval = UTIL_ERROR_BAD_HANDLE;
  }

  return (retval);
}


void CCTK_FCALL CCTK_FNAME (CCTK_InterpLocalUniform)
                           (int *ierror,
                            const int *N_dims,
                            const int *operator_handle,
                            const int *param_table_handle,
                            /***** coordinate system *****/
                            const CCTK_REAL coord_origin[],
                            const CCTK_REAL coord_delta[],
                            /***** interpolation points *****/
                            const int *N_interp_points,
                            const int *interp_coords_type_code,
                            const void *const interp_coords[],
                            /***** input arrays *****/
                            const int *N_input_arrays,
                            const CCTK_INT input_array_dims[],
                            const CCTK_INT input_array_type_codes[],
                            const void *const input_arrays[],
                            /***** output arrays *****/
                            const int *N_output_arrays,
                            const CCTK_INT output_array_type_codes[],
                            void *const output_arrays[])
{
  *ierror = CCTK_InterpLocalUniform (*N_dims, *operator_handle,
                                     *param_table_handle,
                                     coord_origin, coord_delta,
                                     *N_interp_points, *interp_coords_type_code,
                                     interp_coords,
                                     *N_input_arrays, input_array_dims,
                                     input_array_type_codes, input_arrays,
                                     *N_output_arrays, output_array_type_codes,
                                     output_arrays);
}
