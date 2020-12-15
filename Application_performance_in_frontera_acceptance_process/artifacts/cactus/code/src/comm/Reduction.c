 /*@@
   @file      Reduction.c
   @date      1999/05/13
   @author    Gabrielle Allen
   @desc
              This file contains routines to deal with registering and
              using functions providing reduction operations.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_Groups.h"
#include "cctk_Reduction.h"
#include "cctk_WarnLevel.h"
#include "cctk_ActiveThorns.h"

#include "StoreHandledData.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(comm_Reduction_c);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/* prototypes for external C routines are declared in header cctk_Reduction.h
   here only follow the fortran wrapper prototypes */

void CCTK_FCALL CCTK_FNAME(CCTK_ReductionHandle)
     (int *handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(CCTK_ReductionArrayHandle)
     (int *operation_handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(CCTK_Reduce)
     (int *fortranreturn,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const int *num_out_vals,
      const int *type_out_vals,
      void *out_vals,
      const int *num_in_fields,
      ... );
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceArray)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const int *num_out_vals,
      const int *type_out_vals,
      void *out_vals,
      const int *num_dims,
      const int *num_in_arrays,
      const int *type_in_arrays,
      ... );

/* new local array reduction API */
void CCTK_FCALL CCTK_FNAME(CCTK_LocalArrayReductionHandle)
     (int *handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocalArrays)
     (int *fortran_return,
      int *N_dims, int *operation_handle, 
      int *param_table_handle,   int *N_input_arrays,
      const CCTK_INT input_array_dims[], 
      const CCTK_INT input_array_type_codes[],
      const void *const input_arrays[],
      int *M_output_values,
      const CCTK_INT output_number_type_codes[],
      void *const output_values[]);

/* new gridarray reduction API */
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceGridArrays)
     (int *fortranreturn,
      const cGH **GH,
      int *dest_proc,
      int *local_reduce_handle,
      int *param_table_handle,
      int *N_input_arrays,
      const CCTK_INT input_array_variable_indices[],
      int *M_output_values,
      const CCTK_INT output_value_type_codes[],
      void* const output_values[]);

/* new reduce arrays globally API */
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceArraysGlobally)
     (int *fortranreturn,
      const cGH **GH,
      int *dest_proc,
      int *local_reduce_handle,
      int *param_table_handle,
      int *N_input_arrays,
      const void * const input_arrays[],
      int *input_dims,
      const CCTK_INT input_array_dims[],
      const CCTK_INT input_array_type_codes[],      
      int *M_output_values,
      const CCTK_INT output_value_type_codes[],
      void* const output_values[]);
            
/* FIXME: OLD INTERFACE */
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocalScalar)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_scalar,
      void *out_scalar,
      const int *data_type);
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocScalar)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_scalar,
      void *out_scalar,
      const int *data_type);
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocalArray1D)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_array1d,
      void *out_array1d,
      const int *num_in_array1d,
      const int *data_type);
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocArrayToArray1D)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_array1d,
      void *out_array1d,
      const int *num_in_array1d,
      const int *data_type);
void CCTK_FCALL  CCTK_FNAME(CCTK_ReduceLocArrayToArray2D)
     (int  *fortran_return, const cGH **GH,
      const int  *proc,
      const int  *operation_handle,
      const void *in_array2d,
      void *out_array2d,
      const int  *xsize, const int *ysize,
      const int  *data_type);
void CCTK_FCALL  CCTK_FNAME(CCTK_ReduceLocArrayToArray3D)
     (int  *fortran_return, const cGH **GH,
      const int  *proc,
      const int  *operation_handle,
      const void *in_array3d,
      void *out_array3d,
      const int  *xsize, const int *ysize, const int *zsize,
      const int  *data_type);


/********************************************************************
 ********************    Internal Typedefs   ************************
 ********************************************************************/
/* structure holding the routines for a registered reduction operator */
typedef struct
{
  const char *implementation;
  const char *name;
  cReduceOperator reduce_operator;
} t_reduce_operator;

/* structure holding a function pointer to an array reduction operator */
typedef struct
{
  int (*function) (REDUCTION_ARRAY_OPERATOR_REGISTER_ARGLIST);
} t_reduction_array_op;

/* structure holding the routines for a registered local array reduction operator */
typedef struct
{
  const char *implementation;
  const char *name;
  cLocalArrayReduceOperator reduce_operator;
} t_local_array_reduce_operator;

/* structure holding a function pointer to an array reduction operator */
typedef struct
{
  int (*function) (REDUCTION_LOCAL_ARRAY_OPERATOR_REGISTER_ARGLIST);
} t_reduction_local_array_op;

/* structure holding the routines for a registered grid array reduction operator */
typedef struct
{
  const char *implementation;
  cGridArrayReduceOperator reduce_operator;
} t_grid_array_reduce_operator;

/* structure holding a function pointer to a GA reduction operator */
typedef struct
{
  int (*function) (REDUCTION_GRID_ARRAY_OPERATOR_REGISTER_ARGLIST);
} t_reduction_grid_array_op;

/* structure holding the routines for a registered global array reduction operator */
typedef struct
{
  const char *implementation;
  const char *name;
  cReduceArraysGloballyOperator reduce_operator;
} t_reduce_arrays_globally_operator;

/* structure holding a function pointer to a global array reduction operator */
typedef struct
{
  int (*function) (REDUCTION_ARRAYS_GLOBALLY_OPERATOR_REGISTER_ARGLIST);
} t_reduction_arrays_globally_op;

/********************************************************************
 ********************    Static Variables   *************************
 ********************************************************************/
static cHandledData *ReductionOperators = NULL;
static int num_reductions = 0;
static cHandledData *ReductionArrayOperators = NULL;
static int num_reductions_array = 0;

static cHandledData *LocalArrayReductionOperators = NULL;
static int num_local_array_reductions = 0;

static cHandledData *GridArrayReductionOperators = NULL;
static int num_GA_reductions = 0;

static cGridArrayReduceOperator GA_reduc = NULL;
static char global[] ="c";

static cHandledData *ReductionArraysGloballyOperators = NULL;
static int num_reductions_arrays_globally = 0;

static cReduceArraysGloballyOperator ArraysGlobally_reduc = NULL;
static char reduction_arrays_globally_global[] ="Y";

 /*@@
   @routine    CCTKi_RegisterReductionOperator
   @date       April 28 1999
   @author     Gabrielle Allen
   @desc
   Registers "function" as a reduction operator called "name"
   @enddesc
   @var     function
   @vdesc   Routine containing reduction operator
   @vtype   (void (*))
   @vio
   @endvar
   @var     name
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTKi_RegisterReductionOperator(const char *thorn,
                                    cReduceOperator operator,
                                    const char *name)
{
  int handle;
  void *tmp;
  t_reduce_operator *reduce_operator;


  /* Check that the method hasn't already been registered */
  handle = Util_GetHandle(ReductionOperators, name, &tmp);
  reduce_operator = tmp;
  if(handle < 0)
  {
    reduce_operator = malloc (sizeof (t_reduce_operator));
    if (reduce_operator)
    {
      reduce_operator->implementation = CCTK_ThornImplementation(thorn);
      reduce_operator->name = name;
      reduce_operator->reduce_operator = operator;
      handle = Util_NewHandle(&ReductionOperators, name, reduce_operator);

      /* Remember how many reduction operators there are */
      num_reductions++;
    }
  }
  else
  {
    /* Reduction operator with this name already exists. */
    CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
              "CCTK_RegisterReductionOperator: Reduction operator "
              "with this name already exists");
    handle = -1;
  }

  return handle;
}


 /*@@
   @routine    CCTK_ReductionHandle
   @date       April 28 1999
   @author     Gabrielle Allen
   @desc
   Returns the handle of a given reduction operator
   @enddesc
   @var     reduction
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTK_ReductionHandle(const char *reduction)
{
  int handle;


  handle = Util_GetHandle(ReductionOperators, reduction, NULL);
  if (handle < 0)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_ReductionHandle: No handle: '%d' found for reduction operator "
               "'%s'", handle, reduction);
  }

  return handle;
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReductionHandle)(int *handle, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(reduction)
  *handle = CCTK_ReductionHandle(reduction);
  free(reduction);
}


 /*@@
   @routine    CCTK_Reduce
   @date       April 28 1999
   @author     Gabrielle Allen
   @desc
               Generic routine for doing a reduction operation on a set of
               Cactus variables.
   @enddesc
   @var     GH
   @vdesc   pointer to the grid hierarchy
   @vtype   cGH *
   @vio     in
   @endvar
   @var     proc
   @vdesc   processor that receives the result of the reduction operation
            (a negative value means that all processors get the result)
   @vtype   int
   @vio     in
   @endvar
   @var     operation_handle
   @vdesc   the handle specifying the reduction operator
   @vtype   int
   @vio     in
   @endvar
   @var     num_out_vals
   @vdesc   number of elements in the reduction output
   @vtype   int
   @vio     in
   @endvar
   @var     type_out_vals
   @vdesc   datatype of the output values
   @vtype   int
   @vio     in
   @endvar
   @var     out_vals
   @vdesc   pointer to buffer holding the output values
   @vtype   void *
   @vio     in
   @endvar
   @var     num_in_fields
   @vdesc   number of input fields passed in the variable argument list
   @vtype   int
   @vio     in
   @endvar
   @var     <...>
   @vdesc   list of variables indices of input fields
   @vtype   int
   @vio     in
   @endvar
@@*/
int CCTK_Reduce(const cGH *GH,
                int proc,
                int operation_handle,
                int num_out_vals,
                int type_out_vals,
                void *out_vals,
                int num_in_fields,
                ... )
{
  va_list indices;
  int i;
  int retval;
  int *in_fields;
  t_reduce_operator *operator;

  /* Get the pointer to the reduction operator */
  if (operation_handle < 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_Reduce: Invalid handle passed to CCTK_Reduce");
    retval = -1;
  }
  else
  {
    operator = Util_GetHandledData(ReductionOperators,operation_handle);

    if (!operator)
    {
      CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
                "CCTK_Reduce: Reduction operation is not registered"
                "and cannot be called");
      retval = -1;
    }
    else
    {
      /* Fill in the array of variable indices from the variable
         argument list */
      in_fields = malloc(num_in_fields*sizeof(int));
      va_start(indices, num_in_fields);
      for (i=0; i<num_in_fields; i++)
      {
        in_fields[i] = va_arg(indices,int);
      }
      va_end(indices);

      retval = operator->reduce_operator (GH, proc, num_out_vals,
                                          type_out_vals, out_vals,
                                          num_in_fields, in_fields);

      free(in_fields);
    }
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_Reduce)
     (int *fortranreturn,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const int *num_out_vals,
      const int *type_out_vals,
      void *out_vals,
      const int *num_in_fields,
      ... )
{
  va_list indices;
  int retval;
  int i;
  int *in_fields;
  t_reduce_operator *operator;

  /* initialize return code to indicate an error */
  *fortranreturn = -1;

  if (*operation_handle < 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_Reduce: Invalid handle passed to CCTK_Reduce");
    retval = -1;
  }
  else
  {
    /* Get the pointer to the reduction operator */
    operator = Util_GetHandledData(ReductionOperators,*operation_handle);

    if (!operator)
    {
      CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
                "CCTK_Reduce: Reduction operation is not registered"
                " and cannot be called");
      retval = -1;
    }
    else
    {

      /* Fill in the array of variable indices from the variable
         argument list */
      in_fields = malloc (*num_in_fields * sizeof (int));
      va_start(indices, num_in_fields);
      for (i=0; i<*num_in_fields; i++)
      {
        in_fields[i] = *va_arg(indices,int *);
      }
      va_end(indices);

      retval = operator->reduce_operator (*GH, *proc, *num_out_vals,
                                          *type_out_vals, out_vals,
                                          *num_in_fields,in_fields);

      free(in_fields);
    }
  }
  *fortranreturn = retval;
}


 /*@@
   @routine CCTK_RegisterReductionArrayOperator
   @date    Aug 19 1999
   @author  Thomas Radke
   @desc
            Registers "function" as a array reduction operator called "name"
   @enddesc
   @var     function
   @vdesc   Routine containing reduction operator
   @vtype   (int (*))
   @vio
   @endvar
   @var     name
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTK_RegisterReductionArrayOperator
         (int (*function)(REDUCTION_ARRAY_OPERATOR_REGISTER_ARGLIST),
         const char *name)
{
  int handle;
  t_reduction_array_op *data;

  /* Check that the method hasn't already been registered */
  handle = Util_GetHandle(ReductionArrayOperators, name, NULL);

  if(handle < 0)
  {
    data = malloc (sizeof (*data));
    data->function = function;

    /* Get a handle for it. */
    handle = Util_NewHandle(&ReductionArrayOperators, name, data);

    /* Remember how many reduction operators there are */
    num_reductions_array++;

  }
  else
  {
    /* Reduction operator with this name already exists. */
    CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
              "CCTK_RegisterReductionArrayOperator: "
              "Array reduction operator with this name already exists");
    handle = -1;
  }

  return handle;
}


 /*@@
   @routine CCTK_ReductionArrayHandle
   @date    Aug 19 1999
   @author  Thomas Radke
   @desc
            Returns the handle of a given array reduction operator
   @enddesc
   @var     reduction
   @vdesc   String containing name of array reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTK_ReductionArrayHandle(const char *reduction)
{

  int handle;

  handle = Util_GetHandle(ReductionArrayOperators, reduction, NULL);
  if (handle < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_ReductionArrayHandle: "
                "No handle found for array reduction operator '%s'",
                reduction);
  }

  return handle;
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReductionArrayHandle)
     (int *operation_handle, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(reduction)
  *operation_handle = CCTK_ReductionArrayHandle(reduction);
  free(reduction);
}


 /*@@
   @routine    CCTK_ReduceArray
   @date       Aug 19 1999
   @author     Thomas Radke
   @desc
               Generic routine for doing a reduction operation on a set of
               arrays.
   @enddesc
   @var     GH
   @vdesc   pointer to the grid hierarchy
   @vtype   const cGH *
   @vio     in
   @endvar
   @var     proc
   @vdesc   processor that receives the result of the reduction operation
            (a negative value means that all processors get the result)
   @vtype   int
   @vio     in
   @endvar
   @var     operation_handle
   @vdesc   the handle specifying the reduction operator
   @vtype   int
   @vio     in
   @endvar
   @var     num_out_vals
   @vdesc   number of elements in the reduction output
   @vtype   int
   @vio     in
   @endvar
   @var     type_out_vals
   @vdesc   datatype of the reduction output
   @vtype   int
   @vio     in
   @endvar
   @var     out_vals
   @vdesc   pointer to buffer holding the output values
   @vtype   void *
   @vio     in
   @endvar
   @var     num_dims
   @vdesc   number of dimensions of input arrays
   @vtype   int
   @vio     in
   @endvar
   @var     num_in_fields
   @vdesc   number of input fields passed in the variable argument list
   @vtype   int
   @vio     in
   @endvar
   @var     type_in_arrays
   @vdesc   datatype of the input arrays
   @vtype   int
   @vio     in
   @endvar
   @var     <...>
   @vdesc   list of dimensions of input arrays
   @vtype   int
   @vio     in
   @endvar
   @var     <...>
   @vdesc   list of pointers to input arrays
   @vtype   void *
   @vio     in
   @endvar
@@*/
int CCTK_ReduceArray(const cGH *GH,
                     int proc,
                     int operation_handle,
                     int num_out_vals,
                     int type_out_vals,
                     void *out_vals,
                     int num_dims,
                     int num_in_arrays,
                     int type_in_arrays,
                     ... )
{
  va_list indices;
  int i;
  int *dims;
  const void **in_arrays;
  t_reduction_array_op *data;


  /* Get the pointer to the reduction operator */
  if (operation_handle < 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceArray: Invalid handle passed to CCTK_ReduceArray");
    return (-1);
  }

  data = Util_GetHandledData(ReductionArrayOperators,operation_handle);
  if (! data)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceArray: Array reduction operation is not registered "
                 "and cannot be called");
    return (-1);
  }

  /* allocate memory for dims and input array pointers */
  dims = malloc (num_dims * sizeof (int));
  in_arrays = malloc (num_in_arrays * sizeof (void *));

  /* Fill in the arrays of dims and input array pointers
     from the variable argument list */

  va_start(indices, type_in_arrays);

  for (i = 0; i < num_dims; i++)
  {
    dims[i] = va_arg (indices, int);
  }
  for (i = 0; i < num_in_arrays; i++)
  {
    in_arrays[i] = va_arg (indices, void *);
  }

  va_end(indices);

  data->function (GH, proc, num_dims, dims,
                  num_in_arrays, in_arrays, type_in_arrays,
                  num_out_vals, out_vals, type_out_vals);

  free (in_arrays);
  free (dims);

  return (0);
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReduceArray)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const int *num_out_vals,
      const int *type_out_vals,
      void *out_vals,
      const int *num_dims,
      const int *num_in_arrays,
      const int *type_in_arrays,
      ... )
{

  va_list varargs;
  int i;
  int *dims;
  const void **in_arrays;
  t_reduction_array_op *data;


  /* initialize return code to indicate an error */
  *fortran_return = -1;

  /* Get the pointer to the reduction operator */
  if (*operation_handle < 0)
  {
    CCTK_Warn (3,__LINE__,__FILE__,"Cactus",
               "CCTK_ReduceArray: Invalid handle passed to CCTK_ReduceArray");
    return;
  }

  data = Util_GetHandledData (ReductionArrayOperators, *operation_handle);
  if (! data)
  {
    CCTK_Warn (3,__LINE__,__FILE__,"Cactus",
               "CCTK_ReduceArray: Array reduction operation is not registered "
               "and cannot be called");
    return;
  }

  /* allocate memory for dims and input array pointers */
  dims = malloc (*num_dims * sizeof (int));
  in_arrays = malloc (*num_in_arrays * sizeof (void *));

  /* Fill in the arrays of dims and input array pointers
     from the variable argument list */

  va_start (varargs, type_in_arrays);

  for (i = 0; i < *num_dims; i++)
  {
    dims[i] = *va_arg (varargs, int *);
  }
  for (i = 0; i < *num_in_arrays; i++)
  {
    in_arrays[i] = va_arg (varargs, void *);
  }

  va_end (varargs);

  *fortran_return = data->function (*GH, *proc, *num_dims, dims,
                                    *num_in_arrays, in_arrays, *type_in_arrays,
                                    *num_out_vals, out_vals, *type_out_vals);
  free (in_arrays);
  free (dims);
}


 /*@@
   @routine    CCTK_ReduceLocalScalar
   @date       Aug 19 1999
   @author     Thomas Radke
   @desc    Wrapper function for reduction of a single scalar
   @enddesc
   @var     GH
   @vdesc   pointer to the grid hierarchy
   @vtype   const cGH *
   @vio     in
   @endvar
   @var     proc
   @vdesc   processor that receives the result of the reduction operation
            (a negative value means that all processors get the result)
   @vtype   int
   @vio     in
   @endvar
   @var     operation_handle
   @vdesc   the handle specifying the reduction operator
   @vtype   int
   @vio     in
   @endvar
   @var     in_scalar
   @vdesc   pointer to input scalar
   @vtype   void *
   @vio     in
   @endvar
   @var     out_scalar
   @vdesc   pointer to output scalar
   @vtype   void *
   @vio     in
   @endvar
   @var     data_type
   @vdesc   datatype for both input and output scalar
   @vtype   int
   @vio     in
   @endvar
@@*/

/*** FIXME: OLD INTERFACE gerd ***/
int CCTK_ReduceLocalScalar (const cGH *GH, int proc, int operation_handle,
                            const void *in_scalar, void *out_scalar, int data_type)
{
  return (CCTK_ReduceArray (GH, proc, operation_handle,
                            1, data_type, out_scalar,
                            1, 1, data_type, 1, in_scalar));
}

/*** FIXME: OLD INTERFACE gerd ***/
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocalScalar)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_scalar,
      void *out_scalar,
      const int *data_type)
{
  *fortran_return = CCTK_ReduceArray (*GH, *proc, *operation_handle,
                                      1, *data_type, out_scalar,
                                      1, 1, *data_type, 1, in_scalar);
}


int CCTK_ReduceLocScalar (const cGH *GH, int proc, int operation_handle,
                            const void *in_scalar, void *out_scalar, int data_type)
{
  return (CCTK_ReduceArray (GH, proc, operation_handle,
                            1, data_type, out_scalar,
                            1, 1, data_type, 1, in_scalar));
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocScalar)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_scalar,
      void *out_scalar,
      const int *data_type)
{
  *fortran_return = CCTK_ReduceArray (*GH, *proc, *operation_handle,
                                      1, *data_type, out_scalar,
                                      1, 1, *data_type, 1, in_scalar);
}


 /*@@
   @routine    CCTK_ReduceLocArrayToArray1D
   @date       Thu Oct 14 12:10:01 1999
   @author     Gerd Lanfermann
   @desc
        Interface to the migthy CCTK_Reduce for
        reduction of local 1D arrays to local arrays
        (element by element).
   @enddesc
@@*/

/*** FIXME: OLD INTERFACE gerd ***/
int CCTK_ReduceLocalArray1D (const cGH *GH, int proc, int operation_handle,
                             const void *in_array1d, void *out_array1d,
                             int num_in_array1d,
                             int data_type)
{
  return (CCTK_ReduceArray (GH, proc, operation_handle,
                            num_in_array1d, data_type, out_array1d,
                            1, 1, data_type, num_in_array1d, in_array1d));
}

int CCTK_ReduceLocArrayToArray1D(const cGH *GH, int proc, int operation_handle,
                                 const void *in_array1d, void *out_array1d,
                                 int num_in_array1d,
                                 int data_type)
{
  return (CCTK_ReduceArray (GH, proc, operation_handle,
                            num_in_array1d, data_type, out_array1d,
                            1, 1, data_type, num_in_array1d, in_array1d));
}

/*** FIXME: OLD INTERFACE gerd ***/
void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocalArray1D)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_array1d,
      void *out_array1d,
      const int *num_in_array1d,
      const int *data_type)
{
  *fortran_return = CCTK_ReduceArray (*GH, *proc, *operation_handle,
                                      *num_in_array1d, *data_type, out_array1d,
                                      1, 1, *data_type, *num_in_array1d,
                                      in_array1d);
}


/*@@
   @routine    CCTK_ReduceLocArrayToArray1D
   @date       Sat Nov 27 22:52:10 1999
   @author     Gerd Lanfermann
   @desc
       Interface for the reduction of local 1d arrays
       to the mighty CCCTK_reduce interface.
   @enddesc
@@*/


void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocArrayToArray1D)
     (int *fortran_return,
      const cGH **GH,
      const int *proc,
      const int *operation_handle,
      const void *in_array1d,
      void *out_array1d,
      const int *num_in_array1d,
      const int *data_type)
{
  *fortran_return = CCTK_ReduceArray (*GH, *proc, *operation_handle,
                                      *num_in_array1d, *data_type, out_array1d,
                                      1, 1, *data_type, *num_in_array1d,
                                      in_array1d);
}

/*@@
   @routine    CCTK_ReduceLocArrayToArray2D
   @date       Sat Nov 27 22:52:10 1999
   @author     Gerd Lanfermann
   @desc
       Interface for the reduction of local 2d arrays
       to the mighty CCCTK_reduce interface.
   @enddesc
@@*/



int CCTK_ReduceLocArrayToArray2D(const cGH *GH, int proc, int operation_handle,
                                 const void *in_array2d, void *out_array2d,
                                 int xsize, int ysize,
                                 int data_type)
{
  int  lin_size= xsize*ysize;
  return (CCTK_ReduceArray (GH, proc, operation_handle,
                            lin_size,
                            data_type, out_array2d,
                            2, 1, data_type,
                            xsize,ysize, in_array2d));
}

void CCTK_FCALL  CCTK_FNAME(CCTK_ReduceLocArrayToArray2D)
     (int  *fortran_return, const cGH **GH,
      const int  *proc,
      const int  *operation_handle,
      const void *in_array2d,
      void *out_array2d,
      const int  *xsize, const int *ysize,
      const int  *data_type)
{
  int lin_size = (*xsize)*(*ysize);
  *fortran_return =  CCTK_ReduceArray (*GH, *proc, *operation_handle,
                                      lin_size,
                                      *data_type, out_array2d,
                                      2, 1, *data_type,
                                      *xsize, *ysize,
                                      in_array2d);
}

/*@@
   @routine    CCTK_ReduceLocArrayToArray1D
   @date       Sat Nov 27 22:52:10 1999
   @author     Gerd Lanfermann
   @desc
       Interface for the reduction of local 3d arrays
       to 3d arrays.
   @enddesc
@@*/

int CCTK_ReduceLocArrayToArray3D(const cGH *GH, int proc, int operation_handle,
                                 const void *in_array3d, void *out_array3d,
                                 int xsize, int  ysize, int zsize,
                                 int data_type)
{

  int lin_size =  xsize*ysize*zsize;
  return (CCTK_ReduceArray (GH, proc, operation_handle,
                            lin_size,
                            data_type, out_array3d,
                            3, 1, data_type,
                            xsize,ysize,zsize,
                            in_array3d));
}

void CCTK_FCALL  CCTK_FNAME(CCTK_ReduceLocArrayToArray3D)
     (int  *fortran_return, const cGH **GH,
      const int  *proc,
      const int  *operation_handle,
      const void *in_array3d,
      void *out_array3d,
      const int  *xsize, const int *ysize, const int *zsize,
      const int  *data_type)
{
  int lin_size =  (*xsize)*(*ysize)*(*zsize);
  *fortran_return =  CCTK_ReduceArray (*GH, *proc, *operation_handle,
                                       lin_size,
                                       *data_type, out_array3d,
                                       3, 1, *data_type,
                                       *xsize,*ysize,*zsize,
                                       in_array3d);
}

 /* new local arrays api routines */
 /*@@
   @routine    CCTKi_RegisterLocalArrayReductionOperator
   @date       
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
   Registers "function" as a reduction operator called "name"
   @enddesc
   @var     function
   @vdesc   Routine containing reduction operator
   @vtype   (void (*))
   @vio
   @endvar
   @var     name
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTKi_RegisterLocalArrayReductionOperator(const char *thorn,
                                    cLocalArrayReduceOperator operator,
                                    const char *name)
{
  int handle;
  t_local_array_reduce_operator *reduce_operator;


  /* Check that the method hasn't already been registered */
  handle = Util_GetHandle(LocalArrayReductionOperators, name,
                          (void **) &reduce_operator);
  if(handle < 0)
  {
    reduce_operator = malloc (sizeof (t_reduce_operator));
    if (reduce_operator)
    {
      reduce_operator->implementation = CCTK_ThornImplementation(thorn);
      reduce_operator->name = name;
      reduce_operator->reduce_operator = operator;
      handle = Util_NewHandle(&LocalArrayReductionOperators, name, reduce_operator);

      /* Remember how many reduction operators there are */
      num_local_array_reductions++;
    }
  }
  else
  {
    /* Reduction operator with this name already exists. */
    CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
              "CCTK_RegisterLocalArrayReductionOperator: Reduction operator "
              "with this name already exists");
    handle = -1;
  }

  return handle;
}


 /*@@
   @routine    CCTK_LocalArrayReductionHandle
   @date       
   @author     Gabrielle Allen Yaakoub El Khamra
   @desc
   Returns the handle of a given local array reduction operator
   @enddesc
   @var     reduction
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTK_LocalArrayReductionHandle(const char *reduction)
{
  int handle;


  handle = Util_GetHandle(LocalArrayReductionOperators, reduction, NULL);
  if (handle < 0)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "CCTK_LocalArrayReductionHandle: No handle: '%d' found for reduction operator "
               "'%s'", handle, reduction);
  }

  return handle;
}

void CCTK_FCALL CCTK_FNAME(CCTK_LocalArrayReductionHandle)(int *handle, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(reduction)
  *handle = CCTK_LocalArrayReductionHandle(reduction);
  free(reduction);
}


 /*@@
   @routine    CCTK_ReduceLocalArrays
   @date       
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
               Generic routine for doing a reduction operation on a set of
               Cactus variables.
   @enddesc
   @var     output_values
   @vdesc   array of output values
   @vtype   void * const []
   @vio     in
   @endvar
   @var     N_dims
   @vdesc   dimension of the input arrays
   @vtype   int
   @vio     in
   @endvar
   @var     local_reduce_handle
   @vdesc   handle that identifies the reduction operator
   @vtype   int
   @vio     in
   @endvar
   @var     param_table_handle
   @vdesc   handle for a table containing additional parameters
   @vtype   int
   @vio     in
   @endvar
   @var     N_input_arrays
   @vdesc   number of input arrays
   @vtype   int
   @vio     in
   @endvar
   @var     input_array_sizes
   @vdesc   sizes of the input arrays
   @vtype   const CCTK_INT []
   @vio     in
   @endvar
   @var     input_array_type_codes
   @vdesc   types of the input arrays
   @vtype   const CCTK_INT []
   @vio     in
   @endvar
   @var     input_arrays
   @vdesc   input arrays
   @vtype   const void *const []
   @vio     in
   @endvar
   @var     M_output_values
   @vdesc   number of output values
   @vtype   int
   @vio     in
   @endvar
   @var     output_value_type_codes
   @vdesc   types of the output values
   @vtype   const CCTK_INT []
   @vio     in
   @endvar
   @var     output_values
   @vdesc   pointers to the output values
   @vtype   void *const []
   @vio     out
   @endvar
   @returntype int
   @returndesc
            negative for errors
   @endreturndesc
@@*/
int CCTK_ReduceLocalArrays(int N_dims,
                           int local_reduce_handle,
                           int param_table_handle,
                           int N_input_arrays,
                           const CCTK_INT input_array_sizes[],
                           const CCTK_INT input_array_type_codes[],
                           const void *const input_arrays[],
                           int M_output_values,
                           const CCTK_INT output_value_type_codes[],
                           void *const output_values[])
{
  int retval;
  t_local_array_reduce_operator *operator;

  /* Get the pointer to the reduction operator */
  if (local_reduce_handle < 0)
  {
    CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
              "CCTK_LocalArraysReduce: Invalid handle passed to CCTK_ReduceLocalArrays");
    retval = -1;
  }
  else
  {
    operator = Util_GetHandledData(LocalArrayReductionOperators,
                                   local_reduce_handle);

    if (!operator)
    {
      CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
                "CCTK_ReduceLocalArrays: Reduction operation is not registered"
                "and cannot be called");
      retval = -1;
    }
    else
    {
      retval = operator->reduce_operator (N_dims, local_reduce_handle, 
                          param_table_handle, N_input_arrays,
                          input_array_sizes, input_array_type_codes,
                          input_arrays, M_output_values,
                          output_value_type_codes, output_values);
    }
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReduceLocalArrays)
     (int *fortranreturn,int * N_dims, int *operation_handle, 
      int *param_table_handle, int * N_input_arrays,
      const CCTK_INT input_array_dims[], 
      const CCTK_INT input_array_type_codes[],
      const void *const input_arrays[],
      int * M_output_numbers, const CCTK_INT output_number_type_codes[],
      void *const output_numbers[])
{
  int retval;
  t_local_array_reduce_operator *operator;

  /* initialize return code to indicate an error */
  *fortranreturn = -1;

  if (*operation_handle < 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceLocalArrays: Invalid handle passed to CCTK_ReduceLocalArrays");
    retval = -1;
  }
  else
  {
    /* Get the pointer to the reduction operator */
    operator = Util_GetHandledData(LocalArrayReductionOperators,*operation_handle);

    if (!operator)
    {
      CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
                "CCTK_ReduceLocalArrays: Reduction operation is not registered"
                " and cannot be called");
      retval = -1;
    }
    else
    {
      retval = operator->reduce_operator (*N_dims, *operation_handle, 
                          *param_table_handle, *N_input_arrays,
                          input_array_dims, input_array_type_codes,
                          input_arrays, *M_output_numbers,
                          output_number_type_codes, output_numbers);
    }
  }
  *fortranreturn = retval;
}

 /*@@
   @routine    CCTK_NumReduceOperators
   @date       Mon Oct 22 2001
   @author     Gabrielle Allen
   @desc
               The number of reduction operators registered
   @enddesc
   @returntype int
   @returndesc
               number of reduction operators
   @endreturndesc
@@*/

int CCTK_NumReduceOperators()
{
  return num_reductions;
}

 /*@@
   @routine    CCTK_NumReduceOperators
   @date       
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
               The number of reduction operators registered
   @enddesc
   @returntype int
   @returndesc
               number of reduction operators
   @endreturndesc
@@*/

int CCTK_NumLocalArrayReduceOperators(void)
{
  return num_local_array_reductions;
}

void CCTK_FCALL CCTK_FNAME(CCTK_NumLocalArrayReduceOperators)(int *fortranreturn);
void CCTK_FCALL CCTK_FNAME(CCTK_NumLocalArrayReduceOperators)(int *fortranreturn)
{
  *fortranreturn = CCTK_NumLocalArrayReduceOperators();
}

 /*@@
   @routine    CCTK_ReduceOperatorImplementation
   @date       Mon Oct 22 2001
   @author     Gabrielle Allen
   @desc
   Provide the implementation which provides an reduction operator
   @enddesc
   @returntype int
   @returndesc
               Implementation which supplied the interpolation operator
   @endreturndesc
@@*/
const char *CCTK_ReduceOperatorImplementation(int handle)
{
  t_reduce_operator *operator;


  operator = Util_GetHandledData (ReductionOperators, handle);

  return (operator ? operator->implementation : NULL);
}

 /*@@
   @routine    CCTK_LocalArrayReduceOperatorImplementation
   @date       
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
   Provide the implementation which provides an local array reduction operator
   @enddesc
   @returntype int
   @returndesc
               Implementation which supplied the reduction operator
   @endreturndesc
@@*/
const char *CCTK_LocalArrayReduceOperatorImplementation(int handle)
{
  t_local_array_reduce_operator *operator;


  operator = Util_GetHandledData (LocalArrayReductionOperators, handle);

  return (operator ? operator->implementation : NULL);
}



 /*@@
   @routine    CCTK_ReduceOperator
   @date       December 27 2001
   @author     Gabrielle Allen
   @desc
               Returns the name of a reduction operator
   @enddesc
   @var        handle
   @vdesc      Handle for reduction operator
   @vtype      int
   @vio        in
   @endvar

   @returntype const char *
   @returndesc
   The name of the reduction operator, or NULL if the handle
   is invalid
   @endreturndesc
@@*/
const char *CCTK_ReduceOperator (int handle)
{
  const char *name=NULL;
  t_reduce_operator *operator;

  if (handle < 0)
  {
    CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                "CCTK_ReduceOperator: Handle %d invalid", handle);
  }
  else
  {
    operator = Util_GetHandledData (ReductionOperators, handle);
    if (operator)
    {
      name = operator->name;
    }
    else
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTK_ReduceOperator: Handle %d invalid", handle);
    }
  }

  return name;
}
 /*@@
   @routine    CCTK_LocalArrayReduceOperator
   @date       
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
               Returns the name of a reduction operator
   @enddesc
   @var        handle
   @vdesc      Handle for reduction operator
   @vtype      int
   @vio        in
   @endvar

   @returntype const char *
   @returndesc
   The name of the reduction operator, or NULL if the handle
   is invalid
   @endreturndesc
@@*/
const char *CCTK_LocalArrayReduceOperator (int handle)
{
  const char *name=NULL;
  t_local_array_reduce_operator *operator;

  if (handle < 0)
  {
    CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                "CCTK_LocalArrayReduceOperator: Handle %d invalid", handle);
  }
  else
  {
    operator = Util_GetHandledData (LocalArrayReductionOperators, handle);
    if (operator)
    {
      name = operator->name;
    }
    else
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTK_LocalArrayReduceOperator: Handle %d invalid", handle);
    }
  }

  return name;
}

/*@@
   @routine    CCTKi_RegisterGridArrayReductionOperator
   @date       Mon Aug 30 11:27:43 2004
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
   Registers "function" as a reduction operator called "name"
   @enddesc
   @var     function
   @vdesc   Routine containing reduction operator
   @vtype   (void (*))
   @vio
   @endvar
   @var     name
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTKi_RegisterGridArrayReductionOperator(const char *thorn, 
                                             cGridArrayReduceOperator  operator)
{
  int handle;
  t_grid_array_reduce_operator *reduce_operator;

  /* Check that there is no other registered GA reduction */
  if(num_GA_reductions == 0)
  {
    reduce_operator = (t_grid_array_reduce_operator *)malloc (sizeof (t_grid_array_reduce_operator));
    reduce_operator->implementation = CCTK_ThornImplementation(thorn);
    reduce_operator->reduce_operator = operator;
    GA_reduc = operator;
    handle = Util_NewHandle(&GridArrayReductionOperators, global, reduce_operator);
    /* Remember how many reduction operators there are */
    num_GA_reductions++;
  }
  else
  {
    /* Reduction operator with this name already exists. */
    CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
              "CCTK_RegisterGridArrayReductionOperator: Reduction operator "
              "already exists");
    handle = -1;
  }
  return handle;
}

 /*@@
   @routine    CCTK_ReduceGridArrays
   @date       Mon Aug 30 11:27:43 2004
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
               Generic routine for doing a reduction operation on a set of
               Cactus variables.
   @enddesc
   @var     GH
   @vdesc   pointer to the grid hierarchy
   @vtype   cGH *
   @vio     in
   @endvar
   @var     dest_proc
   @vdesc   the number of the processor to which we want to reduce (-1) for all-reduce
   @vtype   int
   @vio     in
   @endvar
   @var     local_reduce_handle
   @vdesc   the handle specifying the reduction operator
   @vtype   int
   @vio     in
   @endvar
   @var     param_table_handle
   @vdesc   the parameter table handle
   @vtype   int
   @vio     in
   @endvar
   @var     N_input_arrays
   @vdesc   number of elements in the reduction input
   @vtype   int
   @vio     in
   @endvar
   @var     input_array_variable_indices
   @vdesc   input arrays
   @vtype   const CCTK_INT []
   @vio     in
   @endvar
   @var     M_output_values
   @vdesc   number of output values
   @vtype   int
   @vio     in
   @endvar
   @var     output_value_type_codes
   @vdesc   array containing output value type codes
   @vtype   const CCTK_IN []
   @vio     in
   @endvar
   @var     output_values
   @vdesc   array of output values
   @vtype   void * const []
   @vio     in
   @endvar
   @returntype int
   @returndesc
            negative for errors
   @endreturndesc
@@*/
int CCTK_ReduceGridArrays(const cGH *GH,
                          int dest_proc,
                          int local_reduce_handle,
                          int param_table_handle,
                          int N_input_arrays,
                          const CCTK_INT input_array_variable_indices[],
                          int M_output_values,
                          const CCTK_INT output_value_type_codes[],
                          void* const output_values[])
{
  int retval;

  /* Get the pointer to the reduction operator */
  if (num_GA_reductions == 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceGridArrays: no grid array reduction registered");
    retval = -1;
  }
  else
  {
    retval = GA_reduc (GH,
                       dest_proc,
                       local_reduce_handle, param_table_handle,
                       N_input_arrays, input_array_variable_indices,
                       M_output_values, output_value_type_codes,
                       output_values);
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReduceGridArrays)
     (int *fortranreturn,
      const cGH **GH,
      int *dest_proc,
      int *local_reduce_handle,
      int *param_table_handle,
      int *N_input_arrays,
      const CCTK_INT input_array_variable_indices[],
      int *M_output_values,
      const CCTK_INT output_value_type_codes[],
      void* const output_values[])
{
  int retval;

  /* Get the pointer to the reduction operator */
  if (num_GA_reductions == 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceGridArrays: no grid array reduction registered");
    retval = -1;
  }
  else
  {
    retval = GA_reduc (*GH,
                       *dest_proc,
                       *local_reduce_handle, *param_table_handle,
                       *N_input_arrays, input_array_variable_indices,
                       *M_output_values, output_value_type_codes,
                       output_values);
  }
  *fortranreturn = retval;
}

 /*@@
   @routine    CCTK_NumGridArrayReductionOperators
   @date       Mon Aug 30 11:27:43 2004
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
               The number of reduction operators registered
   @enddesc
   @returntype int
   @returndesc
               number of reduction operators
   @endreturndesc
@@*/

int CCTK_NumGridArrayReductionOperators(void)
{
  return num_GA_reductions;
}

 /*@@
   @routine    CCTK_GAReductionOperator
   @date       Mon Aug 30 11:27:43 2004
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
               Returns the name of the reduction operator
   @enddesc
   @var        handle
   @vdesc      Handle for reduction operator
   @vtype      int
   @vio        in
   @endvar

   @returntype const char *
   @returndesc
   The name of the reduction operator, or NULL if the handle
   is invalid
   @endreturndesc
@@*/
const char *CCTK_GridArrayReductionOperator (void)
{
  const char *thorn=NULL;
  t_grid_array_reduce_operator *operator;

  if (num_GA_reductions == 0)
  {
    CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                "CCTK_GAReductionOperator: no GA reduction operator");
  }
  else
  {
    Util_GetHandle (GridArrayReductionOperators, global, (void **)&operator);
    if (operator)
    {
      thorn = operator->implementation;
    }
    else
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTK_GAReductionOperator: no GA reduction operator");
    }
  }

  return thorn;
}

/*@@
   @routine    CCTKi_RegisterReduceArraysGloballyOperator
   @date       Mon Aug 30 11:27:43 2004
   @author     Gabrielle Allen, Yaakoub El Khamra
   @desc
   Registers "function" as a reduction operator called "name"
   @enddesc
   @var     function
   @vdesc   Routine containing reduction operator
   @vtype   (void (*))
   @vio
   @endvar
   @var     name
   @vdesc   String containing name of reduction operator
   @vtype   const char *
   @vio     in
   @endvar
@@*/
int CCTKi_RegisterReduceArraysGloballyOperator(const char *thorn, 
                                             cReduceArraysGloballyOperator  operator)
{
  int handle;
  t_reduce_arrays_globally_operator *reduce_operator;

  /* Check that there is no other registered GA reduction */
  if(num_reductions_arrays_globally == 0)
  {
    reduce_operator = (t_reduce_arrays_globally_operator *)malloc (sizeof (t_reduce_arrays_globally_operator));
    reduce_operator->implementation = CCTK_ThornImplementation(thorn);
    reduce_operator->reduce_operator = operator;
    ArraysGlobally_reduc = operator;
    handle = Util_NewHandle(&ReductionArraysGloballyOperators, reduction_arrays_globally_global, reduce_operator);
    /* Remember how many reduction operators there are */
    num_reductions_arrays_globally++;
  }
  else
  {
    /* Reduction operator with this name already exists. */
    CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
              "CCTKi_RegisterReduceArraysGloballyOperator: Reduction operator "
              "already exists");
    handle = -1;
  }
  return handle;
}

 /*@@
   @routine    CCTK_ReduceArraysGlobally
   @date       
   @author     Yaakoub El Khamra
   @desc
               Generic routine for doing a reduction operation on a set of
               arrays in a global manner
   @enddesc
   @var     GH
   @vdesc   pointer to the grid hierarchy
   @vtype   cGH *
   @vio     in
   @endvar
   @var     dest_proc
   @vdesc   the number of the processor to which we want to reduce (-1) for all-reduce
   @vtype   int
   @vio     in
   @endvar
   @var     local_reduce_handle
   @vdesc   the handle specifying the reduction operator
   @vtype   int
   @vio     in
   @endvar
   @var     param_table_handle
   @vdesc   the parameter table handle
   @vtype   int
   @vio     in
   @endvar
   @var     N_input_arrays
   @vdesc   number of elements in the reduction input
   @vtype   int
   @vio     in
   @endvar
   @var     input_arrays
   @vdesc   input arrays
   @vtype   const void *const []
   @vio     in
   @endvar
   @var     input_dims
   @vdesc   number of dimensions
   @vtype   int
   @vio     in
   @endvar
   @var     input_array_sizes
   @vdesc   sizes of the input arrays
   @vtype   const CCTK_INT []
   @vio     in
   @endvar
   @var     input_array_type_codes
   @vdesc   types of the input arrays
   @vtype   const CCTK_INT []
   @vio     in
   @endvar
   @var     M_output_values
   @vdesc   number of output values
   @vtype   int
   @vio     in
   @endvar
   @var     output_value_type_codes
   @vdesc   array containing output value type codes
   @vtype   const CCTK_IN []
   @vio     in
   @endvar
   @var     output_values
   @vdesc   array of output values
   @vtype   void * const []
   @vio     in
   @endvar
   @returntype int
   @returndesc
            negative for errors
   @endreturndesc
@@*/
int CCTK_ReduceArraysGlobally(const cGH *GH,
                          int dest_proc,
                          int local_reduce_handle,
                          int param_table_handle,
                          int N_input_arrays,
                          const void * const input_arrays[],
                          int input_dims,
                          const CCTK_INT input_array_dims[],
                          const CCTK_INT input_array_type_codes[],
                          int M_output_values,
                          const CCTK_INT output_value_type_codes[],
                          void* const output_values[])
{
  int retval;

  /* Get the pointer to the reduction operator */
  if (num_reductions_arrays_globally == 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceArraysGlobally: no grid array reduction registered");
    retval = -1;
  }
  else
  {
    retval = ArraysGlobally_reduc (GH,
                       dest_proc,
                       local_reduce_handle, param_table_handle,
                       N_input_arrays, input_arrays, input_dims,
                       input_array_dims, input_array_type_codes,
                       M_output_values, output_value_type_codes,
                       output_values);
  }
  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_ReduceArraysGlobally)
     (int *fortranreturn,
      const cGH **GH,
      int *dest_proc,
      int *local_reduce_handle,
      int *param_table_handle,
      int *N_input_arrays,
      const void * const input_arrays[],
      int *input_dims,
      const CCTK_INT input_array_dims[],
      const CCTK_INT input_array_type_codes[],
      int *M_output_values,
      const CCTK_INT output_value_type_codes[],
      void* const output_values[])
{
  int retval;

  /* Get the pointer to the reduction operator */
  if (num_reductions_arrays_globally == 0)
  {
    CCTK_Warn(3,__LINE__,__FILE__,"Cactus",
              "CCTK_ReduceArraysGlobally: no grid array reduction registered");
    retval = -1;
  }
  else
  {
    retval = ArraysGlobally_reduc (*GH,
                       *dest_proc,
                       *local_reduce_handle, *param_table_handle,
                       *N_input_arrays, input_arrays, *input_dims,
                       input_array_dims, input_array_type_codes,
                       *M_output_values, output_value_type_codes,
                       output_values);
  }
  *fortranreturn = retval;
}

 /*@@
   @routine    CCTK_NumReductionArraysGloballyOperators
   @date       
   @author     Yaakoub El Khamra
   @desc
               The number of pointwie reduction operators registered
   @enddesc
   @returntype int
   @returndesc
               number of reduction operators
   @endreturndesc
@@*/

int CCTK_NumReductionArraysGloballyOperators(void)
{
  return num_reductions_arrays_globally;
}

 /*@@
   @routine    CCTK_ReductionArraysGloballyOperator
   @date       
   @author     Yaakoub El Khamra
   @desc
               Returns the name of the reduction operator
   @enddesc
   @var        handle
   @vdesc      Handle for reduction operator
   @vtype      int
   @vio        in
   @endvar

   @returntype const char *
   @returndesc
   The name of the reduction operator, or NULL if the handle
   is invalid
   @endreturndesc
@@*/
const char *CCTK_ReductionArraysGloballyOperator (void)
{
  const char *thorn=NULL;
  t_reduce_arrays_globally_operator *operator;

  if (num_reductions_arrays_globally == 0)
  {
    CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                "CCTK_ReductionArraysGloballyOperator: no reduction of arrays globally operator");
  }
  else
  {
    Util_GetHandle (ReductionArraysGloballyOperators, global, (void **)&operator);
    if (operator)
    {
      thorn = operator->implementation;
    }
    else
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTK_ReductionArraysGloballyOperator: no reduction of arrays globally operator");
    }
  }

  return thorn;
}
