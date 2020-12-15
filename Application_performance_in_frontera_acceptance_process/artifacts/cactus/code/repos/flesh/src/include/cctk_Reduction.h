 /*@@
   @header    cctk_Reduction.h
   @date
   @author    Gabrielle Allen
   @desc
              Header file for using reduction operators
   @enddesc
   @version   $Header$
 @@*/

#ifndef _CCTK_REDUCTION_H_
#define _CCTK_REDUCTION_H_ 1

#define REDUCTION_OPERATOR_REGISTER_ARGLIST  \
          const cGH *arg_GH, \
          int arg_proc, \
          int arg_num_outvals, \
          int arg_outtype, \
          void *arg_outvals, \
          int arg_num_invars, \
          const int arg_varlist []


#define REDUCTION_ARRAY_OPERATOR_REGISTER_ARGLIST \
  const cGH *arg_GH, \
  int arg_proc, \
  int arg_nDims, \
  const int arg_dims [], \
  int arg_nArrays, \
  const void *const arg_inArrays [], \
  int arg_inType, \
  int arg_nOutVals, \
  void *arg_outVals, \
  int arg_outType

#define REDUCTION_LOCAL_ARRAY_OPERATOR_REGISTER_ARGLIST \
  int N_dims, \
  int operator_handle, \
  int param_table_handle, \
  int N_input_arrays, \
  const CCTK_INT input_array_dims[], \
  const CCTK_INT input_array_type_codes[], \
  const void *const input_arrays[], \
  int M_output_numbers, \
  const CCTK_INT output_number_type_codes[], \
  void *const output_numbers[]

#define REDUCTION_GRID_ARRAY_OPERATOR_REGISTER_ARGLIST \
  const cGH *GH,  \
  int dest_proc,  \
  int local_reduce_handle, \
  int param_table_handle,  \
  int N_input_arrays,  \
  const CCTK_INT input_array_variable_indices[],  \
  int M_output_values,  \
  const CCTK_INT output_value_type_codes[],  \
  void* const output_values[]

#define REDUCTION_ARRAYS_GLOBALLY_OPERATOR_REGISTER_ARGLIST \
  const cGH *GH,  \
  int dest_proc,  \
  int local_reduce_handle, \
  int param_table_handle,  \
  int N_input_arrays,  \
  const void *const input_arrays[], \
  int N_dims,\
  const CCTK_INT input_array_dims[], \
  const CCTK_INT input_array_type_codes[], \
  int M_output_values,  \
  const CCTK_INT output_value_type_codes[],  \
  void* const output_values[]
    
#ifdef __cplusplus
extern "C"
{
#endif

/* prototype for reduction operator routine */
typedef int (*cReduceOperator) (const cGH *GH,
                                int arg_proc,
                                int arg_num_outvals,
                                int arg_outtype,
                                void *arg_outvals,
                                int arg_num_invars,
                                const int arg_varlist[]);

/* prototype for local array reduction operator routine */
typedef int (*cLocalArrayReduceOperator) (int N_dims, int operator_handle, 
                          int param_table_handle,   int N_input_arrays,
                          const CCTK_INT input_array_dims[], 
                          const CCTK_INT input_array_type_codes[],
                          const void *const input_arrays[],
                          int M_output_numbers,
                          const CCTK_INT output_number_type_codes[],
                          void *const output_numbers[]);

/* prototype for GA reduction operator routine */
typedef int (*cGridArrayReduceOperator) (const cGH *GH,
                                         int dest_proc,
                                         int local_reduce_handle,
                                         int param_table_handle,
                                         int N_input_arrays,
                                         const CCTK_INT input_array_variable_indices[],
                                         int M_output_values,
                                         const CCTK_INT output_value_type_codes[],
                                         void* const output_values[]);

/* prototype for global array reduction operator routine */
typedef int (*cReduceArraysGloballyOperator) (const cGH *GH,
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
                                         void* const output_values[]);

int CCTK_Reduce(const cGH *GH,
                int proc,
                int operation_handle,
                int num_out_vals,
                int type_out_vals,
                void *out_vals,
                int num_in_fields, ...);

int CCTK_ReductionHandle(const char *reduction);

#define CCTK_RegisterReductionOperator(a,b) \
        CCTKi_RegisterReductionOperator(CCTK_THORNSTRING,a,b)

int CCTKi_RegisterReductionOperator(const char *thorn,
                                    cReduceOperator operatorGV,
                                    const char *name);

int CCTK_ReductionArrayHandle(const char *reduction);

int CCTK_RegisterReductionArrayOperator(
         int (*function)(REDUCTION_ARRAY_OPERATOR_REGISTER_ARGLIST),
         const char *name);

const char *CCTK_ReduceOperatorImplementation(int handle);

const char *CCTK_ReduceOperator (int handle);

int CCTK_NumReduceOperators(void);

/* new local array reduction API */
int CCTK_ReduceLocalArrays(int N_dims,
                           int local_reduce_handle,
                           int param_table_handle,
                           int N_input_arrays,
                           const CCTK_INT input_array_sizes[],
                           const CCTK_INT input_array_type_codes[],
                           const void *const input_arrays[],
                           int M_output_values,
                           const CCTK_INT output_value_type_codes[],
                           void *const output_values[]);

int CCTK_LocalArrayReductionHandle(const char *reduction);

#define CCTK_RegisterLocalArrayReductionOperator(a,b) \
        CCTKi_RegisterLocalArrayReductionOperator(CCTK_THORNSTRING,a,b)

int CCTKi_RegisterLocalArrayReductionOperator(const char *thorn,
                           cLocalArrayReduceOperator operatorGV,
                                              const char *name);

int CCTK_RegisterReductionLocalArrayOperator(
         int (*function)(REDUCTION_LOCAL_ARRAY_OPERATOR_REGISTER_ARGLIST),
         const char *name);

const char *CCTK_LocalArrayReduceOperatorImplementation(int handle);

const char *CCTK_LocalArrayReduceOperator (int handle);
int CCTK_NumLocalArrayReduceOperators(void);

/* new GA reduction API */
int CCTK_ReduceGridArrays(const cGH *GH,
                          int dest_proc,
                          int local_reduce_handle,
                          int param_table_handle,
                          int N_input_arrays,
                          const CCTK_INT input_array_variable_indices[],
                          int M_output_values,
                          const CCTK_INT output_value_type_codes[],
                          void* const output_values[]);

#define CCTK_RegisterGridArrayReductionOperator(a) \
        CCTKi_RegisterGridArrayReductionOperator(CCTK_THORNSTRING,a)

int CCTKi_RegisterGridArrayReductionOperator(const char *thorn, cGridArrayReduceOperator 
        operatorGV);

const char *CCTK_GridArrayReductionOperator(void);
int CCTK_NumGridArrayReductionOperators(void);

/* new pointwise reduction API */
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
                          void* const output_values[]);

#define CCTK_RegisterReduceArraysGloballyOperator(a) \
        CCTKi_RegisterReduceArraysGloballyOperator(CCTK_THORNSTRING,a)

int CCTKi_RegisterReduceArraysGloballyOperator(const char *thorn, cReduceArraysGloballyOperator 
        operatorGV);

const char *CCTK_ReductionArraysGloballyOperator(void);
int CCTK_NumReductionArraysGloballyOperators(void);

/* FIXME: old interface - should go */
int CCTK_ReduceLocalScalar (const cGH *GH, int proc, int operation_handle,
                            const void *inScalar, void *outScalar, int dataType);


/* FIXME: old interface - should go */
int CCTK_ReduceLocalArray1D (const cGH *GH, int proc, int operation_handle,
                            const void *in_array1d, void *out_array1d,
                            int num_in_array1d, int data_type);

int CCTK_ReduceLocScalar(const cGH *GH, int proc, int operation_handle,
                         const void *in_scalar, void *out_scalar, int data_type);

int CCTK_ReduceLocArrayToArray1D(const cGH *GH, int proc, int operation_handle,
                                 const void *in_array1d, void *out_array1d,
                                 int num_in_array1d,
                                 int data_type);

int CCTK_ReduceLocArrayToArray2D(const cGH *GH, int proc, int operation_handle,
                                 const void *in_array2d, void *out_array2d,
                                 int xsize, int ysize,
                                 int data_type);

int CCTK_ReduceLocArrayToArray3D(const cGH *GH, int proc, int operation_handle,
                                 const void *in_array3d, void *out_array3d,
                                 int xsize, int  ysize, int zsize,
                                 int data_type);


int CCTK_ReduceArray(const cGH *GH,
                     int proc,
                     int operation_handle,
                     int num_out_vals,
                     int type_out_vals,
                     void *out_vals,
                     int num_dims,
                     int num_in_arrays,
                     int type_in_arrays,
                     ... );
                          


#ifdef __cplusplus
}
#endif

#endif /* _CCTK_REDUCTION_H_ */
