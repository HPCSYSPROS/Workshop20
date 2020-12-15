#include "cctk.h"

module cctk_Reduction
  implicit none
  
  interface

     subroutine CCTK_ReductionHandle (ierr, reduction)
       implicit none
       integer      ierr
       character(*) reduction
     end subroutine CCTK_ReductionHandle

     subroutine CCTK_ReductionArrayHandle (ierr, reduction)
       implicit none
       integer      ierr
       character(*) reduction
     end subroutine CCTK_ReductionArrayHandle
     
     ! Name is too long
!!$     subroutine CCTK_ReduceOperatorImplementation &
!!$          (implementation, implementation_length, handle)
!!$       implicit none
!!$       character(*) implementation
!!$       integer      implementation_length
!!$       integer      handle
!!$     end subroutine CCTK_ReduceOperatorImplementation
     
     subroutine CCTK_ReduceOperator (operator, operator_length, handle)
       implicit none
       character(*) operator
       integer      operator_length
       integer      handle
     end subroutine CCTK_ReduceOperator
     
     subroutine CCTK_NumReduceOperators (num_operators)
       implicit none
       integer num_operators
     end subroutine CCTK_NumReduceOperators
     
     subroutine CCTK_ReduceLocalArrays (ierr, &
          N_dims, &
          operator_handle, &
          param_table_handle, &
          N_input_arrays, &
          input_array_dims, &
          input_array_type_codes, &
          input_arrays, &
          M_output_numbers, &
          output_number_type_codes, &
          output_numbers)
       implicit none
       integer               ierr
       integer               N_dims
       integer               operator_handle
       integer               param_table_handle
       integer               N_input_arrays
       CCTK_INT              input_array_dims(N_input_arrays)
       CCTK_INT              input_array_type_codes(N_input_arrays)
       CCTK_POINTER_TO_CONST input_arrays(N_input_arrays)
       integer               M_output_numbers
       CCTK_INT              output_number_type_codes(M_output_numbers)
       CCTK_POINTER          output_numbers(M_output_numbers)
     end subroutine CCTK_ReduceLocalArrays
     
     subroutine CCTK_LocalArrayReductionHandle (ierr, reduction)
       implicit none
       integer      ierr
       character(*) reduction
     end subroutine CCTK_LocalArrayReductionHandle
     
     ! Name is too long
!!$     subroutine CCTK_LocalArrayReduceOperatorImplementation &
!!$          (implementation, implementation_length, handle)
!!$       implicit none
!!$       character(*) implementation
!!$       integer      implementation_length
!!$       integer      handle
!!$     end subroutine CCTK_LocalArrayReduceOperatorImplementation
     
     subroutine CCTK_LocalArrayReduceOperator &
          (operator, operator_length, handle)
       implicit none
       character(*) operator
       integer      operator_length
       integer      handle
     end subroutine CCTK_LocalArrayReduceOperator
     
     ! Name is too long
!!$     subroutine CCTK_NumLocalArrayReduceOperators (num_operators)
!!$       implicit none
!!$       integer num_operators
!!$     end subroutine CCTK_NumLocalArrayReduceOperators
     
     subroutine CCTK_ReduceGridArrays (ierr, cctkGH, &
          dest_proc, &
          local_reduce_handle, &
          param_table_handle, &
          N_input_arrays, &
          input_array_variable_indices, &
          M_output_values, &
          output_value_type_codes, &
          output_values)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST cctkGH
       integer               dest_proc
       integer               local_reduce_handle
       integer               param_table_handle
       integer               N_input_arrays
       CCTK_INT              input_array_variable_indices(N_input_arrays)
       integer               M_output_values
       CCTK_INT              output_value_type_codes(M_output_values)
       CCTK_POINTER          output_values(M_output_values)
     end subroutine CCTK_ReduceGridArrays
     
     subroutine CCTK_GridArrayReductionOperator (operator, operator_length)
       implicit none
       character(*) operator
       integer      operator_length
     end subroutine CCTK_GridArrayReductionOperator
     
     ! Name is too long
!!$     subroutine CCTK_NumGridArrayReductionOperators (num_operators)
!!$       implicit none
!!$       integer num_operators
!!$     end subroutine CCTK_NumGridArrayReductionOperators
     
!!$     subroutine CCTK_ReduceLocalScalar (ierr, &
!!$          cctkGH,  proc, operation_handle, &
!!$          inScalar, outScalar, dataType)
!!$       implicit none
!!$       integer               ierr
!!$       CCTK_POINTER_TO_CONST cctkGH
!!$       integer               proc
!!$       integer               operation_handle
!!$       CCTK_POINTER_TO_CONST inScalar
!!$       CCTK_POINTER          outScalar
!!$       integer               dataType
!!$     end subroutine CCTK_ReduceLocalScalar
     
!!$     subroutine CCTK_ReduceLocalArray1D (ierr, &
!!$          cctkGH,  proc, operation_handle, &
!!$          in_array1d, out_array1d, num_in_array1d, data_type)
!!$       implicit none
!!$       integer               ierr
!!$       CCTK_POINTER_TO_CONST cctkGH
!!$       integer               proc
!!$       integer               operation_handle
!!$       CCTK_POINTER_TO_CONST in_array1d
!!$       CCTK_POINTER          out_array1d
!!$       integer               num_in_array1d
!!$       integer               data_type
!!$     end subroutine CCTK_ReduceLocalArray1D
     
!!$     subroutine CCTK_ReduceLocScalar (ierr, &
!!$          cctkGH, proc, operation_handle, &
!!$          in_scalar, out_scalar, data_type)
!!$       implicit none
!!$       integer               ierr
!!$       CCTK_POINTER_TO_CONST cctkGH
!!$       integer               proc
!!$       integer               operation_handle
!!$       CCTK_POINTER_TO_CONST in_scalar
!!$       CCTK_POINTER          out_scalar
!!$       integer               data_type
!!$     end subroutine CCTK_ReduceLocScalar
     
!!$     subroutine CCTK_ReduceLocArrayToArray1D (ierr, &
!!$          cctkGH, proc, operation_handle, &
!!$          in_array1d, out_array1d, num_in_array1d, data_type)
!!$       implicit none
!!$       integer               ierr
!!$       CCTK_POINTER_TO_CONST cctkGH
!!$       integer               proc
!!$       integer               operation_handle
!!$       CCTK_POINTER_TO_CONST in_array1d
!!$       CCTK_POINTER          out_array1d
!!$       integer               num_in_array1d
!!$       integer               data_type
!!$     end subroutine CCTK_ReduceLocArrayToArray1D
     
!!$     subroutine CCTK_ReduceLocArrayToArray2D (ierr, &
!!$          cctkGH, proc, operation_handle, &
!!$          in_array2d, out_array2d, xsize, ysize, data_type)
!!$       implicit none
!!$       integer               ierr
!!$       CCTK_POINTER_TO_CONST cctkGH
!!$       integer               proc
!!$       integer               operation_handle
!!$       CCTK_POINTER_TO_CONST in_array2d
!!$       CCTK_POINTER          out_array2d
!!$       integer               xsize
!!$       integer               ysize
!!$       integer               data_type
!!$     end subroutine CCTK_ReduceLocArrayToArray2D
     
!!$     subroutine CCTK_ReduceLocArrayToArray3D (ierr, &
!!$          cctkGH, proc, operation_handle, &
!!$          in_array3d, out_array3d, xsize, ysize, zsize, data_type)
!!$       implicit none
!!$       integer               ierr
!!$       CCTK_POINTER_TO_CONST cctkGH
!!$       integer               proc
!!$       integer               operation_handle
!!$       CCTK_POINTER_TO_CONST in_array3d
!!$       CCTK_POINTER          out_array3d
!!$       integer               xsize
!!$       integer               ysize
!!$       integer               zsize
!!$       integer               data_type
!!$     end subroutine CCTK_ReduceLocArrayToArray3D
     
  end interface
  
  external CCTK_Reduce
  external CCTK_Reduce1
  external CCTK_ReduceLocalScalar
  external CCTK_ReduceLocalArray1D
  external CCTK_ReduceLocScalar
  external CCTK_ReduceLocArrayToArray1D
  external CCTK_ReduceLocArrayToArray2D
  external CCTK_ReduceLocArrayToArray3D
  external CCTK_ReduceArray
  
end module cctk_Reduction
