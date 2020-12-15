#include "cctk.h"

module cctk_Interp
  implicit none

  interface

     subroutine CCTK_InterpHandle (handle, name)
       implicit none
       integer      handle
       character(*) name
     end subroutine CCTK_InterpHandle

#if 0
C    no Fortran 90 interface for following routines
C    their names exceed the maximum allowed length of 31 characters
     subroutine CCTK_InterpRegisterOpLocalUniform (operator_ptr, operator_name, thorn)
       implicit none
       CCTK_FPOINTER operator_ptr
       character(*)  operator_name
       character(*)  thorn
     end subroutine CCTK_InterpRegisterOpLocalUniform

     subroutine CCTK_InterpOperatorImplementation (nchars, imp, handle)
       implicit none
       integer      nchars
       character(*) imp
       integer      handle
     end subroutine CCTK_InterpOperatorImplementation
#endif

     subroutine CCTK_InterpOperator (nchars, operator, handle)
       implicit none
       integer      nchars
       character(*) operator
       integer      handle
     end subroutine CCTK_InterpOperator

     subroutine CCTK_NumInterpOperators (num)
       implicit none
       integer num
     end subroutine CCTK_NumInterpOperators

     ! CCTK_InterpGV is declared below
     ! CCTK_InterpLocal is declared below

     subroutine CCTK_InterpLocalUniform &
          (ierr,                        &
           N_dims,                      &
           operator_handle,             &
           param_table_handle,          &
           coord_origin,                &
           coord_delta,                 &
           N_interp_points,             &
             interp_coords_type_code,   &
             interp_coords,             &
           N_input_arrays,              &
             input_array_dims,          &
             input_array_type_codes,    &
             input_arrays,              &
           N_output_arrays,             &
             output_array_type_codes,   &
             output_arrays)
       implicit none
       integer               ierr
       integer               N_dims
       integer               operator_handle
       integer               param_table_handle
       CCTK_REAL             coord_origin(N_dims)
       CCTK_REAL             coord_delta(N_dims)
       integer               N_interp_points
       integer                 interp_coords_type_code
       CCTK_POINTER_TO_CONST   interp_coords(N_dims)
       integer               N_input_arrays
       CCTK_INT                input_array_dims(N_dims)
       CCTK_INT                input_array_type_codes(N_input_arrays)
       CCTK_POINTER_TO_CONST   input_arrays(N_input_arrays)
       integer               N_output_arrays
       CCTK_INT                output_array_type_codes(N_output_arrays)
       CCTK_POINTER            output_arrays(N_output_arrays)
     end subroutine CCTK_InterpLocalUniform

  end interface

end module cctk_Interp
