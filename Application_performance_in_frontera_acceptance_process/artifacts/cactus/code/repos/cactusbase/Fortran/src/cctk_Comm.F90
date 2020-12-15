#include "cctk.h"

module cctk_Comm
  implicit none

  interface

     subroutine CCTK_SyncGroup (ierr, GH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          group
     end subroutine CCTK_SyncGroup

     subroutine CCTK_Barrier (ierr, GH)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
     end subroutine CCTK_Barrier

! This routine has been made a function instead of a subroutine
! and is declared in cctk.h
!     subroutine CCTK_MyProc (ierr, GH)
!       implicit none
!       integer               ierr
!       CCTK_POINTER_TO_CONST GH
!     end subroutine CCTK_MyProc

! This routine has been made a function instead of a subroutine
! and is declared in cctk.h
!     subroutine CCTK_nProcs (ierr, GH)
!       implicit none
!       integer               ierr
!       CCTK_POINTER_TO_CONST GH
!     end subroutine CCTK_nProcs

     subroutine CCTK_ParallelInit (ierr, GH)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
     end subroutine CCTK_ParallelInit

     subroutine CCTK_Exit (ierr, GH, retval)
       implicit none
       integer      ierr
       CCTK_POINTER GH
       integer      retval
     end subroutine CCTK_Exit

     subroutine CCTK_Abort (ierr, GH, retval)
       implicit none
       integer      ierr
       CCTK_POINTER GH
       integer      retval
     end subroutine CCTK_Abort

     subroutine CCTK_SetupGH (GH, config, convergence_level)
       implicit none
       CCTK_POINTER GH
       CCTK_POINTER config
       integer      convergence_level
     end subroutine CCTK_SetupGH

! CCTK_ArrayGroupSizeB is not available from Fortran

     subroutine CCTK_QueryGroupStorageB (ierr, GH, group, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
       character(*)          groupname
     end subroutine CCTK_QueryGroupStorageB

     subroutine CCTK_GroupDynamicData (ierr, GH, group, data)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
       CCTK_POINTER          data
     end subroutine CCTK_GroupDynamicData

     subroutine CCTK_GroupStorageIncrease &
          (ierr, GH, n_groups, groups, timelevels, status)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               n_groups
       integer               groups(n_groups)
       integer               timelevels(n_groups)
       integer               status(n_groups)
     end subroutine CCTK_GroupStorageIncrease

     subroutine CCTK_GroupStorageDecrease &
          (ierr, GH, n_groups, groups, timelevels, status)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               n_groups
       integer               groups(n_groups)
       integer               timelevels(n_groups)
       integer               status(n_groups)
     end subroutine CCTK_GroupStorageDecrease

     subroutine CCTK_InterpGridArrays           &
          (ierr,                                &
           GH,                                  &
           N_dims,                              &
           local_interp_handle,                 &
           param_table_handle,                  &
           coord_system_handle,                 &
           N_interp_points,                     &
             interp_coords_type,                &
             interp_coords,                     &
           N_input_arrays,                      &
             input_array_indices,               &
           N_output_arrays,                     &
             output_array_types,                &
             output_arrays)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               N_dims
       integer               local_interp_handle
       integer               param_table_handle
       integer               coord_system_handle
       integer               N_interp_points
       integer                 interp_coords_type
       CCTK_POINTER_TO_CONST   interp_coords(N_dims)
       integer               N_input_arrays
       CCTK_INT                input_array_indices(N_input_arrays)
       integer               N_output_arrays
       CCTK_INT                output_array_types(N_output_arrays)
       CCTK_POINTER            output_arrays(N_output_arrays)
     end subroutine CCTK_InterpGridArrays

     subroutine CCTK_QueryGroupStorage (ierr, GH, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          groupname
     end subroutine CCTK_QueryGroupStorage

     subroutine CCTK_QueryGroupStorageI (ierr, GH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
     end subroutine CCTK_QueryGroupStorageI

! CCTK_ArrayGroupSize is not available from Fortran
! CCTK_ArrayGroupSizeI is not available from Fortran

  end interface

end module cctk_Comm
