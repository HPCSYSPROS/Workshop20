#include "cctk.h"

module cctk_Sync
  implicit none

  interface

     subroutine CCTK_SyncGroupI (ierr, cctkGH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST cctkGH
       integer               group
     end subroutine CCTK_SyncGroupI
     
     subroutine CCTK_SyncGroupWithVar (ierr, cctkGH, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST cctkGH
       character(*)          varname
     end subroutine CCTK_SyncGroupWithVar
     
     subroutine CCTK_SyncGroupWithVarI (ierr, cctkGH, var)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST cctkGH
       integer               var
     end subroutine CCTK_SyncGroupWithVarI
     
     subroutine CCTK_SyncGroups (ierr, cctkGH, n_groups, groups)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST cctkGH
       integer               n_groups
       integer               groups(n_groups)
     end subroutine CCTK_SyncGroups
     
  end interface
  
end module cctk_Sync
