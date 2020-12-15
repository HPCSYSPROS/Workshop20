#include "cctk.h"

module cctk_ActiveThorns
  implicit none

  interface

! This routine has been made a function instead of a subroutine
! and is declared in cctk.h
!     subroutine CCTK_IsThornActive (ierr, name)
!       implicit none
!       integer      ierr
!       character(*) name
!     end subroutine CCTK_IsThornActive

     subroutine CCTK_ThornHandle (handle, name)
       implicit none
       CCTK_POINTER_TO_CONST handle
       character(*) name
     end subroutine CCTK_ThornHandle

     subroutine CCTK_IsThornActiveH (ierr, handle)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST handle
     end subroutine CCTK_IsThornActiveH

     subroutine CCTK_IsThornCompiled (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_IsThornCompiled

     subroutine CCTK_IsImplementationActive (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_IsImplementationActive

     subroutine CCTK_IsImplementationCompiled (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_IsImplementationCompiled

     subroutine CCTK_ActivatingThorn (thorn, imp)
       implicit none
       character(*) thorn
       character(*) imp
     end subroutine CCTK_ActivatingThorn

     subroutine CCTK_ImpThornList (thornlist, imp)
       implicit none
       CCTK_POINTER thornlist
       character(*) imp
     end subroutine CCTK_ImpThornList

     subroutine CCTK_ThornImplementation (imp, thorn)
       implicit none
       character(*) imp
       character(*) thorn
     end subroutine CCTK_ThornImplementation

     subroutine CCTK_ImplementationThorn (thorn, imp)
       implicit none
       character(*) thorn
       character(*) imp
     end subroutine CCTK_ImplementationThorn

     subroutine CCTK_ImplementationRequires (implist, imp)
       implicit none
       CCTK_POINTER implist
       character(*) imp
     end subroutine CCTK_ImplementationRequires

     subroutine CCTK_NumCompiledThorns (num)
       implicit none
       integer num
     end subroutine CCTK_NumCompiledThorns

     subroutine CCTK_NumCompiledImplementations (num)
       implicit none
       integer num
     end subroutine CCTK_NumCompiledImplementations

     subroutine CCTK_CompiledThorn (thorn_name, thorn_index)
       implicit none
       character(*) thorn_name
       integer      thorn_index
     end subroutine CCTK_CompiledThorn

     subroutine CCTK_CompiledImplementation (imp_name, imp_index)
       implicit none
       character(*) imp_name
       integer      imp_index
     end subroutine CCTK_CompiledImplementation

  end interface

end module cctk_ActiveThorns
