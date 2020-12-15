#include "cctk.h"

module cctk_GHExtensions
  implicit none

  interface

     subroutine CCTK_RegisterGHExtension (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_RegisterGHExtension

     subroutine CCTK_UnregisterGHExtension (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_UnregisterGHExtension

     subroutine CCTK_RegisterGHExtensionSetupGH (ierr, handle, func)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER func
     end subroutine CCTK_RegisterGHExtensionSetupGH

     subroutine CCTK_RegisterGHExtensionInitGH (ierr, handle, func)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER func
     end subroutine CCTK_RegisterGHExtensionInitGH

#if 0
!    no Fortran 90 interface for following routines
!    their names exceed the maximum allowed length of 31 characters
     subroutine CCTK_RegisterGHExtensionScheduleTraverseGH (ierr, handle, func)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER func
     end subroutine CCTK_RegisterGHExtensionScheduleTraverseGH
#endif

     subroutine CCTK_GHExtensionHandle (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_GHExtensionHandle

     subroutine CCTK_GHExtension (extension, GH, name)
       implicit none
       CCTK_POINTER          extension
       CCTK_POINTER_TO_CONST GH
       character(*)          name
     end subroutine CCTK_GHExtension

  end interface

end module cctk_GHExtensions
