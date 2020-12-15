#include "cctk.h"

module cctk_FortranWrappers
  implicit none

  interface

     subroutine CCTK_RegisterFortranWrapper (ierr, name, function)
       implicit none
       integer       ierr
       character(*)  name
       CCTK_FPOINTER function
     end subroutine CCTK_RegisterFortranWrapper

     subroutine CCTK_FortranWrapper (function, name)
       implicit none
       CCTK_FPOINTER function
       character(*)  name
     end subroutine CCTK_FortranWrapper

  end interface

end module cctk_FortranWrappers
