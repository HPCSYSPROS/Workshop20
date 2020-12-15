#include "cctk.h"

module cctk_Flesh
  implicit none

  interface

     subroutine CCTK_Traverse (ierr, GH, where)
       implicit none
       integer      ierr
       CCTK_POINTER GH
       character(*) where
     end subroutine CCTK_Traverse
     
  end interface
  
end module cctk_Flesh
