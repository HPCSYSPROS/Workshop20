#include "cctk.h"

module cctk_File
  implicit none
  
  interface
     
     subroutine CCTK_CreateDirectory (ierr, mode, pathname)
       implicit none
       integer      ierr
       integer      mode
       character(*) pathname
     end subroutine CCTK_CreateDirectory
     
  end interface
  
end module cctk_File
