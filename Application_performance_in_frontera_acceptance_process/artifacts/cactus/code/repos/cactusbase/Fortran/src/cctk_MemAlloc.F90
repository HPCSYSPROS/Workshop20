#include "cctk.h"

module cctk_MemAlloc
  implicit none

  interface

     subroutine CCTK_Malloc (ptr, size, line, file)
       implicit none
       CCTK_POINTER ptr
       integer      size
       integer      line
       character(*) file
     end subroutine CCTK_Malloc

     subroutine CCTK_Free (ptr)
       implicit none
       CCTK_POINTER ptr
     end subroutine CCTK_Free

     !  CCTK_TotalMemory is already declared in cctk_Malloc1
     
  end interface
  
end module cctk_MemAlloc
