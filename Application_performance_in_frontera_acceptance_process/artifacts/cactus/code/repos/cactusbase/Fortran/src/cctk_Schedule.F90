#include "cctk.h"

module cctk_Schedule
  implicit none

  interface

     subroutine CCTK_SchedulePrint (ierr, where)
       implicit none
       integer      ierr
       character(*) where
     end subroutine CCTK_SchedulePrint

     subroutine CCTK_SchedulePrintTimes (ierr, where)
       implicit none
       integer      ierr
       character(*) where
     end subroutine CCTK_SchedulePrintTimes
     
  end interface
  
end module cctk_Schedule
