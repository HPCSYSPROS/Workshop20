#include "cctk.h"

! This is the name of a subroutine, so the module must be called differently.
module cctk_Malloc1
  implicit none

  interface

     subroutine CCTK_MemStat
       implicit none
     end subroutine CCTK_MemStat

     subroutine CCTK_TotalMemory (amount)
       implicit none
       integer amount
     end subroutine CCTK_TotalMemory

     subroutine CCTK_MemTicketCash (diff, this_ticket)
       implicit none
       integer diff
       integer this_ticket
     end subroutine CCTK_MemTicketCash

     subroutine CCTK_MemTicketRequest (ticket)
       implicit none
       integer ticket
     end subroutine CCTK_MemTicketRequest

     subroutine CCTK_MemTicketDelete (ierr, this_ticket)
       implicit none
       integer ierr
       integer this_ticket
     end subroutine CCTK_MemTicketDelete

  end interface
  
end module cctk_Malloc1
