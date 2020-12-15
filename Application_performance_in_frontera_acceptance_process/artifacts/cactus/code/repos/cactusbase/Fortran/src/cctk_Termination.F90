#include "cctk.h"

module cctk_Termination
  implicit none

  interface

     subroutine CCTK_TerminationReached (ires, cctkGH)
       implicit none
       integer               ires
       CCTK_POINTER_TO_CONST cctkGH
     end subroutine CCTK_TerminationReached
     
     subroutine CCTK_TerminateNext (cctkGH)
       implicit none
       CCTK_POINTER_TO_CONST cctkGH
     end subroutine CCTK_TerminateNext
     
  end interface
  
end module cctk_Termination
