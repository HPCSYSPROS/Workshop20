#include "cctk.h"

module cctk_Main
  implicit none

  interface

     subroutine CCTK_Initialise (ierr, config)
       implicit none
       integer      ierr
       CCTK_POINTER config
     end subroutine CCTK_Initialise

     subroutine CCTK_Evolve (ierr, config)
       implicit none
       integer      ierr
       CCTK_POINTER config
     end subroutine CCTK_Evolve

     subroutine CCTK_Shutdown (ierr, config)
       implicit none
       integer      ierr
       CCTK_POINTER config
     end subroutine CCTK_Shutdown

     subroutine CCTK_MainLoopIndex (main_loop_index)
       implicit none
       integer main_loop_index
     end subroutine CCTK_MainLoopIndex

     subroutine CCTK_SetMainLoopIndex (ierr, main_loop_index)
       implicit none
       integer ierr
       integer main_loop_index
     end subroutine CCTK_SetMainLoopIndex

  end interface
  
end module cctk_Main
