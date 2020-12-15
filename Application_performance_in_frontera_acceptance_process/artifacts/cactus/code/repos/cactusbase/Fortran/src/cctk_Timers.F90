#include "cctk.h"

module cctk_Timers
  implicit none

  interface

     subroutine CCTK_NumTimers (num_timers)
       implicit none
       integer num_timers
     end subroutine CCTK_NumTimers

     subroutine CCTK_NumClocks (num_clocks)
       implicit none
       integer num_clocks
     end subroutine CCTK_NumClocks
     
     subroutine CCTK_TimerName (timer_name, timer_length, timer_handle)
       implicit none
       character(*) timer_name
       integer      timer_length
       integer      timer_handle
     end subroutine CCTK_TimerName
     
     subroutine CCTK_ClockName (clock_name, clock_length, clock_handle)
       implicit none
       character(*) clock_name
       integer      clock_length
       integer      clock_handle
     end subroutine CCTK_ClockName
     
     subroutine CCTK_ClockHandle (handle, nclock)
       implicit none
       integer      handle
       character(*) nclock
     end subroutine CCTK_ClockHandle
     
     subroutine CCTK_TimerCreate (handle, name)
       implicit none
       integer      handle
       character(*) name
     end subroutine CCTK_TimerCreate
     
     subroutine CCTK_TimerCreateI (handle)
       implicit none
       integer handle
     end subroutine CCTK_TimerCreateI
     
     subroutine CCTK_TimerDestroy (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_TimerDestroy
     
     subroutine CCTK_TimerDestroyI (ierr, handle)
       implicit none
       integer ierr
       integer handle
     end subroutine CCTK_TimerDestroyI
     
     subroutine CCTK_TimerStart (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_TimerStart
     
     subroutine CCTK_TimerStartI (ierr, handle)
       implicit none
       integer ierr
       integer handle
     end subroutine CCTK_TimerStartI
     
     subroutine CCTK_TimerStop (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_TimerStop
     
     subroutine CCTK_TimerStopI (ierr, handle)
       implicit none
       integer ierr
       integer handle
     end subroutine CCTK_TimerStopI
     
     subroutine CCTK_TimerReset (ierr, name)
       implicit none
       integer      ierr
       character(*) name
     end subroutine CCTK_TimerReset
     
     subroutine CCTK_TimerResetI (ierr, handle)
       implicit none
       integer ierr
       integer handle
     end subroutine CCTK_TimerResetI
     
     subroutine CCTK_TimerPrintData (ierr, name, nclock)
       implicit none
       integer      ierr
       character(*) name
       character(*) nclock
     end subroutine CCTK_TimerPrintData
     
     subroutine CCTK_TimerPrintDataI (ierr, handle, nclock)
       implicit none
       integer      ierr
       integer      handle
       character(*) nclock
     end subroutine CCTK_TimerPrintDataI
     
  end interface
  
end module cctk_Timers
