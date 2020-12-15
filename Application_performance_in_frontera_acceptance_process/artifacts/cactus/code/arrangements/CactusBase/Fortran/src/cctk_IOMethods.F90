#include "cctk.h"

module cctk_IOMethods
  implicit none

  interface

     subroutine CCTKi_RegisterIOMethod (handle, thorn, name)
       implicit none
       integer handle
       character(*) thorn
       character(*) name
     end subroutine CCTKi_RegisterIOMethod

     subroutine CCTK_RegisterIOMethodOutputGH (ierr, handle, OutputGH)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER OutputGH
     end subroutine CCTK_RegisterIOMethodOutputGH

#if 0
C    no Fortran 90 interface for following routines
C    their names exceed the maximum allowed length of 31 characters
     subroutine CCTK_RegisterIOMethodTimeToOutput (ierr, handle, TimeToOutput)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER TimeToOutput
     end subroutine CCTK_RegisterIOMethodTimeToOutput

     subroutine CCTK_RegisterIOMethodTriggerOutput (ierr, handle, TriggerOutput)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER TriggerOutput
     end subroutine CCTK_RegisterIOMethodTriggerOutput

     subroutine CCTK_RegisterIOMethodOutputVarAs (ierr, handle, OutputVarAs)
       implicit none
       integer       ierr
       integer       handle
       CCTK_FPOINTER OutputVarAs
     end subroutine CCTK_RegisterIOMethodOutputVarAs
#endif

     subroutine CCTK_IOMethodImplementation (nchars, imp, handle)
       implicit none
       integer      nchars
       character(*) imp
       integer      handle
     end subroutine CCTK_IOMethodImplementation

     subroutine CCTK_IOMethod (nchars, method, handle)
       implicit none
       integer      nchars
       character(*) method
       integer      handle
     end subroutine CCTK_IOMethod

     subroutine CCTK_NumIOMethods (num)
       implicit none
       integer num
     end subroutine CCTK_NumIOMethods

  end interface

end module cctk_IOMethods
