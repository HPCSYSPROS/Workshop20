#include "cctk.h"

module cctk_IO
  implicit none

  interface

     subroutine CCTK_OutputGH (ierr, GH)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
     end subroutine CCTK_OutputGH

     subroutine CCTK_OutputVarAs (ierr, GH, var, alias)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          var
       character(*)          alias
     end subroutine CCTK_OutputVarAs

     subroutine CCTK_OutputVar (ierr, GH, var)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          var
     end subroutine CCTK_OutputVar

     subroutine CCTK_OutputVarByMethod (ierr, GH, var, method)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          var
       character(*)          method
     end subroutine CCTK_OutputVarByMethod
     
  end interface
  
end module cctk_IO
