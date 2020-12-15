#include "cctk.h"

module cctk_Parameter
  implicit none

  interface

     subroutine CCTK_ParameterLevel (level)
       implicit none
       integer level
     end subroutine CCTK_ParameterLevel

     subroutine CCTK_ParameterSet (ierr, name, thorn, value)
       implicit none
       integer      ierr
       character(*) name
       character(*) thorn
       character(*) value
     end subroutine CCTK_ParameterSet

     subroutine CCTK_ParameterGet (param, name, thorn, type)
       implicit none
       CCTK_POINTER_TO_CONST param
       character(*)          name
       character(*)          thorn
       integer               type
     end subroutine CCTK_ParameterGet

     subroutine CCTK_ParameterValString (nchars, param, thorn, value)
       implicit none
       integer      nchars
       character(*) param
       character(*) thorn
       character(*) value
     end subroutine CCTK_ParameterValString
  
     subroutine CCTK_ParameterWalk (ierr, first, origin, pfullname, pdata)
       implicit none
       integer      ierr
       integer      first
       character(*) origin
       CCTK_POINTER pfullname
       CCTK_POINTER pdata
     end subroutine CCTK_ParameterWalk

     subroutine CCTK_ParamtereDate (data, name, thorn)
       implicit none
       CCTK_POINTER_TO_CONST data
       character(*)          name
       character(*)          thorn
     end subroutine CCTK_ParamtereDate

     subroutine CCTK_ParameterQueryTimesSet (ierr, name, thorn)
       implicit none
       integer ierr
       character(*) name
       character(*) thorn
     end subroutine CCTK_ParameterQueryTimesSet

  end interface

end module cctk_Parameter
