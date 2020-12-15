#include "cctk.h"

module cctk_WarnLevel
  implicit none

  interface

     subroutine CCTK_Warn (level, line, file, thorn, message)
       implicit none
       integer      level
       integer      line
       character(*) file
       character(*) thorn
       character(*) message
     end subroutine CCTK_Warn

     subroutine CCTK_Error (line, file, thorn, message)
       implicit none
       integer      line
       character(*) file
       character(*) thorn
       character(*) message
     end subroutine CCTK_Error

     subroutine CCTK_ParamWarn (thorn, message)
       implicit none
       character(*) thorn
       character(*) message
     end subroutine CCTK_ParamWarn

     subroutine CCTK_Info (thorn, message)
       implicit none
       character(*) thorn
       character(*) message
     end subroutine CCTK_Info
     
  end interface
  
  ! Do not declare these, because some compilers insist on a
  ! definition once they have seen an external declaration:
!!$  external CCTK_VWarn
!!$  external CCTK_VError
!!$  external CCTK_VParamWarn
!!$  external CCTK_VInfo
  
end module cctk_WarnLevel
