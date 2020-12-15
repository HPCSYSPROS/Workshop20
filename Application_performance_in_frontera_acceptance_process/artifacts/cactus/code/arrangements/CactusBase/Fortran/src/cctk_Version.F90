#include "cctk.h"

module cctk_Version
  implicit none

  interface
     
     subroutine CCTK_FullVersion (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTK_FullVersion
     
     subroutine CCTK_MajorVersion (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTK_MajorVersion
     
     subroutine CCTK_MinorVersion (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTK_MinorVersion
     
     subroutine CCTK_MicroVersion (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTK_MicroVersion
     
     subroutine CCTKi_DateStamp (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTKi_DateStamp
     
     subroutine CCTK_CompileTime (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTK_CompileTime
     
     subroutine CCTKi_CompileDate (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTKi_CompileDate
     
     subroutine CCTKi_CompileUser (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTKi_CompileUser
     
     subroutine CCTKi_RunUser (string, string_length)
       implicit none
       character(*) string
       integer      string_length
     end subroutine CCTKi_RunUser
     
  end interface
  
end module cctk_Version
