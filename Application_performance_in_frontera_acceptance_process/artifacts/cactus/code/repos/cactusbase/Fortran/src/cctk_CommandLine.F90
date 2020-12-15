#include "cctk.h"

module cctk_CommandLine
  implicit none

  interface

! This subroutine has been renamed
! because it conflicts with the name of this module
     subroutine CCTK_GetCommandLine (argc, outargv)
       implicit none
       integer      argc
       CCTK_POINTER outargv
     end subroutine CCTK_GetCommandLine

     subroutine CCTK_ParameterFileName (ierr, filename, filenamelen)
       implicit none
       integer      ierr
       character(*) filename
       integer      filenamelen
     end subroutine CCTK_ParameterFileName

  end interface

end module cctk_CommandLine
