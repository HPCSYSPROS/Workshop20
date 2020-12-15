#include "cctk.h"

module cctk_GroupsOnGH
  implicit none

  interface

     subroutine CCTK_VarDataPtr (ptr, GH, timelevel, fullvarname)
       implicit none
       CCTK_POINTER          ptr
       CCTK_POINTER_TO_CONST GH
       integer               timelevel
       character(*)          fullvarname
     end subroutine CCTK_VarDataPtr

     subroutine CCTK_VarDataPtrI (ptr, GH, timelevel, varindex)
       implicit none
       CCTK_POINTER          ptr
       CCTK_POINTER_TO_CONST GH
       integer               timelevel
       integer               varindex
     end subroutine CCTK_VarDataPtrI

     subroutine CCTK_VarDataPtrB (ptr, GH, timelevel, varindex, fullvarname)
       implicit none
       CCTK_POINTER          ptr
       CCTK_POINTER_TO_CONST GH
       integer               timelevel
       integer               varindex
       character(*)          fullvarname
     end subroutine CCTK_VarDataPtrB

     subroutine CCTK_DisableGroupStorageI (ierr, GH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
     end subroutine CCTK_DisableGroupStorageI

     subroutine CCTK_DisableGroupCommI (ierr, GH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
     end subroutine CCTK_DisableGroupCommI

     subroutine CCTK_EnableGroupStorageI (ierr, GH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
     end subroutine CCTK_EnableGroupStorageI

     subroutine CCTK_EnableGroupCommI (ierr, GH, group)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               group
     end subroutine CCTK_EnableGroupCommI

     subroutine CCTK_GrouplbndGN (ierr, GH, dim, lbnd, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lbnd(dim)
       character(*)          groupname
     end subroutine CCTK_GrouplbndGN

     subroutine CCTK_GrouplbndVN (ierr, GH, dim, lbnd, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lbnd(dim)
       character(*)          varname
     end subroutine CCTK_GrouplbndVN

     subroutine CCTK_GrouplbndGI (ierr, GH, dim, lbnd, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lbnd(dim)
       integer               groupindex
     end subroutine CCTK_GrouplbndGI

     subroutine CCTK_GrouplbndVI (ierr, GH, dim, lbnd, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lbnd(dim)
       integer               varindex
     end subroutine CCTK_GrouplbndVI

     subroutine CCTK_GroupubndGN (ierr, GH, dim, ubnd, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               ubnd(dim)
       character(*)          groupname
     end subroutine CCTK_GroupubndGN

     subroutine CCTK_GroupubndVN (ierr, GH, dim, ubnd, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               ubnd(dim)
       character(*)          varname
     end subroutine CCTK_GroupubndVN

     subroutine CCTK_GroupubndGI (ierr, GH, dim, ubnd, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               ubnd(dim)
       integer               groupindex
     end subroutine CCTK_GroupubndGI

     subroutine CCTK_GroupubndVI (ierr, GH, dim, ubnd, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               ubnd(dim)
       integer               varindex
     end subroutine CCTK_GroupubndVI

     subroutine CCTK_GrouplshGN (ierr, GH, dim, lsh, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lsh(dim)
       character(*)          groupname
     end subroutine CCTK_GrouplshGN

     subroutine CCTK_GrouplshVN (ierr, GH, dim, lsh, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lsh(dim)
       character(*)          varname
     end subroutine CCTK_GrouplshVN

     subroutine CCTK_GrouplshGI (ierr, GH, dim, lsh, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lsh(dim)
       integer               groupindex
     end subroutine CCTK_GrouplshGI

     subroutine CCTK_GrouplshVI (ierr, GH, dim, lsh, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               lsh(dim)
       integer               varindex
     end subroutine CCTK_GrouplshVI

     subroutine CCTK_GroupgshGN (ierr, GH, dim, gsh, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               gsh(dim)
       character(*)          groupname
     end subroutine CCTK_GroupgshGN

     subroutine CCTK_GroupgshVN (ierr, GH, dim, gsh, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               gsh(dim)
       character(*)          varname
     end subroutine CCTK_GroupgshVN

     subroutine CCTK_GroupgshGI (ierr, GH, dim, gsh, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               gsh(dim)
       integer               groupindex
     end subroutine CCTK_GroupgshGI

     subroutine CCTK_GroupgshVI (ierr, GH, dim, gsh, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               gsh(dim)
       integer               varindex
     end subroutine CCTK_GroupgshVI

     subroutine CCTK_GroupbboxGN (ierr, GH, size, bbox, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               size
       integer               bbox(size)
       character(*)          groupname
     end subroutine CCTK_GroupbboxGN

     subroutine CCTK_GroupbboxVN (ierr, GH, size, bbox, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               size
       integer               bbox(size)
       character(*)          varname
     end subroutine CCTK_GroupbboxVN

     subroutine CCTK_GroupbboxGI (ierr, GH, size, bbox, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               size
       integer               bbox(size)
       integer               groupindex
     end subroutine CCTK_GroupbboxGI

     subroutine CCTK_GroupbboxVI (ierr, GH, size, bbox, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               size
       integer               bbox(size)
       integer               varindex
     end subroutine CCTK_GroupbboxVI

     subroutine CCTK_GroupnghostzonesGN (ierr, GH, dim, nghostzones, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               nghostzones(dim)
       character(*)          groupname
     end subroutine CCTK_GroupnghostzonesGN

     subroutine CCTK_GroupnghostzonesVN (ierr, GH, dim, nghostzones, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               nghostzones(dim)
       character(*)          varname
     end subroutine CCTK_GroupnghostzonesVN

     subroutine CCTK_GroupnghostzonesGI (ierr, GH, dim, nghostzones, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               nghostzones(dim)
       integer               groupindex
     end subroutine CCTK_GroupnghostzonesGI

     subroutine CCTK_GroupnghostzonesVI (ierr, GH, dim, nghostzones, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               dim
       integer               nghostzones(dim)
       integer               varindex
     end subroutine CCTK_GroupnghostzonesVI

     subroutine CCTK_ActiveTimeLevels (ierr, GH, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          groupname
     end subroutine CCTK_ActiveTimeLevels

     subroutine CCTK_ActiveTimeLevelsGN (ierr, GH, groupname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          groupname
     end subroutine CCTK_ActiveTimeLevelsGN

     subroutine CCTK_ActiveTimeLevelsGI (ierr, GH, groupindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               groupindex
     end subroutine CCTK_ActiveTimeLevelsGI

     subroutine CCTK_ActiveTimeLevelsVN (ierr, GH, varname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       character(*)          varname
     end subroutine CCTK_ActiveTimeLevelsVN

     subroutine CCTK_ActiveTimeLevelsVI (ierr, GH, varindex)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               varindex
     end subroutine CCTK_ActiveTimeLevelsVI

  end interface
  
end module cctk_GroupsOnGH
