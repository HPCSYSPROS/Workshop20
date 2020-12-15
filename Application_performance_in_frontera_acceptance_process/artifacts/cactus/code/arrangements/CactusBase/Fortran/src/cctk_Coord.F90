#include "cctk.h"

module cctk_Coord
  implicit none

  interface

     subroutine CCTK_CoordDir (dir, name, systemname)
       implicit none
       integer      dir
       character(*) name
       character(*) systemname
     end subroutine CCTK_CoordDir

     subroutine CCTK_CoordIndex (ierr, dir, name, systemname)
       implicit none
       integer      ierr
       integer      dir
       character(*) name
       character(*) systemname
     end subroutine CCTK_CoordIndex

     subroutine CCTK_CoordName (name, namelen, dir, systemname)
       implicit none
       character(*) name
       integer      namelen
       integer      dir
       character(*) systemname
     end subroutine CCTK_CoordName

     subroutine CCTK_CoordRange &
          (ierr, GH, coord_lower, coord_upper, &
          coord_dir, coord_name, system_name)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       CCTK_REAL             coord_lower
       CCTK_REAL             coord_upper
       integer               coord_dir
       character(*)          coord_name
       character(*)          system_name
     end subroutine CCTK_CoordRange

     subroutine CCTK_CoordRangePhysIndex &
          (ierr, GH, lower, upper, coorddir, coordname, systemname)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       integer               lower
       integer               upper
       integer               coorddir
       character(*)          coordname
       character(*)          systemname
     end subroutine CCTK_CoordRangePhysIndex

     subroutine CCTK_CoordRegisterData (ierr, dir, gv, name, systemname)
       implicit none
       integer      ierr
       integer      dir
       character(*) gv
       character(*) name
       character(*) systemname
     end subroutine CCTK_CoordRegisterData

     subroutine CCTK_CoordSystemDim (dim, systemname)
       implicit none
       integer      dim
       character(*) systemname
     end subroutine CCTK_CoordSystemDim

     subroutine CCTK_CoordSystemHandle (handle, systemname)
       implicit none
       integer      handle
       character(*) systemname
     end subroutine CCTK_CoordSystemHandle

     subroutine CCTK_CoordSystemName (name, namelen, handle)
       implicit none
       character(*) name
       integer      namelen
       integer      handle
     end subroutine CCTK_CoordSystemName

     subroutine CCTK_CoordLocalRange &
          (GH, lower, upper, coord_dir, coord_name, system_name)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       CCTK_REAL             lower
       CCTK_REAL             upper
       integer               coord_dir
       character(*)          coord_name
       character(*)          system_name
     end subroutine CCTK_CoordLocalRange

     subroutine CCTK_CoordRegisterRange &
          (GH, coord_min, coord_max, coord_dir, coord_name, system_name)
       implicit none
       integer               ierr
       CCTK_POINTER_TO_CONST GH
       CCTK_REAL             coord_min
       CCTK_REAL             coord_max
       integer               coord_dir
       character(*)          coord_name
       character(*)          system_name
     end subroutine CCTK_CoordRegisterRange

! This routine is commented out because its name is too long.
! Therefore it cannot be called from Fortran anyway.
!     subroutine CCTK_CoordRegisterRangePhysIndex &
!          (ierr, GH, coord_min, coord_max, coord_dir, coord_name, system_name)
!       implicit none
!       integer               ierr
!       CCTK_POINTER_TO_CONST GH
!       integer               coord_min
!       integer               coord_max
!       integer               coord_dir
!       character(*)          coord_name
!       character(*)          system_name
!     end subroutine CCTK_CoordRegisterRangePhysIndex

     subroutine CCTK_NumCoordSystems (num)
       implicit none
       integer num
     end subroutine CCTK_NumCoordSystems

     subroutine CCTK_CoordSystemImplementation (imp, implen, handle)
       implicit none
       character(*) imp
       integer      implen
       integer      handle
     end subroutine CCTK_CoordSystemImplementation

  end interface

end module cctk_Coord
