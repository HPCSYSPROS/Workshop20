! Check if there are enough grid cells and ghost zones to compute the finite
! differences.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine SBP_CheckGridSizes(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  integer :: i
  integer, dimension(3) :: bc
  integer, dimension(3,3) :: np
  character(len=128) :: info_message
  character(len=1), dimension(3) :: direction = (/ 'x', 'y', 'z' /)

  bc(1) = cctk_bbox(1)+cctk_bbox(2)+1    ! 0 processor boundaries -> bc = 3
  bc(2) = cctk_bbox(3)+cctk_bbox(4)+1    ! 1 processor boundary -> bc = 2
  bc(3) = cctk_bbox(5)+cctk_bbox(6)+1    ! 2 processor boundaries -> bc = 1

  if ( ( use_dissipation > 0 ) .and. &
       ( CCTK_EQUALS(dissipation_type,"Kreiss-Oliger" ) ) .and. &
       ( any (cctk_nghostzones < order / 2 + 1) ) .and. &
       ( CCTK_nProcs(cctkGH) > 1 ) ) then
    write(info_message,'(a22,i1,a13,i2,a32)') 'You need ghostsize >= ', &
                                    order/2+1, ' to run with ', order+2, &
                                    ' order Kreiss-Oliger dissipation'
    call CCTK_WARN(0,info_message) 
  end if

  if ( any (cctk_nghostzones < order / 2) .and. CCTK_nProcs(cctkGH) > 1 )  then
    write(info_message,'(a22,i1,a13,i1,a25)') 'You need ghostsize >= ', &
                                    order/2, ' to run with ', order, &
                                    ' order finite differences'
    call CCTK_WARN(0,info_message)
  end if

  if ( CCTK_EQUALS(norm_type,"Diagonal") ) then
    np(1,:) = 2*cctk_nghostzones
    np(2,:) = 3*order/2
    np(3,:) = 2*order
  else
    np(1,:) = 2*cctk_nghostzones
    np(2,:) = 3*order/2 + 1
    np(3,:) = 2*(order+1)
  end if
    
  do i = 1, 3 
    if ( cctk_lsh(i) < np(bc(i),i) ) then
      write(info_message,'(a1,a35,i1,a25)') &

         direction(i), '-direction is to small to run with ', order, &
         ' order finite differences'
      ! This is a warning only since this thorn may be active, but
      ! never called to calculate finite differences
      call CCTK_WARN(1,info_message)
    end if
  end do
end subroutine SBP_CheckGridSizes
