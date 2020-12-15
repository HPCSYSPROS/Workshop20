#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff2_2_1 ( nsize, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(2), save :: a
  CCTK_REAL, dimension(3,2), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_2_2_1 ( a, q )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1,1) = q(1,1); dd(2,1) = q(2,1); dd(3,1) = q(3,1)
    imin(1) = 1; imax(1) = 3
    dd(1,2) = q(1,2); dd(2,2) = q(2,2); dd(3,2) = q(3,2)
    imin(2) = 1; imax(2) = 3
    il = 3
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-2,nsize-1) = q(3,2); dd(nsize-1,nsize-1) = q(2,2)
    dd(nsize,nsize-1) = q(1,2)
    imin(nsize-1) = nsize-2; imax(nsize-1) = nsize
    dd(nsize-2,nsize) = q(3,1); dd(nsize-1,nsize) = q(2,1)
    dd(nsize,nsize) = q(1,1)
    imin(nsize) = nsize-2; imax(nsize) = nsize
    ir = nsize - 2
  end if
  do i = il, ir
    dd(i-1,i) = a(2); dd(i:i+1,i) = a(1:2)
    imin(i) = i-1; imax(i) = i+1
  end do

end subroutine set_coeff2_2_1
