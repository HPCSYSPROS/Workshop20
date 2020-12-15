#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff_6 ( nsize, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(3), save :: a
  CCTK_REAL, dimension(5,3), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_6 ( a, q )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1:5,1) = q(1:5,1)
    imin(1) = 1; imax(1) = 1

    dd(1:5,2) = q(1:5,2)
    imin(2) = 1; imax(2) = 3

    dd(1:5,3) = q(1:5,3)
    imin(3) = 1; imax(3) = 5

    il = 4
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-4:nsize,nsize-2) = -q(5:1:-1,3)
    imin(nsize-2) = nsize-4; imax(nsize-2) = nsize

    dd(nsize-4:nsize,nsize-1) = -q(5:1:-1,2)
    imin(nsize-1) = nsize-2; imax(nsize-1) = nsize

    dd(nsize-4:nsize,nsize) = -q(5:1:-1,1)
    imin(nsize) = nsize; imax(nsize) = nsize

    ir = nsize - 3
  end if
  do i = il, ir
    dd(i-3:i-1,i) = -a(3:1:-1);
    dd(i+1:i+3,i) = a(1:3)
    imin(i) = i-3; imax(i) = i+3
  end do
end subroutine set_coeff_6
