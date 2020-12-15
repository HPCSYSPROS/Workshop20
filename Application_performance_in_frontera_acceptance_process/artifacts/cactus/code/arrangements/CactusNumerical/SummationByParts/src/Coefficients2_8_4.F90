#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff2_8_4 ( nsize, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(5), save :: a
  CCTK_REAL, dimension(12,8), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_2_8_4 ( a, q )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1:7,1) = q(1:7,1)
    imin(1) = 1; imax(1) = 7

    dd(1:8,2) = q(1:8,2)
    imin(2) = 1; imax(2) = 8

    dd(1:8,3) = q(1:8,3)
    imin(3) = 1; imax(3) = 8

    dd(1:8,4) = q(1:8,4)
    imin(4) = 1; imax(4) = 8

    dd(1:9,5) = q(1:9,5)
    imin(5) = 1; imax(5) = 9

    dd(1:10,6) = q(1:10,6)
    imin(6) = 1; imax(6) = 10

    dd(1:11,7) = q(1:11,7)
    imin(7) = 1; imax(7) = 11

    dd(2:12,8) = q(2:12,8)
    imin(8) = 2; imax(8) = 12

    il = 9
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-11:nsize-1,nsize-7) = q(12:2:-1,8)
    imin(nsize-7) = nsize-11; imax(nsize-7) = nsize-1

    dd(nsize-10:nsize,nsize-6) = q(11:1:-1,7)
    imin(nsize-6) = nsize-10; imax(nsize-6) = nsize

    dd(nsize-9:nsize,nsize-5) = q(10:1:-1,6)
    imin(nsize-5) = nsize-9; imax(nsize-5) = nsize

    dd(nsize-8:nsize,nsize-4) = q(9:1:-1,5)
    imin(nsize-4) = nsize-8; imax(nsize-4) = nsize

    dd(nsize-7:nsize,nsize-3) = q(8:1:-1,4)
    imin(nsize-3) = nsize-7; imax(nsize-3) = nsize

    dd(nsize-7:nsize,nsize-2) = q(8:1:-1,3)
    imin(nsize-2) = nsize-7; imax(nsize-2) = nsize

    dd(nsize-7:nsize,nsize-1) = q(8:1:-1,2)
    imin(nsize-1) = nsize-7; imax(nsize-1) = nsize

    dd(nsize-6:nsize,nsize) = q(7:1:-1,1)
    imin(nsize) = nsize-6; imax(nsize) = nsize

    ir = nsize - 8
  end if
  do i = il, ir
    dd(i-4:i-1,i) = a(5:2:-1); dd(i:i+4,i) = a(1:5)
    imin(i) = i-4; imax(i) = i+4
  end do
end subroutine set_coeff2_8_4
