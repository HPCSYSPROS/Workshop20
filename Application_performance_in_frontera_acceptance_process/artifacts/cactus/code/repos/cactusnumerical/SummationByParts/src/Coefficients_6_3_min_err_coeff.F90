#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff_6_3_opt ( nsize, loc_order, bb, gsize, imin, imax, dd )

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
  CCTK_REAL, dimension(9,6), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_6_3_opt ( a, q )
    first = .false. 
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1:6,1) = q(1:6,1)
    imin(1) = 1; imax(1) = 6

    dd(1:6,2) = q(1:6,2)
    imin(2) = 1; imax(2) = 6

    dd(1:6,3) = q(1:6,3)
    imin(3) = 1; imax(3) = 6

    dd(1:7,4) = q(1:7,4)
    imin(4) = 1; imax(4) = 7

    dd(1:8,5) = q(1:8,5)
    imin(5) = 1; imax(5) = 8

    dd(1:9,6) = q(1:9,6)
    imin(6) = 1; imax(6) = 9

    il = 7
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-8:nsize,nsize-5) = -q(9:1:-1,6)
    imin(nsize-5) = nsize-8; imax(nsize-5) = nsize

    dd(nsize-7:nsize,nsize-4) = -q(8:1:-1,5)
    imin(nsize-4) = nsize-7; imax(nsize-4) = nsize

    dd(nsize-6:nsize,nsize-3) = -q(7:1:-1,4)
    imin(nsize-3) = nsize-6; imax(nsize-3) = nsize

    dd(nsize-5:nsize,nsize-2) = -q(6:1:-1,3)
    imin(nsize-2) = nsize-5; imax(nsize-2) = nsize

    dd(nsize-5:nsize,nsize-1) = -q(6:1:-1,2)
    imin(nsize-1) = nsize-5; imax(nsize-1) = nsize

    dd(nsize-5:nsize,nsize) = -q(6:1:-1,1)
    imin(nsize) = nsize-5; imax(nsize) = nsize

    ir = nsize - 6
  end if
  do i = il, ir
    dd(i-3:i-1,i) = -a(3:1:-1); dd(i+1:i+3,i) = a(1:3)
    imin(i) = i-3; imax(i) = i+3
  end do
end subroutine set_coeff_6_3_opt

subroutine set_coeff_up_6_3_opt ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, dir, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(-3:3), save :: a1, a2
  CCTK_REAL, dimension(9,6), save :: q1, q2

    CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_6_3 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_6_3_opt ( a1, q1, a2, q2 )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    if ( dir == -1 ) then
      dd(1:6,1) = q1(1:6,1)
      imin(1) = 1; imax(1) = 6

      dd(1:6,2) = q1(1:6,2)
      imin(2) = 1; imax(2) = 6

      dd(1:6,3) = q1(1:6,3)
      imin(3) = 1; imax(3) = 6

      dd(1:6,4) = q1(1:6,4)
      imin(4) = 1; imax(4) = 6

      dd(1:7,5) = q1(1:7,5)
      imin(5) = 1; imax(5) = 7

      dd(1:8,6) = q1(1:8,6)
      imin(6) = 1; imax(6) = 8
    else
      dd(1:6,1) = q2(1:6,1)
      imin(1) = 1; imax(1) = 6

      dd(1:6,2) = q2(1:6,2)
      imin(2) = 1; imax(2) = 6

      dd(1:6,3) = q2(1:6,3)
      imin(3) = 1; imax(3) = 6

      dd(1:7,4) = q2(1:7,4)
      imin(4) = 1; imax(4) = 7

      dd(1:8,5) = q2(1:8,5)
      imin(5) = 1; imax(5) = 8

      dd(1:9,6) = q2(1:9,6)
      imin(6) = 1; imax(6) = 9
    end if
    il = 7
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize-8:nsize,nsize-5) = -q2(9:1:-1,6)
      imin(nsize-5) = nsize-8; imax(nsize-5) = nsize

      dd(nsize-7:nsize,nsize-4) = -q2(8:1:-1,5)
      imin(nsize-4) = nsize-7; imax(nsize-4) = nsize

      dd(nsize-6:nsize,nsize-3) = -q2(7:1:-1,4)
      imin(nsize-3) = nsize-6; imax(nsize-3) = nsize

      dd(nsize-5:nsize,nsize-2) = -q2(6:1:-1,3)
      imin(nsize-2) = nsize-5; imax(nsize-2) = nsize

      dd(nsize-5:nsize,nsize-1) = -q2(6:1:-1,2)
      imin(nsize-1) = nsize-5; imax(nsize-1) = nsize

      dd(nsize-5:nsize,nsize) = -q2(6:1:-1,1)
      imin(nsize) = nsize-5; imax(nsize) = nsize
    else
      dd(nsize-7:nsize,nsize-5) = -q1(8:1:-1,6)
      imin(nsize-5) = nsize-7; imax(nsize-5) = nsize

      dd(nsize-6:nsize,nsize-4) = -q1(7:1:-1,5)
      imin(nsize-4) = nsize-6; imax(nsize-4) = nsize

      dd(nsize-5:nsize,nsize-3) = -q1(6:1:-1,4)
      imin(nsize-3) = nsize-5; imax(nsize-3) = nsize

      dd(nsize-5:nsize,nsize-2) = -q1(6:1:-1,3)
      imin(nsize-2) = nsize-5; imax(nsize-2) = nsize

      dd(nsize-5:nsize,nsize-1) = -q1(6:1:-1,2)
      imin(nsize-1) = nsize-5; imax(nsize-1) = nsize

      dd(nsize-5:nsize,nsize) = -q1(6:1:-1,1)
      imin(nsize) = nsize-5; imax(nsize) = nsize
    end if
    ir = nsize - 6
  end if

  do i = il, ir
    if ( dir == -1 ) then
      dd(i-3:i+2,i) = a1(-3:2)
      imin(i) = i-3; imax(i) = i+2
    else
      dd(i-2:i+3,i) = a2(-2:3)
      imin(i) = i-2; imax(i) = i+3
    end if
  end do
end subroutine set_coeff_up_6_3_opt

subroutine set_coeff_up_6 ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, dir, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(-4:4), save :: a1, a2
  CCTK_REAL, dimension(5,4), save :: q1, q2

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_2 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_6 ( a1, q1, a2, q2 )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    if ( dir == -1 ) then
      dd(1,1) = q1(1,1);
      imin(1) = 1; imax(1) = 1

      dd(1:2,2) = q1(1:2,2);
      imin(2) = 1; imax(2) = 2

      dd(1:3,3) = q1(1:3,3);
      imin(3) = 1; imax(3) = 3

      dd(1:5,4) = q1(1:5,4);
      imin(4) = 1; imax(4) = 5
    else
      dd(1,1) = q2(1,1);
      imin(1) = 1; imax(1) = 1

      dd(1:5,2) = q2(1:5,4);
      imin(2) = 1; imax(2) = 5

      dd(1:7,3) = a2(-2:4);
      imin(3) = 1; imax(3) = 7

      dd(2:8,4) = a2(-2:4);
      imin(4) = 2; imax(4) = 8
    end if
    il = 5
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize,nsize) = q1(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-4:nsize,nsize-1) = q1(1:5,4)
      imin(nsize-1) = nsize-4; imax(nsize-1) = nsize

      dd(nsize-6:nsize,nsize-2) = a1(-4:2)
      imin(nsize-2) = nsize-6; imax(nsize-2) = nsize

      dd(nsize-7:nsize-1,nsize-3) = a1(-4:2)
      imin(nsize-3) = nsize-7; imax(nsize-3) = nsize-1
    else
      dd(nsize,nsize) = q2(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-1:nsize,nsize-1) = q2(4:5,2)
      imin(nsize-1) = nsize-1; imax(nsize-1) = nsize

      dd(nsize-2:nsize,nsize-2) = q2(3:5,3)
      imin(nsize-2) = nsize-2; imax(nsize-2) = nsize

      dd(nsize-4:nsize,nsize-3) = q2(1:5,4)
      imin(nsize-3) = nsize-4; imax(nsize-3) = nsize
    end if
    ir = nsize - 4
  end if
  do i = il, ir
    if ( dir == -1 ) then
      dd(i-4:i+2,i) = a1(-4:2)
      imin(i) = i-4; imax(i) = i+2
    else
      dd(i-2:i+4,i) = a2(-2:4)
      imin(i) = i-2; imax(i) = i+4
    end if
  end do
end subroutine set_coeff_up_6
