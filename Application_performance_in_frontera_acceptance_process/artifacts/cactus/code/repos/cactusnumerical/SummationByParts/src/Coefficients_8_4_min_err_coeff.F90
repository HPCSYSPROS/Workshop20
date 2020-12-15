#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff_8_4_opt ( nsize, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(4), save :: a
  CCTK_REAL, dimension(12,8), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_8_4_opt ( a, q )
    first = .false. 
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1:8,1) = q(1:8,1)
    imin(1) = 1; imax(1) = 8

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

    dd(1:12,8) = q(1:12,8)
    imin(8) = 1; imax(8) = 12

    il = 9
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-11:nsize,nsize-7) = -q(12:1:-1,8)
    imin(nsize-7) = nsize-11; imax(nsize-7) = nsize

    dd(nsize-10:nsize,nsize-6) = -q(11:1:-1,7)
    imin(nsize-6) = nsize-10; imax(nsize-6) = nsize

    dd(nsize-9:nsize,nsize-5) = -q(10:1:-1,6)
    imin(nsize-5) = nsize-9; imax(nsize-5) = nsize

    dd(nsize-8:nsize,nsize-4) = -q(9:1:-1,5)
    imin(nsize-4) = nsize-8; imax(nsize-4) = nsize

    dd(nsize-7:nsize,nsize-3) = -q(8:1:-1,4)
    imin(nsize-3) = nsize-7; imax(nsize-3) = nsize

    dd(nsize-7:nsize,nsize-2) = -q(8:1:-1,3)
    imin(nsize-2) = nsize-7; imax(nsize-2) = nsize

    dd(nsize-7:nsize,nsize-1) = -q(8:1:-1,2)
    imin(nsize-1) = nsize-7; imax(nsize-1) = nsize

    dd(nsize-7:nsize,nsize) = -q(8:1:-1,1)
    imin(nsize) = nsize-7; imax(nsize) = nsize

    ir = nsize - 8
  end if
  do i = il, ir
    dd(i-4:i-1,i) = -a(4:1:-1); dd(i+1:i+4,i) = a(1:4)
    imin(i) = i-4; imax(i) = i+4
  end do 
end subroutine set_coeff_8_4_opt

subroutine set_coeff_up_8_4_opt ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

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
  CCTK_REAL, dimension(12,8), save :: q1, q2

    CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_8_4_opt called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_8_4_opt ( a1, q1, a2, q2 )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    if ( dir == -1 ) then
      dd(1:8,1) = q1(1:8,1)
      imin(1) = 1; imax(1) = 8

      dd(1:8,2) = q1(1:8,2)
      imin(2) = 1; imax(2) = 8

      dd(1:8,3) = q1(1:8,3)
      imin(3) = 1; imax(3) = 8

      dd(1:8,4) = q1(1:8,4)
      imin(4) = 1; imax(4) = 8

      dd(1:8,5) = q1(1:8,5)
      imin(5) = 1; imax(5) = 8

      dd(1:9,6) = q1(1:9,6)
      imin(6) = 1; imax(6) = 9

      dd(1:10,7) = q1(1:10,7)
      imin(7) = 1; imax(7) = 10

      dd(1:11,8) = q1(1:11,8)
      imin(8) = 1; imax(8) = 11
    else
      dd(1:8,1) = q2(1:8,1)
      imin(1) = 1; imax(1) = 8

      dd(1:8,2) = q2(1:8,2)
      imin(2) = 1; imax(2) = 8

      dd(1:8,3) = q2(1:8,3)
      imin(3) = 1; imax(3) = 8

      dd(1:8,4) = q2(1:8,4)
      imin(4) = 1; imax(4) = 8

      dd(1:9,5) = q2(1:9,5)
      imin(5) = 1; imax(5) = 9

      dd(1:10,6) = q2(1:10,6)
      imin(6) = 1; imax(6) = 10

      dd(1:11,7) = q2(1:11,7)
      imin(7) = 1; imax(7) = 11

      dd(1:12,8) = q2(1:12,8)
      imin(8) = 1; imax(8) = 12
    end if
    il = 9
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
   if ( dir == -1 ) then
      dd(nsize-11:nsize,nsize-7) = -q2(12:1:-1,8)
      imin(nsize-7) = nsize-11; imax(nsize-7) = nsize

      dd(nsize-10:nsize,nsize-6) = -q2(11:1:-1,7)
      imin(nsize-6) = nsize-10; imax(nsize-6) = nsize

      dd(nsize-9:nsize,nsize-5) = -q2(10:1:-1,6)
      imin(nsize-5) = nsize-9; imax(nsize-5) = nsize

      dd(nsize-8:nsize,nsize-4) = -q2(9:1:-1,5)
      imin(nsize-4) = nsize-8; imax(nsize-4) = nsize

      dd(nsize-7:nsize,nsize-3) = -q2(8:1:-1,4)
      imin(nsize-3) = nsize-7; imax(nsize-3) = nsize

      dd(nsize-7:nsize,nsize-2) = -q2(8:1:-1,3)
      imin(nsize-2) = nsize-7; imax(nsize-2) = nsize

      dd(nsize-7:nsize,nsize-1) = -q2(8:1:-1,2)
      imin(nsize-1) = nsize-7; imax(nsize-1) = nsize

      dd(nsize-7:nsize,nsize) = -q2(8:1:-1,1)
      imin(nsize) = nsize-7; imax(nsize) = nsize
    else
      dd(nsize-10:nsize,nsize-7) = -q1(11:1:-1,8)
      imin(nsize-7) = nsize-10; imax(nsize-7) = nsize

      dd(nsize-9:nsize,nsize-6) = -q1(10:1:-1,7)
      imin(nsize-6) = nsize-9; imax(nsize-6) = nsize

      dd(nsize-8:nsize,nsize-5) = -q1(9:1:-1,6)
      imin(nsize-5) = nsize-8; imax(nsize-5) = nsize

      dd(nsize-7:nsize,nsize-4) = -q1(8:1:-1,5)
      imin(nsize-4) = nsize-7; imax(nsize-4) = nsize

      dd(nsize-7:nsize,nsize-3) = -q1(8:1:-1,4)
      imin(nsize-3) = nsize-7; imax(nsize-3) = nsize

      dd(nsize-7:nsize,nsize-2) = -q1(8:1:-1,3)
      imin(nsize-2) = nsize-7; imax(nsize-2) = nsize

      dd(nsize-7:nsize,nsize-1) = -q1(8:1:-1,2)
      imin(nsize-1) = nsize-7; imax(nsize-1) = nsize

      dd(nsize-7:nsize,nsize) = -q1(8:1:-1,1)
      imin(nsize) = nsize-7; imax(nsize) = nsize
    end if
    ir = nsize - 8
  end if

  do i = il, ir
    if ( dir == -1 ) then
      dd(i-4:i+3,i) = a1(-4:3)
      imin(i) = i-4; imax(i) = i+3
    else
      dd(i-3:i+4,i) = a2(-3:4)
      imin(i) = i-3; imax(i) = i+4
    end if
  end do
end subroutine set_coeff_up_8_4_opt

subroutine set_coeff_up_8 ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, dir, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(-5:5), save :: a1, a2
  CCTK_REAL, dimension(7,5), save :: q1, q2

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_8 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_8 ( a1, q1, a2, q2 )
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

      dd(1:7,5) = q1(1:7,5);
      imin(5) = 1; imax(5) = 7
    else
      dd(1,1) = q2(1,1);
      imin(1) = 1; imax(1) = 1

      dd(1:5,2) = q2(3:7,4);
      imin(2) = 1; imax(2) = 5

      dd(1:7,3) = q2(1:7,5);
      imin(3) = 1; imax(3) = 7

      dd(1:9,4) = a2(-3:5);
      imin(4) = 1; imax(4) = 9

      dd(2:10,5) = a2(-3:5);
      imin(5) = 2; imax(5) = 10
    end if
    il = 6
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize,nsize) = q1(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-4:nsize,nsize-1) = q1(1:5,4)
      imin(nsize-1) = nsize-4; imax(nsize-1) = nsize

      dd(nsize-6:nsize,nsize-2) = q1(1:7,5)
      imin(nsize-2) = nsize-6; imax(nsize-2) = nsize

      dd(nsize-8:nsize,nsize-3) = a1(-5:3)
      imin(nsize-3) = nsize-8; imax(nsize-3) = nsize

      dd(nsize-9:nsize-1,nsize-4) = a1(-5:3)
      imin(nsize-4) = nsize-9; imax(nsize-4) = nsize-1
    else
      dd(nsize,nsize) = q2(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-1:nsize,nsize-1) = q2(6:7,2)
      imin(nsize-1) = nsize-1; imax(nsize-1) = nsize

      dd(nsize-2:nsize,nsize-2) = q2(5:7,3)
      imin(nsize-2) = nsize-2; imax(nsize-2) = nsize

      dd(nsize-4:nsize,nsize-3) = q2(3:7,4)
      imin(nsize-3) = nsize-4; imax(nsize-3) = nsize

      dd(nsize-6:nsize,nsize-4) = q2(1:7,5)
      imin(nsize-4) = nsize-6; imax(nsize-4) = nsize
    end if
    ir = nsize - 5
  end if
  do i = il, ir
    if ( dir == -1 ) then
      dd(i-5:i+3,i) = a1(-5:3)
      imin(i) = i-5; imax(i) = i+3
    else
      dd(i-3:i+5,i) = a2(-3:5)
      imin(i) = i-3; imax(i) = i+5
    end if
  end do
end subroutine set_coeff_up_8
