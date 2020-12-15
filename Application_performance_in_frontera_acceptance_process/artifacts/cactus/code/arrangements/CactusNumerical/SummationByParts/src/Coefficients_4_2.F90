#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff_4_2 ( nsize, loc_order, bb, gsize, imin, imax, dd )

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
  CCTK_REAL, dimension(6,4), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_4_2 ( a, q )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1:4,1) = q(1:4,1)
    imin(1) = 1; imax(1) = 4

    dd(1:3,2) = q(1:3,2)
    imin(2) = 1; imax(2) = 3

    dd(1:5,3) = q(1:5,3)
    imin(3) = 1; imax(3) = 5

    dd(1:6,4) = q(1:6,4)
    imin(4) = 1; imax(4) = 6

    il = 5
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-5:nsize,nsize-3) = -q(6:1:-1,4)
    imin(nsize-3) = nsize-5; imax(nsize-3) = nsize

    dd(nsize-4:nsize,nsize-2) = -q(5:1:-1,3)
    imin(nsize-2) = nsize-4; imax(nsize-2) = nsize

    dd(nsize-2:nsize,nsize-1) = -q(3:1:-1,2)
    imin(nsize-1) = nsize-2; imax(nsize-1) = nsize

    dd(nsize-3:nsize,nsize) = -q(4:1:-1,1)
    imin(nsize) = nsize-3; imax(nsize) = nsize
    ir = nsize - 4
  end if
!$omp parallel do private(i)
  do i = il, ir
    dd(i-2:i-1,i) = -a(2:1:-1); dd(i+1:i+2,i) = a(1:2)
    imin(i) = i-2; imax(i) = i+2
  end do
end subroutine set_coeff_4_2

subroutine set_coeff_up_4_2 ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, dir, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(-2:2), save :: a1, a2
  CCTK_REAL, dimension(6,4), save :: q1, q2

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_4_2 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_4_2 ( a1, q1, a2, q2 )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    if ( dir == -1 ) then
      dd(1:4,1) = q1(1:4,1)
      imin(1) = 1; imax(1) = 4

      dd(1:4,2) = q1(1:4,2)
      imin(2) = 1; imax(2) = 4

      dd(1:4,3) = q1(1:4,3)
      imin(3) = 1; imax(3) = 4

      dd(1:5,4) = q1(1:5,4)
      imin(4) = 1; imax(4) = 5
    else
      dd(1:4,1) = q2(1:4,1)
      imin(1) = 1; imax(1) = 4

      dd(1:4,2) = q2(1:4,2)
      imin(2) = 1; imax(2) = 4

      dd(2:5,3) = q2(2:5,3)
      imin(3) = 2; imax(3) = 5

      dd(1:6,4) = q2(1:6,4)
      imin(4) = 1; imax(4) = 6
    end if
    il = 5
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize-5:nsize,nsize-3) = -q2(6:1:-1,4)
      imin(nsize-3) = nsize-5; imax(nsize-3) = nsize

      dd(nsize-4:nsize-1,nsize-2) = -q2(5:2:-1,3)
      imin(nsize-2) = nsize-4; imax(nsize-2) = nsize-1

      dd(nsize-3:nsize,nsize-1) = -q2(4:1:-1,2)
      imin(nsize-1) = nsize-3; imax(nsize-1) = nsize

      dd(nsize-3:nsize,nsize) = -q2(4:1:-1,1)
      imin(nsize) = nsize-3; imax(nsize) = nsize
    else
      dd(nsize-4:nsize,nsize-3) = -q1(5:1:-1,4)
      imin(nsize-3) = nsize-4; imax(nsize-3) = nsize

      dd(nsize-3:nsize,nsize-2) = -q1(4:1:-1,3)
      imin(nsize-2) = nsize-3; imax(nsize-2) = nsize

      dd(nsize-3:nsize,nsize-1) = -q1(4:1:-1,2)
      imin(nsize-1) = nsize-3; imax(nsize-1) = nsize

      dd(nsize-3:nsize,nsize) = -q1(4:1:-1,1)
      imin(nsize) = nsize-3; imax(nsize) = nsize
    end if
    ir = nsize - 4
  end if
  do i = il, ir
    if ( dir == -1 ) then
      dd(i-2:i+1,i) = a1(-2:1)
      imin(i) = i-2; imax(i) = i+1
    else
      dd(i-1:i+2,i) = a2(-1:2)
      imin(i) = i-1; imax(i) = i+2
    end if
  end do
end subroutine set_coeff_up_4_2

subroutine set_coeff_up_4 ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

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
  CCTK_REAL, dimension(4,3), save :: q1, q2

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_2 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_4 ( a1, q1, a2, q2 )
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
    else
      dd(1,1) = q2(1,1);
      imin(1) = 1; imax(1) = 1

      dd(1:5,2) = a2(-1:3);
      imin(2) = 1; imax(2) = 5

      dd(2:6,3) = a2(-1:3);
      imin(3) = 2; imax(3) = 6
    end if
    il = 4
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize,nsize) = q1(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-4:nsize,nsize-1) = a1(-3:1)
      imin(nsize-1) = nsize-4; imax(nsize-1) = nsize

      dd(nsize-5:nsize-1,nsize-2) = a1(-3:1)
      imin(nsize-2) = nsize-5; imax(nsize-2) = nsize-1
    else
      dd(nsize,nsize) = q2(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-1:nsize,nsize-1) = q2(3:4,2)
      imin(nsize-1) = nsize-1; imax(nsize-1) = nsize

      dd(nsize-2:nsize,nsize-2) = q2(2:4,3)
      imin(nsize-2) = nsize-2; imax(nsize-2) = nsize
    end if
    ir = nsize - 3
  end if
  do i = il, ir
    if ( dir == -1 ) then
      dd(i-3:i+1,i) = a1(-3:1)
      imin(i) = i-3; imax(i) = i+1
    else
      dd(i-1:i+3,i) = a2(-1:3)
      imin(i) = i-1; imax(i) = i+3
    end if
  end do
end subroutine set_coeff_up_4
