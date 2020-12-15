#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff_2_1 ( nsize, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(1), save :: a
  CCTK_REAL, dimension(3,2), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_2_1 ( a, q )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1,1) = q(1,1); dd(2,1) = q(2,1)
    imin(1) = 1; imax(1) = 2
    dd(1,2) = q(1,2); dd(3,2) = q(3,2)
    imin(2) = 1; imax(2) = 3
    il = 3
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize-2,nsize-1) = -q(3,2); dd(nsize,nsize-1) = -q(1,2)
    imin(nsize-1) = nsize-2; imax(nsize-1) = nsize
    dd(nsize-1,nsize) = -q(2,1); dd(nsize,nsize) = -q(1,1)
    imin(nsize) = nsize-1; imax(nsize) = nsize
    ir = nsize - 2
  end if
!$omp parallel do private(i)
  do i = il, ir
    dd(i-1,i) = -a(1); dd(i+1,i) = a(1)
    imin(i) = i-1; imax(i) = i+1
  end do

end subroutine set_coeff_2_1

subroutine set_coeff_up_2_1 ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: nsize, dir, loc_order
  CCTK_INT, dimension(2), intent(IN) :: bb
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, dimension(nsize), intent(OUT) :: imin, imax
  CCTK_REAL, dimension(nsize,nsize), intent(OUT) :: dd

  CCTK_REAL, dimension(-1:1), save :: a1, a2
  CCTK_REAL, dimension(3,2), save :: q1, q2

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_2_1 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_2_1 ( a1, q1, a2, q2 )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    if ( dir == -1 ) then
      dd(1,1) = q1(1,1); dd(2,1) = q1(2,1)
      imin(1) = 1; imax(1) = 2
      dd(1,2) = q1(1,2); dd(2,2) = q1(2,2)
      imin(2) = 1; imax(2) = 2
    else
      dd(1,1) = q2(1,1); dd(2,1) = q2(2,1)
      imin(1) = 1; imax(1) = 2
      dd(2,2) = q2(2,2); dd(3,2) = q2(3,2)
      imin(2) = 2; imax(2) = 3
    end if
    il = 3
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize-2,nsize-1) = -q2(3,2); dd(nsize-1,nsize-1) = -q2(2,2)
      imin(nsize-1) = nsize-2; imax(nsize-1) = nsize-1
      dd(nsize-1,nsize) = -q2(2,1); dd(nsize,nsize) = -q2(1,1)
      imin(nsize) = nsize-1; imax(nsize) = nsize
    else
      dd(nsize-1,nsize-1) = -q1(2,2); dd(nsize,nsize-1) = -q1(1,2)
      imin(nsize-1) = nsize-1; imax(nsize-1) = nsize
      dd(nsize-1,nsize) = -q1(2,1); dd(nsize,nsize) = -q1(1,1)
      imin(nsize) = nsize-1; imax(nsize) = nsize
    end if
    ir = nsize - 2
  end if
  do i = il, ir
    if ( dir == -1 ) then
      dd(i-1:i,i) = a1(-1:0)
      imin(i) = i-1; imax(i) = i
    else
      dd(i:i+1,i) = a2(0:1)
      imin(i) = i; imax(i) = i+1
    end if
  end do
end subroutine set_coeff_up_2_1

subroutine set_coeff_up_2 ( nsize, dir, loc_order, bb, gsize, imin, imax, dd )

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
  CCTK_REAL, dimension(3,2), save :: q1, q2

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( .not. ( dir == -1 .or. dir == 1 ) ) then
    call CCTK_WARN(0, 'Internal error: set_coeff_up_2 called with invalid dir value')
  end if

  if ( first ) then
    call coeffs_up_2 ( a1, q1, a2, q2 )
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
    else
      dd(1,1) = q2(1,1);
      imin(1) = 1; imax(1) = 1

      dd(2:4,2) = a2(0:2);
      imin(2) = 2; imax(2) = 4
    end if
    il = 3
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    if ( dir == -1 ) then
      dd(nsize,nsize) = q1(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-3:nsize-1,nsize-1) = a1(-2:0)
      imin(nsize-1) = nsize-3; imax(nsize-1) = nsize-1
    else
      dd(nsize,nsize) = q2(1,1)
      imin(nsize) = nsize; imax(nsize) = nsize

      dd(nsize-1:nsize,nsize-1) = q2(2:3,2)
      imin(nsize-1) = nsize-1; imax(nsize-1) = nsize
    end if
    ir = nsize - 2
  end if
  do i = il, ir
    if ( dir == -1 ) then
      dd(i-2:i,i) = a1(-2:0)
      imin(i) = i-2; imax(i) = i
    else
      dd(i:i+2,i) = a2(0:2)
      imin(i) = i; imax(i) = i+2
    end if
  end do
end subroutine set_coeff_up_2
