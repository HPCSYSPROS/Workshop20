#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine set_coeff2_2 ( nsize, loc_order, bb, gsize, imin, imax, dd )

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
  CCTK_REAL, dimension(1,1), save :: q 

  CCTK_INT :: i, il, ir

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_2_2 ( a, q )
    first = .false.
  end if

  dd = zero
  imin = 0
  imax = -1

  if ( bb(1) == 0 ) then
    il = 1 + gsize
  else
    dd(1,1) = q(1,1);
    imin(1) = 1; imax(1) = 1
    il = 2
  end if
  if ( bb(2) == 0 ) then
    ir = nsize - gsize
  else
    dd(nsize,nsize) = q(1,1)
    imin(nsize) = nsize; imax(nsize) = nsize
    ir = nsize - 1
  end if
  do i = il, ir
    dd(i-1,i) = a(2); dd(i:i+1,i) = a(1:2)
    imin(i) = i-1; imax(i) = i+1
  end do

end subroutine set_coeff2_2
