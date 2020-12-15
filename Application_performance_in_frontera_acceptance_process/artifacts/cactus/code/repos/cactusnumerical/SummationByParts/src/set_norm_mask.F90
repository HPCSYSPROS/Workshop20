! $Header$

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
      
subroutine SBP_SetNormMask (CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  integer ::  onesided(6)
  CCTK_INT :: np
  CCTK_REAL, dimension(:), allocatable :: mask_x, mask_y, mask_z

! Note: The first number is twice the value from the paper, since Carpet
! Multiplies by 1/2 at the boundary when doing the sum reduction.
  CCTK_REAL, dimension(2), parameter :: bmask_2 = (/ 1.0_wp, 1.0_wp /) 
  CCTK_REAL, dimension(3), parameter :: bmask_3 = (/ 1.0_wp, 1.0_wp, 1.0_wp /) 
  CCTK_REAL, dimension(4), parameter :: bmask_4 = (/ 34.0_wp/48.0_wp, &
                                                     59.0_wp/48.0_wp, &
                                                     43.0_wp/48.0_wp, &
                                                     49.0_wp/48.0_wp /)
  CCTK_REAL, dimension(6), parameter :: bmask_6 = (/ 27298.0_wp/43200._wp, &
                                                     12013.0_wp/8640._wp, &
                                                     2711.0_wp/4320.0_wp, &
                                                     5359.0_wp/4320.0_wp, &
                                                     7877.0_wp/8640.0_wp, &
                                                     43801.0_wp/43200.0_wp /)
  CCTK_REAL, dimension(8), parameter :: bmask_8 = (/ 2996278.0_wp/5080320.0_wp,&
                                                     1107307.0_wp/725760.0_wp, &
                                                     20761.0_wp/80640.0_wp, &
                                                     1304999.0_wp/725760.0_wp, &
                                                     299527.0_wp/725760.0_wp, &
                                                     103097.0_wp/80640.0_wp, &
                                                     670091.0_wp/725760.0_wp, &
                                                     5127739.0_wp/5080320.0_wp/)
                                                     
  CCTK_REAL, dimension(8) :: bmask
  CCTK_INT :: ni, nj, nk
  CCTK_INT :: i, j, k

  ni = cctk_lsh(1); nj = cctk_lsh(2); nk = cctk_lsh(3)

  call SBP_determine_onesided_stencil (cctkGH, onesided)


  if (any (onesided/=0)) then
     allocate ( mask_x(ni), mask_y(nj), mask_z(nk) )
    mask_x = 1.0d0; mask_y = 1.0d0; mask_z = 1.0d0

    if ( CCTK_EQUALS(norm_type,'Diagonal') ) then
      select case (order)
      case (2)
        bmask(1:2) = bmask_2
      case (4)
        bmask(1:4) = bmask_4
      case(6)
        bmask(1:6) = bmask_6
      case(8)
        bmask(1:8) = bmask_8
      case default
        call CCTK_WARN (0, "Unknown stencil specified")
      end select 
      np = order
    else
      select case (order)
      case (4,6)
        bmask(1:3) = bmask_3
        np = 3
      case default
        call CCTK_WARN (0, "Unknown stencil specified")
      end select
    end if

    if (onesided(1)/=0) mask_x(1:np) = bmask(1:np)
    if (onesided(2)/=0) mask_x(ni:ni-np+1:-1) = bmask(1:np)
    if (onesided(3)/=0) mask_y(1:np) = bmask(1:np)
    if (onesided(4)/=0) mask_y(nj:nj-np+1:-1) = bmask(1:np)
    if (onesided(5)/=0) mask_z(1:np) = bmask(1:np)
    if (onesided(6)/=0) mask_z(nk:nk-np+1:-1) = bmask(1:np)

!$omp parallel workshare
    forall (i=1:ni, j=1:nj, k=1:nk)
       nmask(i,j,k) = mask_x(i) * mask_y(j) * mask_z(k)
    end forall
!$omp end parallel workshare

  else

!$omp parallel workshare
    nmask = 1.0_wp
!$omp end parallel workshare

  end if

end subroutine SBP_SetNormMask
