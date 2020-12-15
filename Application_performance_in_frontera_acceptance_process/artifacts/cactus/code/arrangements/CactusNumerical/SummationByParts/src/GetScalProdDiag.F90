! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
subroutine SBP_GetScalProdDiag ( cctkGH, dir, nsize, sigmad )

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_POINTER_TO_CONST, intent(IN) :: cctkGH
  CCTK_INT, intent(IN) :: dir
  CCTK_INT, intent(IN) :: nsize
  CCTK_REAL, dimension(nsize), intent(OUT) :: sigmad

  CCTK_REAL, parameter :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  integer, dimension(6) :: onesided
  CCTK_INT, dimension(6) :: offset
  CCTK_INT :: ol, or

  CCTK_REAL, dimension(2), parameter :: bmask_2 = (/ 0.5_wp, 1.0_wp /) 
  CCTK_REAL, dimension(4), parameter :: bmask_4 = (/ 17.0_wp/48.0_wp, &
                                                     59.0_wp/48.0_wp, &
                                                     43.0_wp/48.0_wp, &
                                                     49.0_wp/48.0_wp /)
  CCTK_REAL, dimension(6), parameter :: bmask_6 = (/ 13649.0_wp/43200._wp, &
                                                     12013.0_wp/8640._wp, &
                                                     2711.0_wp/4320.0_wp, &
                                                     5359.0_wp/4320.0_wp, &
                                                     7877.0_wp/8640.0_wp, &
                                                     43801.0_wp/43200.0_wp /)
  CCTK_REAL, dimension(8), parameter :: bmask_8 = (/ 1498139.0_wp/5080320.0_wp,&
                                                     1107307.0_wp/725760.0_wp, &
                                                     20761.0_wp/80640.0_wp, &
                                                     1304999.0_wp/725760.0_wp, &
                                                     299527.0_wp/725760.0_wp, &
                                                     103097.0_wp/80640.0_wp, &
                                                     670091.0_wp/725760.0_wp, &
                                                     5127739.0_wp/5080320.0_wp/)
                                                     
  if ( dir < 0 .or. dir > 2 ) then
    call CCTK_WARN(0, 'Error: dir is outside the legal range')
  end if

  call SBP_determine_onesided_stencil (cctkGH, onesided)
  
  call get_shiftout ( cctkGH, offset )

  ol = offset(dir*2+1)
  or = offset(dir*2+2)

  sigmad(1:ol) = zero
  sigmad(nsize+1-or:nsize) = zero
  sigmad(1+ol:nsize-or) = 1.0_wp

  if ( onesided(dir*2+1) == 1 ) then
    select case (order)
    case (2)
      sigmad(1+ol:2+ol) = bmask_2
    case (4)
      sigmad(1+ol:4+ol) = bmask_4
    case (6)
      sigmad(1+ol:6+ol) = bmask_6
    case (8)
      sigmad(1+ol:8+ol) = bmask_8
    end select
  end if
  if ( onesided(dir*2+2) == 1 ) then
    select case (order)
    case (2)
      sigmad(nsize-or:nsize-1-or:-1) = bmask_2
    case (4)
      sigmad(nsize-or:nsize-3-or:-1) = bmask_4
    case (6)
      sigmad(nsize-or:nsize-5-or:-1) = bmask_6
    case (8)
      sigmad(nsize-or:nsize-7-or:-1) = bmask_8
    end select
  end if

end subroutine SBP_GetScalProdDiag
