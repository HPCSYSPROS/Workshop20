! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_2_1 (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

  implicit none

  DECLARE_CCTK_PARAMETERS

  integer ::  ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(in) :: var
  CCTK_REAL, dimension(ni,nj,nk), intent(inout) :: rhs
  CCTK_INT, dimension(6), intent(in) :: bb
  CCTK_INT, dimension(3), intent(in) :: gsize
  CCTK_INT, dimension(6), intent(in) :: offset
  CCTK_REAL, dimension(3), intent(in) :: delta
  CCTK_REAL, intent(in) :: epsilon 

  CCTK_REAL :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  CCTK_REAL, dimension(3,2) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  if ( scale_with_h > 0 ) then
    idel = epsilon / ( 4 * delta(1) )
  else
    idel = epsilon / 4
  end if

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) ) * idel

    il = 3 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) ) * idel

    ir = or - 2
  end if
!$omp parallel workshare
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -2.0_wp * var(il:ir,:,:) + &
                             ( var(il-1:ir-1,:,:) + &
                               var(il+1:ir+1,:,:) ) ) * idel

!$omp end parallel workshare
  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 4 * delta(2) )
    else
      idel = epsilon / 4
    end if
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) ) * idel
  
      jl = 3 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) ) * idel
  
      jr = or - 2
    end if
!$omp parallel workshare
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -2.0_wp * var(:,jl:jr,:) + &
                               ( var(:,jl-1:jr-1,:) + &
                                 var(:,jl+1:jr+1,:) ) ) * idel
!$omp end parallel workshare
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 4 * delta(3) )
    else
      idel = epsilon / 4
    end if
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) ) * idel
  
      kl = 3 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) ) * idel
  
      kr = or - 2
    end if
!$omp parallel workshare
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -2.0_wp * var(:,:,kl:kr) + &
                               ( var(:,:,kl-1:kr-1) + &
                                 var(:,:,kl+1:kr+1) ) ) * idel
!$omp end parallel workshare
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(3,2), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -2.0_wp
    a(2,1) = 2.0_wp
    a(3,1) = 0.0_wp
    a(1,2) = 1.0_wp
    a(2,2) = -2.0_wp
    a(3,2) = 1.0_wp

  end subroutine set_coeff

end subroutine dissipation_2_1
