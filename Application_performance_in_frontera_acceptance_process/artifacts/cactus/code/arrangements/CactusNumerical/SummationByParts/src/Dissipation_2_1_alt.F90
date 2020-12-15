! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_2_1_alt (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

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
  CCTK_REAL, dimension(4,2) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  idel = epsilon / ( 16 * delta(1) )

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)

!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) ) * idel
!$OMP END PARALLEL WORKSHARE

    il = 3 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)

!$OMP PARALLEL WORKSHARE
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) ) * idel
!$OMP END PARALLEL WORKSHARE

    ir = or - 2
  end if

!$OMP PARALLEL WORKSHARE
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -6.0_wp * var(il:ir,:,:) + &
                      4.0_wp * ( var(il-1:ir-1,:,:) + &
                                 var(il+1:ir+1,:,:) ) - &
                               ( var(il-2:ir-2,:,:) + &
                                 var(il+2:ir+2,:,:) ) ) * idel
!$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 16 * delta(2) )
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)

!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      jl = 3 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)

!$OMP PARALLEL WORKSHARE
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 2
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -6.0_wp * var(:,jl:jr,:) + &
                        4.0_wp * ( var(:,jl-1:jr-1,:) + &
                                   var(:,jl+1:jr+1,:) ) - &
                                 ( var(:,jl-2:jr-2,:) + &
                                   var(:,jl+2:jr+2,:) ) ) * idel
!$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 16 * delta(3) )
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      kl = 3 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      kr = or - 2
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -6.0_wp * var(:,:,kl:kr) + &
                        4.0_wp * ( var(:,:,kl-1:kr-1) + &
                                   var(:,:,kl+1:kr+1) ) - &
                                 ( var(:,:,kl-2:kr-2) + &
                                   var(:,:,kl+2:kr+2) ) ) * idel
!$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(4,2), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -2.0_wp
    a(2,1) = 4.0_wp
    a(3,1) = -2.0_wp
    a(4,1) = 0.0_wp
    a(1,2) = 2.0_wp
    a(2,2) = -5.0_wp
    a(3,2) = 4.0_wp
    a(4,2) = -1.0_wp

  end subroutine set_coeff

end subroutine dissipation_2_1_alt

subroutine dissipation_2_1_delta (var, ni, nj, nk, bb, gsize, offset, &
                                  dx, dy, dz, epsilon, rhs)

  implicit none

  DECLARE_CCTK_PARAMETERS

  integer ::  ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(in) :: var
  CCTK_REAL, dimension(ni,nj,nk), intent(inout) :: rhs
  CCTK_INT, dimension(6), intent(in) :: bb
  CCTK_INT, dimension(3), intent(in) :: gsize
  CCTK_INT, dimension(6), intent(in) :: offset
  CCTK_REAL, dimension(ni,nj,nk), intent(in) :: dx, dy, dz
  CCTK_REAL, intent(in) :: epsilon 

  CCTK_REAL :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  CCTK_REAL, dimension(4,2) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  idel = epsilon / 16

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)

!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) ) * idel / dx(1+ol,:,:)
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) ) * &
                 idel / dx(2+ol,:,:)
!$OMP END PARALLEL WORKSHARE

    il = 3 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)

!$OMP PARALLEL WORKSHARE
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) ) * &
                 idel / dx(or-1,:,:)
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) ) * idel / dx(or,:,:)
!$OMP END PARALLEL WORKSHARE

    ir = or - 2
  end if

!$OMP PARALLEL WORKSHARE
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -6.0_wp * var(il:ir,:,:) + &
                      4.0_wp * ( var(il-1:ir-1,:,:) + &
                                 var(il+1:ir+1,:,:) ) - &
                               ( var(il-2:ir-2,:,:) + &
                                 var(il+2:ir+2,:,:) ) ) * &
                   idel / dx(il:ir,:,:)
!$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 16
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)

!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) ) * idel / dy(:,1+ol,:)
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) ) * &
                   idel /dy(:,2+ol,:)
!$OMP END PARALLEL WORKSHARE
  
      jl = 3 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)

!$OMP PARALLEL WORKSHARE
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) ) * &
                   idel / dy(:,or-1,:)
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) ) * idel / dy(:,or,:)
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 2
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -6.0_wp * var(:,jl:jr,:) + &
                        4.0_wp * ( var(:,jl-1:jr-1,:) + &
                                   var(:,jl+1:jr+1,:) ) - &
                                 ( var(:,jl-2:jr-2,:) + &
                                   var(:,jl+2:jr+2,:) ) ) * &
                     idel / dy(:,jl:jr,:)
!$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 16
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) ) * idel / dz (:,:,1+ol)
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) ) * &
                   idel / dz(:,:,2+ol)
!$OMP END PARALLEL WORKSHARE
  
      kl = 3 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) ) * &
                   idel / dz(:,:,or-1)
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) ) * idel / dz(:,:,or)
!$OMP END PARALLEL WORKSHARE
  
      kr = or - 2
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -6.0_wp * var(:,:,kl:kr) + &
                        4.0_wp * ( var(:,:,kl-1:kr-1) + &
                                   var(:,:,kl+1:kr+1) ) - &
                                 ( var(:,:,kl-2:kr-2) + &
                                   var(:,:,kl+2:kr+2) ) ) * &
                     idel / dz(:,:,kl:kr)
!$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(4,2), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -2.0_wp
    a(2,1) = 4.0_wp
    a(3,1) = -2.0_wp
    a(4,1) = 0.0_wp
    a(1,2) = 2.0_wp
    a(2,2) = -5.0_wp
    a(3,2) = 4.0_wp
    a(4,2) = -1.0_wp

  end subroutine set_coeff

end subroutine dissipation_2_1_delta
