! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_8_4_alt (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

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
  CCTK_REAL, dimension(13,8) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or, i, j, k

  call set_coeff ( a )

  idel = epsilon / ( 1024 * delta(1) )

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)

!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) + &
                   a(5,1) * var(5+ol,:,:) + a(6,1) * var(6+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) + a(6,2) * var(6+ol,:,:) + &
                   a(7,2) * var(7+ol,:,:) ) * idel
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) + &
                   a(7,3) * var(7+ol,:,:) + a(8,3) * var(8+ol,:,:) ) * idel
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) + a(8,4) * var(8+ol,:,:) + &
                   a(9,4) * var(9+ol,:,:) ) * idel
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(1,5) * var(1+ol,:,:) + a(2,5) * var(2+ol,:,:) + &
                   a(3,5) * var(3+ol,:,:) + a(4,5) * var(4+ol,:,:) + &
                   a(5,5) * var(5+ol,:,:) + a(6,5) * var(6+ol,:,:) + &
                   a(7,5) * var(7+ol,:,:) + a(8,5) * var(8+ol,:,:) + &
                   a(9,5) * var(9+ol,:,:) + a(10,5) * var(10+ol,:,:) ) * idel
    rhs(6+ol,:,:) = rhs(6+ol,:,:) + &
                 ( a(1,6) * var(1+ol,:,:) + a(2,6) * var(2+ol,:,:) + &
                   a(3,6) * var(3+ol,:,:) + a(4,6) * var(4+ol,:,:) + &
                   a(5,6) * var(5+ol,:,:) + a(6,6) * var(6+ol,:,:) + &
                   a(7,6) * var(7+ol,:,:) + a(8,6) * var(8+ol,:,:) + &
                   a(9,6) * var(9+ol,:,:) + a(10,6) * var(10+ol,:,:) + &
                   a(11,6) * var(11+ol,:,:) ) * idel
    rhs(7+ol,:,:) = rhs(7+ol,:,:) + &
                 ( a(2,7) * var(2+ol,:,:) + a(3,7) * var(3+ol,:,:) + &
                   a(4,7) * var(4+ol,:,:) + a(5,7) * var(5+ol,:,:) + &
                   a(6,7) * var(6+ol,:,:) + a(7,7) * var(7+ol,:,:) + &
                   a(8,7) * var(8+ol,:,:) + a(9,7) * var(9+ol,:,:) + &
                   a(10,7) * var(10+ol,:,:) + a(11,7) * var(11+ol,:,:) + &
                   a(12,7) * var(12+ol,:,:) ) * idel
    rhs(8+ol,:,:) = rhs(8+ol,:,:) + &
                 ( a(3,8) * var(3+ol,:,:) + a(4,8) * var(4+ol,:,:) + &
                   a(5,8) * var(5+ol,:,:) + a(6,8) * var(6+ol,:,:) + &
                   a(7,8) * var(7+ol,:,:) + a(8,8) * var(8+ol,:,:) + &
                   a(9,8) * var(9+ol,:,:) + a(10,8) * var(10+ol,:,:) + &
                   a(11,8) * var(11+ol,:,:) + a(12,8) * var(12+ol,:,:) + &
                   a(13,8) * var(13+ol,:,:) ) * idel
!$OMP END PARALLEL WORKSHARE

    il = 9 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)

!$OMP PARALLEL WORKSHARE
    rhs(or-7,:,:) = rhs(or-7,:,:) + &
                 ( a(3,8) * var(or-2,:,:) + a(4,8) * var(or-3,:,:) + &
                   a(5,8) * var(or-4,:,:) + a(6,8) * var(or-5,:,:) + &
                   a(7,8) * var(or-6,:,:) + a(8,8) * var(or-7,:,:) + &
                   a(9,8) * var(or-8,:,:) + a(10,8) * var(or-9,:,:) + &
                   a(11,8) * var(or-10,:,:) + a(12,8) * var(or-11,:,:) + &
                   a(13,8) * var(or-12,:,:) ) * idel
    rhs(or-6,:,:) = rhs(or-6,:,:) + &
                 ( a(2,7) * var(or-1,:,:) + a(3,7) * var(or-2,:,:) + &
                   a(4,7) * var(or-3,:,:) + a(5,7) * var(or-4,:,:) + &
                   a(6,7) * var(or-5,:,:) + a(7,7) * var(or-6,:,:) + &
                   a(8,7) * var(or-7,:,:) + a(9,7) * var(or-8,:,:) + &
                   a(10,7) * var(or-9,:,:) + a(11,7) * var(or-10,:,:) + &
                   a(12,7) * var(or-11,:,:) ) * idel
    rhs(or-5,:,:) = rhs(or-5,:,:) + &
                 ( a(1,6) * var(or,:,:) + a(2,6) * var(or-1,:,:) + &
                   a(3,6) * var(or-2,:,:) + a(4,6) * var(or-3,:,:) + &
                   a(5,6) * var(or-4,:,:) + a(6,6) * var(or-5,:,:) + &
                   a(7,6) * var(or-6,:,:) + a(8,6) * var(or-7,:,:) + &
                   a(9,6) * var(or-8,:,:) + a(10,6) * var(or-9,:,:) + &
                   a(11,6) * var(or-10,:,:) ) * idel
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                 ( a(1,5) * var(or,:,:) + a(2,5) * var(or-1,:,:) + &
                   a(3,5) * var(or-2,:,:) + a(4,5) * var(or-3,:,:) + &
                   a(5,5) * var(or-4,:,:) + a(6,5) * var(or-5,:,:) + &
                   a(7,5) * var(or-6,:,:) + a(8,5) * var(or-7,:,:) + &
                   a(9,5) * var(or-8,:,:) + a(10,5) * var(or-9,:,:) ) * idel
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                 ( a(1,4) * var(or,:,:) + a(2,4) * var(or-1,:,:) + &
                   a(3,4) * var(or-2,:,:) + a(4,4) * var(or-3,:,:) + &
                   a(5,4) * var(or-4,:,:) + a(6,4) * var(or-5,:,:) + &
                   a(7,4) * var(or-6,:,:) + a(8,4) * var(or-7,:,:) + &
                   a(9,4) * var(or-8,:,:) ) * idel
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                 ( a(1,3) * var(or,:,:) + a(2,3) * var(or-1,:,:) + &
                   a(3,3) * var(or-2,:,:) + a(4,3) * var(or-3,:,:) + &
                   a(5,3) * var(or-4,:,:) + a(6,3) * var(or-5,:,:) + &
                   a(7,3) * var(or-6,:,:) + a(8,3) * var(or-7,:,:) ) * idel
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) + &
                   a(5,2) * var(or-4,:,:) + a(6,2) * var(or-5,:,:) + &
                   a(7,2) * var(or-6,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) + a(4,1) * var(or-3,:,:) + &
                   a(5,1) * var(or-4,:,:) + a(6,1) * var(or-5,:,:) ) * idel
!$OMP END PARALLEL WORKSHARE

    ir = or - 8
  end if

!! !$OMP PARALLEL WORKSHARE
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) SHARED(var,rhs,il,ir,nj,nk,idel)
  do k = 1, nk
     do j = 1, nj
       do i = il, ir
          rhs(i,j,k) = rhs(i,j,k) + &
                   ( -252.0_wp * var(i,j,k) + &
                      210.0_wp * ( var(i-1,j,k) + &
                                   var(i+1,j,k) ) - &
                      120.0_wp * ( var(i-2,j,k) + &
                                   var(i+2,j,k) ) + &
                       45.0_wp * ( var(i-3,j,k) + &
                                   var(i+3,j,k) ) - &
                       10.0_wp * ( var(i-4,j,k) + &
                                   var(i+4,j,k) ) + &
                                 ( var(i-5,j,k) + &
                                   var(i+5,j,k) ) ) * idel
       end do
     end do
   end do
!$OMP END PARALLEL DO
!! !$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 1024 * delta(2) )
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3);

!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) + &
                     a(5,1) * var(:,5+ol,:) + a(6,1) * var(:,6+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) + a(6,2) * var(:,6+ol,:) + &
                     a(7,2) * var(:,7+ol,:) ) * idel
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) + &
                     a(7,3) * var(:,7+ol,:) + a(8,3) * var(:,8+ol,:) ) * idel
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) + a(8,4) * var(:,8+ol,:) + &
                     a(9,4) * var(:,9+ol,:) ) * idel
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(1,5) * var(:,1+ol,:) + a(2,5) * var(:,2+ol,:) + &
                     a(3,5) * var(:,3+ol,:) + a(4,5) * var(:,4+ol,:) + &
                     a(5,5) * var(:,5+ol,:) + a(6,5) * var(:,6+ol,:) + &
                     a(7,5) * var(:,7+ol,:) + a(8,5) * var(:,8+ol,:) + &
                     a(9,5) * var(:,9+ol,:) + a(10,5) * var(:,10+ol,:) ) * idel
      rhs(:,6+ol,:) = rhs(:,6+ol,:) + &
                   ( a(1,6) * var(:,1+ol,:) + a(2,6) * var(:,2+ol,:) + &
                     a(3,6) * var(:,3+ol,:) + a(4,6) * var(:,4+ol,:) + &
                     a(5,6) * var(:,5+ol,:) + a(6,6) * var(:,6+ol,:) + &
                     a(7,6) * var(:,7+ol,:) + a(8,6) * var(:,8+ol,:) + &
                     a(9,6) * var(:,9+ol,:) + a(10,6) * var(:,10+ol,:) + &
                     a(11,6) * var(:,11+ol,:) ) * idel
      rhs(:,7+ol,:) = rhs(:,7+ol,:) + &
                   ( a(2,7) * var(:,2+ol,:) + a(3,7) * var(:,3+ol,:) + &
                     a(4,7) * var(:,4+ol,:) + a(5,7) * var(:,5+ol,:) + &
                     a(6,7) * var(:,6+ol,:) + a(7,7) * var(:,7+ol,:) + &
                     a(8,7) * var(:,8+ol,:) + a(9,7) * var(:,9+ol,:) + &
                     a(10,7) * var(:,10+ol,:) + a(11,7) * var(:,11+ol,:) + &
                     a(12,7) * var(:,12+ol,:) ) * idel
      rhs(:,8+ol,:) = rhs(:,8+ol,:) + &
                   ( a(3,8) * var(:,3+ol,:) + a(4,8) * var(:,4+ol,:) + &
                     a(5,8) * var(:,5+ol,:) + a(6,8) * var(:,6+ol,:) + &
                     a(7,8) * var(:,7+ol,:) + a(8,8) * var(:,8+ol,:) + &
                     a(9,8) * var(:,9+ol,:) + a(10,8) * var(:,10+ol,:) + &
                     a(11,8) * var(:,11+ol,:) + a(12,8) * var(:,12+ol,:) + &
                     a(13,8) * var(:,13+ol,:) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      jl = 9 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)

!$OMP PARALLEL WORKSHARE
      rhs(:,or-7,:) = rhs(:,or-7,:) + &
                   ( a(3,8) * var(:,or-2,:) + a(4,8) * var(:,or-3,:) + &
                     a(5,8) * var(:,or-4,:) + a(6,8) * var(:,or-5,:) + &
                     a(7,8) * var(:,or-6,:) + a(8,8) * var(:,or-7,:) + &
                     a(9,8) * var(:,or-8,:) + a(10,8) * var(:,or-9,:) + &
                     a(11,8) * var(:,or-10,:) + a(12,8) * var(:,or-11,:) + &
                     a(13,8) * var(:,or-12,:) ) * idel
      rhs(:,or-6,:) = rhs(:,or-6,:) + &
                   ( a(2,7) * var(:,or-1,:) + a(3,7) * var(:,or-2,:) + &
                     a(4,7) * var(:,or-3,:) + a(5,7) * var(:,or-4,:) + &
                     a(6,7) * var(:,or-5,:) + a(7,7) * var(:,or-6,:) + &
                     a(8,7) * var(:,or-7,:) + a(9,7) * var(:,or-8,:) + &
                     a(10,7) * var(:,or-9,:) + a(11,7) * var(:,or-10,:) + &
                     a(12,7) * var(:,or-11,:) ) * idel
      rhs(:,or-5,:) = rhs(:,or-5,:) + &
                   ( a(1,6) * var(:,or,:) + a(2,6) * var(:,or-1,:) + &
                     a(3,6) * var(:,or-2,:) + a(4,6) * var(:,or-3,:) + &
                     a(5,6) * var(:,or-4,:) + a(6,6) * var(:,or-5,:) + &
                     a(7,6) * var(:,or-6,:) + a(8,6) * var(:,or-7,:) + &
                     a(9,6) * var(:,or-8,:) + a(10,6) * var(:,or-9,:) + &
                     a(11,6) * var(:,or-10,:) ) * idel
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                   ( a(1,5) * var(:,or,:) + a(2,5) * var(:,or-1,:) + &
                     a(3,5) * var(:,or-2,:) + a(4,5) * var(:,or-3,:) + &
                     a(5,5) * var(:,or-4,:) + a(6,5) * var(:,or-5,:) + &
                     a(7,5) * var(:,or-6,:) + a(8,5) * var(:,or-7,:) + &
                     a(9,5) * var(:,or-8,:) + a(10,5) * var(:,or-9,:) ) * idel
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                   ( a(1,4) * var(:,or,:) + a(2,4) * var(:,or-1,:) + &
                     a(3,4) * var(:,or-2,:) + a(4,4) * var(:,or-3,:) + &
                     a(5,4) * var(:,or-4,:) + a(6,4) * var(:,or-5,:) + &
                     a(7,4) * var(:,or-6,:) + a(8,4) * var(:,or-7,:) + &
                     a(9,4) * var(:,or-8,:) ) * idel
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                   ( a(1,3) * var(:,or,:) + a(2,3) * var(:,or-1,:) + &
                     a(3,3) * var(:,or-2,:) + a(4,3) * var(:,or-3,:) + &
                     a(5,3) * var(:,or-4,:) + a(6,3) * var(:,or-5,:) + &
                     a(7,3) * var(:,or-6,:) + a(8,3) * var(:,or-7,:) ) * idel
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) + &
                     a(5,2) * var(:,or-4,:) + a(6,2) * var(:,or-5,:) + &
                     a(7,2) * var(:,or-6,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) + a(4,1) * var(:,or-3,:) + &
                     a(5,1) * var(:,or-4,:) + a(6,1) * var(:,or-5,:) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 8
    end if

!! !$OMP PARALLEL WORKSHARE
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) SHARED(var,rhs,ni,jl,jr,nk,idel)
  do k = 1, nk
     do j = jl, jr
       do i = 1, ni
          rhs(i,j,k) = rhs(i,j,k) + &
                     ( -252.0_wp * var(i,j,k) + &
                        210.0_wp * ( var(i,j-1,k) + &
                                     var(i,j+1,k) ) - &
                        120.0_wp * ( var(i,j-2,k) + &
                                     var(i,j+2,k) ) + &
                         45.0_wp * ( var(i,j-3,k) + &
                                     var(i,j+3,k) ) - &
                         10.0_wp * ( var(i,j-4,k) + &
                                     var(i,j+4,k) ) + &
                                   ( var(i,j-5,k) + &
                                     var(i,j+5,k) ) ) * idel
       end do
     end do
   end do
!$OMP END PARALLEL DO
!! !$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 1024 * delta(3) )
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5);

!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) + &
                     a(5,1) * var(:,:,5+ol) + a(6,1) * var(:,:,6+ol ) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) + a(6,2) * var(:,:,6+ol) + &
                     a(7,2) * var(:,:,7+ol ) ) * idel
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) + &
                     a(7,3) * var(:,:,7+ol) + a(8,3) * var(:,:,8+ol ) ) * idel
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) + a(8,4) * var(:,:,8+ol) + &
                     a(9,4) * var(:,:,9+ol ) ) * idel
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(1,5) * var(:,:,1+ol) + a(2,5) * var(:,:,2+ol) + &
                     a(3,5) * var(:,:,3+ol) + a(4,5) * var(:,:,4+ol) + &
                     a(5,5) * var(:,:,5+ol) + a(6,5) * var(:,:,6+ol) + &
                     a(7,5) * var(:,:,7+ol) + a(8,5) * var(:,:,8+ol) + &
                     a(9,5) * var(:,:,9+ol) + a(10,5) * var(:,:,10+ol ) ) * idel
      rhs(:,:,6+ol) = rhs(:,:,6+ol) + &
                   ( a(1,6) * var(:,:,1+ol) + a(2,6) * var(:,:,2+ol) + &
                     a(3,6) * var(:,:,3+ol) + a(4,6) * var(:,:,4+ol) + &
                     a(5,6) * var(:,:,5+ol) + a(6,6) * var(:,:,6+ol) + &
                     a(7,6) * var(:,:,7+ol) + a(8,6) * var(:,:,8+ol) + &
                     a(9,6) * var(:,:,9+ol) + a(10,6) * var(:,:,10+ol) + &
                     a(11,6) * var(:,:,11+ol ) ) * idel
      rhs(:,:,7+ol) = rhs(:,:,7+ol) + &
                   ( a(2,7) * var(:,:,2+ol) + a(3,7) * var(:,:,3+ol) + &
                     a(4,7) * var(:,:,4+ol) + a(5,7) * var(:,:,5+ol) + &
                     a(6,7) * var(:,:,6+ol) + a(7,7) * var(:,:,7+ol) + &
                     a(8,7) * var(:,:,8+ol) + a(9,7) * var(:,:,9+ol) + &
                     a(10,7) * var(:,:,10+ol) + a(11,7) * var(:,:,11+ol) + &
                     a(12,7) * var(:,:,12+ol ) ) * idel
      rhs(:,:,8+ol) = rhs(:,:,8+ol) + &
                   ( a(3,8) * var(:,:,3+ol) + a(4,8) * var(:,:,4+ol) + &
                     a(5,8) * var(:,:,5+ol) + a(6,8) * var(:,:,6+ol) + &
                     a(7,8) * var(:,:,7+ol) + a(8,8) * var(:,:,8+ol) + &
                     a(9,8) * var(:,:,9+ol) + a(10,8) * var(:,:,10+ol) + &
                     a(11,8) * var(:,:,11+ol) + a(12,8) * var(:,:,12+ol) + &
                     a(13,8) * var(:,:,13+ol ) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      kl = 9 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-7) = rhs(:,:,or-7) + &
                   ( a(3,8) * var(:,:,or-2) + a(4,8) * var(:,:,or-3) + &
                     a(5,8) * var(:,:,or-4) + a(6,8) * var(:,:,or-5) + &
                     a(7,8) * var(:,:,or-6) + a(8,8) * var(:,:,or-7) + &
                     a(9,8) * var(:,:,or-8) + a(10,8) * var(:,:,or-9) + &
                     a(11,8) * var(:,:,or-10) + a(12,8) * var(:,:,or-11) + &
                     a(13,8) * var(:,:,or-12) ) * idel
      rhs(:,:,or-6) = rhs(:,:,or-6) + &
                   ( a(2,7) * var(:,:,or-1) + a(3,7) * var(:,:,or-2) + &
                     a(4,7) * var(:,:,or-3) + a(5,7) * var(:,:,or-4) + &
                     a(6,7) * var(:,:,or-5) + a(7,7) * var(:,:,or-6) + &
                     a(8,7) * var(:,:,or-7) + a(9,7) * var(:,:,or-8) + &
                     a(10,7) * var(:,:,or-9) + a(11,7) * var(:,:,or-10) + &
                     a(12,7) * var(:,:,or-11) ) * idel
      rhs(:,:,or-5) = rhs(:,:,or-5) + &
                   ( a(1,6) * var(:,:,or) + a(2,6) * var(:,:,or-1) + &
                     a(3,6) * var(:,:,or-2) + a(4,6) * var(:,:,or-3) + &
                     a(5,6) * var(:,:,or-4) + a(6,6) * var(:,:,or-5) + &
                     a(7,6) * var(:,:,or-6) + a(8,6) * var(:,:,or-7) + &
                     a(9,6) * var(:,:,or-8) + a(10,6) * var(:,:,or-9) + &
                     a(11,6) * var(:,:,or-10) ) * idel
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                   ( a(1,5) * var(:,:,or) + a(2,5) * var(:,:,or-1) + &
                     a(3,5) * var(:,:,or-2) + a(4,5) * var(:,:,or-3) + &
                     a(5,5) * var(:,:,or-4) + a(6,5) * var(:,:,or-5) + &
                     a(7,5) * var(:,:,or-6) + a(8,5) * var(:,:,or-7) + &
                     a(9,5) * var(:,:,or-8) + a(10,5) * var(:,:,or-9) ) * idel
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                   ( a(1,4) * var(:,:,or) + a(2,4) * var(:,:,or-1) + &
                     a(3,4) * var(:,:,or-2) + a(4,4) * var(:,:,or-3) + &
                     a(5,4) * var(:,:,or-4) + a(6,4) * var(:,:,or-5) + &
                     a(7,4) * var(:,:,or-6) + a(8,4) * var(:,:,or-7) + &
                     a(9,4) * var(:,:,or-8) ) * idel
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                   ( a(1,3) * var(:,:,or) + a(2,3) * var(:,:,or-1) + &
                     a(3,3) * var(:,:,or-2) + a(4,3) * var(:,:,or-3) + &
                     a(5,3) * var(:,:,or-4) + a(6,3) * var(:,:,or-5) + &
                     a(7,3) * var(:,:,or-6) + a(8,3) * var(:,:,or-7) ) * idel
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) + &
                     a(5,2) * var(:,:,or-4) + a(6,2) * var(:,:,or-5) + &
                     a(7,2) * var(:,:,or-6) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) + a(4,1) * var(:,:,or-3) + &
                     a(5,1) * var(:,:,or-4) + a(6,1) * var(:,:,or-5) ) * idel
!$OMP END PARALLEL WORKSHARE

      kr = or - 8
    end if

!! !$OMP PARALLEL WORKSHARE
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) SHARED(var,rhs,ni,nj,kl,kr,idel)
  do k = kl, kr
     do j = 1, nj
       do i = 1, ni
          rhs(i,j,k) = rhs(i,j,k) + &
                     ( -252.0_wp * var(i,j,k) + &
                        210.0_wp * ( var(i,j,k-1) + &
                                     var(i,j,k+1) ) - &
                        120.0_wp * ( var(i,j,k-2) + &
                                     var(i,j,k+2) ) + &
                         45.0_wp * ( var(i,j,k-3) + &
                                     var(i,j,k+3) ) - &
                         10.0_wp * ( var(i,j,k-4) + &
                                     var(i,j,k+4) ) + &
                                   ( var(i,j,k-5) + &
                                     var(i,j,k+5) ) ) * idel
       end do
     end do
   end do
!$OMP END PARALLEL DO
!! !$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(13,8), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -3.3910872088637970174997113084967416241083103770745_wp
    a(2,1) = 16.955436044318985087498556542483708120541551885372_wp
    a(3,1) = -33.910872088637970174997113084967416241083103770745_wp
    a(4,1) = 33.910872088637970174997113084967416241083103770745_wp
    a(5,1) = -16.955436044318985087498556542483708120541551885372_wp
    a(6,1) = 3.3910872088637970174997113084967416241083103770745_wp
    a(7,1) = 0.0_wp
    a(8,1) = 0.0_wp
    a(9,1) = 0.0_wp
    a(10,1) = 0.0_wp
    a(11,1) = 0.0_wp
    a(12,1) = 0.0_wp
    a(13,1) = 0.0_wp
    a(1,2) = 3.2771399440263630592058029074141137010783820566473_wp
    a(2,2) = -17.041127708937087907870175118553391245607586694566_wp
    a(3,2) = 36.048539384289993651263831981555250711862202623121_wp
    a(4,2) = -39.325679328316356710469634888969364412940584679768_wp
    a(5,2) = 22.939979608184541414440620351898795907548674396531_wp
    a(6,2) = -6.5542798880527261184116058148282274021567641132947_wp
    a(7,2) = 0.65542798880527261184116058148282274021567641132947_wp
    a(8,2) = 0.0_wp
    a(9,2) = 0.0_wp
    a(10,2) = 0.0_wp
    a(11,2) = 0.0_wp
    a(12,2) = 0.0_wp
    a(13,2) = 0.0_wp
    a(1,3) = -38.842059631038967294446317614758441308222147295410_wp
    a(2,3) = 213.63132797071432011945474688117142719522181012475_wp
    a(3,3) = -489.40995135109098791002360194595636048359905592216_wp
    a(4,3) = 602.05192428110399306391792302875584027744328307885_wp
    a(5,3) = -427.26265594142864023890949376234285439044362024951_wp
    a(6,3) = 174.78926833967535282500842926641298588699966282934_wp
    a(7,3) = -38.842059631038967294446317614758441308222147295410_wp
    a(8,3) = 3.8842059631038967294446317614758441308222147295410_wp
    a(9,3) = 0.0_wp
    a(10,3) = 0.0_wp
    a(11,3) = 0.0_wp
    a(12,3) = 0.0_wp
    a(13,3) = 0.0_wp
    a(1,4) = 5.5613835719414344378807953109542612676331552744485_wp
    a(2,4) = -33.368301431648606627284771865725567605798931646691_wp
    a(3,4) = 86.201445365092233787152327319791049648313906753952_wp
    a(4,4) = -125.68726872587641829610597402756630464850930920254_wp
    a(5,4) = 114.00836322479940597655630387456235598647968312619_wp
    a(6,4) = -66.736602863297213254569543731451135211597863293382_wp
    a(7,4) = 25.026226073736454970463578899294175704349198735018_wp
    a(8,4) = -5.5613835719414344378807953109542612676331552744485_wp
    a(9,4) = 0.55613835719414344378807953109542612676331552744485_wp
    a(10,4) = 0.0_wp
    a(11,4) = 0.0_wp
    a(12,4) = 0.0_wp
    a(13,4) = 0.0_wp
    a(1,5) = -12.115101476661536355654081268132755978592914829047_wp
    a(2,5) = 84.805710336630754489578568876929291850150403803330_wp
    a(3,5) = -266.53223248655379982438978789892063152904412623904_wp
    a(4,5) = 496.71916054312299058181733199344299512230950799093_wp
    a(5,5) = -608.17809412840912505383487966026435012536432441817_wp
    a(6,5) = 508.83426201978452693747141326157575110090242281998_wp
    a(7,5) = -290.76243543987687253569795043518614348622995589713_wp
    a(8,5) = 109.03591328995382720088673141319480380733623346142_wp
    a(9,5) = -24.230202953323072711308162536265511957185829658094_wp
    a(10,5) = 2.4230202953323072711308162536265511957185829658094_wp
    a(11,5) = 0.0_wp
    a(12,5) = 0.0_wp
    a(13,5) = 0.0_wp
    a(1,6) = 0.78217600900123184961735065035840034142603567514089_wp
    a(2,6) = -7.8217600900123184961735065035840034142603567514089_wp
    a(3,6) = 35.197920405055433232780779266128015364171605381340_wp
    a(4,6) = -93.861121080147821954082078043008040971124281016906_wp
    a(5,6) = 164.25696189025868841964363657526407169946749177959_wp
    a(6,6) = -197.10835426831042610357236389031688603936099013550_wp
    a(7,6) = 164.25696189025868841964363657526407169946749177959_wp
    a(8,6) = -93.861121080147821954082078043008040971124281016906_wp
    a(9,6) = 35.197920405055433232780779266128015364171605381340_wp
    a(10,6) = -7.8217600900123184961735065035840034142603567514089_wp
    a(11,6) = 0.78217600900123184961735065035840034142603567514089_wp
    a(12,6) = 0.0_wp
    a(13,6) = 0.0_wp
    a(1,7) = 0.0_wp
    a(2,7) = 1.0830767761393601764536458481012280421614377748694_wp
    a(3,7) = -10.830767761393601764536458481012280421614377748694_wp
    a(4,7) = 48.738454926271207940414063164555261897264699869122_wp
    a(5,7) = -129.96921313672322117443750177214736505937253298433_wp
    a(6,7) = 227.44612298926563705526562810125788885390193272257_wp
    a(7,7) = -272.93534758711876446631875372150946662468231926708_wp
    a(8,7) = 227.44612298926563705526562810125788885390193272257_wp
    a(9,7) = -129.96921313672322117443750177214736505937253298433_wp
    a(10,7) = 48.738454926271207940414063164555261897264699869122_wp
    a(11,7) = -10.830767761393601764536458481012280421614377748694_wp
    a(12,7) = 1.0830767761393601764536458481012280421614377748694_wp
    a(13,7) = 0.0_wp
    a(1,8) = 0.0_wp
    a(2,8) = 0.0_wp
    a(3,8) = 0.99075245444434671889501396229410272246695863420506_wp
    a(4,8) = -9.9075245444434671889501396229410272246695863420506_wp
    a(5,8) = 44.583860449995602350275628303234622511013138539228_wp
    a(6,8) = -118.89029453332160626740167547529232669603503610461_wp
    a(7,8) = 208.05801543331281096795293208176157171806131318306_wp
    a(8,8) = -249.66961851997537316154351849811388606167357581967_wp
    a(9,8) = 208.05801543331281096795293208176157171806131318306_wp
    a(10,8) = -118.89029453332160626740167547529232669603503610461_wp
    a(11,8) = 44.583860449995602350275628303234622511013138539228_wp
    a(12,8) = -9.9075245444434671889501396229410272246695863420506_wp
    a(13,8) = 0.99075245444434671889501396229410272246695863420506_wp

  end subroutine set_coeff

end subroutine dissipation_8_4_alt

subroutine dissipation_8_4_delta (var, ni, nj, nk, bb, gsize, offset, &
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
  CCTK_REAL, dimension(13,8) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  idel = epsilon / 1024

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)

!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) + &
                   a(5,1) * var(5+ol,:,:) + a(6,1) * var(6+ol,:,:) ) * &
                 idel / dx(1+ol,:,:)
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) + a(6,2) * var(6+ol,:,:) + &
                   a(7,2) * var(7+ol,:,:) ) * idel / dx(2+ol,:,:)
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) + &
                   a(7,3) * var(7+ol,:,:) + a(8,3) * var(8+ol,:,:) ) * &
                 idel / dx(3+ol,:,:)
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) + a(8,4) * var(8+ol,:,:) + &
                   a(9,4) * var(9+ol,:,:) ) * idel / dx(4+ol,:,:)
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(1,5) * var(1+ol,:,:) + a(2,5) * var(2+ol,:,:) + &
                   a(3,5) * var(3+ol,:,:) + a(4,5) * var(4+ol,:,:) + &
                   a(5,5) * var(5+ol,:,:) + a(6,5) * var(6+ol,:,:) + &
                   a(7,5) * var(7+ol,:,:) + a(8,5) * var(8+ol,:,:) + &
                   a(9,5) * var(9+ol,:,:) + a(10,5) * var(10+ol,:,:) ) * &
                 idel / dx(5+ol,:,:)
    rhs(6+ol,:,:) = rhs(6+ol,:,:) + &
                 ( a(1,6) * var(1+ol,:,:) + a(2,6) * var(2+ol,:,:) + &
                   a(3,6) * var(3+ol,:,:) + a(4,6) * var(4+ol,:,:) + &
                   a(5,6) * var(5+ol,:,:) + a(6,6) * var(6+ol,:,:) + &
                   a(7,6) * var(7+ol,:,:) + a(8,6) * var(8+ol,:,:) + &
                   a(9,6) * var(9+ol,:,:) + a(10,6) * var(10+ol,:,:) + &
                   a(11,6) * var(11+ol,:,:) ) * idel / dx(6+ol,:,:)
    rhs(7+ol,:,:) = rhs(7+ol,:,:) + &
                 ( a(2,7) * var(2+ol,:,:) + a(3,7) * var(3+ol,:,:) + &
                   a(4,7) * var(4+ol,:,:) + a(5,7) * var(5+ol,:,:) + &
                   a(6,7) * var(6+ol,:,:) + a(7,7) * var(7+ol,:,:) + &
                   a(8,7) * var(8+ol,:,:) + a(9,7) * var(9+ol,:,:) + &
                   a(10,7) * var(10+ol,:,:) + a(11,7) * var(11+ol,:,:) + &
                   a(12,7) * var(12+ol,:,:) ) * idel / dx(7+ol,:,:)
    rhs(8+ol,:,:) = rhs(8+ol,:,:) + &
                 ( a(3,8) * var(3+ol,:,:) + a(4,8) * var(4+ol,:,:) + &
                   a(5,8) * var(5+ol,:,:) + a(6,8) * var(6+ol,:,:) + &
                   a(7,8) * var(7+ol,:,:) + a(8,8) * var(8+ol,:,:) + &
                   a(9,8) * var(9+ol,:,:) + a(10,8) * var(10+ol,:,:) + &
                   a(11,8) * var(11+ol,:,:) + a(12,8) * var(12+ol,:,:) + &
                   a(13,8) * var(13+ol,:,:) ) * idel / dx(8+ol,:,:)
!$OMP END PARALLEL WORKSHARE

    il = 9 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)

!$OMP PARALLEL WORKSHARE
    rhs(or-7,:,:) = rhs(or-7,:,:) + &
                 ( a(3,8) * var(or-2,:,:) + a(4,8) * var(or-3,:,:) + &
                   a(5,8) * var(or-4,:,:) + a(6,8) * var(or-5,:,:) + &
                   a(7,8) * var(or-6,:,:) + a(8,8) * var(or-7,:,:) + &
                   a(9,8) * var(or-8,:,:) + a(10,8) * var(or-9,:,:) + &
                   a(11,8) * var(or-10,:,:) + a(12,8) * var(or-11,:,:) + &
                   a(13,8) * var(or-12,:,:) ) * idel / dx(or-7,:,:)
    rhs(or-6,:,:) = rhs(or-6,:,:) + &
                 ( a(2,7) * var(or-1,:,:) + a(3,7) * var(or-2,:,:) + &
                   a(4,7) * var(or-3,:,:) + a(5,7) * var(or-4,:,:) + &
                   a(6,7) * var(or-5,:,:) + a(7,7) * var(or-6,:,:) + &
                   a(8,7) * var(or-7,:,:) + a(9,7) * var(or-8,:,:) + &
                   a(10,7) * var(or-9,:,:) + a(11,7) * var(or-10,:,:) + &
                   a(12,7) * var(or-11,:,:) ) * idel / dx(or-6,:,:)
    rhs(or-5,:,:) = rhs(or-5,:,:) + &
                 ( a(1,6) * var(or,:,:) + a(2,6) * var(or-1,:,:) + &
                   a(3,6) * var(or-2,:,:) + a(4,6) * var(or-3,:,:) + &
                   a(5,6) * var(or-4,:,:) + a(6,6) * var(or-5,:,:) + &
                   a(7,6) * var(or-6,:,:) + a(8,6) * var(or-7,:,:) + &
                   a(9,6) * var(or-8,:,:) + a(10,6) * var(or-9,:,:) + &
                   a(11,6) * var(or-10,:,:) ) * idel / dx(or-5,:,:)
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                 ( a(1,5) * var(or,:,:) + a(2,5) * var(or-1,:,:) + &
                   a(3,5) * var(or-2,:,:) + a(4,5) * var(or-3,:,:) + &
                   a(5,5) * var(or-4,:,:) + a(6,5) * var(or-5,:,:) + &
                   a(7,5) * var(or-6,:,:) + a(8,5) * var(or-7,:,:) + &
                   a(9,5) * var(or-8,:,:) + a(10,5) * var(or-9,:,:) ) * &
                 idel / dx(or-4,:,:)
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                 ( a(1,4) * var(or,:,:) + a(2,4) * var(or-1,:,:) + &
                   a(3,4) * var(or-2,:,:) + a(4,4) * var(or-3,:,:) + &
                   a(5,4) * var(or-4,:,:) + a(6,4) * var(or-5,:,:) + &
                   a(7,4) * var(or-6,:,:) + a(8,4) * var(or-7,:,:) + &
                   a(9,4) * var(or-8,:,:) ) * idel / dx(or-3,:,:)
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                 ( a(1,3) * var(or,:,:) + a(2,3) * var(or-1,:,:) + &
                   a(3,3) * var(or-2,:,:) + a(4,3) * var(or-3,:,:) + &
                   a(5,3) * var(or-4,:,:) + a(6,3) * var(or-5,:,:) + &
                   a(7,3) * var(or-6,:,:) + a(8,3) * var(or-7,:,:) ) * &
                 idel / dx(or-2,:,:)
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) + &
                   a(5,2) * var(or-4,:,:) + a(6,2) * var(or-5,:,:) + &
                   a(7,2) * var(or-6,:,:) ) * idel / dx(or-1,:,:)
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) + a(4,1) * var(or-3,:,:) + &
                   a(5,1) * var(or-4,:,:) + a(6,1) * var(or-5,:,:) ) * &
                 idel / dx(or,:,:)
!$OMP END PARALLEL WORKSHARE

    ir = or - 8
  end if

!$OMP PARALLEL WORKSHARE
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -252.0_wp * var(il:ir,:,:) + &
                      210.0_wp * ( var(il-1:ir-1,:,:) + &
                                   var(il+1:ir+1,:,:) ) - &
                      120.0_wp * ( var(il-2:ir-2,:,:) + &
                                   var(il+2:ir+2,:,:) ) + &
                       45.0_wp * ( var(il-3:ir-3,:,:) + &
                                   var(il+3:ir+3,:,:) ) - &
                       10.0_wp * ( var(il-4:ir-4,:,:) + &
                                   var(il+4:ir+4,:,:) ) + &
                                 ( var(il-5:ir-5,:,:) + &
                                   var(il+5:ir+5,:,:) ) ) * idel / dx(il:ir,:,:)
!$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 1024
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3);

!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) + &
                     a(5,1) * var(:,5+ol,:) + a(6,1) * var(:,6+ol,:) ) * &
                   idel / dy(:,1+ol,:)
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) + a(6,2) * var(:,6+ol,:) + &
                     a(7,2) * var(:,7+ol,:) ) * idel / dy(:,2+ol,:)
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) + &
                     a(7,3) * var(:,7+ol,:) + a(8,3) * var(:,8+ol,:) ) * &
                   idel / dy(:,3+ol,:)
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) + a(8,4) * var(:,8+ol,:) + &
                     a(9,4) * var(:,9+ol,:) ) * idel / dy(:,4+ol,:)
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(1,5) * var(:,1+ol,:) + a(2,5) * var(:,2+ol,:) + &
                     a(3,5) * var(:,3+ol,:) + a(4,5) * var(:,4+ol,:) + &
                     a(5,5) * var(:,5+ol,:) + a(6,5) * var(:,6+ol,:) + &
                     a(7,5) * var(:,7+ol,:) + a(8,5) * var(:,8+ol,:) + &
                     a(9,5) * var(:,9+ol,:) + a(10,5) * var(:,10+ol,:) ) * &
                   idel / dy(:,5+ol,:)
      rhs(:,6+ol,:) = rhs(:,6+ol,:) + &
                   ( a(1,6) * var(:,1+ol,:) + a(2,6) * var(:,2+ol,:) + &
                     a(3,6) * var(:,3+ol,:) + a(4,6) * var(:,4+ol,:) + &
                     a(5,6) * var(:,5+ol,:) + a(6,6) * var(:,6+ol,:) + &
                     a(7,6) * var(:,7+ol,:) + a(8,6) * var(:,8+ol,:) + &
                     a(9,6) * var(:,9+ol,:) + a(10,6) * var(:,10+ol,:) + &
                     a(11,6) * var(:,11+ol,:) ) * idel / dy(:,6+ol,:)
      rhs(:,7+ol,:) = rhs(:,7+ol,:) + &
                   ( a(2,7) * var(:,2+ol,:) + a(3,7) * var(:,3+ol,:) + &
                     a(4,7) * var(:,4+ol,:) + a(5,7) * var(:,5+ol,:) + &
                     a(6,7) * var(:,6+ol,:) + a(7,7) * var(:,7+ol,:) + &
                     a(8,7) * var(:,8+ol,:) + a(9,7) * var(:,9+ol,:) + &
                     a(10,7) * var(:,10+ol,:) + a(11,7) * var(:,11+ol,:) + &
                     a(12,7) * var(:,12+ol,:) ) * idel / dy(:,7+ol,:)
      rhs(:,8+ol,:) = rhs(:,8+ol,:) + &
                   ( a(3,8) * var(:,3+ol,:) + a(4,8) * var(:,4+ol,:) + &
                     a(5,8) * var(:,5+ol,:) + a(6,8) * var(:,6+ol,:) + &
                     a(7,8) * var(:,7+ol,:) + a(8,8) * var(:,8+ol,:) + &
                     a(9,8) * var(:,9+ol,:) + a(10,8) * var(:,10+ol,:) + &
                     a(11,8) * var(:,11+ol,:) + a(12,8) * var(:,12+ol,:) + &
                     a(13,8) * var(:,13+ol,:) ) * idel / dy(:,8+ol,:)
!$OMP END PARALLEL WORKSHARE
  
      jl = 9 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)

!$OMP PARALLEL WORKSHARE
      rhs(:,or-7,:) = rhs(:,or-7,:) + &
                   ( a(3,8) * var(:,or-2,:) + a(4,8) * var(:,or-3,:) + &
                     a(5,8) * var(:,or-4,:) + a(6,8) * var(:,or-5,:) + &
                     a(7,8) * var(:,or-6,:) + a(8,8) * var(:,or-7,:) + &
                     a(9,8) * var(:,or-8,:) + a(10,8) * var(:,or-9,:) + &
                     a(11,8) * var(:,or-10,:) + a(12,8) * var(:,or-11,:) + &
                     a(13,8) * var(:,or-12,:) ) * idel / dy(:,or-7,:)
      rhs(:,or-6,:) = rhs(:,or-6,:) + &
                   ( a(2,7) * var(:,or-1,:) + a(3,7) * var(:,or-2,:) + &
                     a(4,7) * var(:,or-3,:) + a(5,7) * var(:,or-4,:) + &
                     a(6,7) * var(:,or-5,:) + a(7,7) * var(:,or-6,:) + &
                     a(8,7) * var(:,or-7,:) + a(9,7) * var(:,or-8,:) + &
                     a(10,7) * var(:,or-9,:) + a(11,7) * var(:,or-10,:) + &
                     a(12,7) * var(:,or-11,:) ) * idel / dy(:,or-6,:)
      rhs(:,or-5,:) = rhs(:,or-5,:) + &
                   ( a(1,6) * var(:,or,:) + a(2,6) * var(:,or-1,:) + &
                     a(3,6) * var(:,or-2,:) + a(4,6) * var(:,or-3,:) + &
                     a(5,6) * var(:,or-4,:) + a(6,6) * var(:,or-5,:) + &
                     a(7,6) * var(:,or-6,:) + a(8,6) * var(:,or-7,:) + &
                     a(9,6) * var(:,or-8,:) + a(10,6) * var(:,or-9,:) + &
                     a(11,6) * var(:,or-10,:) ) * idel / dy(:,or-5,:)
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                   ( a(1,5) * var(:,or,:) + a(2,5) * var(:,or-1,:) + &
                     a(3,5) * var(:,or-2,:) + a(4,5) * var(:,or-3,:) + &
                     a(5,5) * var(:,or-4,:) + a(6,5) * var(:,or-5,:) + &
                     a(7,5) * var(:,or-6,:) + a(8,5) * var(:,or-7,:) + &
                     a(9,5) * var(:,or-8,:) + a(10,5) * var(:,or-9,:) ) * &
                   idel / dy(:,or-4,:)
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                   ( a(1,4) * var(:,or,:) + a(2,4) * var(:,or-1,:) + &
                     a(3,4) * var(:,or-2,:) + a(4,4) * var(:,or-3,:) + &
                     a(5,4) * var(:,or-4,:) + a(6,4) * var(:,or-5,:) + &
                     a(7,4) * var(:,or-6,:) + a(8,4) * var(:,or-7,:) + &
                     a(9,4) * var(:,or-8,:) ) * idel / dy(:,or-3,:)
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                   ( a(1,3) * var(:,or,:) + a(2,3) * var(:,or-1,:) + &
                     a(3,3) * var(:,or-2,:) + a(4,3) * var(:,or-3,:) + &
                     a(5,3) * var(:,or-4,:) + a(6,3) * var(:,or-5,:) + &
                     a(7,3) * var(:,or-6,:) + a(8,3) * var(:,or-7,:) ) * &
                   idel / dy(:,or-2,:)
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) + &
                     a(5,2) * var(:,or-4,:) + a(6,2) * var(:,or-5,:) + &
                     a(7,2) * var(:,or-6,:) ) * idel / dy(:,or-1,:)
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) + a(4,1) * var(:,or-3,:) + &
                     a(5,1) * var(:,or-4,:) + a(6,1) * var(:,or-5,:) ) * &
                   idel / dy(:,or,:)
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 8
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -252.0_wp * var(:,jl:jr,:) + &
                        210.0_wp * ( var(:,jl-1:jr-1,:) + &
                                     var(:,jl+1:jr+1,:) ) - &
                        120.0_wp * ( var(:,jl-2:jr-2,:) + &
                                     var(:,jl+2:jr+2,:) ) + &
                         45.0_wp * ( var(:,jl-3:jr-3,:) + &
                                     var(:,jl+3:jr+3,:) ) - &
                         10.0_wp * ( var(:,jl-4:jr-4,:) + &
                                     var(:,jl+4:jr+4,:) ) + &
                                   ( var(:,jl-5:jr-5,:) + &
                                     var(:,jl+5:jr+5,:) ) ) * &
                     idel / dy(:,jl:jr,:)
!$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 1024
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5);

!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) + &
                     a(5,1) * var(:,:,5+ol) + a(6,1) * var(:,:,6+ol ) ) * &
                   idel / dz(:,:,1+ol)
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) + a(6,2) * var(:,:,6+ol) + &
                     a(7,2) * var(:,:,7+ol ) ) * idel / dz(:,:,2+ol)
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) + &
                     a(7,3) * var(:,:,7+ol) + a(8,3) * var(:,:,8+ol ) ) * &
                   idel / dz(:,:,3+ol)
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) + a(8,4) * var(:,:,8+ol) + &
                     a(9,4) * var(:,:,9+ol ) ) * idel / dz(:,:,4+ol)
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(1,5) * var(:,:,1+ol) + a(2,5) * var(:,:,2+ol) + &
                     a(3,5) * var(:,:,3+ol) + a(4,5) * var(:,:,4+ol) + &
                     a(5,5) * var(:,:,5+ol) + a(6,5) * var(:,:,6+ol) + &
                     a(7,5) * var(:,:,7+ol) + a(8,5) * var(:,:,8+ol) + &
                     a(9,5) * var(:,:,9+ol) + a(10,5) * var(:,:,10+ol ) ) * &
                   idel / dz(:,:,5+ol)
      rhs(:,:,6+ol) = rhs(:,:,6+ol) + &
                   ( a(1,6) * var(:,:,1+ol) + a(2,6) * var(:,:,2+ol) + &
                     a(3,6) * var(:,:,3+ol) + a(4,6) * var(:,:,4+ol) + &
                     a(5,6) * var(:,:,5+ol) + a(6,6) * var(:,:,6+ol) + &
                     a(7,6) * var(:,:,7+ol) + a(8,6) * var(:,:,8+ol) + &
                     a(9,6) * var(:,:,9+ol) + a(10,6) * var(:,:,10+ol) + &
                     a(11,6) * var(:,:,11+ol ) ) * idel / dz(:,:,6+ol)
      rhs(:,:,7+ol) = rhs(:,:,7+ol) + &
                   ( a(2,7) * var(:,:,2+ol) + a(3,7) * var(:,:,3+ol) + &
                     a(4,7) * var(:,:,4+ol) + a(5,7) * var(:,:,5+ol) + &
                     a(6,7) * var(:,:,6+ol) + a(7,7) * var(:,:,7+ol) + &
                     a(8,7) * var(:,:,8+ol) + a(9,7) * var(:,:,9+ol) + &
                     a(10,7) * var(:,:,10+ol) + a(11,7) * var(:,:,11+ol) + &
                     a(12,7) * var(:,:,12+ol ) ) * idel / dz(:,:,7+ol)
      rhs(:,:,8+ol) = rhs(:,:,8+ol) + &
                   ( a(3,8) * var(:,:,3+ol) + a(4,8) * var(:,:,4+ol) + &
                     a(5,8) * var(:,:,5+ol) + a(6,8) * var(:,:,6+ol) + &
                     a(7,8) * var(:,:,7+ol) + a(8,8) * var(:,:,8+ol) + &
                     a(9,8) * var(:,:,9+ol) + a(10,8) * var(:,:,10+ol) + &
                     a(11,8) * var(:,:,11+ol) + a(12,8) * var(:,:,12+ol) + &
                     a(13,8) * var(:,:,13+ol ) ) * idel / dz(:,:,8+ol)
!$OMP END PARALLEL WORKSHARE
  
      kl = 9 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-7) = rhs(:,:,or-7) + &
                   ( a(3,8) * var(:,:,or-2) + a(4,8) * var(:,:,or-3) + &
                     a(5,8) * var(:,:,or-4) + a(6,8) * var(:,:,or-5) + &
                     a(7,8) * var(:,:,or-6) + a(8,8) * var(:,:,or-7) + &
                     a(9,8) * var(:,:,or-8) + a(10,8) * var(:,:,or-9) + &
                     a(11,8) * var(:,:,or-10) + a(12,8) * var(:,:,or-11) + &
                     a(13,8) * var(:,:,or-12) ) * idel / dz(:,:,or-7)
      rhs(:,:,or-6) = rhs(:,:,or-6) + &
                   ( a(2,7) * var(:,:,or-1) + a(3,7) * var(:,:,or-2) + &
                     a(4,7) * var(:,:,or-3) + a(5,7) * var(:,:,or-4) + &
                     a(6,7) * var(:,:,or-5) + a(7,7) * var(:,:,or-6) + &
                     a(8,7) * var(:,:,or-7) + a(9,7) * var(:,:,or-8) + &
                     a(10,7) * var(:,:,or-9) + a(11,7) * var(:,:,or-10) + &
                     a(12,7) * var(:,:,or-11) ) * idel / dz(:,:,or-6)
      rhs(:,:,or-5) = rhs(:,:,or-5) + &
                   ( a(1,6) * var(:,:,or) + a(2,6) * var(:,:,or-1) + &
                     a(3,6) * var(:,:,or-2) + a(4,6) * var(:,:,or-3) + &
                     a(5,6) * var(:,:,or-4) + a(6,6) * var(:,:,or-5) + &
                     a(7,6) * var(:,:,or-6) + a(8,6) * var(:,:,or-7) + &
                     a(9,6) * var(:,:,or-8) + a(10,6) * var(:,:,or-9) + &
                     a(11,6) * var(:,:,or-10) ) * idel / dz(:,:,or-5)
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                   ( a(1,5) * var(:,:,or) + a(2,5) * var(:,:,or-1) + &
                     a(3,5) * var(:,:,or-2) + a(4,5) * var(:,:,or-3) + &
                     a(5,5) * var(:,:,or-4) + a(6,5) * var(:,:,or-5) + &
                     a(7,5) * var(:,:,or-6) + a(8,5) * var(:,:,or-7) + &
                     a(9,5) * var(:,:,or-8) + a(10,5) * var(:,:,or-9) ) * &
                   idel / dz(:,:,or-4)
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                   ( a(1,4) * var(:,:,or) + a(2,4) * var(:,:,or-1) + &
                     a(3,4) * var(:,:,or-2) + a(4,4) * var(:,:,or-3) + &
                     a(5,4) * var(:,:,or-4) + a(6,4) * var(:,:,or-5) + &
                     a(7,4) * var(:,:,or-6) + a(8,4) * var(:,:,or-7) + &
                     a(9,4) * var(:,:,or-8) ) * idel / dz(:,:,or-3)
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                   ( a(1,3) * var(:,:,or) + a(2,3) * var(:,:,or-1) + &
                     a(3,3) * var(:,:,or-2) + a(4,3) * var(:,:,or-3) + &
                     a(5,3) * var(:,:,or-4) + a(6,3) * var(:,:,or-5) + &
                     a(7,3) * var(:,:,or-6) + a(8,3) * var(:,:,or-7) ) * &
                   idel / dz(:,:,or-2)
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) + &
                     a(5,2) * var(:,:,or-4) + a(6,2) * var(:,:,or-5) + &
                     a(7,2) * var(:,:,or-6) ) * idel / dz(:,:,or-1)
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) + a(4,1) * var(:,:,or-3) + &
                     a(5,1) * var(:,:,or-4) + a(6,1) * var(:,:,or-5) ) * &
                   idel / dz(:,:,or)
!$OMP END PARALLEL WORKSHARE

      kr = or - 8
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -252.0_wp * var(:,:,kl:kr) + &
                        210.0_wp * ( var(:,:,kl-1:kr-1) + &
                                     var(:,:,kl+1:kr+1) ) - &
                        120.0_wp * ( var(:,:,kl-2:kr-2) + &
                                     var(:,:,kl+2:kr+2) ) + &
                         45.0_wp * ( var(:,:,kl-3:kr-3) + &
                                     var(:,:,kl+3:kr+3) ) - &
                         10.0_wp * ( var(:,:,kl-4:kr-4) + &
                                     var(:,:,kl+4:kr+4) ) + &
                                   ( var(:,:,kl-5:kr-5) + &
                                     var(:,:,kl+5:kr+5) ) ) * &
                     idel / dz(:,:,kl:kr)
!$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(13,8), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -3.3910872088637970174997113084967416241083103770745_wp
    a(2,1) = 16.955436044318985087498556542483708120541551885372_wp
    a(3,1) = -33.910872088637970174997113084967416241083103770745_wp
    a(4,1) = 33.910872088637970174997113084967416241083103770745_wp
    a(5,1) = -16.955436044318985087498556542483708120541551885372_wp
    a(6,1) = 3.3910872088637970174997113084967416241083103770745_wp
    a(7,1) = 0.0_wp
    a(8,1) = 0.0_wp
    a(9,1) = 0.0_wp
    a(10,1) = 0.0_wp
    a(11,1) = 0.0_wp
    a(12,1) = 0.0_wp
    a(13,1) = 0.0_wp
    a(1,2) = 3.2771399440263630592058029074141137010783820566473_wp
    a(2,2) = -17.041127708937087907870175118553391245607586694566_wp
    a(3,2) = 36.048539384289993651263831981555250711862202623121_wp
    a(4,2) = -39.325679328316356710469634888969364412940584679768_wp
    a(5,2) = 22.939979608184541414440620351898795907548674396531_wp
    a(6,2) = -6.5542798880527261184116058148282274021567641132947_wp
    a(7,2) = 0.65542798880527261184116058148282274021567641132947_wp
    a(8,2) = 0.0_wp
    a(9,2) = 0.0_wp
    a(10,2) = 0.0_wp
    a(11,2) = 0.0_wp
    a(12,2) = 0.0_wp
    a(13,2) = 0.0_wp
    a(1,3) = -38.842059631038967294446317614758441308222147295410_wp
    a(2,3) = 213.63132797071432011945474688117142719522181012475_wp
    a(3,3) = -489.40995135109098791002360194595636048359905592216_wp
    a(4,3) = 602.05192428110399306391792302875584027744328307885_wp
    a(5,3) = -427.26265594142864023890949376234285439044362024951_wp
    a(6,3) = 174.78926833967535282500842926641298588699966282934_wp
    a(7,3) = -38.842059631038967294446317614758441308222147295410_wp
    a(8,3) = 3.8842059631038967294446317614758441308222147295410_wp
    a(9,3) = 0.0_wp
    a(10,3) = 0.0_wp
    a(11,3) = 0.0_wp
    a(12,3) = 0.0_wp
    a(13,3) = 0.0_wp
    a(1,4) = 5.5613835719414344378807953109542612676331552744485_wp
    a(2,4) = -33.368301431648606627284771865725567605798931646691_wp
    a(3,4) = 86.201445365092233787152327319791049648313906753952_wp
    a(4,4) = -125.68726872587641829610597402756630464850930920254_wp
    a(5,4) = 114.00836322479940597655630387456235598647968312619_wp
    a(6,4) = -66.736602863297213254569543731451135211597863293382_wp
    a(7,4) = 25.026226073736454970463578899294175704349198735018_wp
    a(8,4) = -5.5613835719414344378807953109542612676331552744485_wp
    a(9,4) = 0.55613835719414344378807953109542612676331552744485_wp
    a(10,4) = 0.0_wp
    a(11,4) = 0.0_wp
    a(12,4) = 0.0_wp
    a(13,4) = 0.0_wp
    a(1,5) = -12.115101476661536355654081268132755978592914829047_wp
    a(2,5) = 84.805710336630754489578568876929291850150403803330_wp
    a(3,5) = -266.53223248655379982438978789892063152904412623904_wp
    a(4,5) = 496.71916054312299058181733199344299512230950799093_wp
    a(5,5) = -608.17809412840912505383487966026435012536432441817_wp
    a(6,5) = 508.83426201978452693747141326157575110090242281998_wp
    a(7,5) = -290.76243543987687253569795043518614348622995589713_wp
    a(8,5) = 109.03591328995382720088673141319480380733623346142_wp
    a(9,5) = -24.230202953323072711308162536265511957185829658094_wp
    a(10,5) = 2.4230202953323072711308162536265511957185829658094_wp
    a(11,5) = 0.0_wp
    a(12,5) = 0.0_wp
    a(13,5) = 0.0_wp
    a(1,6) = 0.78217600900123184961735065035840034142603567514089_wp
    a(2,6) = -7.8217600900123184961735065035840034142603567514089_wp
    a(3,6) = 35.197920405055433232780779266128015364171605381340_wp
    a(4,6) = -93.861121080147821954082078043008040971124281016906_wp
    a(5,6) = 164.25696189025868841964363657526407169946749177959_wp
    a(6,6) = -197.10835426831042610357236389031688603936099013550_wp
    a(7,6) = 164.25696189025868841964363657526407169946749177959_wp
    a(8,6) = -93.861121080147821954082078043008040971124281016906_wp
    a(9,6) = 35.197920405055433232780779266128015364171605381340_wp
    a(10,6) = -7.8217600900123184961735065035840034142603567514089_wp
    a(11,6) = 0.78217600900123184961735065035840034142603567514089_wp
    a(12,6) = 0.0_wp
    a(13,6) = 0.0_wp
    a(1,7) = 0.0_wp
    a(2,7) = 1.0830767761393601764536458481012280421614377748694_wp
    a(3,7) = -10.830767761393601764536458481012280421614377748694_wp
    a(4,7) = 48.738454926271207940414063164555261897264699869122_wp
    a(5,7) = -129.96921313672322117443750177214736505937253298433_wp
    a(6,7) = 227.44612298926563705526562810125788885390193272257_wp
    a(7,7) = -272.93534758711876446631875372150946662468231926708_wp
    a(8,7) = 227.44612298926563705526562810125788885390193272257_wp
    a(9,7) = -129.96921313672322117443750177214736505937253298433_wp
    a(10,7) = 48.738454926271207940414063164555261897264699869122_wp
    a(11,7) = -10.830767761393601764536458481012280421614377748694_wp
    a(12,7) = 1.0830767761393601764536458481012280421614377748694_wp
    a(13,7) = 0.0_wp
    a(1,8) = 0.0_wp
    a(2,8) = 0.0_wp
    a(3,8) = 0.99075245444434671889501396229410272246695863420506_wp
    a(4,8) = -9.9075245444434671889501396229410272246695863420506_wp
    a(5,8) = 44.583860449995602350275628303234622511013138539228_wp
    a(6,8) = -118.89029453332160626740167547529232669603503610461_wp
    a(7,8) = 208.05801543331281096795293208176157171806131318306_wp
    a(8,8) = -249.66961851997537316154351849811388606167357581967_wp
    a(9,8) = 208.05801543331281096795293208176157171806131318306_wp
    a(10,8) = -118.89029453332160626740167547529232669603503610461_wp
    a(11,8) = 44.583860449995602350275628303234622511013138539228_wp
    a(12,8) = -9.9075245444434671889501396229410272246695863420506_wp
    a(13,8) = 0.99075245444434671889501396229410272246695863420506_wp

  end subroutine set_coeff

end subroutine dissipation_8_4_delta
