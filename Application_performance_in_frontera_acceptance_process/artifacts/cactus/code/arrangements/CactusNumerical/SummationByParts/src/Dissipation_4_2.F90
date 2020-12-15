! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_4_2 (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

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
  CCTK_REAL, dimension(6,4) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or
  CCTK_INT :: i, j, k

  call set_coeff ( a )

  if ( scale_with_h > 0 ) then
    idel = epsilon / ( 16 * delta(1) )
  else
    idel = epsilon / 16
  end if

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)
!$omp parallel do private(j,k)
    do k=1,nk
    do j=1,nj
    rhs(1+ol,j,k) = rhs(1+ol,j,k) + &
                 ( a(1,1) * var(1+ol,j,k) + a(2,1) * var(2+ol,j,k) + &
                   a(3,1) * var(3+ol,j,k) ) * idel
    rhs(2+ol,j,k) = rhs(2+ol,j,k) + &
                 ( a(1,2) * var(1+ol,j,k) + a(2,2) * var(2+ol,j,k) + &
                   a(3,2) * var(3+ol,j,k) + a(4,2) * var(4+ol,j,k) ) * idel
    rhs(3+ol,j,k) = rhs(3+ol,j,k) + &
                 ( a(1,3) * var(1+ol,j,k) + a(2,3) * var(2+ol,j,k) + &
                   a(3,3) * var(3+ol,j,k) + a(4,3) * var(4+ol,j,k) + &
                   a(5,3) * var(5+ol,j,k) ) * idel
    rhs(4+ol,j,k) = rhs(4+ol,j,k) + &
                 ( a(2,4) * var(2+ol,j,k) + a(3,4) * var(3+ol,j,k) + &
                   a(4,4) * var(4+ol,j,k) + a(5,4) * var(5+ol,j,k) + &
                   a(6,4) * var(6+ol,j,k) ) * idel
    end do
    end do
!$omp end parallel do

    il = 5 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)
!$omp parallel do private(j,k)
    do k=1,nk
    do j=1,nj
    rhs(or-3,j,k) = rhs(or-3,j,k) + &
                 ( a(2,4) * var(or-1,j,k) + a(3,4) * var(or-2,j,k) + &
                   a(4,4) * var(or-3,j,k) + a(5,4) * var(or-4,j,k) + &
                   a(6,4) * var(or-5,j,k) ) * idel
    rhs(or-2,j,k) = rhs(or-2,j,k) + &
                 ( a(1,3) * var(or,j,k) + a(2,3) * var(or-1,j,k) + &
                   a(3,3) * var(or-2,j,k) + a(4,3) * var(or-3,j,k) + &
                   a(5,3) * var(or-4,j,k) ) * idel
    rhs(or-1,j,k) = rhs(or-1,j,k) + &
                 ( a(1,2) * var(or,j,k) + a(2,2) * var(or-1,j,k) + &
                   a(3,2) * var(or-2,j,k) + a(4,2) * var(or-3,j,k) ) * idel
    rhs(or,j,k)   = rhs(or,j,k) + &
                 ( a(1,1) * var(or,j,k) + a(2,1) * var(or-1,j,k) + &
                   a(3,1) * var(or-2,j,k) ) * idel
    end do
    end do
!$omp end parallel do

    ir = or - 4
  end if
!$omp parallel do private(i,j,k)
  do k=1,nk
  do j=1,nj
  do i=il,ir
  rhs(i,j,k) = rhs(i,j,k) + &
                   ( -6.0_wp * var(i,j,k) + &
                      4.0_wp * ( var(i-1,j,k) + &
                                 var(i+1,j,k) ) - &
                               ( var(i-2,j,k) + &
                                 var(i+2,j,k) ) ) * idel
  end do
  end do
  end do
!$omp end parallel do

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 16 * delta(2) )
    else
      idel = epsilon / 16
    end if
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)
!$omp parallel do private(i,k)
      do k=1,nk
      do i=1,ni
      rhs(i,1+ol,k) = rhs(i,1+ol,k) + &
                   ( a(1,1) * var(i,1+ol,k) + a(2,1) * var(i,2+ol,k) + &
                     a(3,1) * var(i,3+ol,k) ) * idel
      rhs(i,2+ol,k) = rhs(i,2+ol,k) + &
                   ( a(1,2) * var(i,1+ol,k) + a(2,2) * var(i,2+ol,k) + &
                     a(3,2) * var(i,3+ol,k) + a(4,2) * var(i,4+ol,k) ) * idel
      rhs(i,3+ol,k) = rhs(i,3+ol,k) + &
                   ( a(1,3) * var(i,1+ol,k) + a(2,3) * var(i,2+ol,k) + &
                     a(3,3) * var(i,3+ol,k) + a(4,3) * var(i,4+ol,k) + &
                     a(5,3) * var(i,5+ol,k) ) * idel
      rhs(i,4+ol,k) = rhs(i,4+ol,k) + &
                   ( a(2,4) * var(i,2+ol,k) + a(3,4) * var(i,3+ol,k) + &
                     a(4,4) * var(i,4+ol,k) + a(5,4) * var(i,5+ol,k) + &
                     a(6,4) * var(i,6+ol,k) ) * idel
      end do
      end do
!$omp end parallel do
  
      jl = 5 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)
!$omp parallel do private(i,k)
      do k=1,nk
      do i=1,ni
      rhs(i,or-3,k) = rhs(i,or-3,k) + &
                   ( a(2,4) * var(i,or-1,k) + a(3,4) * var(i,or-2,k) + &
                     a(4,4) * var(i,or-3,k) + a(5,4) * var(i,or-4,k) + &
                     a(6,4) * var(i,or-5,k) ) * idel
      rhs(i,or-2,k) = rhs(i,or-2,k) + &
                   ( a(1,3) * var(i,or,k) + a(2,3) * var(i,or-1,k) + &
                     a(3,3) * var(i,or-2,k) + a(4,3) * var(i,or-3,k) + &
                     a(5,3) * var(i,or-4,k) ) * idel
      rhs(i,or-1,k) = rhs(i,or-1,k) + &
                   ( a(1,2) * var(i,or,k) + a(2,2) * var(i,or-1,k) + &
                     a(3,2) * var(i,or-2,k) + a(4,2) * var(i,or-3,k) ) * idel
      rhs(i,or,k)   = rhs(i,or,k) + &
                   ( a(1,1) * var(i,or,k) + a(2,1) * var(i,or-1,k) + &
                     a(3,1) * var(i,or-2,k) ) * idel
      end do
      end do
!$omp end parallel do
  
      jr = or - 4
    end if
!$omp parallel do private(i,j,k)
    do k=1,nk
    do j=jl,jr
    do i=1,ni
    rhs(i,j,k) = rhs(i,j,k) + &
                     ( -6.0_wp * var(i,j,k) + &
                        4.0_wp * ( var(i,j-1,k) + &
                                   var(i,j+1,k) ) - &
                                 ( var(i,j-2,k) + &
                                   var(i,j+2,k) ) ) * idel
    end do
    end do
    end do
!$omp end parallel do
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 16 * delta(3) )
    else
      idel = epsilon / 16
    end if
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)
!$omp parallel do private(i,j)
      do j=1,nj
      do i=1,ni
      rhs(i,j,1+ol) = rhs(i,j,1+ol) + &
                   ( a(1,1) * var(i,j,1+ol) + a(2,1) * var(i,j,2+ol) + &
                     a(3,1) * var(i,j,3+ol) ) * idel
      rhs(i,j,2+ol) = rhs(i,j,2+ol) + &
                   ( a(1,2) * var(i,j,1+ol) + a(2,2) * var(i,j,2+ol) + &
                     a(3,2) * var(i,j,3+ol) + a(4,2) * var(i,j,4+ol) ) * idel
      rhs(i,j,3+ol) = rhs(i,j,3+ol) + &
                   ( a(1,3) * var(i,j,1+ol) + a(2,3) * var(i,j,2+ol) + &
                     a(3,3) * var(i,j,3+ol) + a(4,3) * var(i,j,4+ol) + &
                     a(5,3) * var(i,j,5+ol) ) * idel
      rhs(i,j,4+ol) = rhs(i,j,4+ol) + &
                   ( a(2,4) * var(i,j,2+ol) + a(3,4) * var(i,j,3+ol) + &
                     a(4,4) * var(i,j,4+ol) + a(5,4) * var(i,j,5+ol) + &
                     a(6,4) * var(i,j,6+ol) ) * idel
      end do
      end do
!$omp end parallel do
  
      kl = 5 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)
!$omp parallel do private(i,j)
      do j=1,nj
      do i=1,ni
      rhs(i,j,or-3) = rhs(i,j,or-3) + &
                   ( a(2,4) * var(i,j,or-1) + a(3,4) * var(i,j,or-2) + &
                     a(4,4) * var(i,j,or-3) + a(5,4) * var(i,j,or-4) + &
                     a(6,4) * var(i,j,or-5) ) * idel
      rhs(i,j,or-2) = rhs(i,j,or-2) + &
                   ( a(1,3) * var(i,j,or) + a(2,3) * var(i,j,or-1) + &
                     a(3,3) * var(i,j,or-2) + a(4,3) * var(i,j,or-3) + &
                     a(5,3) * var(i,j,or-4) ) * idel
      rhs(i,j,or-1) = rhs(i,j,or-1) + &
                   ( a(1,2) * var(i,j,or) + a(2,2) * var(i,j,or-1) + &
                     a(3,2) * var(i,j,or-2) + a(4,2) * var(i,j,or-3) ) * idel
      rhs(i,j,or)   = rhs(i,j,or) + &
                   ( a(1,1) * var(i,j,or) + a(2,1) * var(i,j,or-1) + &
                     a(3,1) * var(i,j,or-2) ) * idel
      end do
      end do
!$omp end parallel do
  
      kr = or - 4
    end if
!$omp parallel do private(i,j,k)
    do k=kl,kr
    do j=1,nj
    do i=1,ni
    rhs(i,j,k) = rhs(i,j,k) + &
                     ( -6.0_wp * var(i,j,k) + &
                        4.0_wp * ( var(i,j,k-1) + &
                                   var(i,j,k+1) ) - &
                                 ( var(i,j,k-2) + &
                                   var(i,j,k+2) ) ) * idel
    end do
    end do
    end do
!$omp end parallel do
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(6,4), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -2.8235294117647058823529411764705882352941176470588_wp
    a(2,1) = 5.6470588235294117647058823529411764705882352941176_wp
    a(3,1) = -2.8235294117647058823529411764705882352941176470588_wp
    a(4,1) = 0.0_wp
    a(5,1) = 0.0_wp
    a(6,1) = 0.0_wp
    a(1,2) = 1.6271186440677966101694915254237288135593220338983_wp
    a(2,2) = -4.0677966101694915254237288135593220338983050847458_wp
    a(3,2) = 3.2542372881355932203389830508474576271186440677966_wp
    a(4,2) = -0.81355932203389830508474576271186440677966101694915_wp
    a(5,2) = 0.0_wp
    a(6,2) = 0.0_wp
    a(1,3) = -1.1162790697674418604651162790697674418604651162791_wp
    a(2,3) = 4.4651162790697674418604651162790697674418604651163_wp
    a(3,3) = -6.6976744186046511627906976744186046511627906976744_wp
    a(4,3) = 4.4651162790697674418604651162790697674418604651163_wp
    a(5,3) = -1.1162790697674418604651162790697674418604651162791_wp
    a(6,3) = 0.0_wp
    a(1,4) = 0.0_wp
    a(2,4) = -0.97959183673469387755102040816326530612244897959184_wp
    a(3,4) = 3.9183673469387755102040816326530612244897959183673_wp
    a(4,4) = -5.8775510204081632653061224489795918367346938775510_wp
    a(5,4) = 3.9183673469387755102040816326530612244897959183673_wp
    a(6,4) = -0.97959183673469387755102040816326530612244897959184_wp

  end subroutine set_coeff

end subroutine dissipation_4_2
