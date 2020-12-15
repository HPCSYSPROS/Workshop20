! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_4_2_alt (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

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
  CCTK_REAL, dimension(7,4) :: a
  CCTK_REAL :: idel

  CCTK_INT :: imin, imax, jmin, jmax, kmin, kmax
  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or
  CCTK_INT :: i, j, k

  call set_coeff ( a )

  idel = epsilon / ( 64 * delta(1) )

  ! Determine extent of non-ghost region
  imin =  1; if (bb(1)==0) imin = imin + gsize(1)
  jmin =  1; if (bb(3)==0) jmin = jmin + gsize(2)
  kmin =  1; if (bb(5)==0) kmin = kmin + gsize(3)
  imax = ni; if (bb(2)==0) imax = imax - gsize(1)
  jmax = nj; if (bb(4)==0) jmax = jmax - gsize(2)
  kmax = nk; if (bb(6)==0) kmax = kmax - gsize(3)

  if ( bb(1) == 0 ) then
    il = imin
  else
    ol = offset(1)
!$omp parallel do private(j,k)
    do k=kmin,kmax
    do j=jmin,jmax
    rhs(1+ol,j,k) = rhs(1+ol,j,k) + &
                 ( a(1,1) * var(1+ol,j,k) + a(2,1) * var(2+ol,j,k) + &
                   a(3,1) * var(3+ol,j,k) + a(4,1) * var(4+ol,j,k) ) * idel
    rhs(2+ol,j,k) = rhs(2+ol,j,k) + &
                 ( a(1,2) * var(1+ol,j,k) + a(2,2) * var(2+ol,j,k) + &
                   a(3,2) * var(3+ol,j,k) + a(4,2) * var(4+ol,j,k) + &
                   a(5,2) * var(5+ol,j,k) ) * idel
    rhs(3+ol,j,k) = rhs(3+ol,j,k) + &
                 ( a(1,3) * var(1+ol,j,k) + a(2,3) * var(2+ol,j,k) + &
                   a(3,3) * var(3+ol,j,k) + a(4,3) * var(4+ol,j,k) + &
                   a(5,3) * var(5+ol,j,k) + a(6,3) * var(6+ol,j,k) ) * idel
    rhs(4+ol,j,k) = rhs(4+ol,j,k) + &
                 ( a(1,4) * var(1+ol,j,k) + a(2,4) * var(2+ol,j,k) + &
                   a(3,4) * var(3+ol,j,k) + a(4,4) * var(4+ol,j,k) + &
                   a(5,4) * var(5+ol,j,k) + a(6,4) * var(6+ol,j,k) + &
                   a(7,4) * var(7+ol,j,k) ) * idel
    end do
    end do
!$omp end parallel do

    il = 5 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = imax
  else
    or = ni - offset(2)
!$omp parallel do private(j,k)
    do k=kmin,kmax
    do j=jmin,jmax
    rhs(or-3,j,k) = rhs(or-3,j,k) + &
                 ( a(1,4) * var(or,j,k) + a(2,4) * var(or-1,j,k) + &
                   a(3,4) * var(or-2,j,k) + a(4,4) * var(or-3,j,k) + &
                   a(5,4) * var(or-4,j,k) + a(6,4) * var(or-5,j,k) + &
                   a(7,4) * var(or-6,j,k) ) * idel
    rhs(or-2,j,k) = rhs(or-2,j,k) + &
                 ( a(1,3) * var(or,j,k) + a(2,3) * var(or-1,j,k) + &
                   a(3,3) * var(or-2,j,k) + a(4,3) * var(or-3,j,k) + &
                   a(5,3) * var(or-4,j,k) + a(6,3) * var(or-5,j,k) ) * idel
    rhs(or-1,j,k) = rhs(or-1,j,k) + &
                 ( a(1,2) * var(or,j,k) + a(2,2) * var(or-1,j,k) + &
                   a(3,2) * var(or-2,j,k) + a(4,2) * var(or-3,j,k) + &
                   a(5,2) * var(or-4,j,k) ) * idel
    rhs(or,j,k)   = rhs(or,j,k) + &
                 ( a(1,1) * var(or,j,k) + a(2,1) * var(or-1,j,k) + &
                   a(3,1) * var(or-2,j,k) + a(4,1) * var(or-3,j,k) ) * idel
    end do
    end do
!$omp end parallel do

    ir = or - 4
  end if
!$omp parallel do private(i,j,k)
  do k=kmin,kmax
  do j=jmin,jmax
  do i=il,ir
  rhs(i,j,k) = rhs(i,j,k) + &
                   ( -20.0_wp * var(i,j,k) + &
                      15.0_wp * ( var(i-1,j,k) + &
                                  var(i+1,j,k) ) - &
                       6.0_wp * ( var(i-2,j,k) + &
                                  var(i+2,j,k) ) + &
                                ( var(i-3,j,k) + &
                                  var(i+3,j,k) ) ) * idel
  end do
  end do
  end do
!$omp end parallel do

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 64 * delta(2) )
  
    if ( bb(3) == 0 ) then
      jl = jmin
    else
      ol = offset(3)
!$omp parallel do private(i,k)
      do k=kmin,kmax
      do i=imin,imax
      rhs(i,1+ol,k) = rhs(i,1+ol,k) + &
                   ( a(1,1) * var(i,1+ol,k) + a(2,1) * var(i,2+ol,k) + &
                     a(3,1) * var(i,3+ol,k) + a(4,1) * var(i,4+ol,k) ) * idel
      rhs(i,2+ol,k) = rhs(i,2+ol,k) + &
                   ( a(1,2) * var(i,1+ol,k) + a(2,2) * var(i,2+ol,k) + &
                     a(3,2) * var(i,3+ol,k) + a(4,2) * var(i,4+ol,k) + &
                     a(5,2) * var(i,5+ol,k) ) * idel
      rhs(i,3+ol,k) = rhs(i,3+ol,k) + &
                   ( a(1,3) * var(i,1+ol,k) + a(2,3) * var(i,2+ol,k) + &
                     a(3,3) * var(i,3+ol,k) + a(4,3) * var(i,4+ol,k) + &
                     a(5,3) * var(i,5+ol,k) + a(6,3) * var(i,6+ol,k) ) * idel
      rhs(i,4+ol,k) = rhs(i,4+ol,k) + &
                   ( a(1,4) * var(i,1+ol,k) + a(2,4) * var(i,2+ol,k) + &
                     a(3,4) * var(i,3+ol,k) + a(4,4) * var(i,4+ol,k) + &
                     a(5,4) * var(i,5+ol,k) + a(6,4) * var(i,6+ol,k) + &
                     a(7,4) * var(i,7+ol,k) ) * idel
      end do
      end do
!$omp end parallel do
  
      jl = 5 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = jmax
    else
      or = nj - offset(4)
!$omp parallel do private(i,k)
      do k=kmin,kmax
      do i=imin,imax
      rhs(i,or-3,k) = rhs(i,or-3,k) + &
                   ( a(1,4) * var(i,or,k) + a(2,4) * var(i,or-1,k) + &
                     a(3,4) * var(i,or-2,k) + a(4,4) * var(i,or-3,k) + &
                     a(5,4) * var(i,or-4,k) + a(6,4) * var(i,or-5,k) + &
                     a(7,4) * var(i,or-6,k) ) * idel
      rhs(i,or-2,k) = rhs(i,or-2,k) + &
                   ( a(1,3) * var(i,or,k) + a(2,3) * var(i,or-1,k) + &
                     a(3,3) * var(i,or-2,k) + a(4,3) * var(i,or-3,k) + &
                     a(5,3) * var(i,or-4,k) + a(6,3) * var(i,or-5,k) ) * idel
      rhs(i,or-1,k) = rhs(i,or-1,k) + &
                   ( a(1,2) * var(i,or,k) + a(2,2) * var(i,or-1,k) + &
                     a(3,2) * var(i,or-2,k) + a(4,2) * var(i,or-3,k) + &
                     a(5,2) * var(i,or-4,k) ) * idel
      rhs(i,or,k)   = rhs(i,or,k) + &
                   ( a(1,1) * var(i,or,k) + a(2,1) * var(i,or-1,k) + &
                     a(3,1) * var(i,or-2,k) + a(4,1) * var(i,or-3,k) ) * idel
      end do
      end do
!$omp end parallel do
  
      jr = or - 4
    end if
!$omp parallel do private(i,j,k)
    do k=kmin,kmax
    do j=jl,jr
    do i=imin,imax
    rhs(i,j,k) = rhs(i,j,k) + &
                     ( -20.0_wp * var(i,j,k) + &
                        15.0_wp * ( var(i,j-1,k) + &
                                    var(i,j+1,k) ) - &
                         6.0_wp * ( var(i,j-2,k) + &
                                    var(i,j+2,k) ) + &
                                  ( var(i,j-3,k) + &
                                    var(i,j+3,k) ) ) * idel
    end do
    end do
    end do
!$omp end parallel do
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 64 * delta(3) )
  
    if ( bb(5) == 0 ) then
      kl = kmin
    else
      ol = offset(5)
!$omp parallel do private(i,j)
      do j=jmin,jmax
      do i=imin,imax
      rhs(i,j,1+ol) = rhs(i,j,1+ol) + &
                   ( a(1,1) * var(i,j,1+ol) + a(2,1) * var(i,j,2+ol) + &
                     a(3,1) * var(i,j,3+ol) + a(4,1) * var(i,j,4+ol) ) * idel
      rhs(i,j,2+ol) = rhs(i,j,2+ol) + &
                   ( a(1,2) * var(i,j,1+ol) + a(2,2) * var(i,j,2+ol) + &
                     a(3,2) * var(i,j,3+ol) + a(4,2) * var(i,j,4+ol) + &
                     a(5,2) * var(i,j,5+ol) ) * idel
      rhs(i,j,3+ol) = rhs(i,j,3+ol) + &
                   ( a(1,3) * var(i,j,1+ol) + a(2,3) * var(i,j,2+ol) + &
                     a(3,3) * var(i,j,3+ol) + a(4,3) * var(i,j,4+ol) + &
                     a(5,3) * var(i,j,5+ol) + a(6,3) * var(i,j,6+ol) ) * idel
      rhs(i,j,4+ol) = rhs(i,j,4+ol) + &
                   ( a(1,4) * var(i,j,1+ol) + a(2,4) * var(i,j,2+ol) + &
                     a(3,4) * var(i,j,3+ol) + a(4,4) * var(i,j,4+ol) + &
                     a(5,4) * var(i,j,5+ol) + a(6,4) * var(i,j,6+ol) + &
                     a(7,4) * var(i,j,7+ol) ) * idel
      end do
      end do
!$omp end parallel do
  
      kl = 5 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = kmax
    else
      or = nk - offset(6)
!$omp parallel do private(i,j)
      do j=jmin,jmax
      do i=imin,imax
      rhs(i,j,or-3) = rhs(i,j,or-3) + &
                   ( a(1,4) * var(i,j,or) + a(2,4) * var(i,j,or-1) + &
                     a(3,4) * var(i,j,or-2) + a(4,4) * var(i,j,or-3) + &
                     a(5,4) * var(i,j,or-4) + a(6,4) * var(i,j,or-5) + &
                     a(7,4) * var(i,j,or-6) ) * idel
      rhs(i,j,or-2) = rhs(i,j,or-2) + &
                   ( a(1,3) * var(i,j,or) + a(2,3) * var(i,j,or-1) + &
                     a(3,3) * var(i,j,or-2) + a(4,3) * var(i,j,or-3) + &
                     a(5,3) * var(i,j,or-4) + a(6,3) * var(i,j,or-5) ) * idel
      rhs(i,j,or-1) = rhs(i,j,or-1) + &
                   ( a(1,2) * var(i,j,or) + a(2,2) * var(i,j,or-1) + &
                     a(3,2) * var(i,j,or-2) + a(4,2) * var(i,j,or-3) + &
                     a(5,2) * var(i,j,or-4) ) * idel
      rhs(i,j,or)   = rhs(i,j,or) + &
                   ( a(1,1) * var(i,j,or) + a(2,1) * var(i,j,or-1) + &
                     a(3,1) * var(i,j,or-2) + a(4,1) * var(i,j,or-3) ) * idel
      end do
      end do
!$omp end parallel do
  
      kr = or - 4
    end if
!$omp parallel do private(i,j,k)
    do k=kl,kr
    do j=jmin,jmax
    do i=imin,imax
    rhs(i,j,k) = rhs(i,j,k) + &
                     ( -20.0_wp * var(i,j,k) + &
                        15.0_wp * ( var(i,j,k-1) + &
                                    var(i,j,k+1) ) - &
                         6.0_wp * ( var(i,j,k-2) + &
                                    var(i,j,k+2) ) + &
                                  ( var(i,j,k-3) + &
                                    var(i,j,k+3) ) ) * idel
    end do
    end do
    end do
!$omp end parallel do
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(7,4), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -2.8235294117647058823529411764705882352941176470588_wp
    a(2,1) = 8.4705882352941176470588235294117647058823529411765_wp
    a(3,1) = -8.4705882352941176470588235294117647058823529411765_wp
    a(4,1) = 2.8235294117647058823529411764705882352941176470588_wp
    a(5,1) = 0.0_wp
    a(6,1) = 0.0_wp
    a(7,1) = 0.0_wp
    a(1,2) = 2.4406779661016949152542372881355932203389830508475_wp
    a(2,2) = -8.1355932203389830508474576271186440677966101694915_wp
    a(3,2) = 9.7627118644067796610169491525423728813559322033898_wp
    a(4,2) = -4.8813559322033898305084745762711864406779661016949_wp
    a(5,2) = 0.81355932203389830508474576271186440677966101694915_wp
    a(6,2) = 0.0_wp
    a(7,2) = 0.0_wp
    a(1,3) = -3.3488372093023255813953488372093023255813953488372_wp
    a(2,3) = 13.395348837209302325581395348837209302325581395349_wp
    a(3,3) = -21.209302325581395348837209302325581395348837209302_wp
    a(4,3) = 16.744186046511627906976744186046511627906976744186_wp
    a(5,3) = -6.6976744186046511627906976744186046511627906976744_wp
    a(6,3) = 1.1162790697674418604651162790697674418604651162791_wp
    a(7,3) = 0.0_wp
    a(1,4) = 0.97959183673469387755102040816326530612244897959184_wp
    a(2,4) = -5.8775510204081632653061224489795918367346938775510_wp
    a(3,4) = 14.693877551020408163265306122448979591836734693878_wp
    a(4,4) = -19.591836734693877551020408163265306122448979591837_wp
    a(5,4) = 14.693877551020408163265306122448979591836734693878_wp
    a(6,4) = -5.8775510204081632653061224489795918367346938775510_wp
    a(7,4) = 0.97959183673469387755102040816326530612244897959184_wp

  end subroutine set_coeff

end subroutine dissipation_4_2_alt

subroutine dissipation_4_2_delta (var, ni, nj, nk, bb, gsize, offset, &
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
  CCTK_REAL, dimension(7,4) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  idel = epsilon / 64

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)
!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) ) * &
                 idel / dx(1+ol,:,:)
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) ) * idel / dx(2+ol,:,:)
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) ) * &
                 idel / dx(3+ol,:,:)
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) ) * idel / dx(4+ol,:,:)
!$OMP END PARALLEL WORKSHARE

    il = 5 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)
!$OMP PARALLEL WORKSHARE
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                 ( a(1,4) * var(or,:,:) + a(2,4) * var(or-1,:,:) + &
                   a(3,4) * var(or-2,:,:) + a(4,4) * var(or-3,:,:) + &
                   a(5,4) * var(or-4,:,:) + a(6,4) * var(or-5,:,:) + &
                   a(7,4) * var(or-6,:,:) ) * idel / dx(or-3,:,:)
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                 ( a(1,3) * var(or,:,:) + a(2,3) * var(or-1,:,:) + &
                   a(3,3) * var(or-2,:,:) + a(4,3) * var(or-3,:,:) + &
                   a(5,3) * var(or-4,:,:) + a(6,3) * var(or-5,:,:) ) * &
                 idel / dx(or-2,:,:)
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) + &
                   a(5,2) * var(or-4,:,:) ) * idel / dx(or-1,:,:)
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) + a(4,1) * var(or-3,:,:) ) * &
                 idel / dx(or,:,:)
!$OMP END PARALLEL WORKSHARE

    ir = or - 4
  end if
!$OMP PARALLEL WORKSHARE
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -20.0_wp * var(il:ir,:,:) + &
                      15.0_wp * ( var(il-1:ir-1,:,:) + &
                                  var(il+1:ir+1,:,:) ) - &
                       6.0_wp * ( var(il-2:ir-2,:,:) + &
                                  var(il+2:ir+2,:,:) ) + &
                                ( var(il-3:ir-3,:,:) + &
                                  var(il+3:ir+3,:,:) ) ) * idel / dx(il:ir,:,:)
!$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 64
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)
!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) ) * &
                   idel / dy(:,1+ol,:)
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) ) * idel / dy(:,2+ol,:)
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) ) * &
                   idel / dy(:,3+ol,:)
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) ) * idel / dy(:,4+ol,:)
!$OMP END PARALLEL WORKSHARE
  
      jl = 5 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)
!$OMP PARALLEL WORKSHARE
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                   ( a(1,4) * var(:,or,:) + a(2,4) * var(:,or-1,:) + &
                     a(3,4) * var(:,or-2,:) + a(4,4) * var(:,or-3,:) + &
                     a(5,4) * var(:,or-4,:) + a(6,4) * var(:,or-5,:) + &
                     a(7,4) * var(:,or-6,:) ) * idel / dy(:,or-3,:)
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                   ( a(1,3) * var(:,or,:) + a(2,3) * var(:,or-1,:) + &
                     a(3,3) * var(:,or-2,:) + a(4,3) * var(:,or-3,:) + &
                     a(5,3) * var(:,or-4,:) + a(6,3) * var(:,or-5,:) ) * &
                   idel / dy(:,or-2,:)
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) + &
                     a(5,2) * var(:,or-4,:) ) * idel / dy(:,or-1,:)
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) + a(4,1) * var(:,or-3,:) ) * &
                   idel / dy(:,or,:)
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 4
    end if
!$OMP PARALLEL WORKSHARE
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -20.0_wp * var(:,jl:jr,:) + &
                        15.0_wp * ( var(:,jl-1:jr-1,:) + &
                                    var(:,jl+1:jr+1,:) ) - &
                         6.0_wp * ( var(:,jl-2:jr-2,:) + &
                                    var(:,jl+2:jr+2,:) ) + &
                                  ( var(:,jl-3:jr-3,:) + &
                                    var(:,jl+3:jr+3,:) ) ) * &
                     idel / dy(:,jl:jr,:)
!$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 64
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)
!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) ) * &
                   idel / dz(:,:,1+ol)
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) ) * idel / dz(:,:,2+ol)
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) ) * &
                   idel / dz(:,:,3+ol)
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) ) * idel / dz(:,:,4+ol)
!$OMP END PARALLEL WORKSHARE
  
      kl = 5 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)
!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                   ( a(1,4) * var(:,:,or) + a(2,4) * var(:,:,or-1) + &
                     a(3,4) * var(:,:,or-2) + a(4,4) * var(:,:,or-3) + &
                     a(5,4) * var(:,:,or-4) + a(6,4) * var(:,:,or-5) + &
                     a(7,4) * var(:,:,or-6) ) * idel / dz(:,:,or-3)
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                   ( a(1,3) * var(:,:,or) + a(2,3) * var(:,:,or-1) + &
                     a(3,3) * var(:,:,or-2) + a(4,3) * var(:,:,or-3) + &
                     a(5,3) * var(:,:,or-4) + a(6,3) * var(:,:,or-5) ) * &
                   idel / dz(:,:,or-2)
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) + &
                     a(5,2) * var(:,:,or-4) ) * idel / dz(:,:,or-1)
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) + a(4,1) * var(:,:,or-3) ) * &
                   idel / dz(:,:,or)
!$OMP END PARALLEL WORKSHARE
  
      kr = or - 4
    end if
!$OMP PARALLEL WORKSHARE
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -20.0_wp * var(:,:,kl:kr) + &
                        15.0_wp * ( var(:,:,kl-1:kr-1) + &
                                    var(:,:,kl+1:kr+1) ) - &
                         6.0_wp * ( var(:,:,kl-2:kr-2) + &
                                    var(:,:,kl+2:kr+2) ) + &
                                  ( var(:,:,kl-3:kr-3) + &
                                    var(:,:,kl+3:kr+3) ) ) * &
                     idel / dz(:,:,kl:kr)
!$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(7,4), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -2.8235294117647058823529411764705882352941176470588_wp
    a(2,1) = 8.4705882352941176470588235294117647058823529411765_wp
    a(3,1) = -8.4705882352941176470588235294117647058823529411765_wp
    a(4,1) = 2.8235294117647058823529411764705882352941176470588_wp
    a(5,1) = 0.0_wp
    a(6,1) = 0.0_wp
    a(7,1) = 0.0_wp
    a(1,2) = 2.4406779661016949152542372881355932203389830508475_wp
    a(2,2) = -8.1355932203389830508474576271186440677966101694915_wp
    a(3,2) = 9.7627118644067796610169491525423728813559322033898_wp
    a(4,2) = -4.8813559322033898305084745762711864406779661016949_wp
    a(5,2) = 0.81355932203389830508474576271186440677966101694915_wp
    a(6,2) = 0.0_wp
    a(7,2) = 0.0_wp
    a(1,3) = -3.3488372093023255813953488372093023255813953488372_wp
    a(2,3) = 13.395348837209302325581395348837209302325581395349_wp
    a(3,3) = -21.209302325581395348837209302325581395348837209302_wp
    a(4,3) = 16.744186046511627906976744186046511627906976744186_wp
    a(5,3) = -6.6976744186046511627906976744186046511627906976744_wp
    a(6,3) = 1.1162790697674418604651162790697674418604651162791_wp
    a(7,3) = 0.0_wp
    a(1,4) = 0.97959183673469387755102040816326530612244897959184_wp
    a(2,4) = -5.8775510204081632653061224489795918367346938775510_wp
    a(3,4) = 14.693877551020408163265306122448979591836734693878_wp
    a(4,4) = -19.591836734693877551020408163265306122448979591837_wp
    a(5,4) = 14.693877551020408163265306122448979591836734693878_wp
    a(6,4) = -5.8775510204081632653061224489795918367346938775510_wp
    a(7,4) = 0.97959183673469387755102040816326530612244897959184_wp

  end subroutine set_coeff

end subroutine dissipation_4_2_delta
