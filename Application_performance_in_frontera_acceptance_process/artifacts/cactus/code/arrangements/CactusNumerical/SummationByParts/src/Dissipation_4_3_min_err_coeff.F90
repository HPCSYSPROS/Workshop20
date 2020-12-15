! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine dissipation_4_3_opt (var, lsh, gsh, lbnd, bb, gsize, offset, &
                                delta, epsilon, dfl, npatches, patch, rhs)

  use dissipation_coeff

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, dimension(3), intent(in) ::  lsh, gsh, lbnd
  CCTK_INT :: ni, nj, nk
  CCTK_REAL, dimension(lsh(1),lsh(2),lsh(3)), intent(in) :: var
  CCTK_REAL, dimension(lsh(1),lsh(2),lsh(3)), intent(inout) :: rhs
  CCTK_INT, dimension(6), intent(in) :: bb
  CCTK_INT, dimension(3), intent(in) :: gsize
  CCTK_INT, dimension(6), intent(in) :: offset
  CCTK_REAL, dimension(3), intent(in) :: delta
  CCTK_REAL, intent(in) :: epsilon
  CCTK_REAL, dimension(3), intent(in) :: dfl
  CCTK_INT, intent(in) :: npatches, patch

  CCTK_REAL :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  integer :: i, j, k
  CCTK_REAL, dimension(:,:), allocatable :: atmp, d, b, h
  CCTK_REAL, dimension(:,:), pointer :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or
  CCTK_INT :: nii, njj, nkk

  ni = lsh(1); nj = lsh(2); nk = lsh(3)

  if ( first ) then
    allocate ( savedx(npatches), savedy(npatches), savedz(npatches) )
    savedx = .false.; savedy = .false.; savedz = .false.
    allocate ( xcoeff(npatches), ycoeff(npatches), zcoeff(npatches) )
    first = .false.
  end if

  if ( .not. savedx(patch) ) then
    nii = ni
    if ( bb(1) /= 0 ) then
      nii = nii - offset(1)
    endif
    if ( bb(2) /= 0 ) then
      nii = nii - offset(2)
    endif
    allocate ( atmp(nii,nii), d(nii,nii), b(nii,nii), h(nii,nii) )
    atmp = zero; d = zero; b = zero; h = zero

    call set_dmatrix ( d, bb(1:2) )

    call set_bmatrix ( b, bb(1:2), lsh(1), gsh(1)-offset(1)-offset(2), lbnd(1), delta(1), dfl(1) )
  
    call set_hmatrix ( h, bb(1:2) )


    atmp = -transpose ( matmul ( h, matmul ( transpose(d), matmul ( b, d ) ) ) )

    allocate ( xcoeff(patch)%coeff(nii,nii) )

    xcoeff(patch)%coeff = atmp

    savedx(patch) = .true.

    deallocate ( atmp, d, b, h )

  end if

  a => xcoeff(patch)%coeff

  if ( scale_with_h > 0 ) then
    idel = epsilon / ( 16 * delta(1) )
  else
    idel = epsilon / 16
  end if

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) + a(6,2) * var(6+ol,:,:) + &
                   a(7,2) * var(7+ol,:,:) ) * idel
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) + &
                   a(7,3) * var(7+ol,:,:) ) * idel
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) ) * idel
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(1,5) * var(1+ol,:,:) + a(2,5) * var(2+ol,:,:) + &
                   a(3,5) * var(3+ol,:,:) + a(4,5) * var(4+ol,:,:) + &
                   a(5,5) * var(5+ol,:,:) + a(6,5) * var(6+ol,:,:) + &
                   a(7,5) * var(7+ol,:,:) ) * idel
    il = 6 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                    ( a(nii-6,nii-4) * var(or-6,:,:) + &
                      a(nii-5,nii-4) * var(or-5,:,:) + &
                      a(nii-4,nii-4) * var(or-4,:,:) + &
                      a(nii-3,nii-4) * var(or-3,:,:) + &
                      a(nii-2,nii-4) * var(or-2,:,:) + &
                      a(nii-1,nii-4) * var(or-1,:,:) + &
                      a(nii,nii-4) * var(or,:,:) ) * idel
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                    ( a(nii-6,nii-3) * var(or-6,:,:) + &
                      a(nii-5,nii-3) * var(or-5,:,:) + &
                      a(nii-4,nii-3) * var(or-4,:,:) + &
                      a(nii-3,nii-3) * var(or-3,:,:) + &
                      a(nii-2,nii-3) * var(or-2,:,:) + &
                      a(nii-1,nii-3) * var(or-1,:,:) + &
                      a(nii,nii-3) * var(or,:,:) ) * idel
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                    ( a(nii-6,nii-2) * var(or-6,:,:) + &
                      a(nii-5,nii-2) * var(or-5,:,:) + &
                      a(nii-4,nii-2) * var(or-4,:,:) + &
                      a(nii-3,nii-2) * var(or-3,:,:) + &
                      a(nii-2,nii-2) * var(or-2,:,:) + &
                      a(nii-1,nii-2) * var(or-1,:,:) + &
                      a(nii,nii-2) * var(or,:,:) ) * idel
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                    ( a(nii-6,nii-1) * var(or-6,:,:) + &
                      a(nii-5,nii-1) * var(or-5,:,:) + &
                      a(nii-4,nii-1) * var(or-4,:,:) + &
                      a(nii-3,nii-1) * var(or-3,:,:) + &
                      a(nii-2,nii-1) * var(or-2,:,:) + &
                      a(nii-1,nii-1) * var(or-1,:,:) + &
                      a(nii,nii-1) * var(or,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                    ( a(nii-2,nii) * var(or-2,:,:) + &
                      a(nii-1,nii) * var(or-1,:,:) + &
                      a(nii,nii) * var(or,:,:) ) * idel
    ir = or - 5
  end if
  if (bb(1) /= 0) then
    do i = il, ir
        rhs(i,:,:) = rhs(i,:,:) + &
                       ( a(i-ol-2,i-ol) * var(i-2,:,:) + &
                         a(i-ol-1,i-ol) * var(i-1,:,:) + &
                         a(i-ol,i-ol) * var(i,:,:) + &
                         a(i-ol+1,i-ol) * var(i+1,:,:) + &
                         a(i-ol+2,i-ol) * var(i+2,:,:) ) * idel
    end do
  else
    do i = il, ir
        rhs(i,:,:) = rhs(i,:,:) + &
                       ( a(i-2,i) * var(i-2,:,:) + &
                         a(i-1,i) * var(i-1,:,:) + &
                         a(i,i) * var(i,:,:) + &
                         a(i+1,i) * var(i+1,:,:) + &
                         a(i+2,i) * var(i+2,:,:) ) * idel
    end do
  endif

  if ( zero_derivs_y == 0 ) then
    if ( .not. savedy(patch) ) then
      njj = nj
      if ( bb(3) /= 0 ) then
        njj = njj - offset(3)
      endif
      if ( bb(4) /= 0 ) then
        njj = njj - offset(4)
      endif
      allocate ( atmp(njj,njj), d(njj,njj), b(njj,njj), h(njj,njj) )
      atmp = zero; d = zero; b = zero; h = zero
  
      call set_dmatrix ( d, bb(3:4) )
  
      call set_bmatrix ( b, bb(3:4), lsh(2), gsh(2)-offset(3)-offset(4), lbnd(2), delta(2), dfl(2) )
    
      call set_hmatrix ( h, bb(3:4) )
  
      atmp = -transpose ( matmul ( h, matmul ( transpose(d), matmul ( b, d ) ) ) )
  
      allocate ( ycoeff(patch)%coeff(njj,njj) )

      ycoeff(patch)%coeff = atmp

      savedy(patch) = .true.

      deallocate ( atmp, d, b, h )

    end if

    a => ycoeff(patch)%coeff

    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 16 * delta(2) )
    else
      idel = epsilon / 16
    end if
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) + a(6,2) * var(:,6+ol,:) + &
                     a(7,2) * var(:,7+ol,:) ) * idel
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) + &
                     a(7,3) * var(:,7+ol,:) ) * idel
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) ) * idel
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(1,5) * var(:,1+ol,:) + a(2,5) * var(:,2+ol,:) + &
                     a(3,5) * var(:,3+ol,:) + a(4,5) * var(:,4+ol,:) + &
                     a(5,5) * var(:,5+ol,:) + a(6,5) * var(:,6+ol,:) + &
                     a(7,5) * var(:,7+ol,:) ) * idel
      jl = 6 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                      ( a(njj-6,njj-4) * var(:,or-6,:) + &
                        a(njj-5,njj-4) * var(:,or-5,:) + &
                        a(njj-4,njj-4) * var(:,or-4,:) + &
                        a(njj-3,njj-4) * var(:,or-3,:) + &
                        a(njj-2,njj-4) * var(:,or-2,:) + &
                        a(njj-1,njj-4) * var(:,or-1,:) + &
                        a(njj,njj-4) * var(:,or,:) ) * idel
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                      ( a(njj-6,njj-3) * var(:,or-6,:) + &
                        a(njj-5,njj-3) * var(:,or-5,:) + &
                        a(njj-4,njj-3) * var(:,or-4,:) + &
                        a(njj-3,njj-3) * var(:,or-3,:) + &
                        a(njj-2,njj-3) * var(:,or-2,:) + &
                        a(njj-1,njj-3) * var(:,or-1,:) + &
                        a(njj,njj-3) * var(:,or,:) ) * idel
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                      ( a(njj-6,njj-2) * var(:,or-6,:) + &
                        a(njj-5,njj-2) * var(:,or-5,:) + &
                        a(njj-4,njj-2) * var(:,or-4,:) + &
                        a(njj-3,njj-2) * var(:,or-3,:) + &
                        a(njj-2,njj-2) * var(:,or-2,:) + &
                        a(njj-1,njj-2) * var(:,or-1,:) + &
                        a(njj,njj-2) * var(:,or,:) ) * idel
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                      ( a(njj-6,njj-1) * var(:,or-6,:) + &
                        a(njj-5,njj-1) * var(:,or-5,:) + &
                        a(njj-4,njj-1) * var(:,or-4,:) + &
                        a(njj-3,njj-1) * var(:,or-3,:) + &
                        a(njj-2,njj-1) * var(:,or-2,:) + &
                        a(njj-1,njj-1) * var(:,or-1,:) + &
                        a(njj,njj-1) * var(:,or,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                      ( a(njj-2,njj) * var(:,or-2,:) + &
                        a(njj-1,njj) * var(:,or-1,:) + &
                        a(njj,njj) * var(:,or,:) ) * idel
      jr = or - 5
    end if
  
    if ( bb(3) /= 0 ) then
      do j = jl, jr
         rhs(:,j,:) = rhs(:,j,:) + &
                        ( a(j-ol-2,j-ol) * var(:,j-2,:) + &
                          a(j-ol-1,j-ol) * var(:,j-1,:) + &
                          a(j-ol,j-ol) * var(:,j,:) + &
                          a(j-ol+1,j-ol) * var(:,j+1,:) + &
                          a(j-ol+2,j-ol) * var(:,j+2,:) ) * idel
      end do
    else
      do j = jl, jr
         rhs(:,j,:) = rhs(:,j,:) + &
                        ( a(j-2,j) * var(:,j-2,:) + &
                          a(j-1,j) * var(:,j-1,:) + &
                          a(j,j) * var(:,j,:) + &
                          a(j+1,j) * var(:,j+1,:) + &
                          a(j+2,j) * var(:,j+2,:) ) * idel
      end do
    endif
  
  end if

  if ( zero_derivs_z == 0 ) then
    if ( .not. savedz(patch) ) then
      nkk = nk
      if ( bb(5) /= 0 ) then
        nkk = nkk - offset(5)
      endif
      if ( bb(6) /= 0 ) then
        nkk = nkk - offset(6)
      endif
      allocate ( atmp(nkk,nkk), d(nkk,nkk), b(nkk,nkk), h(nkk,nkk) )
      atmp = zero; d = zero; b = zero; h = zero
  
      call set_dmatrix ( d, bb(5:6) )
  
      call set_bmatrix ( b, bb(5:6), lsh(3), gsh(3)-offset(5)-offset(6), lbnd(3), delta(3), dfl(3) )
  
      call set_hmatrix ( h, bb(5:6) )
  
      atmp = -transpose ( matmul ( h, matmul ( transpose(d), matmul ( b, d ) ) ) )
  
      allocate ( zcoeff(patch)%coeff(nkk,nkk) )

      zcoeff(patch)%coeff = atmp

      savedz(patch) = .true.

      deallocate ( atmp, d, b, h )

    end if

    a => zcoeff(patch)%coeff

    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 16 * delta(3) )
    else
      idel = epsilon / 16
    end if
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) + a(6,2) * var(:,:,6+ol) + &
                     a(7,2) * var(:,:,7+ol) ) * idel
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) + &
                     a(7,3) * var(:,:,7+ol) ) * idel
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) ) * idel
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(1,5) * var(:,:,1+ol) + a(2,5) * var(:,:,2+ol) + &
                     a(3,5) * var(:,:,3+ol) + a(4,5) * var(:,:,4+ol) + &
                     a(5,5) * var(:,:,5+ol) + a(6,5) * var(:,:,6+ol) + &
                     a(7,5) * var(:,:,7+ol) ) * idel
      kl = 6 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                      ( a(nkk-6,nkk-4) * var(:,:,or-6) + &
                        a(nkk-5,nkk-4) * var(:,:,or-5) + &
                        a(nkk-4,nkk-4) * var(:,:,or-4) + &
                        a(nkk-3,nkk-4) * var(:,:,or-3) + &
                        a(nkk-2,nkk-4) * var(:,:,or-2) + &
                        a(nkk-1,nkk-4) * var(:,:,or-1) + &
                        a(nkk,nkk-4) * var(:,:,or) ) * idel
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                      ( a(nkk-6,nkk-3) * var(:,:,or-6) + &
                        a(nkk-5,nkk-3) * var(:,:,or-5) + &
                        a(nkk-4,nkk-3) * var(:,:,or-4) + &
                        a(nkk-3,nkk-3) * var(:,:,or-3) + &
                        a(nkk-2,nkk-3) * var(:,:,or-2) + &
                        a(nkk-1,nkk-3) * var(:,:,or-1) + &
                        a(nkk,nkk-3) * var(:,:,or) ) * idel
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                      ( a(nkk-6,nkk-2) * var(:,:,or-6) + &
                        a(nkk-5,nkk-2) * var(:,:,or-5) + &
                        a(nkk-4,nkk-2) * var(:,:,or-4) + &
                        a(nkk-3,nkk-2) * var(:,:,or-3) + &
                        a(nkk-2,nkk-2) * var(:,:,or-2) + &
                        a(nkk-1,nkk-2) * var(:,:,or-1) + &
                        a(nkk,nkk-2) * var(:,:,or) ) * idel
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                      ( a(nkk-6,nkk-1) * var(:,:,or-6) + &
                        a(nkk-5,nkk-1) * var(:,:,or-5) + &
                        a(nkk-4,nkk-1) * var(:,:,or-4) + &
                        a(nkk-3,nkk-1) * var(:,:,or-3) + &
                        a(nkk-2,nkk-1) * var(:,:,or-2) + &
                        a(nkk-1,nkk-1) * var(:,:,or-1) + &
                        a(nkk,nkk-1) * var(:,:,or) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                      ( a(nkk-2,nkk) * var(:,:,or-2) + &
                        a(nkk-1,nkk) * var(:,:,or-1) + &
                        a(nkk,nkk) * var(:,:,or) ) * idel
      kr = or - 5
    end if
  
    if ( bb(5) /= 0) then
      do k = kl, kr
         rhs(:,:,k) = rhs(:,:,k) + &
                        ( a(k-ol-2,k-ol) * var(:,:,k-2) + &
                          a(k-ol-1,k-ol) * var(:,:,k-1) + &
                          a(k-ol,k-ol) * var(:,:,k) + &
                          a(k-ol+1,k-ol) * var(:,:,k+1) + &
                          a(k-ol+2,k-ol) * var(:,:,k+2) ) * idel
      end do
    else
      do k = kl, kr
         rhs(:,:,k) = rhs(:,:,k) + &
                        ( a(k-2,k) * var(:,:,k-2) + &
                          a(k-1,k) * var(:,:,k-1) + &
                          a(k,k) * var(:,:,k) + &
                          a(k+1,k) * var(:,:,k+1) + &
                          a(k+2,k) * var(:,:,k+2) ) * idel
      end do
    endif
  
  end if

contains

  subroutine set_dmatrix ( d, bb )

    implicit none

    CCTK_REAL, dimension(:,:), intent(out) :: d
    CCTK_INT, dimension(2), intent(in) :: bb
    CCTK_INT :: n
    CCTK_REAL, dimension(3), save :: ac = (/ 1.0_wp, -2.0_wp, 1.0_wp /)
    CCTK_INT :: i

    n = size(d,1)
    if ( bb(1) == 0 ) then
      d(1,1:2) = ac(2:3)
    else
      d(1,1:3) = ac
    end if
    if ( bb(2) == 0 ) then
      d(n,n-1:n) = ac(1:2)
    else
      d(n,n-2:n) = ac
    end if
    do i = 2, n-1
      d(i,i-1:i+1) = ac
    end do
  end subroutine set_dmatrix

  subroutine set_bmatrix ( b, bb, lsh, gsh, lbnd, h, dfl )

    implicit none

    CCTK_REAL, dimension(:,:), intent(out) :: b
    CCTK_INT, dimension(2), intent(in) :: bb
    CCTK_INT, intent(in) :: lsh, gsh, lbnd
    CCTK_REAL, intent(in) :: h, dfl
    CCTK_INT :: n, i
    CCTK_REAL :: xp

    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    n = size(b,1)
    xp = ( gsh - 1 ) * dfl - 1
    do i = 1, n
      b(i,i) = bf ( real(i+lbnd-1,wp), xp, real(gsh-1,wp), h )
    end do

  end subroutine set_bmatrix
    
  subroutine set_hmatrix ( h, bb )

    implicit none

    CCTK_REAL, dimension(:,:), intent(out) :: h
    CCTK_INT, dimension(2), intent(in) :: bb
    CCTK_INT :: n, i

    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    n = size(h,1)
    do i = 1, n
      h(i,i) = 1.0_wp
    end do
    if ( bb(1) /= 0 ) then
      h(1:5,1:5) = sigma ( )
    end if
    if ( bb(2) /= 0 ) then
      h(n:n-4:-1,n:n-4:-1) = sigma ( )
    end if

  end subroutine set_hmatrix
    
  function bf ( x, xp, xmax, h )

    implicit none

    CCTK_REAL :: bf
    CCTK_REAL, intent(in) :: x, xp, xmax, h

    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    if ( x < xp ) then
      bf = ( 1.0_wp - h ) / xp * x + h
    else if ( x > xmax - xp ) then
      bf = ( h - 1.0_wp ) / xp * ( x + xp - xmax ) + 1
    else
      bf = 1.0_wp
    end if

  end function bf

  function sigma ( )

    implicit none

    CCTK_REAL, dimension(5,5) :: sigma

    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    sigma(1,1) = 4.186595370392226897362216859769846226369_wp
    sigma(2,1) = 0.0_wp
    sigma(3,1) = 0.0_wp
    sigma(4,1) = 0.0_wp
    sigma(5,1) = 0.0_wp
    sigma(1,2) = 0.0_wp
    sigma(2,2) = 0.6725191921225620731888714836983116420871_wp
    sigma(3,2) = 0.3613418181134949259370502966736306984367_wp
    sigma(4,2) = -0.2021316117293899791481674539631879662707_wp
    sigma(5,2) = 0.03455320708729270824077678274955265350304_wp
    sigma(1,3) = 0.0_wp
    sigma(2,3) = 0.3613418181134949259370502966736306984367_wp
    sigma(3,3) = 0.7206133711630147057720442098623847362950_wp
    sigma(4,3) = 0.1376472340546569368321616389764958792591_wp
    sigma(5,3) = -0.04136405531324488624637892257286207044784_wp
    sigma(1,4) = 0.0_wp
    sigma(2,4) = -0.2021316117293899791481674539631879662707_wp
    sigma(3,4) = 0.1376472340546569368321616389764958792591_wp
    sigma(4,4) = 0.9578653607931026822074133441449909689509_wp
    sigma(5,4) = 0.02069353627247161734563597102894256809696_wp
    sigma(1,5) = 0.0_wp
    sigma(2,5) = 0.03455320708729270824077678274955265350304_wp
    sigma(3,5) = -0.04136405531324488624637892257286207044784_wp
    sigma(4,5) = 0.02069353627247161734563597102894256809696_wp
    sigma(5,5) = 0.9908272703370861473007798925906968380654_wp

  end function sigma

end subroutine dissipation_4_3_opt
