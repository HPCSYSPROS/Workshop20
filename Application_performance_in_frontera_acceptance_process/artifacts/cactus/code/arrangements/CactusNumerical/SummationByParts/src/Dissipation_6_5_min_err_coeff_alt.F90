! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine dissipation_6_5_alt (var, lsh, gsh, lbnd, bb, gsize, offset, &
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

    call set_bmatrix ( b, bb(1:2), lsh(1), gsh(1)-offset(1)-offset(2), lbnd(1), delta(1), dfl(1) )

    call set_hmatrix ( h, bb(1:2) )

    call set_dmatrix ( d, bb(1:2) )

    atmp = -transpose ( matmul ( h, matmul ( transpose(d), matmul ( b, d ) ) ) )

    allocate ( xcoeff(patch)%coeff(nii,nii) )

    xcoeff(patch)%coeff = atmp

    savedx(patch) = .true.

    deallocate ( atmp, d, b, h )

  end if

  a => xcoeff(patch)%coeff

  idel = epsilon /  ( 256.0_wp * delta(1) )

  if ( bb(1) == 0 ) then

    il = 1 + gsize(1)

  else
    ol = offset(1)

    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) + &
                   a(5,1) * var(5+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) + a(6,2) * var(6+ol,:,:) + &
                   a(7,2) * var(7+ol,:,:) + a(8,2) * var(8+ol,:,:) + &
                   a(9,2) * var(9+ol,:,:) + a(10,2) * var(10+ol,:,:) + &
                   a(11,2) * var(11+ol,:,:) ) * idel
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) + &
                   a(7,3) * var(7+ol,:,:) + a(8,3) * var(8+ol,:,:) + &
                   a(9,3) * var(9+ol,:,:) + a(10,3) * var(10+ol,:,:) + &
                   a(11,3) * var(11+ol,:,:) ) * idel
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) + a(8,4) * var(8+ol,:,:) + &
                   a(9,4) * var(9+ol,:,:) + a(10,4) * var(10+ol,:,:) + &
                   a(11,4) * var(11+ol,:,:) ) * idel
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(1,5) * var(1+ol,:,:) + a(2,5) * var(2+ol,:,:) + &
                   a(3,5) * var(3+ol,:,:) + a(4,5) * var(4+ol,:,:) + &
                   a(5,5) * var(5+ol,:,:) + a(6,5) * var(6+ol,:,:) + &
                   a(7,5) * var(7+ol,:,:) + a(8,5) * var(8+ol,:,:) + &
                   a(9,5) * var(9+ol,:,:) + a(10,5) * var(10+ol,:,:) + &
                   a(11,5) * var(11+ol,:,:) ) * idel
    rhs(6+ol,:,:) = rhs(6+ol,:,:) + &
                 ( a(1,6) * var(1+ol,:,:) + a(2,6) * var(2+ol,:,:) + &
                   a(3,6) * var(3+ol,:,:) + a(4,6) * var(4+ol,:,:) + &
                   a(5,6) * var(5+ol,:,:) + a(6,6) * var(6+ol,:,:) + &
                   a(7,6) * var(7+ol,:,:) + a(8,6) * var(8+ol,:,:) + &
                   a(9,6) * var(9+ol,:,:) + a(10,6) * var(10+ol,:,:) + &
                   a(11,6) * var(11+ol,:,:) ) * idel
    rhs(7+ol,:,:) = rhs(7+ol,:,:) + &
                 ( a(1,7) * var(1+ol,:,:) + a(2,7) * var(2+ol,:,:) + &
                   a(3,7) * var(3+ol,:,:) + a(4,7) * var(4+ol,:,:) + &
                   a(5,7) * var(5+ol,:,:) + a(6,7) * var(6+ol,:,:) + &
                   a(7,7) * var(7+ol,:,:) + a(8,7) * var(8+ol,:,:) + &
                   a(9,7) * var(9+ol,:,:) + a(10,7) * var(10+ol,:,:) + &
                   a(11,7) * var(11+ol,:,:) ) * idel

    il = 8 + ol

  end if

  if ( bb(2) == 0 ) then

    ir = ni - gsize(1)

  else
    or = ni - offset(2)

    rhs(or-6,:,:) = rhs(or-6,:,:) + &
                    ( a(nii-10,nii-6) * var(or-10,:,:) + &
                      a(nii-9,nii-6) * var(or-9,:,:) + &
                      a(nii-8,nii-6) * var(or-8,:,:) + &
                      a(nii-7,nii-6) * var(or-7,:,:) + &
                      a(nii-6,nii-6) * var(or-6,:,:) + &
                      a(nii-5,nii-6) * var(or-5,:,:) + &
                      a(nii-4,nii-6) * var(or-4,:,:) + &
                      a(nii-3,nii-6) * var(or-3,:,:) + &
                      a(nii-2,nii-6) * var(or-2,:,:) + &
                      a(nii-1,nii-6) * var(or-1,:,:) + &
                      a(nii,nii-6) * var(or,:,:) ) * idel
    rhs(or-5,:,:) = rhs(or-5,:,:) + &
                    ( a(nii-10,nii-5) * var(or-10,:,:) + &
                      a(nii-9,nii-5) * var(or-9,:,:) + &
                      a(nii-8,nii-5) * var(or-8,:,:) + &
                      a(nii-7,nii-5) * var(or-7,:,:) + &
                      a(nii-6,nii-5) * var(or-6,:,:) + &
                      a(nii-5,nii-5) * var(or-5,:,:) + &
                      a(nii-4,nii-5) * var(or-4,:,:) + &
                      a(nii-3,nii-5) * var(or-3,:,:) + &
                      a(nii-2,nii-5) * var(or-2,:,:) + &
                      a(nii-1,nii-5) * var(or-1,:,:) + &
                      a(nii,nii-5) * var(or,:,:) ) * idel
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                    ( a(nii-10,nii-4) * var(or-10,:,:) + &
                      a(nii-9,nii-4) * var(or-9,:,:) + &
                      a(nii-8,nii-4) * var(or-8,:,:) + &
                      a(nii-7,nii-4) * var(or-7,:,:) + &
                      a(nii-6,nii-4) * var(or-6,:,:) + &
                      a(nii-5,nii-4) * var(or-5,:,:) + &
                      a(nii-4,nii-4) * var(or-4,:,:) + &
                      a(nii-3,nii-4) * var(or-3,:,:) + &
                      a(nii-2,nii-4) * var(or-2,:,:) + &
                      a(nii-1,nii-4) * var(or-1,:,:) + &
                      a(nii,nii-4) * var(or,:,:) ) * idel
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                    ( a(nii-10,nii-3) * var(or-10,:,:) + &
                      a(nii-9,nii-3) * var(or-9,:,:) + &
                      a(nii-8,nii-3) * var(or-8,:,:) + &
                      a(nii-7,nii-3) * var(or-7,:,:) + &
                      a(nii-6,nii-3) * var(or-6,:,:) + &
                      a(nii-5,nii-3) * var(or-5,:,:) + &
                      a(nii-4,nii-3) * var(or-4,:,:) + &
                      a(nii-3,nii-3) * var(or-3,:,:) + &
                      a(nii-2,nii-3) * var(or-2,:,:) + &
                      a(nii-1,nii-3) * var(or-1,:,:) + &
                      a(nii,nii-3) * var(or,:,:) ) * idel
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                    ( a(nii-10,nii-2) * var(or-10,:,:) + &
                      a(nii-9,nii-2) * var(or-9,:,:) + &
                      a(nii-8,nii-2) * var(or-8,:,:) + &
                      a(nii-7,nii-2) * var(or-7,:,:) + &
                      a(nii-6,nii-2) * var(or-6,:,:) + &
                      a(nii-5,nii-2) * var(or-5,:,:) + &
                      a(nii-4,nii-2) * var(or-4,:,:) + &
                      a(nii-3,nii-2) * var(or-3,:,:) + &
                      a(nii-2,nii-2) * var(or-2,:,:) + &
                      a(nii-1,nii-2) * var(or-1,:,:) + &
                      a(nii,nii-2) * var(or,:,:) ) * idel
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                    ( a(nii-10,nii-1) * var(or-10,:,:) + &
                      a(nii-9,nii-1) * var(or-9,:,:) + &
                      a(nii-8,nii-1) * var(or-8,:,:) + &
                      a(nii-7,nii-1) * var(or-7,:,:) + &
                      a(nii-6,nii-1) * var(or-6,:,:) + &
                      a(nii-5,nii-1) * var(or-5,:,:) + &
                      a(nii-4,nii-1) * var(or-4,:,:) + &
                      a(nii-3,nii-1) * var(or-3,:,:) + &
                      a(nii-2,nii-1) * var(or-2,:,:) + &
                      a(nii-1,nii-1) * var(or-1,:,:) + &
                      a(nii,nii-1) * var(or,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                    ( a(nii-4,nii) * var(or-4,:,:) + &
                      a(nii-3,nii) * var(or-3,:,:) + &
                      a(nii-2,nii) * var(or-2,:,:) + &
                      a(nii-1,nii) * var(or-1,:,:) + &
                      a(nii,nii) * var(or,:,:) ) * idel

    ir = or - 7

  end if

  if ( bb(1) /= 0 ) then
     do i = il, ir
        rhs(i,:,:) = rhs(i,:,:) + &
                       ( a(i-ol-4,i-ol) * var(i-4,:,:) + &
                         a(i-ol-3,i-ol) * var(i-3,:,:) + &
                         a(i-ol-2,i-ol) * var(i-2,:,:) + &
                         a(i-ol-1,i-ol) * var(i-1,:,:) + &
                         a(i-ol,i-ol) * var(i,:,:) + &
                         a(i-ol+1,i-ol) * var(i+1,:,:) + &
                         a(i-ol+2,i-ol) * var(i+2,:,:) + &
                         a(i-ol+3,i-ol) * var(i+3,:,:) + &
                         a(i-ol+4,i-ol) * var(i+4,:,:) ) * idel
     end do

  else

     do i = il, ir
        rhs(i,:,:) = rhs(i,:,:) + &
                       ( a(i-4,i) * var(i-4,:,:) + &
                         a(i-3,i) * var(i-3,:,:) + &
                         a(i-2,i) * var(i-2,:,:) + &
                         a(i-1,i) * var(i-1,:,:) + &
                         a(i,i) * var(i,:,:) + &
                         a(i+1,i) * var(i+1,:,:) + &
                         a(i+2,i) * var(i+2,:,:) + &
                         a(i+3,i) * var(i+3,:,:) + &
                         a(i+4,i) * var(i+4,:,:) ) * idel
     end do

  end if ! bb(1) /= 0


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
  
      call set_bmatrix ( b, bb(3:4), lsh(2), gsh(2)-offset(3)-offset(4), lbnd(2), delta(2), dfl(2) )
    
      call set_hmatrix ( h, bb(3:4) )
  
      call set_dmatrix ( d, bb(3:4) )
  
      atmp = -transpose ( matmul ( h, matmul ( transpose(d), matmul ( b, d ) ) ) )

      allocate ( ycoeff(patch)%coeff(njj,njj) )

      ycoeff(patch)%coeff = atmp
 
      savedy(patch) = .true.

      deallocate ( atmp, d, b, h )

    end if
 
    a => ycoeff(patch)%coeff
 
    idel = epsilon /  ( 256.0_wp * delta(2) )
   
    if ( bb(3) == 0 ) then
  
      jl = 1 + gsize(2)
  
    else
      ol = offset(3)

      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) + &
                     a(5,1) * var(:,5+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) + a(6,2) * var(:,6+ol,:) + &
                     a(7,2) * var(:,7+ol,:) + a(8,2) * var(:,8+ol,:) + &
                     a(9,2) * var(:,9+ol,:) + a(10,2) * var(:,10+ol,:) + &
                     a(11,2) * var(:,11+ol,:) ) * idel
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) + &
                     a(7,3) * var(:,7+ol,:) + a(8,3) * var(:,8+ol,:) + &
                     a(9,3) * var(:,9+ol,:) + a(10,3) * var(:,10+ol,:) + &
                     a(11,3) * var(:,11+ol,:) ) * idel
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) + a(8,4) * var(:,8+ol,:) + &
                     a(9,4) * var(:,9+ol,:) + a(10,4) * var(:,10+ol,:) + &
                     a(11,4) * var(:,11+ol,:) ) * idel
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(1,5) * var(:,1+ol,:) + a(2,5) * var(:,2+ol,:) + &
                     a(3,5) * var(:,3+ol,:) + a(4,5) * var(:,4+ol,:) + &
                     a(5,5) * var(:,5+ol,:) + a(6,5) * var(:,6+ol,:) + &
                     a(7,5) * var(:,7+ol,:) + a(8,5) * var(:,8+ol,:) + &
                     a(9,5) * var(:,9+ol,:) + a(10,5) * var(:,10+ol,:) + &
                     a(11,5) * var(:,11+ol,:) ) * idel
      rhs(:,6+ol,:) = rhs(:,6+ol,:) + &
                   ( a(1,6) * var(:,1+ol,:) + a(2,6) * var(:,2+ol,:) + &
                     a(3,6) * var(:,3+ol,:) + a(4,6) * var(:,4+ol,:) + &
                     a(5,6) * var(:,5+ol,:) + a(6,6) * var(:,6+ol,:) + &
                     a(7,6) * var(:,7+ol,:) + a(8,6) * var(:,8+ol,:) + &
                     a(9,6) * var(:,9+ol,:) + a(10,6) * var(:,10+ol,:) + &
                     a(11,6) * var(:,11+ol,:) ) * idel
      rhs(:,7+ol,:) = rhs(:,7+ol,:) + &
                   ( a(1,7) * var(:,1+ol,:) + a(2,7) * var(:,2+ol,:) + &
                     a(3,7) * var(:,3+ol,:) + a(4,7) * var(:,4+ol,:) + &
                     a(5,7) * var(:,5+ol,:) + a(6,7) * var(:,6+ol,:) + &
                     a(7,7) * var(:,7+ol,:) + a(8,7) * var(:,8+ol,:) + &
                     a(9,7) * var(:,9+ol,:) + a(10,7) * var(:,10+ol,:) + &
                     a(11,7) * var(:,11+ol,:) ) * idel
  
      jl = 8 + ol
  
    end if
  
    if ( bb(4) == 0 ) then
  
      jr = nj - gsize(2)
  
    else
      or = nj - offset(4)
  
      rhs(:,or-6,:) = rhs(:,or-6,:) + &
                      ( a(njj-10,njj-6) * var(:,or-10,:) + &
                        a(njj-9,njj-6) * var(:,or-9,:) + &
                        a(njj-8,njj-6) * var(:,or-8,:) + &
                        a(njj-7,njj-6) * var(:,or-7,:) + &
                        a(njj-6,njj-6) * var(:,or-6,:) + &
                        a(njj-5,njj-6) * var(:,or-5,:) + &
                        a(njj-4,njj-6) * var(:,or-4,:) + &
                        a(njj-3,njj-6) * var(:,or-3,:) + &
                        a(njj-2,njj-6) * var(:,or-2,:) + &
                        a(njj-1,njj-6) * var(:,or-1,:) + &
                        a(njj,njj-6) * var(:,or,:) ) * idel
      rhs(:,or-5,:) = rhs(:,or-5,:) + &
                      ( a(njj-10,njj-5) * var(:,or-10,:) + &
                        a(njj-9,njj-5) * var(:,or-9,:) + &
                        a(njj-8,njj-5) * var(:,or-8,:) + &
                        a(njj-7,njj-5) * var(:,or-7,:) + &
                        a(njj-6,njj-5) * var(:,or-6,:) + &
                        a(njj-5,njj-5) * var(:,or-5,:) + &
                        a(njj-4,njj-5) * var(:,or-4,:) + &
                        a(njj-3,njj-5) * var(:,or-3,:) + &
                        a(njj-2,njj-5) * var(:,or-2,:) + &
                        a(njj-1,njj-5) * var(:,or-1,:) + &
                        a(njj,njj-5) * var(:,or,:) ) * idel
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                      ( a(njj-10,njj-4) * var(:,or-10,:) + &
                        a(njj-9,njj-4) * var(:,or-9,:) + &
                        a(njj-8,njj-4) * var(:,or-8,:) + &
                        a(njj-7,njj-4) * var(:,or-7,:) + &
                        a(njj-6,njj-4) * var(:,or-6,:) + &
                        a(njj-5,njj-4) * var(:,or-5,:) + &
                        a(njj-4,njj-4) * var(:,or-4,:) + &
                        a(njj-3,njj-4) * var(:,or-3,:) + &
                        a(njj-2,njj-4) * var(:,or-2,:) + &
                        a(njj-1,njj-4) * var(:,or-1,:) + &
                        a(njj,njj-4) * var(:,or,:) ) * idel
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                      ( a(njj-10,njj-3) * var(:,or-10,:) + &
                        a(njj-9,njj-3) * var(:,or-9,:) + &
                        a(njj-8,njj-3) * var(:,or-8,:) + &
                        a(njj-7,njj-3) * var(:,or-7,:) + &
                        a(njj-6,njj-3) * var(:,or-6,:) + &
                        a(njj-5,njj-3) * var(:,or-5,:) + &
                        a(njj-4,njj-3) * var(:,or-4,:) + &
                        a(njj-3,njj-3) * var(:,or-3,:) + &
                        a(njj-2,njj-3) * var(:,or-2,:) + &
                        a(njj-1,njj-3) * var(:,or-1,:) + &
                        a(njj,njj-3) * var(:,or,:) ) * idel
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                      ( a(njj-10,njj-2) * var(:,or-10,:) + &
                        a(njj-9,njj-2) * var(:,or-9,:) + &
                        a(njj-8,njj-2) * var(:,or-8,:) + &
                        a(njj-7,njj-2) * var(:,or-7,:) + &
                        a(njj-6,njj-2) * var(:,or-6,:) + &
                        a(njj-5,njj-2) * var(:,or-5,:) + &
                        a(njj-4,njj-2) * var(:,or-4,:) + &
                        a(njj-3,njj-2) * var(:,or-3,:) + &
                        a(njj-2,njj-2) * var(:,or-2,:) + &
                        a(njj-1,njj-2) * var(:,or-1,:) + &
                        a(njj,njj-2) * var(:,or,:) ) * idel
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                      ( a(njj-10,njj-1) * var(:,or-10,:) + &
                        a(njj-9,njj-1) * var(:,or-9,:) + &
                        a(njj-8,njj-1) * var(:,or-8,:) + &
                        a(njj-7,njj-1) * var(:,or-7,:) + &
                        a(njj-6,njj-1) * var(:,or-6,:) + &
                        a(njj-5,njj-1) * var(:,or-5,:) + &
                        a(njj-4,njj-1) * var(:,or-4,:) + &
                        a(njj-3,njj-1) * var(:,or-3,:) + &
                        a(njj-2,njj-1) * var(:,or-2,:) + &
                        a(njj-1,njj-1) * var(:,or-1,:) + &
                        a(njj,njj-1) * var(:,or,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                      ( a(njj-4,njj) * var(:,or-4,:) + &
                        a(njj-3,njj) * var(:,or-3,:) + &
                        a(njj-2,njj) * var(:,or-2,:) + &
                        a(njj-1,njj) * var(:,or-1,:) + &
                        a(njj,njj) * var(:,or,:) ) * idel
  
      jr = or - 7
   
    end if

    do j = jl, jr

      rhs(:,j,:) = rhs(:,j,:) + &
                     ( a(j-ol-4,j-ol) * var(:,j-4,:) + &
                       a(j-ol-3,j-ol) * var(:,j-3,:) + &
                       a(j-ol-2,j-ol) * var(:,j-2,:) + &
                       a(j-ol-1,j-ol) * var(:,j-1,:) + &
                       a(j-ol,j-ol) * var(:,j,:) + &
                       a(j-ol+1,j-ol) * var(:,j+1,:) + &
                       a(j-ol+2,j-ol) * var(:,j+2,:) + &
                       a(j-ol+3,j-ol) * var(:,j+3,:) + &
                       a(j-ol+4,j-ol) * var(:,j+4,:) ) * idel
  
    end do
  
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
  
      call set_bmatrix ( b, bb(5:6), lsh(3), gsh(3)-offset(5)-offset(6), lbnd(3), delta(3), dfl(3) )

      call set_hmatrix ( h, bb(5:6) )

      call set_dmatrix ( d, bb(5:6) )

      atmp = -transpose ( matmul ( h, matmul ( transpose(d), matmul ( b, d ) ) ) )

      allocate ( zcoeff(patch)%coeff(nkk,nkk) )

      zcoeff(patch)%coeff = atmp

      savedz(patch) = .true.

      deallocate ( atmp, d, b, h )

    end if
  
    a => zcoeff(patch)%coeff

    idel = epsilon /  ( 256.0_wp * delta(3) )
  
    if ( bb(5) == 0 ) then
  
      kl = 1 + gsize(3)
  
    else
      ol = offset(5)
  
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) + &
                     a(5,1) * var(:,:,5+ol) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) + a(6,2) * var(:,:,6+ol) + &
                     a(7,2) * var(:,:,7+ol) + a(8,2) * var(:,:,8+ol) + &
                     a(9,2) * var(:,:,9+ol) + a(10,2) * var(:,:,10+ol) + &
                     a(11,2) * var(:,:,11+ol) ) * idel
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) + &
                     a(7,3) * var(:,:,7+ol) + a(8,3) * var(:,:,8+ol) + &
                     a(9,3) * var(:,:,9+ol) + a(10,3) * var(:,:,10+ol) + &
                     a(11,3) * var(:,:,11+ol) ) * idel
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) + a(8,4) * var(:,:,8+ol) + &
                     a(9,4) * var(:,:,9+ol) + a(10,4) * var(:,:,10+ol) + &
                     a(11,4) * var(:,:,11+ol) ) * idel
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(1,5) * var(:,:,1+ol) + a(2,5) * var(:,:,2+ol) + &
                     a(3,5) * var(:,:,3+ol) + a(4,5) * var(:,:,4+ol) + &
                     a(5,5) * var(:,:,5+ol) + a(6,5) * var(:,:,6+ol) + &
                     a(7,5) * var(:,:,7+ol) + a(8,5) * var(:,:,8+ol) + &
                     a(9,5) * var(:,:,9+ol) + a(10,5) * var(:,:,10+ol) + &
                     a(11,5) * var(:,:,11+ol) ) * idel
      rhs(:,:,6+ol) = rhs(:,:,6+ol) + &
                   ( a(1,6) * var(:,:,1+ol) + a(2,6) * var(:,:,2+ol) + &
                     a(3,6) * var(:,:,3+ol) + a(4,6) * var(:,:,4+ol) + &
                     a(5,6) * var(:,:,5+ol) + a(6,6) * var(:,:,6+ol) + &
                     a(7,6) * var(:,:,7+ol) + a(8,6) * var(:,:,8+ol) + &
                     a(9,6) * var(:,:,9+ol) + a(10,6) * var(:,:,10+ol) + &
                     a(11,6) * var(:,:,11+ol) ) * idel
      rhs(:,:,7+ol) = rhs(:,:,7+ol) + &
                   ( a(1,7) * var(:,:,1+ol) + a(2,7) * var(:,:,2+ol) + &
                     a(3,7) * var(:,:,3+ol) + a(4,7) * var(:,:,4+ol) + &
                     a(5,7) * var(:,:,5+ol) + a(6,7) * var(:,:,6+ol) + &
                     a(7,7) * var(:,:,7+ol) + a(8,7) * var(:,:,8+ol) + &
                     a(9,7) * var(:,:,9+ol) + a(10,7) * var(:,:,10+ol) + &
                     a(11,7) * var(:,:,11+ol) ) * idel
  
      kl = 8 + ol
  
    end if
  
    if ( bb(6) == 0 ) then
  
      kr = nk - gsize(3)
  
    else
      or = nk - offset(6)
  
      rhs(:,:,or-6) = rhs(:,:,or-6) + &
                      ( a(nkk-10,nkk-6) * var(:,:,or-10) + &
                        a(nkk-9,nkk-6) * var(:,:,or-9) + &
                        a(nkk-8,nkk-6) * var(:,:,or-8) + &
                        a(nkk-7,nkk-6) * var(:,:,or-7) + &
                        a(nkk-6,nkk-6) * var(:,:,or-6) + &
                        a(nkk-5,nkk-6) * var(:,:,or-5) + &
                        a(nkk-4,nkk-6) * var(:,:,or-4) + &
                        a(nkk-3,nkk-6) * var(:,:,or-3) + &
                        a(nkk-2,nkk-6) * var(:,:,or-2) + &
                        a(nkk-1,nkk-6) * var(:,:,or-1) + &
                        a(nkk,nkk-6) * var(:,:,or) ) * idel
      rhs(:,:,or-5) = rhs(:,:,or-5) + &
                      ( a(nkk-10,nkk-5) * var(:,:,or-10) + &
                        a(nkk-9,nkk-5) * var(:,:,or-9) + &
                        a(nkk-8,nkk-5) * var(:,:,or-8) + &
                        a(nkk-7,nkk-5) * var(:,:,or-7) + &
                        a(nkk-6,nkk-5) * var(:,:,or-6) + &
                        a(nkk-5,nkk-5) * var(:,:,or-5) + &
                        a(nkk-4,nkk-5) * var(:,:,or-4) + &
                        a(nkk-3,nkk-5) * var(:,:,or-3) + &
                        a(nkk-2,nkk-5) * var(:,:,or-2) + &
                        a(nkk-1,nkk-5) * var(:,:,or-1) + &
                        a(nkk,nkk-5) * var(:,:,or) ) * idel
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                      ( a(nkk-10,nkk-4) * var(:,:,or-10) + &
                        a(nkk-9,nkk-4) * var(:,:,or-9) + &
                        a(nkk-8,nkk-4) * var(:,:,or-8) + &
                        a(nkk-7,nkk-4) * var(:,:,or-7) + &
                        a(nkk-6,nkk-4) * var(:,:,or-6) + &
                        a(nkk-5,nkk-4) * var(:,:,or-5) + &
                        a(nkk-4,nkk-4) * var(:,:,or-4) + &
                        a(nkk-3,nkk-4) * var(:,:,or-3) + &
                        a(nkk-2,nkk-4) * var(:,:,or-2) + &
                        a(nkk-1,nkk-4) * var(:,:,or-1) + &
                        a(nkk,nkk-4) * var(:,:,or) ) * idel
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                      ( a(nkk-10,nkk-3) * var(:,:,or-10) + &
                        a(nkk-9,nkk-3) * var(:,:,or-9) + &
                        a(nkk-8,nkk-3) * var(:,:,or-8) + &
                        a(nkk-7,nkk-3) * var(:,:,or-7) + &
                        a(nkk-6,nkk-3) * var(:,:,or-6) + &
                        a(nkk-5,nkk-3) * var(:,:,or-5) + &
                        a(nkk-4,nkk-3) * var(:,:,or-4) + &
                        a(nkk-3,nkk-3) * var(:,:,or-3) + &
                        a(nkk-2,nkk-3) * var(:,:,or-2) + &
                        a(nkk-1,nkk-3) * var(:,:,or-1) + &
                        a(nkk,nkk-3) * var(:,:,or) ) * idel
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                      ( a(nkk-10,nkk-2) * var(:,:,or-10) + &
                        a(nkk-9,nkk-2) * var(:,:,or-9) + &
                        a(nkk-8,nkk-2) * var(:,:,or-8) + &
                        a(nkk-7,nkk-2) * var(:,:,or-7) + &
                        a(nkk-6,nkk-2) * var(:,:,or-6) + &
                        a(nkk-5,nkk-2) * var(:,:,or-5) + &
                        a(nkk-4,nkk-2) * var(:,:,or-4) + &
                        a(nkk-3,nkk-2) * var(:,:,or-3) + &
                        a(nkk-2,nkk-2) * var(:,:,or-2) + &
                        a(nkk-1,nkk-2) * var(:,:,or-1) + &
                        a(nkk,nkk-2) * var(:,:,or) ) * idel
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                      ( a(nkk-10,nkk-1) * var(:,:,or-10) + &
                        a(nkk-9,nkk-1) * var(:,:,or-9) + &
                        a(nkk-8,nkk-1) * var(:,:,or-8) + &
                        a(nkk-7,nkk-1) * var(:,:,or-7) + &
                        a(nkk-6,nkk-1) * var(:,:,or-6) + &
                        a(nkk-5,nkk-1) * var(:,:,or-5) + &
                        a(nkk-4,nkk-1) * var(:,:,or-4) + &
                        a(nkk-3,nkk-1) * var(:,:,or-3) + &
                        a(nkk-2,nkk-1) * var(:,:,or-2) + &
                        a(nkk-1,nkk-1) * var(:,:,or-1) + &
                        a(nkk,nkk-1) * var(:,:,or) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                      ( a(nkk-4,nkk) * var(:,:,or-4) + &
                        a(nkk-3,nkk) * var(:,:,or-3) + &
                        a(nkk-2,nkk) * var(:,:,or-2) + &
                        a(nkk-1,nkk) * var(:,:,or-1) + &
                        a(nkk,nkk) * var(:,:,or) ) * idel
  
      kr = or - 7
  
    end if
  
    do k = kl, kr
  
      if ( bb(5) /= 0) then
         rhs(:,:,k) = rhs(:,:,k) + &
                     ( a(k-ol-4,k-ol) * var(:,:,k-4) + &
                       a(k-ol-3,k-ol) * var(:,:,k-3) + &
                       a(k-ol-2,k-ol) * var(:,:,k-2) + &
                       a(k-ol-1,k-ol) * var(:,:,k-1) + &
                       a(k-ol,k-ol) * var(:,:,k) + &
                       a(k-ol+1,k-ol) * var(:,:,k+1) + &
                       a(k-ol+2,k-ol) * var(:,:,k+2) + &
                       a(k-ol+3,k-ol) * var(:,:,k+3) + &
                       a(k-ol+4,k-ol) * var(:,:,k+4) ) * idel
      else
         rhs(:,:,k) = rhs(:,:,k) + &
                     ( a(k-4,k) * var(:,:,k-4) + &
                       a(k-3,k) * var(:,:,k-3) + &
                       a(k-2,k) * var(:,:,k-2) + &
                       a(k-1,k) * var(:,:,k-1) + &
                       a(k,k) * var(:,:,k) + &
                       a(k+1,k) * var(:,:,k+1) + &
                       a(k+2,k) * var(:,:,k+2) + &
                       a(k+3,k) * var(:,:,k+3) + &
                       a(k+4,k) * var(:,:,k+4) ) * idel
      endif
  
    end do
  
  end if

contains

  subroutine set_dmatrix ( d, bb )

    implicit none

    CCTK_REAL, dimension(:,:), intent(out) :: d
    CCTK_INT, dimension(2), intent(in) :: bb
    CCTK_INT :: n
    CCTK_REAL, dimension(5), save :: ac = (/ 1.0_wp, -4.0_wp, 6.0_wp, &
                                                     -4.0_wp, 1.0_wp /)
    CCTK_INT :: i

    d = zero

    n = size(d,1)

    if ( bb(1) == 0 ) then
      d(1,1:3) = ac(3:5)
      d(2,1:4) = ac(2:5)
    else
      d(1,1:5) = ac
      d(2,1:5) = ac
    end if
    if ( bb(2) == 0 ) then
      d(n-1,n-3:n) = ac(1:4)
      d(n,n-2:n) = ac(1:3)
    else
      d(n-1,n-4:n) = ac
      d(n,n-4:n) = ac
    end if
    do i = 3, n-2
      d(i,i-2:i+2) = ac
    end do
  end subroutine set_dmatrix


  subroutine set_bmatrix ( b, bb, lsh, gsh, lbnd, h, dfl )

    implicit none

    CCTK_REAL, dimension(:,:), intent(inout) :: b
    CCTK_INT, dimension(2), intent(in) :: bb
    CCTK_INT, intent(in) :: lsh, gsh, lbnd
    CCTK_REAL, intent(in) :: h, dfl
    CCTK_INT :: n
    CCTK_INT :: i
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

    h = zero

    do i = 1, n
      h(i,i) = 1.0_wp
    end do

    if ( bb(1) /= 0 ) then
      h(1:7,1:7) = sigma ( )
    end if

    if ( bb(2) /= 0 ) then
      h(n:n-6:-1,n:n-6:-1) = sigma ( )
    end if

  end subroutine set_hmatrix
    

  function bf ( x, xp, xmax, h )

    implicit none

    CCTK_REAL :: bf
    CCTK_REAL, intent(in) :: x, xp, xmax, h
    CCTK_REAL :: h2, h21, xp2, xp3

    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    if ( x < xp ) then

      h2 = h**2
      h21 = 1.0_wp - h2
      xp2 = 1.0_wp / xp**2
      xp3 = xp2 / xp
      bf = ( -2.0_wp * h21 * xp3 * x + 3.0_wp * h21 * xp2 ) * x**2 + h2

    else if ( x > xmax - xp ) then

      h2 = h**2
      h21 = 1.0_wp - h2
      xp2 = 1.0_wp / xp**2
      xp3 = xp2 / xp
      bf = ( 2.0_wp * h21 * xp3 * ( x - xmax ) + &
             3.0_wp * h21 * xp2 ) * ( x - xmax )**2 + h2

    else

      bf = 1.0_wp

    end if

  end function bf


  function sigma ( )

    implicit none

    CCTK_REAL, dimension(7,7) :: sigma

    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    sigma(1,1) = 4.930709842221048047321555312222552006914_wp
    sigma(2,1) = 0.0_wp
    sigma(3,1) = 0.0_wp
    sigma(4,1) = 0.0_wp
    sigma(5,1) = 0.0_wp
    sigma(6,1) = 0.0_wp
    sigma(7,1) = 0.0_wp
    sigma(1,2) = 0.0_wp
    sigma(2,2) = 1.499128158967098993672808967791430710830_wp
    sigma(3,2) = 0.8363320637711490458766332755492551638364_wp
    sigma(4,2) = -0.2634194189384563273844251938374368097488_wp
    sigma(5,2) = -0.1281142262970477369698649073120970227895_wp
    sigma(6,2) = -0.1830121923043950334684235233648902467425_wp
    sigma(7,2) = 0.01275879120879511857581469566228951245742_wp
    sigma(1,3) = 0.0_wp
    sigma(2,3) = 0.8363320637711490458766332755492551638364_wp
    sigma(3,3) = 0.8785380112113613457715293301885586949116_wp
    sigma(4,3) = 0.1454039706555408804118493444561963110723_wp
    sigma(5,3) = -0.03346859376328961536059718572568704595897_wp
    sigma(6,3) = -0.1759548942026651677201193466207997400242_wp
    sigma(7,3) = 0.02020404709761823823384626194775259449799_wp
    sigma(1,4) = 0.0_wp
    sigma(2,4) = -0.2634194189384563273844251938374368097488_wp
    sigma(3,4) = 0.1454039706555408804118493444561963110723_wp
    sigma(4,4) = 0.5690707861055695191099807759535994201984_wp
    sigma(5,4) = 0.3374169937858825026860519890688362591356_wp
    sigma(6,4) = -0.03076801586199402080957339555604653302670_wp
    sigma(7,4) = 0.004077760901760668929595481886645691229036_wp
    sigma(1,5) = 0.0_wp
    sigma(2,5) = -0.1281142262970477369698649073120970227895_wp
    sigma(3,5) = -0.03346859376328961536059718572568704595897_wp
    sigma(4,5) = 0.3374169937858825026860519890688362591356_wp
    sigma(5,5) = 0.5001980842834785426891480539161939116986_wp
    sigma(6,5) = 0.2904403943485980286998891564955850119600_wp
    sigma(7,5) = -0.05655813410599246578174268632557910837116_wp
    sigma(1,6) = 0.0_wp
    sigma(2,6) = -0.1830121923043950334684235233648902467425_wp
    sigma(3,6) = -0.1759548942026651677201193466207997400242_wp
    sigma(4,6) = -0.03076801586199402080957339555604653302670_wp
    sigma(5,6) = 0.2904403943485980286998891564955850119600_wp
    sigma(6,6) = 0.8640467497799918355533313312208192675134_wp
    sigma(7,6) = 0.03704582714017916999084939064121697271831_wp
    sigma(1,7) = 0.0_wp
    sigma(2,7) = 0.01275879120879511857581469566228951245742_wp
    sigma(3,7) = 0.02020404709761823823384626194775259449799_wp
    sigma(4,7) = 0.004077760901760668929595481886645691229036_wp
    sigma(5,7) = -0.05655813410599246578174268632557910837116_wp
    sigma(6,7) = 0.03704582714017916999084939064121697271831_wp
    sigma(7,7) = 0.9909401061257062767072573096147576992963_wp

  end function sigma

end subroutine dissipation_6_5_alt
