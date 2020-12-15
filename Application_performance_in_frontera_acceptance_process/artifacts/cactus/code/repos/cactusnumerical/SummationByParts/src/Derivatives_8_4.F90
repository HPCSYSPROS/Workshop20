#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine deriv_gf_8_4 ( var, ni, nj, nk, dir, bb, gsize, offset, delta, dvar )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(IN) :: var
  CCTK_INT, intent(IN) :: dir
  CCTK_INT, intent(IN) :: bb(2)
  CCTK_INT, intent(IN) :: gsize
  CCTK_INT, intent(IN) :: offset(2)
  CCTK_REAL, intent(IN) :: delta
  CCTK_REAL, dimension(ni,nj,nk), intent(OUT) :: dvar

  CCTK_REAL, dimension(4), save :: a
  CCTK_REAL, dimension(12,8), save :: q 
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_8_4 ( a, q )
    first = .false. 
  end if

  idel = 1.0_wp / delta

  if (gsize < 4) call CCTK_WARN (0, "not enough ghostzones")

  direction: select case (dir)
  case (0) direction
    if ( bb(1) == 0 ) then
      il = 1 + gsize
    else
      ol = offset(1)
      dvar(1+ol,:,:) = ( q(1,1) * var(1+ol,:,:) + q(2,1) * var(2+ol,:,:) + &
                         q(3,1) * var(3+ol,:,:) + q(4,1) * var(4+ol,:,:) + &
                         q(5,1) * var(5+ol,:,:) + q(6,1) * var(6+ol,:,:) + &
                         q(7,1) * var(7+ol,:,:) + q(8,1) * var(8+ol,:,:) ) * idel
      dvar(2+ol,:,:) = ( q(1,2) * var(1+ol,:,:) + q(3,2) * var(3+ol,:,:) + &
                         q(4,2) * var(4+ol,:,:) + q(5,2) * var(5+ol,:,:) + &
                         q(6,2) * var(6+ol,:,:) + q(7,2) * var(7+ol,:,:) + &
                         q(8,2) * var(8+ol,:,:) ) * idel
      dvar(3+ol,:,:) = ( q(1,3) * var(1+ol,:,:) + q(2,3) * var(2+ol,:,:) + &
                         q(4,3) * var(4+ol,:,:) + q(5,3) * var(5+ol,:,:) + &
                         q(6,3) * var(6+ol,:,:) + q(7,3) * var(7+ol,:,:) + &
                         q(8,3) * var(8+ol,:,:) ) * idel
      dvar(4+ol,:,:) = ( q(1,4) * var(1+ol,:,:) + q(2,4) * var(2+ol,:,:) + &
                         q(3,4) * var(3+ol,:,:) + q(5,4) * var(5+ol,:,:) + &
                         q(6,4) * var(6+ol,:,:) + q(7,4) * var(7+ol,:,:) + &
                         q(8,4) * var(8+ol,:,:) ) * idel
      dvar(5+ol,:,:) = ( q(1,5) * var(1+ol,:,:) + q(2,5) * var(2+ol,:,:) + &
                         q(3,5) * var(3+ol,:,:) + q(4,5) * var(4+ol,:,:) + &
                         q(6,5) * var(6+ol,:,:) + q(7,5) * var(7+ol,:,:) + &
                         q(8,5) * var(8+ol,:,:) + q(9,5) * var(9+ol,:,:) ) * idel
      dvar(6+ol,:,:) = ( q(1,6) * var(1+ol,:,:) + q(2,6) * var(2+ol,:,:) + &
                         q(3,6) * var(3+ol,:,:) + q(4,6) * var(4+ol,:,:) + &
                         q(5,6) * var(5+ol,:,:) + q(7,6) * var(7+ol,:,:) + &
                         q(8,6) * var(8+ol,:,:) + q(9,6) * var(9+ol,:,:) + &
                         q(10,6) * var(10+ol,:,:) ) * idel
      dvar(7+ol,:,:) = ( q(1,7) * var(1+ol,:,:) + q(2,7) * var(2+ol,:,:) + &
                         q(3,7) * var(3+ol,:,:) + q(4,7) * var(4+ol,:,:) + &
                         q(5,7) * var(5+ol,:,:) + q(6,7) * var(6+ol,:,:) + &
                         q(8,7) * var(8+ol,:,:) + q(9,7) * var(9+ol,:,:) + &
                         q(10,7) * var(10+ol,:,:) + q(11,7) * var(11+ol,:,:) ) * idel
      dvar(8+ol,:,:) = ( q(1,8) * var(1+ol,:,:) + q(2,8) * var(2+ol,:,:) + &
                         q(3,8) * var(3+ol,:,:) + q(4,8) * var(4+ol,:,:) + &
                         q(5,8) * var(5+ol,:,:) + q(6,8) * var(6+ol,:,:) + &
                         q(7,8) * var(7+ol,:,:) + q(9,8) * var(9+ol,:,:) + &
                         q(10,8) * var(10+ol,:,:) + q(11,8) * var(11+ol,:,:) + &
                         q(12,8) * var(12+ol,:,:) ) * idel
      il = 9 + ol
    end if
    if ( bb(2) == 0 ) then
      ir = ni - gsize
    else
      or = ni - offset(2)
      dvar(or,:,:) =   - ( q(1,1) * var(or,:,:) + q(2,1) * var(or-1,:,:) + &
                           q(3,1) * var(or-2,:,:) + q(4,1) * var(or-3,:,:) + &
                           q(5,1) * var(or-4,:,:) + q(6,1) * var(or-5,:,:) + &
                           q(7,1) * var(or-6,:,:) + &
                           q(8,1) * var(or-7,:,:) ) * idel
      dvar(or-1,:,:) = - ( q(1,2) * var(or,:,:) + q(3,2) * var(or-2,:,:) + &
                           q(4,2) * var(or-3,:,:) + q(5,2) * var(or-4,:,:) + &
                           q(6,2) * var(or-5,:,:) + q(7,2) * var(or-6,:,:) + &
                           q(8,2) * var(or-7,:,:) ) * idel
      dvar(or-2,:,:) = - ( q(1,3) * var(or,:,:) + q(2,3) * var(or-1,:,:) + &
                           q(4,3) * var(or-3,:,:) + q(5,3) * var(or-4,:,:) + &
                           q(6,3) * var(or-5,:,:) + q(7,3) * var(or-6,:,:) + &
                           q(8,3) * var(or-7,:,:) ) * idel
      dvar(or-3,:,:) = - ( q(1,4) * var(or,:,:) + q(2,4) * var(or-1,:,:) + &
                           q(3,4) * var(or-2,:,:) + q(5,4) * var(or-4,:,:) + &
                           q(6,4) * var(or-5,:,:) + q(7,4) * var(or-6,:,:) + &
                           q(8,4) * var(or-7,:,:) ) * idel
      dvar(or-4,:,:) = - ( q(1,5) * var(or,:,:) + q(2,5) * var(or-1,:,:) + &
                           q(3,5) * var(or-2,:,:) + q(4,5) * var(or-3,:,:) + &
                           q(6,5) * var(or-5,:,:) + q(7,5) * var(or-6,:,:) + &
                           q(8,5) * var(or-7,:,:) + &
                           q(9,5) * var(or-8,:,:) ) * idel
      dvar(or-5,:,:) = - ( q(1,6) * var(or,:,:) + q(2,6) * var(or-1,:,:) + &
                           q(3,6) * var(or-2,:,:) + q(4,6) * var(or-3,:,:) + &
                           q(5,6) * var(or-4,:,:) + q(7,6) * var(or-6,:,:) + &
                           q(8,6) * var(or-7,:,:) + q(9,6) * var(or-8,:,:) + &
                           q(10,6) * var(or-9,:,:) ) * idel
      dvar(or-6,:,:) = - ( q(1,7) * var(or,:,:) + q(2,7) * var(or-1,:,:) + &
                           q(3,7) * var(or-2,:,:) + q(4,7) * var(or-3,:,:) + &
                           q(5,7) * var(or-4,:,:) + q(6,7) * var(or-5,:,:) + &
                           q(8,7) * var(or-7,:,:) + q(9,7) * var(or-8,:,:) + &
                           q(10,7) * var(or-9,:,:) + &
                           q(11,7) * var(or-10,:,:) ) * idel
      dvar(or-7,:,:) = - ( q(1,8) * var(or,:,:) + q(2,8) * var(or-1,:,:) + &
                           q(3,8) * var(or-2,:,:) + q(4,8) * var(or-3,:,:) + &
                           q(5,8) * var(or-4,:,:) + q(6,8) * var(or-5,:,:) + &
                           q(7,8) * var(or-6,:,:) + q(9,8) * var(or-8,:,:) + &
                           q(10,8) * var(or-9,:,:) + &
                           q(11,8) * var(or-10,:,:) + &
                           q(12,8) * var(or-11,:,:) ) * idel
      ir = or - 8
    end if
    if (il > ir+1) call CCTK_WARN (0, "domain too small")
    dvar(il:ir,:,:) = ( a(1) * ( var(il+1:ir+1,:,:) - &
                                 var(il-1:ir-1,:,:) ) + &
                        a(2) * ( var(il+2:ir+2,:,:) - &
                                 var(il-2:ir-2,:,:) ) + &
                        a(3) * ( var(il+3:ir+3,:,:) - &
                                 var(il-3:ir-3,:,:) ) + &
                        a(4) * ( var(il+4:ir+4,:,:) - &
                                 var(il-4:ir-4,:,:) ) ) * idel
  case (1) direction
    if ( zero_derivs_y /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        jl = 1 + gsize
      else
        ol = offset(1)
        dvar(:,1+ol,:) = ( q(1,1) * var(:,1+ol,:) + q(2,1) * var(:,2+ol,:) + &
                           q(3,1) * var(:,3+ol,:) + q(4,1) * var(:,4+ol,:) + &
                           q(5,1) * var(:,5+ol,:) + q(6,1) * var(:,6+ol,:) + &
                           q(7,1) * var(:,7+ol,:) + q(8,1) * var(:,8+ol,:) ) * idel
        dvar(:,2+ol,:) = ( q(1,2) * var(:,1+ol,:) + q(3,2) * var(:,3+ol,:) + &
                           q(4,2) * var(:,4+ol,:) + q(5,2) * var(:,5+ol,:) + &
                           q(6,2) * var(:,6+ol,:) + q(7,2) * var(:,7+ol,:) + &
                           q(8,2) * var(:,8+ol,:) ) * idel
        dvar(:,3+ol,:) = ( q(1,3) * var(:,1+ol,:) + q(2,3) * var(:,2+ol,:) + &
                           q(4,3) * var(:,4+ol,:) + q(5,3) * var(:,5+ol,:) + &
                           q(6,3) * var(:,6+ol,:) + q(7,3) * var(:,7+ol,:) + &
                           q(8,3) * var(:,8+ol,:) ) * idel
        dvar(:,4+ol,:) = ( q(1,4) * var(:,1+ol,:) + q(2,4) * var(:,2+ol,:) + &
                           q(3,4) * var(:,3+ol,:) + q(5,4) * var(:,5+ol,:) + &
                           q(6,4) * var(:,6+ol,:) + q(7,4) * var(:,7+ol,:) + &
                           q(8,4) * var(:,8+ol,:) ) * idel
        dvar(:,5+ol,:) = ( q(1,5) * var(:,1+ol,:) + q(2,5) * var(:,2+ol,:) + &
                           q(3,5) * var(:,3+ol,:) + q(4,5) * var(:,4+ol,:) + &
                           q(6,5) * var(:,6+ol,:) + q(7,5) * var(:,7+ol,:) + &
                           q(8,5) * var(:,8+ol,:) + q(9,5) * var(:,9+ol,:) ) * idel
        dvar(:,6+ol,:) = ( q(1,6) * var(:,1+ol,:) + q(2,6) * var(:,2+ol,:) + &
                           q(3,6) * var(:,3+ol,:) + q(4,6) * var(:,4+ol,:) + &
                           q(5,6) * var(:,5+ol,:) + q(7,6) * var(:,7+ol,:) + &
                           q(8,6) * var(:,8+ol,:) + q(9,6) * var(:,9+ol,:) + &
                           q(10,6) * var(:,10+ol,:) ) * idel
        dvar(:,7+ol,:) = ( q(1,7) * var(:,1+ol,:) + q(2,7) * var(:,2+ol,:) + &
                           q(3,7) * var(:,3+ol,:) + q(4,7) * var(:,4+ol,:) + &
                           q(5,7) * var(:,5+ol,:) + q(6,7) * var(:,6+ol,:) + &
                           q(8,7) * var(:,8+ol,:) + q(9,7) * var(:,9+ol,:) + &
                           q(10,7) * var(:,10+ol,:) + q(11,7) * var(:,11+ol,:) ) * idel
        dvar(:,8+ol,:) = ( q(1,8) * var(:,1+ol,:) + q(2,8) * var(:,2+ol,:) + &
                           q(3,8) * var(:,3+ol,:) + q(4,8) * var(:,4+ol,:) + &
                           q(5,8) * var(:,5+ol,:) + q(6,8) * var(:,6+ol,:) + &
                           q(7,8) * var(:,7+ol,:) + q(9,8) * var(:,9+ol,:) + &
                           q(10,8) * var(:,10+ol,:) + q(11,8) * var(:,11+ol,:) + &
                           q(12,8) * var(:,12+ol,:) ) * idel
        jl = 9 + ol
      end if
      if ( bb(2) == 0 ) then
        jr = nj - gsize
      else
        or = nj - offset(2)
        dvar(:,or,:) =   - ( q(1,1) * var(:,or,:) + q(2,1) * var(:,or-1,:) + &
                             q(3,1) * var(:,or-2,:) + q(4,1) * var(:,or-3,:) + &
                             q(5,1) * var(:,or-4,:) + q(6,1) * var(:,or-5,:) + &
                             q(7,1) * var(:,or-6,:) + &
                             q(8,1) * var(:,or-7,:) ) * idel
        dvar(:,or-1,:) = - ( q(1,2) * var(:,or,:) + q(3,2) * var(:,or-2,:) + &
                             q(4,2) * var(:,or-3,:) + q(5,2) * var(:,or-4,:) + &
                             q(6,2) * var(:,or-5,:) + q(7,2) * var(:,or-6,:) + &
                             q(8,2) * var(:,or-7,:) ) * idel
        dvar(:,or-2,:) = - ( q(1,3) * var(:,or,:) + q(2,3) * var(:,or-1,:) + &
                             q(4,3) * var(:,or-3,:) + q(5,3) * var(:,or-4,:) + &
                             q(6,3) * var(:,or-5,:) + q(7,3) * var(:,or-6,:) + &
                             q(8,3) * var(:,or-7,:) ) * idel
        dvar(:,or-3,:) = - ( q(1,4) * var(:,or,:) + q(2,4) * var(:,or-1,:) + &
                             q(3,4) * var(:,or-2,:) + q(5,4) * var(:,or-4,:) + &
                             q(6,4) * var(:,or-5,:) + q(7,4) * var(:,or-6,:) + &
                             q(8,4) * var(:,or-7,:) ) * idel
        dvar(:,or-4,:) = - ( q(1,5) * var(:,or,:) + q(2,5) * var(:,or-1,:) + &
                             q(3,5) * var(:,or-2,:) + q(4,5) * var(:,or-3,:) + &
                             q(6,5) * var(:,or-5,:) + q(7,5) * var(:,or-6,:) + &
                             q(8,5) * var(:,or-7,:) + &
                             q(9,5) * var(:,or-8,:) ) * idel
        dvar(:,or-5,:) = - ( q(1,6) * var(:,or,:) + q(2,6) * var(:,or-1,:) + &
                             q(3,6) * var(:,or-2,:) + q(4,6) * var(:,or-3,:) + &
                             q(5,6) * var(:,or-4,:) + q(7,6) * var(:,or-6,:) + &
                             q(8,6) * var(:,or-7,:) + q(9,6) * var(:,or-8,:) + &
                             q(10,6) * var(:,or-9,:) ) * idel
        dvar(:,or-6,:) = - ( q(1,7) * var(:,or,:) + q(2,7) * var(:,or-1,:) + &
                             q(3,7) * var(:,or-2,:) + q(4,7) * var(:,or-3,:) + &
                             q(5,7) * var(:,or-4,:) + q(6,7) * var(:,or-5,:) + &
                             q(8,7) * var(:,or-7,:) + q(9,7) * var(:,or-8,:) + &
                             q(10,7) * var(:,or-9,:) + &
                             q(11,7) * var(:,or-10,:) ) * idel
        dvar(:,or-7,:) = - ( q(1,8) * var(:,or,:) + q(2,8) * var(:,or-1,:) + &
                             q(3,8) * var(:,or-2,:) + q(4,8) * var(:,or-3,:) + &
                             q(5,8) * var(:,or-4,:) + q(6,8) * var(:,or-5,:) + &
                             q(7,8) * var(:,or-6,:) + q(9,8) * var(:,or-8,:) + &
                             q(10,8) * var(:,or-9,:) + &
                             q(11,8) * var(:,or-10,:) + &
                             q(12,8) * var(:,or-11,:) ) * idel
        jr = or - 8
      end if
      if (jl > jr+1) call CCTK_WARN (0, "domain too small")
      dvar(:,jl:jr,:) = ( a(1) * ( var(:,jl+1:jr+1,:) - &
                                   var(:,jl-1:jr-1,:) ) + &
                          a(2) * ( var(:,jl+2:jr+2,:) - &
                                   var(:,jl-2:jr-2,:) ) + &
                          a(3) * ( var(:,jl+3:jr+3,:) - &
                                   var(:,jl-3:jr-3,:) ) + &
                          a(4) * ( var(:,jl+4:jr+4,:) - &
                                   var(:,jl-4:jr-4,:) ) ) * idel
    end if
  case (2) direction
    if ( zero_derivs_z /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        kl = 1 + gsize
      else
        ol = offset(1)
        dvar(:,:,1+ol) = ( q(1,1) * var(:,:,1+ol) + q(2,1) * var(:,:,2+ol) + &
                           q(3,1) * var(:,:,3+ol) + q(4,1) * var(:,:,4+ol) + &
                           q(5,1) * var(:,:,5+ol) + q(6,1) * var(:,:,6+ol) + &
                           q(7,1) * var(:,:,7+ol) + q(8,1) * var(:,:,8+ol) ) * idel
        dvar(:,:,2+ol) = ( q(1,2) * var(:,:,1+ol) + q(3,2) * var(:,:,3+ol) + &
                           q(4,2) * var(:,:,4+ol) + q(5,2) * var(:,:,5+ol) + &
                           q(6,2) * var(:,:,6+ol) + q(7,2) * var(:,:,7+ol) + &
                           q(8,2) * var(:,:,8+ol) ) * idel
        dvar(:,:,3+ol) = ( q(1,3) * var(:,:,1+ol) + q(2,3) * var(:,:,2+ol) + &
                           q(4,3) * var(:,:,4+ol) + q(5,3) * var(:,:,5+ol) + &
                           q(6,3) * var(:,:,6+ol) + q(7,3) * var(:,:,7+ol) + &
                           q(8,3) * var(:,:,8+ol) ) * idel
        dvar(:,:,4+ol) = ( q(1,4) * var(:,:,1+ol) + q(2,4) * var(:,:,2+ol) + &
                           q(3,4) * var(:,:,3+ol) + q(5,4) * var(:,:,5+ol) + &
                           q(6,4) * var(:,:,6+ol) + q(7,4) * var(:,:,7+ol) + &
                           q(8,4) * var(:,:,8+ol) ) * idel
        dvar(:,:,5+ol) = ( q(1,5) * var(:,:,1+ol) + q(2,5) * var(:,:,2+ol) + &
                           q(3,5) * var(:,:,3+ol) + q(4,5) * var(:,:,4+ol) + &
                           q(6,5) * var(:,:,6+ol) + q(7,5) * var(:,:,7+ol) + &
                           q(8,5) * var(:,:,8+ol) + q(9,5) * var(:,:,9+ol) ) * idel
        dvar(:,:,6+ol) = ( q(1,6) * var(:,:,1+ol) + q(2,6) * var(:,:,2+ol) + &
                           q(3,6) * var(:,:,3+ol) + q(4,6) * var(:,:,4+ol) + &
                           q(5,6) * var(:,:,5+ol) + q(7,6) * var(:,:,7+ol) + &
                           q(8,6) * var(:,:,8+ol) + q(9,6) * var(:,:,9+ol) + &
                           q(10,6) * var(:,:,10+ol) ) * idel
        dvar(:,:,7+ol) = ( q(1,7) * var(:,:,1+ol) + q(2,7) * var(:,:,2+ol) + &
                           q(3,7) * var(:,:,3+ol) + q(4,7) * var(:,:,4+ol) + &
                           q(5,7) * var(:,:,5+ol) + q(6,7) * var(:,:,6+ol) + &
                           q(8,7) * var(:,:,8+ol) + q(9,7) * var(:,:,9+ol) + &
                           q(10,7) * var(:,:,10+ol) + q(11,7) * var(:,:,11+ol) ) * idel
        dvar(:,:,8+ol) = ( q(1,8) * var(:,:,1+ol) + q(2,8) * var(:,:,2+ol) + &
                           q(3,8) * var(:,:,3+ol) + q(4,8) * var(:,:,4+ol) + &
                           q(5,8) * var(:,:,5+ol) + q(6,8) * var(:,:,6+ol) + &
                           q(7,8) * var(:,:,7+ol) + q(9,8) * var(:,:,9+ol) + &
                           q(10,8) * var(:,:,10+ol) + q(11,8) * var(:,:,11+ol) + &
                           q(12,8) * var(:,:,12+ol) ) * idel
        kl = 9 + ol
      end if
      if ( bb(2) == 0 ) then
        kr = nk - gsize
      else 
        or = nk - offset(2)
        dvar(:,:,or) =   - ( q(1,1) * var(:,:,or) + q(2,1) * var(:,:,or-1) + &
                             q(3,1) * var(:,:,or-2) + q(4,1) * var(:,:,or-3) + &
                             q(5,1) * var(:,:,or-4) + q(6,1) * var(:,:,or-5) + &
                             q(7,1) * var(:,:,or-6) + &
                             q(8,1) * var(:,:,or-7) ) * idel
        dvar(:,:,or-1) = - ( q(1,2) * var(:,:,or) + q(3,2) * var(:,:,or-2) + &
                             q(4,2) * var(:,:,or-3) + q(5,2) * var(:,:,or-4) + &
                             q(6,2) * var(:,:,or-5) + q(7,2) * var(:,:,or-6) + &
                             q(8,2) * var(:,:,or-7) ) * idel
        dvar(:,:,or-2) = - ( q(1,3) * var(:,:,or) + q(2,3) * var(:,:,or-1) + &
                             q(4,3) * var(:,:,or-3) + q(5,3) * var(:,:,or-4) + &
                             q(6,3) * var(:,:,or-5) + q(7,3) * var(:,:,or-6) + &
                             q(8,3) * var(:,:,or-7) ) * idel
        dvar(:,:,or-3) = - ( q(1,4) * var(:,:,or) + q(2,4) * var(:,:,or-1) + &
                             q(3,4) * var(:,:,or-2) + q(5,4) * var(:,:,or-4) + &
                             q(6,4) * var(:,:,or-5) + q(7,4) * var(:,:,or-6) + &
                             q(8,4) * var(:,:,or-7) ) * idel
        dvar(:,:,or-4) = - ( q(1,5) * var(:,:,or) + q(2,5) * var(:,:,or-1) + &
                             q(3,5) * var(:,:,or-2) + q(4,5) * var(:,:,or-3) + &
                             q(6,5) * var(:,:,or-5) + q(7,5) * var(:,:,or-6) + &
                             q(8,5) * var(:,:,or-7) + &
                             q(9,5) * var(:,:,or-8) ) * idel
        dvar(:,:,or-5) = - ( q(1,6) * var(:,:,or) + q(2,6) * var(:,:,or-1) + &
                             q(3,6) * var(:,:,or-2) + q(4,6) * var(:,:,or-3) + &
                             q(5,6) * var(:,:,or-4) + q(7,6) * var(:,:,or-6) + &
                             q(8,6) * var(:,:,or-7) + q(9,6) * var(:,:,or-8) + &
                             q(10,6) * var(:,:,or-9) ) * idel
        dvar(:,:,or-6) = - ( q(1,7) * var(:,:,or) + q(2,7) * var(:,:,or-1) + &
                             q(3,7) * var(:,:,or-2) + q(4,7) * var(:,:,or-3) + &
                             q(5,7) * var(:,:,or-4) + q(6,7) * var(:,:,or-5) + &
                             q(8,7) * var(:,:,or-7) + q(9,7) * var(:,:,or-8) + &
                             q(10,7) * var(:,:,or-9) + &
                             q(11,7) * var(:,:,or-10) ) * idel
        dvar(:,:,or-7) = - ( q(1,8) * var(:,:,or) + q(2,8) * var(:,:,or-1) + &
                             q(3,8) * var(:,:,or-2) + q(4,8) * var(:,:,or-3) + &
                             q(5,8) * var(:,:,or-4) + q(6,8) * var(:,:,or-5) + &
                             q(7,8) * var(:,:,or-6) + q(9,8) * var(:,:,or-8) + &
                             q(10,8) * var(:,:,or-9) + &
                             q(11,8) * var(:,:,or-10) + &
                             q(12,8) * var(:,:,or-11) ) * idel
        kr = or - 8
      end if
      if (kl > kr+1) call CCTK_WARN (0, "domain too small")
      dvar(:,:,kl:kr) = ( a(1) * ( var(:,:,kl+1:kr+1) - &
                                   var(:,:,kl-1:kr-1) ) + &
                          a(2) * ( var(:,:,kl+2:kr+2) - &
                                   var(:,:,kl-2:kr-2) ) + &
                          a(3) * ( var(:,:,kl+3:kr+3) - &
                                   var(:,:,kl-3:kr-3) ) + &
                          a(4) * ( var(:,:,kl+4:kr+4) - &
                                   var(:,:,kl-4:kr-4) ) ) * idel
    end if
  end select direction
end subroutine deriv_gf_8_4
