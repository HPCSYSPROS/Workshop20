#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine deriv_gf_8_4_opt ( var, ni, nj, nk, dir, bb, gsize, offset, delta, dvar )

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
    call coeffs_1_8_4_opt ( a, q )
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
end subroutine deriv_gf_8_4_opt

subroutine up_deriv_gf_8_4_opt ( var, ni, nj, nk, dir, bb, gsize, delta, up, dvar )

  use All_Coeffs_mod

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(IN) :: ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(IN) :: var
  CCTK_INT, intent(IN) :: dir
  CCTK_INT, intent(IN) :: bb(2)
  CCTK_INT, intent(IN) :: gsize
  CCTK_REAL, intent(IN) :: delta
  CCTK_REAL, dimension(ni,nj,nk), intent(IN) :: up
  CCTK_REAL, dimension(ni,nj,nk), intent(OUT) :: dvar

  CCTK_REAL, dimension(-4:4), save :: a1, a2
  CCTK_REAL, dimension(12,8), save :: q1, q2
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_up_8_4_opt ( a1, q1, a2, q2 )
    first = .false. 
  end if

  idel = 1.0_wp / delta

  if (gsize < 4) call CCTK_WARN (0, "not enough ghostzones")

  direction: select case (dir)
  case (0) direction
    if ( bb(1) == 0 ) then
      il = 1 + gsize
    else
      where ( up(1,:,:) < zero )
        dvar(1,:,:) = ( q1(1,1) * var(1,:,:) + q1(2,1) * var(2,:,:) + &
                        q1(3,1) * var(3,:,:) + q1(4,1) * var(4,:,:) + &
                        q1(5,1) * var(5,:,:) + q1(6,1) * var(6,:,:) + &
                        q1(7,1) * var(7,:,:) + q1(8,1) * var(8,:,:) ) * idel
      elsewhere
        dvar(1,:,:) = ( q2(1,1) * var(1,:,:) + q2(2,1) * var(2,:,:) + &
                        q2(3,1) * var(3,:,:) + q2(4,1) * var(4,:,:) + &
                        q2(5,1) * var(5,:,:) + q2(6,1) * var(6,:,:) + &
                        q2(7,1) * var(7,:,:) + q2(8,1) * var(8,:,:) ) * idel
      end where
      where ( up(2,:,:) < zero )
        dvar(2,:,:) = ( q1(1,2) * var(1,:,:) + q1(2,2) * var(2,:,:) + &
                        q1(3,2) * var(3,:,:) + q1(4,2) * var(4,:,:) + &
                        q1(5,2) * var(5,:,:) + q1(6,2) * var(6,:,:) + &
                        q1(7,2) * var(7,:,:) + q1(8,2) * var(8,:,:) ) * idel
      elsewhere
        dvar(2,:,:) = ( q2(1,2) * var(1,:,:) + q2(2,2) * var(2,:,:) + &
                        q2(3,2) * var(3,:,:) + q2(4,2) * var(4,:,:) + &
                        q2(5,2) * var(5,:,:) + q2(6,2) * var(6,:,:) + &
                        q2(7,2) * var(7,:,:) + q2(8,2) * var(8,:,:) ) * idel
      end where
      where ( up(3,:,:) < zero )
        dvar(3,:,:) = ( q1(1,3) * var(1,:,:) + q1(2,3) * var(2,:,:) + &
                        q1(3,3) * var(3,:,:) + q1(4,3) * var(4,:,:) + &
                        q1(5,3) * var(5,:,:) + q1(6,3) * var(6,:,:) + &
                        q1(7,3) * var(7,:,:) + q1(8,3) * var(8,:,:) ) * idel
      elsewhere
        dvar(3,:,:) = ( q2(1,3) * var(1,:,:) + q2(2,3) * var(2,:,:) + &
                        q2(3,3) * var(3,:,:) + q2(4,3) * var(4,:,:) + &
                        q2(5,3) * var(5,:,:) + q2(6,3) * var(6,:,:) + &
                        q2(7,3) * var(7,:,:) + q2(8,3) * var(8,:,:) ) * idel
      end where
      where ( up(4,:,:) < zero )
        dvar(4,:,:) = ( q1(1,4) * var(1,:,:) + q1(2,4) * var(2,:,:) + &
                        q1(3,4) * var(3,:,:) + q1(4,4) * var(4,:,:) + &
                        q1(5,4) * var(5,:,:) + q1(6,4) * var(6,:,:) + &
                        q1(7,4) * var(7,:,:) + q1(8,4) * var(8,:,:) ) * idel
      elsewhere
        dvar(4,:,:) = ( q2(1,4) * var(1,:,:) + q2(2,4) * var(2,:,:) + &
                        q2(3,4) * var(3,:,:) + q2(4,4) * var(4,:,:) + &
                        q2(5,4) * var(5,:,:) + q2(6,4) * var(6,:,:) + &
                        q2(7,4) * var(7,:,:) + q2(8,4) * var(8,:,:) ) * idel
      end where
      where ( up(5,:,:) < zero )
        dvar(5,:,:) = ( q1(1,5) * var(1,:,:) + q1(2,5) * var(2,:,:) + &
                        q1(3,5) * var(3,:,:) + q1(4,5) * var(4,:,:) + &
                        q1(5,5) * var(5,:,:) + q1(6,5) * var(6,:,:) + &
                        q1(7,5) * var(7,:,:) + q1(8,5) * var(8,:,:) + &
                        q1(9,5) * var(9,:,:) ) * idel
      elsewhere
        dvar(5,:,:) = ( q2(1,5) * var(1,:,:) + q2(2,5) * var(2,:,:) + &
                        q2(3,5) * var(3,:,:) + q2(4,5) * var(4,:,:) + &
                        q2(5,5) * var(5,:,:) + q2(6,5) * var(6,:,:) + &
                        q2(7,5) * var(7,:,:) + q2(8,5) * var(8,:,:) + &
                        q2(9,5) * var(9,:,:) ) * idel
      end where
      where ( up(6,:,:) < zero )
        dvar(6,:,:) = ( q1(1,6) * var(1,:,:) + q1(2,6) * var(2,:,:) + &
                        q1(3,6) * var(3,:,:) + q1(4,6) * var(4,:,:) + &
                        q1(5,6) * var(5,:,:) + q1(6,6) * var(6,:,:) + &
                        q1(7,6) * var(7,:,:) + q1(8,6) * var(8,:,:) + &
                        q1(9,6) * var(9,:,:) + q1(10,6) * var(10,:,:) ) * idel
      elsewhere
        dvar(6,:,:) = ( q2(1,6) * var(1,:,:) + q2(2,6) * var(2,:,:) + &
                        q2(3,6) * var(3,:,:) + q2(4,6) * var(4,:,:) + &
                        q2(5,6) * var(5,:,:) + q2(6,6) * var(6,:,:) + &
                        q2(7,6) * var(7,:,:) + q2(8,6) * var(8,:,:) + &
                        q2(9,6) * var(9,:,:) + q2(10,6) * var(10,:,:) ) * idel
      end where
      where ( up(7,:,:) < zero )
        dvar(7,:,:) = ( q1(1,7) * var(1,:,:) + q1(2,7) * var(2,:,:) + &
                        q1(3,7) * var(3,:,:) + q1(4,7) * var(4,:,:) + &
                        q1(5,7) * var(5,:,:) + q1(6,7) * var(6,:,:) + &
                        q1(7,7) * var(7,:,:) + q1(8,7) * var(8,:,:) + &
                        q1(9,7) * var(9,:,:) + q1(10,7) * var(10,:,:) + &
                        q1(11,7) * var(11,:,:) ) * idel
      elsewhere
        dvar(7,:,:) = ( q2(1,7) * var(1,:,:) + q2(2,7) * var(2,:,:) + &
                        q2(3,7) * var(3,:,:) + q2(4,7) * var(4,:,:) + &
                        q2(5,7) * var(5,:,:) + q2(6,7) * var(6,:,:) + &
                        q2(7,7) * var(7,:,:) + q2(8,7) * var(8,:,:) + &
                        q2(9,7) * var(9,:,:) + q2(10,7) * var(10,:,:) + &
                        q2(11,7) * var(11,:,:) ) * idel
      end where
      where ( up(8,:,:) < zero )
        dvar(8,:,:) = ( q1(1,8) * var(1,:,:) + q1(2,8) * var(2,:,:) + &
                        q1(3,8) * var(3,:,:) + q1(4,8) * var(4,:,:) + &
                        q1(5,8) * var(5,:,:) + q1(6,8) * var(6,:,:) + &
                        q1(7,8) * var(7,:,:) + q1(8,8) * var(8,:,:) + &
                        q1(9,8) * var(9,:,:) + q1(10,8) * var(10,:,:) + &
                        q1(11,8) * var(11,:,:) + &
                        q1(12,8) * var(12,:,:) ) * idel
      elsewhere
        dvar(8,:,:) = ( q2(1,8) * var(1,:,:) + q2(2,8) * var(2,:,:) + &
                        q2(3,8) * var(3,:,:) + q2(4,8) * var(4,:,:) + &
                        q2(5,8) * var(5,:,:) + q2(6,8) * var(6,:,:) + &
                        q2(7,8) * var(7,:,:) + q2(8,8) * var(8,:,:) + &
                        q2(9,8) * var(9,:,:) + q2(10,8) * var(10,:,:) + &
                        q2(11,8) * var(11,:,:) + &
                        q2(12,8) * var(12,:,:) ) * idel
      end where
      il = 9
    end if
    if ( bb(2) == 0 ) then
      ir = ni - gsize
    else
      where ( up(ni,:,:) < zero )
        dvar(ni,:,:) =   - ( q2(1,1) * var(ni,:,:) + &
                             q2(2,1) * var(ni-1,:,:) + &
                             q2(3,1) * var(ni-2,:,:) + &
                             q2(4,1) * var(ni-3,:,:) + &
                             q2(5,1) * var(ni-4,:,:) + &
                             q2(6,1) * var(ni-5,:,:) + &
                             q2(7,1) * var(ni-6,:,:) + &
                             q2(8,1) * var(ni-7,:,:) ) * idel
      elsewhere
        dvar(ni,:,:) =   - ( q1(1,1) * var(ni,:,:) + &
                             q1(2,1) * var(ni-1,:,:) + &
                             q1(3,1) * var(ni-2,:,:) + &
                             q1(4,1) * var(ni-3,:,:) + &
                             q1(5,1) * var(ni-4,:,:) + &
                             q1(6,1) * var(ni-5,:,:) + &
                             q1(7,1) * var(ni-6,:,:) + &
                             q1(8,1) * var(ni-7,:,:) ) * idel
      end where
      where ( up(ni-1,:,:) < zero )
        dvar(ni-1,:,:) = - ( q2(1,2) * var(ni,:,:) + &
                             q2(2,2) * var(ni-1,:,:) + &
                             q2(3,2) * var(ni-2,:,:) + &
                             q2(4,2) * var(ni-3,:,:) + &
                             q2(5,2) * var(ni-4,:,:) + &
                             q2(6,2) * var(ni-5,:,:) + &
                             q2(7,2) * var(ni-6,:,:) + &
                             q2(8,2) * var(ni-7,:,:) ) * idel
      elsewhere
        dvar(ni-1,:,:) = - ( q1(1,2) * var(ni,:,:) + &
                             q1(2,2) * var(ni-1,:,:) + &
                             q1(3,2) * var(ni-2,:,:) + &
                             q1(4,2) * var(ni-3,:,:) + &
                             q1(5,2) * var(ni-4,:,:) + &
                             q1(6,2) * var(ni-5,:,:) + &
                             q1(7,2) * var(ni-6,:,:) + &
                             q1(8,2) * var(ni-7,:,:) ) * idel
      end where
      where ( up(ni-2,:,:) < zero )
        dvar(ni-2,:,:) = - ( q2(1,3) * var(ni,:,:) + &
                             q2(2,3) * var(ni-1,:,:) + &
                             q2(3,3) * var(ni-2,:,:) + &
                             q2(4,3) * var(ni-3,:,:) + &
                             q2(5,3) * var(ni-4,:,:) + &
                             q2(6,3) * var(ni-5,:,:) + &
                             q2(7,3) * var(ni-6,:,:) + &
                             q2(8,3) * var(ni-7,:,:) ) * idel
      elsewhere
        dvar(ni-2,:,:) = - ( q1(1,3) * var(ni,:,:) + &
                             q1(2,3) * var(ni-1,:,:) + &
                             q1(3,3) * var(ni-2,:,:) + &
                             q1(4,3) * var(ni-3,:,:) + &
                             q1(5,3) * var(ni-4,:,:) + &
                             q1(6,3) * var(ni-5,:,:) + &
                             q1(7,3) * var(ni-6,:,:) + &
                             q1(8,3) * var(ni-7,:,:) ) * idel
      end where
      where ( up(ni-3,:,:) < zero )
        dvar(ni-3,:,:) = - ( q2(1,4) * var(ni,:,:) + &
                             q2(2,4) * var(ni-1,:,:) + &
                             q2(3,4) * var(ni-2,:,:) + &
                             q2(4,4) * var(ni-3,:,:) + &
                             q2(5,4) * var(ni-4,:,:) + &
                             q2(6,4) * var(ni-5,:,:) + &
                             q2(7,4) * var(ni-6,:,:) + &
                             q2(8,4) * var(ni-7,:,:) ) * idel
      elsewhere
        dvar(ni-3,:,:) = - ( q1(1,4) * var(ni,:,:) + &
                             q1(2,4) * var(ni-1,:,:) + &
                             q1(3,4) * var(ni-2,:,:) + &
                             q1(4,4) * var(ni-3,:,:) + &
                             q1(5,4) * var(ni-4,:,:) + &
                             q1(6,4) * var(ni-5,:,:) + &
                             q1(7,4) * var(ni-6,:,:) + &
                             q1(8,4) * var(ni-7,:,:) ) * idel
      end where
      where ( up(ni-4,:,:) < zero )
        dvar(ni-4,:,:) = - ( q2(1,5) * var(ni,:,:) + &
                             q2(2,5) * var(ni-1,:,:) + &
                             q2(3,5) * var(ni-2,:,:) + &
                             q2(4,5) * var(ni-3,:,:) + &
                             q2(5,5) * var(ni-4,:,:) + &
                             q2(6,5) * var(ni-5,:,:) + &
                             q2(7,5) * var(ni-6,:,:) + &
                             q2(8,5) * var(ni-7,:,:) + &
                             q2(9,5) * var(ni-8,:,:) ) * idel
      elsewhere
        dvar(ni-4,:,:) = - ( q1(1,5) * var(ni,:,:) + &
                             q1(2,5) * var(ni-1,:,:) + &
                             q1(3,5) * var(ni-2,:,:) + &
                             q1(4,5) * var(ni-3,:,:) + &
                             q1(5,5) * var(ni-4,:,:) + &
                             q1(6,5) * var(ni-5,:,:) + &
                             q1(7,5) * var(ni-6,:,:) + &
                             q1(8,5) * var(ni-7,:,:) + &
                             q1(9,5) * var(ni-8,:,:) ) * idel
      end where
      where ( up(ni-5,:,:) < zero )
        dvar(ni-5,:,:) = - ( q2(1,6) * var(ni,:,:) + &
                             q2(2,6) * var(ni-1,:,:) + &
                             q2(3,6) * var(ni-2,:,:) + &
                             q2(4,6) * var(ni-3,:,:) + &
                             q2(5,6) * var(ni-4,:,:) + &
                             q2(6,6) * var(ni-5,:,:) + &
                             q2(7,6) * var(ni-6,:,:) + &
                             q2(8,6) * var(ni-7,:,:) + &
                             q2(9,6) * var(ni-8,:,:) + &
                             q2(10,6) * var(ni-9,:,:) ) * idel
      elsewhere
        dvar(ni-5,:,:) = - ( q1(1,6) * var(ni,:,:) + &
                             q1(2,6) * var(ni-1,:,:) + &
                             q1(3,6) * var(ni-2,:,:) + &
                             q1(4,6) * var(ni-3,:,:) + &
                             q1(5,6) * var(ni-4,:,:) + &
                             q1(6,6) * var(ni-5,:,:) + &
                             q1(7,6) * var(ni-6,:,:) + &
                             q1(8,6) * var(ni-7,:,:) + &
                             q1(9,6) * var(ni-8,:,:) + &
                             q1(10,6) * var(ni-9,:,:) ) * idel
      end where
      where ( up(ni-6,:,:) < zero )
        dvar(ni-6,:,:) = - ( q2(1,7) * var(ni,:,:) + &
                             q2(2,7) * var(ni-1,:,:) + &
                             q2(3,7) * var(ni-2,:,:) + &
                             q2(4,7) * var(ni-3,:,:) + &
                             q2(5,7) * var(ni-4,:,:) + &
                             q2(6,7) * var(ni-5,:,:) + &
                             q2(7,7) * var(ni-6,:,:) + &
                             q2(8,7) * var(ni-7,:,:) + &
                             q2(9,7) * var(ni-8,:,:) + &
                             q2(10,7) * var(ni-9,:,:) + &
                             q2(11,7) * var(ni-10,:,:) ) * idel
      elsewhere
        dvar(ni-6,:,:) = - ( q1(1,7) * var(ni,:,:) + &
                             q1(2,7) * var(ni-1,:,:) + &
                             q1(3,7) * var(ni-2,:,:) + &
                             q1(4,7) * var(ni-3,:,:) + &
                             q1(5,7) * var(ni-4,:,:) + &
                             q1(6,7) * var(ni-5,:,:) + &
                             q1(7,7) * var(ni-6,:,:) + &
                             q1(8,7) * var(ni-7,:,:) + &
                             q1(9,7) * var(ni-8,:,:) + &
                             q1(10,7) * var(ni-9,:,:) + &
                             q1(11,7) * var(ni-10,:,:) ) * idel
      end where
      where ( up(ni-7,:,:) < zero )
        dvar(ni-7,:,:) = - ( q2(1,8) * var(ni,:,:) + &
                             q2(2,8) * var(ni-1,:,:) + &
                             q2(3,8) * var(ni-2,:,:) + &
                             q2(4,8) * var(ni-3,:,:) + &
                             q2(5,8) * var(ni-4,:,:) + &
                             q2(6,8) * var(ni-5,:,:) + &
                             q2(7,8) * var(ni-6,:,:) + &
                             q2(8,8) * var(ni-7,:,:) + &
                             q2(9,8) * var(ni-8,:,:) + &
                             q2(10,8) * var(ni-9,:,:) + &
                             q2(11,8) * var(ni-10,:,:) + &
                             q2(12,8) * var(ni-11,:,:) ) * idel
      elsewhere
        dvar(ni-7,:,:) = - ( q1(1,8) * var(ni,:,:) + &
                             q1(2,8) * var(ni-1,:,:) + &
                             q1(3,8) * var(ni-2,:,:) + &
                             q1(4,8) * var(ni-3,:,:) + &
                             q1(5,8) * var(ni-4,:,:) + &
                             q1(6,8) * var(ni-5,:,:) + &
                             q1(7,8) * var(ni-6,:,:) + &
                             q1(8,8) * var(ni-7,:,:) + &
                             q1(9,8) * var(ni-8,:,:) + &
                             q1(10,8) * var(ni-9,:,:) + &
                             q1(11,8) * var(ni-10,:,:) + &
                             q1(12,8) * var(ni-11,:,:) ) * idel
      end where
      ir = ni - 8
    end if
    if (il > ir+1) call CCTK_WARN (0, "domain too small")
    where ( up(il:ir,:,:) < zero )
      dvar(il:ir,:,:) = ( a1(-4) * var(il-4:ir-4,:,:) + &
                          a1(-3) * var(il-3:ir-3,:,:) + &
                          a1(-2) * var(il-2:ir-2,:,:) + &
                          a1(-1) * var(il-1:ir-1,:,:) + &
                          a1(0) * var(il:ir,:,:) + &
                          a1(1) * var(il+1:ir+1,:,:) + &
                          a1(2) * var(il+2:ir+2,:,:) + &
                          a1(3) * var(il+3:ir+3,:,:) + &
                          a1(4) * var(il+4:ir+4,:,:) ) * idel
    elsewhere
      dvar(il:ir,:,:) = ( a2(-4) * var(il-4:ir-4,:,:) + &
                          a2(-3) * var(il-3:ir-3,:,:) + &
                          a2(-2) * var(il-2:ir-2,:,:) + &
                          a2(-1) * var(il-1:ir-1,:,:) + &
                          a2(0) * var(il:ir,:,:) + &
                          a2(1) * var(il+1:ir+1,:,:) + &
                          a2(2) * var(il+2:ir+2,:,:) + &
                          a2(3) * var(il+3:ir+3,:,:) + &
                          a2(4) * var(il+4:ir+4,:,:) ) * idel
    end where
  case (1) direction
    if ( zero_derivs_y /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        jl = 1 + gsize
      else
        where ( up(:,1,:) < zero )
          dvar(:,1,:) = ( q1(1,1) * var(:,1,:) + q1(2,1) * var(:,2,:) + &
                          q1(3,1) * var(:,3,:) + q1(4,1) * var(:,4,:) + &
                          q1(5,1) * var(:,5,:) + q1(6,1) * var(:,6,:) + &
                          q1(7,1) * var(:,7,:) + q1(8,1) * var(:,8,:) ) * idel
        elsewhere
          dvar(:,1,:) = ( q2(1,1) * var(:,1,:) + q2(2,1) * var(:,2,:) + &
                          q2(3,1) * var(:,3,:) + q2(4,1) * var(:,4,:) + &
                          q2(5,1) * var(:,5,:) + q2(6,1) * var(:,6,:) + &
                          q2(7,1) * var(:,7,:) + q2(8,1) * var(:,8,:) ) * idel
        end where
        where ( up(:,2,:) < zero )
          dvar(:,2,:) = ( q1(1,2) * var(:,1,:) + q1(2,2) * var(:,2,:) + &
                          q1(3,2) * var(:,3,:) + q1(4,2) * var(:,4,:) + &
                          q1(5,2) * var(:,5,:) + q1(6,2) * var(:,6,:) + &
                          q1(7,2) * var(:,7,:) + q1(8,2) * var(:,8,:) ) * idel
        elsewhere
          dvar(:,2,:) = ( q2(1,2) * var(:,1,:) + q2(2,2) * var(:,2,:) + &
                          q2(3,2) * var(:,3,:) + q2(4,2) * var(:,4,:) + &
                          q2(5,2) * var(:,5,:) + q2(6,2) * var(:,6,:) + &
                          q2(7,2) * var(:,7,:) + q2(8,2) * var(:,8,:) ) * idel
        end where
        where ( up(:,3,:) < zero )
          dvar(:,3,:) = ( q1(1,3) * var(:,1,:) + q1(2,3) * var(:,2,:) + &
                          q1(3,3) * var(:,3,:) + q1(4,3) * var(:,4,:) + &
                          q1(5,3) * var(:,5,:) + q1(6,3) * var(:,6,:) + &
                          q1(7,3) * var(:,7,:) + q1(8,3) * var(:,8,:) ) * idel
        elsewhere
          dvar(:,3,:) = ( q2(1,3) * var(:,1,:) + q2(2,3) * var(:,2,:) + &
                          q2(3,3) * var(:,3,:) + q2(4,3) * var(:,4,:) + &
                          q2(5,3) * var(:,5,:) + q2(6,3) * var(:,6,:) + &
                          q2(7,3) * var(:,7,:) + q2(8,3) * var(:,8,:) ) * idel
        end where
        where ( up(:,4,:) < zero )
          dvar(:,4,:) = ( q1(1,4) * var(:,1,:) + q1(2,4) * var(:,2,:) + &
                          q1(3,4) * var(:,3,:) + q1(4,4) * var(:,4,:) + &
                          q1(5,4) * var(:,5,:) + q1(6,4) * var(:,6,:) + &
                          q1(7,4) * var(:,7,:) + q1(8,4) * var(:,8,:) ) * idel
        elsewhere
          dvar(:,4,:) = ( q2(1,4) * var(:,1,:) + q2(2,4) * var(:,2,:) + &
                          q2(3,4) * var(:,3,:) + q2(4,4) * var(:,4,:) + &
                          q2(5,4) * var(:,5,:) + q2(6,4) * var(:,6,:) + &
                          q2(7,4) * var(:,7,:) + q2(8,4) * var(:,8,:) ) * idel
        end where
        where ( up(:,5,:) < zero )
          dvar(:,5,:) = ( q1(1,5) * var(:,1,:) + q1(2,5) * var(:,2,:) + &
                          q1(3,5) * var(:,3,:) + q1(4,5) * var(:,4,:) + &
                          q1(5,5) * var(:,5,:) + q1(6,5) * var(:,6,:) + &
                          q1(7,5) * var(:,7,:) + q1(8,5) * var(:,8,:) + &
                          q1(9,5) * var(:,9,:) ) * idel
        elsewhere
          dvar(:,5,:) = ( q2(1,5) * var(:,1,:) + q2(2,5) * var(:,2,:) + &
                          q2(3,5) * var(:,3,:) + q2(4,5) * var(:,4,:) + &
                          q2(5,5) * var(:,5,:) + q2(6,5) * var(:,6,:) + &
                          q2(7,5) * var(:,7,:) + q2(8,5) * var(:,8,:) + &
                          q2(9,5) * var(:,9,:) ) * idel
        end where
        where ( up(:,6,:) < zero )
          dvar(:,6,:) = ( q1(1,6) * var(:,1,:) + q1(2,6) * var(:,2,:) + &
                          q1(3,6) * var(:,3,:) + q1(4,6) * var(:,4,:) + &
                          q1(5,6) * var(:,5,:) + q1(6,6) * var(:,6,:) + &
                          q1(7,6) * var(:,7,:) + q1(8,6) * var(:,8,:) + &
                          q1(9,6) * var(:,9,:) + q1(10,6) * var(:,10,:) ) * idel
        elsewhere
          dvar(:,6,:) = ( q2(1,6) * var(:,1,:) + q2(2,6) * var(:,2,:) + &
                          q2(3,6) * var(:,3,:) + q2(4,6) * var(:,4,:) + &
                          q2(5,6) * var(:,5,:) + q2(6,6) * var(:,6,:) + &
                          q2(7,6) * var(:,7,:) + q2(8,6) * var(:,8,:) + &
                          q2(9,6) * var(:,9,:) + q2(10,6) * var(:,10,:) ) * idel
        end where
        where ( up(:,7,:) < zero )
          dvar(:,7,:) = ( q1(1,7) * var(:,1,:) + q1(2,7) * var(:,2,:) + &
                          q1(3,7) * var(:,3,:) + q1(4,7) * var(:,4,:) + &
                          q1(5,7) * var(:,5,:) + q1(6,7) * var(:,6,:) + &
                          q1(7,7) * var(:,7,:) + q1(8,7) * var(:,8,:) + &
                          q1(9,7) * var(:,9,:) + q1(10,7) * var(:,10,:) + &
                          q1(11,7) * var(:,11,:) ) * idel
        elsewhere
          dvar(:,7,:) = ( q2(1,7) * var(:,1,:) + q2(2,7) * var(:,2,:) + &
                          q2(3,7) * var(:,3,:) + q2(4,7) * var(:,4,:) + &
                          q2(5,7) * var(:,5,:) + q2(6,7) * var(:,6,:) + &
                          q2(7,7) * var(:,7,:) + q2(8,7) * var(:,8,:) + &
                          q2(9,7) * var(:,9,:) + q2(10,7) * var(:,10,:) + &
                          q2(11,7) * var(:,11,:) ) * idel
        end where
        where ( up(:,8,:) < zero )
          dvar(:,8,:) = ( q1(1,8) * var(:,1,:) + q1(2,8) * var(:,2,:) + &
                          q1(3,8) * var(:,3,:) + q1(4,8) * var(:,4,:) + &
                          q1(5,8) * var(:,5,:) + q1(6,8) * var(:,6,:) + &
                          q1(7,8) * var(:,7,:) + q1(8,8) * var(:,8,:) + &
                          q1(9,8) * var(:,9,:) + q1(10,8) * var(:,10,:) + &
                          q1(11,8) * var(:,11,:) + &
                          q1(12,8) * var(:,12,:) ) * idel
        elsewhere
          dvar(:,8,:) = ( q2(1,8) * var(:,1,:) + q2(2,8) * var(:,2,:) + &
                          q2(3,8) * var(:,3,:) + q2(4,8) * var(:,4,:) + &
                          q2(5,8) * var(:,5,:) + q2(6,8) * var(:,6,:) + &
                          q2(7,8) * var(:,7,:) + q2(8,8) * var(:,8,:) + &
                          q2(9,8) * var(:,9,:) + q2(10,8) * var(:,10,:) + &
                          q2(11,8) * var(:,11,:) + &
                          q2(12,8) * var(:,12,:) ) * idel
        end where
        jl = 9
      end if
      if ( bb(2) == 0 ) then
        jr = nj - gsize
      else
        where ( up(:,nj,:) < zero )
          dvar(:,nj,:) =   - ( q2(1,1) * var(:,nj,:) + &
                               q2(2,1) * var(:,nj-1,:) + &
                               q2(3,1) * var(:,nj-2,:) + &
                               q2(4,1) * var(:,nj-3,:) + &
                               q2(5,1) * var(:,nj-4,:) + &
                               q2(6,1) * var(:,nj-5,:) + &
                               q2(7,1) * var(:,nj-6,:) + &
                               q2(8,1) * var(:,nj-7,:) ) * idel
        elsewhere
          dvar(:,nj,:) =   - ( q1(1,1) * var(:,nj,:) + &
                               q1(2,1) * var(:,nj-1,:) + &
                               q1(3,1) * var(:,nj-2,:) + &
                               q1(4,1) * var(:,nj-3,:) + &
                               q1(5,1) * var(:,nj-4,:) + &
                               q1(6,1) * var(:,nj-5,:) + &
                               q1(7,1) * var(:,nj-6,:) + &
                               q1(8,1) * var(:,nj-7,:) ) * idel
        end where
        where ( up(:,nj-1,:) < zero )
          dvar(:,nj-1,:) = - ( q2(1,2) * var(:,nj,:) + &
                               q2(2,2) * var(:,nj-1,:) + &
                               q2(3,2) * var(:,nj-2,:) + &
                               q2(4,2) * var(:,nj-3,:) + &
                               q2(5,2) * var(:,nj-4,:) + &
                               q2(6,2) * var(:,nj-5,:) + &
                               q2(7,2) * var(:,nj-6,:) + &
                               q2(8,2) * var(:,nj-7,:) ) * idel
        elsewhere
          dvar(:,nj-1,:) = - ( q1(1,2) * var(:,nj,:) + &
                               q1(2,2) * var(:,nj-1,:) + &
                               q1(3,2) * var(:,nj-2,:) + &
                               q1(4,2) * var(:,nj-3,:) + &
                               q1(5,2) * var(:,nj-4,:) + &
                               q1(6,2) * var(:,nj-5,:) + &
                               q1(7,2) * var(:,nj-6,:) + &
                               q1(8,2) * var(:,nj-7,:) ) * idel
        end where
        where ( up(:,nj-2,:) < zero )
          dvar(:,nj-2,:) = - ( q2(1,3) * var(:,nj,:) + &
                               q2(2,3) * var(:,nj-1,:) + &
                               q2(3,3) * var(:,nj-2,:) + &
                               q2(4,3) * var(:,nj-3,:) + &
                               q2(5,3) * var(:,nj-4,:) + &
                               q2(6,3) * var(:,nj-5,:) + &
                               q2(7,3) * var(:,nj-6,:) + &
                               q2(8,3) * var(:,nj-7,:) ) * idel
        elsewhere
          dvar(:,nj-2,:) = - ( q1(1,3) * var(:,nj,:) + &
                               q1(2,3) * var(:,nj-1,:) + &
                               q1(3,3) * var(:,nj-2,:) + &
                               q1(4,3) * var(:,nj-3,:) + &
                               q1(5,3) * var(:,nj-4,:) + &
                               q1(6,3) * var(:,nj-5,:) + &
                               q1(7,3) * var(:,nj-6,:) + &
                               q1(8,3) * var(:,nj-7,:) ) * idel
        end where
        where ( up(:,nj-3,:) < zero )
          dvar(:,nj-3,:) = - ( q2(1,4) * var(:,nj,:) + &
                               q2(2,4) * var(:,nj-1,:) + &
                               q2(3,4) * var(:,nj-2,:) + &
                               q2(4,4) * var(:,nj-3,:) + &
                               q2(5,4) * var(:,nj-4,:) + &
                               q2(6,4) * var(:,nj-5,:) + &
                               q2(7,4) * var(:,nj-6,:) + &
                               q2(8,4) * var(:,nj-7,:) ) * idel
        elsewhere
          dvar(:,nj-3,:) = - ( q1(1,4) * var(:,nj,:) + &
                               q1(2,4) * var(:,nj-1,:) + &
                               q1(3,4) * var(:,nj-2,:) + &
                               q1(4,4) * var(:,nj-3,:) + &
                               q1(5,4) * var(:,nj-4,:) + &
                               q1(6,4) * var(:,nj-5,:) + &
                               q1(7,4) * var(:,nj-6,:) + &
                               q1(8,4) * var(:,nj-7,:) ) * idel
        end where
        where ( up(:,nj-4,:) < zero )
          dvar(:,nj-4,:) = - ( q2(1,5) * var(:,nj,:) + &
                               q2(2,5) * var(:,nj-1,:) + &
                               q2(3,5) * var(:,nj-2,:) + &
                               q2(4,5) * var(:,nj-3,:) + &
                               q2(5,5) * var(:,nj-4,:) + &
                               q2(6,5) * var(:,nj-5,:) + &
                               q2(7,5) * var(:,nj-6,:) + &
                               q2(8,5) * var(:,nj-7,:) + &
                               q2(9,5) * var(:,nj-8,:) ) * idel
        elsewhere
          dvar(:,nj-4,:) = - ( q1(1,5) * var(:,nj,:) + &
                               q1(2,5) * var(:,nj-1,:) + &
                               q1(3,5) * var(:,nj-2,:) + &
                               q1(4,5) * var(:,nj-3,:) + &
                               q1(5,5) * var(:,nj-4,:) + &
                               q1(6,5) * var(:,nj-5,:) + &
                               q1(7,5) * var(:,nj-6,:) + &
                               q1(8,5) * var(:,nj-7,:) + &
                               q1(9,5) * var(:,nj-8,:) ) * idel
        end where
        where ( up(:,nj-5,:) < zero )
          dvar(:,nj-5,:) = - ( q2(1,6) * var(:,nj,:) + &
                               q2(2,6) * var(:,nj-1,:) + &
                               q2(3,6) * var(:,nj-2,:) + &
                               q2(4,6) * var(:,nj-3,:) + &
                               q2(5,6) * var(:,nj-4,:) + &
                               q2(6,6) * var(:,nj-5,:) + &
                               q2(7,6) * var(:,nj-6,:) + &
                               q2(8,6) * var(:,nj-7,:) + &
                               q2(9,6) * var(:,nj-8,:) + &
                               q2(10,6) * var(:,nj-9,:) ) * idel
        elsewhere
          dvar(:,nj-5,:) = - ( q1(1,6) * var(:,nj,:) + &
                               q1(2,6) * var(:,nj-1,:) + &
                               q1(3,6) * var(:,nj-2,:) + &
                               q1(4,6) * var(:,nj-3,:) + &
                               q1(5,6) * var(:,nj-4,:) + &
                               q1(6,6) * var(:,nj-5,:) + &
                               q1(7,6) * var(:,nj-6,:) + &
                               q1(8,6) * var(:,nj-7,:) + &
                               q1(9,6) * var(:,nj-8,:) + &
                               q1(10,6) * var(:,nj-9,:) ) * idel
        end where
        where ( up(:,nj-6,:) < zero )
          dvar(:,nj-6,:) = - ( q2(1,7) * var(:,nj,:) + &
                               q2(2,7) * var(:,nj-1,:) + &
                               q2(3,7) * var(:,nj-2,:) + &
                               q2(4,7) * var(:,nj-3,:) + &
                               q2(5,7) * var(:,nj-4,:) + &
                               q2(6,7) * var(:,nj-5,:) + &
                               q2(7,7) * var(:,nj-6,:) + &
                               q2(8,7) * var(:,nj-7,:) + &
                               q2(9,7) * var(:,nj-8,:) + &
                               q2(10,7) * var(:,nj-9,:) + &
                               q2(11,7) * var(:,nj-10,:) ) * idel
        elsewhere
          dvar(:,nj-6,:) = - ( q1(1,7) * var(:,nj,:) + &
                               q1(2,7) * var(:,nj-1,:) + &
                               q1(3,7) * var(:,nj-2,:) + &
                               q1(4,7) * var(:,nj-3,:) + &
                               q1(5,7) * var(:,nj-4,:) + &
                               q1(6,7) * var(:,nj-5,:) + &
                               q1(7,7) * var(:,nj-6,:) + &
                               q1(8,7) * var(:,nj-7,:) + &
                               q1(9,7) * var(:,nj-8,:) + &
                               q1(10,7) * var(:,nj-9,:) + &
                               q1(11,7) * var(:,nj-10,:) ) * idel
        end where
        where ( up(:,nj-7,:) < zero )
          dvar(:,nj-7,:) = - ( q2(1,8) * var(:,nj,:) + &
                               q2(2,8) * var(:,nj-1,:) + &
                               q2(3,8) * var(:,nj-2,:) + &
                               q2(4,8) * var(:,nj-3,:) + &
                               q2(5,8) * var(:,nj-4,:) + &
                               q2(6,8) * var(:,nj-5,:) + &
                               q2(7,8) * var(:,nj-6,:) + &
                               q2(8,8) * var(:,nj-7,:) + &
                               q2(9,8) * var(:,nj-8,:) + &
                               q2(10,8) * var(:,nj-9,:) + &
                               q2(11,8) * var(:,nj-10,:) + &
                               q2(12,8) * var(:,nj-11,:) ) * idel
        elsewhere
          dvar(:,nj-7,:) = - ( q1(1,8) * var(:,nj,:) + &
                               q1(2,8) * var(:,nj-1,:) + &
                               q1(3,8) * var(:,nj-2,:) + &
                               q1(4,8) * var(:,nj-3,:) + &
                               q1(5,8) * var(:,nj-4,:) + &
                               q1(6,8) * var(:,nj-5,:) + &
                               q1(7,8) * var(:,nj-6,:) + &
                               q1(8,8) * var(:,nj-7,:) + &
                               q1(9,8) * var(:,nj-8,:) + &
                               q1(10,8) * var(:,nj-9,:) + &
                               q1(11,8) * var(:,nj-10,:) + &
                               q1(12,8) * var(:,nj-11,:) ) * idel
        end where
        jr = nj - 8
      end if
      if (jl > jr+1) call CCTK_WARN (0, "domain too small")
      where ( up(:,jl:jr,:) < zero )
        dvar(:,jl:jr,:) = ( a1(-4) * var(:,jl-4:jr-4,:) + &
                            a1(-3) * var(:,jl-3:jr-3,:) + &
                            a1(-2) * var(:,jl-2:jr-2,:) + &
                            a1(-1) * var(:,jl-1:jr-1,:) + &
                            a1(0)  * var(:,jl:jr,:) + &
                            a1(1)  * var(:,jl+1:jr+1,:) + &
                            a1(2)  * var(:,jl+2:jr+2,:) + &
                            a1(3)  * var(:,jl+3:jr+3,:) + &
                            a1(4)  * var(:,jl+4:jr+4,:) ) * idel
      elsewhere
        dvar(:,jl:jr,:) = ( a2(-4) * var(:,jl-4:jr-4,:) + &
                            a2(-3) * var(:,jl-3:jr-3,:) + &
                            a2(-2) * var(:,jl-2:jr-2,:) + &
                            a2(-1) * var(:,jl-1:jr-1,:) + &
                            a2(0)  * var(:,jl:jr,:) + &
                            a2(1)  * var(:,jl+1:jr+1,:) + &
                            a2(2)  * var(:,jl+2:jr+2,:) + &
                            a2(3)  * var(:,jl+3:jr+3,:) + &
                            a2(4)  * var(:,jl+4:jr+4,:) ) * idel
      end where
    end if
  case (2) direction
    if ( zero_derivs_z /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        kl = 1 + gsize
      else
        where ( up(:,:,1) < zero )
          dvar(:,:,1) = ( q2(1,1) * var(:,:,1) + q2(2,1) * var(:,:,2) + &
                          q2(3,1) * var(:,:,3) + q2(4,1) * var(:,:,4) + &
                          q2(5,1) * var(:,:,5) + q2(6,1) * var(:,:,6) + &
                          q2(7,1) * var(:,:,7) + q2(8,1) * var(:,:,8) ) * idel
        elsewhere
          dvar(:,:,1) = ( q2(1,1) * var(:,:,1) + q2(2,1) * var(:,:,2) + &
                          q2(3,1) * var(:,:,3) + q2(4,1) * var(:,:,4) + &
                          q2(5,1) * var(:,:,5) + q2(6,1) * var(:,:,6) + &
                          q2(7,1) * var(:,:,7) + q2(8,1) * var(:,:,8) ) * idel
        end where
        where ( up(:,:,2) < zero )
          dvar(:,:,2) = ( q1(1,2) * var(:,:,1) + q1(2,2) * var(:,:,2) + &
                          q1(3,2) * var(:,:,3) + q1(4,2) * var(:,:,4) + &
                          q1(5,2) * var(:,:,5) + q1(6,2) * var(:,:,6) + &
                          q1(7,2) * var(:,:,7) + q1(8,2) * var(:,:,8) ) * idel
        elsewhere
          dvar(:,:,2) = ( q2(1,2) * var(:,:,1) + q2(2,2) * var(:,:,2) + &
                          q2(3,2) * var(:,:,3) + q2(4,2) * var(:,:,4) + &
                          q2(5,2) * var(:,:,5) + q2(6,2) * var(:,:,6) + &
                          q2(7,2) * var(:,:,7) + q2(8,2) * var(:,:,8) ) * idel
        end where
        where ( up(:,:,3) < zero )
          dvar(:,:,3) = ( q1(1,3) * var(:,:,1) + q1(2,3) * var(:,:,2) + &
                          q1(3,3) * var(:,:,3) + q1(4,3) * var(:,:,4) + &
                          q1(5,3) * var(:,:,5) + q1(6,3) * var(:,:,6) + &
                          q1(7,3) * var(:,:,7) + q1(8,3) * var(:,:,8) ) * idel
        elsewhere
          dvar(:,:,3) = ( q2(1,3) * var(:,:,1) + q2(2,3) * var(:,:,2) + &
                          q2(3,3) * var(:,:,3) + q2(4,3) * var(:,:,4) + &
                          q2(5,3) * var(:,:,5) + q2(6,3) * var(:,:,6) + &
                          q2(7,3) * var(:,:,7) + q2(8,3) * var(:,:,8) ) * idel
        end where
        where ( up(:,:,4) < zero )
          dvar(:,:,4) = ( q1(1,4) * var(:,:,1) + q1(2,4) * var(:,:,2) + &
                          q1(3,4) * var(:,:,3) + q1(4,4) * var(:,:,4) + &
                          q1(5,4) * var(:,:,5) + q1(6,4) * var(:,:,6) + &
                          q1(7,4) * var(:,:,7) + q1(8,4) * var(:,:,8) ) * idel
        elsewhere
          dvar(:,:,4) = ( q2(1,4) * var(:,:,1) + q2(2,4) * var(:,:,2) + &
                          q2(3,4) * var(:,:,3) + q2(4,4) * var(:,:,4) + &
                          q2(5,4) * var(:,:,5) + q2(6,4) * var(:,:,6) + &
                          q2(7,4) * var(:,:,7) + q2(8,4) * var(:,:,8) ) * idel
        end where
        where ( up(:,:,5) < zero )
          dvar(:,:,5) = ( q1(1,5) * var(:,:,1) + q1(2,5) * var(:,:,2) + &
                          q1(3,5) * var(:,:,3) + q1(4,5) * var(:,:,4) + &
                          q1(5,5) * var(:,:,5) + q1(6,5) * var(:,:,6) + &
                          q1(7,5) * var(:,:,7) + q1(8,5) * var(:,:,8) + &
                          q1(9,5) * var(:,:,9) ) * idel
        elsewhere
          dvar(:,:,5) = ( q2(1,5) * var(:,:,1) + q2(2,5) * var(:,:,2) + &
                          q2(3,5) * var(:,:,3) + q2(4,5) * var(:,:,4) + &
                          q2(5,5) * var(:,:,5) + q2(6,5) * var(:,:,6) + &
                          q2(7,5) * var(:,:,7) + q2(8,5) * var(:,:,8) + &
                          q2(9,5) * var(:,:,9) ) * idel
        end where
        where ( up(:,:,6) < zero )
          dvar(:,:,6) = ( q1(1,6) * var(:,:,1) + q1(2,6) * var(:,:,2) + &
                          q1(3,6) * var(:,:,3) + q1(4,6) * var(:,:,4) + &
                          q1(5,6) * var(:,:,5) + q1(6,6) * var(:,:,6) + &
                          q1(7,6) * var(:,:,7) + q1(8,6) * var(:,:,8) + &
                          q1(9,6) * var(:,:,9) + q1(10,6) * var(:,:,10) ) * idel
        elsewhere
          dvar(:,:,6) = ( q2(1,6) * var(:,:,1) + q2(2,6) * var(:,:,2) + &
                          q2(3,6) * var(:,:,3) + q2(4,6) * var(:,:,4) + &
                          q2(5,6) * var(:,:,5) + q2(6,6) * var(:,:,6) + &
                          q2(7,6) * var(:,:,7) + q2(8,6) * var(:,:,8) + &
                          q2(9,6) * var(:,:,9) + q2(10,6) * var(:,:,10) ) * idel
        end where
        where ( up(:,:,7) < zero )
          dvar(:,:,7) = ( q1(1,7) * var(:,:,1) + q1(2,7) * var(:,:,2) + &
                          q1(3,7) * var(:,:,3) + q1(4,7) * var(:,:,4) + &
                          q1(5,7) * var(:,:,5) + q1(6,7) * var(:,:,6) + &
                          q1(7,7) * var(:,:,7) + q1(8,7) * var(:,:,8) + &
                          q1(9,7) * var(:,:,9) + q1(10,7) * var(:,:,10) + &
                          q1(11,7) * var(:,:,11) ) * idel
        elsewhere
          dvar(:,:,7) = ( q2(1,7) * var(:,:,1) + q2(2,7) * var(:,:,2) + &
                          q2(3,7) * var(:,:,3) + q2(4,7) * var(:,:,4) + &
                          q2(5,7) * var(:,:,5) + q2(6,7) * var(:,:,6) + &
                          q2(7,7) * var(:,:,7) + q2(8,7) * var(:,:,8) + &
                          q2(9,7) * var(:,:,9) + q2(10,7) * var(:,:,10) + &
                          q2(11,7) * var(:,:,11) ) * idel
        end where
        where ( up(:,:,8) < zero )
          dvar(:,:,8) = ( q1(1,8) * var(:,:,1) + q1(2,8) * var(:,:,2) + &
                          q1(3,8) * var(:,:,3) + q1(4,8) * var(:,:,4) + &
                          q1(5,8) * var(:,:,5) + q1(6,8) * var(:,:,6) + &
                          q1(7,8) * var(:,:,7) + q1(8,8) * var(:,:,8) + &
                          q1(9,8) * var(:,:,9) + q1(10,8) * var(:,:,10) + &
                          q1(11,8) * var(:,:,11) + &
                          q1(12,8) * var(:,:,12) ) * idel
        elsewhere
          dvar(:,:,8) = ( q2(1,8) * var(:,:,1) + q2(2,8) * var(:,:,2) + &
                          q2(3,8) * var(:,:,3) + q2(4,8) * var(:,:,4) + &
                          q2(5,8) * var(:,:,5) + q2(6,8) * var(:,:,6) + &
                          q2(7,8) * var(:,:,7) + q2(8,8) * var(:,:,8) + &
                          q2(9,8) * var(:,:,9) + q2(10,8) * var(:,:,10) + &
                          q2(11,8) * var(:,:,11) + &
                          q2(12,8) * var(:,:,12) ) * idel
        end where
        kl = 9
      end if
      if ( bb(2) == 0 ) then
        kr = nk - gsize
      else 
        where ( up(:,:,nk) < zero )
          dvar(:,:,nk) =   - ( q2(1,1) * var(:,:,nk) + &
                               q2(2,1) * var(:,:,nk-1) + &
                               q2(3,1) * var(:,:,nk-2) + &
                               q2(4,1) * var(:,:,nk-3) + &
                               q2(5,1) * var(:,:,nk-4) + &
                               q2(6,1) * var(:,:,nk-5) + &
                               q2(7,1) * var(:,:,nk-6) + &
                               q2(8,1) * var(:,:,nk-7) ) * idel
        elsewhere
          dvar(:,:,nk) =   - ( q1(1,1) * var(:,:,nk) + &
                               q1(2,1) * var(:,:,nk-1) + &
                               q1(3,1) * var(:,:,nk-2) + &
                               q1(4,1) * var(:,:,nk-3) + &
                               q1(5,1) * var(:,:,nk-4) + &
                               q1(6,1) * var(:,:,nk-5) + &
                               q1(7,1) * var(:,:,nk-6) + &
                               q1(8,1) * var(:,:,nk-7) ) * idel
        end where
        where ( up(:,:,nk-1) < zero )
          dvar(:,:,nk-1) = - ( q2(1,2) * var(:,:,nk) + &
                               q2(2,2) * var(:,:,nk-1) + &
                               q2(3,2) * var(:,:,nk-2) + &
                               q2(4,2) * var(:,:,nk-3) + &
                               q2(5,2) * var(:,:,nk-4) + &
                               q2(6,2) * var(:,:,nk-5) + &
                               q2(7,2) * var(:,:,nk-6) + &
                               q2(8,2) * var(:,:,nk-7) ) * idel
        elsewhere
          dvar(:,:,nk-1) = - ( q1(1,2) * var(:,:,nk) + &
                               q1(2,2) * var(:,:,nk-1) + &
                               q1(3,2) * var(:,:,nk-2) + &
                               q1(4,2) * var(:,:,nk-3) + &
                               q1(5,2) * var(:,:,nk-4) + &
                               q1(6,2) * var(:,:,nk-5) + &
                               q1(7,2) * var(:,:,nk-6) + &
                               q1(8,2) * var(:,:,nk-7) ) * idel
        end where
        where ( up(:,:,nk-2) < zero )
          dvar(:,:,nk-2) = - ( q2(1,3) * var(:,:,nk) + &
                               q2(2,3) * var(:,:,nk-1) + &
                               q2(3,3) * var(:,:,nk-2) + &
                               q2(4,3) * var(:,:,nk-3) + &
                               q2(5,3) * var(:,:,nk-4) + &
                               q2(6,3) * var(:,:,nk-5) + &
                               q2(7,3) * var(:,:,nk-6) + &
                               q2(8,3) * var(:,:,nk-7) ) * idel
        elsewhere
          dvar(:,:,nk-2) = - ( q1(1,3) * var(:,:,nk) + &
                               q1(2,3) * var(:,:,nk-1) + &
                               q1(3,3) * var(:,:,nk-2) + &
                               q1(4,3) * var(:,:,nk-3) + &
                               q1(5,3) * var(:,:,nk-4) + &
                               q1(6,3) * var(:,:,nk-5) + &
                               q1(7,3) * var(:,:,nk-6) + &
                               q1(8,3) * var(:,:,nk-7) ) * idel
        end where
        where ( up(:,:,nk-3) < zero )
          dvar(:,:,nk-3) = - ( q2(1,4) * var(:,:,nk) + &
                               q2(2,4) * var(:,:,nk-1) + &
                               q2(3,4) * var(:,:,nk-2) + &
                               q2(4,4) * var(:,:,nk-3) + &
                               q2(5,4) * var(:,:,nk-4) + &
                               q2(6,4) * var(:,:,nk-5) + &
                               q2(7,4) * var(:,:,nk-6) + &
                               q2(8,4) * var(:,:,nk-7) ) * idel
        elsewhere
          dvar(:,:,nk-3) = - ( q1(1,4) * var(:,:,nk) + &
                               q1(2,4) * var(:,:,nk-1) + &
                               q1(3,4) * var(:,:,nk-2) + &
                               q1(4,4) * var(:,:,nk-3) + &
                               q1(5,4) * var(:,:,nk-4) + &
                               q1(6,4) * var(:,:,nk-5) + &
                               q1(7,4) * var(:,:,nk-6) + &
                               q1(8,4) * var(:,:,nk-7) ) * idel
        end where
        where ( up(:,:,nk-4) < zero )
          dvar(:,:,nk-4) = - ( q2(1,5) * var(:,:,nk) + &
                               q2(2,5) * var(:,:,nk-1) + &
                               q2(3,5) * var(:,:,nk-2) + &
                               q2(4,5) * var(:,:,nk-3) + &
                               q2(5,5) * var(:,:,nk-4) + &
                               q2(6,5) * var(:,:,nk-5) + &
                               q2(7,5) * var(:,:,nk-6) + &
                               q2(8,5) * var(:,:,nk-7) + &
                               q2(9,5) * var(:,:,nk-8) ) * idel
        elsewhere
          dvar(:,:,nk-4) = - ( q1(1,5) * var(:,:,nk) + &
                               q1(2,5) * var(:,:,nk-1) + &
                               q1(3,5) * var(:,:,nk-2) + &
                               q1(4,5) * var(:,:,nk-3) + &
                               q1(5,5) * var(:,:,nk-4) + &
                               q1(6,5) * var(:,:,nk-5) + &
                               q1(7,5) * var(:,:,nk-6) + &
                               q1(8,5) * var(:,:,nk-7) + &
                               q1(9,5) * var(:,:,nk-8) ) * idel
        end where
        where ( up(:,:,nk-5) < zero )
          dvar(:,:,nk-5) = - ( q2(1,6) * var(:,:,nk) + &
                               q2(2,6) * var(:,:,nk-1) + &
                               q2(3,6) * var(:,:,nk-2) + &
                               q2(4,6) * var(:,:,nk-3) + &
                               q2(5,6) * var(:,:,nk-4) + &
                               q2(6,6) * var(:,:,nk-5) + &
                               q2(7,6) * var(:,:,nk-6) + &
                               q2(8,6) * var(:,:,nk-7) + &
                               q2(9,6) * var(:,:,nk-8) + &
                               q2(10,6) * var(:,:,nk-9) ) * idel
        elsewhere
          dvar(:,:,nk-5) = - ( q1(1,6) * var(:,:,nk) + &
                               q1(2,6) * var(:,:,nk-1) + &
                               q1(3,6) * var(:,:,nk-2) + &
                               q1(4,6) * var(:,:,nk-3) + &
                               q1(5,6) * var(:,:,nk-4) + &
                               q1(6,6) * var(:,:,nk-5) + &
                               q1(7,6) * var(:,:,nk-6) + &
                               q1(8,6) * var(:,:,nk-7) + &
                               q1(9,6) * var(:,:,nk-8) + &
                               q1(10,6) * var(:,:,nk-9) ) * idel
        end where
        where ( up(:,:,nk-6) < zero )
          dvar(:,:,nk-6) = - ( q2(1,7) * var(:,:,nk) + &
                               q2(2,7) * var(:,:,nk-1) + &
                               q2(3,7) * var(:,:,nk-2) + &
                               q2(4,7) * var(:,:,nk-3) + &
                               q2(5,7) * var(:,:,nk-4) + &
                               q2(6,7) * var(:,:,nk-5) + &
                               q2(7,7) * var(:,:,nk-6) + &
                               q2(8,7) * var(:,:,nk-7) + &
                               q2(9,7) * var(:,:,nk-8) + &
                               q2(10,7) * var(:,:,nk-9) + &
                               q2(11,7) * var(:,:,nk-10) ) * idel
        elsewhere
          dvar(:,:,nk-6) = - ( q1(1,7) * var(:,:,nk) + &
                               q1(2,7) * var(:,:,nk-1) + &
                               q1(3,7) * var(:,:,nk-2) + &
                               q1(4,7) * var(:,:,nk-3) + &
                               q1(5,7) * var(:,:,nk-4) + &
                               q1(6,7) * var(:,:,nk-5) + &
                               q1(7,7) * var(:,:,nk-6) + &
                               q1(8,7) * var(:,:,nk-7) + &
                               q1(9,7) * var(:,:,nk-8) + &
                               q1(10,7) * var(:,:,nk-9) + &
                               q1(11,7) * var(:,:,nk-10) ) * idel
        end where
        where ( up(:,:,nk-7) < zero )
          dvar(:,:,nk-7) = - ( q2(1,8) * var(:,:,nk) + &
                               q2(2,8) * var(:,:,nk-1) + &
                               q2(3,8) * var(:,:,nk-2) + &
                               q2(4,8) * var(:,:,nk-3) + &
                               q2(5,8) * var(:,:,nk-4) + &
                               q2(6,8) * var(:,:,nk-5) + &
                               q2(7,8) * var(:,:,nk-6) + &
                               q2(8,8) * var(:,:,nk-7) + &
                               q2(9,8) * var(:,:,nk-8) + &
                               q2(10,8) * var(:,:,nk-9) + &
                               q2(11,8) * var(:,:,nk-10) + &
                               q2(12,8) * var(:,:,nk-11) ) * idel
        elsewhere
          dvar(:,:,nk-7) = - ( q1(1,8) * var(:,:,nk) + &
                               q1(2,8) * var(:,:,nk-1) + &
                               q1(3,8) * var(:,:,nk-2) + &
                               q1(4,8) * var(:,:,nk-3) + &
                               q1(5,8) * var(:,:,nk-4) + &
                               q1(6,8) * var(:,:,nk-5) + &
                               q1(7,8) * var(:,:,nk-6) + &
                               q1(8,8) * var(:,:,nk-7) + &
                               q1(9,8) * var(:,:,nk-8) + &
                               q1(10,8) * var(:,:,nk-9) + &
                               q1(11,8) * var(:,:,nk-10) + &
                               q1(12,8) * var(:,:,nk-11) ) * idel
        end where
        kr = nk - 8
      end if
      if (kl > kr+1) call CCTK_WARN (0, "domain too small")
      where ( up(:,:,kl:kr) < zero )
        dvar(:,:,kl:kr) = ( a1(-4) * var(:,:,kl-4:kr-4) + &
                            a1(-3) * var(:,:,kl-3:kr-3) + &
                            a1(-2) * var(:,:,kl-2:kr-2) + &
                            a1(-1) * var(:,:,kl-1:kr-1) + &
                            a1(0)  * var(:,:,kl:kr) + &
                            a1(1)  * var(:,:,kl+1:kr+1) + &
                            a1(2)  * var(:,:,kl+2:kr+2) + &
                            a1(3)  * var(:,:,kl+3:kr+3) + &
                            a1(4)  * var(:,:,kl+4:kr+4) ) * idel
      elsewhere
        dvar(:,:,kl:kr) = ( a2(-4) * var(:,:,kl-4:kr-4) + &
                            a2(-3) * var(:,:,kl-3:kr-3) + &
                            a2(-2) * var(:,:,kl-2:kr-2) + &
                            a2(-1) * var(:,:,kl-1:kr-1) + &
                            a2(0)  * var(:,:,kl:kr) + &
                            a2(1)  * var(:,:,kl+1:kr+1) + &
                            a2(2)  * var(:,:,kl+2:kr+2) + &
                            a2(3)  * var(:,:,kl+3:kr+3) + &
                            a2(4)  * var(:,:,kl+4:kr+4) ) * idel
      end where
    end if
  end select direction
end subroutine up_deriv_gf_8_4_opt
