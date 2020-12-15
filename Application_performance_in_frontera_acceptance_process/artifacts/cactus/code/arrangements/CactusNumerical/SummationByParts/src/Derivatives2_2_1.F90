#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine deriv2_gf_2_1 ( var, ni, nj, nk, dir, bb, gsize, delta, dvar2 )

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
  CCTK_REAL, dimension(ni,nj,nk), intent(OUT) :: dvar2

  CCTK_REAL, dimension(2), save :: a
  CCTK_REAL, dimension(3,2), save :: q 
  CCTK_REAL :: idel2

  CCTK_INT :: il, ir, jl, jr, kl, kr

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_2_2_1 ( a, q )
    first = .false.
  end if

  idel2 = 1.0_wp / delta**2

  direction: select case (dir)
  case (0) direction
    if ( bb(1) == 0 ) then
      il = 1 + gsize
    else
      dvar2(1,:,:) = ( q(1,1) * var(1,:,:) + q(2,1) * var(2,:,:) + &
                       q(3,1) * var(3,:,:) ) * idel2
      dvar2(2,:,:) = ( q(1,2) * var(1,:,:) + q(2,2) * var(2,:,:) + &
                       q(3,2) * var(3,:,:) ) * idel2
      il = 3
    end if
    if ( bb(2) == 0 ) then
      ir = ni - gsize
    else
      dvar2(ni-1,:,:) = ( q(1,2) * var(ni,:,:) + &
                          q(2,2) * var(ni-1,:,:) + &
                          q(3,2) * var(ni-2,:,:) ) * idel2
      dvar2(ni,:,:) = ( q(1,1) * var(ni,:,:) + &
                        q(2,1) * var(ni-1,:,:) + &
                        q(3,1) * var(ni-2,:,:) ) * idel2
      ir = ni - 2
    end if
    dvar2(il:ir,:,:) = ( a(2) * ( var(il+1:ir+1,:,:) + &
                                  var(il-1:ir-1,:,:) ) + &
                         a(1) * var(il:ir,:,:) ) * idel2
  case (1) direction
    if ( zero_derivs_y /= 0 ) then
      dvar2 = zero
    else
      if ( bb(1) == 0 ) then
        jl = 1 + gsize
      else
        dvar2(:,1,:) = ( q(1,1) * var(:,1,:) + q(2,1) * var(:,2,:) + &
                         q(3,1) * var(:,3,:) ) * idel2
        dvar2(:,2,:) = ( q(1,2) * var(:,1,:) + q(2,2) * var(:,2,:) + &
                         q(3,2) * var(:,3,:) ) * idel2
        jl = 3
      end if
      if ( bb(2) == 0 ) then
        jr = nj - gsize
      else
        dvar2(:,nj-1,:) = ( q(1,2) * var(:,nj,:) + &
                            q(2,2) * var(:,nj-1,:) + &
                            q(3,2) * var(:,nj-2,:) ) * idel2
        dvar2(:,nj,:) = ( q(1,1) * var(:,nj,:) + &
                          q(2,1) * var(:,nj-1,:) + &
                          q(3,1) * var(:,nj-2,:) ) * idel2
        jr = nj - 2
      end if
      dvar2(:,jl:jr,:) = ( a(2) * ( var(:,jl+1:jr+1,:) + &
                                    var(:,jl-1:jr-1,:) ) + &
                           a(1) * var(:,jl:jr,:) ) * idel2
    end if
  case (2) direction
    if ( zero_derivs_z /= 0 ) then
      dvar2 = zero
    else
      if ( bb(1) == 0 ) then
        kl = 1 + gsize
      else
        dvar2(:,:,1) = ( q(1,1) * var(:,:,1) + q(2,1) * var(:,:,2) + &
                         q(3,1) * var(:,:,3) ) * idel2
        dvar2(:,:,2) = ( q(1,2) * var(:,:,1) + q(2,2) * var(:,:,2) + &
                         q(3,2) * var(:,:,3) ) * idel2
        kl = 3
      end if
      if ( bb(2) == 0 ) then
        kr = nk - gsize
      else 
        dvar2(:,:,nk-1) = ( q(1,2) * var(:,:,nk) + &
                            q(2,2) * var(:,:,nk-1) + &
                            q(3,2) * var(:,:,nk-2) ) * idel2
        dvar2(:,:,nk) = ( q(1,1) * var(:,:,nk) + &
                          q(2,1) * var(:,:,nk-1) + &
                          q(3,1) * var(:,:,nk-2) ) * idel2
        kr = nk - 2
      end if
      dvar2(:,:,kl:kr) = ( a(2) * ( var(:,:,kl+1:kr+1) +  &
                                    var(:,:,kl-1:kr-1) ) + &
                           a(1) * var(:,:,kl:kr) ) * idel2
    end if
  end select direction
end subroutine deriv2_gf_2_1
