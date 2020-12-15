#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine deriv_gf_2_1 ( var, ni, nj, nk, dir, bb, gsize, offset, delta, dvar )

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

  CCTK_REAL, dimension(1), save :: a
  CCTK_REAL, dimension(3,2), save :: q 
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_1_2_1 ( a, q )
    first = .false.
  end if

  idel = 1.0_wp / delta

  if (gsize < 1) call CCTK_WARN (0, "not enough ghostzones")

  direction: select case (dir)
  case (0) direction
    if ( bb(1) == 0 ) then
      il = 1 + gsize
    else
      ol = offset(1)
!$omp parallel workshare
      dvar(1+ol,:,:) = ( q(1,1) * var(1+ol,:,:) + q(2,1) * var(2+ol,:,:) ) * idel
      dvar(2+ol,:,:) = ( q(1,2) * var(1+ol,:,:) + q(3,2) * var(3+ol,:,:) ) * idel
!$omp end parallel workshare
      il = 3 + ol
    end if
    if ( bb(2) == 0 ) then
      ir = ni - gsize
    else
      or = ni - offset(2)
!$omp parallel workshare
      dvar(or-1,:,:) = - ( q(1,2) * var(or,:,:) + &
                           q(3,2) * var(or-2,:,:) ) * idel
      dvar(or,:,:) = - ( q(1,1) * var(or,:,:) + &
                         q(2,1) * var(or-1,:,:) ) * idel
!$omp end parallel workshare
      ir = or - 2
    end if
    if (il > ir+1) call CCTK_WARN (0, "domain too small")
!$omp parallel workshare
    dvar(il:ir,:,:) = a(1) * ( var(il+1:ir+1,:,:) - var(il-1:ir-1,:,:) ) * idel
!$omp end parallel workshare
  case (1) direction
    if ( zero_derivs_y /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        jl = 1 + gsize
      else
        ol = offset(1)
!$omp parallel workshare
        dvar(:,1+ol,:) = ( q(1,1) * var(:,1+ol,:) + q(2,1) * var(:,2+ol,:) ) * idel
        dvar(:,2+ol,:) = ( q(1,2) * var(:,1+ol,:) + q(3,2) * var(:,3+ol,:) ) * idel
!$omp end parallel workshare
        jl = 3 + ol
      end if
      if ( bb(2) == 0 ) then
        jr = nj - gsize
      else
        or = nj - offset(2)
!$omp parallel workshare
        dvar(:,or-1,:) = - ( q(1,2) * var(:,or,:) + &
                             q(3,2) * var(:,or-2,:) ) * idel
        dvar(:,or,:) = - ( q(1,1) * var(:,or,:) + &
                           q(2,1) * var(:,or-1,:) ) * idel
!$omp end parallel workshare
        jr = or - 2
      end if
      if (jl > jr+1) call CCTK_WARN (0, "domain too small")
!$omp parallel workshare
      dvar(:,jl:jr,:) = a(1) * ( var(:,jl+1:jr+1,:) - &
                                 var(:,jl-1:jr-1,:) ) * idel
!$omp end parallel workshare
    end if
  case (2) direction
    if ( zero_derivs_z /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        kl = 1 + gsize
      else
        ol = offset(1)
!$omp parallel workshare
        dvar(:,:,1+ol) = ( q(1,1) * var(:,:,1+ol) + q(2,1) * var(:,:,2+ol) ) * idel
        dvar(:,:,2+ol) = ( q(1,2) * var(:,:,1+ol) + q(3,2) * var(:,:,3+ol) ) * idel
!$omp end parallel workshare
        kl = 3 + ol
      end if
      if ( bb(2) == 0 ) then
        kr = nk - gsize
      else 
        or = nk - offset(2)
!$omp parallel workshare
        dvar(:,:,or-1) = - ( q(1,2) * var(:,:,or) + &
                             q(3,2) * var(:,:,or-2) ) * idel
        dvar(:,:,or) = - ( q(1,1) * var(:,:,or) + &
                           q(2,1) * var(:,:,or-1) ) * idel
!$omp end parallel workshare
        kr = or - 2
      end if
      if (kl > kr+1) call CCTK_WARN (0, "domain too small")
!$omp parallel workshare
      dvar(:,:,kl:kr) = a(1) * ( var(:,:,kl+1:kr+1) -  &
                                 var(:,:,kl-1:kr-1) ) * idel
!$omp end parallel workshare
    end if
  end select direction
end subroutine deriv_gf_2_1

subroutine up_deriv_gf_2_1 ( var, ni, nj, nk, dir, bb, gsize, delta, up, dvar )

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

  CCTK_REAL, dimension(-1:1), save :: a1, a2
  CCTK_REAL, dimension(3,2), save :: q1, q2 
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr

  logical, save :: first = .true.

  if ( first ) then
    call coeffs_up_2_1 ( a1, q1, a2, q2 )
    first = .false.
  end if

  idel = 1.0_wp / delta

  if (gsize < 1) call CCTK_WARN (0, "not enough ghostzones")

  direction: select case (dir)
  case (0) direction
    if ( bb(1) == 0 ) then
      il = 1 + gsize
    else
      where ( up(1,:,:) < zero )
        dvar(1,:,:) = ( q1(1,1) * var(1,:,:) + q1(2,1) * var(2,:,:) ) * idel
      elsewhere
        dvar(1,:,:) = ( q2(1,1) * var(1,:,:) + q2(2,1) * var(2,:,:) ) * idel
      end where
      where ( up(2,:,:) < zero ) 
        dvar(2,:,:) = ( q1(1,2) * var(1,:,:) + q1(2,2) * var(2,:,:) + &
                        q1(3,2) * var(3,:,:) ) * idel
      elsewhere
        dvar(2,:,:) = ( q2(1,2) * var(1,:,:) + q2(2,2) * var(2,:,:) + &
                        q2(3,2) * var(3,:,:) ) * idel
      end where
      il = 3
    end if
    if ( bb(2) == 0 ) then
      ir = ni - gsize
    else
      where ( up(ni-1,:,:) < zero )
        dvar(ni-1,:,:) = - ( q2(1,2) * var(ni,:,:) + &
                             q2(2,2) * var(ni-1,:,:) + &
                             q2(3,2) * var(ni-2,:,:) ) * idel
      elsewhere
        dvar(ni-1,:,:) = - ( q1(1,2) * var(ni,:,:) + &
                             q1(2,2) * var(ni-1,:,:) + &
                             q1(3,2) * var(ni-2,:,:) ) * idel
      end where
      where ( up(ni,:,:) < zero )
        dvar(ni,:,:) = - ( q2(1,1) * var(ni,:,:) + &
                           q2(2,1) * var(ni-1,:,:) ) * idel
      elsewhere
        dvar(ni,:,:) = - ( q1(1,1) * var(ni,:,:) + &
                           q1(2,1) * var(ni-1,:,:) ) * idel
      end where
      ir = ni - 2
    end if
    if (il > ir+1) call CCTK_WARN (0, "domain too small")
    where ( up(il:ir,:,:) < zero ) 
      dvar(il:ir,:,:) = ( a1(-1) * var(il-1:ir-1,:,:) + &
                          a1(0)  * var(il:ir,:,:) + &
                          a1(1)  * var(il+1:ir+1,:,:) ) * idel
    elsewhere
      dvar(il:ir,:,:) = ( a2(-1) * var(il-1:ir-1,:,:) + &
                          a2(0)  * var(il:ir,:,:) + &
                          a2(1)  * var(il+1:ir+1,:,:) ) * idel
    end where
  case (1) direction
    if ( zero_derivs_y /= 0 ) then
      dvar = zero
    else
      if ( bb(1) == 0 ) then
        jl = 1 + gsize
      else
        where ( up(:,1,:) < zero )
          dvar(:,1,:) = ( q1(1,1) * var(:,1,:) + q1(2,1) * var(:,2,:) ) * idel
        elsewhere
          dvar(:,1,:) = ( q2(1,1) * var(:,1,:) + q2(2,1) * var(:,2,:) ) * idel
        end where
        where ( up(:,2,:) < zero )
          dvar(:,2,:) = ( q1(1,2) * var(:,1,:) + q1(2,2) * var(:,2,:) + &
                          q1(3,2) * var(:,3,:) ) * idel
        elsewhere
          dvar(:,2,:) = ( q2(1,2) * var(:,1,:) + q2(2,2) * var(:,2,:) + &
                          q2(3,2) * var(:,3,:) ) * idel
        end where
        jl = 3
      end if
      if ( bb(2) == 0 ) then
        jr = nj - gsize
      else
        where ( up(:,nj-1,:) < zero )
          dvar(:,nj-1,:) = - ( q2(1,2) * var(:,nj,:) + &
                               q2(2,2) * var(:,nj-1,:) + &
                               q2(3,2) * var(:,nj-2,:) ) * idel
        elsewhere
          dvar(:,nj-1,:) = - ( q1(1,2) * var(:,nj,:) + &
                               q1(2,2) * var(:,nj-1,:) + &
                               q1(3,2) * var(:,nj-2,:) ) * idel
        end where
        where ( up(:,nj,:) < zero )
          dvar(:,nj,:) = - ( q2(1,1) * var(:,nj,:) + &
                             q2(2,1) * var(:,nj-1,:) ) * idel
        elsewhere
          dvar(:,nj,:) = - ( q1(1,1) * var(:,nj,:) + &
                             q1(2,1) * var(:,nj-1,:) ) * idel
        end where
        jr = nj - 2
      end if
      if (jl > jr+1) call CCTK_WARN (0, "domain too small")
      where ( up(:,jl:jr,:) < zero ) 
        dvar(:,jl:jr,:) = ( a1(-1) * var(:,jl-1:jr-1,:) + &
                            a1(0)  * var(:,jl:jr,:) + &
                            a1(1)  * var(:,jl+1:jr+1,:) ) * idel
      elsewhere
        dvar(:,jl:jr,:) = ( a2(-1) * var(:,jl-1:jr-1,:) + &
                            a2(0)  * var(:,jl:jr,:) + &
                            a2(1)  * var(:,jl+1:jr+1,:) ) * idel
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
          dvar(:,:,1) = ( q1(1,1) * var(:,:,1) + q1(2,1) * var(:,:,2) ) * idel
        elsewhere
          dvar(:,:,1) = ( q2(1,1) * var(:,:,1) + q2(2,1) * var(:,:,2) ) * idel
        end where
        where ( up(:,:,2) < zero )  
          dvar(:,:,2) = ( q1(1,2) * var(:,:,1) + q1(2,2) * var(:,:,2) + &
                          q1(3,2) * var(:,:,3) ) * idel
        elsewhere
          dvar(:,:,2) = ( q2(1,2) * var(:,:,1) + q2(2,2) * var(:,:,2) + &
                          q2(3,2) * var(:,:,3) ) * idel
        end where
        kl = 3
      end if
      if ( bb(2) == 0 ) then
        kr = nk - gsize
      else 
        where ( up(:,:,nk-1) < zero )
          dvar(:,:,nk-1) = - ( q2(1,2) * var(:,:,nk) + &
                               q2(2,2) * var(:,:,nk-1) + &
                               q2(3,2) * var(:,:,nk-2) ) * idel
        elsewhere
          dvar(:,:,nk-1) = - ( q1(1,2) * var(:,:,nk) + &
                               q1(2,2) * var(:,:,nk-1) + &
                               q1(3,2) * var(:,:,nk-2) ) * idel
        end where
        where ( up(:,:,nk) < zero )
          dvar(:,:,nk) = - ( q2(1,1) * var(:,:,nk) + &
                             q2(2,1) * var(:,:,nk-1) ) * idel
        elsewhere
          dvar(:,:,nk) = - ( q1(1,1) * var(:,:,nk) + &
                             q1(2,1) * var(:,:,nk-1) ) * idel
        end where
        kr = nk - 2
      end if
      if (kl > kr+1) call CCTK_WARN (0, "domain too small")
      where ( up(:,:,kl:kr) < zero )
        dvar(:,:,kl:kr) = ( a1(-1) * var(:,:,kl-1:kr-1) + &
                            a1(0)  * var(:,:,kl:kr) + &
                            a1(1)  * var(:,:,kl+1:kr+1) ) * idel
      elsewhere
        dvar(:,:,kl:kr) = ( a2(-1) * var(:,:,kl-1:kr-1) + &
                            a2(0)  * var(:,:,kl:kr) + &
                            a2(1)  * var(:,:,kl+1:kr+1) ) * idel
      end where
    end if
  end select direction
end subroutine up_deriv_gf_2_1
