#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "GRHydro_Macros.h"


subroutine Con2Prim_pt_bisect(cctk_iteration,ii,jj,kk,&
     handle, dens, sx, sy, sz, tau, rho, velx, vely, &
     velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, sdetg, x, y, z, r, epsnegative, GRHydro_rho_min, pmin, epsmin, &
     GRHydro_reflevel, GRHydro_C2P_failed, local_perc_ptol, myproc)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: pmax0 = 1.0d0

  CCTK_REAL dens, sx, sy, sz, tau, rho, velx, vely, velz, epsilon, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, sdetg, isdetg, w_lorentz, x, &
       y, z, r, GRHydro_rho_min, local_perc_ptol
  
  CCTK_REAL :: s2, f, df, vlowx, vlowy, vlowz
  CCTK_INT  :: count, handle, GRHydro_reflevel, myproc
  CCTK_REAL :: GRHydro_C2P_failed
  CCTK_REAL :: udens, usx, usy, usz, utau, pold, pnew, &
            dummy1,dummy2,pmin, epsmin

  CCTK_REAL :: p1,p2,tmp1,tmp2,tmp3
  CCTK_REAL :: pin
  CCTK_REAL :: pmax
  logical :: done

  CCTK_INT cctk_iteration,ii,jj,kk
  character(len=200) warnline
  logical epsnegative

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xtemp,xye,tmp,plow
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars


!!$  Undensitize the variables 

  isdetg = 1.0d0/sdetg
  udens = dens * isdetg
  usx = sx * isdetg
  usy = sy * isdetg
  usz = sz * isdetg
  utau = tau * isdetg
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.0d0*usx*usz*uxz + 2.0d0*usy*usz*uyz  

  pin = press
!!$  Set initial guess for pressure:
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
                      rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)
  !pold = max(1.d-10,xpress)
  ! This is the lowest admissible pressure, otherwise we are off the physical regime
  ! Sometimes, we may end up with bad initial guesses. In that case we
  ! need to start off with something reasonable at least
  plow = max(pmin, sqrt(s2) - utau - udens)
  
  ! Start out with some reasonable initial guess
  pold = max(plow+1.d-10,xpress)
  pmax = min(1.0d6*pold,pmax0)

!!$  Check that the variables have a chance of being physical
  if( (utau + pold + udens)**2 - s2 .le. 0.0d0) then
    if (c2p_reset_pressure .ne. 0) then
      pold = sqrt(s2 + c2p_reset_pressure_to_value) - utau - udens
    else 
      GRHydro_C2P_failed = 1
    endif
  endif
  
!!$  Calculate rho and epsilon 

!define temporary variables to speed up
  rho = udens * sqrt( (utau + pold + udens)**2 - s2)/(utau + pold + udens)
  w_lorentz = (utau + pold + udens) / sqrt( (utau + pold + udens)**2 - s2)
  epsilon = (sqrt( (utau + pold + udens)**2 - s2) - pold * w_lorentz - &
       udens)/udens
  
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
                      rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)

  tmp1 = (xpress - pold)  !f(c)
  tmp2 = (xpress - plow)  !f(a)
  tmp3 = (xpress - pmax)  !f(b)

  done = .false.

  ! make sure we are bracketing the root
  if(tmp2*tmp3 .gt. 0.0d0) then
     ! In this case, we are not bracketing the root.
     ! Based on what we do with pold / plow / pmin above,
     ! the only thing that can be is that pmax is too low.
     ! -> make pmax bigger:
     pmax = pmax0
  endif


  if (tmp1*tmp2 < 0.0d0) then
     pnew = (plow + pold)/2.0d0
     p1 = plow
     p2 = pold
  else
     pnew = (pold + pmax)/2.0d0
     p1 = pold
     p2 = pmax
  endif

!!$Find the root
  count = 0

  do while(.not.done.or.count.lt.GRHydro_countmin)

    count = count + 1

    if (count > GRHydro_countmax*2) then

      !$OMP CRITICAL
      call CCTK_WARN(1, 'count > 2 * GRHydro_countmax! ')
      write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a28,i8)') 'cctk_iteration:', cctk_iteration
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3i5,3g16.7)') 'ijk, xyz location: ',ii,jj,kk,x,y,z
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'radius: ',r
      call CCTK_WARN(1,warnline)
      write(warnline,*) "uxx", uxx
      call CCTK_WARN(1,warnline)
      write(warnline,*) "uxy", uxy
      call CCTK_WARN(1,warnline)
      write(warnline,*) "uxz", uxz
      call CCTK_WARN(1,warnline)
      write(warnline,*) "uyy", uyy
      call CCTK_WARN(1,warnline)
      write(warnline,*) "uyz", uyz
      call CCTK_WARN(1,warnline)
      write(warnline,*) "uzz", uzz
      call CCTK_WARN(1,warnline)
      write(warnline,*) "scon", sx, sy, sz
      call CCTK_WARN(1,warnline)
      write(warnline,*) "dens,tau", dens, tau
      call CCTK_WARN(1,warnline)
      write(warnline,*) "rho,eps", rho,epsilon
      call CCTK_WARN(1,warnline)
      write(warnline,*) "sdetg", sdetg
      call CCTK_WARN(1,warnline)
      write(warnline,'(a21,1P10E15.6)') 'udens, utau, pnew, s2', udens,utau,pnew,s2
      call CCTK_WARN(1,warnline)
      write(warnline,'(1P10E15.6)') (abs(pnew - pold)/abs(pnew))
      call CCTK_WARN(1,warnline)
      ! We will also accept the root if we can't do any better, 
      ! namely if p2-p1 < 1.0d-14 (round off)
!      if(abs(p2-p1).le.1.0d-14*max(p2,p1)) then
!          if(myproc.eq.9) then
      call CCTK_WARN(0,"Not accepting root (cannot do better)")
      !$OMP END CRITICAL

      ! for safety, let's set the point to atmosphere
      SET_ATMO_MIN(rho, GRHydro_rho_min, r)
      udens = rho
      dens = sdetg * rho
      pnew = pmin
      epsilon = epsmin
      ! w_lorentz=1, so the expression for utau reduces to:
      utau  = rho + rho*epsmin - udens
      sx = 0.d0
      sy = 0.d0
      sz = 0.d0
      s2 = 0.d0
      usx = 0.d0
      usy = 0.d0
      usz = 0.d0
      w_lorentz = 1.d0
      goto 51
    end if

    
!!$    Recalculate primitive variables and function
    tmp = (utau + pnew + udens)**2 - s2
    if (tmp .le. 0.0d0) then
       call CCTK_WARN(0,"Bad bisect c2p!")
    endif

    rho = udens * sqrt( tmp)/(utau + pnew + udens)
    w_lorentz = (utau + pnew + udens) / sqrt( tmp)
    epsilon = (sqrt( tmp) - pnew * w_lorentz - &
         udens)/udens

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,xtemp,xye,xpress,keyerr,anyerr)
    
    tmp1 = xpress - pnew
    tmp2 = xpress - p1
    tmp3 = xpress - p2

    if (abs(tmp1/pnew).lt.local_perc_ptol) then
       done = .true.
    else if (tmp2*tmp3 > 0.0d0) then
       ! In this case we are no longer bracketing the loop.
       ! This can happen due to round-off error at critical points, e.g.,
       ! values of tmp that are very small.
       ! We will accept the root as long as it is within an order
       ! of magnitude of where we want to be.
       if(abs(tmp1/pnew).lt.10.0d0*local_perc_ptol) then
          done = .true.
          !$OMP CRITICAL
          write(warnline,"(A7,1P10E22.14)") "BISECT:", abs(tmp1/pnew)
          call CCTK_WARN(1,warnline)
          write(warnline,*) "Not bracketing root, but accepting result, since < 10 x GRHydro_perc_ptol"
          call CCTK_WARN(1,warnline)
          !$OMP END CRITICAL
       endif
    endif
    
    if(.not.done) then
       if (tmp1*tmp2 < 0.0d0) then
          p1 = p1
          p2 = pnew
          pnew = (p1 + p2)*0.5d0
       else
          p1 = pnew
          p2 = p2
          pnew = (p1 + p2)*0.5d0
       endif
    endif

  enddo
  
!!$  Calculate primitive variables from root

  !if (rho .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
  IF_BELOW_ATMO(rho, GRHydro_rho_min, GRHydro_atmo_tolerance, r) then
    SET_ATMO_MIN(rho, GRHydro_rho_min, r) !GRHydro_rho_min
    udens = rho
    dens = sdetg * rho
!    epsilon = (sqrt( (utau + pnew + udens)**2) - pnew -  udens)/udens
    epsilon = epsmin
    ! w_lorentz=1, so the expression for utau reduces to:
    utau  = rho + rho*epsmin - udens
    sx = 0.d0
    sy = 0.d0
    sz = 0.d0
    s2 = 0.d0
    usx = 0.d0
    usy = 0.d0
    usz = 0.d0
    w_lorentz = 1.d0
  end if

51 press = pnew
  vlowx = usx / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowy = usy / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowz = usz / ( (rho + rho*epsilon + press) * w_lorentz**2)
  velx = uxx * vlowx + uxy * vlowy + uxz * vlowz
  vely = uxy * vlowx + uyy * vlowy + uyz * vlowz
  velz = uxz * vlowx + uyz * vlowy + uzz * vlowz
  
    
!!$If all else fails, use the polytropic EoS

  if(epsilon .lt. 0.0d0) then
    epsnegative = .true.
  endif

end subroutine Con2Prim_pt_bisect


