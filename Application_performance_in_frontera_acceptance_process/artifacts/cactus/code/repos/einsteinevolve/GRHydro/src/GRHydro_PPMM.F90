 /*@@
   @file      GRHydro_PPMM.F90
   @date      Sun Feb 10 16:53:29 2002
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke, Toni Font, Luca Baiotti, Frank Loeffler
   @desc 
   Routines to do PPM reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

 /*@@
   @routine    SimplePPM_1dM
   @date       Thu Feb 14 19:08:52 2002
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke, Toni Font
   @desc 
   The simple PPM reconstruction routine that applies along
   each one dimensional slice.

   @enddesc 
   @calls     
   @calledby   
   @history 
   Written in frustration when IH couldn''t get Toni''s original code 
   to work.
   @endhistory 

@@*/

subroutine SimplePPM_1dM(handle,poly,nx,dx,&
     rho,velx,vely,velz,Bvcx,Bvcy,Bvcz,psidc,eps,press,&
     rhominus,velxminus,velyminus,velzminus,Bvcxminus,Bvcyminus,Bvczminus,psidcminus,epsminus,&
     rhoplus,velxplus,velyplus,velzplus,Bvcxplus,Bvcyplus,Bvczplus,psidcplus,epsplus,dc_flag,&
     trivial_rp, hydro_excision_mask,&
     gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
     dir, ni, nj, nrx, nry, nrz, ev_l, ev_r, xw)

  USE GRHydro_Scalars
  USE GRHydro_EigenproblemM

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: handle,poly,nx
  CCTK_REAL :: dx
  CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps
  CCTK_REAL, dimension(nx) :: Bvcx,Bvcy,Bvcz,psidc
  CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus
  CCTK_REAL, dimension(nx) :: Bvcxminus,Bvcyminus,Bvczminus,psidcminus
  CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus
  CCTK_REAL, dimension(nx) :: Bvcxplus,Bvcyplus,Bvczplus,psidcplus

  CCTK_INT :: i,s,dc_flag
  CCTK_REAL, dimension(nx) :: drho,dvelx,dvely,dvelz,deps
  CCTK_REAL, dimension(nx) :: dBvcx,dBvcy,dBvcz,dpsidc
  CCTK_REAL, dimension(nx) :: dmrho,dmvelx,dmvely,dmvelz,dmeps
  CCTK_REAL, dimension(nx) :: dmBvcx,dmBvcy,dmBvcz,dmpsidc
  CCTK_REAL, dimension(nx) :: press,dpress,d2rho,tilde_flatten
  CCTK_REAL :: dpress2,dvel,w,flatten,eta,etatilde

  logical, dimension(nx) :: trivial_rp

  CCTK_INT, dimension(nx) :: hydro_excision_mask

  CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
                              beta, alp, w_lorentz
  CCTK_INT :: dir, ni, nj, nrx, nry, nrz
  CCTK_REAL, dimension(nrx, nry, nrz) :: ev_l, ev_r, xw

  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det
  CCTK_REAL, dimension(5) :: lam
  CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz
  CCTK_REAL, dimension(nx) :: xwind, l_ev_l, l_ev_r

  logical :: cond

  CCTK_REAL, dimension(nx) :: vx,vy,vz
  
  
  !!$ reconstruct w^i = w_lorentz*v^i to ensure slower than light
  !!$ speeds?
  if (reconstruct_Wv.ne.0) then
    ! all velocity like quantities are now w_lorentz*velocity. We will convert
    ! back to ordinary velocity at the end, using w_lorentz = sqrt(1 + g_{ij}
    ! w^i w^j). 
    vx = w_lorentz*velx
    vy = w_lorentz*vely
    vz = w_lorentz*velz
  else
    vx = velx
    vy = vely
    vz = velz
  end if


!!$  Initially, all the Riemann problems will be trivial

trivial_rp = .true.

!!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
!!$  This is the expression for an even grid.

  do i = 2, nx - 1
    drho(i) = 0.5d0 * (rho(i+1) - rho(i-1))
    dvelx(i) = 0.5d0 * (vx(i+1) - vx(i-1))
    dvely(i) = 0.5d0 * (vy(i+1) - vy(i-1))
    dvelz(i) = 0.5d0 * (vz(i+1) - vz(i-1))
    dBvcx(i) = 0.5d0 * (Bvcx(i+1) - Bvcx(i-1))
    dBvcy(i) = 0.5d0 * (Bvcy(i+1) - Bvcy(i-1))
    dBvcz(i) = 0.5d0 * (Bvcz(i+1) - Bvcz(i-1))
    if(dc_flag.ne.0)dpsidc(i) = 0.5d0 * (psidc(i+1) - psidc(i-1))

    dpress(i) = press(i+1) - press(i-1)
    d2rho(i) = (rho(i+1) - 2.d0 * rho(i) + rho(i-1))! / 6.d0 / dx / dx
    ! since we use d2rho only for the condition d2rho(i+1)*d2rhoi(i-1)<0 
    ! the denominator is not necessary
  end do
  if (poly .eq. 0) then
    do i = 2, nx - 1
      deps(i) = 0.5d0 * (eps(i+1) - eps(i-1))
    end do
  end if

!!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

  do i = 2, nx - 1
#define STEEP(x,dx,dmx)                                              \
    if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
      dmx(i) = sign(1.d0, dx(i)) *                                   \
           min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
                           2.d0 * abs(x(i+1) - x(i)))              &&\
    else                                                           &&\
      dmx(i) = 0.d0                                                &&\
    end if
    STEEP(rho, drho, dmrho)
    STEEP(vx, dvelx, dmvelx)
    STEEP(vy, dvely, dmvely)
    STEEP(vz, dvelz, dmvelz)
    STEEP(Bvcx, dBvcx, dmBvcx)
    STEEP(Bvcy, dBvcy, dmBvcy)
    STEEP(Bvcz, dBvcz, dmBvcz)
    if(dc_flag.ne.0) then
       STEEP(psidc,dpsidc,dmpsidc)
    endif

  end do
  if (poly .eq. 0) then
    do i = 2, nx - 1
      STEEP(eps, deps, dmeps)
    end do
  end if

!!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

  do i = 2, nx-2
    rhoplus(i) = 0.5d0 * (rho(i) + rho(i+1)) + &
         (dmrho(i) - dmrho(i+1)) / 6.d0
    rhominus(i+1) = rhoplus(i)
    velxplus(i) = 0.5d0 * (vx(i) + vx(i+1)) + &
         (dmvelx(i) - dmvelx(i+1)) / 6.d0
    velxminus(i+1) = velxplus(i)
    velyplus(i) = 0.5d0 * (vy(i) + vy(i+1)) + &
         (dmvely(i) - dmvely(i+1)) / 6.d0
    velyminus(i+1) = velyplus(i)
    velzplus(i) = 0.5d0 * (vz(i) + vz(i+1)) + &
         (dmvelz(i) - dmvelz(i+1)) / 6.d0
    velzminus(i+1) = velzplus(i)

    Bvcxplus(i) = 0.5d0 * (Bvcx(i) + Bvcx(i+1)) + &
         (dmBvcx(i) - dmBvcx(i+1)) / 6.d0
    Bvcxminus(i+1) = Bvcxplus(i)
    Bvcyplus(i) = 0.5d0 * (Bvcy(i) + Bvcy(i+1)) + &
         (dmBvcy(i) - dmBvcy(i+1)) / 6.d0
    Bvcyminus(i+1) = Bvcyplus(i)
    Bvczplus(i) = 0.5d0 * (Bvcz(i) + Bvcz(i+1)) + &
         (dmBvcz(i) - dmBvcz(i+1)) / 6.d0
    Bvczminus(i+1) = Bvczplus(i)
    if(dc_flag.ne.0) then
       psidcplus(i) = 0.5d0 * (psidc(i) + psidc(i+1)) + &
            (dmpsidc(i) - dmpsidc(i+1)) / 6.d0
       psidcminus(i+1) = psidcplus(i)
    endif
  end do
  if (poly .eq. 0) then
    do i = 2, nx-2
      epsplus(i) = 0.5d0 * (eps(i) + eps(i+1)) + &
           (dmeps(i) - dmeps(i+1)) / 6.d0
      epsminus(i+1) = epsplus(i)
    end do
  end if

!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which mat be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even. 
!!$Note also that we don''t have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).

  if (ppm_detect .ne. 0) then

    do i = 3, nx - 2
      if ( (d2rho(i+1)*d2rho(i-1) < 0.d0).and.(abs(rho(i+1)-rho(i-1)) - &
           ppm_epsilon_shock * min(abs(rho(i+1)), abs(rho(i-1))) > 0.d0) ) then
        etatilde = (rho(i-2) - rho(i+2) + 4.d0 * drho(i)) / (drho(i) * 12.d0)
      else
        etatilde = 0.d0
      end if
      eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
      if (ppm_k0 * abs(drho(i)) * min(press(i-1),press(i+1)) < &
           abs(dpress(i)) * min(rho(i-1), rho(i+1))) then
        eta = 0.d0
      end if
      if (eta > 0.d0) then
        trivial_rp(i-1) = .false.
        trivial_rp(i) = .false.
      end if
      rhominus(i) = rhominus(i) * (1.d0 - eta) + &
           (rho(i-1) + 0.5d0 * dmrho(i-1)) * eta
      rhoplus(i) = rhoplus(i) * (1.d0 - eta) + &
           (rho(i+1) - 0.5d0 * dmrho(i+1)) * eta
    end do

  end if

  !!$ mppm
#define D_UPW(x) (0.5d0 * (x(i) + x(i+1)))
#define LEFT1(x)  (13.d0*x(i+1)-5.d0*x(i+2)+x(i+3)+3.d0*x(i  ))/12.d0
#define RIGHT1(x) (13.d0*x(i  )-5.d0*x(i-1)+x(i-2)+3.d0*x(i+1))/12.d0
  if (ppm_mppm .gt. 0) then
    l_ev_l=0.d0
    l_ev_r=0.d0
    xwind=0.d0
    do i=3, nx - 3
      agxx = 0.5d0*( gxx(i) + gxx(i+1) )
      agxy = 0.5d0*( gxy(i) + gxy(i+1) )
      agxz = 0.5d0*( gxz(i) + gxz(i+1) )
      agyy = 0.5d0*( gyy(i) + gyy(i+1) )
      agyz = 0.5d0*( gyz(i) + gyz(i+1) )
      agzz = 0.5d0*( gzz(i) + gzz(i+1) )
      det = SPATIAL_DETERMINANT(agxx, agxy, agxz, \
                              agyy, agyz, agzz)
      call UpperMetric (uxx, uxy, uxz, uyy, uyz, uzz, &
                        det, agxx, agxy, agxz, agyy, agyz, agzz)
      call eigenvaluesM(handle,&
                 D_UPW(rho), D_UPW(velx), D_UPW(vely), D_UPW(velz), &
                 D_UPW(eps), D_UPW(press), D_UPW(w_lorentz),&
                 D_UPW(Bvcx), D_UPW(Bvcy), D_UPW(Bvcz), lam, &
                 agxx, agxy, agxz, agyy, agyz, agzz, &
                 uxx, D_UPW(alp), D_UPW(beta))
      l_ev_l(i)=lam(1)
      l_ev_r(i)=lam(5)
      xwind(i) = (lam(1) + lam(5)) / (abs(lam(1)) + abs(lam(5)))
      xwind(i) = min(1.d0, max(-1.d0, xwind(i)))
#define LEFTPLUS(x,xplus)    xplus(i)   =       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define LEFTMINUS(x,xminus)  xminus(i+1)=       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define RIGHTPLUS(x,xplus)   xplus(i)   =       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define RIGHTMINUS(x,xminus) xminus(i+1)=       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define CHECK(x,xc) if (x(i+1) .gt. x(i)) then && xc=min(x(i+1),max(x(i),xc)) && else && xc=min(x(i),max(x(i+1),xc)) && endif
!!$      xwind(i)=0.d0
      if (xwind(i) .lt. 0.0d0) then
        LEFTPLUS(rho, rhoplus)
        LEFTMINUS(rho, rhominus)
        LEFTPLUS(velx, velxplus)
        LEFTMINUS(velx, velxminus)
        LEFTPLUS(vely, velyplus)
        LEFTMINUS(vely, velyminus)
        LEFTPLUS(velz, velzplus)
        LEFTMINUS(velz, velzminus)
        LEFTPLUS(Bvcx, Bvcxplus)
        LEFTMINUS(Bvcx, Bvcxminus)
        LEFTPLUS(Bvcy, Bvcyplus)
        LEFTMINUS(Bvcy, Bvcyminus)
        LEFTPLUS(Bvcz, Bvczplus)
        LEFTMINUS(Bvcz, Bvczminus)
        if(dc_flag.ne.0) then
           LEFTPLUS(psidc, psidcplus)
           LEFTMINUS(psidc, psidcminus)
        endif
        if (poly .eq. 0) then
          LEFTPLUS(eps, epsplus)
          LEFTMINUS(eps, epsminus)
        end if
      else
        RIGHTPLUS(rho, rhoplus)
        RIGHTMINUS(rho, rhominus)
        RIGHTPLUS(velx, velxplus)
        RIGHTMINUS(velx, velxminus)
        RIGHTPLUS(vely, velyplus)
        RIGHTMINUS(vely, velyminus)
        RIGHTPLUS(velz, velzplus)
        RIGHTMINUS(velz, velzminus)
        RIGHTPLUS(Bvcx, Bvcxplus)
        RIGHTMINUS(Bvcx, Bvcxminus)
        RIGHTPLUS(Bvcy, Bvcyplus)
        RIGHTMINUS(Bvcy, Bvcyminus)
        RIGHTPLUS(Bvcz, Bvczplus)
        RIGHTMINUS(Bvcz, Bvczminus)
        if(dc_flag.ne.0) then
           RIGHTPLUS(psidc, psidcplus)
           RIGHTMINUS(psidc, psidcminus)
        endif
        if (poly .eq. 0) then
          RIGHTPLUS(eps, epsplus)
          RIGHTMINUS(eps, epsminus)
        end if
      end if
      CHECK(rho, rhoplus(i))
      CHECK(rho, rhominus(i+1))
      CHECK(velx, velxplus(i))
      CHECK(velx, velxminus(i+1))
      CHECK(vely, velyplus(i))
      CHECK(vely, velyminus(i+1))
      CHECK(velz, velzplus(i))
      CHECK(velz, velzminus(i+1))
      CHECK(Bvcx, Bvcxplus(i))
      CHECK(Bvcx, Bvcxminus(i+1))
      CHECK(Bvcy, Bvcyplus(i))
      CHECK(Bvcy, Bvcyminus(i+1))
      CHECK(Bvcz, Bvczplus(i))
      CHECK(Bvcz, Bvczminus(i+1))
      if(dc_flag.ne.0) then
         CHECK(psidc, psidcplus(i))
         CHECK(psidc, psidcminus(i+1))
      endif
      if (poly .eq. 0) then
        CHECK(eps, epsplus(i))
        CHECK(eps, epsminus(i+1))
      end if
!!$      if ((dir .eq. 1) .and. (ni .eq. 4) .and. (nj .eq. 4)) then
!!$        write (*,*) rhoplus(i), rhominus(i+1)
!!$      end if
    end do
    !!$ mppm debug output
    if (ppm_mppm_debug_eigenvalues .gt. 0) then
      if (dir .eq. 1) then
        ev_l(:,ni,nj) = l_ev_l
        ev_r(:,ni,nj) = l_ev_r
        xw(:,ni,nj) = xwind
      else if (dir .eq. 2) then
        ev_l(ni,:,nj) = l_ev_l
        ev_r(ni,:,nj) = l_ev_r
        xw(ni,:,nj) = xwind
      else if (dir .eq. 3) then
        ev_l(ni,nj,:) = l_ev_l
        ev_r(ni,nj,:) = l_ev_r
        xw(ni,nj,:) = xwind
      else
        write (*,*) "flux direction not 1 to 3 ?"
      end if
    end if
  end if

!!$  Zone flattening. See appendix of C&W, p. 197-8.

  do i = 3, nx - 2
    dpress2 = press(i+2) - press(i-2)
    dvel = vx(i+1) - vx(i-1)
    if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
         (dvel < 0.d0) ) then
      w = 1.d0
    else
      w = 0.d0
    end if
    if (abs(dpress2) < ppm_small) then
      tilde_flatten(i) = 1.d0
    else
      tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
           (dpress(i) / dpress2 - ppm_omega1)))
    end if
  end do



  if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
    do i = 3, nx - 2
      flatten = tilde_flatten(i)
      if (abs(1.d0 - flatten) > 0.d0) then
        trivial_rp(i-1) = .false.
        trivial_rp(i) = .false.
      end if
      rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
      rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
      velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
      velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
      velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
      velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
      velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
      velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
      Bvcxplus(i) = flatten * Bvcxplus(i) + (1.d0 - flatten) * Bvcx(i)
      Bvcxminus(i) = flatten * Bvcxminus(i) + (1.d0 - flatten) * Bvcx(i)
      Bvcyplus(i) = flatten * Bvcyplus(i) + (1.d0 - flatten) * Bvcy(i)
      Bvcyminus(i) = flatten * Bvcyminus(i) + (1.d0 - flatten) * Bvcy(i)
      Bvczplus(i) = flatten * Bvczplus(i) + (1.d0 - flatten) * Bvcz(i)
      Bvczminus(i) = flatten * Bvczminus(i) + (1.d0 - flatten) * Bvcz(i)
      if(dc_flag.ne.0) then
         psidcplus(i) = flatten * psidcplus(i) + (1.d0 - flatten) * psidc(i)
         psidcminus(i) = flatten * psidcminus(i) + (1.d0 - flatten) * psidc(i)
      endif
      if (poly .eq. 0) then
        epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
        epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
      end if
    end do
  else  !!$ Really implement C&W, page 197; which requires stencil 4.
    do i = 4, nx - 3
      s=int(sign(1.d0, -dpress(i)))
      flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
      if (abs(1.d0 - flatten) > 0.d0) then
        trivial_rp(i-1) = .false.
        trivial_rp(i) = .false.
      end if
      rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
      rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
      velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
      velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
      velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
      velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
      velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
      velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
      Bvcxplus(i) = flatten * Bvcxplus(i) + (1.d0 - flatten) * Bvcx(i)
      Bvcxminus(i) = flatten * Bvcxminus(i) + (1.d0 - flatten) * Bvcx(i)
      Bvcyplus(i) = flatten * Bvcyplus(i) + (1.d0 - flatten) * Bvcy(i)
      Bvcyminus(i) = flatten * Bvcyminus(i) + (1.d0 - flatten) * Bvcy(i)
      Bvczplus(i) = flatten * Bvczplus(i) + (1.d0 - flatten) * Bvcz(i)
      Bvczminus(i) = flatten * Bvczminus(i) + (1.d0 - flatten) * Bvcz(i)
      if(dc_flag.ne.0) then
         psidcplus(i) = flatten * psidcplus(i) + (1.d0 - flatten) * psidc(i)
         psidcminus(i) = flatten * psidcminus(i) + (1.d0 - flatten) * psidc(i)
      endif
      if (poly .eq. 0) then
        epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
        epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
      end if
    end do
  end if


!!$ Monotonicity. See (1.10) of C&W.

do i = GRHydro_stencil, nx - GRHydro_stencil + 1
#define MON(xminus,x,xplus)                                       \
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) )     \
        .and. ((xplus(i)-x(i))*(x(i)-xminus(i)) .le. 0.d0)) then&&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
      xminus(i) = x(i)                                          &&\
      xplus(i) = x(i)                                           &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) >                         \
                (xplus(i) - xminus(i))**2) then                 &&\
      xminus(i) = 3.d0 * x(i) - 2.d0 * xplus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) <                         \
               -(xplus(i) - xminus(i))**2) then                 &&\
      xplus(i) = 3.d0 * x(i) - 2.d0 * xminus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if                                                      &&\
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) ) ) then     &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if

    MON(rhominus,rho,rhoplus)
    MON(velxminus,vx,velxplus)
    MON(velyminus,vy,velyplus)
    MON(velzminus,vz,velzplus)
    MON(Bvcxminus,Bvcx,Bvcxplus)
    MON(Bvcyminus,Bvcy,Bvcyplus)
    MON(Bvczminus,Bvcz,Bvczplus)
    if(dc_flag.ne.0) then
       MON(psidcminus,psidc,psidcplus)
    endif
  end do
  if (poly .eq. 0) then
    do i = GRHydro_stencil, nx - GRHydro_stencil + 1
      MON(epsminus,eps,epsplus)
    end do
  end if

  if (check_for_trivial_rp .eq. 0) then
    trivial_rp = .false.
  end if

  !!$ excision
  do i = 1, nx
    if (GRHydro_enable_internal_excision /= 0 .and. &
        (hydro_excision_mask(i) .ne. 0)) then
      if (i .gt. 1) then
        trivial_rp(i-1)=.true.
      end if
      trivial_rp(i)=.true.
    else
      !!$ Do not optimize cond away by combining the 'if's. Fortran does not
      !!$  have to follow the order of sub-expressions given here and might
      !!$  access outside the array range
      cond = .false.
      if (i .gt. 1 .and. GRHydro_enable_internal_excision /= 0) then
        cond = hydro_excision_mask(i-1) .ne. 0
      end if
      if (cond) then
        rhominus(i)=rho(i)
        rhoplus(i)=rho(i)
        velxminus(i)=vx(i)
        velxplus(i)=vx(i)
        velyminus(i)=vy(i)
        velyplus(i)=vy(i)
        velzminus(i)=vz(i)
        velzplus(i)=vz(i)
        Bvcxminus(i)=Bvcx(i)  
        Bvcxplus(i)=Bvcx(i)
        Bvcyminus(i)=Bvcy(i)
        Bvcyplus(i)=Bvcy(i)
        Bvczminus(i)=Bvcz(i)
        Bvczplus(i)=Bvcz(i)
        rhominus(i-1)=rho(i)
        rhoplus(i-1)=rho(i)
        velxminus(i-1)=vx(i)
        velxplus(i-1)=vx(i)
        velyminus(i-1)=vy(i)
        velyplus(i-1)=vy(i)
        velzminus(i-1)=vz(i)
        velzplus(i-1)=vz(i)
        Bvcxminus(i-1)=Bvcx(i)
        Bvcxplus(i-1)=Bvcx(i)
        Bvcyminus(i-1)=Bvcy(i)
        Bvcyplus(i-1)=Bvcy(i)
        Bvczminus(i-1)=Bvcz(i)
        Bvczplus(i-1)=Bvcz(i)
        if (dc_flag.ne.0) then
          psidcminus(i)=psidc(i)
          psidcplus(i)=psidc(i)
          psidcminus(i-1)=psidc(i)
          psidcplus(i-1)=psidc(i)
        end if
       if (poly .eq. 0) then
          epsminus(i)=eps(i)
          epsplus(i)=eps(i)
          epsminus(i-1)=eps(i)
          epsplus(i-1)=eps(i)
        end if
      else
        cond = .false.
        if ((i.gt.2) .and. (i.lt.nx) .and. GRHydro_enable_internal_excision /= 0) then
          cond = (ppm_mppm .eq. 0) .and. (hydro_excision_mask(i-2) .ne. 0)
        end if
        if (cond) then
          call PPM_TVD(rho(i-1), rho(i), rho(i+1), rhominus(i), rhoplus(i))
          call PPM_TVD(vx(i-1), vx(i), vx(i+1), velxminus(i), velxplus(i))
          call PPM_TVD(vy(i-1), vy(i), vy(i+1), velyminus(i), velyplus(i))
          call PPM_TVD(vz(i-1), vz(i), vz(i+1), velzminus(i), velzplus(i))
          call PPM_TVD(Bvcx(i-1), Bvcx(i), Bvcx(i+1), Bvcxminus(i), Bvcxplus(i))
          call PPM_TVD(Bvcy(i-1), Bvcy(i), Bvcy(i+1), Bvcyminus(i), Bvcyplus(i))
          call PPM_TVD(Bvcz(i-1), Bvcz(i), Bvcz(i+1), Bvczminus(i), Bvczplus(i))
          if(dc_flag.ne.0) then
             call PPM_TVD(psidc(i-1), psidc(i), psidc(i+1), psidcminus(i), psidcplus(i))
          endif
          if (poly .eq. 0) then
            call PPM_TVD(eps(i-1), eps(i), eps(i+1), epsminus(i), epsplus(i))
          end if
        end if
      end if
      cond = .false.
      if (i .lt. nx .and. GRHydro_enable_internal_excision /= 0) then
        cond = hydro_excision_mask(i+1) .ne. 0
      end if
      if (cond) then
        rhominus(i)=rho(i)
        rhoplus(i)=rho(i)
        velxminus(i)=vx(i)
        velxplus(i)=vx(i)
        velyminus(i)=vy(i)
        velyplus(i)=vy(i)
        velzminus(i)=vz(i)
        velzplus(i)=vz(i)
        Bvcxminus(i)=Bvcx(i)
        Bvcxplus(i)=Bvcx(i)
        Bvcyminus(i)=Bvcy(i)
        Bvcyplus(i)=Bvcy(i)
        Bvczminus(i)=Bvcz(i)
        Bvczplus(i)=Bvcz(i)
        rhominus(i+1)=rho(i)
        rhoplus(i+1)=rho(i)
        velxminus(i+1)=vx(i)
        velxplus(i+1)=vx(i)
        velyminus(i+1)=vy(i)
        velyplus(i+1)=vy(i)
        velzminus(i+1)=vz(i)
        velzplus(i+1)=vz(i)
        Bvcxminus(i+1)=Bvcx(i)
        Bvcxplus(i+1)=Bvcx(i)
        Bvcyminus(i+1)=Bvcy(i)
        Bvcyplus(i+1)=Bvcy(i)
        Bvczminus(i+1)=Bvcz(i)
        Bvczplus(i+1)=Bvcz(i)
        if (dc_flag.ne.0) then
          psidcminus(i)=psidc(i)
          psidcplus(i)=psidc(i)
          psidcminus(i+1)=psidc(i)
          psidcplus(i+1)=psidc(i)
        endif
        if (poly .eq. 0) then
          epsminus(i)=eps(i)
          epsplus(i)=eps(i)
          epsminus(i+1)=eps(i)
          epsplus(i+1)=eps(i)
        endif
      else
        cond = .false.
        if ((i.lt.nx-1) .and. (i.gt.1) .and. GRHydro_enable_internal_excision /= 0) then
          cond = (ppm_mppm .eq. 0) .and. (hydro_excision_mask(i+2) .ne. 0)
        end if
        if (cond) then
          call PPM_TVD(rho(i-1), rho(i), rho(i+1), rhominus(i), rhoplus(i))
          call PPM_TVD(vx(i-1), vx(i), vx(i+1), velxminus(i), velxplus(i))
          call PPM_TVD(vy(i-1), vy(i), vy(i+1), velyminus(i), velyplus(i))
          call PPM_TVD(vz(i-1), vz(i), vz(i+1), velzminus(i), velzplus(i))
          call PPM_TVD(Bvcx(i-1), Bvcx(i), Bvcx(i+1), Bvcxminus(i), Bvcxplus(i))
          call PPM_TVD(Bvcy(i-1), Bvcy(i), Bvcy(i+1), Bvcyminus(i), Bvcyplus(i))
          call PPM_TVD(Bvcz(i-1), Bvcz(i), Bvcz(i+1), Bvczminus(i), Bvczplus(i))
          if(dc_flag.ne.0) then
             call PPM_TVD(Bvcz(i-1), Bvcz(i), Bvcz(i+1), Bvczminus(i), Bvczplus(i))
          endif
          if (poly .eq. 0) then
            call PPM_TVD(eps(i-1), eps(i), eps(i+1), epsminus(i), epsplus(i))
          end if
        end if
      end if
    end if
  end do

  !!$ transform back to vel if Wv was used!
  if (reconstruct_Wv.ne.0) then
    do i = grhydro_stencil, nx - grhydro_stencil + 1
      ! divide out the Loretnz factor obtained from w_lorentz =
      ! sqrt(1+g_{ij} w^i w^j) for both the
      ! plus and minus quantities this should by construction ensure
      ! that any Lorentz factor calculated from them later on is
      ! physical (ie. > 1.d0)
      agxx = 0.5d0*( gxx(i) + gxx(i-1) )
      agxy = 0.5d0*( gxy(i) + gxy(i-1) )
      agxz = 0.5d0*( gxz(i) + gxz(i-1) )
      agyy = 0.5d0*( gyy(i) + gyy(i-1) )
      agyz = 0.5d0*( gyz(i) + gyz(i-1) )
      agzz = 0.5d0*( gzz(i) + gzz(i-1) )
      w = sqrt( 1.d0 + agxx*velxminus(i)*velxminus(i) + &
                       agyy*velyminus(i)*velyminus(i) &
              + agzz*velzminus(i)*velzminus(i) + &
                2.d0*agxy*velxminus(i)*velyminus(i) &
              + 2.d0*agxz*velxminus(i)*velzminus(i) + &
                2.d0*agyz*velyminus(i)*velzminus(i) )
      velxminus(i) = velxminus(i)/w
      velyminus(i) = velyminus(i)/w
      velzminus(i) = velzminus(i)/w
      
      agxx = 0.5d0*( gxx(i) + gxx(i+1) )
      agxy = 0.5d0*( gxy(i) + gxy(i+1) )
      agxz = 0.5d0*( gxz(i) + gxz(i+1) )
      agyy = 0.5d0*( gyy(i) + gyy(i+1) )
      agyz = 0.5d0*( gyz(i) + gyz(i+1) )
      agzz = 0.5d0*( gzz(i) + gzz(i+1) )
      w = sqrt( 1.d0 + agxx*velxplus(i)*velxplus(i) + &
                agyy*velyplus(i)*velyplus(i) &
              + agzz*velzplus(i)*velzplus(i) + &
               2.d0*agxy*velxplus(i)*velyplus(i) &
              + 2.d0*agxz*velxplus(i)*velzplus(i) + &
               2.d0*agyz*velyplus(i)*velzplus(i) )
      velxplus(i) = velxplus(i)/w
      velyplus(i) = velyplus(i)/w
      velzplus(i) = velzplus(i)/w
    end do
  end if
  
  
return

end subroutine SimplePPM_1dM


