 /*@@
   @file      GRHydro_PPM.F90
   @date      Sun Feb 10 16:53:29 2002
   @author    Ian Hawke, Toni Font, Luca Baiotti, Frank Loeffler
   @desc 
   Routines to do PPM reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

subroutine PPM_TVD(origm, orig, origp, bextm, bextp)
  CCTK_REAL :: origm, orig, origp, bextm, bextp
  CCTK_REAL :: dloc, dupw, delta

  dupw = orig - origm
  dloc = origp - orig
  if (dupw*dloc < 0.d0) then
    delta=0.d0
  else if (abs(dupw) < abs(dloc)) then
    delta=dupw
  else
    delta=dloc
  end if
  bextm = orig - 0.5d0 * delta
  bextp = orig + 0.5d0 * delta
end subroutine PPM_TVD

 /*@@
   @routine    SimplePPM_1d
   @date       Thu Feb 14 19:08:52 2002
   @author     Ian Hawke, Toni Font, Christian Reisswig
   @desc 
   The simple PPM reconstruction routine that applies along
   each one dimensional slice.

   @enddesc 
   @calls     
   @calledby   
   @history 
   Written in frustration when IH couldn''t get Toni''s original code 
   to work.
   Later extended to enhanced PPM scheme by CR.
   @endhistory 

@@*/

#define SpaceMask_CheckStateBitsF90_1D(mask,i,type_bits,state_bits) \
  (iand(mask((i)),(type_bits)).eq.(state_bits))


subroutine SimplePPM_1d(apply_enhanced_ppm,handle,poly,&
     nx,dx,rho,velx,vely,velz,eps,press,rhominus,&
     velxminus,velyminus,velzminus,epsminus,rhoplus,velxplus,velyplus,&
     velzplus,epsplus,trivial_rp, hydro_excision_mask,&
     gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
     dir, ni, nj, nrx, nry, nrz, ev_l, ev_r, xw)

  USE GRHydro_Scalars
  USE GRHydro_Eigenproblem

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, parameter :: one = 1

  logical :: apply_enhanced_ppm

  CCTK_INT :: handle,poly,nx
  CCTK_REAL :: dx
  CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps
  CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus
  CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus

  CCTK_INT :: i,s
  CCTK_REAL, dimension(nx) :: drho,dvelx,dvely,dvelz,deps
  CCTK_REAL, dimension(nx) :: dmrho,dmvelx,dmvely,dmvelz,dmeps
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

  CCTK_REAL :: D2a, D2aL, D2aR, D2aC, D2aLim, rhi, daplus, daminus, D3a, D3aLL, D3aL, D3aR, D2aLL, D2aRR, D3aMin, D3aMax

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


#define STEEP(x,dx,dmx)                                              \
         if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
            dmx(i) = sign(one, dx(i)) *                                   \
               min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
                                 2.d0 * abs(x(i+1) - x(i)))              &&\
         else                                                           &&\
            dmx(i) = 0.d0                                                &&\
         end if


  if (.not.apply_enhanced_ppm) then
      !! This is the original PPM algorithm by Colella & Woodward 1984.

      !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
      !!$  This is the expression for an even grid.

      do i = 2, nx - 1
         drho(i) = 0.5d0 * (rho(i+1) - rho(i-1))
         dvelx(i) = 0.5d0 * (vx(i+1) - vx(i-1))
         dvely(i) = 0.5d0 * (vy(i+1) - vy(i-1))
         dvelz(i) = 0.5d0 * (vz(i+1) - vz(i-1))
         dpress(i) = press(i+1) - press(i-1)
         d2rho(i) = ((rho(i+1) + rho(i-1)) - 2.d0 * rho(i))! / 6.d0 / dx / dx
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
         STEEP(rho, drho, dmrho)
         STEEP(vx, dvelx, dmvelx)
         STEEP(vy, dvely, dmvely)
         STEEP(vz, dvelz, dmvelz)
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
      end do
      if (poly .eq. 0) then
         do i = 2, nx-2
            epsplus(i) = 0.5d0 * (eps(i) + eps(i+1)) + &
               (dmeps(i) - dmeps(i+1)) / 6.d0
            epsminus(i+1) = epsplus(i)
         end do
      end if
  else
!! This is the modified PPM algorithm by Colella & Sekora 2008 and McCorquodale & Colella 2011.
!! This uses a better limiter based on second derivatives that preserves
!! accuracy at local extrema. It also uses a higher-order interpolation polynomial.

      !!$  Initial boundary states (sixth order accurate). See (17) of Colella and Sekora 2008, p.7071
#define APPROX_AT_CELL_INTERFACE_STENCIL4(a, ah)  \
      ah = 37.0d0/60.0d0*(a(i)+a(i+1)) - 2.0d0/15.0d0*(a(i-1)+a(i+2)) + 1.0d0/60.0d0*(a(i-2)+a(i+3))
      
      !!$  Initial boundary states (4th order accurate). See (16) of Colella and Sekora 2008, p.7071
#define APPROX_AT_CELL_INTERFACE(a, ah)  \
      ah = 7.0d0/12.0d0*(a(i)+a(i+1)) - 1.0d0/12.0d0*(a(i-1)+a(i+2))
      
#define LIMIT(a,ah,C, alim) \
      if ((min(a(i),a(i+1)) .le. ah) .and. (ah .le. max(a(i),a(i+1)))) then       &&\
         alim = ah   &&\
      else &&\
         D2a  = 3.0d0 * ((a(i)   + a(i+1)) - 2.0d0*ah    )                                 &&\
         D2aL =         ((a(i-1) + a(i+1)) - 2.0d0*a(i)  )                                 &&\
         D2aR =         ((a(i)   + a(i+2)) - 2.0d0*a(i+1))                                 &&\
         D2aLim = sign(one, D2a)*min(C*abs(D2aL), C*abs(D2aR), abs(D2a))/3.0d0             &&\
         if (D2a*D2aR .ge. 0 .and. D2a*D2aL .ge. 0) then                     &&\
            alim = 0.5d0*(a(i)+a(i+1)) - D2aLim &&\
         else                                                                             &&\
            alim = 0.5d0*(a(i)+a(i+1))                                                    &&\
         end if                                                                           &&\
      endif

#define LIMIT_EPS(a,ah,C, alim) \
      if ((min(a(i),a(i+1)) .le. ah) .and. (ah .le. max(a(i),a(i+1)))) then       &&\
         alim = ah   &&\
      else &&\
         D2a  = 3.0d0 * ((a(i)   + a(i+1)) - 2.0d0*ah    )                                 &&\
         D2aL =         ((a(i-1) + a(i+1)) - 2.0d0*a(i)  )                                 &&\
         D2aR =         ((a(i)   + a(i+2)) - 2.0d0*a(i+1))                                 &&\
         D2aLim = sign(one, D2a)*min(C*abs(D2aL), C*abs(D2aR), abs(D2a))/3.0d0             &&\
         if (D2a*D2aR .ge. 0 .and. D2a*D2aL .ge. 0 .and. abs(D2aLim) .lt. abs(0.5d0*(a(i)+a(i+1)))) then                     &&\
            alim = 0.5d0*(a(i)+a(i+1)) - D2aLim &&\
         else                                                                             &&\
            alim = 0.5d0*(a(i)+a(i+1))                                                    &&\
         end if                                                                           &&\
      endif



   if (PPM3) then
      
      !! We initialize "plus" \equiv a_j+1/2 with (16) via APPROX_AT_CELL_INTERFACE, 
      !! then checking for (13) of Colella & Sekora 2008 and applying
      !! (18) and (19) if (13) is not satisfied. This is done with LIMIT.
      do i = 2, nx-2
         APPROX_AT_CELL_INTERFACE(rho, rhoplus(i))
         LIMIT(rho, rhoplus(i), enhanced_ppm_C2, rhoplus(i))
         rhominus(i+1) = rhoplus(i)
         
         APPROX_AT_CELL_INTERFACE(vx, velxplus(i))
         LIMIT(vx, velxplus(i), enhanced_ppm_C2, velxplus(i))
         velxminus(i+1) = velxplus(i)
         
         APPROX_AT_CELL_INTERFACE(vy, velyplus(i))
         LIMIT(vy, velyplus(i), enhanced_ppm_C2, velyplus(i))
         velyminus(i+1) = velyplus(i)
         
         APPROX_AT_CELL_INTERFACE(vz, velzplus(i))
         LIMIT(vz, velzplus(i), enhanced_ppm_C2, velzplus(i))
         velzminus(i+1) = velzplus(i)
      end do
      if (poly .eq. 0) then
         do i = 2, nx-2
            APPROX_AT_CELL_INTERFACE(eps, epsplus(i))
            LIMIT_EPS(eps, epsplus(i), enhanced_ppm_C2, epsplus(i))
            epsminus(i+1) = epsplus(i)
         end do
      end if
      
    else

      !! Same as above but for 4 stencil points using (17) of Colella & Sekora 2008 as
      !! initial states.
      do i = 3, nx-3
         APPROX_AT_CELL_INTERFACE_STENCIL4(rho, rhoplus(i))
         LIMIT(rho, rhoplus(i), enhanced_ppm_C2, rhoplus(i))
         rhominus(i+1) = rhoplus(i)
         
         APPROX_AT_CELL_INTERFACE_STENCIL4(vx, velxplus(i))
         LIMIT(vx, velxplus(i), enhanced_ppm_C2, velxplus(i))
         velxminus(i+1) = velxplus(i)
         
         APPROX_AT_CELL_INTERFACE_STENCIL4(vy, velyplus(i))
         LIMIT(vy, velyplus(i), enhanced_ppm_C2, velyplus(i))
         velyminus(i+1) = velyplus(i)
         
         APPROX_AT_CELL_INTERFACE_STENCIL4(vz, velzplus(i))
         LIMIT(vz, velzplus(i), enhanced_ppm_C2, velzplus(i))
         velzminus(i+1) = velzplus(i)
      end do
      if (poly .eq. 0) then
         do i = 3, nx-3
            APPROX_AT_CELL_INTERFACE_STENCIL4(eps, epsplus(i))
            LIMIT(eps, epsplus(i), enhanced_ppm_C2, epsplus(i))
            epsminus(i+1) = epsplus(i)
         end do
      end if
    endif
    
    !! Finally compute pressure gradient needed for flattening and shock detection
    do i = 2, nx-1
       dpress(i) = press(i+1) - press(i-1)
    end do
    
  endif


!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which mat be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even. 
!!$Note also that we don''t have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).

  if (.not.apply_enhanced_ppm) then 
   ! Only for 1984 PPM scheme!
   if (ppm_detect .ne. 0) then

      do i = 3, nx - 2
         if ( (d2rho(i+1)*d2rho(i-1) < 0.d0).and.(abs(rho(i+1)-rho(i-1)) - &
            ppm_epsilon_shock * min(abs(rho(i+1)), abs(rho(i-1))) > 0.d0) ) then
         etatilde = ((rho(i-2) - rho(i+2)) + 4.d0 * drho(i)) / (drho(i) * 12.d0)
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
         call eigenvalues(handle,&
                  D_UPW(rho), D_UPW(velx), D_UPW(vely), D_UPW(velz), &
                  D_UPW(eps), D_UPW(w_lorentz), lam, &
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


   if (.not.apply_enhanced_ppm) then
      ! In 1984 PPM, flattening is applied before constraining parabolic profiles.
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
            if (poly .eq. 0) then
               epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
               epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
            end if
         end do
      else  !!$ Really implement C&W, page 197; which requires stencil 4.
         do i = 4, nx - 3
            s=int(sign(one, -dpress(i)))
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
            if (poly .eq. 0) then
               epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
               epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
            end if
         end do
      end if
   endif

!!$ Monotonicity. See (1.10) of C&W.
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


!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34
!! This requires 4 stencil points.
#define MON_WITH_LOCAL_EXTREMUM_STENCIL4(aminus, a, aplus, C, C3) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC =  (a(i-1) + a(i+1)) - 2.0d0*a(i)                               &&\
         D2aL =  (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                             &&\
         D2aLL = (a(i-3) + a(i-1)) - 2.0d0*a(i-2)                             &&\
         D2aR =  (a(i)   + a(i+2)) - 2.0d0*a(i+1)                             &&\
         D2aRR = (a(i+1) + a(i+3)) - 2.0d0*a(i+2)                             &&\
         D3a = D2aR - D2aC   &&\
         D3aL = D2aC - D2aL  &&\
         D3aR = D2aRR - D2aR &&\
         D3aLL = D2aL - D2aLL &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            D3aMin = min(D3aLL, D3aL, D3a, D3aR)   &&\
            D3aMax = max(D3aLL, D3aL, D3a, D3aR)   &&\
            if (C3 * max(abs(D3aMin), abs(D3aMax)) .le. D3aMax-D3aMin) then &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) + daplus * rhi     &&\
                     aminus(i) = a(i) - daminus * rhi    &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus) &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
                  endif &&\
            endif &&\
         endif &&\
         trivial_rp(i-1) = .false.                                 &&\
         trivial_rp(i) = .false.                                   &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
         trivial_rp(i-1) = .false.  &&\
         trivial_rp(i) = .false.    &&\
      endif

!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This contains an additional limiter in case the correction becomes larger than
!! the corrected value. This is to avoid negative values for epsilon.
!! This requires 4 stencil points.
#define MON_WITH_LOCAL_EXTREMUM_STENCIL4_EPS(aminus, a, aplus, C, C3) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC =  (a(i-1) + a(i+1)) - 2.0d0*a(i)                               &&\
         D2aL =  (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                             &&\
         D2aLL = (a(i-3) + a(i-1)) - 2.0d0*a(i-2)                             &&\
         D2aR =  (a(i)   + a(i+2)) - 2.0d0*a(i+1)                             &&\
         D2aRR = (a(i+1) + a(i+3)) - 2.0d0*a(i+2)                             &&\
         D3a = D2aR - D2aC   &&\
         D3aL = D2aC - D2aL  &&\
         D3aR = D2aRR - D2aR &&\
         D3aLL = D2aL - D2aLL &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            D3aMin = min(D3aLL, D3aL, D3a, D3aR)   &&\
            D3aMax = max(D3aLL, D3aL, D3a, D3aR)   &&\
            if (C3 * max(abs(D3aMin), abs(D3aMax)) .le. D3aMax-D3aMin) then &&\
               if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) + daplus * rhi     &&\
                     aminus(i) = a(i) - daminus * rhi    &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)  &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
                  endif &&\
               else &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) &&\
                     aminus(i) = a(i) &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) &&\
                  endif &&\
               endif &&\
            endif &&\
         endif &&\
         trivial_rp(i-1) = .false.                                 &&\
         trivial_rp(i) = .false.                                   &&\
      else                                                         &&\
         if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
         else &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i)               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i)               &&\
            end if                                                    &&\
         endif &&\
         trivial_rp(i-1) = .false.  &&\
         trivial_rp(i) = .false.    &&\
      endif



!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
#define MON_WITH_LOCAL_EXTREMUM(aminus, a, aplus, C) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC = (a(i-1) + a(i+1)) - 2.0d0*a(i)                                &&\
         D2aL = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                              &&\
         D2aR = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                              &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) + daplus * rhi     &&\
                  aminus(i) = a(i) - daminus * rhi    &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
               endif &&\
         endif &&\
         trivial_rp(i-1) = .false.                                 &&\
         trivial_rp(i) = .false.                                   &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
         trivial_rp(i-1) = .false.  &&\
         trivial_rp(i) = .false.    &&\
      endif   
      

!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
!! This contains an additional limiter in case the correction becomes larger than
!! the corrected value. This is to avoid negative values for epsilon.
#define MON_WITH_LOCAL_EXTREMUM_EPS(aminus, a, aplus, C) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  && \
         D2aC = (a(i-1) + a(i+1)) - 2.0d0*a(i)                                &&\
         D2aL = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                              &&\
         D2aR = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                              &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) + daplus * rhi     &&\
                  aminus(i) = a(i) - daminus * rhi    &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
               endif &&\
            else &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) &&\
                  aminus(i) = a(i) &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) &&\
               endif &&\
            endif &&\
         endif &&\
         trivial_rp(i-1) = .false.                                 &&\
         trivial_rp(i) = .false.                                   &&\
      else                                                         &&\
         if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
         else &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i)               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i)               &&\
            end if                                                    &&\
         endif &&\
         trivial_rp(i-1) = .false.  &&\
         trivial_rp(i) = .false.    &&\
      endif   


   if (apply_enhanced_ppm) then
      !! Constrain parabolic profiles, PPM 2011/2008
      if (PPM3) then
            do i = 3, nx - 2
                  MON_WITH_LOCAL_EXTREMUM(rhominus,rho,rhoplus, enhanced_ppm_C2)
                  MON_WITH_LOCAL_EXTREMUM(velxminus,vx,velxplus, enhanced_ppm_C2)
                  MON_WITH_LOCAL_EXTREMUM(velyminus,vy,velyplus, enhanced_ppm_C2)
                  MON_WITH_LOCAL_EXTREMUM(velzminus,vz,velzplus, enhanced_ppm_C2)
            end do
            if (poly .eq. 0) then
                  do i = 3, nx - 2
                     MON_WITH_LOCAL_EXTREMUM_EPS(epsminus,eps,epsplus, enhanced_ppm_C2)
                  end do
            end if
            
      else
            do i = 4, nx - 3
                  MON_WITH_LOCAL_EXTREMUM_STENCIL4(rhominus,rho,rhoplus, enhanced_ppm_C2, enhanced_ppm_C3)
                  MON_WITH_LOCAL_EXTREMUM_STENCIL4(velxminus,vx,velxplus, enhanced_ppm_C2, enhanced_ppm_C3)
                  MON_WITH_LOCAL_EXTREMUM_STENCIL4(velyminus,vy,velyplus, enhanced_ppm_C2, enhanced_ppm_C3)
                  MON_WITH_LOCAL_EXTREMUM_STENCIL4(velzminus,vz,velzplus, enhanced_ppm_C2, enhanced_ppm_C3)
            end do
            if (poly .eq. 0) then
                  do i = 4, nx - 3
                     MON_WITH_LOCAL_EXTREMUM_STENCIL4_EPS(epsminus,eps,epsplus, enhanced_ppm_C2, enhanced_ppm_C3)
                  end do
            endif
      end if 
      
      ! Apply flattening after constraining the parabolic profiles
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
            if (poly .eq. 0) then
               epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
               epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
            end if
         end do
      else  !!$ Really implement C&W, page 197; which requires stencil 4.
         do i = 4, nx - 3
            s=int(sign(one, -dpress(i)))
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
            if (poly .eq. 0) then
               epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
               epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
            end if
         end do
      end if
    
  endif

  if (.not.apply_enhanced_ppm) then
     !! Constrain parabolic profiles, PPM 1984
     do i = GRHydro_stencil, nx - GRHydro_stencil + 1
         ! original Colella&Woodward monotonicity preservation
         MON(rhominus,rho,rhoplus)
         MON(velxminus,vx,velxplus)
         MON(velyminus,vy,velyplus)
         MON(velzminus,vz,velzplus)
     end do
     if (poly .eq. 0) then
         do i = GRHydro_stencil, nx - GRHydro_stencil + 1
            MON(epsminus,eps,epsplus)
         end do
     end if
  endif

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
        rhominus(i-1)=rho(i)
        rhoplus(i-1)=rho(i)
        velxminus(i-1)=vx(i)
        velxplus(i-1)=vx(i)
        velyminus(i-1)=vy(i)
        velyplus(i-1)=vy(i)
        velzminus(i-1)=vz(i)
        velzplus(i-1)=vz(i)
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
        rhominus(i+1)=rho(i)
        rhoplus(i+1)=rho(i)
        velxminus(i+1)=vx(i)
        velxplus(i+1)=vx(i)
        velyminus(i+1)=vy(i)
        velyplus(i+1)=vy(i)
        velzminus(i+1)=vz(i)
        velzplus(i+1)=vz(i)
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

end subroutine SimplePPM_1d

!!!! routine for doing PPM to temperature
subroutine SimplePPM_temperature_1d(apply_enhanced_ppm,&     
     nx,dx,velx,temperature,press,&
     tempminus,&
     tempplus, hydro_excision_mask)

  USE GRHydro_Scalars
  USE GRHydro_Eigenproblem

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  logical :: apply_enhanced_ppm 

  CCTK_REAL, parameter :: one = 1.0d0

  CCTK_INT :: nx
  CCTK_REAL :: dx
  CCTK_REAL, dimension(nx) :: temperature, velx
  CCTK_REAL, dimension(nx) :: tempminus
  CCTK_REAL, dimension(nx) :: tempplus

  CCTK_INT :: i,s
  CCTK_REAL, dimension(nx) :: dtemp
  CCTK_REAL, dimension(nx) :: dmtemp
  CCTK_REAL, dimension(nx) :: press,dpress,tilde_flatten
  CCTK_REAL :: dpress2,dvel,w,flatten

  CCTK_INT, dimension(nx) :: hydro_excision_mask


  CCTK_REAL :: D2a, D2aL, D2aR, D2aC, D2aLim, rhi, daplus, daminus, D3a, D3aLL, D3aL, D3aR, D2aLL, D2aRR, D3aMin, D3aMax

  logical :: cond
  
#define STEEP(x,dx,dmx)                                              \
         if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
            dmx(i) = sign(one, dx(i)) *                                   \
               min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
                                 2.d0 * abs(x(i+1) - x(i)))              &&\
         else                                                           &&\
            dmx(i) = 0.d0                                                &&\
         end if


  if (.not.apply_enhanced_ppm) then
      !! This is the original PPM algorithm by Colella & Woodward 1984.

      !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
      !!$  This is the expression for an even grid.

     do i = 2, nx - 1
        dpress(i) = press(i+1) - press(i-1)
        ! the denominator is not necessary
        dtemp(i) = 0.5d0 * (temperature(i+1) - temperature(i-1))
     end do

      !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

     do i = 2, nx - 1
        STEEP(temperature, dtemp, dmtemp)
     enddo
      

!!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

     do i = 2, nx-2
        tempplus(i) = 0.5d0 * (temperature(i) + temperature(i+1)) + &
             (dmtemp(i) - dmtemp(i+1)) / 6.d0
        tempminus(i+1) = tempplus(i)
     enddo
  else
     !! This is the modified PPM algorithm by Colella & Sekora 2008 and McCorquodale & Colella 2011.
     !! This uses a better limiter based on second derivatives that preserves
     !! accuracy at local extrema. It also uses a higher-order interpolation polynomial.

      !!$  Initial boundary states (sixth order accurate). See (17) of Colella and Sekora 2008, p.7071
#define APPROX_AT_CELL_INTERFACE_STENCIL4T(a, ah)  \
      ah = 37.0d0/60.0d0*(a(i)+a(i+1)) - 2.0d0/15.0d0*(a(i-1)+a(i+2)) + 1.0d0/60.0d0*(a(i-2)+a(i+3))
      
      !!$  Initial boundary states (4th order accurate). See (16) of Colella and Sekora 2008, p.7071
#define APPROX_AT_CELL_INTERFACET(a, ah)  \
      ah = 7.0d0/12.0d0*(a(i)+a(i+1)) - 1.0d0/12.0d0*(a(i-1)+a(i+2))
      
#define LIMITT(a,ah,C, alim) \
      if ((min(a(i),a(i+1)) .le. ah) .and. (ah .le. max(a(i),a(i+1)))) then       &&\
         alim = ah   &&\
      else &&\
         D2a  = 3.0d0 * ((a(i)   + a(i+1)) - 2.0d0*ah    )                                     &&\
         D2aL =         ((a(i-1) + a(i+1)) - 2.0d0*a(i)  )                                     &&\
         D2aR =         ((a(i)   + a(i+2)) - 2.0d0*a(i+1))                                     &&\
         D2aLim = sign(one, D2a)*min(C*abs(D2aL), C*abs(D2aR), abs(D2a))/3.0d0             &&\
         if (D2a*D2aR .ge. 0 .and. D2a*D2aL .ge. 0) then                     &&\
            alim = 0.5d0*(a(i)+a(i+1)) - D2aLim &&\
         else                                                                             &&\
            alim = 0.5d0*(a(i)+a(i+1))                                                    &&\
         end if                                                                           &&\
      endif

#define LIMITT_EPS(a,ah,C, alim) \
      if ((min(a(i),a(i+1)) .le. ah) .and. (ah .le. max(a(i),a(i+1)))) then       &&\
         alim = ah   &&\
      else &&\
         D2a  = 3.0d0 * ((a(i)   + a(i+1)) - 2.0d0*ah    )                                     &&\
         D2aL =         ((a(i-1) + a(i+1)) - 2.0d0*a(i)  )                                     &&\
         D2aR =         ((a(i)   + a(i+2)) - 2.0d0*a(i+1))                                     &&\
         D2aLim = sign(one, D2a)*min(C*abs(D2aL), C*abs(D2aR), abs(D2a))/3.0d0             &&\
         if (D2a*D2aR .ge. 0 .and. D2a*D2aL .ge. 0 .and. abs(D2aLim) .lt. abs(0.5d0*(a(i)+a(i+1)))) then                     &&\
            alim = 0.5d0*(a(i)+a(i+1)) - D2aLim &&\
         else                                                                             &&\
            alim = 0.5d0*(a(i)+a(i+1))                                                    &&\
         end if&&\
      endif



      if (PPM3) then
      
      !! We initialize "plus" \equiv a_j+1/2 with (16) via APPROX_AT_CELL_INTERFACE, 
      !! then checking for (13) of Colella & Sekora 2008 and applying
      !! (18) and (19) if (13) is not satisfied. This is done with LIMIT.
         do i = 2, nx-2

            APPROX_AT_CELL_INTERFACET(temperature, tempplus(i))
            LIMITT_EPS(temperature, tempplus(i), enhanced_ppm_C2, tempplus(i))
            tempminus(i+1) = tempplus(i)

         enddo
      
      else

         !! Same as above but for 4 stencil points using (17) of Colella & Sekora 2008 as
         !! initial states.
         do i = 3, nx-3
         
            APPROX_AT_CELL_INTERFACE_STENCIL4T(temperature, tempplus(i))
            LIMITT_EPS(temperature, tempplus(i), enhanced_ppm_C2, tempplus(i))
            tempminus(i+1) = tempplus(i)
            
         enddo
        
      endif

      !! Finally compute pressure gradient needed for flattening and shock detection
      do i = 2, nx-1
         dpress(i) = press(i+1) - press(i-1)
      end do
      
   endif

!!$  Zone flattening. See appendix of C&W, p. 197-8.
  do i = 3, nx - 2
    dpress2 = press(i+2) - press(i-2)
    dvel = velx(i+1) - velx(i-1)
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


   if (.not.apply_enhanced_ppm) then
      ! In 1984 PPM, flattening is applied before constraining parabolic profiles.
      if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
         do i = 3, nx - 2
            flatten = tilde_flatten(i)
            tempplus(i) = flatten * tempplus(i) + (1.d0 - flatten) * temperature(i)
            tempminus(i) = flatten * tempminus(i) + (1.d0 - flatten) * temperature(i)
         end do
      else  !!$ Really implement C&W, page 197; which requires stencil 4.
         do i = 4, nx - 3
            s=int(sign(one, -dpress(i)))
            flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
            tempplus(i) = flatten * tempplus(i) + (1.d0 - flatten) * temperature(i)
            tempminus(i) = flatten * tempminus(i) + (1.d0 - flatten) * temperature(i)
         end do
      end if
   endif

!!$ Monotonicity. See (1.10) of C&W.
#define MONT(xminus,x,xplus)                                       \
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) )     \
        .and. ((xplus(i)-x(i))*(x(i)-xminus(i)) .le. 0.d0)) then&&\
      xminus(i) = x(i)                                          &&\
      xplus(i) = x(i)                                           &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) >                         \
                (xplus(i) - xminus(i))**2) then                 &&\
      xminus(i) = 3.d0 * x(i) - 2.d0 * xplus(i)                 &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) <                         \
               -(xplus(i) - xminus(i))**2) then                 &&\
      xplus(i) = 3.d0 * x(i) - 2.d0 * xminus(i)                 &&\
    end if                                                      &&\


!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34
!! This requires 4 stencil points.
#define MONT_WITH_LOCAL_EXTREMUM_STENCIL4(aminus, a, aplus, C, C3) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC  = (a(i-1) + a(i+1)) - 2.0d0*a(i)                               &&\
         D2aL  = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                             &&\
         D2aLL = (a(i-3) + a(i-1)) - 2.0d0*a(i-2)                             &&\
         D2aR  = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                             &&\
         D2aRR = (a(i+1) + a(i+3)) - 2.0d0*a(i+2)                             &&\
         D3a = D2aR - D2aC   &&\
         D3aL = D2aC - D2aL  &&\
         D3aR = D2aRR - D2aR &&\
         D3aLL = D2aL - D2aLL &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            D3aMin = min(D3aLL, D3aL, D3a, D3aR)   &&\
            D3aMax = max(D3aLL, D3aL, D3a, D3aR)   &&\
            if (C3 * max(abs(D3aMin), abs(D3aMax)) .le. D3aMax-D3aMin) then &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) + daplus * rhi     &&\
                     aminus(i) = a(i) - daminus * rhi    &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
                  endif &&\
            endif &&\
         endif &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
      endif

!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This contains an additional limiter in case the correction becomes larger than
!! the corrected value. This is to avoid negative values for epsilon.
!! This requires 4 stencil points.
#define MONT_WITH_LOCAL_EXTREMUM_STENCIL4_EPS(aminus, a, aplus, C, C3) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC  = (a(i-1) + a(i+1)) - 2.0d0*a(i)                               &&\
         D2aL  = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                             &&\
         D2aLL = (a(i-3) + a(i-1)) - 2.0d0*a(i-2)                             &&\
         D2aR  = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                             &&\
         D2aRR = (a(i+1) + a(i+3)) - 2.0d0*a(i+2)                             &&\
         D3a = D2aR - D2aC   &&\
         D3aL = D2aC - D2aL  &&\
         D3aR = D2aRR - D2aR &&\
         D3aLL = D2aL - D2aLL &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            D3aMin = min(D3aLL, D3aL, D3a, D3aR)   &&\
            D3aMax = max(D3aLL, D3aL, D3a, D3aR)   &&\
            if (C3 * max(abs(D3aMin), abs(D3aMax)) .le. D3aMax-D3aMin) then &&\
               if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) + daplus * rhi     &&\
                     aminus(i) = a(i) - daminus * rhi    &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
                  endif &&\
               else &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) &&\
                     aminus(i) = a(i) &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) &&\
                  endif &&\
               endif &&\
            endif &&\
         endif &&\
      else                                                         &&\
         if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
         else &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i)               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i)               &&\
            end if                                                    &&\
         endif &&\
      endif


!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
#define MONT_WITH_LOCAL_EXTREMUM(aminus, a, aplus, C) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC = (a(i-1) + a(i+1)) - 2.0d0*a(i)                                &&\
         D2aL = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                              &&\
         D2aR = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                              &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) + daplus * rhi     &&\
                  aminus(i) = a(i) - daminus * rhi    &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
               endif &&\
         endif &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
      endif   
      

!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
!! This contains an additional limiter in case the correction becomes larger than
!! the corrected value. This is to avoid negative values for epsilon.
#define MONT_WITH_LOCAL_EXTREMUM_EPS(aminus, a, aplus, C) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC = (a(i-1) + a(i+1)) - 2.0d0*a(i)                                &&\
         D2aL = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                              &&\
         D2aR = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                              &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) + daplus * rhi     &&\
                  aminus(i) = a(i) - daminus * rhi    &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
               endif &&\
            else &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) &&\
                  aminus(i) = a(i) &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) &&\
               endif &&\
            endif &&\
         endif &&\
      else                                                         &&\
         if (abs(daplus) .le. abs(a(i)) .and. abs(daminus) .le. abs(a(i))) then &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
         else &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i)               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i)               &&\
            end if                                                    &&\
         endif &&\
      endif   


   if (apply_enhanced_ppm) then
      !! Constrain parabolic profiles, PPM 2011/2008
      if (PPM3) then
            do i = 3, nx - 2
               MONT_WITH_LOCAL_EXTREMUM_EPS(tempminus,temperature,tempplus, enhanced_ppm_C2)
            enddo
      else
            do i = 4, nx - 3
               MONT_WITH_LOCAL_EXTREMUM_STENCIL4_EPS(tempminus,temperature,tempplus, enhanced_ppm_C2, enhanced_ppm_C3)
            enddo
      end if 
      
      ! Apply flattening after constraining the parabolic profiles
      if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
         do i = 3, nx - 2
            flatten = tilde_flatten(i)
            tempplus(i) = flatten * tempplus(i) + (1.d0 - flatten) * temperature(i)
            tempminus(i) = flatten * tempminus(i) + (1.d0 - flatten) * temperature(i)
         end do
      else  !!$ Really implement C&W, page 197; which requires stencil 4.
         do i = 4, nx - 3
            s=int(sign(one, -dpress(i)))
            flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
            tempplus(i) = flatten * tempplus(i) + (1.d0 - flatten) * temperature(i)
            tempminus(i) = flatten * tempminus(i) + (1.d0 - flatten) * temperature(i)
         end do
      end if
    
  endif

  if (.not.apply_enhanced_ppm) then
     !! Constrain parabolic profiles, PPM 1984
     do i = GRHydro_stencil, nx - GRHydro_stencil + 1
         ! original Colella&Woodward monotonicity preservation
        MONT(tempminus,temperature,tempplus)
     enddo
  endif

  !!$ excision
  do i = 1, nx
    if (GRHydro_enable_internal_excision /= 0 .and. &
        (hydro_excision_mask(i) .ne. 0)) then
       ! do nothing
    else
      !!$ Do not optimize cond away by combining the 'if's. Fortran does not
      !!$  have to follow the order of sub-expressions given here and might
      !!$  access outside the array range
      cond = .false.
      if (i .gt. 1 .and. GRHydro_enable_internal_excision /= 0) then
        cond = hydro_excision_mask(i-1) .ne. 0
      end if
      if (cond) then
        tempminus(i) = temperature(i)
        tempplus(i) = temperature(i)
        tempminus(i-1) = temperature(i)
        tempplus(i-1) = temperature(i)
      else
        cond = .false.
        if ((i.gt.2) .and. (i.lt.nx) .and. GRHydro_enable_internal_excision /= 0) then
          cond = (ppm_mppm .eq. 0) .and. (hydro_excision_mask(i-2) .ne. 0)
        end if
        if (cond) then
          call PPM_TVD(temperature(i-1), temperature(i), temperature(i+1), tempminus(i), tempplus(i))
        end if
      end if
      cond = .false.
      if (i .lt. nx .and. GRHydro_enable_internal_excision /= 0) then
        cond = hydro_excision_mask(i+1) .ne. 0
      end if
      if (cond) then
        tempminus(i) = temperature(i)
        tempplus(i) = temperature(i)
        tempminus(i-1) = temperature(i)
        tempplus(i-1) = temperature(i)
      else
        cond = .false.
        if ((i.lt.nx-1) .and. (i.gt.1) .and. GRHydro_enable_internal_excision /= 0) then
          cond = (ppm_mppm .eq. 0) .and. (hydro_excision_mask(i+2) .ne. 0)
        end if
        if (cond) then
          call PPM_TVD(temperature(i-1), temperature(i), temperature(i+1), tempminus(i), tempplus(i))
        end if
      end if
    end if
  end do
  
  
  
return

end subroutine SimplePPM_temperature_1d





subroutine SimplePPM_tracer_1d(apply_enhanced_ppm,&
     nx,dx,rho,velx,vely,velz, &
     tracer,tracerminus,tracerplus,press)

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  logical :: apply_enhanced_ppm

  CCTK_REAL, parameter :: one = 1

  CCTK_INT :: nx
  CCTK_REAL :: dx
  CCTK_REAL, dimension(nx) :: rho,velx,vely,velz
  CCTK_REAL, dimension(nx,number_of_tracers) :: tracer,tracerminus,tracerplus
  CCTK_REAL :: tracerflatomega


  CCTK_INT :: i,s,itracer
  CCTK_REAL, dimension(nx) :: press,dpress,tilde_flatten
  CCTK_REAL, dimension(nx,number_of_tracers) :: dmtracer, dtracer, tracerflat!, d2tracer
  CCTK_REAL :: dpress2,w,flatten,dvel

!!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
!!$  This is the expression for an even grid.


  do i = 2, nx - 1
    dpress(i) = press(i+1) - press(i-1)
  end do

  do itracer=1,number_of_tracers
     do i = 2, nx - 1
        dtracer(i,itracer) = 0.5d0 * (tracer(i+1,itracer) - tracer(i-1,itracer))
!        d2tracer(i,itracer) = (tracer(i+1) - 2.d0 * tracer(i) + tracer(i-1))! / 6.d0 / dx / dx
!    ! since we use d2tracer only for the condition d2tracer(i+1)*d2tracer(i-1)<0 
!    ! the denominator is not necessary
     enddo
  enddo

!!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

  do itracer=1,number_of_tracers
     do i = 2, nx - 1
        if( (tracer(i+1,itracer) - tracer(i,itracer)) * &
             (tracer(i,itracer) - tracer(i-1,itracer)) > 0.0d0 ) then
           dmtracer(i,itracer) = sign(one,dtracer(i,itracer)) * &
                min(abs(dtracer(i,itracer)), 2.0d0 * &
                abs(tracer(i,itracer) - tracer(i-1,itracer)), &
                2.0d0 * abs(tracer(i+1,itracer) - tracer(i,itracer)))
        else
           dmtracer(i,itracer) = 0.0d0
        endif
     end do
  enddo

!!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

    do itracer=1,number_of_tracers
       do i = 2, nx - 2
          tracerplus(i,itracer) = 0.5d0 * (tracer(i,itracer) + tracer(i+1,itracer)) + &
               (dmtracer(i,itracer) - dmtracer(i+1,itracer)) / 6.d0
          tracerminus(i+1,itracer) = tracerplus(i,itracer)
       enddo
    enddo
    

!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which mat be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even. 
!!$Note also that we do not have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).

!!! We might play around with this for the tracer. CURRENTLY TURNED OFF

#if 0
  if (ppm_detect .eq. 1000) then
     do itracer=1,number_of_tracers
        
        do i = 3, nx - 2
           if ( (dtracer(i+1,itracer)*dtracer(i-1,itracer) > 0.d0) & !make sure this is not an extremum
           .and.(abs(tracer(i+1,itracer)-tracer(i-1,itracer)) - & !this is to prevent steepening
                !of relatively small composition jumps
                ppm_epsilon_shock * min(tracer(i+1,itracer), tracer(i-1,itracer)) > 0.d0 )  & 
                .and. & ! the actual criterion from Plewa & Mueller
                 ((tracer(i+1,itracer)-tracer(i-1,itracer)) / &
                 (tracer(i+2,itracer)-tracer(i-2,itracer)) > ppm_omega1 ) ) then

           etatilde = (tracer(i-2,itracer) - tracer(i+2,itracer) + & 
                4.d0 * dtracer(i,itracer)) / (dtracer(i,itracer) * 12.d0)

           write(*,*) "Additional Steepening in Zone: ",i

        else
           etatilde = 0.d0
        end if
        eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
        if (ppm_k0 * abs(dtracer(i,itracer)) * min(press(i-1),press(i+1)) < &
             abs(dpress(i)) * min(tracer(i-1,itracer), tracer(i+1,itracer))) then
           eta = 0.d0
        end if
        tracerminus(i,itracer) = tracerminus(i,itracer) * (1.d0 - eta) + &
             (tracer(i-1,itracer) + 0.5d0 * dmtracer(i-1,itracer)) * eta
        tracerplus(i,itracer) = tracerplus(i,itracer) * (1.d0 - eta) + &
             (tracer(i+1,itracer) - 0.5d0 * dmtracer(i+1,itracer)) * eta
     end do

  enddo

  end if
#endif

!!$  Zone flattening. See appendix of C&W, p. 197-8.

  do i = 3, nx - 2
    dpress2 = press(i+2) - press(i-2)
    dvel = velx(i+1) - velx(i-1)
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

  if (PPM3) then
     do itracer=1,number_of_tracers
        do i = 3, nx - 2
           flatten = tilde_flatten(i)
           tracerplus(i,itracer) = flatten * tracerplus(i,itracer) & 
                + (1.d0 - flatten) * tracer(i,itracer)
           tracerminus(i,itracer) = flatten * tracerminus(i,itracer) & 
                + (1.d0 - flatten) * tracer(i,itracer)
        end do
     enddo
  else  !!$ Really implement C&W, page 197; which requires stencil 4.
     do itracer=1,number_of_tracers
        do i = 4, nx - 3
           s=int(sign(one, -dpress(i)))
           flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
           tracerplus(i,itracer) = flatten * tracerplus(i,itracer) + &
                (1.d0 - flatten) * tracer(i,itracer)
           tracerminus(i,itracer) = flatten * tracerminus(i,itracer) & 
                + (1.d0 - flatten) * tracer(i,itracer)
        end do
     enddo
  end if


!! Additional flattening a la Plewa & Mueller                                                                 

#if 1
  do itracer=1,number_of_tracers
     do i = 2, nx - 1
        if ( ( tracer(i+1,itracer) - tracer(i,itracer) ) * &
           ( tracer(i,itracer) - tracer(i-1,itracer) ) < 0.0d0 ) then
           tracerflat(i,itracer) = 1.0d0
        else
           tracerflat(i,itracer) = 0.0d0
        endif
     enddo
  enddo

  do itracer=1,number_of_tracers
     do i = 3, nx -2

        tracerflatomega = 0.5d0 * max(tracerflat(i-1,itracer),2.0d0*tracerflat(i,itracer), &
             tracerflat(i+1,itracer)) * ppm_omega_tracer

        tracerplus(i,itracer) = tracerflatomega*tracer(i,itracer) + &
             (1.0d0 - tracerflatomega)*tracerplus(i,itracer)

        tracerminus(i,itracer) = tracerflatomega*tracer(i,itracer) + &
             (1.0d0 - tracerflatomega)*tracerminus(i,itracer)

     enddo
  enddo


#endif

!!$ Monotonicity. See (1.10) of C&W.                                                                          


  do itracer=1,number_of_tracers
     do i = GRHydro_stencil, nx - GRHydro_stencil + 1
        if (((tracerplus(i,itracer)-tracer(i,itracer))*      &
           (tracer(i,itracer)-tracerminus(i,itracer)) .le. 0.d0)) then
           tracerminus(i,itracer) = tracer(i,itracer)
           tracerplus(i,itracer) = tracer(i,itracer)
        else if ((tracerplus(i,itracer) - tracerminus(i,itracer)) * (tracer(i,itracer) - 0.5d0 * &
             (tracerplus(i,itracer) + tracerminus(i,itracer))) > &
           (tracerplus(i,itracer) - tracerminus(i,itracer))**2 / 6.d0) then
           tracerminus(i,itracer) = 3.d0 * tracer(i,itracer) - 2.d0 * tracerplus(i,itracer)
        else if ((tracerplus(i,itracer) - tracerminus(i,itracer)) * (tracer(i,itracer) - 0.5d0 * &
             (tracerplus(i,itracer) + tracerminus(i,itracer))) <  &
           -(tracerplus(i,itracer) - tracerminus(i,itracer))**2 / 6.d0 ) then
           tracerplus(i,itracer) = 3.d0 * tracer(i,itracer) - 2.d0 * tracerminus(i,itracer)
        end if
     end do
  enddo



end subroutine SimplePPM_tracer_1d


subroutine SimplePPM_ye_1d(apply_enhanced_ppm,&
     nx,dx,rho,velx,vely,velz, &
     Y_e,Y_e_minus,Y_e_plus,press)

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  logical :: apply_enhanced_ppm

  CCTK_REAL, parameter :: one = 1

  CCTK_INT :: nx
  CCTK_REAL :: dx
  CCTK_REAL, dimension(nx) :: rho,velx,vely,velz
  CCTK_REAL, dimension(nx) :: Y_e,Y_e_minus,Y_e_plus


  CCTK_INT :: i,s
  CCTK_REAL, dimension(nx) :: press,dpress,tilde_flatten
  CCTK_REAL, dimension(nx) :: dmY_e, dY_e
  CCTK_REAL :: dpress2,w,flatten,dvel

  CCTK_REAL :: D2a, D2aL, D2aR, D2aC, D2aLim, rhi, daplus, daminus, D3a, D3aLL, D3aL, D3aR, D2aLL, D2aRR, D3aMin, D3aMax

  do i = 2, nx - 1
       dpress(i) = press(i+1) - press(i-1)
  end do

  
  if (.not.apply_enhanced_PPM) then
  
   !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
   !!$  This is the expression for an even grid.


      do i = 2, nx - 1
            dY_e(i) = 0.5d0 * (Y_e(i+1) - Y_e(i-1))
      !        d2Y_e(i,iY_e) = (Y_e(i+1) - 2.d0 * Y_e(i) + Y_e(i-1))! / 6.d0 / dx / dx
      !    ! since we use d2Y_e only for the condition d2Y_e(i+1)*d2Y_e(i-1)<0 
      !    ! the denominator is not necessary
      enddo


      !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

      do i = 2, nx - 1
         if( (Y_e(i+1) - Y_e(i)) * &
               (Y_e(i) - Y_e(i-1)) > 0.0d0 ) then
            dmY_e(i) = sign(one,dY_e(i)) * &
                  min(abs(dY_e(i)), 2.0d0 * &
                  abs(Y_e(i) - Y_e(i-1)), &
                  2.0d0 * abs(Y_e(i+1) - Y_e(i)))
         else
            dmY_e(i) = 0.0d0
         endif
      end do

      !!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

      do i = 2, nx - 2
         Y_e_plus(i) = 0.5d0 * (Y_e(i) + Y_e(i+1)) + &
               (dmY_e(i) - dmY_e(i+1)) / 6.d0
         Y_e_minus(i+1) = Y_e_plus(i)
      enddo
   else
!! This is the modified PPM algorithm by Colella & Sekora 2008.
!! This uses a better limiter based on second derivatives that preserves
!! accuracy at local extrema. It also uses a higher-order interpolation polynomial.


      !!$  Initial boundary states (sixth order accurate). See (17) of Colella and Sekora, p.7071
#undef APPROX_AT_CELL_INTERFACE_STENCIL4
!!$  Initial boundary states (sixth order accurate). See (17) of Colella and Sekora 2008, p.7071
#define APPROX_AT_CELL_INTERFACE_STENCIL4(a, ah)  \
      ah = 37.0d0/60.0d0*(a(i)+a(i+1)) - 2.0d0/15.0d0*(a(i-1)+a(i+2)) + 1.0d0/60.0d0*(a(i-2)+a(i+3))

#undef APPROX_AT_CELL_INTERFACE
      !!$  Initial boundary states (4th order accurate). See (16) of Colella and Sekora 2008, p.7071
#define APPROX_AT_CELL_INTERFACE(a, ah)  \
      ah = 7.0d0/12.0d0*(a(i)+a(i+1)) - 1.0d0/12.0d0*(a(i-1)+a(i+2))


#undef LIMIT
#define LIMIT(a,ah,C, alim) \
      if ((min(a(i),a(i+1)) .le. ah) .and. (ah .le. max(a(i),a(i+1)))) then       &&\
         alim = ah   &&\
      else &&\
         D2a  = 3.0d0 * ((a(i)   + a(i+1)) - 2.0d0*ah    )                   &&\
         D2aL =         ((a(i-1) + a(i+1)) - 2.0d0*a(i)  )                   &&\
         D2aR =         ((a(i)   + a(i+2)) - 2.0d0*a(i+1))                   &&\
         if (D2a*D2aR .ge. 0 .and. D2a*D2aL .ge. 0) then                     &&\
            alim = 0.5d0*(a(i)+a(i+1)) - sign(one, D2a)*min(C*abs(D2aL), C*abs(D2aR), abs(D2a))/3.0d0 &&\
         else                                                                             &&\
            alim = 0.5d0*(a(i)+a(i+1))                                                    &&\
         end if                                                                           &&\
      endif

      if (PPM3) then
         !! We initialize "plus" \equiv a_j+1/2 with (16) via APPROX_AT_CELL_INTERFACE, 
         !! then checking for (13) of Colella & Sekora 2008 and applying
         !! (18) and (19) if (13) is not satisfied. This is done with LIMIT.
         !! There, we also switch to lower order if we are near atmosphere values.
         do i = 2, nx-2
            APPROX_AT_CELL_INTERFACE(Y_e, Y_e_plus(i))
            LIMIT(Y_e, Y_e_plus(i), enhanced_ppm_C2, Y_e_plus(i))
            Y_e_minus(i+1) = Y_e_plus(i)
         end do
         
      else

         !! Same as above but for 4 stencil points using (17) of Colella & Sekora 2008 as
         !! initial states.
         do i = 3, nx-3
            APPROX_AT_CELL_INTERFACE_STENCIL4(Y_e, Y_e_plus(i))
            LIMIT(Y_e, Y_e_plus(i), enhanced_ppm_C2, Y_e_plus(i))
            Y_e_minus(i+1) = Y_e_plus(i)
         end do
      endif
    
   end if



!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which mat be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even. 
!!$Note also that we do not have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).


!!$  Zone flattening. See appendix of C&W, p. 197-8.

  do i = 3, nx - 2
    dpress2 = press(i+2) - press(i-2)
    dvel = velx(i+1) - velx(i-1)
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

  if (.not.apply_enhanced_ppm) then
      if (PPM3) then
         do i = 3, nx - 2
            flatten = tilde_flatten(i)
            Y_e_plus(i) = flatten * Y_e_plus(i) & 
                  + (1.d0 - flatten) * Y_e(i)
            Y_e_minus(i) = flatten * Y_e_minus(i) & 
                  + (1.d0 - flatten) * Y_e(i)
         end do
      else  !!$ Really implement C&W, page 197; which requires stencil 4.
         do i = 4, nx - 3
            s=int(sign(one, -dpress(i)))
            flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
            Y_e_plus(i) = flatten * Y_e_plus(i) + &
                  (1.d0 - flatten) * Y_e(i)
            Y_e_minus(i) = flatten * Y_e_minus(i) & 
                  + (1.d0 - flatten) * Y_e(i)
         end do
      end if
  endif
  

#undef MON_WITH_LOCAL_EXTREMUM_STENCIL4
!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34
!! This requires 4 stencil points.
#define MON_WITH_LOCAL_EXTREMUM_STENCIL4(aminus, a, aplus, C, C3) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC  = (a(i-1) + a(i+1)) - 2.0d0*a(i)                               &&\
         D2aL  = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                             &&\
         D2aLL = (a(i-3) + a(i-1)) - 2.0d0*a(i-2)                             &&\
         D2aR  = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                             &&\
         D2aRR = (a(i+1) + a(i+3)) - 2.0d0*a(i+2)                             &&\
         D3a = D2aR - D2aC   &&\
         D3aL = D2aC - D2aL  &&\
         D3aR = D2aRR - D2aR &&\
         D3aLL = D2aL - D2aLL &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            D3aMin = min(D3aLL, D3aL, D3a, D3aR)   &&\
            D3aMax = max(D3aLL, D3aL, D3a, D3aR)   &&\
            if (C3 * max(abs(D3aMin), abs(D3aMax)) .le. D3aMax-D3aMin) then &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) + daplus * rhi     &&\
                     aminus(i) = a(i) - daminus * rhi    &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
                  endif &&\
            endif &&\
         endif &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
      endif
            




#undef MON_WITH_LOCAL_EXTREMUM
!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
#define MON_WITH_LOCAL_EXTREMUM(aminus, a, aplus, C) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC = (a(i-1) + a(i+1)) - 2.0d0*a(i)                                &&\
         D2aL = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                              &&\
         D2aR = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                              &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) + daplus * rhi     &&\
                  aminus(i) = a(i) - daminus * rhi    &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
               endif &&\
         endif      &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) - 2.0d0*(aminus(i) - a(i))               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*(aplus(i) - a(i))               &&\
            end if  &&\
      endif





  if (apply_enhanced_ppm) then
      ! Constrain parabolic profiles
      if (PPM3) then
         do i = GRHydro_stencil, nx - GRHydro_stencil + 1
               MON_WITH_LOCAL_EXTREMUM(Y_e_minus,Y_e,Y_e_plus, enhanced_ppm_C2)
         end do
      else
         do i = GRHydro_stencil, nx - GRHydro_stencil + 1
               MON_WITH_LOCAL_EXTREMUM_STENCIL4(Y_e_minus,Y_e,Y_e_plus, enhanced_ppm_C2, enhanced_ppm_C3)
         end do
      endif
      
      ! In PPM 2011/2008, flattening is applied after constraining parabolic profiles.
      if (PPM3) then
         do i = 3, nx - 2
            flatten = tilde_flatten(i)
            Y_e_plus(i) = flatten * Y_e_plus(i) & 
                  + (1.d0 - flatten) * Y_e(i)
            Y_e_minus(i) = flatten * Y_e_minus(i) & 
                  + (1.d0 - flatten) * Y_e(i)
         end do
      else  !!$ Really implement C&W, page 197; which requires stencil 4.
         do i = 4, nx - 3
            s=int(sign(one, -dpress(i)))
            flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
            Y_e_plus(i) = flatten * Y_e_plus(i) + &
                  (1.d0 - flatten) * Y_e(i)
            Y_e_minus(i) = flatten * Y_e_minus(i) & 
                  + (1.d0 - flatten) * Y_e(i)
         end do
      end if
      
  else
      !!$ Monotonicity. See (1.10) of C&W.
      do i = GRHydro_stencil, nx - GRHydro_stencil + 1
         if (((Y_e_plus(i)-Y_e(i))*      &
               (Y_e(i)-Y_e_minus(i)) .le. 0.d0)) then
            Y_e_minus(i) = Y_e(i)
            Y_e_plus(i) = Y_e(i)
         else if ((Y_e_plus(i) - Y_e_minus(i)) * (Y_e(i) - 0.5d0 * &
               (Y_e_plus(i) + Y_e_minus(i))) > &
               (Y_e_plus(i) - Y_e_minus(i))**2 / 6.d0) then
            Y_e_minus(i) = 3.d0 * Y_e(i) - 2.0d0 * Y_e_plus(i)
         else if ((Y_e_plus(i) - Y_e_minus(i)) * (Y_e(i) - 0.5d0 * &
               (Y_e_plus(i) + Y_e_minus(i))) <  &
               -(Y_e_plus(i) - Y_e_minus(i))**2 / 6.0d0 ) then
            Y_e_plus(i) = 3.d0 * Y_e(i) - 2.d0 * Y_e_minus(i)
         end if
      end do
  end if


end subroutine SimplePPM_ye_1d
