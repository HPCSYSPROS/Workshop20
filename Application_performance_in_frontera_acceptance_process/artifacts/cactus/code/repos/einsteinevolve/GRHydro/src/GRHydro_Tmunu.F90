 /*@@
   @file      GRHydro_Tmunu.F90
   @date      Thu Apr 16 19:38:40 2009
   @author    Ian Hawke
   @histpry
    Apr. 2009: Luca Baiotti copied and adapted for the Tmunu-thorn mechanism the original include file
   @desc 
     The calculation of the stress energy tensor.
     The version used here was worked out by Miguel Alcubierre. I
     think it was an extension of the routine from GR3D, written
     by Mark Miller.
     C version added by Ian Hawke.

     Lower components of the stress-energy tensor obtained from
     the hydro variables.  The components are given by:

     T      =  rho h  u   u    +  P  g
      mu nu            mu  nu         mu nu

     where rho is the energy density of the fluid, h the enthalpy
     and P the pressure.  The enthalpy is given in terms of the
     basic variables as:

     h  =  1  +  e  +  P/rho

     with e the internal energy (eps here).

     In the expresion for T_{mu,nu} we also have the four-velocity
     of the fluid given by (v_i is the 3-velocity field):

                                 i
     u  =  W ( - alpha +  v  beta  )
      0                    i

     u  =  W v
      i       i
                                                i  -1/2
     with W the Lorentz factor:   W = ( 1 -  v v  )
                                              i

     and where alpha and beta are the lapse and shift vector.

     Finally, the 4 metric is given by

                   2             i
     g   =  - alpha  + beta  beta
      00                   i

     g   =  beta
      0i        i


     g   =  gamma      (the spatial metric)
      ij        ij


   @enddesc 
 @@*/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)

 subroutine GRHydro_Tmunu(CCTK_ARGUMENTS)

   implicit none

   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   DECLARE_CCTK_FUNCTIONS

   CCTK_REAL velxlow, velylow, velzlow
   CCTK_REAL betaxlow, betaylow, betazlow, beta2
   CCTK_REAL utlow, uxlow, uylow, uzlow
   CCTK_REAL rhoenthalpy
   CCTK_REAL ut,ux,uy,uz,bst,bsx,bsy,bsz,bs2
   CCTK_REAL dampfac
   CCTK_INT i,j,k

!!$ Damping factor
   dampfac = 1.0

   !$OMP PARALLEL DO PRIVATE(i,j,k,velxlow, velylow, velzlow,&
   !$OMP betaxlow, betaylow, betazlow, beta2, utlow, uxlow, uylow, uzlow,&
   !$OMP rhoenthalpy, ut,ux,uy,uz,bst,bsx,bsy,bsz,bs2,dampfac)
   do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)

          velxlow = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k)
          velylow = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k)
          velzlow = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k)

!!$       Calculate lower components and square of shift vector.

          betaxlow = gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k) + gxz(i,j,k)*betaz(i,j,k)
          betaylow = gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k) + gyz(i,j,k)*betaz(i,j,k)
          betazlow = gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k) + gzz(i,j,k)*betaz(i,j,k)
          
          beta2 = betax(i,j,k)*betaxlow + betay(i,j,k)*betaylow + betaz(i,j,k)*betazlow 
           
!!$       Calculate the specific relativistic enthalpy times rho times the
!!$       square of the lorentz factor.

         rhoenthalpy = w_lorentz(i,j,k)**2*(rho(i,j,k)*(1.0d0 + eps(i,j,k)) + press(i,j,k))

!!$       Calculate lower components of 4-velocity (without the Lorent factor).

         utlow = (-alp(i,j,k) + velx(i,j,k)*betaxlow + vely(i,j,k)*betaylow + velz(i,j,k)*betazlow)

         uxlow = velxlow
         uylow = velylow
         uzlow = velzlow


!!$      Initialize damping factor
         dampfac = 1.0
!!$      Apply tanh blending for Tmunu.
         if ((Tmunu_damping_radius_min .gt. 0) .and. (r(i,j,k) .gt. Tmunu_damping_radius_min)) then
            !  0.5 * (1.0 - tanh(4.0*(x-x0)/sigma0))
            if (r(i,j,k) .lt. Tmunu_damping_radius_max) then
               dampfac = 0.5d0 * (1.0d0 - tanh((8.0d0*r(i,j,k)-4.0d0*(Tmunu_damping_radius_max+Tmunu_damping_radius_min))/(Tmunu_damping_radius_max-Tmunu_damping_radius_min)))
            else
               dampfac = 0.0
               continue   ! no need to add anything to Tmunu at the current point (it's zero anyway!)
            endif
         else
            dampfac = 1.0
         endif
         
!!$       Calculate Tmunu (the lower components!).

         eTtt(i,j,k) = eTtt(i,j,k) + dampfac * (rhoenthalpy*utlow**2 + press(i,j,k)*(beta2 - alp(i,j,k)**2))

         eTtx(i,j,k) = eTtx(i,j,k) + dampfac * (rhoenthalpy*utlow*uxlow + press(i,j,k)*betaxlow)
         eTty(i,j,k) = eTty(i,j,k) + dampfac * (rhoenthalpy*utlow*uylow + press(i,j,k)*betaylow)
         eTtz(i,j,k) = eTtz(i,j,k) + dampfac * (rhoenthalpy*utlow*uzlow + press(i,j,k)*betazlow)

         eTxx(i,j,k) = eTxx(i,j,k) + dampfac * (rhoenthalpy*uxlow**2 + press(i,j,k)*gxx(i,j,k))
         eTyy(i,j,k) = eTyy(i,j,k) + dampfac * (rhoenthalpy*uylow**2 + press(i,j,k)*gyy(i,j,k))
         eTzz(i,j,k) = eTzz(i,j,k) + dampfac * (rhoenthalpy*uzlow**2 + press(i,j,k)*gzz(i,j,k))
         
         eTxy(i,j,k) = eTxy(i,j,k) + dampfac * (rhoenthalpy*uxlow*uylow + press(i,j,k)*gxy(i,j,k))
         eTxz(i,j,k) = eTxz(i,j,k) + dampfac * (rhoenthalpy*uxlow*uzlow + press(i,j,k)*gxz(i,j,k))
         eTyz(i,j,k) = eTyz(i,j,k) + dampfac * (rhoenthalpy*uylow*uzlow + press(i,j,k)*gyz(i,j,k))

       end do
     end do
   end do
   !$OMP END PARALLEL DO

   return
 
 end subroutine GRHydro_Tmunu
