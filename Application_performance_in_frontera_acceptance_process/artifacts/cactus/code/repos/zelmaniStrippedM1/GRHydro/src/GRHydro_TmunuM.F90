 /*@@
   @file      GRHydro_Tmunu.F90
   @date      Aug 30, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
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

     T      =  (rho h +b^2) u   u    +  (P+b^2/2)  g      - b   b
      mu nu                   mu  nu                mu nu     mu nu 

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
#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define bcomx(i,j,k) bcom(i,j,k,1)
#define bcomy(i,j,k) bcom(i,j,k,2)
#define bcomz(i,j,k) bcom(i,j,k,3)

 subroutine GRHydro_TmunuM(CCTK_ARGUMENTS)

   implicit none

   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   DECLARE_CCTK_FUNCTIONS

   CCTK_REAL :: velxlow, velylow, velzlow
   CCTK_REAL :: betaxlow, betaylow, betazlow, beta2
   CCTK_REAL :: Bvecxlow,Bvecylow,Bveczlow
   CCTK_REAL :: bdotv,b2,bxlow,bylow,bzlow,btlow,dum1,dum2
   CCTK_REAL :: utlow,rhohstarw2,pstar
   CCTK_REAL :: bdotbeta,vdotbeta
   CCTK_INT :: i,j,k


   !$OMP PARALLEL DO PRIVATE(i,j,k,velxlow, velylow, velzlow,&
   !$OMP Bvecxlow,Bvecylow,Bveczlow, bdotv,dum1,dum2,b2,bxlow,bylow,bzlow,&
   !$OMP betaxlow, betaylow, betazlow, beta2, bdotbeta,vdotbeta,utlow, btlow,&
   !$OMP rhohstarw2,pstar)

   do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)

          ! need separate dum1, dum2 b/c of Fortrans aliasing rules
          call calc_vlow_blow(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
               gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
               velx(i,j,k),vely(i,j,k),velz(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k), &
               velxlow,velylow,velzlow,Bvecxlow,Bvecylow,Bveczlow, &
               bdotv,b2,dum1,dum2,bxlow,bylow,bzlow)
          
!!$       Calculate lower components and square of shift vector.
          
          
          betaxlow = gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k) + gxz(i,j,k)*betaz(i,j,k)
          betaylow = gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k) + gyz(i,j,k)*betaz(i,j,k)
          betazlow = gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k) + gzz(i,j,k)*betaz(i,j,k)
          beta2 = betax(i,j,k)*betaxlow + betay(i,j,k)*betaylow + betaz(i,j,k)*betazlow 
          
          bdotbeta = betaxlow*Bvecx(i,j,k)+betaylow*Bvecy(i,j,k)+betazlow*Bvecz(i,j,k)
          vdotbeta = betaxlow*velx(i,j,k)+betaylow*vely(i,j,k)+betazlow*velz(i,j,k)

!!$   u0 low is missing the w_lorentz factor (see below)!!
          utlow = -1.d0*alp(i,j,k) + vdotbeta
          
          btlow = -1.0d0*w_lorentz(i,j,k)*alp(i,j,k)*bdotv + &
               bdotbeta/w_lorentz(i,j,k) + w_lorentz(i,j,k)*bdotv*vdotbeta

!!$       Calculate the specific relativistic enthalpy times rho + the mag. field contribution times the
!!$       square of the lorentz factor.

         rhohstarw2 = w_lorentz(i,j,k)**2*(rho(i,j,k)*(1.0d0 + eps(i,j,k)) + press(i,j,k) + b2)
         pstar =  press(i,j,k)+0.5d0*b2

!!$       Calculate lower components of 4-velocity (without the Lorent factor).
!!$         uxlow = velxlow
!!$         uylow = velylow
!!$         uzlow = velzlow

!!$       Calculate Tmunu (the lower components!).

         eTtt(i,j,k) = eTtt(i,j,k) + rhohstarw2*utlow**2 + pstar*(beta2 - alp(i,j,k)**2) - btlow**2

         eTtx(i,j,k) = eTtx(i,j,k) + rhohstarw2*utlow*velxlow + pstar*betaxlow - btlow*bxlow
         eTty(i,j,k) = eTty(i,j,k) + rhohstarw2*utlow*velylow + pstar*betaylow - btlow*bylow
         eTtz(i,j,k) = eTtz(i,j,k) + rhohstarw2*utlow*velzlow + pstar*betazlow - btlow*bzlow

         eTxx(i,j,k) = eTxx(i,j,k) + rhohstarw2*velxlow**2 + pstar*gxx(i,j,k) - bxlow**2
         eTyy(i,j,k) = eTyy(i,j,k) + rhohstarw2*velylow**2 + pstar*gyy(i,j,k) - bylow**2
         eTzz(i,j,k) = eTzz(i,j,k) + rhohstarw2*velzlow**2 + pstar*gzz(i,j,k) - bzlow**2
         
         eTxy(i,j,k) = eTxy(i,j,k) + rhohstarw2*velxlow*velylow + pstar*gxy(i,j,k) - bxlow*bylow
         eTxz(i,j,k) = eTxz(i,j,k) + rhohstarw2*velxlow*velzlow + pstar*gxz(i,j,k) - bxlow*bzlow
         eTyz(i,j,k) = eTyz(i,j,k) + rhohstarw2*velylow*velzlow + pstar*gyz(i,j,k) - bylow*bzlow

         if(calculate_bcom .ne. 0) then
           bcom_sq(i,j,k) = b2 
           bcom0(i,j,k) = w_lorentz(i,j,k)*bdotv/alp(i,j,k)
           bcomx(i,j,k) = Bvecx(i,j,k)/w_lorentz(i,j,k) + bcom0(i,j,k)*(alp(i,j,k)*velx(i,j,k)-betax(i,j,k))
           bcomy(i,j,k) = Bvecy(i,j,k)/w_lorentz(i,j,k) + bcom0(i,j,k)*(alp(i,j,k)*vely(i,j,k)-betay(i,j,k))
           bcomz(i,j,k) = Bvecz(i,j,k)/w_lorentz(i,j,k) + bcom0(i,j,k)*(alp(i,j,k)*velz(i,j,k)-betaz(i,j,k))
         endif 

       end do
     end do
   end do
   !$OMP END PARALLEL DO

   return
 
 end subroutine GRHydro_TmunuM
