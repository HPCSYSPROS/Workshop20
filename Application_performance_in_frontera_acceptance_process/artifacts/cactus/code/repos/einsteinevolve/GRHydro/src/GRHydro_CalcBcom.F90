 /*@@
   @file      GRHydro_CalcBcom.F90
   @date      Jan 01, 2012
   @author    Bruno Mundim
   @histpry
    Based on GRHydro_TmunuM.F90 file.
   @desc 
     The calculation of the magnetic pressure and the comoving magnetic
     field.
   @enddesc 
 @@*/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define bcomx(i,j,k) bcom(i,j,k,1)
#define bcomy(i,j,k) bcom(i,j,k,2)
#define bcomz(i,j,k) bcom(i,j,k,3)

 subroutine GRHydro_CalcBcom(CCTK_ARGUMENTS)

   implicit none

   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_FUNCTIONS

   CCTK_REAL :: velxlow, velylow, velzlow
   CCTK_REAL :: Bvecxlow,Bvecylow,Bveczlow
   CCTK_REAL :: bdotv,b2,bxlow,bylow,bzlow,dum1,dum2
   CCTK_INT :: i,j,k
   integer :: timelevels

   call CCTK_ActiveTimeLevels(timelevels, cctkGH, "GRHydro::bcom")
   if(timelevels.eq.0) then
     call CCTK_ERROR("No storage for GRHydro::bcom")
     STOP
   end if

! Only compute if Tmunu is not computed! (Bcom is also computed with Tmunu)
   if (stress_energy_state .ne. 0) then
      return
   endif

   !$OMP PARALLEL DO PRIVATE(i,j,k,velxlow, velylow, velzlow,&
   !$OMP Bvecxlow,Bvecylow,Bveczlow, bdotv,dum1,dum2,b2,bxlow,bylow,bzlow)

   do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)

          call calc_vlow_blow(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
               gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
               velx(i,j,k),vely(i,j,k),velz(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k), &
               velxlow,velylow,velzlow,Bvecxlow,Bvecylow,Bveczlow, &
               bdotv,b2,dum1,dum2,bxlow,bylow,bzlow)

           bcom_sq(i,j,k) = b2 
           bcom0(i,j,k) = w_lorentz(i,j,k)*bdotv/alp(i,j,k)
           bcomx(i,j,k) = Bvecx(i,j,k)/w_lorentz(i,j,k) + bcom0(i,j,k)*(alp(i,j,k)*velx(i,j,k)-betax(i,j,k))
           bcomy(i,j,k) = Bvecy(i,j,k)/w_lorentz(i,j,k) + bcom0(i,j,k)*(alp(i,j,k)*vely(i,j,k)-betay(i,j,k))
           bcomz(i,j,k) = Bvecz(i,j,k)/w_lorentz(i,j,k) + bcom0(i,j,k)*(alp(i,j,k)*velz(i,j,k)-betaz(i,j,k))

       end do
     end do
   end do
   !$OMP END PARALLEL DO

   return
 
 end subroutine GRHydro_CalcBcom
