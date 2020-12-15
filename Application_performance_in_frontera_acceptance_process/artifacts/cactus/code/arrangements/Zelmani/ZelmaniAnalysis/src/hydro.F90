#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine ZelmaniAnalysis_HydroLocalF(CCTK_ARGUMENTS)

  ! use module matdet from TAT/TGRTensor
  use matdet
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  ! local stuff
  integer i,j,k
  logical err
  integer nx, ny, nz
  integer ghostx,ghosty,ghostz
  integer aa,bb,cc,dd

  ! some helpers
  CCTK_REAL, parameter :: pi = 3.14159265358979d0
  CCTK_REAL, parameter :: invpi16 = 1.0d0 / (16.0d0 * pi)
  CCTK_REAL, parameter :: invpi24 = 1.0d0 / (24.0d0 * pi)
  CCTK_REAL dx,dy,dz
  CCTK_REAL w2,sg,g512,detg,rho0,rhoh
  CCTK_REAL gg(3,3),gu(3,3),at(3,3),atu(3,3),atul(3,3)
  CCTK_REAL at2
  
  CCTK_REAL Jx, Jy, Omega, cyl_rad_sq
  CCTK_REAL invalp
  CCTK_REAL ulx, uly, ulz, gtx, gty, gtz
  

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  ghostx = cctk_nghostzones(1)
  ghosty = cctk_nghostzones(2)
  ghostz = cctk_nghostzones(3)
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  do k=ghostz,nz-ghostz
     do j=ghosty,ny-ghosty
        do i=ghostx,nx-ghostx

           ! some definitions
           gg(1,1) = gxx(i,j,k)
           gg(1,2) = gxy(i,j,k)
           gg(1,3) = gxz(i,j,k)
           gg(2,1) = gg(1,2)
           gg(2,2) = gyy(i,j,k)
           gg(2,3) = gxz(i,j,k)
           gg(3,1) = gg(1,3)
           gg(3,2) = gg(2,3)
           gg(3,3) = gzz(i,j,k)

           ! the simple stuff first & helpers
           call calc_symdet3(gg,detg,err)
           g512 = detg**(5.0d0/12.0d0)
           sg = sqrt(detg)
           rhoh = rho(i,j,k) + rho(i,j,k)*eps(i,j,k) + press(i,j,k)
           w2 = w_lorentz(i,j,k) * w_lorentz(i,j,k)
           rho0 = rhoh * w2 - press(i,j,k)
           invalp = 1.0d0/alp(i,j,k)
           
           ! final integrand for M_p
           proper_mass_local(i,j,k) = dens(i,j,k) * (1.0d0 + eps(i,j,k))

           ! final integrand for E_kin
           total_kinetic_energy_local(i,j,k) = (w2 - 1.0d0)  &
                * sg * rho(i,j,k)

           ! final integrand for M_grav
           gravitational_mass_local(i,j,k) = rho0 * g512

           ! final integrand for M_bary and M_grav inside 10^12 g/ccm
           if(rho(i,j,k).gt.1.61930347d-06) then
              gravitational_mass_local_in1e12(i,j,k) = rho0 * g512
              baryonic_mass_local_in1e12(i,j,k)      = dens(i,j,k)
           else
              gravitational_mass_local_in1e12(i,j,k) = 0.0d0
              baryonic_mass_local_in1e12(i,j,k)      = 0.0d0
           endif

           ! angular momentum and rotational energy (matter based)
           cyl_rad_sq = x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k)
           if (cyl_rad_sq > 1.0d-10) then
              Omega = (x(i,j,k)*(alp(i,j,k)*vel(i,j,k,2)-betay(i,j,k)) - &
		   y(i,j,k)*(alp(i,j,k)*vel(i,j,k,1)-betax(i,j,k))) / &
		   cyl_rad_sq
           else
              Omega = 0.0d0
           endif

           gtx = gxx(i,j,k) * betax(i,j,k) + &
                gxy(i,j,k) * betay(i,j,k)  + &
                gxz(i,j,k) * betaz(i,j,k)

           gty = gxy(i,j,k) * betax(i,j,k) + &
                gyy(i,j,k) * betay(i,j,k)  + &
                gyz(i,j,k) * betaz(i,j,k)    

           gtz = gxz(i,j,k) * betax(i,j,k) + &
                gyz(i,j,k) * betay(i,j,k)  + &
                gzz(i,j,k) * betaz(i,j,k)

           ulx = w_lorentz(i,j,k) *  &
                (gtx * invalp +      &
	         gxx(i,j,k) * (vel(i,j,k,1) - betax(i,j,k) * invalp) + &
                 gxy(i,j,k) * (vel(i,j,k,2) - betay(i,j,k) * invalp) + &
                 gxz(i,j,k) * (vel(i,j,k,3) - betaz(i,j,k) * invalp))


           uly = w_lorentz(i,j,k) *  &
                (gty * invalp +  &
                gxy(i,j,k) * (vel(i,j,k,1) - betax(i,j,k) * invalp) + &
                gyy(i,j,k) * (vel(i,j,k,2) - betay(i,j,k) * invalp) + &
                gyz(i,j,k) * (vel(i,j,k,3) - betaz(i,j,k) * invalp))
    
           ulz = w_lorentz(i,j,k) *  &
                (gtz * invalp +  &
                gxz(i,j,k) * (vel(i,j,k,1) - betax(i,j,k) * invalp) + &
                gyz(i,j,k) * (vel(i,j,k,2) - betay(i,j,k) * invalp) + &
                gzz(i,j,k) * (vel(i,j,k,3) - betaz(i,j,k) * invalp))


           Jx = rhoh * w_lorentz(i,j,k) * ulx;
           Jy = rhoh * w_lorentz(i,j,k) * uly;

           ! final integrand for Jz
           angular_momentum_local(i,j,k) = sg * &
                 (x(i,j,k)*Jy - y(i,j,k)*Jx)

           ! final integrand for T
           kinetic_energy_local(i,j,k) = 0.5d0 * Omega * &
                angular_momentum_local(i,j,k)
           
           if(do_adm_full_hydro/=0) then
           ! Gear up to computing the ADM Mass and
           ! angular momentum a la Ho et al. 2002

              ! upper metric
              call ZA_UpperMetric(gg,detg,gu)
              ! A tilde and its upper components
              ! A tilde comes from ML and has
              ! covariant components
              at(1,1) = at11(i,j,k)
              at(1,2) = at12(i,j,k)
              at(1,3) = at13(i,j,k)
              at(2,1) = at12(i,j,k)
              at(2,2) = at22(i,j,k)
              at(2,3) = at23(i,j,k)
              at(3,1) = at13(i,j,k)
              at(3,2) = at23(i,j,k)
              at(3,3) = at33(i,j,k)

              ! first compute atul which is Atilde^j_k
              do aa=1,3
                 do bb=1,3
                    atul(aa,bb) = 0.0d0
                    do cc=1,3
                       atul(aa,bb) = atul(aa,bb) + gu(aa,cc)*at(cc,bb)
                    enddo
                 enddo
              enddo

              ! now the full thing Atilde^jk
              do aa=1,3
                 do bb=1,3
                    atu(aa,bb) = 0.0d0
                    do cc=1,3
                       atu(aa,bb) = atu(aa,bb) + gu(aa,cc)*atul(aa,cc)
                    enddo
                 enddo
              enddo
              
              ! Atilde_ij Atilde^ij
              at2 = 0.0d0
              do aa=1,3
                 do bb=1,3
                    at2 = at2 + at(aa,bb)*atu(aa,bb)
                 enddo
              enddo

           endif

        enddo
     enddo
  enddo


end subroutine ZelmaniAnalysis_HydroLocalF

subroutine ZA_UpperMetric(gg,detg,gu)

  implicit none
  CCTK_REAL gg(3,3),gu(3,3),detg
  CCTK_REAL invdetg 

  invdetg = 1.0d0/detg

  gu(1,1) = (-gg(2,3)*gg(2,3) + gg(2,2)*gg(3,3))*invdetg
  gu(1,2) = ( gg(1,3)*gg(2,3) - gg(1,2)*gg(3,3))*invdetg
  gu(1,3) = (-gg(1,3)*gg(2,2) + gg(1,2)*gg(2,3))*invdetg
  gu(2,1) =   gu(1,2)
  gu(2,2) = (-gg(1,3)*gg(1,3) + gg(1,1)*gg(3,3))*invdetg
  gu(2,3) = ( gg(1,2)*gg(1,3) - gg(1,1)*gg(2,3))*invdetg
  gu(3,1) = gu(1,3)
  gu(3,2) = gu(2,3)
  gu(3,3) = (-gg(1,2)*gg(1,2) + gg(1,1)*gg(2,2))*invdetg

end subroutine ZA_UpperMetric
