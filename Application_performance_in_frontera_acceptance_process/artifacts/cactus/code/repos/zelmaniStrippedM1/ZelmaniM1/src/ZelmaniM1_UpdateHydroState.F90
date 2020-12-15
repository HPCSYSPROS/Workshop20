 /*@@
   @file      ZelmaniM1_UpdateHydroState.F90
   @date      Wed Mar 13 14:18:38 2002
   @author    
   @desc 
   update they hydro variables via a con2prim
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)

!subroutine zm1_UpdateHydroState(CCTK_ARGUMENTS)
!  ! This routine solves for the primitive variables
!  ! after update of tau and scon in radiation-matter coupling.
!  ! The code here is basically a copy of GRHydro's GRHydro_Con2PrimHot.F90.
!
!  implicit none
!  DECLARE_CCTK_ARGUMENTS
!  DECLARE_CCTK_PARAMETERS
!
!  call CCTK_WARN(0,"This subroutine should never be executed.")
!
!end subroutine zm1_UpdateHydroState


subroutine zm1_UpdateEntropy(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  integer :: keytemp
  integer :: keyerr(1)
  integer :: anyerr
  integer :: i,j,k
  real*8  :: dummy(1),rf_precision
  integer :: n,eoskey
  eoskey = 4
  rf_precision = 1.0d-12

  n=1
  keytemp = 1
  keyerr = 0
  anyerr = 0

  ! this is a bit slow, but will be called only when
  ! output is requested, so should be fine
  !$OMP PARALLEL DO PRIVATE(keyerr,anyerr,dummy)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
          
           call EOS_Omni_short(eoskey,keytemp,rf_precision,n,&
                rho(i,j,k),dummy,temperature(i,j,k),y_e(i,j,k),&
                dummy,entropy(i,j,k),dummy,dummy,dummy,dummy,&
                dummy,keyerr,anyerr)

           if(anyerr.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(0,"EOS error in ZelmaniM1 Prim2Con EOS call. Improve error message. Aborting")
              !$OMP END CRITICAL
           endif
           
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO  



end subroutine zm1_UpdateEntropy

#if 1
! this is the old code, which was doing a prim2con after
! a direct update of eps and Y_e, but without momentum update
subroutine zm1_UpdateHydroState(CCTK_ARGUMENTS)

  implicit none

  ! save memory when MP is not used
  ! TARGET as to be before DECLARE_CCTK_ARGUMENTS for gcc 4.1
  TARGET gaa, gab, gac, gbb, gbc, gcc
  TARGET gxx, gxy, gxz, gyy, gyz, gzz
  TARGET lvel, vel
  TARGET gaa_p, gab_p, gac_p, gbb_p, gbc_p, gcc_p
  TARGET gxx_p, gxy_p, gxz_p, gyy_p, gyz_p, gzz_p
  TARGET lvel_p, vel_p
  TARGET gaa_p_p, gab_p_p, gac_p_p, gbb_p_p, gbc_p_p, gcc_p_p
  TARGET gxx_p_p, gxy_p_p, gxz_p_p, gyy_p_p, gyz_p_p, gzz_p_p
  TARGET lvel_p_p, vel_p_p

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  CCTK_REAL :: dummy1, dummy2

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(:,:,:), POINTER :: g11, g12, g13, g22, g23, g33
  CCTK_REAL, DIMENSION(:,:,:,:), POINTER :: vup
  
  !CCTK_REAL, DIMENSION(:,:,:), POINTER :: g11_p, g12_p, g13_p, g22_p, &
  !                                        g23_p, g33_p
  !CCTK_REAL, DIMENSION(:,:,:,:), POINTER :: vup_p
  !CCTK_REAL, DIMENSION(:,:,:), POINTER :: g11_p_p, g12_p_p, g13_p_p, g22_p_p, &
  !                                        g23_p_p, g33_p_p
  !CCTK_REAL, DIMENSION(:,:,:,:), POINTER :: vup_p_p

! begin EOS Omni vars
  integer :: n,keytemp,anyerr
  integer :: eoskey
  integer, allocatable :: keyerr(:,:,:)
  real(8) :: rf_precision
  character(len=256) :: warnline
  n=cctk_lsh(1)*cctk_lsh(2)*cctk_lsh(3);keytemp=0;anyerr=0;
  eoskey = 4
  rf_precision = 1.0d-12
  allocate(keyerr(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)))
! end EOS Omni vars
  
  if (zm1_verbose>0) then 
    call CCTK_INFO('Updating EoS and Conservative Variables in ZelmaniM1.')
  endif
   
  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    g11 => gaa
    g12 => gab
    g13 => gac
    g22 => gbb
    g23 => gbc
    g33 => gcc
    vup => lvel
  else
    g11 => gxx
    g12 => gxy
    g13 => gxz
    g22 => gyy
    g23 => gyz
    g33 => gzz
    vup => vel
  end if

  if(evolve_temper.ne.1) then
     call CCTK_WARN(0,"What's the point of neutrinos if we are not using a nuclear EOS?")
  endif

!$OMP PARALLEL DO PRIVATE(det,i,j,k,anyerr,keyerr, dummy1, dummy2)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
         if (y_e(i,j,k) < GRHydro_Y_e_min) y_e(i,j,k) = GRHydro_Y_e_min      
         if (y_e(i,j,k) > GRHydro_Y_e_max) y_e(i,j,k) = GRHydro_Y_e_max  
      enddo
    enddo
  enddo    
  
  ! first do a grid-function wide call of the EOS
  keytemp = 0
  call EOS_Omni_pressOMP(eoskey,keytemp,rf_precision,n,&
       rho,eps,temperature,y_e,press,keyerr,anyerr)
  
  if(anyerr.ne.0) then
     call CCTK_WARN(0,"EOS error in ZelmaniM1 Prim2Con EOS call. Improve error message. Aborting")
  endif

!$OMP PARALLEL DO PRIVATE(det,i,j,k,anyerr,keyerr, dummy1, dummy2)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
    do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
      do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)
         
         y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)

         call prim2con_np(g11(i,j,k), g12(i,j,k), g13(i,j,k), g22(i,j,k), &
              g23(i,j,k), g33(i,j,k), sdetg(i,j,k), dens(i,j,k), &
              scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), tau(i,j,k),&
              rho(i,j,k), velx(i,j,k), vely(i,j,k), velz(i,j,k), &
              eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))

           !!! NEED TO PERFORM A MOMENTUM UPDATE HERE !!!
        
      end do
    end do
  end do
!$OMP END PARALLEL DO

  deallocate(keyerr)
contains 
   subroutine prim2con_np(gxx, gxy, gxz, gyy, gyz, gzz, sdet, ddens, &
      dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w)
 
   implicit none

   DECLARE_CCTK_PARAMETERS
 
   CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdet
   CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
        deps, dpress, w, vlowx, vlowy, vlowz

   w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
        *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
        *dvely*dvelz))

   vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
   vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
   vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz
   
   ddens = sdet * drho * w
   dsx = sdet * (drho*(1+deps)+dpress)*w*w * vlowx
   dsy = sdet * (drho*(1+deps)+dpress)*w*w * vlowy
   dsz = sdet * (drho*(1+deps)+dpress)*w*w * vlowz
   dtau = sdet * ((drho*(1+deps)+dpress)*w*w - dpress) - ddens

end subroutine prim2con_np



end subroutine zm1_UpdateHydroState
#endif
