 /*@@
   @file      GRHydro_UpdateMaskM.F90
   @date      Sep 2, 2010
   @author    
   @desc 
   Alter the update terms if inside the atmosphere or excision region
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_Macros.h"
#include "SpaceMask.h"

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)
#define velx_p(i,j,k) vup_p(i,j,k,1)
#define vely_p(i,j,k) vup_p(i,j,k,2)
#define velz_p(i,j,k) vup_p(i,j,k,3)
#define velx_p_p(i,j,k) vup_p_p(i,j,k,1)
#define vely_p_p(i,j,k) vup_p_p(i,j,k,2)
#define velz_p_p(i,j,k) vup_p_p(i,j,k,3)

!!$ We don't need to adapt GRHydroUpdateAtmosphereMask, GRHydro_SetupMask, or  
!!$ since we need to evolve Bvec in the atmosphere

!!$ In GRHydro_AtmosphereResetM, we just need to switch the P2C calls to MHD

 /*@@
   @routine    GRHydro_AtmosphereResetM
   @date       Sep 2, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   After MoL has evolved, if a point is supposed to be reset then do so.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_AtmosphereResetM(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
! end EOS Omni vars

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)

! begin EOS Omni vars
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xye=0.0d0;xeps=0.0d0;xtemp=0.0d0
! end EOS Omni vars

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define vel faulty_vel
#define Bvec faulty_Bvec

  if (verbose.eq.1) call CCTK_INFO("Entering AtmosphereReset.")

!$OMP PARALLEL DO PRIVATE(sdet,keytemp,i,j,k,anyerr,keyerr)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
         
        if (atmosphere_mask(i, j, k) .ne. 0) then

          rho(i,j,k) = GRHydro_rho_min
          vup(i,j,k,1) = 0.0d0
          vup(i,j,k,2) = 0.0d0
          vup(i,j,k,3) = 0.0d0
          sdet = sdetg(i,j,k)

          Bprim(i,j,k,:) = 0.0d0
          
          if(evolve_temper.ne.0) then

             ! set the temperature to be relatively low
             temperature(i,j,k) = grhydro_hot_atmo_temp
             y_e(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)
             call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                  i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11(i,j,k),&
                  g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
                  sdet, dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3),&
                  tau(i,j,k),Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
                  rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3),&
                  eps(i,j,k),press(i,j,k),Bprim(i,j,k,1), &
                  Bprim(i,j,k,2), Bprim(i,j,k,3), w_lorentz(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k))
             y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)

          else

            call prim2conpolytypeM(GRHydro_polytrope_handle, &
                 g11(i,j,k), g12(i,j,k), g13(i,j,k), &
                 g22(i,j,k), g23(i,j,k), g33(i,j,k), sdet, &
                 dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                 tau(i,j,k), Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
                 rho(i,j,k), vup(i,j,k,1), vup(i,j,k,2), &
                 vup(i,j,k,3), eps(i,j,k), press(i,j,k), &
                 Bprim(i,j,k,1),Bprim(i,j,k,2),Bprim(i,j,k,3),w_lorentz(i,j,k))
            if (wk_atmosphere .eq. 0) then
              atmosphere_mask(i, j, k) = 0
              atmosphere_mask_real(i, j, k) = 0
            end if

          end if
        endif

      end do
    end do
  end do
!$OMP END PARALLEL DO

!!$  call GRHydro_BoundariesM(CCTK_PASS_FTOF)
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_vel
#undef faulty_Bvec
  
end subroutine GRHydro_AtmosphereResetM

subroutine GRHydro_InitialAtmosphereResetM(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet
  CCTK_REAL :: rho_min

  CCTK_INT :: eos_handle


  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11_p, g12_p, g13_p, g22_p, &
                                          g23_p, g33_p
  pointer (pg11_p,g11_p), (pg12_p,g12_p), (pg13_p,g13_p), (pg22_p,g22_p), &
                                          (pg23_p,g23_p), (pg33_p,g33_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup_p, Bprim_p
  pointer (pvup_p,vup_p), (pBprim_p,Bprim_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11_p_p, g12_p_p, g13_p_p, g22_p_p, &
                                          g23_p_p, g33_p_p
  pointer (pg11_p_p,g11_p_p), (pg12_p_p,g12_p_p), (pg13_p_p,g13_p_p), (pg22_p_p,g22_p_p), &
                                          (pg23_p_p,g23_p_p), (pg33_p_p,g33_p_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup_p_p, Bprim_p_p
  pointer (pvup_p_p,vup_p_p), (pBprim_p_p,Bprim_p_p)

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xye=0.0d0;xeps=0.0d0;xtemp=0.0d0
! end EOS Omni vars

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
  if (timelevels .gt. 1) then
    if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
      pg11_p = loc(gaa_p)
      pg12_p = loc(gab_p)
      pg13_p = loc(gac_p)
      pg22_p = loc(gbb_p)
      pg23_p = loc(gbc_p)
      pg33_p = loc(gcc_p)
      pvup_p = loc(lvel_p)
      pBprim_p = loc(Bvec_p)
    else
      pg11_p = loc(gxx_p)
      pg12_p = loc(gxy_p)
      pg13_p = loc(gxz_p)
      pg22_p = loc(gyy_p)
      pg23_p = loc(gyz_p)
      pg33_p = loc(gzz_p)
      pvup_p = loc(vel_p)
      pBprim_p = loc(Bvec_p)
    end if
  end if 
  if (timelevels .gt. 2) then
    if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
      pg11_p_p = loc(gaa_p_p)
      pg12_p_p = loc(gab_p_p)
      pg13_p_p = loc(gac_p_p)
      pg22_p_p = loc(gbb_p_p)
      pg23_p_p = loc(gbc_p_p)
      pg33_p_p = loc(gcc_p_p)
      pvup_p_p = loc(lvel_p_p)
      pBprim_p_p = loc(lBvec_p_p)
    else
      pg11_p_p = loc(gxx_p_p)
      pg12_p_p = loc(gxy_p_p)
      pg13_p_p = loc(gxz_p_p)
      pg22_p_p = loc(gyy_p_p)
      pg23_p_p = loc(gyz_p_p)
      pg33_p_p = loc(gzz_p_p)
      pvup_p_p = loc(vel_p_p)
      pBprim_p_p = loc(Bvec_p_p)
    end if
  end if 
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define vel faulty_vel
#define Bvec faulty_Bvec
#define gxx_p faulty_gxx_p
#define gxy_p faulty_gxy_p
#define gxz_p faulty_gxz_p
#define gyy_p faulty_gyy_p
#define gyz_p faulty_gyz_p
#define gzz_p faulty_gzz_p
#define vel_p faulty_vel_p
#define Bvec_p faulty_Bvec_p
#define gxx_p_p faulty_gxx_p_p
#define gxy_p_p faulty_gxy_p_p
#define gxz_p_p faulty_gxz_p_p
#define gyy_p_p faulty_gyy_p_p
#define gyz_p_p faulty_gyz_p_p
#define gzz_p_p faulty_gzz_p_p
#define vel_p_p faulty_vel_p_p
#define Bvec_p_p faulty_Bvec_p_p


  eos_handle = GRHydro_polytrope_handle
  
  rho_min = GRHydro_rho_min
  if (initial_atmosphere_factor .gt. 0) then
    rho_min = rho_min * initial_atmosphere_factor
  endif
  
!$OMP PARALLEL DO PRIVATE(sdet,keytemp,i,j,k,anyerr,keyerr)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        if (rho(i,j,k) .le. rho_min .or. &
            (GRHydro_enable_internal_excision .ne. 0 .and. &
             hydro_excision_mask(i,j,k) .ne. 0) )  then
          rho(i,j,k) = rho_min
          vup(i,j,k,1) = 0.0d0
          vup(i,j,k,2) = 0.0d0
          vup(i,j,k,3) = 0.0d0

          sdet = sdetg(i,j,k)

          if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
             temperature(i,j,k) = grhydro_hot_atmo_temp
             y_e(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)
             call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                  i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11(i,j,k),&
                  g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
                  sdet, dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3),&
                  tau(i,j,k),Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
                  rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3),&
                  eps(i,j,k),press(i,j,k),Bprim(i,j,k,1), &
                  Bprim(i,j,k,2), Bprim(i,j,k,3), w_lorentz(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k))
             y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)

          else

          keytemp = 0
          call EOS_Omni_press(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                 rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
          call EOS_Omni_EpsFromPress(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                 rho(i,j,k),xeps,xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)

          call prim2conpolytypeM(eos_handle, &
               g11(i,j,k), g12(i,j,k), g13(i,j,k), &
               g22(i,j,k), g23(i,j,k), g33(i,j,k), sdet, &
               dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
               tau(i,j,k), Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
               rho(i,j,k), vup(i,j,k,1), vup(i,j,k,2), &
               vup(i,j,k,3), eps(i,j,k), press(i,j,k), &
               Bprim(i,j,k,1),Bprim(i,j,k,2),Bprim(i,j,k,3),w_lorentz(i,j,k))
          end if
        end if
        if (timelevels .gt. 1) then
          if (rho_p(i,j,k) .le. rho_min) then
            rho_p(i,j,k) = rho_min
            vup_p(i,j,k,1) = 0.0d0
            vup_p(i,j,k,2) = 0.0d0
            vup_p(i,j,k,3) = 0.0d0

            sdet = sqrt(SPATIAL_DETERMINANT(g11_p(i,j,k), g12_p(i,j,k), g13_p(i,j,k), \
                 g22_p(i,j,k), g23_p(i,j,k), g33_p(i,j,k)))

            if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
               temperature_p(i,j,k) = grhydro_hot_atmo_temp
               y_e_p(i,j,k) = grhydro_hot_atmo_Y_e
               keytemp = 1
               call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p(i,j,k),eps_p(i,j,k),temperature_p(i,j,k),y_e_p(i,j,k),&
                    press_p(i,j,k),keyerr,anyerr)
               call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                    i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11_p(i,j,k),&
                    g12_p(i,j,k),g13_p(i,j,k),g22_p(i,j,k),g23_p(i,j,k),g33_p(i,j,k),&
                    sdet, dens_p(i,j,k),scon_p(i,j,k,1),scon_p(i,j,k,2),scon_p(i,j,k,3),&
                    tau_p(i,j,k),Bcons_p(i,j,k,1),Bcons_p(i,j,k,2),Bcons_p(i,j,k,3),&
                    rho_p(i,j,k),vup_p(i,j,k,1),vup_p(i,j,k,2),vup_p(i,j,k,3),&
                    eps_p(i,j,k),press_p(i,j,k),Bprim_p(i,j,k,1), &
                    Bprim_p(i,j,k,2), Bprim_p(i,j,k,3), w_lorentz_p(i,j,k),&
                    temperature_p(i,j,k),y_e_p(i,j,k))
               y_e_con_p(i,j,k) = dens_p(i,j,k) * y_e_p(i,j,k)

            else

            keytemp = 0
            call EOS_Omni_press(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p(i,j,k),eps_p(i,j,k),xtemp,xye,press_p(i,j,k),keyerr,anyerr)
            call EOS_Omni_EpsFromPress(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p(i,j,k),xeps,xtemp,xye,press_p(i,j,k),eps_p(i,j,k),keyerr,anyerr)

            call prim2conpolytypeM(eos_handle, &
                 g11_p(i,j,k), g12_p(i,j,k), g13_p(i,j,k), &
                 g22_p(i,j,k), g23_p(i,j,k), g33_p(i,j,k), sdet, &
                 dens_p(i,j,k), scon_p(i,j,k,1), scon_p(i,j,k,2), scon_p(i,j,k,3), &
                 tau_p(i,j,k), Bcons_p(i,j,k,1),Bcons_p(i,j,k,2),Bcons_p(i,j,k,3),&
                 rho_p(i,j,k), vup_p(i,j,k,1), vup_p(i,j,k,2), &
                 vup_p(i,j,k,3), eps_p(i,j,k), press_p(i,j,k), &
                 Bprim_p(i,j,k,1),Bprim_p(i,j,k,2),Bprim_p(i,j,k,3),w_lorentz_p(i,j,k))
            end if
          end if
        end if

        if (timelevels .gt. 2) then
          if (rho_p_p(i,j,k) .le. rho_min) then
            rho_p_p(i,j,k) = rho_min
            vup_p_p(i,j,k,1) = 0.0d0
            vup_p_p(i,j,k,2) = 0.0d0
            vup_p_p(i,j,k,3) = 0.0d0

            sdet = sqrt(SPATIAL_DETERMINANT(g11_p_p(i,j,k), g12_p_p(i,j,k), g13_p_p(i,j,k), \
                 g22_p_p(i,j,k), g23_p_p(i,j,k), g33_p_p(i,j,k)))

            if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
               temperature_p_p(i,j,k) = grhydro_hot_atmo_temp
               y_e_p_p(i,j,k) = grhydro_hot_atmo_Y_e
               keytemp = 1
               call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p_p(i,j,k),eps_p_p(i,j,k),temperature_p_p(i,j,k),y_e_p_p(i,j,k),&
                    press_p_p(i,j,k),keyerr,anyerr)
               call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                    i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11_p_p(i,j,k),&
                    g12_p_p(i,j,k),g13_p_p(i,j,k),g22_p_p(i,j,k),g23_p_p(i,j,k),g33_p_p(i,j,k),&
                    sdet, dens_p_p(i,j,k),scon_p_p(i,j,k,1),scon_p_p(i,j,k,2),scon_p_p(i,j,k,3),&
                    tau_p_p(i,j,k),Bcons_p_p(i,j,k,1),Bcons_p_p(i,j,k,2),Bcons_p_p(i,j,k,3),&
                    rho_p_p(i,j,k),vup_p_p(i,j,k,1),vup_p_p(i,j,k,2),vup_p_p(i,j,k,3),&
                    eps_p_p(i,j,k),press_p_p(i,j,k),Bprim_p_p(i,j,k,1), &
                    Bprim_p_p(i,j,k,2), Bprim_p_p(i,j,k,3), w_lorentz_p_p(i,j,k),&
                    temperature_p_p(i,j,k),y_e_p_p(i,j,k))
               y_e_con_p_p(i,j,k) = dens_p_p(i,j,k) * y_e_p_p(i,j,k)

            else

            keytemp = 0
            call EOS_Omni_press(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p_p(i,j,k),eps_p_p(i,j,k),xtemp,xye,press_p_p(i,j,k),keyerr,anyerr)
            call EOS_Omni_EpsFromPress(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p_p(i,j,k),xeps,xtemp,xye,press_p_p(i,j,k),eps_p_p(i,j,k),keyerr,anyerr)

            call prim2conpolytypeM(eos_handle, &
                 g11_p_p(i,j,k), g12_p_p(i,j,k), g13_p_p(i,j,k), &
                 g22_p_p(i,j,k), g23_p_p(i,j,k), g33_p_p(i,j,k), sdet, &
                 dens_p_p(i,j,k), scon_p_p(i,j,k,1), scon_p_p(i,j,k,2), scon_p_p(i,j,k,3), &
                 tau_p_p(i,j,k), Bcons_p_p(i,j,k,1),Bcons_p_p(i,j,k,2),Bcons_p_p(i,j,k,3),&
                 rho_p_p(i,j,k), vup_p_p(i,j,k,1), vup_p_p(i,j,k,2), &
                 vup_p_p(i,j,k,3), eps_p_p(i,j,k), press_p_p(i,j,k), &
                 Bprim_p_p(i,j,k,1),Bprim_p_p(i,j,k,2),Bprim_p_p(i,j,k,3),w_lorentz_p_p(i,j,k))
            end if
          end if
        end if

      end do
    end do
  end do
!$OMP END PARALLEL DO

  write(*,*) "     GRHydro_InitialAtmosphereReset"
!!$  call GRHydro_BoundariesM(CCTK_PASS_FTOF)
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_vel
#undef faulty_Bvec
#undef faulty_gxx_p
#undef faulty_gxy_p
#undef faulty_gxz_p
#undef faulty_gyy_p
#undef faulty_gyz_p
#undef faulty_gzz_p
#undef faulty_vel_p
#undef faulty_Bvec_p
#undef faulty_gxx_p_p
#undef faulty_gxy_p_p
#undef faulty_gxz_p_p
#undef faulty_gyy_p_p
#undef faulty_gyz_p_p
#undef faulty_gzz_p_p
#undef faulty_vel_p_p
#undef faulty_Bvec_p_p

end subroutine GRHydro_InitialAtmosphereResetM


 /*@@
   @routine    GRHydro_AtmosphereResetAM
   @date       Sep 2, 2010
   @author     Tanja Bode, Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   After MoL has evolved, if a point is supposed to be reset then do so.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_AtmosphereResetAM(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
! end EOS Omni vars

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)

! begin EOS Omni vars
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xye=0.0d0;xeps=0.0d0;xtemp=0.0d0
! end EOS Omni vars

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define vel faulty_vel
#define Bvec faulty_Bvec

  if (verbose.eq.1) call CCTK_INFO("Entering AtmosphereReset.")

!$OMP PARALLEL DO PRIVATE(sdet,keytemp,i,j,k,anyerr,keyerr)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
         
        if (atmosphere_mask(i, j, k) .ne. 0) then

          rho(i,j,k) = GRHydro_rho_min
          vup(i,j,k,1) = 0.0d0
          vup(i,j,k,2) = 0.0d0
          vup(i,j,k,3) = 0.0d0
          sdet = sdetg(i,j,k)

          if(evolve_temper.ne.0) then

             ! set the temperature to be relatively low
             temperature(i,j,k) = grhydro_hot_atmo_temp
             y_e(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)
             call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                  i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11(i,j,k),&
                  g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
                  sdet, dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3),&
                  tau(i,j,k),Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
                  rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3),&
                  eps(i,j,k),press(i,j,k),Bprim(i,j,k,1), &
                  Bprim(i,j,k,2), Bprim(i,j,k,3), w_lorentz(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k))
             y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)

          else

            call prim2conpolytypeM(GRHydro_polytrope_handle, &
                 g11(i,j,k), g12(i,j,k), g13(i,j,k), &
                 g22(i,j,k), g23(i,j,k), g33(i,j,k), sdet, &
                 dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                 tau(i,j,k), Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
                 rho(i,j,k), vup(i,j,k,1), vup(i,j,k,2), &
                 vup(i,j,k,3), eps(i,j,k), press(i,j,k), &
                 Bprim(i,j,k,1),Bprim(i,j,k,2),Bprim(i,j,k,3),w_lorentz(i,j,k))
            if (wk_atmosphere .eq. 0) then
              atmosphere_mask(i, j, k) = 0
              atmosphere_mask_real(i, j, k) = 0
            end if

          end if
        endif

      end do
    end do
  end do
!$OMP END PARALLEL DO

!!$  call GRHydro_BoundariesM(CCTK_PASS_FTOF)
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_vel
#undef faulty_Bvec
  
end subroutine GRHydro_AtmosphereResetAM

subroutine GRHydro_InitialAtmosphereResetAM(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet
  CCTK_REAL :: rho_min

  CCTK_INT :: eos_handle


  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11_p, g12_p, g13_p, g22_p, &
                                          g23_p, g33_p
  pointer (pg11_p,g11_p), (pg12_p,g12_p), (pg13_p,g13_p), (pg22_p,g22_p), &
                                          (pg23_p,g23_p), (pg33_p,g33_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup_p, Bprim_p
  pointer (pvup_p,vup_p), (pBprim_p,Bprim_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11_p_p, g12_p_p, g13_p_p, g22_p_p, &
                                          g23_p_p, g33_p_p
  pointer (pg11_p_p,g11_p_p), (pg12_p_p,g12_p_p), (pg13_p_p,g13_p_p), (pg22_p_p,g22_p_p), &
                                          (pg23_p_p,g23_p_p), (pg33_p_p,g33_p_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup_p_p, Bprim_p_p
  pointer (pvup_p_p,vup_p_p), (pBprim_p_p,Bprim_p_p)

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xye=0.0d0;xeps=0.0d0;xtemp=0.0d0
! end EOS Omni vars

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
  if (timelevels .gt. 1) then
    if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
      pg11_p = loc(gaa_p)
      pg12_p = loc(gab_p)
      pg13_p = loc(gac_p)
      pg22_p = loc(gbb_p)
      pg23_p = loc(gbc_p)
      pg33_p = loc(gcc_p)
      pvup_p = loc(lvel_p)
      pBprim_p = loc(Bvec_p)
    else
      pg11_p = loc(gxx_p)
      pg12_p = loc(gxy_p)
      pg13_p = loc(gxz_p)
      pg22_p = loc(gyy_p)
      pg23_p = loc(gyz_p)
      pg33_p = loc(gzz_p)
      pvup_p = loc(vel_p)
      pBprim_p = loc(Bvec_p)
    end if
  end if 
  if (timelevels .gt. 2) then
    if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
      pg11_p_p = loc(gaa_p_p)
      pg12_p_p = loc(gab_p_p)
      pg13_p_p = loc(gac_p_p)
      pg22_p_p = loc(gbb_p_p)
      pg23_p_p = loc(gbc_p_p)
      pg33_p_p = loc(gcc_p_p)
      pvup_p_p = loc(lvel_p_p)
      pBprim_p_p = loc(lBvec_p_p)
    else
      pg11_p_p = loc(gxx_p_p)
      pg12_p_p = loc(gxy_p_p)
      pg13_p_p = loc(gxz_p_p)
      pg22_p_p = loc(gyy_p_p)
      pg23_p_p = loc(gyz_p_p)
      pg33_p_p = loc(gzz_p_p)
      pvup_p_p = loc(vel_p_p)
      pBprim_p_p = loc(Bvec_p_p)
    end if
  end if 
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define vel faulty_vel
#define Bvec faulty_Bvec
#define gxx_p faulty_gxx_p
#define gxy_p faulty_gxy_p
#define gxz_p faulty_gxz_p
#define gyy_p faulty_gyy_p
#define gyz_p faulty_gyz_p
#define gzz_p faulty_gzz_p
#define vel_p faulty_vel_p
#define Bvec_p faulty_Bvec_p
#define gxx_p_p faulty_gxx_p_p
#define gxy_p_p faulty_gxy_p_p
#define gxz_p_p faulty_gxz_p_p
#define gyy_p_p faulty_gyy_p_p
#define gyz_p_p faulty_gyz_p_p
#define gzz_p_p faulty_gzz_p_p
#define vel_p_p faulty_vel_p_p
#define Bvec_p_p faulty_Bvec_p_p


  eos_handle = GRHydro_polytrope_handle
  
  rho_min = GRHydro_rho_min
  if (initial_atmosphere_factor .gt. 0) then
    rho_min = rho_min * initial_atmosphere_factor
  endif
  
!$OMP PARALLEL DO PRIVATE(sdet,keytemp,i,j,k,anyerr,keyerr)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        if (rho(i,j,k) .le. rho_min .or. &
            (GRHydro_enable_internal_excision .ne. 0 .and. &
             hydro_excision_mask(i,j,k) .ne. 0) )  then
          rho(i,j,k) = rho_min
          vup(i,j,k,1) = 0.0d0
          vup(i,j,k,2) = 0.0d0
          vup(i,j,k,3) = 0.0d0

          sdet = sdetg(i,j,k)

          if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
             temperature(i,j,k) = grhydro_hot_atmo_temp
             y_e(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)
             call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                  i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11(i,j,k),&
                  g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
                  sdet, dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3),&
                  tau(i,j,k),Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
                  rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3),&
                  eps(i,j,k),press(i,j,k),Bprim(i,j,k,1), &
                  Bprim(i,j,k,2), Bprim(i,j,k,3), w_lorentz(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k))
             y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)

          else

          keytemp = 0
          call EOS_Omni_press(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                 rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
          call EOS_Omni_EpsFromPress(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                 rho(i,j,k),xeps,xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)

          call prim2conpolytypeM(eos_handle, &
               g11(i,j,k), g12(i,j,k), g13(i,j,k), &
               g22(i,j,k), g23(i,j,k), g33(i,j,k), sdet, &
               dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
               tau(i,j,k), Bcons(i,j,k,1),Bcons(i,j,k,2),Bcons(i,j,k,3),&
               rho(i,j,k), vup(i,j,k,1), vup(i,j,k,2), &
               vup(i,j,k,3), eps(i,j,k), press(i,j,k), &
               Bprim(i,j,k,1),Bprim(i,j,k,2),Bprim(i,j,k,3),w_lorentz(i,j,k))
          end if
        end if
        if (timelevels .gt. 1) then
          if (rho_p(i,j,k) .le. rho_min) then
            rho_p(i,j,k) = rho_min
            vup_p(i,j,k,1) = 0.0d0
            vup_p(i,j,k,2) = 0.0d0
            vup_p(i,j,k,3) = 0.0d0

            sdet = sqrt(SPATIAL_DETERMINANT(g11_p(i,j,k), g12_p(i,j,k), g13_p(i,j,k), \
                 g22_p(i,j,k), g23_p(i,j,k), g33_p(i,j,k)))

            if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
               temperature_p(i,j,k) = grhydro_hot_atmo_temp
               y_e_p(i,j,k) = grhydro_hot_atmo_Y_e
               keytemp = 1
               call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p(i,j,k),eps_p(i,j,k),temperature_p(i,j,k),y_e_p(i,j,k),&
                    press_p(i,j,k),keyerr,anyerr)
               call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                    i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11_p(i,j,k),&
                    g12_p(i,j,k),g13_p(i,j,k),g22_p(i,j,k),g23_p(i,j,k),g33_p(i,j,k),&
                    sdet, dens_p(i,j,k),scon_p(i,j,k,1),scon_p(i,j,k,2),scon_p(i,j,k,3),&
                    tau_p(i,j,k),Bcons_p(i,j,k,1),Bcons_p(i,j,k,2),Bcons_p(i,j,k,3),&
                    rho_p(i,j,k),vup_p(i,j,k,1),vup_p(i,j,k,2),vup_p(i,j,k,3),&
                    eps_p(i,j,k),press_p(i,j,k),Bprim_p(i,j,k,1), &
                    Bprim_p(i,j,k,2), Bprim_p(i,j,k,3), w_lorentz_p(i,j,k),&
                    temperature_p(i,j,k),y_e_p(i,j,k))
               y_e_con_p(i,j,k) = dens_p(i,j,k) * y_e_p(i,j,k)

            else

            keytemp = 0
            call EOS_Omni_press(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p(i,j,k),eps_p(i,j,k),xtemp,xye,press_p(i,j,k),keyerr,anyerr)
            call EOS_Omni_EpsFromPress(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p(i,j,k),xeps,xtemp,xye,press_p(i,j,k),eps_p(i,j,k),keyerr,anyerr)

            call prim2conpolytypeM(eos_handle, &
                 g11_p(i,j,k), g12_p(i,j,k), g13_p(i,j,k), &
                 g22_p(i,j,k), g23_p(i,j,k), g33_p(i,j,k), sdet, &
                 dens_p(i,j,k), scon_p(i,j,k,1), scon_p(i,j,k,2), scon_p(i,j,k,3), &
                 tau_p(i,j,k), Bcons_p(i,j,k,1),Bcons_p(i,j,k,2),Bcons_p(i,j,k,3),&
                 rho_p(i,j,k), vup_p(i,j,k,1), vup_p(i,j,k,2), &
                 vup_p(i,j,k,3), eps_p(i,j,k), press_p(i,j,k), &
                 Bprim_p(i,j,k,1),Bprim_p(i,j,k,2),Bprim_p(i,j,k,3),w_lorentz_p(i,j,k))
            end if
          end if
        end if

        if (timelevels .gt. 2) then
          if (rho_p_p(i,j,k) .le. rho_min) then
            rho_p_p(i,j,k) = rho_min
            vup_p_p(i,j,k,1) = 0.0d0
            vup_p_p(i,j,k,2) = 0.0d0
            vup_p_p(i,j,k,3) = 0.0d0

            sdet = sqrt(SPATIAL_DETERMINANT(g11_p_p(i,j,k), g12_p_p(i,j,k), g13_p_p(i,j,k), \
                 g22_p_p(i,j,k), g23_p_p(i,j,k), g33_p_p(i,j,k)))

            if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
               temperature_p_p(i,j,k) = grhydro_hot_atmo_temp
               y_e_p_p(i,j,k) = grhydro_hot_atmo_Y_e
               keytemp = 1
               call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p_p(i,j,k),eps_p_p(i,j,k),temperature_p_p(i,j,k),y_e_p_p(i,j,k),&
                    press_p_p(i,j,k),keyerr,anyerr)
               call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                    i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),g11_p_p(i,j,k),&
                    g12_p_p(i,j,k),g13_p_p(i,j,k),g22_p_p(i,j,k),g23_p_p(i,j,k),g33_p_p(i,j,k),&
                    sdet, dens_p_p(i,j,k),scon_p_p(i,j,k,1),scon_p_p(i,j,k,2),scon_p_p(i,j,k,3),&
                    tau_p_p(i,j,k),Bcons_p_p(i,j,k,1),Bcons_p_p(i,j,k,2),Bcons_p_p(i,j,k,3),&
                    rho_p_p(i,j,k),vup_p_p(i,j,k,1),vup_p_p(i,j,k,2),vup_p_p(i,j,k,3),&
                    eps_p_p(i,j,k),press_p_p(i,j,k),Bprim_p_p(i,j,k,1), &
                    Bprim_p_p(i,j,k,2), Bprim_p_p(i,j,k,3), w_lorentz_p_p(i,j,k),&
                    temperature_p_p(i,j,k),y_e_p_p(i,j,k))
               y_e_con_p_p(i,j,k) = dens_p_p(i,j,k) * y_e_p_p(i,j,k)

            else

            keytemp = 0
            call EOS_Omni_press(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p_p(i,j,k),eps_p_p(i,j,k),xtemp,xye,press_p_p(i,j,k),keyerr,anyerr)
            call EOS_Omni_EpsFromPress(eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho_p_p(i,j,k),xeps,xtemp,xye,press_p_p(i,j,k),eps_p_p(i,j,k),keyerr,anyerr)

            call prim2conpolytypeM(eos_handle, &
                 g11_p_p(i,j,k), g12_p_p(i,j,k), g13_p_p(i,j,k), &
                 g22_p_p(i,j,k), g23_p_p(i,j,k), g33_p_p(i,j,k), sdet, &
                 dens_p_p(i,j,k), scon_p_p(i,j,k,1), scon_p_p(i,j,k,2), scon_p_p(i,j,k,3), &
                 tau_p_p(i,j,k), Bcons_p_p(i,j,k,1),Bcons_p_p(i,j,k,2),Bcons_p_p(i,j,k,3),&
                 rho_p_p(i,j,k), vup_p_p(i,j,k,1), vup_p_p(i,j,k,2), &
                 vup_p_p(i,j,k,3), eps_p_p(i,j,k), press_p_p(i,j,k), &
                 Bprim_p_p(i,j,k,1),Bprim_p_p(i,j,k,2),Bprim_p_p(i,j,k,3),w_lorentz_p_p(i,j,k))
            end if
          end if
        end if

      end do
    end do
  end do
!$OMP END PARALLEL DO

  write(*,*) "     GRHydro_InitialAtmosphereReset"
!!$  call GRHydro_BoundariesM(CCTK_PASS_FTOF)
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_vel
#undef faulty_Bvec
#undef faulty_gxx_p
#undef faulty_gxy_p
#undef faulty_gxz_p
#undef faulty_gyy_p
#undef faulty_gyz_p
#undef faulty_gzz_p
#undef faulty_vel_p
#undef faulty_Bvec_p
#undef faulty_gxx_p_p
#undef faulty_gxy_p_p
#undef faulty_gxz_p_p
#undef faulty_gyy_p_p
#undef faulty_gyz_p_p
#undef faulty_gzz_p_p
#undef faulty_vel_p_p
#undef faulty_Bvec_p_p

end subroutine GRHydro_InitialAtmosphereResetAM

