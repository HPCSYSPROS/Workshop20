 /*@@
   @file      GRHydro_UpdateMask.F90
   @date      Wed Mar 13 14:18:38 2002
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

subroutine GRHydroUpdateAtmosphereMask(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_REAL :: frac
  logical :: evolve_Bvec
  logical :: evolve_Avec

  frac = CCTK_DELTA_TIME

  evolve_Bvec = .false.
  evolve_Avec = .false.

  if (CCTK_EQUALS(Bvec_evolution_method,"GRHydro")) then
     evolve_Bvec = .true.
  elseif  (CCTK_EQUALS(Bvec_evolution_method,"GRHydro_Avec") .or. CCTK_EQUALS(Bvec_evolution_method,"GRHydro_Avec_centered")) then
     evolve_Avec = .true.
  endif

  if(evolve_temper.ne.1.and.evolve_Y_e.ne.1) then
     !$OMP PARALLEL DO PRIVATE(k,j,i)
     do k = 1+GRHydro_stencil, cctk_lsh(3)-GRHydro_stencil
        do j = 1+GRHydro_stencil, cctk_lsh(2)-GRHydro_stencil
           do i = 1+GRHydro_stencil, cctk_lsh(1)-GRHydro_stencil
              if ( GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0) .or. &
                   (atmosphere_mask(i,j,k) .ne. 0) .or. &
                      (tau(i,j,k) + frac * taurhs(i,j,k) .le. 0.d0) .or. &
                      (dens(i,j,k) + frac * densrhs(i,j,k) .le. 0.d0) ) then
                 densrhs(i,j,k) = 0.0d0
                 srhs(i,j,k,:)   = 0.0d0
                 taurhs(i,j,k)  = 0.0d0
                 
                 if (evolve_Bvec) then
                    Bconsrhs(i,j,k,:) = 0.0d0
                 endif
                 if (evolve_Avec) then
                    Avecrhs(i,j,k,:) = 0.0d0
                 endif
                 
                 ! TODO: Need to set Avecrhs and Aphirhs to zero for centered / or staggered case if vector potential is evolved

                 ! Set real-valued mask! This will be sync'ed and right after syncing translated to
                 ! our standard integer based mask (so that atmosphere_mask is still valid!).
                 atmosphere_mask_real(i,j,k) = 1
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     ! allow negative total energy
     !$OMP PARALLEL DO PRIVATE(k,j,i)
     do k = 1+GRHydro_stencil, cctk_lsh(3)-GRHydro_stencil
        do j = 1+GRHydro_stencil, cctk_lsh(2)-GRHydro_stencil
           do i = 1+GRHydro_stencil, cctk_lsh(1)-GRHydro_stencil
              if ( GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0) .or. &
                   (atmosphere_mask(i,j,k) .ne. 0) .or. &
                      (dens(i,j,k) + frac * densrhs(i,j,k) .le. 0.d0) ) then
                 y_e_con_rhs(i,j,k) = 0.0d0
                 densrhs(i,j,k) = 0.0d0
                 srhs(i,j,k,:)   = 0.0d0
                 taurhs(i,j,k)  = 0.0d0
                 
                 if (evolve_Bvec) then
                    Bconsrhs(i,j,k,:) = 0.0d0
                 endif
                 if (evolve_Avec) then
                    Avecrhs(i,j,k,:) = 0.0d0
                 endif

                 ! TODO: Need to set Avecrhs and Aphirhs to zero for centered / or staggered case if vector potential is evolved

                 ! Set real-valued mask! This will be sync'ed and right after syncing translated to
                 ! our standard integer based mask (so that atmosphere_mask is still valid!).
                 atmosphere_mask_real(i,j,k) = 1
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

end subroutine GRHydroUpdateAtmosphereMask


subroutine GRHydroPostSyncAtmosphereMask(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k

!! This sets the integer atmo mask based on the real-valued (and sync'ed) atmo mask

   !$OMP PARALLEL DO PRIVATE(k,j,i)
   do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
         do i = 1, cctk_lsh(1)
            if ( abs(atmosphere_mask_real(i,j,k)) .gt. 0.5) then
               atmosphere_mask(i,j,k) = 1
            end if
         end do
      end do
   end do
   !$OMP END PARALLEL DO

end subroutine GRHydroPostSyncAtmosphereMask


 /*@@
   @routine    GRHydroCopyIntegerMask
   @date       Wed Jul  4 15:40:16 PDT 2012
   @author     Roland Haas
   @desc 
   Initializes real valued mask with integer valued one.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydroCopyIntegerMask(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k

!! This sets the real atmo mask based on the integer-valued atmo mask

   !$OMP PARALLEL DO PRIVATE(k,j,i)
   do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
         do i = 1, cctk_lsh(1)
            atmosphere_mask_real(i,j,k) = atmosphere_mask(i,j,k)
         end do
      end do
   end do
   !$OMP END PARALLEL DO

end subroutine GRHydroCopyIntegerMask


 /*@@
   @routine    GRHydro_SetupMask
   @date       Thu Jun 20 13:27:28 2002
   @author     Ian Hawke
   @desc 
   Initialize the mask to be zero.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_SetupMask(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  atmosphere_mask = 0
  atmosphere_mask_real = 0
  
  call CCTK_INFO("Setting up the atmosphere mask: all points are not_atmosphere")
  
end subroutine GRHydro_SetupMask

 /*@@
   @routine    GRHydro_InitAtmosMask
   @date       Thu Jun 20 13:27:28 2002
   @author     Ian Hawke
   @desc 
   Initialize the mask based on rho_min. This is used only if wk_atmosphere=yes.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_InitAtmosMask(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i,j,k
  CCTK_REAL :: dummy1, dummy2
  
  !$OMP PARALLEL DO PRIVATE(i,j,k, dummy1,dummy2)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        !if (rho(i,j,k) .le. GRHydro_rho_min) then
        IF_BELOW_ATMO(rho(i,j,k), GRHydro_rho_min, 0.0, r(i,j,k)) then
          atmosphere_mask(i,j,k) = 1
          atmosphere_mask_real(i,j,k) = 1
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  call CCTK_INFO("Setting up the atmosphere mask: points with rho<rho_min are set to be atmosphere")
  
end subroutine GRHydro_InitAtmosMask

 /*@@
   @routine    GRHydro_AtmosphereReset
   @date       Thu Jun 20 13:30:51 2002
   @author     Ian Hawke
   @desc 
   After MoL has evolved, if a point is supposed to be reset then do so.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_AtmosphereReset(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet, dummy1, dummy2


  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

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
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  if (verbose.eq.1) call CCTK_INFO("Entering AtmosphereReset.")

!$OMP PARALLEL DO PRIVATE(sdet,keytemp,i,j,k,anyerr,keyerr, dummy1, dummy2)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        if (atmosphere_mask(i, j, k) .ne. 0) then

          SET_ATMO_MIN(rho(i,j,k), GRHydro_rho_min, r(i,j,k))
          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0
          sdet = sdetg(i,j,k)
          
          if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
             temperature(i,j,k) = grhydro_hot_atmo_temp
             y_e(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)

             call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel,&
                  cctk_iteration,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),&
                  g11(i,j,k),g12(i,j,k),&
                  g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                  sdet,dens(i,j,k),scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                  tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
                  velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k))
             y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)
          else
             call prim2conpolytype(GRHydro_polytrope_handle, &
                  g11(i,j,k), g12(i,j,k), g13(i,j,k), &
                  g22(i,j,k), g23(i,j,k), g33(i,j,k), sdet, &
                  dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                  tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
                  velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))
             if (wk_atmosphere .eq. 0) then
                atmosphere_mask(i, j, k) = 0
                atmosphere_mask_real(i, j, k) = 0
             end if
          endif

        end if

      end do
    end do
  end do
!$OMP END PARALLEL DO


end subroutine GRHydro_AtmosphereReset

subroutine GRHydro_InitialAtmosphereReset(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet
  CCTK_REAL :: rho_min, dummy1, dummy2

  CCTK_INT :: eos_handle


  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11_p, g12_p, g13_p, g22_p, &
                                          g23_p, g33_p
  pointer (pg11_p,g11_p), (pg12_p,g12_p), (pg13_p,g13_p), (pg22_p,g22_p), &
                                          (pg23_p,g23_p), (pg33_p,g33_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup_p
  pointer (pvup_p,vup_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11_p_p, g12_p_p, g13_p_p, g22_p_p, &
                                          g23_p_p, g33_p_p
  pointer (pg11_p_p,g11_p_p), (pg12_p_p,g12_p_p), (pg13_p_p,g13_p_p), (pg22_p_p,g22_p_p), &
                                          (pg23_p_p, g23_p_p), (pg33_p_p, g33_p_p)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup_p_p
  pointer (pvup_p_p,vup_p_p)

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
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
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
    else
      pg11_p = loc(gxx_p)
      pg12_p = loc(gxy_p)
      pg13_p = loc(gxz_p)
      pg22_p = loc(gyy_p)
      pg23_p = loc(gyz_p)
      pg33_p = loc(gzz_p)
      pvup_p = loc(vel_p)
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
    else
      pg11_p_p = loc(gxx_p_p)
      pg12_p_p = loc(gxy_p_p)
      pg13_p_p = loc(gxz_p_p)
      pg22_p_p = loc(gyy_p_p)
      pg23_p_p = loc(gyz_p_p)
      pg33_p_p = loc(gzz_p_p)
      pvup_p_p = loc(vel_p_p)
    end if
  end if 


  eos_handle = GRHydro_polytrope_handle
  
  rho_min = GRHydro_rho_min
  if (initial_atmosphere_factor .gt. 0) then
    rho_min = rho_min * initial_atmosphere_factor
  endif

!$OMP PARALLEL DO PRIVATE(sdet,keytemp,i,j,k,anyerr,keyerr, dummy1, dummy2)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        SET_ATMO_MIN(dummy2, rho_min, r(i,j,k))
        if (rho(i,j,k) .le. dummy2 .or. &
            GRHydro_enable_internal_excision /= 0 .and. &
            hydro_excision_mask(i,j,k) .ne. 0) then 
          SET_ATMO_MIN(rho(i,j,k), dummy2, r(i,j,k))
          sdet = sdetg(i,j,k)

          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0

          if(evolve_temper.ne.0) then
!             ! set the temperature to be relatively low
             temperature(i,j,k) = grhydro_hot_atmo_temp
             y_e(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)

             call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel,&
                  cctk_iteration,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),&
                  g11(i,j,k),g12(i,j,k),&
                  g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                  sdet,dens(i,j,k),scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                  tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
                  velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k))
             y_e_con(i,j,k) = dens(i,j,k) * y_e(i,j,k)
          else
             keytemp = 0
             call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
             call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),xeps,xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
             call prim2conpolytype(eos_handle, &
                  g11(i,j,k), g12(i,j,k), g13(i,j,k), &
                  g22(i,j,k), g23(i,j,k), g33(i,j,k), sdet, &
                  dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                  tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
                  velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))
          endif
        end if
        if (timelevels .gt. 1) then
          if (rho_p(i,j,k) .le. dummy2) then
            SET_ATMO_MIN(rho_p(i,j,k), dummy2, r(i,j,k))
            velx_p(i,j,k) = 0.0d0
            vely_p(i,j,k) = 0.0d0
            velz_p(i,j,k) = 0.0d0

            sdet = sqrt(SPATIAL_DETERMINANT(g11_p(i,j,k), g12_p(i,j,k), g13_p(i,j,k), \
                  g22_p(i,j,k), g23_p(i,j,k), g33_p(i,j,k)))

            if(evolve_temper.ne.0) then
             ! set the temperature to be relatively low
             temperature_p(i,j,k) = grhydro_hot_atmo_temp
             y_e_p(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho_p(i,j,k),eps_p(i,j,k),temperature_p(i,j,k),y_e_p(i,j,k),&
                  press_p(i,j,k),keyerr,anyerr)

             call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel,&
                  cctk_iteration,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),&
                  g11_p(i,j,k),g12_p(i,j,k),&
                  g13_p(i,j,k),g22_p(i,j,k),g23_p(i,j,k),g33_p(i,j,k), &
                  sdet,dens_p(i,j,k),scon_p(i,j,k,1), scon_p(i,j,k,2), scon_p(i,j,k,3), &
                  tau_p(i,j,k), rho_p(i,j,k), velx_p(i,j,k), vely_p(i,j,k), &
                  velz_p(i,j,k), eps_p(i,j,k), press_p(i,j,k), w_lorentz_p(i,j,k),&
                  temperature_p(i,j,k),y_e_p(i,j,k))
             y_e_con_p(i,j,k) = dens_p(i,j,k) * y_e_p(i,j,k)
            else
               keytemp = 0
               call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p(i,j,k),eps_p(i,j,k),xtemp,xye,press_p(i,j,k),keyerr,anyerr)
               call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p(i,j,k),xeps,xtemp,xye,press_p(i,j,k),eps_p(i,j,k),keyerr,anyerr)
               call prim2conpolytype(eos_handle, &
                    g11_p(i,j,k), g12_p(i,j,k), g13_p(i,j,k), &
                    g22_p(i,j,k), g23_p(i,j,k), g33_p(i,j,k), sdet, &
                    dens_p(i,j,k), scon_p(i,j,k,1), scon_p(i,j,k,2), scon_p(i,j,k,3), &
                    tau_p(i,j,k), rho_p(i,j,k), velx_p(i,j,k), vely_p(i,j,k), &
                    velz_p(i,j,k), eps_p(i,j,k), press_p(i,j,k), w_lorentz_p(i,j,k))
            endif

          endif
        end if
        if (timelevels .gt. 2) then
          if (rho_p_p(i,j,k) .le. dummy2) then
            SET_ATMO_MIN(rho_p_p(i,j,k), dummy2, r(i,j,k))
            velx_p_p(i,j,k) = 0.0d0
            vely_p_p(i,j,k) = 0.0d0
            velz_p_p(i,j,k) = 0.0d0
            sdet = sqrt(SPATIAL_DETERMINANT(g11_p_p(i,j,k), g12_p_p(i,j,k), g13_p_p(i,j,k), \
                  g22_p_p(i,j,k), g23_p_p(i,j,k), g33_p_p(i,j,k)))

            if(evolve_temper.ne.0) then
             ! set the temperature to be relatively low
             temperature_p_p(i,j,k) = grhydro_hot_atmo_temp
             y_e_p_p(i,j,k) = grhydro_hot_atmo_Y_e
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho_p_p(i,j,k),eps_p_p(i,j,k),temperature_p_p(i,j,k),y_e_p_p(i,j,k),&
                  press_p_p(i,j,k),keyerr,anyerr)
             call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel,&
                  cctk_iteration,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),&
                  g11_p_p(i,j,k),g12_p_p(i,j,k),&
                  g13_p_p(i,j,k),g22_p_p(i,j,k),g23_p_p(i,j,k),g33_p_p(i,j,k), &
                  sdet,dens_p_p(i,j,k),scon_p_p(i,j,k,1), scon_p_p(i,j,k,2), scon_p_p(i,j,k,3), &
                  tau_p_p(i,j,k), rho_p_p(i,j,k), velx_p_p(i,j,k), vely_p_p(i,j,k), &
                  velz_p_p(i,j,k), eps_p_p(i,j,k), press_p_p(i,j,k), w_lorentz_p_p(i,j,k),&
                  temperature_p_p(i,j,k),y_e_p_p(i,j,k))
             y_e_con_p_p(i,j,k) = dens_p_p(i,j,k) * y_e_p_p(i,j,k)
            else
               keytemp = 0
               call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p_p(i,j,k),eps_p_p(i,j,k),xtemp,xye,press_p_p(i,j,k),keyerr,anyerr)
               call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                    rho_p_p(i,j,k),xeps,xtemp,xye,press_p_p(i,j,k),eps_p_p(i,j,k),keyerr,anyerr)
               call prim2conpolytype(eos_handle, &
                    g11_p_p(i,j,k), g12_p_p(i,j,k), g13_p_p(i,j,k), &
                    g22_p_p(i,j,k), g23_p_p(i,j,k), g33_p_p(i,j,k), sdet, &
                    dens_p_p(i,j,k), scon_p_p(i,j,k,1), scon_p_p(i,j,k,2), scon_p_p(i,j,k,3), &
                    tau_p_p(i,j,k), rho_p_p(i,j,k), velx_p_p(i,j,k), vely_p_p(i,j,k), &
                    velz_p_p(i,j,k), eps_p_p(i,j,k), press_p_p(i,j,k), w_lorentz_p_p(i,j,k))
            endif
          endif
        endif

      end do
    end do
  end do
!$OMP END PARALLEL DO

end subroutine GRHydro_InitialAtmosphereReset

