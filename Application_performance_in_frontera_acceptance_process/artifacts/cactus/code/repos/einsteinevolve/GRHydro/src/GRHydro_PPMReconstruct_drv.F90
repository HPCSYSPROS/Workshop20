 /*@@
   @file      GRHydro_PPMReconstruct_drv.F90
   @date      Tue Jul 19 13:22:03 EDT 2011
   @author    Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   Driver routine to perform the PPM reconstructions.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)


 /*@@
   @routine    GRHydro_PPMReconstruct_drv
   @date       Tue Jul 19 13:24:34 EDT 2011
   @author     Luca Baiotti, Ian Hawke, Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   A driver routine to do PPM reconstructions. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_PPMReconstruct_drv(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates

  integer :: nx, ny, nz, i, j, k

  logical, dimension(:,:,:), allocatable :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz

  CCTK_INT :: ierr

  logical :: apply_enhanced_ppm

  allocate(trivial_rp(cctk_ash(1),cctk_ash(2),cctk_ash(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_ERROR("Allocation problems with trivial_rp")
    STOP
  end if
    
  call SpaceMask_GetTypeBits(type_bitsx, "Hydro_RiemannProblemX")
  call SpaceMask_GetStateBits(trivialx, "Hydro_RiemannProblemX", &
       &"trivial")
  call SpaceMask_GetStateBits(not_trivialx, "Hydro_RiemannProblemX", &
       &"not_trivial")
  call SpaceMask_GetTypeBits(type_bitsy, "Hydro_RiemannProblemY")
  call SpaceMask_GetStateBits(trivialy, "Hydro_RiemannProblemY", &
       &"trivial")
  call SpaceMask_GetStateBits(not_trivialy, "Hydro_RiemannProblemY", &
       &"not_trivial")
  call SpaceMask_GetTypeBits(type_bitsz, "Hydro_RiemannProblemZ")
  call SpaceMask_GetStateBits(trivialz, "Hydro_RiemannProblemZ", &
       &"trivial")
  call SpaceMask_GetStateBits(not_trivialz, "Hydro_RiemannProblemZ", &
       &"not_trivial")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  ! if use_enhanced_ppm, allow old PPM on one level
  if (GRHydro_oppm_reflevel .eq. (-1) .or. &
       GRHydro_reflevel .ne. GRHydro_oppm_reflevel) then
     apply_enhanced_ppm = use_enhanced_ppm .ne. 0
  else
     apply_enhanced_ppm = .false.
  end if


!!$ PPM starts:
    if (flux_direction == 1) then
         !$OMP PARALLEL DO PRIVATE(i, j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
           do j = GRHydro_stencil, ny - GRHydro_stencil + 1
           if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
             call SimplePPM_1d(apply_enhanced_ppm,&
                  GRHydro_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                  rho(:,j,k),vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3),eps(:,j,k),&
                  press(:,j,k),rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                  velzminus(:,j,k),epsminus(:,j,k),rhoplus(:,j,k),&
                  velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),epsplus(:,j,k),&
                  trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                  gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                  gzz(:,j,k), betax(:,j,k), alp(:,j,k),&
                  w_lorentz(:,j,k), &
                  1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                  GRHydro_mppm_eigenvalue_x_right, &
                  GRHydro_mppm_xwind)
           else
             call SimplePPM_1d(apply_enhanced_ppm,&
                  GRHydro_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                  rho(:,j,k),lvel(:,j,k,1),lvel(:,j,k,2),lvel(:,j,k,3),eps(:,j,k),&
                  press(:,j,k),rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                  velzminus(:,j,k),epsminus(:,j,k),rhoplus(:,j,k),&
                  velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),epsplus(:,j,k),&
                  trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                  gaa(:,j,k), gab(:,j,k), gac(:,j,k), gbb(:,j,k), gbc(:,j,k), &
                  gcc(:,j,k), betaa(:,j,k), alp(:,j,k),&
                  w_lorentz(:,j,k), &
                  1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                  GRHydro_mppm_eigenvalue_x_right, &
                  GRHydro_mppm_xwind)
           end if
             do i = 1, nx
                if (trivial_rp(i,j,k)) then
                  SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, trivialx)
                else
                  SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, not_trivialx)
                end if
             end do
           end do
         end do
         !$OMP END PARALLEL DO
      if(evolve_temper.ne.0.and.reconstruct_temper.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_temperature_1d(apply_enhanced_ppm,&
                      nx,CCTK_DELTA_SPACE(1),vel(:,j,k,1), & 
                      temperature(:,j,k),press(:,j,k),&
                      tempminus(:,j,k),tempplus(:,j,k), hydro_excision_mask)
               else
                   call SimplePPM_temperature_1d(apply_enhanced_ppm,&
                      nx,CCTK_DELTA_SPACE(1),lvel(:,j,k,1), & 
                      temperature(:,j,k),press(:,j,k),&
                      tempminus(:,j,k),tempplus(:,j,k), hydro_excision_mask)
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      endif
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                      nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3), &
                      tracer(:,j,k,:),tracerminus(:,j,k,:),tracerplus(:,j,k,:), &
                      press(:,j,k))
               else
                 call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                      nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      lvel(:,j,k,1),lvel(:,j,k,2),lvel(:,j,k,3), &
                      tracer(:,j,k,:),tracerminus(:,j,k,:),tracerplus(:,j,k,:), &
                      press(:,j,k))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_ye_1d(apply_enhanced_ppm,&
                      nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3), &
                      Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k), &
                      press(:,j,k))
               else
                 call SimplePPM_ye_1d(apply_enhanced_ppm,&
                      nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      lvel(:,j,k,1),lvel(:,j,k,2),lvel(:,j,k,3), &
                      Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k), &
                      press(:,j,k))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if

    else if (flux_direction == 2) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, nz - GRHydro_stencil + 1
          do j = GRHydro_stencil, nx - GRHydro_stencil + 1
            if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
              call SimplePPM_1d(apply_enhanced_ppm,&
                   GRHydro_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                   rho(j,:,k),vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1),eps(j,:,k),&
                   press(j,:,k),rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                   velxminus(j,:,k),epsminus(j,:,k),rhoplus(j,:,k),&
                   velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),epsplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                   gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                   gxx(j,:,k), betay(j,:,k), alp(j,:,k),&
                   w_lorentz(j,:,k), &
                   2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                   GRHydro_mppm_eigenvalue_y_right, &
                   GRHydro_mppm_xwind)
            else
              call SimplePPM_1d(apply_enhanced_ppm,&
                   GRHydro_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                   rho(j,:,k),lvel(j,:,k,2),lvel(j,:,k,3),lvel(j,:,k,1),eps(j,:,k),&
                   press(j,:,k),rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                   velxminus(j,:,k),epsminus(j,:,k),rhoplus(j,:,k),&
                   velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),epsplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                   gbb(j,:,k), gbc(j,:,k), gab(j,:,k), gcc(j,:,k), gac(j,:,k), &
                   gaa(j,:,k), betab(j,:,k), alp(j,:,k),&
                   w_lorentz(j,:,k), &
                   2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                   GRHydro_mppm_eigenvalue_y_right, &
                   GRHydro_mppm_xwind)
            end if
            do i = 1, ny
              if (trivial_rp(j,i,k)) then
                 SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, trivialy)
              else
                 SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, not_trivialy)
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      if(evolve_temper.ne.0.and.reconstruct_temper.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_temperature_1d(apply_enhanced_ppm,&
                      ny,CCTK_DELTA_SPACE(2),vel(j,:,k,2), & 
                      temperature(j,:,k),press(j,:,k),&
                      tempminus(j,:,k),tempplus(j,:,k), hydro_excision_mask)
               else
                 call SimplePPM_temperature_1d(apply_enhanced_ppm,&
                      ny,CCTK_DELTA_SPACE(2),lvel(j,:,k,2), & 
                      temperature(j,:,k),press(j,:,k),&
                      tempminus(j,:,k),tempplus(j,:,k), hydro_excision_mask)
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      endif
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                     ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                     vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1), &
                     tracer(j,:,k,:),tracerminus(j,:,k,:),tracerplus(j,:,k,:), &
                     press(j,:,k))
              else
                call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                     ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                     lvel(j,:,k,2),lvel(j,:,k,3),lvel(j,:,k,1), &
                     tracer(j,:,k,:),tracerminus(j,:,k,:),tracerplus(j,:,k,:), &
                     press(j,:,k))
              end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_ye_1d(apply_enhanced_ppm,&
                      ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                      vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1), &
                      Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k), &
                      press(j,:,k))
               else
                 call SimplePPM_ye_1d(apply_enhanced_ppm,&
                      ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                      lvel(j,:,k,2),lvel(j,:,k,3),lvel(j,:,k,1), &
                      Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k), &
                      press(j,:,k))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if

    else if (flux_direction == 3) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, ny - GRHydro_stencil + 1
          do j = GRHydro_stencil, nx - GRHydro_stencil + 1
            if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
              call SimplePPM_1d(apply_enhanced_ppm,&
                   GRHydro_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                   rho(j,k,:),vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2),eps(j,k,:),&
                   press(j,k,:),rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
                   velyminus(j,k,:),epsminus(j,k,:),rhoplus(j,k,:),&
                   velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),epsplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                   gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                   gyy(j,k,:), betaz(j,k,:), alp(j,k,:),&
                   w_lorentz(j,k,:), &
                   3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                   GRHydro_mppm_eigenvalue_z_right, &
                   GRHydro_mppm_xwind)
            else
              call SimplePPM_1d(apply_enhanced_ppm,&
                   GRHydro_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                   rho(j,k,:),lvel(j,k,:,3),lvel(j,k,:,1),lvel(j,k,:,2),eps(j,k,:),&
                   press(j,k,:),rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
                   velyminus(j,k,:),epsminus(j,k,:),rhoplus(j,k,:),&
                   velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),epsplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                   gcc(j,k,:), gac(j,k,:), gbc(j,k,:), gaa(j,k,:), gab(j,k,:), &
                   gbb(j,k,:), betac(j,k,:), alp(j,k,:),&
                   w_lorentz(j,k,:), &
                   3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                   GRHydro_mppm_eigenvalue_z_right, &
                   GRHydro_mppm_xwind)
            end if
            do i = 1, nz
              if (trivial_rp(j,k,i)) then
                SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, trivialz)
              else
                SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, not_trivialz)
              end if
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      if(evolve_temper.ne.0.and.reconstruct_temper.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_temperature_1d(apply_enhanced_ppm,&
                      nz,CCTK_DELTA_SPACE(3),vel(j,k,:,3), & 
                      temperature(j,k,:),press(j,k,:),&
                      tempminus(j,k,:),tempplus(j,k,:), hydro_excision_mask)
               else
                 call SimplePPM_temperature_1d(apply_enhanced_ppm,&
                      nz,CCTK_DELTA_SPACE(3),lvel(j,k,:,3), & 
                      temperature(j,k,:),press(j,k,:),&
                      tempminus(j,k,:),tempplus(j,k,:), hydro_excision_mask)
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      endif
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1

               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                      nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2), &
                      tracer(j,k,:,:),tracerminus(j,k,:,:),tracerplus(j,k,:,:), &
                      press(j,k,:))
               else
                 call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                      nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      lvel(j,k,:,3),lvel(j,k,:,1),lvel(j,k,:,2), &
                      tracer(j,k,:,:),tracerminus(j,k,:,:),tracerplus(j,k,:,:), &
                      press(j,k,:))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_ye_1d(apply_enhanced_ppm,&
                      nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2), &
                      Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:), &
                      press(j,k,:))
               else
                 call SimplePPM_ye_1d(apply_enhanced_ppm,&
                      nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      lvel(j,k,:,3),lvel(j,k,:,1),lvel(j,k,:,2), &
                      Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:), &
                      press(j,k,:))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if

    else
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if
!!$ PPM ends.

  deallocate(trivial_rp)

end subroutine GRHydro_PPMReconstruct_drv
