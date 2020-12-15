 /*@@
   @file      GRHydro_PPMMReconstruct_drv.F90
   @date      Wed Jul 27 15:17:03 EDT 2011
   @author    Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   Driver routine to perform the magnetic version of PPM reconstruction.
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
#define Bconsx(i,j,k) Bcons(i,j,k,1)
#define Bconsy(i,j,k) Bcons(i,j,k,2)
#define Bconsz(i,j,k) Bcons(i,j,k,3)


 /*@@
   @routine    GRHydro_PPMMReconstruct_drv
   @date       Wed Jul 27 15:17:55 EDT 2011
   @author     Luca Baiotti, Ian Hawke, Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   A driver routine to do magnetic version of PPM reconstructions. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_PPMMReconstruct_drv(CCTK_ARGUMENTS)
  
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

  CCTK_REAL, dimension(:,:,:),allocatable :: &
       dum, dump, dumm

  CCTK_INT :: ierr

  allocate(trivial_rp(cctk_ash(1),cctk_ash(2),cctk_ash(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_ERROR("Allocation problems with trivial_rp")
    STOP
  end if
  
!!$ The dum variables are used as dummies if MHD on but divergence cleaning off
  allocate(dum(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       dump(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       dumm(cctk_ash(1),cctk_ash(2),cctk_ash(3)),STAT=ierr)

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

  !initialize trivial_rp to false                                               
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           trivial_rp(i,j,k) = .false.
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO 

!!$ PPM starts:
    if (flux_direction == 1) then
            ! constraint transport needs to be able to average fluxes in the directions
            ! other that flux_direction, which in turn need the primitives on interfaces
            !$OMP PARALLEL DO PRIVATE(i, j, k)
            do k = GRHydro_stencil, nz - GRHydro_stencil + 1 + transport_constraints
              do j = GRHydro_stencil, ny - GRHydro_stencil + 1 + transport_constraints
               if(clean_divergence.ne.0) then
                 if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                   call SimplePPM_1dM(GRHydro_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                        rho(:,j,k),vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3),&
                        Bvec(:,j,k,1),Bvec(:,j,k,2),Bvec(:,j,k,3),psidc(:,j,k),eps(:,j,k),press(:,j,k),&
                        rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                        Bvecxminus(:,j,k),Bvecyminus(:,j,k),Bveczminus(:,j,k),psidcminus(:,j,k),epsminus(:,j,k),&
                        rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                        Bvecxplus(:,j,k),Bvecyplus(:,j,k),Bveczplus(:,j,k),psidcplus(:,j,k),epsplus(:,j,k),&
                        1,trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                        gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                        gzz(:,j,k), betax(:,j,k), alp(:,j,k),&
                        w_lorentz(:,j,k), &
                        1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                        GRHydro_mppm_eigenvalue_x_right, &
                        GRHydro_mppm_xwind)
                 else
                   call SimplePPM_1dM(GRHydro_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                        rho(:,j,k),lvel(:,j,k,1),lvel(:,j,k,2),lvel(:,j,k,3),&
                        lBvec(:,j,k,1),lBvec(:,j,k,2),lBvec(:,j,k,3),psidc(:,j,k),eps(:,j,k),press(:,j,k),&
                        rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                        Bvecxminus(:,j,k),Bvecyminus(:,j,k),Bveczminus(:,j,k),psidcminus(:,j,k),epsminus(:,j,k),&
                        rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                        Bvecxplus(:,j,k),Bvecyplus(:,j,k),Bveczplus(:,j,k),psidcplus(:,j,k),epsplus(:,j,k),&
                        1,trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                        gaa(:,j,k), gab(:,j,k), gac(:,j,k), gbb(:,j,k), gbc(:,j,k), &
                        gcc(:,j,k), betaa(:,j,k), alp(:,j,k),&
                        w_lorentz(:,j,k), &
                        1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                        GRHydro_mppm_eigenvalue_x_right, &
                        GRHydro_mppm_xwind)
                 end if
               else  !clean_divergence
                 if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                   call SimplePPM_1dM(GRHydro_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                        rho(:,j,k),vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3),&
                        Bvec(:,j,k,1),Bvec(:,j,k,2),Bvec(:,j,k,3),dum(:,j,k),eps(:,j,k),press(:,j,k),&
                        rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                        Bvecxminus(:,j,k),Bvecyminus(:,j,k),Bveczminus(:,j,k),dumm(:,j,k),epsminus(:,j,k),&
                        rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                        Bvecxplus(:,j,k),Bvecyplus(:,j,k),Bveczplus(:,j,k),dump(:,j,k),epsplus(:,j,k),&
                        0,trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                        gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                        gzz(:,j,k), betax(:,j,k), alp(:,j,k),&
                        w_lorentz(:,j,k), &
                        1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                        GRHydro_mppm_eigenvalue_x_right, &
                        GRHydro_mppm_xwind)
                 else
                   call SimplePPM_1dM(GRHydro_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                        rho(:,j,k),lvel(:,j,k,1),lvel(:,j,k,2),lvel(:,j,k,3),&
                        lBvec(:,j,k,1),lBvec(:,j,k,2),lBvec(:,j,k,3),dum(:,j,k),eps(:,j,k),press(:,j,k),&
                        rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                        Bvecxminus(:,j,k),Bvecyminus(:,j,k),Bveczminus(:,j,k),dumm(:,j,k),epsminus(:,j,k),&
                        rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                        Bvecxplus(:,j,k),Bvecyplus(:,j,k),Bveczplus(:,j,k),dump(:,j,k),epsplus(:,j,k),&
                        0,trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                        gaa(:,j,k), gab(:,j,k), gac(:,j,k), gbb(:,j,k), gbc(:,j,k), &
                        gcc(:,j,k), betaa(:,j,k), alp(:,j,k),&
                        w_lorentz(:,j,k), &
                        1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                        GRHydro_mppm_eigenvalue_x_right, &
                        GRHydro_mppm_xwind)
                 end if
               endif  !clean_divergence
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

      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_tracer_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3), &
                      tracer(:,j,k,:),tracerminus(:,j,k,:),tracerplus(:,j,k,:), &
                      press(:,j,k))
               else
                 call SimplePPM_tracer_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
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
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1 + transport_constraints
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1 + transport_constraints
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_ye_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3), &
                      Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k), &
                      press(:,j,k))
               else
                 call SimplePPM_ye_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                      lvel(:,j,k,1),lvel(:,j,k,2),lvel(:,j,k,3), &
                      Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k), &
                      press(:,j,k))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if

    else if (flux_direction == 2) then
          ! constraint transport needs to be able to average fluxes in the directions
          ! other that flux_direction, which in turn need the primitives on interfaces
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1 + transport_constraints
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1 + transport_constraints
             if(clean_divergence.ne.0) then
              if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                call SimplePPM_1dM(GRHydro_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                     rho(j,:,k),vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1),&
                     Bvec(j,:,k,2),Bvec(j,:,k,3),Bvec(j,:,k,1),psidc(j,:,k),eps(j,:,k),press(j,:,k),&
                     rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                     Bvecyminus(j,:,k),Bveczminus(j,:,k),Bvecxminus(j,:,k),psidcminus(j,:,k),epsminus(j,:,k),&
                     rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                     Bvecyplus(j,:,k),Bveczplus(j,:,k),Bvecxplus(j,:,k),psidcplus(j,:,k),epsplus(j,:,k),&
                     1,trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                     gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                     gxx(j,:,k), betay(j,:,k), alp(j,:,k),&
                     w_lorentz(j,:,k), &
                     2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                     GRHydro_mppm_eigenvalue_y_right, &
                     GRHydro_mppm_xwind)
              else
                call SimplePPM_1dM(GRHydro_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                     rho(j,:,k),lvel(j,:,k,2),lvel(j,:,k,3),lvel(j,:,k,1),&
                     lBvec(j,:,k,2),lBvec(j,:,k,3),lBvec(j,:,k,1),psidc(j,:,k),eps(j,:,k),press(j,:,k),&
                     rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                     Bvecyminus(j,:,k),Bveczminus(j,:,k),Bvecxminus(j,:,k),psidcminus(j,:,k),epsminus(j,:,k),&
                     rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                     Bvecyplus(j,:,k),Bveczplus(j,:,k),Bvecxplus(j,:,k),psidcplus(j,:,k),epsplus(j,:,k),&
                     1,trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                     gbb(j,:,k), gbc(j,:,k), gab(j,:,k), gcc(j,:,k), gac(j,:,k), &
                     gaa(j,:,k), betab(j,:,k), alp(j,:,k),&
                     w_lorentz(j,:,k), &
                     2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                     GRHydro_mppm_eigenvalue_y_right, &
                     GRHydro_mppm_xwind)
              end if
             else  !clean_divergence
              if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                call SimplePPM_1dM(GRHydro_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                     rho(j,:,k),vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1),&
                     Bvec(j,:,k,2),Bvec(j,:,k,3),Bvec(j,:,k,1),dum(j,:,k),eps(j,:,k),press(j,:,k),&
                     rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                     Bvecyminus(j,:,k),Bveczminus(j,:,k),Bvecxminus(j,:,k),dumm(j,:,k),epsminus(j,:,k),&
                     rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                     Bvecyplus(j,:,k),Bveczplus(j,:,k),Bvecxplus(j,:,k),dump(j,:,k),epsplus(j,:,k),&
                     0,trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                     gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                     gxx(j,:,k), betay(j,:,k), alp(j,:,k),&
                     w_lorentz(j,:,k), &
                     2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                     GRHydro_mppm_eigenvalue_y_right, &
                     GRHydro_mppm_xwind)
              else
                call SimplePPM_1dM(GRHydro_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                     rho(j,:,k),lvel(j,:,k,2),lvel(j,:,k,3),lvel(j,:,k,1),&
                     lBvec(j,:,k,2),lBvec(j,:,k,3),lBvec(j,:,k,1),dum(j,:,k),eps(j,:,k),press(j,:,k),&
                     rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                     Bvecyminus(j,:,k),Bveczminus(j,:,k),Bvecxminus(j,:,k),dumm(j,:,k),epsminus(j,:,k),&
                     rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                     Bvecyplus(j,:,k),Bveczplus(j,:,k),Bvecxplus(j,:,k),dump(j,:,k),epsplus(j,:,k),&
                     0,trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                     gbb(j,:,k), gbc(j,:,k), gab(j,:,k), gcc(j,:,k), gac(j,:,k), &
                     gaa(j,:,k), betab(j,:,k), alp(j,:,k),&
                     w_lorentz(j,:,k), &
                     2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                     GRHydro_mppm_eigenvalue_y_right, &
                     GRHydro_mppm_xwind)
              end if
             endif  !clean_divergence
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
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                call SimplePPM_tracer_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                     vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1), &
                     tracer(j,:,k,:),tracerminus(j,:,k,:),tracerplus(j,:,k,:), &
                     press(j,:,k))
              else
                call SimplePPM_tracer_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
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
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1 + transport_constraints
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1 + transport_constraints
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_ye_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                      vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1), &
                      Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k), &
                      press(j,:,k))
               else
                 call SimplePPM_ye_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                      lvel(j,:,k,2),lvel(j,:,k,3),lvel(j,:,k,1), &
                      Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k), &
                      press(j,:,k))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      end if

    else if (flux_direction == 3) then
          ! constraint transport needs to be able to average fluxes in the directions
          ! other that flux_direction, which in turn need the primitives on interfaces
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, ny - GRHydro_stencil + 1 + transport_constraints
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1 + transport_constraints
             if(clean_divergence.ne.0) then
              if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                call SimplePPM_1dM(GRHydro_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                     rho(j,k,:),vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2),&
                     Bvec(j,k,:,3),Bvec(j,k,:,1),Bvec(j,k,:,2),psidc(j,k,:),eps(j,k,:),press(j,k,:),&
                     rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                     Bveczminus(j,k,:),Bvecxminus(j,k,:),Bvecyminus(j,k,:),psidcminus(j,k,:),epsminus(j,k,:),&
                     rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                     Bveczplus(j,k,:),Bvecxplus(j,k,:),Bvecyplus(j,k,:),psidcplus(j,k,:),epsplus(j,k,:),&
                     1,trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                     gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                     gyy(j,k,:), betaz(j,k,:), alp(j,k,:),&
                     w_lorentz(j,k,:), &
                     3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                     GRHydro_mppm_eigenvalue_z_right, &
                     GRHydro_mppm_xwind)
              else
                call SimplePPM_1dM(GRHydro_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                     rho(j,k,:),lvel(j,k,:,3),lvel(j,k,:,1),lvel(j,k,:,2),&
                     lBvec(j,k,:,3),lBvec(j,k,:,1),lBvec(j,k,:,2),psidc(j,k,:),eps(j,k,:),press(j,k,:),&
                     rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                     Bveczminus(j,k,:),Bvecxminus(j,k,:),Bvecyminus(j,k,:),psidcminus(j,k,:),epsminus(j,k,:),&
                     rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                     Bveczplus(j,k,:),Bvecxplus(j,k,:),Bvecyplus(j,k,:),psidcplus(j,k,:),epsplus(j,k,:),&
                     1,trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                     gcc(j,k,:), gac(j,k,:), gbc(j,k,:), gaa(j,k,:), gab(j,k,:), &
                     gbb(j,k,:), betac(j,k,:), alp(j,k,:),&
                     w_lorentz(j,k,:), &
                     3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                     GRHydro_mppm_eigenvalue_z_right, &
                     GRHydro_mppm_xwind)
              end if
             else  !clean_divergence
              if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                call SimplePPM_1dM(GRHydro_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                     rho(j,k,:),vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2),&
                     Bvec(j,k,:,3),Bvec(j,k,:,1),Bvec(j,k,:,2),dum(j,k,:),eps(j,k,:),press(j,k,:),&
                     rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                     Bveczminus(j,k,:),Bvecxminus(j,k,:),Bvecyminus(j,k,:),dumm(j,k,:),epsminus(j,k,:),&
                     rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                     Bveczplus(j,k,:),Bvecxplus(j,k,:),Bvecyplus(j,k,:),dump(j,k,:),epsplus(j,k,:),&
                     0,trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                     gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                     gyy(j,k,:), betaz(j,k,:), alp(j,k,:),&
                     w_lorentz(j,k,:), &
                     3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                     GRHydro_mppm_eigenvalue_z_right, &
                     GRHydro_mppm_xwind)
              else
                call SimplePPM_1dM(GRHydro_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                     rho(j,k,:),lvel(j,k,:,3),lvel(j,k,:,1),lvel(j,k,:,2),&
                     lBvec(j,k,:,3),lBvec(j,k,:,1),lBvec(j,k,:,2),dum(j,k,:),eps(j,k,:),press(j,k,:),&
                     rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                     Bveczminus(j,k,:),Bvecxminus(j,k,:),Bvecyminus(j,k,:),dumm(j,k,:),epsminus(j,k,:),&
                     rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                     Bveczplus(j,k,:),Bvecxplus(j,k,:),Bvecyplus(j,k,:),dump(j,k,:),epsplus(j,k,:),&
                     0,trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                     gcc(j,k,:), gac(j,k,:), gbc(j,k,:), gaa(j,k,:), gab(j,k,:), &
                     gbb(j,k,:), betac(j,k,:), alp(j,k,:),&
                     w_lorentz(j,k,:), &
                     3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                     GRHydro_mppm_eigenvalue_z_right, &
                     GRHydro_mppm_xwind)
              end if
             endif  !clean_divergence
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
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1

               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_tracer_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2), &
                      tracer(j,k,:,:),tracerminus(j,k,:,:),tracerplus(j,k,:,:), &
                      press(j,k,:))
               else
                 call SimplePPM_tracer_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
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
          do k = GRHydro_stencil, ny - GRHydro_stencil + 1 + transport_constraints
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1 + transport_constraints
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call SimplePPM_ye_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2), &
                      Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:), &
                      press(j,k,:))
               else
                 call SimplePPM_ye_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
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
  deallocate(dum,dump,dumm)

end subroutine GRHydro_PPMMReconstruct_drv

