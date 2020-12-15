 /*@@
   @file      GRHydro_ReconstructPoly.F90
   @date      Sat Jan 26 02:13:25 2002
   @author    
   @desc 
   Wrapper routine to perform the reconstruction for polytropes.
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
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)
#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define Bconsx(i,j,k) Bcons(i,j,k,1)
#define Bconsy(i,j,k) Bcons(i,j,k,2)
#define Bconsz(i,j,k) Bcons(i,j,k,3)


 /*@@
   @routine    ReconstructionPolytype
   @date       Tue Mar 19 11:36:55 2002
   @author     Ian Hawke  
   @desc 
   If using a polytropic type EOS, we do not have to reconstruct all the 
   variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine ReconstructionPolytype(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  logical :: apply_enhanced_ppm

  integer :: nx, ny, nz, i, j, k, itracer

  logical, dimension(:,:,:), allocatable :: trivial_rp
!!$  logical, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz

  CCTK_REAL, dimension(:,:,:),allocatable :: &
       &lbetax, lbetay, lbetaz, dum, dump, dumm
!!$  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: &
!!$       &lbetax, lbetay, lbetaz

  CCTK_INT :: ierr

  CCTK_REAL :: local_min_tracer

  CCTK_INT :: GRHydro_UseGeneralCoordinates

  allocate(trivial_rp(cctk_ash(1),cctk_ash(2),cctk_ash(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_ERROR("Allocation problems with trivial_rp")
    STOP
  end if
  
  allocate(lbetax(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       lbetay(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       lbetaz(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       dum(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       dump(cctk_ash(1),cctk_ash(2),cctk_ash(3)),&
       dumm(cctk_ash(1),cctk_ash(2),cctk_ash(3)),STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_ERROR("Allocation problems with lbeta")
    STOP
  end if
  
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    call CCTK_ERROR("MP not yet supported in GRHydro_ReconstructionPolytype")
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

  ! if use_enhanced_ppm, allow old PPM on one level
  if (GRHydro_oppm_reflevel .eq. (-1) .or. &
       GRHydro_reflevel .ne. GRHydro_oppm_reflevel) then
     apply_enhanced_ppm = use_enhanced_ppm .ne. 0
  else
     apply_enhanced_ppm = .false.
  end if

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

!!$  Currently only option is reconstruction on primitive variables.
!!$  Should change this.

  ! Initialize plus and minus states
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           lbetax(i,j,k) = betax(i,j,k)
           lbetay(i,j,k) = betay(i,j,k)
           lbetaz(i,j,k) = betaz(i,j,k)
           trivial_rp(i,j,k) = .false.
           ! must initialize rho and eps plus minus
           ! to cell average in order to have sane values
           ! in the boundary zones for EOS call
           rhoplus(i,j,k) = rho(i,j,k)
           rhominus(i,j,k)= rho(i,j,k)
           epsplus(i,j,k) = eps(i,j,k)
           epsminus(i,j,k) = eps(i,j,k)
           velxplus(i,j,k) = 0.0d0
           velxminus(i,j,k) = 0.0d0
           velyplus(i,j,k) = 0.0d0
           velyminus(i,j,k) = 0.0d0
           velzplus(i,j,k) = 0.0d0
           velzminus(i,j,k) = 0.0d0

           if(evolve_mhd.ne.0) then
              Bvecxplus(i,j,k) = 0.0d0
              Bvecxminus(i,j,k) = 0.0d0
              Bvecyplus(i,j,k) = 0.0d0
              Bvecyminus(i,j,k) = 0.0d0
              Bveczplus(i,j,k) = 0.0d0
              Bveczminus(i,j,k) = 0.0d0
              if(clean_divergence.ne.0) then
                 psidcplus(i,j,k) = 0.0d0
                 psidcminus(i,j,k) = 0.0d0
              endif
           endif
           if (evolve_entropy .ne. 0) then
              entropyplus(i,j,k) = 0.0d0
              entropyminus(i,j,k) = 0.0d0
           endif

           if (evolve_tracer .ne. 0) then
              tracerplus(i,j,k,:) = 0.0d0
              tracerminus(i,j,k,:) = 0.0d0
           endif

           if (evolve_Y_e .ne. 0) then
              ! set this to the cell center values
              ! to make sure we have good Y_e even in
              ! the boundary region (for full GF EOS calls)
              Y_e_plus(i,j,k) = Y_e(i,j,k)
              Y_e_minus(i,j,k) = Y_e(i,j,k)
           endif
          if(evolve_temper .ne. 0) then
              ! set this to cell center value to have
              ! good initial guesses at interfaces
              ! in case we don't reconstruct temp
              tempplus(i,j,k) = temperature(i,j,k)
              tempminus(i,j,k) = temperature(i,j,k)
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO             



  if (CCTK_EQUALS(recon_method,"tvd")) then

     if (evolve_tracer .ne. 0) then
        do itracer=1,number_of_tracers
           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                tracer(:,:,:,itracer), tracerplus(:,:,:,itracer), &
                tracerminus(:,:,:,itracer), &
                trivial_rp, hydro_excision_mask)
        enddo
     end if
     
     if (evolve_Y_e .ne. 0) then
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Y_e(:,:,:), Y_e_plus(:,:,:), &
             Y_e_minus(:,:,:), &
             trivial_rp, hydro_excision_mask)
     endif

    if (CCTK_EQUALS(recon_vars,"primitive")) then
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           rho, rhoplus, rhominus, trivial_rp, hydro_excision_mask)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           vel(:,:,:,1), velxplus, velxminus, trivial_rp, hydro_excision_mask)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           vel(:,:,:,2), velyplus, velyminus, trivial_rp, hydro_excision_mask)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           vel(:,:,:,3), velzplus, velzminus, trivial_rp, hydro_excision_mask)
      if(evolve_mhd.ne.0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bvec(:,:,:,1), Bvecxplus, Bvecxminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bvec(:,:,:,2), Bvecyplus, Bvecyminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bvec(:,:,:,3), Bveczplus, Bveczminus, trivial_rp, hydro_excision_mask)
      endif
    else if (CCTK_EQUALS(recon_vars,"conservative")) then
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           dens, densplus, densminus, trivial_rp, hydro_excision_mask)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           scon(:,:,:,1), sxplus, sxminus, trivial_rp, hydro_excision_mask)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           scon(:,:,:,2), syplus, syminus, trivial_rp, hydro_excision_mask)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           scon(:,:,:,3), szplus, szminus, trivial_rp, hydro_excision_mask)
      if(evolve_mhd.ne.0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bcons(:,:,:,1), Bconsxplus, Bconsxminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bcons(:,:,:,2), Bconsyplus, Bconsyminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bcons(:,:,:,3), Bconszplus, Bconszminus, trivial_rp, hydro_excision_mask)
      endif
    else
      call CCTK_ERROR("Variable type to reconstruct not recognized.")
      STOP
    end if

    if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
       call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
            psidc, psidcplus, psidcminus, trivial_rp, hydro_excision_mask)
    endif

    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (trivial_rp(i,j,k)) then
            if (flux_direction == 1) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, trivialx)
            else if (flux_direction == 2) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsy, trivialy)
            else if (flux_direction == 3) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsz, trivialz)
            end if
          else
            if (flux_direction == 1) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, not_trivialx)
            else if (flux_direction == 2) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsy, not_trivialy)
            else if (flux_direction == 3) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsz, not_trivialz)
            end if
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  else if (CCTK_EQUALS(recon_method,"ppm")) then

    if (flux_direction == 1) then
      if(evolve_mhd.ne.0) then
        if(clean_divergence.ne.0) then
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
              call SimplePPM_1dM(GRHydro_eos_handle,1,nx,CCTK_DELTA_SPACE(1),&
                   rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),&
                   Bvecx(:,j,k),Bvecy(:,j,k),Bvecz(:,j,k),psidc(:,j,k),eps(:,j,k),press(:,j,k),&
                   rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                   Bvecxminus(:,j,k),Bvecyminus(:,j,k),Bveczminus(:,j,k),psidcminus(:,j,k),epsminus(:,j,k),&
                   rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                   Bvecxplus(:,j,k),Bvecyplus(:,j,k),Bveczplus(:,j,k),psidcplus(:,j,k),epsplus(:,j,k),&
                   1,trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                   gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                   gzz(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                   w_lorentz(:,j,k), &
                   1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                   GRHydro_mppm_eigenvalue_x_right, &
                   GRHydro_mppm_xwind)
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
        else  !clean_divergence
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
              call SimplePPM_1dM(GRHydro_eos_handle,1,nx,CCTK_DELTA_SPACE(1),&
                   rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),&
                   Bvecx(:,j,k),Bvecy(:,j,k),Bvecz(:,j,k),dum(:,j,k),eps(:,j,k),press(:,j,k),&
                   rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                   Bvecxminus(:,j,k),Bvecyminus(:,j,k),Bveczminus(:,j,k),dumm(:,j,k),epsminus(:,j,k),&
                   rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                   Bvecxplus(:,j,k),Bvecyplus(:,j,k),Bveczplus(:,j,k),dump(:,j,k),epsplus(:,j,k),&
                   0,trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                   gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                   gzz(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                   w_lorentz(:,j,k), &
                   1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                   GRHydro_mppm_eigenvalue_x_right, &
                   GRHydro_mppm_xwind)
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
        endif  !clean_divergence
      else  !evolve_mhd
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, nz - GRHydro_stencil + 1
          do j = GRHydro_stencil, ny - GRHydro_stencil + 1
            call SimplePPM_1d(apply_enhanced_ppm,&
                 GRHydro_eos_handle,1,nx,CCTK_DELTA_SPACE(1),&
                 rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),eps(:,j,k),&
                 press(:,j,k),rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                 velzminus(:,j,k),epsminus(:,j,k),rhoplus(:,j,k),&
                 velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),epsplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k),&
                 gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                 gzz(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                 w_lorentz(:,j,k), &
                 1, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_x_left, &
                 GRHydro_mppm_eigenvalue_x_right, &
                 GRHydro_mppm_xwind)
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
      endif  !evolve_mhd
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
               call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                    nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                    velx(:,j,k),vely(:,j,k),velz(:,j,k), &
                    tracer(:,j,k,:),tracerminus(:,j,k,:),tracerplus(:,j,k,:), &
                    press(:,j,k))
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, ny - GRHydro_stencil + 1
               call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                    nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                    velx(:,j,k),vely(:,j,k),velz(:,j,k), &
                    Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k), &
                    press(:,j,k))
            end do
         end do
         !$OMP END PARALLEL DO
      end if

    else if (flux_direction == 2) then
      if(evolve_mhd.ne.0) then
        if(clean_divergence.ne.0) then
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              call SimplePPM_1dM(GRHydro_eos_handle,1,ny,CCTK_DELTA_SPACE(2),&
                   rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),&
                   Bvecy(j,:,k),Bvecz(j,:,k),Bvecx(j,:,k),psidc(j,:,k),eps(j,:,k),press(j,:,k),&
                   rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                   Bvecyminus(j,:,k),Bveczminus(j,:,k),Bvecxminus(j,:,k),psidcminus(j,:,k),epsminus(j,:,k),&
                   rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                   Bvecyplus(j,:,k),Bveczplus(j,:,k),Bvecxplus(j,:,k),psidcplus(j,:,k),epsplus(j,:,k),&
                   1,trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                   gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                   gxx(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                   w_lorentz(j,:,k), &
                   2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                   GRHydro_mppm_eigenvalue_y_right, &
                   GRHydro_mppm_xwind)
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
        else  !clean_divergence
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              call SimplePPM_1dM(GRHydro_eos_handle,1,ny,CCTK_DELTA_SPACE(2),&
                   rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),&
                   Bvecy(j,:,k),Bvecz(j,:,k),Bvecx(j,:,k),dum(j,:,k),eps(j,:,k),press(j,:,k),&
                   rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                   Bvecyminus(j,:,k),Bveczminus(j,:,k),Bvecxminus(j,:,k),dumm(j,:,k),epsminus(j,:,k),&
                   rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                   Bvecyplus(j,:,k),Bveczplus(j,:,k),Bvecxplus(j,:,k),dump(j,:,k),epsplus(j,:,k),&
                   0,trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                   gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                   gxx(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                   w_lorentz(j,:,k), &
                   2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                   GRHydro_mppm_eigenvalue_y_right, &
                   GRHydro_mppm_xwind)
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
        endif  !clean_divergence
      else  !evolve_mhd
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, nz - GRHydro_stencil + 1
          do j = GRHydro_stencil, nx - GRHydro_stencil + 1
            call SimplePPM_1d(apply_enhanced_ppm,&
                 GRHydro_eos_handle,1,ny,CCTK_DELTA_SPACE(2),&
                 rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),eps(j,:,k),&
                 press(j,:,k),rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                 velxminus(j,:,k),epsminus(j,:,k),rhoplus(j,:,k),&
                 velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),epsplus(j,:,k),&
                 trivial_rp(j,:,k), hydro_excision_mask(j,:,k),&
                 gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                 gxx(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                 w_lorentz(j,:,k), &
                 2, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_y_left, &
                 GRHydro_mppm_eigenvalue_y_right, &
                 GRHydro_mppm_xwind)
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
      endif  !evolve_mhd
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                   ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                   vely(j,:,k),velz(j,:,k),velx(j,:,k), &
                   tracer(j,:,k,:),tracerminus(j,:,k,:),tracerplus(j,:,k,:), &
                   press(j,:,k))
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, nz - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
               call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                    ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                    velx(j,:,k),vely(j,:,k),velz(j,:,k), &
                    Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k), &
                    press(j,:,k))
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      
    else if (flux_direction == 3) then
      if(evolve_mhd.ne.0) then
        if(clean_divergence.ne.0) then
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              call SimplePPM_1dM(GRHydro_eos_handle,1,nz,CCTK_DELTA_SPACE(3),&
                   rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),&
                   Bvecz(j,k,:),Bvecx(j,k,:),Bvecy(j,k,:),psidc(j,k,:),eps(j,k,:),press(j,k,:),&
                   rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                   Bveczminus(j,k,:),Bvecxminus(j,k,:),Bvecyminus(j,k,:),psidcminus(j,k,:),epsminus(j,k,:),&
                   rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                   Bveczplus(j,k,:),Bvecxplus(j,k,:),Bvecyplus(j,k,:),psidcplus(j,k,:),epsplus(j,k,:),&
                   1,trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                   gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                   gyy(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
                   w_lorentz(j,k,:), &
                   3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                   GRHydro_mppm_eigenvalue_z_right, &
                   GRHydro_mppm_xwind)
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
        else  !clean_divergence
          !$OMP PARALLEL DO PRIVATE(i, j, k)
          do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
              call SimplePPM_1dM(GRHydro_eos_handle,1,nz,CCTK_DELTA_SPACE(3),&
                   rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),&
                   Bvecz(j,k,:),Bvecx(j,k,:),Bvecy(j,k,:),dum(:,j,k),eps(j,k,:),press(j,k,:),&
                   rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                   Bveczminus(j,k,:),Bvecxminus(j,k,:),Bvecyminus(j,k,:),dumm(j,k,:),epsminus(j,k,:),&
                   rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                   Bveczplus(j,k,:),Bvecxplus(j,k,:),Bvecyplus(j,k,:),dump(j,k,:),epsplus(j,k,:),&
                   0,trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                   gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                   gyy(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
                   w_lorentz(j,k,:), &
                   3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                   GRHydro_mppm_eigenvalue_z_right, &
                   GRHydro_mppm_xwind)
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
        endif  !clean_divergence
      else  !evolve_mhd
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, ny - GRHydro_stencil + 1
          do j = GRHydro_stencil, nx - GRHydro_stencil + 1
            call SimplePPM_1d(apply_enhanced_ppm,&
                 GRHydro_eos_handle,1,nz,CCTK_DELTA_SPACE(3),&
                 rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),eps(j,k,:),&
                 press(j,k,:),rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
                 velyminus(j,k,:),epsminus(j,k,:),rhoplus(j,k,:),&
                 velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),epsplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:),&
                 gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                 gyy(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
                 w_lorentz(j,k,:), &
                 3, j, k, nx, ny, nz, GRHydro_mppm_eigenvalue_z_left, &
                 GRHydro_mppm_eigenvalue_z_right, &
                 GRHydro_mppm_xwind)
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
      endif  !evolve_mhd
      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1

               call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                    nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                    velz(j,k,:),velx(j,k,:),vely(j,k,:), &
                    tracer(j,k,:,:),tracerminus(j,k,:,:),tracerplus(j,k,:,:), &
                    press(j,k,:))
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
         do k = GRHydro_stencil, ny - GRHydro_stencil + 1
            do j = GRHydro_stencil, nx - GRHydro_stencil + 1
               call SimplePPM_tracer_1d(apply_enhanced_ppm,&
                    nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                    velx(j,k,:),vely(j,k,:),velz(j,k,:), &
                    Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:), &
                    press(j,k,:))
            end do
         end do
         !$OMP END PARALLEL DO
      end if
    else
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if
  else if (CCTK_EQUALS(recon_method,"eno")) then

     if (evolve_tracer .ne. 0) then
        do itracer=1,number_of_tracers
           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                tracer(:,:,:,itracer), tracerplus(:,:,:,itracer), &
                tracerminus(:,:,:,itracer), trivial_rp, &
                hydro_excision_mask)
        enddo
     end if

     if (evolve_Y_e .ne. 0) then
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Y_e(:,:,:), Y_e_plus(:,:,:), &
             Y_e_minus(:,:,:), &
             trivial_rp, hydro_excision_mask)
     endif

    if (flux_direction == 1) then
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1
        do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 rho(:,j,k),rhominus(:,j,k),rhoplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 velx(:,j,k),velxminus(:,j,k),velxplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 vely(:,j,k),velyminus(:,j,k),velyplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 velz(:,j,k),velzminus(:,j,k),velzplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if(evolve_mhd.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    Bvecx(:,j,k),Bvecxminus(:,j,k),Bvecxplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    Bvecy(:,j,k),Bvecyminus(:,j,k),Bvecyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    Bvecz(:,j,k),Bveczminus(:,j,k),Bveczplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 dens(:,j,k),densminus(:,j,k),densplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 sx(:,j,k),sxminus(:,j,k),sxplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 sy(:,j,k),syminus(:,j,k),syplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 sz(:,j,k),szminus(:,j,k),szplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if(evolve_mhd.ne.0) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      Bconsx(:,j,k),Bconsxminus(:,j,k),Bconsxplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      Bconsy(:,j,k),Bconsyminus(:,j,k),Bconsyplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      Bconsz(:,j,k),Bconszminus(:,j,k),Bconszplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
          else
            !$OMP CRITICAL
            call CCTK_ERROR("Variable type to reconstruct not recognized.")
            STOP
            !$OMP END CRITICAL
          end if

          if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
             call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                  psidc(:,j,k),psidcminus(:,j,k),psidcplus(:,j,k),&
                  trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
          endif

          do i = 1, cctk_lsh(1)
            if (trivial_rp(i,j,k)) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, trivialx)
            else
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, not_trivialx)
            end if
          end do
        end do
       end do
        !$OMP END PARALLEL DO
     else if (flux_direction == 2) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1
           do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1
              if (CCTK_EQUALS(recon_vars,"primitive")) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      rho(j,:,k),rhominus(j,:,k),rhoplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      velx(j,:,k),velxminus(j,:,k),velxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      vely(j,:,k),velyminus(j,:,k),velyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      velz(j,:,k),velzminus(j,:,k),velzplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 if(evolve_mhd.ne.0) then
                    call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                         Bvecx(j,:,k),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                         Bvecy(j,:,k),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                         Bvecz(j,:,k),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 end if

              else if (CCTK_EQUALS(recon_vars,"conservative")) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      dens(j,:,k),densminus(j,:,k),densplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      sx(j,:,k),sxminus(j,:,k),sxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      sy(j,:,k),syminus(j,:,k),syplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      sz(j,:,k),szminus(j,:,k),szplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 if(evolve_mhd.ne.0) then
                    call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                         Bconsx(j,:,k),Bconsxminus(j,:,k),Bconsxplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                         Bconsy(j,:,k),Bconsyminus(j,:,k),Bconsyplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                         Bconsz(j,:,k),Bconszminus(j,:,k),Bconszplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 end if
              else
                 !$OMP CRITICAL
                 call CCTK_ERROR("Variable type to reconstruct not recognized.")
                 STOP
                 !$OMP END CRITICAL
              end if

              if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      psidc(j,:,k),psidcminus(j,:,k),psidcplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              endif

              do i = 1, cctk_lsh(2)
                 if (trivial_rp(j,i,k)) then
                    SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, trivialy)
                 else
                    SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, not_trivialy)
                 end if
              end do
           end do
        end do
      !$OMP END PARALLEL DO
    else if (flux_direction == 3) then
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1
        do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 rho(j,k,:),rhominus(j,k,:),rhoplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 velx(j,k,:),velxminus(j,k,:),velxplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 vely(j,k,:),velyminus(j,k,:),velyplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 velz(j,k,:),velzminus(j,k,:),velzplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if(evolve_mhd.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Bvecx(j,k,:),Bvecxminus(j,k,:),Bvecxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Bvecy(j,k,:),Bvecyminus(j,k,:),Bvecyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Bvecz(j,k,:),Bveczminus(j,k,:),Bveczplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 dens(j,k,:),densminus(j,k,:),densplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 sx(j,k,:),sxminus(j,k,:),sxplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 sy(j,k,:),syminus(j,k,:),syplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 sz(j,k,:),szminus(j,k,:),szplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if(evolve_mhd.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Bconsx(j,k,:),Bconsxminus(j,k,:),Bconsxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Bconsy(j,k,:),Bconsyminus(j,k,:),Bconsyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Bconsz(j,k,:),Bconszminus(j,k,:),Bconszplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
          else
            !$OMP CRITICAL
            call CCTK_ERROR("Variable type to reconstruct not recognized.")
            STOP
            !$OMP END CRITICAL
          end if

          if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
             call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                  psidc(j,k,:),psidcminus(j,k,:),psidcplus(j,k,:),&
                  trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
          endif

          do i = 1, cctk_lsh(3)
            if (trivial_rp(j,k,i)) then
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, trivialz)
            else
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, not_trivialz)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    else
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if
  else if (CCTK_EQUALS(recon_method,"weno")) then

     if (evolve_tracer .ne. 0) then
        do itracer=1,number_of_tracers
           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                tracer(:,:,:,itracer), tracerplus(:,:,:,itracer), &
                tracerminus(:,:,:,itracer), trivial_rp, &
                hydro_excision_mask)
        enddo
     end if

     if (evolve_Y_e .ne. 0) then
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Y_e(:,:,:), Y_e_plus(:,:,:), &
             Y_e_minus(:,:,:), &
             trivial_rp, hydro_excision_mask)
     endif

    if (flux_direction == 1) then
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1
        do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 rho(:,j,k),rhominus(:,j,k),rhoplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 velx(:,j,k),velxminus(:,j,k),velxplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 vely(:,j,k),velyminus(:,j,k),velyplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 velz(:,j,k),velzminus(:,j,k),velzplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if(evolve_mhd.ne.0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    Bvecx(:,j,k),Bvecxminus(:,j,k),Bvecxplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    Bvecy(:,j,k),Bvecyminus(:,j,k),Bvecyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    Bvecz(:,j,k),Bveczminus(:,j,k),Bveczplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 dens(:,j,k),densminus(:,j,k),densplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 sx(:,j,k),sxminus(:,j,k),sxplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 sy(:,j,k),syminus(:,j,k),syplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 sz(:,j,k),szminus(:,j,k),szplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if(evolve_mhd.ne.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      Bconsx(:,j,k),Bconsxminus(:,j,k),Bconsxplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      Bconsy(:,j,k),Bconsyminus(:,j,k),Bconsyplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      Bconsz(:,j,k),Bconszminus(:,j,k),Bconszplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
          else
            !$OMP CRITICAL
            call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
            !$OMP END CRITICAL
          end if

          if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
             call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                  psidc(:,j,k),psidcminus(:,j,k),psidcplus(:,j,k),&
                  trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
          endif

          do i = 1, cctk_lsh(1)
            if (trivial_rp(i,j,k)) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, trivialx)
            else
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, not_trivialx)
            end if
          end do
        end do
       end do
        !$OMP END PARALLEL DO
     else if (flux_direction == 2) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1
           do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1
              if (CCTK_EQUALS(recon_vars,"primitive")) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      rho(j,:,k),rhominus(j,:,k),rhoplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      velx(j,:,k),velxminus(j,:,k),velxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      vely(j,:,k),velyminus(j,:,k),velyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      velz(j,:,k),velzminus(j,:,k),velzplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 if(evolve_mhd.ne.0) then
                    call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                         Bvecx(j,:,k),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                         Bvecy(j,:,k),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                         Bvecz(j,:,k),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 end if

              else if (CCTK_EQUALS(recon_vars,"conservative")) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      dens(j,:,k),densminus(j,:,k),densplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      sx(j,:,k),sxminus(j,:,k),sxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      sy(j,:,k),syminus(j,:,k),syplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      sz(j,:,k),szminus(j,:,k),szplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 if(evolve_mhd.ne.0) then
                    call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                         Bconsx(j,:,k),Bconsxminus(j,:,k),Bconsxplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                         Bconsy(j,:,k),Bconsyminus(j,:,k),Bconsyplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                    call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                         Bconsz(j,:,k),Bconszminus(j,:,k),Bconszplus(j,:,k),&
                         trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 end if
              else
                 !$OMP CRITICAL
                 call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
                 !$OMP END CRITICAL
              end if

              if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      psidc(j,:,k),psidcminus(j,:,k),psidcplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              endif

              do i = 1, cctk_lsh(2)
                 if (trivial_rp(j,i,k)) then
                    SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, trivialy)
                 else
                    SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, not_trivialy)
                 end if
              end do
           end do
        end do
      !$OMP END PARALLEL DO
    else if (flux_direction == 3) then
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1
        do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 rho(j,k,:),rhominus(j,k,:),rhoplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 velx(j,k,:),velxminus(j,k,:),velxplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 vely(j,k,:),velyminus(j,k,:),velyplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 velz(j,k,:),velzminus(j,k,:),velzplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if(evolve_mhd.ne.0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Bvecx(j,k,:),Bvecxminus(j,k,:),Bvecxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Bvecy(j,k,:),Bvecyminus(j,k,:),Bvecyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Bvecz(j,k,:),Bveczminus(j,k,:),Bveczplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 dens(j,k,:),densminus(j,k,:),densplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 sx(j,k,:),sxminus(j,k,:),sxplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 sy(j,k,:),syminus(j,k,:),syplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 sz(j,k,:),szminus(j,k,:),szplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if(evolve_mhd.ne.0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Bconsx(j,k,:),Bconsxminus(j,k,:),Bconsxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Bconsy(j,k,:),Bconsyminus(j,k,:),Bconsyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Bconsz(j,k,:),Bconszminus(j,k,:),Bconszplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
          else
            !$OMP CRITICAL
            call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
            !$OMP END CRITICAL
          end if

          if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
             call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                  psidc(j,k,:),psidcminus(j,k,:),psidcplus(j,k,:),&
                  trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
          endif

          do i = 1, cctk_lsh(3)
            if (trivial_rp(j,k,i)) then
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, trivialz)
            else
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, not_trivialz)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    else
      call CCTK_WARN(0, "Flux direction not x,y,z")
    end if
  else
    call CCTK_ERROR("Reconstruction method not recognized!")
    STOP
  end if

  deallocate(trivial_rp)
  deallocate(lbetax, lbetay, lbetaz)
  deallocate(dum, dump, dumm)

  !$OMP PARALLEL WORKSHARE
  where ( (rhoplus < GRHydro_rho_min).or.(rhominus < GRHydro_rho_min).or.&
          (epsplus < 0.d0).or.(epsminus < 0.d0) )
    rhoplus = rho
    rhominus = rho
    velxplus = vel(:,:,:,1)
    velxminus = vel(:,:,:,1)
    velyplus = vel(:,:,:,2)
    velyminus = vel(:,:,:,2)
    velzplus = vel(:,:,:,3)
    velzminus = vel(:,:,:,3)
    epsplus = eps
    epsminus = eps
  end where
  !$OMP END PARALLEL WORKSHARE

  if (evolve_tracer .ne. 0) then
    if (use_min_tracer .ne. 0) then
      local_min_tracer = min_tracer
    else
      local_min_tracer = 0d0
    end if
   
    !$OMP PARALLEL WORKSHARE
    where( (tracerplus  .le. local_min_tracer).or.&
           (tracerminus .le. local_min_tracer) )
      tracerplus = tracer
      tracerminus = tracer
    end where
    !$OMP END PARALLEL WORKSHARE
    ! Call the conserved tracer routine in any case because (accord. to
    ! Christian Ott) this is the only way this works
    call Prim2ConservativeTracer(CCTK_PASS_FTOF)
  endif

  if (CCTK_EQUALS(recon_vars,"primitive").or.&
       CCTK_EQUALS(recon_method,"ppm")) then
     if(evolve_mhd.ne.0) then
        call Prim2ConservativePolytypeM(CCTK_PASS_FTOF)
     else
        call Prim2ConservativePolytype(CCTK_PASS_FTOF)
     endif
  else if (CCTK_EQUALS(recon_vars,"conservative")) then
     if(evolve_mhd.ne.0) then
        call Con2PrimBoundsPolytype(CCTK_PASS_FTOF)
     else
        call Con2PrimBoundsPolytype(CCTK_PASS_FTOF)
     endif
  else
     call CCTK_ERROR("Variable type to reconstruct not recognized.")
     STOP
  end if
  
  return
  
end subroutine ReconstructionPolytype

