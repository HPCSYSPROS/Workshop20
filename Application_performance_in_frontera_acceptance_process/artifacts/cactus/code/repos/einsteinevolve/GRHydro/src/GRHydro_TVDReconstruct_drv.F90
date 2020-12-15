 /*@@
   @file      GRHydro_TVDReconstruct_drv.F90
   @date      Tue Jul 19 13:22:03 EDT 2011
   @author    Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   Driver routine to perform the TVD reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)
#define Bvecx(i,j,k) Bprim(i,j,k,1)
#define Bvecy(i,j,k) Bprim(i,j,k,2)
#define Bvecz(i,j,k) Bprim(i,j,k,3)
#define Bconsx(i,j,k) Bcons(i,j,k,1)
#define Bconsy(i,j,k) Bcons(i,j,k,2)
#define Bconsz(i,j,k) Bcons(i,j,k,3)


 /*@@
   @routine    GRHydro_TVDReconstruct_drv
   @date       Tue Jul 19 13:24:34 EDT 2011
   @author     Luca Baiotti, Ian Hawke, Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   A driver routine to do TVD reconstruction. Currently just does
   TVD on the primitive variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_TVDReconstruct_drv(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    call TVDReconstruct_drv(gaa, gab, gac, gbb, gbc, gcc, lvel, lBvec)
  else
    call TVDReconstruct_drv(gxx, gxy, gxz, gyy, gyz, gzz, vel, Bvec)
  end if

contains
subroutine TVDReconstruct_drv(g11, g12, g13, g22, g23, g33, vup, Bprim)

  implicit none

  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim

  integer :: nx, ny, nz, i, j, k, itracer

  logical, dimension(:,:,:), allocatable :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz

  CCTK_INT :: ierr

  ! variables used when reconstruct_Wv == true
  CCTK_REAL :: inv_w                                 ! 1/Lorentz factor
  CCTK_REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: wvel ! vel*w_lorentz
  CCTK_REAL :: agxx,agxy,agxz,agyy,agyz,agzz         ! metric components

  allocate(trivial_rp(cctk_ash(1),cctk_ash(2),cctk_ash(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_ERROR("Allocation problems with trivial_rp")
    STOP
  end if

  !!$ reconstruct w^i = w_lorentz*v^i to ensure slower than light
  !!$ speeds?
  if (reconstruct_Wv .ne. 0) then
    ! all velocity like quantities are now w_lorentz*velocity. We will convert
    ! back to ordinary velocity at the end, using w_lorentz = sqrt(1 + g_{ij}
    ! w^i w^j).
    allocate(wvel(cctk_ash(1),cctk_ash(2),cctk_ash(3),3),STAT=ierr)
    if (ierr .ne. 0) then
      call CCTK_ERROR("Allocation problems with wvel")
      STOP
    end if
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=1,cctk_lsh(3)
      do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
          wvel(i,j,k,1) = vup(i,j,k,1) * w_lorentz(i,j,k)
          wvel(i,j,k,2) = vup(i,j,k,2) * w_lorentz(i,j,k)
          wvel(i,j,k,3) = vup(i,j,k,3) * w_lorentz(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
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

!!$ TVD starts:
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

      if (reconstruct_Wv.ne.0) then
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             wvel(:,:,:,1), velxplus, velxminus, trivial_rp, hydro_excision_mask)
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             wvel(:,:,:,2), velyplus, velyminus, trivial_rp, hydro_excision_mask)
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             wvel(:,:,:,3), velzplus, velzminus, trivial_rp, hydro_excision_mask)
        ! divide out the Loretnz factor obtained from w_lorentz =
        ! sqrt(1+g_{ij} w^i w^j) for both the
        ! plus and minus quantities this should by construction ensure
        ! that any Lorentz factor calculated from them later on is
        ! physical (ie. > 1.d0)
        ! constraint transport needs to be able to average fluxes in the directions
        ! other that flux_direction, which in turn need the primitives on interfaces
        !$OMP PARALLEL DO PRIVATE(i,j,k,inv_w,agxx,agxy,agxz,agyy,agyz,agzz)
        do k = GRHydro_stencil, nz - GRHydro_stencil + 1 + transport_constraints*(1-zoffset)
          do j = GRHydro_stencil, ny - GRHydro_stencil + 1 + transport_constraints*(1-yoffset)
            do i = GRHydro_stencil, nx - GRHydro_stencil + 1 + transport_constraints*(1-xoffset)
              agxx = 0.5d0*( g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset) )
              agxy = 0.5d0*( g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset) )
              agxz = 0.5d0*( g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset) )
              agyy = 0.5d0*( g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset) )
              agyz = 0.5d0*( g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset) )
              agzz = 0.5d0*( g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset) )
              inv_w = 1d0/sqrt( 1.d0 + agxx*velxminus(i,j,k)*velxminus(i,j,k) &
                                     + agyy*velyminus(i,j,k)*velyminus(i,j,k) &
                                     + agzz*velzminus(i,j,k)*velzminus(i,j,k) &
                                + 2.d0*agxy*velxminus(i,j,k)*velyminus(i,j,k) &
                                + 2.d0*agxz*velxminus(i,j,k)*velzminus(i,j,k) &
                                + 2.d0*agyz*velyminus(i,j,k)*velzminus(i,j,k) )
              velxminus(i,j,k) = velxminus(i,j,k)*inv_w
              velyminus(i,j,k) = velyminus(i,j,k)*inv_w
              velzminus(i,j,k) = velzminus(i,j,k)*inv_w
              
              agxx = 0.5d0*( g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset) )
              agxy = 0.5d0*( g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset) )
              agxz = 0.5d0*( g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset) )
              agyy = 0.5d0*( g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset) )
              agyz = 0.5d0*( g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset) )
              agzz = 0.5d0*( g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset) )
              inv_w = 1d0/sqrt( 1.d0 + agxx*velxplus(i,j,k)*velxplus(i,j,k) &
                                     + agyy*velyplus(i,j,k)*velyplus(i,j,k) &
                                     + agzz*velzplus(i,j,k)*velzplus(i,j,k) &
                                + 2.d0*agxy*velxplus(i,j,k)*velyplus(i,j,k) &
                                + 2.d0*agxz*velxplus(i,j,k)*velzplus(i,j,k) &
                                + 2.d0*agyz*velyplus(i,j,k)*velzplus(i,j,k) )
              velxplus(i,j,k) = velxplus(i,j,k)*inv_w
              velyplus(i,j,k) = velyplus(i,j,k)*inv_w
              velzplus(i,j,k) = velzplus(i,j,k)*inv_w
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      else
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             vup(:,:,:,1), velxplus, velxminus, trivial_rp, hydro_excision_mask)
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             vup(:,:,:,2), velyplus, velyminus, trivial_rp, hydro_excision_mask)
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             vup(:,:,:,3), velzplus, velzminus, trivial_rp, hydro_excision_mask)
      end if
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             eps, epsplus, epsminus, trivial_rp, hydro_excision_mask)
      if(evolve_mhd.ne.0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bprim(:,:,:,1), Bvecxplus, Bvecxminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bprim(:,:,:,2), Bvecyplus, Bvecyminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bprim(:,:,:,3), Bveczplus, Bveczminus, trivial_rp, hydro_excision_mask)
      endif
      if (evolve_entropy .ne. 0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              entropy, entropyplus, entropyminus, trivial_rp, hydro_excision_mask)
      endif
      if (evolve_temper .ne. 0 .and. reconstruct_temper .ne. 0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              temperature, tempplus, tempminus, trivial_rp, hydro_excision_mask)
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
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           tau, tauplus, tauminus, trivial_rp, hydro_excision_mask)
      if(evolve_mhd.ne.0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bcons(:,:,:,1), Bconsxplus, Bconsxminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bcons(:,:,:,2), Bconsyplus, Bconsyminus, trivial_rp, hydro_excision_mask)
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              Bcons(:,:,:,3), Bconszplus, Bconszminus, trivial_rp, hydro_excision_mask)
      endif
      if (evolve_entropy .ne. 0) then
         call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
              entropycons, entropyconsplus, entropyconsminus, trivial_rp, hydro_excision_mask)
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
!!$ TVD ends.

  deallocate(trivial_rp)
  ! Fortran 90 will automatically deallocate wvel for us

end subroutine
end subroutine GRHydro_TVDReconstruct_drv


