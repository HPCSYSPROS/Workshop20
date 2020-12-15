 /*@@
   @file      GRHydro_ENOReconstruct_drv.F90
   @date      Tue Jul 19 13:22:03 EDT 2011
   @author    Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   Driver routine to perform the ENO reconstruction.
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
   @routine    GRHydro_ENOReconstruct_drv
   @date       Tue Jul 19 13:24:34 EDT 2011
   @author     Luca Baiotti, Ian Hawke, Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   A driver routine to do ENO reconstruction. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_ENOReconstruct_drv(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates

  integer :: nx, ny, nz, i, j, k, itracer

  logical, dimension(:,:,:), allocatable :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz

  CCTK_REAL, dimension(:,:,:),allocatable :: &
       &dum, dump, dumm

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

  if (ierr .ne. 0) then
    call CCTK_ERROR("Allocation problems with dum dump dumm")
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

!!$ ENO starts:
     if (evolve_tracer .ne. 0) then
        do itracer=1,number_of_tracers
           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                tracer(:,:,:,itracer), tracerplus(:,:,:,itracer), &
                tracerminus(:,:,:,itracer), trivial_rp, &
                hydro_excision_mask)
        enddo
     end if

    if (flux_direction == 1) then
      ! constraint transport needs to be able to average fluxes in the directions
      ! other that flux_direction, which in turn need the primitives on interfaces
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1 + transport_constraints
        do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1 + transport_constraints
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 rho(:,j,k),rhominus(:,j,k),rhoplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                   vel(:,j,k,1),velxminus(:,j,k),velxplus(:,j,k),&
                   trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                   vel(:,j,k,2),velyminus(:,j,k),velyplus(:,j,k),&
                   trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                   vel(:,j,k,3),velzminus(:,j,k),velzplus(:,j,k),&
                   trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            else
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                   lvel(:,j,k,1),velxminus(:,j,k),velxplus(:,j,k),&
                   trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                   lvel(:,j,k,2),velyminus(:,j,k),velyplus(:,j,k),&
                   trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                   lvel(:,j,k,3),velzminus(:,j,k),velzplus(:,j,k),&
                   trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            end if
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 eps(:,j,k),epsminus(:,j,k),epsplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if(evolve_mhd.ne.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      Bvec(:,j,k,1),Bvecxminus(:,j,k),Bvecxplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      Bvec(:,j,k,2),Bvecyminus(:,j,k),Bvecyplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      Bvec(:,j,k,3),Bveczminus(:,j,k),Bveczplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               else
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      lBvec(:,j,k,1),Bvecxminus(:,j,k),Bvecxplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      lBvec(:,j,k,2),Bvecyminus(:,j,k),Bvecyplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                      lBvec(:,j,k,3),Bveczminus(:,j,k),Bveczplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               end if
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    entropy(:,j,k),entropyminus(:,j,k),entropyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
            if (evolve_Y_e.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
            if (evolve_temper.ne.0.and.reconstruct_temper.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    temperature(:,j,k),tempminus(:,j,k),tempplus(:,j,k),&
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
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 tau(:,j,k),tauminus(:,j,k),tauplus(:,j,k),&
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
            if (evolve_entropy .ne. 0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(1),&
                    entropycons(:,j,k),entropyconsminus(:,j,k),entropyconsplus(:,j,k),&
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
      ! constraint transport needs to be able to average fluxes in the directions
      ! other that flux_direction, which in turn need the primitives on interfaces
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1 + transport_constraints
        do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1 + transport_constraints
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 rho(j,:,k),rhominus(j,:,k),rhoplus(j,:,k),&
                 trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                   vel(j,:,k,1),velxminus(j,:,k),velxplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                   vel(j,:,k,2),velyminus(j,:,k),velyplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                   vel(j,:,k,3),velzminus(j,:,k),velzplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            else
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                   lvel(j,:,k,1),velxminus(j,:,k),velxplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                   lvel(j,:,k,2),velyminus(j,:,k),velyplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                   lvel(j,:,k,3),velzminus(j,:,k),velzplus(j,:,k),&
                   trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            end if
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 eps(j,:,k),epsminus(j,:,k),epsplus(j,:,k),&
                 trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            if(evolve_mhd.ne.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      Bvec(j,:,k,1),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      Bvec(j,:,k,2),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      Bvec(j,:,k,3),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
               else
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      lBvec(j,:,k,1),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      lBvec(j,:,k,2),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                      lBvec(j,:,k,3),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
               end if
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                    entropy(j,:,k),entropyminus(j,:,k),entropyplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
            if (evolve_Y_e.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                    Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
            if (evolve_temper.ne.0.and.reconstruct_temper.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                    temperature(j,:,k),tempminus(j,:,k),tempplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
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
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 tau(j,:,k),tauminus(j,:,k),tauplus(j,:,k),&
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
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(2),&
                    entropycons(j,:,k),entropyconsminus(j,:,k),entropyconsplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
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
      ! constraint transport needs to be able to average fluxes in the directions
      ! other that flux_direction, which in turn need the primitives on interfaces
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1 + transport_constraints
        do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1 + transport_constraints
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 rho(j,k,:),rhominus(j,k,:),rhoplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                   vel(j,k,:,1),velxminus(j,k,:),velxplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                   vel(j,k,:,2),velyminus(j,k,:),velyplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                   vel(j,k,:,3),velzminus(j,k,:),velzplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            else
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                   lvel(j,k,:,1),velxminus(j,k,:),velxplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                   lvel(j,k,:,2),velyminus(j,k,:),velyplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
              call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                   lvel(j,k,:,3),velzminus(j,k,:),velzplus(j,k,:),&
                   trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            end if
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 eps(j,k,:),epsminus(j,k,:),epsplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if(evolve_mhd.ne.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                      Bvec(j,k,:,1),Bvecxminus(j,k,:),Bvecxplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                      Bvec(j,k,:,2),Bvecyminus(j,k,:),Bvecyplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                      Bvec(j,k,:,3),Bveczminus(j,k,:),Bveczplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               else
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                      lBvec(j,k,:,1),Bvecxminus(j,k,:),Bvecxplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                      lBvec(j,k,:,2),Bvecyminus(j,k,:),Bvecyplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                      lBvec(j,k,:,3),Bveczminus(j,k,:),Bveczplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               end if
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    entropy(j,k,:),entropyminus(j,k,:),entropyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
            if (evolve_Y_e.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
            if (evolve_temper.ne.0.and.reconstruct_temper.ne.0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    temperature(j,k,:),tempminus(j,k,:),tempplus(j,k,:),&
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
            call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 tau(j,k,:),tauminus(j,k,:),tauplus(j,k,:),&
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
            if (evolve_entropy .ne. 0) then
               call GRHydro_ENOReconstruct1d(eno_order,cctk_lsh(3),&
                    entropycons(j,k,:),entropyconsminus(j,k,:),entropyconsplus(j,k,:),&
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
!!$ ENO ends.

  deallocate(trivial_rp)
  deallocate(dum,dump,dumm)

end subroutine GRHydro_ENOReconstruct_drv

