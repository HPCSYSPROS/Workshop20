 /*@@
   @file      GRHydro_WENOReconstruct_drv.F90
   @date      Fri Jan 3 2013
   @author    Christian Reisswig, Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   Driver routine to perform the WENO reconstruction.
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


!!$ transform back to vel if Wv was used
subroutine undo_Wv(nx,velxminus,velyminus,velzminus,&
                      velxplus,velyplus,velzplus,&
                      gxx,gxy,gxz,gyy,gyz,gzz)

   implicit none

   DECLARE_CCTK_PARAMETERS
   DECLARE_CCTK_FUNCTIONS
   
   CCTK_INT :: i,nx
   CCTK_REAL :: w, agxx, agxy, agxz, agyy, agyz, agzz
   CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz
   CCTK_REAL, dimension(nx) :: velxminus,velyminus,velzminus
   CCTK_REAL, dimension(nx) :: velxplus,velyplus,velzplus
                      
   do i = grhydro_stencil, nx - grhydro_stencil + 1
      ! divide out the Loretnz factor obtained from w_lorentz =
      ! sqrt(1+g_{ij} w^i w^j) for both the
      ! plus and minus quantities this should by construction ensure
      ! that any Lorentz factor calculated from them later on is
      ! physical (ie. > 1.d0)
      agxx = 0.5d0*( gxx(i) + gxx(i-1) )
      agxy = 0.5d0*( gxy(i) + gxy(i-1) )
      agxz = 0.5d0*( gxz(i) + gxz(i-1) )
      agyy = 0.5d0*( gyy(i) + gyy(i-1) )
      agyz = 0.5d0*( gyz(i) + gyz(i-1) )
      agzz = 0.5d0*( gzz(i) + gzz(i-1) )
      w = sqrt( 1.d0 + agxx*velxminus(i)*velxminus(i) + &
                        agyy*velyminus(i)*velyminus(i) &
               + agzz*velzminus(i)*velzminus(i) + &
                  2.d0*agxy*velxminus(i)*velyminus(i) &
               + 2.d0*agxz*velxminus(i)*velzminus(i) + &
                  2.d0*agyz*velyminus(i)*velzminus(i) )
      velxminus(i) = velxminus(i)/w
      velyminus(i) = velyminus(i)/w
      velzminus(i) = velzminus(i)/w
      
      agxx = 0.5d0*( gxx(i) + gxx(i+1) )
      agxy = 0.5d0*( gxy(i) + gxy(i+1) )
      agxz = 0.5d0*( gxz(i) + gxz(i+1) )
      agyy = 0.5d0*( gyy(i) + gyy(i+1) )
      agyz = 0.5d0*( gyz(i) + gyz(i+1) )
      agzz = 0.5d0*( gzz(i) + gzz(i+1) )
      w = sqrt( 1.d0 + agxx*velxplus(i)*velxplus(i) + &
                  agyy*velyplus(i)*velyplus(i) &
               + agzz*velzplus(i)*velzplus(i) + &
               2.d0*agxy*velxplus(i)*velyplus(i) &
               + 2.d0*agxz*velxplus(i)*velzplus(i) + &
               2.d0*agyz*velyplus(i)*velzplus(i) )
      velxplus(i) = velxplus(i)/w
      velyplus(i) = velyplus(i)/w
      velzplus(i) = velzplus(i)/w
   end do
end subroutine undo_Wv

  

 /*@@
   @routine    GRHydro_WENOReconstruct_drv
   @date       Fri Jan 3 2013
   @author     Christian Reisswig, Luca Baiotti, Ian Hawke, Bruno C. Mundim, Joshua Faber, Christian D. Ott
   @desc 
   A driver routine to do WENO reconstruction. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_WENOReconstruct_drv(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

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

  ! 1D arrays in x, y and z direction for reconstructing Wv instead of v
!   CCTK_REAL, dimension(nx) :: vxX,vyX,vzX
!   CCTK_REAL, dimension(ny) :: vxY,vyY,vzY
!   CCTK_REAL, dimension(nz) :: vxZ,vyZ,vzZ
  
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

!!$ WENO starts:
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
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 rho(:,j,k),rhominus(:,j,k),rhoplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if (reconstruct_Wv.eq.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    vel(:,j,k,1),velxminus(:,j,k),velxplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    vel(:,j,k,2),velyminus(:,j,k),velyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    vel(:,j,k,3),velzminus(:,j,k),velzplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    lvel(:,j,k,1),velxminus(:,j,k),velxplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    lvel(:,j,k,2),velyminus(:,j,k),velyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    lvel(:,j,k,3),velzminus(:,j,k),velzplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               end if
            else
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    w_lorentz(:,j,k)*vel(:,j,k,1),velxminus(:,j,k),velxplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    w_lorentz(:,j,k)*vel(:,j,k,2),velyminus(:,j,k),velyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    w_lorentz(:,j,k)*vel(:,j,k,3),velzminus(:,j,k),velzplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call undo_Wv(int(nx,ik), velxminus(:,j,k), velyminus(:,j,k), velzminus(:,j,k),&
                                  velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                                  gxx(:,j,k),gxy(:,j,k),gxz(:,j,k),gyy(:,j,k),gyz(:,j,k),gzz(:,j,k))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    w_lorentz(:,j,k)*lvel(:,j,k,1),velxminus(:,j,k),velxplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    w_lorentz(:,j,k)*lvel(:,j,k,2),velyminus(:,j,k),velyplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    w_lorentz(:,j,k)*lvel(:,j,k,3),velzminus(:,j,k),velzplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call undo_Wv(int(nx,ik), velxminus(:,j,k), velyminus(:,j,k), velzminus(:,j,k),&
                                  velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                                  gaa(:,j,k),gab(:,j,k),gac(:,j,k),gbb(:,j,k),gbc(:,j,k),gcc(:,j,k))
               end if
            endif
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 eps(:,j,k),epsminus(:,j,k),epsplus(:,j,k),&
                 trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            if(evolve_Y_e.ne.0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
            if(evolve_temper.ne.0.and.reconstruct_temper.ne.0) then 
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    temperature(:,j,k),tempminus(:,j,k),tempplus(:,j,k),&
                    trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
            endif
            if(evolve_mhd.ne.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      Bvec(:,j,k,1),Bvecxminus(:,j,k),Bvecxplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      Bvec(:,j,k,2),Bvecyminus(:,j,k),Bvecyplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      Bvec(:,j,k,3),Bveczminus(:,j,k),Bveczplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      lBvec(:,j,k,1),Bvecxminus(:,j,k),Bvecxplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      lBvec(:,j,k,2),Bvecyminus(:,j,k),Bvecyplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                      lBvec(:,j,k,3),Bveczminus(:,j,k),Bveczplus(:,j,k),&
                      trivial_rp(:,j,k), hydro_excision_mask(:,j,k))
               end if
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                    entropy(:,j,k),entropyminus(:,j,k),entropyplus(:,j,k),&
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
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
                 tau(:,j,k),tauminus(:,j,k),tauplus(:,j,k),&
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
            if (evolve_entropy .ne. 0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(1),&
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
      ! constraint transport needs to be able to average fluxes in the directions
      ! other that flux_direction, which in turn need the primitives on interfaces
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1 + transport_constraints
        do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1 + transport_constraints
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                 rho(j,:,k),rhominus(j,:,k),rhoplus(j,:,k),&
                 trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            if (reconstruct_Wv.eq.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    vel(j,:,k,1),velxminus(j,:,k),velxplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    vel(j,:,k,2),velyminus(j,:,k),velyplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    vel(j,:,k,3),velzminus(j,:,k),velzplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    lvel(j,:,k,1),velxminus(j,:,k),velxplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    lvel(j,:,k,2),velyminus(j,:,k),velyplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    lvel(j,:,k,3),velzminus(j,:,k),velzplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
               end if
            else
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    w_lorentz(j,:,k)*vel(j,:,k,1),velxminus(j,:,k),velxplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    w_lorentz(j,:,k)*vel(j,:,k,2),velyminus(j,:,k),velyplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    w_lorentz(j,:,k)*vel(j,:,k,3),velzminus(j,:,k),velzplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call undo_Wv(int(ny,ik), velyminus(j,:,k), velzminus(j,:,k), velxminus(j,:,k),&
                                  velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                                  gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), gxx(j,:,k))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    w_lorentz(j,:,k)*lvel(j,:,k,1),velxminus(j,:,k),velxplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    w_lorentz(j,:,k)*lvel(j,:,k,2),velyminus(j,:,k),velyplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    w_lorentz(j,:,k)*lvel(j,:,k,3),velzminus(j,:,k),velzplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call undo_Wv(int(ny,ik), velyminus(j,:,k), velzminus(j,:,k), velxminus(j,:,k),&
                                  velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                                  gbb(j,:,k), gbc(j,:,k), gab(j,:,k), gcc(j,:,k), gac(j,:,k), gaa(j,:,k))
               end if
            endif
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                 eps(j,:,k),epsminus(j,:,k),epsplus(j,:,k),&
                 trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            if(evolve_Y_e.ne.0) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                 Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k),&
                 trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
            if(evolve_temper.ne.0.and.reconstruct_temper.ne.0) then 
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    temperature(j,:,k),tempminus(j,:,k),tempplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
            if(evolve_mhd.ne.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      Bvec(j,:,k,1),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      Bvec(j,:,k,2),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      Bvec(j,:,k,3),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      lBvec(j,:,k,1),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      lBvec(j,:,k,2),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      lBvec(j,:,k,3),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
              end if
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    entropy(j,:,k),entropyminus(j,:,k),entropyplus(j,:,k),&
                    trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
            endif
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
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                 tau(j,:,k),tauminus(j,:,k),tauplus(j,:,k),&
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
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
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
      ! constraint transport needs to be able to average fluxes in the directions
      ! other that flux_direction, which in turn need the primitives on interfaces
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1 + transport_constraints
        do j = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1 + transport_constraints
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 rho(j,k,:),rhominus(j,k,:),rhoplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if (reconstruct_Wv.eq.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    vel(j,k,:,1),velxminus(j,k,:),velxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    vel(j,k,:,2),velyminus(j,k,:),velyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    vel(j,k,:,3),velzminus(j,k,:),velzplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      lBvec(j,:,k,1),Bvecxminus(j,:,k),Bvecxplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      lBvec(j,:,k,2),Bvecyminus(j,:,k),Bvecyplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                      lBvec(j,:,k,3),Bveczminus(j,:,k),Bveczplus(j,:,k),&
                      trivial_rp(j,:,k), hydro_excision_mask(j,:,k))
               end if
            else
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*vel(j,k,:,1),velxminus(j,k,:),velxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*vel(j,k,:,2),velyminus(j,k,:),velyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*vel(j,k,:,3),velzminus(j,k,:),velzplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call undo_Wv(int(nz,ik), velzminus(j,k,:), velxminus(j,k,:), velyminus(j,k,:),&
                                  velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                                  gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:),gyy(j,k,:))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*lvel(j,k,:,1),velxminus(j,k,:),velxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*lvel(j,k,:,2),velyminus(j,k,:),velyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*lvel(j,k,:,3),velzminus(j,k,:),velzplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call undo_Wv(int(nz,ik), velzminus(j,k,:), velxminus(j,k,:), velyminus(j,k,:),&
                                  velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                                  gcc(j,k,:), gac(j,k,:), gbc(j,k,:), gaa(j,k,:), gab(j,k,:),gbb(j,k,:))
               end if
            endif
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 eps(j,k,:),epsminus(j,k,:),epsplus(j,k,:),&
                 trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            if(evolve_Y_e.ne.0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
            if(evolve_temper.ne.0.and.reconstruct_temper.ne.0) then 
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    temperature(j,k,:),tempminus(j,k,:),tempplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
            endif
            if(evolve_mhd.ne.0) then
               if (GRHydro_UseGeneralCoordinates(cctkGH).eq.0) then
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                      Bvec(j,k,:,1),Bvecxminus(j,k,:),Bvecxplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                      Bvec(j,k,:,2),Bvecyminus(j,k,:),Bvecyplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                      Bvec(j,k,:,3),Bveczminus(j,k,:),Bveczplus(j,k,:),&
                      trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               else
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*lvel(j,k,:,1),velxminus(j,k,:),velxplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*lvel(j,k,:,2),velyminus(j,k,:),velyplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
                 call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                    w_lorentz(j,k,:)*lvel(j,k,:,3),velzminus(j,k,:),velzplus(j,k,:),&
                    trivial_rp(j,k,:), hydro_excision_mask(j,k,:))
               end if
            endif
            if (evolve_entropy .ne. 0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(2),&
                    entropy(j,k,:),entropyminus(j,k,:),entropyplus(j,k,:),&
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
            call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
                 tau(j,k,:),tauminus(j,k,:),tauplus(j,k,:),&
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
            if (evolve_entropy .ne. 0) then
               call GRHydro_WENOReconstruct1d(WENO_order,cctk_lsh(3),&
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
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if
!!$ WENO ends.


  deallocate(trivial_rp)
  deallocate(dum,dump,dumm)


end subroutine GRHydro_WENOReconstruct_drv

