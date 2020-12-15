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

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)
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

  integer :: nx, ny, nz, i, j, k

  logical, dimension(:,:,:), allocatable :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz

  CCTK_INT :: ierr

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: beta1, beta2, beta3
  pointer (pbeta1,beta1), (pbeta2,beta2), (pbeta3,beta3)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  logical :: apply_enhanced_ppm

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pbeta1 = loc(betaa)
    pbeta2 = loc(betab)
    pbeta3 = loc(betac)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pbeta1 = loc(betax)
    pbeta2 = loc(betay)
    pbeta3 = loc(betaz)
    pvup = loc(vel)
  end if

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

           trivial_rp(:,j,k) = .false.
           call PPMtyc(nx,CCTK_DELTA_SPACE(1),rho(:,j,k),velx(:,j,k),&
                vely(:,j,k),velz(:,j,k),eps(:,j,k),temperature(:,j,k),Y_e(:,j,k),&
                press(:,j,k),rhominus(:,j,k),&
                velxminus(:,j,k),velyminus(:,j,k),velzminus(:,j,k),&
                epsminus(:,j,k),tempminus(:,j,k), y_e_minus(:,j,k),&
                rhoplus(:,j,k),&
                velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),&
                epsplus(:,j,k),tempplus(:,j,k),y_e_plus(:,j,k))
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
      else if (flux_direction == 2) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, nz - GRHydro_stencil + 1
          do j = GRHydro_stencil, nx - GRHydro_stencil + 1
           trivial_rp(j,:,k) = .false.
           call PPMtyc(ny,CCTK_DELTA_SPACE(2),rho(j,:,k),&
                vely(j,:,k),velz(j,:,k),velx(j,:,k),&
                eps(j,:,k),temperature(j,:,k),Y_e(j,:,k), &
                press(j,:,k),rhominus(j,:,k),&
                velyminus(j,:,k),velzminus(j,:,k),velxminus(j,:,k),&
                epsminus(j,:,k),tempminus(j,:,k), y_e_minus(j,:,k),rhoplus(j,:,k),&
                velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),&
                epsplus(j,:,k),tempplus(j,:,k),y_e_plus(j,:,k))
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

    else if (flux_direction == 3) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = GRHydro_stencil, ny - GRHydro_stencil + 1
          do j = GRHydro_stencil, nx - GRHydro_stencil + 1
           trivial_rp(j,k,:) = .false.
           call PPMtyc(nz,CCTK_DELTA_SPACE(3),rho(j,k,:),&
                velz(j,k,:),velx(j,k,:),vely(j,k,:),&
                eps(j,k,:),temperature(j,k,:),Y_e(j,k,:), press(j,k,:),rhominus(j,k,:),&
                velzminus(j,k,:),velxminus(j,k,:),velyminus(j,k,:),&
                epsminus(j,k,:),tempminus(j,k,:),y_e_minus(j,k,:),rhoplus(j,k,:),&
                velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),&
                epsplus(j,k,:),tempplus(j,k,:),y_e_plus(j,k,:))
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
    else
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if
!!$ PPM ends.

  deallocate(trivial_rp)

end subroutine GRHydro_PPMReconstruct_drv

