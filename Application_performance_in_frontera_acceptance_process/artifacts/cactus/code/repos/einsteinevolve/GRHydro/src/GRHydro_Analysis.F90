
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

subroutine GRHydro_CalcDivB(CCTK_ARGUMENTS)
    
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k,itracer
  CCTK_REAL :: idx, idy, idz, Bcons_l1, Bcons_r1, Bcons_l2, Bcons_r2, Bcons_l3, Bcons_r3
  logical :: do_Avec, do_Bcons

  ! for boundary conditions
  integer :: ierr
  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1

  idx = 1.d0 / CCTK_DELTA_SPACE(1)
  idy = 1.d0 / CCTK_DELTA_SPACE(2)
  idz = 1.d0 / CCTK_DELTA_SPACE(3)

  do_Avec = CCTK_EQUALS(Bvec_evolution_method,"GRHydro_Avec")
  do_Bcons = CCTK_EQUALS(Bvec_evolution_method,"GRHydro")

  divB = 0d0

  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,Bcons_l1,Bcons_r1,Bcons_l2,Bcons_r2,Bcons_l3,Bcons_r3)
  do k = 2, cctk_lsh(3)-1
    do j = 2, cctk_lsh(2)-1
      do i = 2, cctk_lsh(1)-1

        if(transport_constraints.ne.0) then
          ! edge based divergence (see WhiskyMHD & Bruno's thesis, Eq. 7.27)
          divB(i,j,k) = &
            0.25d0*(Bcons(i+1,j,k,1)-Bcons(i,j,k,1)+ &
                    Bcons(i+1,j+1,k,1)-Bcons(i,j+1,k,1)+ &
                    Bcons(i+1,j,k+1,1)-Bcons(i,j,k+1,1)+ &
                    Bcons(i+1,j+1,k+1,1)-Bcons(i,j+1,k+1,1))*idx + &
            0.25d0*(Bcons(i,j+1,k,2)-Bcons(i,j,k,2)+ &
                    Bcons(i+1,j+1,k,2)-Bcons(i+1,j,k,2)+ &
                    Bcons(i,j+1,k+1,2)-Bcons(i,j,k+1,2)+ &
                    Bcons(i+1,j+1,k+1,2)-Bcons(i+1,j,k+1,2))*idy + &
            0.25d0*(Bcons(i,j,k+1,3)-Bcons(i,j,k,3)+ &
                    Bcons(i+1,j,k+1,3)-Bcons(i+1,j,k,3)+ &
                    Bcons(i,j+1,k+1,3)-Bcons(i,j+1,k,3)+ &
                    Bcons(i+1,j+1,k+1,3)-Bcons(i+1,j+1,k,3))*idz
        else if(do_Avec) then
          Bcons_l1 = 0.5d0 * (sdetg(i,j,k)*Bvec(i,j,k,1) + &
                   sdetg(i-1,j,k)*Bvec(i-1,j,k,1))
          Bcons_l2 = 0.5d0 * (sdetg(i,j,k)*Bvec(i,j,k,2) + &
                   sdetg(i,j-1,k)*Bvec(i,j-1,k,2))
          Bcons_l3 = 0.5d0 * (sdetg(i,j,k)*Bvec(i,j,k,3) + &
                   sdetg(i,j,k-1)*Bvec(i,j,k-1,3))
          Bcons_r1 = 0.5d0 * (sdetg(i,j,k)*Bvec(i,j,k,1) + &
                   sdetg(i+1,j,k)*Bvec(i+1,j,k,1))
          Bcons_r2 = 0.5d0 * (sdetg(i,j,k)*Bvec(i,j,k,2) + &
                   sdetg(i,j+1,k)*Bvec(i,j+1,k,2))
          Bcons_r3 = 0.5d0 * (sdetg(i,j,k)*Bvec(i,j,k,3) + &
                   sdetg(i,j,k+1)*Bvec(i,j,k+1,3))

          divB(i,j,k) = (Bcons_l1 - Bcons_r1) * idx + (Bcons_l2 - Bcons_r2) * idy + (Bcons_l3 - Bcons_r3) * idz

        else if(do_Bcons) then
          Bcons_l1 = 0.5d0 * (Bcons(i,j,k,1) + &
                   Bcons(i-1,j,k,1))
          Bcons_l2 = 0.5d0 * (Bcons(i,j,k,2) + &
                   Bcons(i,j-1,k,2))
          Bcons_l3 = 0.5d0 * (Bcons(i,j,k,3) + &
                   Bcons(i,j,k-1,3))
          Bcons_r1 = 0.5d0 * (Bcons(i,j,k,1) + &
                   Bcons(i+1,j,k,1))
          Bcons_r2 = 0.5d0 * (Bcons(i,j,k,2) + &
                   Bcons(i,j+1,k,2))
          Bcons_r3 = 0.5d0 * (Bcons(i,j,k,3) + &
                   Bcons(i,j,k+1,3))

          divB(i,j,k) = (Bcons_l1 - Bcons_r1) * idx + (Bcons_l2 - Bcons_r2) * idy + (Bcons_l3 - Bcons_r3) * idz

        else
          call CCTK_ERROR("Internal error. Do now know how to compute divB")
          STOP

        endif

      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
 
!!$Flat boundaries if required
  if (CCTK_EQUALS(bound,"flat")) then
    ierr = Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::divB", "Flat")
  endif

!!$None boundaries if required
  if (CCTK_EQUALS(bound,"none")) then
    ierr = Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::divB", "None")
  endif

!!$Scalar boundaries if required
  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_ERROR("Until somebody uses this I see no reason to support it")
    STOP
  end if

  if (ierr < 0) then
    call CCTK_ERROR("problems with applying the chosen boundary condition")
    STOP
  end if

end subroutine GRHydro_CalcDivB 

