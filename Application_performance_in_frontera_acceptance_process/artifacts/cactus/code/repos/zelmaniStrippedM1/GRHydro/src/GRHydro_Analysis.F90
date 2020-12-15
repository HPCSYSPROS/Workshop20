
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

subroutine GRHydro_Analysis_Init(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_REAL :: det, sdet

  divB = 0.0d0

  !$OMP PARALLEL DO PRIVATE(i,j,k,det,sdet)
  do k = 1, cctk_lsh(3) ! we need to compute Evec on all faces/edges where the fluxes are defined
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           sdet = sdetg(i,j,k)
           Bcons(i,j,k,1) = sdet * Bvec(i,j,k,1)
           Bcons(i,j,k,2) = sdet * Bvec(i,j,k,2)
           Bcons(i,j,k,3) = sdet * Bvec(i,j,k,3)
        end do
     end do
  end do
  !$OMP END PARALLEL DO


end subroutine GRHydro_Analysis_Init

subroutine GRHydro_CalcDivB(CCTK_ARGUMENTS)
    
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k,itracer
  CCTK_REAL :: idx, idy, idz, Bcons_l1, Bcons_r1, Bcons_l2, Bcons_r2, Bcons_l3, Bcons_r3

  idx = 1.d0 / CCTK_DELTA_SPACE(1)
  idy = 1.d0 / CCTK_DELTA_SPACE(2)
  idz = 1.d0 / CCTK_DELTA_SPACE(3)

      !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,Bcons_l1,Bcons_r1,Bcons_l2,Bcons_r2,Bcons_l3,Bcons_r3)
      do k = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(3) - GRHydro_stencil ! we need to compute Evec on all faces/edges where the fluxes are defined
        do j = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(2) - GRHydro_stencil
          do i = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(1) - GRHydro_stencil

           if(evolve_mhd.ne.0) then
               if(track_divB.ne.0) then
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
                 else
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

                 endif
               endif
            endif
    
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
 
end subroutine GRHydro_CalcDivB 
