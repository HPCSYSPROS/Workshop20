  /*@@
   @file      GRHydro_BvecfromAvec 
   @date      Aug 31, 2010
   @author    Tanja Bode
   @desc 
   Calculate B^i (at cell center) from Avec (at edges)
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define Avecx(i,j,k) Avec(i,j,k,1)
#define Avecy(i,j,k) Avec(i,j,k,2)
#define Avecz(i,j,k) Avec(i,j,k,3)

! Second order f.d. for volume-centered quantity from volume-centered quantity
#define DIFF_X_2(q) (0.5d0 * (q(i+1,j,k) - q(i-1,j,k)) * idx)
#define DIFF_Y_2(q) (0.5d0 * (q(i,j+1,k) - q(i,j-1,k)) * idy)
#define DIFF_Z_2(q) (0.5d0 * (q(i,j,k+1) - q(i,j,k-1)) * idz)

! Fourth order f.d. for volume-centered quantity from volume-centered quantity
#define DIFF_X_4(q) ((- q(i+2,j,k) + 8.d0 * q(i+1,j,k) \
                      - 8.d0 * q(i-1,j,k) + q(i-2,j,k)) / 12.d0 * idx)
#define DIFF_Y_4(q) ((- q(i,j+2,k) + 8.d0 * q(i,j+1,k) \
                      - 8.d0 * q(i,j-1,k) + q(i,j-2,k)) / 12.d0 * idy)
#define DIFF_Z_4(q) ((- q(i,j,k+2) + 8.d0 * q(i,j,k+1) \
                      - 8.d0 * q(i,j,k-1) + q(i,j,k-2)) / 12.d0 * idz)

 /*@@
   @routine   GRHydro_BvecfromAvec
   @date      Oct 24
   @author    Tanja Bode
   @desc 
   Calculate B^i (at cell center) from Avec (at edges)
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_BvecfromAvec(CCTK_ARGUMENTS)

  implicit none

  CCTK_REAL :: Az_y, Ay_z, Ax_z, Az_x, Ay_x, Ax_y
  CCTK_REAL :: sdet, isdet
  CCTK_REAL :: dx, dy, dz, idx, idy, idz
  CCTK_INT  :: i, j, k, nx, ny, nz, GRHydro_UseGeneralCoordinates
  CCTK_INT  :: local_spatial_order

  logical, allocatable, dimension (:,:,:) :: force_spatial_second_order

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
     call CCTK_ERROR("Bvec from Avec only in Cartesian at the moment.");
     STOP
  end if

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  idx = 1.d0/dx
  idy = 1.d0/dy
  idz = 1.d0/dz
 
  allocate (force_spatial_second_order(nx,ny,nz))
  force_spatial_second_order = .FALSE.
  
  if (sources_spatial_order > 2) then
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1 + GRHydro_stencil, nz - GRHydro_stencil
      do j = 1 + GRHydro_stencil, ny - GRHydro_stencil
        do i = 1 + GRHydro_stencil, nx - GRHydro_stencil
          if ((i < 3).or.(i > cctk_lsh(1) - 2).or. &
               (j < 3).or.(j > cctk_lsh(2) - 2).or. &
               (k < 3).or.(k > cctk_lsh(3) - 2) ) then
            force_spatial_second_order(i,j,k) = .TRUE.
          else if ( use_mask > 0 ) then
            if (minval(emask(i-2:i+2,j-2:j+2,k-2:k+2)) < 0.75d0) then
              force_spatial_second_order(i,j,k) = .TRUE.
            end if
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if
  
 
  call CCTK_WARN(1,"Bvec from Avec start."); 
  !$OMP PARALLEL DO PRIVATE( i, j, k, sdet, isdet, local_spatial_order, &
  !$OMP                      Az_y, Ay_z, Ax_z, Az_x, Ay_x, Ax_y)
  do k = GRHydro_stencil, cctk_lsh(3)-GRHydro_stencil+1
    do j = GRHydro_stencil, cctk_lsh(2)-GRHydro_stencil+1
      do i = GRHydro_stencil, cctk_lsh(1)-GRHydro_stencil+1

          local_spatial_order = sources_spatial_order
          if (force_spatial_second_order(i,j,k)) then
            local_spatial_order = 2
          end if

          if (local_spatial_order .eq. 2) then
             Az_y = DIFF_Y_2(Avecz)
             Ay_z = DIFF_Z_2(Avecy)
             Ax_z = DIFF_Z_2(Avecx)
             Az_x = DIFF_X_2(Avecz)
             Ay_x = DIFF_X_2(Avecy)
             Ax_y = DIFF_Y_2(Avecx)
          else
             Az_y = DIFF_Y_4(Avecz)
             Ay_z = DIFF_Z_4(Avecy)
             Ax_z = DIFF_Z_4(Avecx)
             Az_x = DIFF_X_4(Avecz)
             Ay_x = DIFF_X_4(Avecy)
             Ax_y = DIFF_Y_4(Avecx)
          end if

          sdet = sdetg(i,j,k)
          isdet = 1.d0/sdet

          Bvecx(i,j,k) = isdet*( Az_y - Ay_z )
          Bvecy(i,j,k) = isdet*( Ax_z - Az_x )
          Bvecz(i,j,k) = isdet*( Ay_x - Ax_y )

      end do
    end do
  end do
  !$OMP END PARALLEL DO
  call CCTK_WARN(1,"Bvec from Avec end."); 


end subroutine GRHydro_BvecfromAvec
