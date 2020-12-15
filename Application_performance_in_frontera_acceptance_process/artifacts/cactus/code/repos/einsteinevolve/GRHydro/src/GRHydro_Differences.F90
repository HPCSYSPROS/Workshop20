#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

 /*@@
   @routine    GRHydro_DiffRho
   @date       Fri Feb 18 10:02:31 2005
   @author     Ian Hawke
   @desc 
   Compute first differences in rho
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_DiffRho(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: idx, idy, idz

  CCTK_REAL :: dr_eps = 1.d-2

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

!!$  idx =  1.d0 / CCTK_DELTA_SPACE(1)
!!$  idy =  1.d0 / CCTK_DELTA_SPACE(2)
!!$  idz =  1.d0 / CCTK_DELTA_SPACE(3)

  idx = 1.d0
  idy = 1.d0
  idz = 1.d0

  DiffRho = 0.d0

  if (CCTK_EQUALS(gradient_method, "First diff")) then

!!$        Simple first differences

    do k = 2, nz - 1
      do j = 2, ny - 1
        do i = 2, nx - 1
          
          DiffRho(i,j,k) = 0.5d0 * ( &
                                   idx * abs(rho(i+1,j,k) - rho(i-1,j,k)) + &
                                   idy * abs(rho(i,j+1,k) - rho(i,j-1,k)) + &
                                   idz * abs(rho(i,j,k+1) - rho(i,j,k-1)) ) / &
                           rho(i,j,k)    
          
        end do
      end do
    end do
    
  else if (CCTK_EQUALS(gradient_method, "Curvature")) then

    do k = 2, nz - 1
      do j = 2, ny - 1
        do i = 2, nx - 1                 

!!$        Gradient as in Paramesh / Flash

          DiffRho(i,j,k) = &
               sqrt ( &
!!$             xx
             (abs(rho(i+1,j,k) - 2.d0 * rho(i,j,k) + rho(i-1,j,k)) / &
             ( abs(rho(i+1,j,k) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i-1,j,k)) + &
               dr_eps * (abs(rho(i+1,j,k)) + 2.d0 * abs(rho(i,j,k)) + &
                         abs(rho(i-1,j,k))) ))**2 + &
!!$                         xy
             (abs(rho(i+1,j+1,k) - rho(i+1,j-1,k) - &
                 rho(i-1,j+1,k) + rho(i-1,j-1,k)) / &
             ( abs(rho(i+1,j,k) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i-1,j,k)) + &
               dr_eps * (abs(rho(i+1,j+1,k)) + abs(rho(i+1,j-1,k)) + &
                         abs(rho(i-1,j+1,k)) + abs(rho(i-1,j-1,k))) ))**2 + &
!!$                         xz
             (abs(rho(i+1,j,k+1) - rho(i+1,j,k-1) - &
                 rho(i-1,j,k+1) + rho(i-1,j,k-1)) / &
             ( abs(rho(i+1,j,k) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i-1,j,k)) + &
               dr_eps * (abs(rho(i+1,j,k+1)) + abs(rho(i+1,j,k-1)) + &
                         abs(rho(i-1,j,k+1)) + abs(rho(i-1,j,k-1))) ))**2 + &
!!$             yx
             (abs(rho(i+1,j+1,k) - rho(i+1,j-1,k) - &
                 rho(i-1,j+1,k) + rho(i-1,j-1,k)) / &
             ( abs(rho(i,j+1,k) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i,j-1,k)) + &
               dr_eps * (abs(rho(i+1,j+1,k)) + abs(rho(i+1,j-1,k)) + &
                         abs(rho(i-1,j+1,k)) + abs(rho(i-1,j-1,k))) ))**2 + &
!!$                         yy
             (abs(rho(i,j+1,k) - 2.d0 * rho(i,j,k) + rho(i,j-1,k)) / &
             ( abs(rho(i,j+1,k) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i,j-1,k)) + &
               dr_eps * (abs(rho(i,j+1,k)) + 2.d0 * abs(rho(i,j,k)) + &
                         abs(rho(i,j-1,k))) ))**2 + &
!!$                         yz
             (abs(rho(i,j+1,k+1) - rho(i,j+1,k-1) - &
                 rho(i,j-1,k+1) + rho(i,j-1,k-1)) / &
             ( abs(rho(i,j+1,k) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i,j-1,k)) + &
               dr_eps * (abs(rho(i,j+1,k+1)) + abs(rho(i,j+1,k-1)) + &
                         abs(rho(i,j-1,k+1)) + abs(rho(i,j-1,k-1))) ))**2 + &
!!$             zx
             (abs(rho(i+1,j,k+1) - rho(i+1,j,k-1) - &
                 rho(i-1,j,k+1) + rho(i-1,j,k-1)) / &
             ( abs(rho(i,j,k+1) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i,j,k-1)) + &
               dr_eps * (abs(rho(i+1,j,k+1)) + abs(rho(i+1,j,k-1)) + &
                         abs(rho(i-1,j,k+1)) + abs(rho(i-1,j,k-1))) ))**2 + &
!!$                         zy
             (abs(rho(i,j+1,k+1) - rho(i,j-1,k+1) - &
                 rho(i,j+1,k-1) + rho(i,j-1,k-1)) / &
             ( abs(rho(i,j,k+1) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i,j,k-1)) + &
               dr_eps * (abs(rho(i,j+1,k+1)) + abs(rho(i,j-1,k+1)) + &
                         abs(rho(i,j+1,k-1)) + abs(rho(i,j-1,k-1))) ))**2 + &
!!$                         zz
             (abs(rho(i,j,k+1) - 2.d0 * rho(i,j,k) + rho(i,j,k-1)) / &
             ( abs(rho(i,j,k+1) - rho(i,j,k)) + &
               abs(rho(i,j,k) - rho(i,j,k-1)) + &
               dr_eps * (abs(rho(i,j,k+1)) + 2.d0 * abs(rho(i,j,k)) + &
                         abs(rho(i,j,k-1))) ))**2 &
          )
          
        end do
      end do
    end do
    
  else if (CCTK_EQUALS(gradient_method, "Rho weighted")) then

    DiffRho = &
         rho*(CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3))**2

  else
    call CCTK_ERROR("Value of gradient_method not recognized")
    STOP
  end if

end subroutine GRHydro_DiffRho
