/*@@
  @file      GRHydro_Shift.F90
  @date      Fri Mar 14 17:36:37 2003
  @author    Ian Hawke
  @desc 
  Set the shift so that it is approximately comoving with the
  velocity.
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydro_ComovingShift
   @date       Fri Mar 14 17:37:30 2003
   @author     Ian Hawke
   @desc 
   Implements the comoving shift
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_ComovingShift(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: ierr, Reduction_Handle, VarIndex
  CCTK_REAL :: rho_max
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: attenuate

  call CCTK_ReductionHandle(Reduction_Handle, "maximum")
  call CCTK_VarIndex(VarIndex, "HydroBase::rho")

  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, rho_max, 1, VarIndex)
  if (ierr .ne. 0) then
    call CCTK_ERROR("Failed to find rho_max")
    STOP
  end if

!!$  Hard wire in nuclear density

  if (CCTK_EQUALS(comoving_attenuate,"sqrt")) then
    attenuate = sqrt(rho / rho_max)
  else if (CCTK_EQUALS(comoving_attenuate,"tanh")) then
    attenuate = 0.5d0*(tanh(comoving_tanh_factor*&
         (rho/rho_max - comoving_tanh_offset)) + 1.d0 )
  else
    attenuate = 1.d0
  end if
  
  if (rho_max < 0.000324d0) then
    
    if (CCTK_EQUALS(comoving_v_method,"projected")) then
      
      where (rho > 10.d0 * GRHydro_rho_min)
        
        betax = -comoving_factor*x*alp*(x*vel(:,:,:,1) + y*vel(:,:,:,2) + z*vel(:,:,:,3)) / &
             (r**2 + 1.d-10) * attenuate
        betay = -comoving_factor*y*alp*(x*vel(:,:,:,1) + y*vel(:,:,:,2) + z*vel(:,:,:,3)) / &
             (r**2 + 1.d-10) * attenuate
        betaz = -comoving_factor*z*alp*(x*vel(:,:,:,1) + y*vel(:,:,:,2) + z*vel(:,:,:,3)) / &
             (r**2 + 1.d-10) * attenuate
        
      elsewhere
        
        betax = 0.d0
        betay = 0.d0
        betaz = 0.d0
        
      end where
      
    else if (CCTK_EQUALS(comoving_v_method,"components")) then
      
      where (rho > 10.d0 * GRHydro_rho_min)
        
        betax = comoving_factor * alp * x / (r + 1.d-10) * sqrt( (&
             (vel_p(:,:,:,1) - betax_p / alp_p)**2 + &
             (vel_p(:,:,:,2) - betay_p / alp_p)**2 + &
             (vel_p(:,:,:,3) - betaz_p / alp_p)**2) ) * attenuate
        betay = comoving_factor * alp * y / (r + 1.d-10) * sqrt( (&
             (vel_p(:,:,:,1) - betax_p / alp_p)**2 + &
             (vel_p(:,:,:,2) - betay_p / alp_p)**2 + &
             (vel_p(:,:,:,3) - betaz_p / alp_p)**2) ) * attenuate
        betaz = comoving_factor * alp * z / (r + 1.d-10) * sqrt( (&
             (vel_p(:,:,:,1) - betax_p / alp_p)**2 + &
             (vel_p(:,:,:,2) - betay_p / alp_p)**2 + &
             (vel_p(:,:,:,3) - betaz_p / alp_p)**2) ) * attenuate
        
      elsewhere
        
        betax = 0.d0
        betay = 0.d0
        betaz = 0.d0
        
      end where
      
    end if
    
  else
    
    betax = 0.d0
    betay = 0.d0
    betaz = 0.d0
    
  end if
  
end subroutine GRHydro_ComovingShift

 /*@@
   @routine    GRHydro_SetUpCoords
   @date       Fri Mar 14 19:11:13 2003
   @author     Ian Hawke
   @desc 
   Initializes the coordinates that are evolved with the
   comoving shift
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_SetUpCoords(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS

  GRHydro_x = x
  GRHydro_y = y
  GRHydro_z = z

end subroutine GRHydro_SetUpCoords

 /*@@
   @routine    GRHydro_EvolveCoords
   @date       Wed Mar 19 12:04:20 2003
   @author     Ian Hawke
   @desc 
   Evolve the coordinate system
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_EvolveCoords(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS

  GRHydro_x_rhs = -betax
  GRHydro_y_rhs = -betay
  GRHydro_z_rhs = -betaz

end subroutine GRHydro_EvolveCoords
