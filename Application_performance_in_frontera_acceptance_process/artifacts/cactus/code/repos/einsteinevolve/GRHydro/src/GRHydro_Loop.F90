 /*@@
   @file      GRHydro_Loop.F90
   @date      Sat Jan 26 01:50:14 2002
   @author    
   @desc 
   Routines controlling loop counters and direction offsets
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydroStartLoop
   @date       Sat Jan 26 01:50:46 2002
   @author     Ian Hawke
   @desc 
   Set up the counters before the loop
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydroStartLoop(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  
  flux_direction = 3
  
  xoffset = 0
  yoffset = 0
  zoffset = 1

end subroutine GRHydroStartLoop

 /*@@
   @routine    Advance the counters whilst in the loop
   @date       Sat Jan 26 01:51:29 2002
   @author     Ian Hawke
   @desc 
   Just increments the counter and resets the directions.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydroAdvanceLoop(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  
  flux_direction = flux_direction-1
  
  if (flux_direction .eq. 2) then
    
    xoffset = 0
    yoffset = 1
    zoffset = 0
    
  else if (flux_direction .eq. 1) then
    
    xoffset = 1
    yoffset = 0
    zoffset = 0
    
  else
    
    xoffset = -1000000
    yoffset = -1000000
    zoffset = -1000000
    
  end if
  
end subroutine GRHydroAdvanceLoop
