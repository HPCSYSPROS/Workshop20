 /*@@
   @file      GRHydro_RoeAverage.F90
   @date      Sat Jan 26 01:54:09 2002
   @author    
   @desc 
   Calculates the Roe average of two states.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine roeaverage    
   @date       Sat Jan 26 01:54:46 2002
   @author     Ian Hawke
   @desc 
   Calculate the Roe average of two states. Currently just 
   arithmetic averaging.
   @enddesc 
   @calls     
   @calledby   
   @history
   
   @endhistory 

@@*/

subroutine roeaverage(l, r, ave)
  
  implicit none
  
  CCTK_REAL, dimension(5) :: l, r, ave
  
  ave = 0.5d0 * (l + r)
  
end subroutine roeaverage
  
