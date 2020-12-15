 /*@@
   @file      GRHydro_DivergenceClean.F90
   @date      Nov 24, 2010
   @author    
   @desc 
   Routines controlling divergence cleaning
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydro_InitDivergenceClean
   @date       Nov 24, 2010
   @author     Joshua Faber
   @desc 
   Set Psi=0 initially for use with divergence cleaning
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_InitDivergenceClean(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  psidc=0.0

end subroutine GRHydro_InitDivergenceClean

  

