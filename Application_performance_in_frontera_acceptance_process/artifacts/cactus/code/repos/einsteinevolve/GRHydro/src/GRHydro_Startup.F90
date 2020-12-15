 /*@@
   @file      GRHydro_Startup.F90
   @date      Sun Feb 10 00:02:52 2002
   @author    Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

 /*@@
   @routine    GRHydro_Startup
   @date       Sun Feb 10 00:03:09 2002
   @author     Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

integer function GRHydro_Startup()

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: ierr

  call CCTK_RegisterBanner(ierr, "GRHydro: relativistic hydrodynamics, no ice.")
  
  GRHydro_Startup = 0

end function GRHydro_Startup
