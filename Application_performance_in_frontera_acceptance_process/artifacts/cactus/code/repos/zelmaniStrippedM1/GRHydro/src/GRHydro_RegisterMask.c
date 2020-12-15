 /*@@
   @file      GRHydro_RegisterMask.c
   @date      Sun Jan 26 01:55:25 2003
   @author    Ian Hawke
   @desc 
   Routines to register states with SpaceMask.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "SpaceMask.h"

#include <stdio.h>
#include <stdlib.h>

 /*@@
   @routine    GRHydro_RegisterMask
   @date       Sun Jan 26 01:56:06 2003
   @author     
   @desc 
   Register the different mask states with the SpaceMask thorn.
   
   At the moment, the recognized states and values are

   Hydro_RiemannProblem (trivial, not_trivial)

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int GRHydro_RegisterMask(void)
{
  
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;

  const char *rp_list[2] = {"trivial","not_trivial"};
  
  ierr = SpaceMask_RegisterType("Hydro_RiemannProblemX", 2, rp_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the x Riemann Problem with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_RiemannProblemY", 2, rp_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the y Riemann Problem with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_RiemannProblemZ", 2, rp_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the z Riemann Problem with the mask!");
  }

  return 0;
}

