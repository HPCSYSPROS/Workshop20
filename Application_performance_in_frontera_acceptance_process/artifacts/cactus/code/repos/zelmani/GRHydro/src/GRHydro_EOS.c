 /*@@
   @file      GRHydro_EOS.c
   @date      Wed Feb  6 18:25:33 2002
   @author    
   @desc 
   Sets the EOS handle number for use by
   all the GRHydro routines
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <stdlib.h>

 /*@@
   @routine    GRHydro_EOSHandle
   @date       Wed Feb  6 18:28:01 2002
   @author     Ian Hawke
   @desc 
   Sets the EOS handle number
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void GRHydro_EOSHandle(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Trying to get EOS handles");
  EOS_Omni_GetHandle("2D_Polytrope");
  CCTK_INFO("Trying to get EOS handles");
  if (!(*GRHydro_eos_handle       = EOS_Omni_GetHandle(GRHydro_eos_table)))
    CCTK_WARN(0, "Cannot get EOS handle, aborting.");
  if (!(*GRHydro_polytrope_handle = EOS_Omni_GetHandle("2D_Polytrope")))
    CCTK_WARN(0, "Cannot get polytrope EOS handle, aborting.");
  CCTK_VInfo(CCTK_THORNSTRING, "GRHydro will use the %s equation of state.",
                               GRHydro_eos_table);

}
