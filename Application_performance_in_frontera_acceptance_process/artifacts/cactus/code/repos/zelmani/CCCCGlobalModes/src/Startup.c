#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


int CCCCGlobalModes_Startup(void) 
{
  const char *banner = "CCCCGlobalModes: Performing Global Mode Analysis on Rho";
  
  CCTK_RegisterBanner(banner);

  return 0;
  
}
