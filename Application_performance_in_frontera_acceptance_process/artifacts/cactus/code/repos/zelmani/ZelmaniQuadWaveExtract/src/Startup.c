#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


int ZelmaniQuadWaveExtract_Startup(void) 
{
  const char *banner = "ZelmaniQuadWaveExtract: Performing beautiful quadruople wave extraction";
  
  CCTK_RegisterBanner(banner);

  return 0;
  
}
