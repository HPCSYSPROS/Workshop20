/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_BSSN";
  CCTK_RegisterBanner(banner);
  return 0;
}
