/* $Header$ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void dissipation_basegrid (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int d;
  
  if (cctk_dim != 3) {
    CCTK_WARN (0, "This thorn supports only dim=3");
  }
  
  for (d=0; d<cctk_dim; ++d) {
    if (cctk_nghostzones[d] < (order+1)/2) {
      CCTK_WARN (0, "This thorn requires at least (order+1)/2 ghost zones");
    }
  }

  for (d=0;d<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];d++) {
    epsdisA[d]=epsdis;
    epsdisB[d]=epsdis2;
  }
}
