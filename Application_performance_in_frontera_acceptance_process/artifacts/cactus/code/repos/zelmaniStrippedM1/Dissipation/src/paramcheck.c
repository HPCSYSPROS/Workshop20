/* $Header$ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void
dissipation_paramcheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int i;
  int want_horizon;
  
  if (extra_dissipation_in_horizons)
  {
    want_horizon = 0;
    for (i=0; i<100; ++i)
    {
      want_horizon = want_horizon | (horizon_number[i] >= 0);
    }
    
    if (want_horizon && ! CCTK_IsFunctionAliased ("HorizonRadiusInDirection"))
    {
      CCTK_PARAMWARN ("The aliased function \"HorizonRadiusInDirection\" must be defined when the parameter \"extra_dissipation_in_horizons\" is set and one of the sources is AHFinderDirect");
    }
  }
  if (respect_emask && !use_mask)
  {
    CCTK_PARAMWARN ("You can only respect the emask if it is used (use_mask)");
  }
}
