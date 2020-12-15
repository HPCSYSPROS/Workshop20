#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void CoreCollapseControl_ParamCheck(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals("corecollapsecontrol::bounce_criterion", "entropy") &&
     CCTK_ActiveTimeLevels(cctkGH, "HydroBase::entropy") <= 0) {
    CCTK_PARAMWARN("You must not set corecollapsecontrol::bounce_criterion = entropy if not actually using a hot EOS");
  }

}
