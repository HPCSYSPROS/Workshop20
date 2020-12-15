// GRHydro_LastMoLPostStep.c
//
// Compute is this is the last MoL PostStep call. Code taken from Christian
// Reisswig's rejected MoL changes.
//
// Roland Haas
// Sun Jun  3 17:35:53 PDT 2012

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void GRHydro_SetLastMoLPostStep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  const CCTK_INT *MoL_Intermediate_Step = 
    CCTK_VarDataPtr(cctkGH,0,"MoL::MoL_Intermediate_Step");
  if(NULL == MoL_Intermediate_Step)
  {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Could not get data pointer for MoL::MoL_Intermediate_Step");
  }

  // If counter becomes zero, the only thing left to do is to call PostStep!
  // NOTE: MoL's counter is uninitialised during initial data setup!
  *InLastMoLPostStep = cctk_iteration == 0 || *MoL_Intermediate_Step == 0;
}

void GRHydro_ClearLastMoLPostStep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  *InLastMoLPostStep = 0;
}
