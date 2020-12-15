#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



! Set the stress energy (storage) state grid scalar from the parameter.

subroutine TmunuBase_SetStressEnergyState (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  stress_energy_state = stress_energy_storage
  if (support_old_CalcTmunu_mechanism /= 0) then
     stress_energy_2_state = stress_energy_state
  else 
     stress_energy_2_state = 0
  end if
end subroutine TmunuBase_SetStressEnergyState
