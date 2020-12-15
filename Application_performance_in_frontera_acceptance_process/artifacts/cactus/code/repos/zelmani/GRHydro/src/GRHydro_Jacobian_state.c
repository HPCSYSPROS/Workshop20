#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"




/**
   Checks states for Jacobians
*/
void
GRHydro_check_Jacobian_state (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  static CCTK_INT idxJacobian_state = -1, idxiJacobian_state = -1;

  if (!CCTK_IsImplementationActive("Coordinates"))
     return; /* Multipatch infrastructure not active */

  if (idxJacobian_state == -1)
  {
     idxJacobian_state = CCTK_VarIndex("Coordinates::jacobian_state");
     assert(idxJacobian_state >= 0);
  }

  if (idxiJacobian_state == -1)
  {
     idxiJacobian_state = CCTK_VarIndex("Coordinates::inverse_jacobian_state");
     assert(idxiJacobian_state >= 0);
  }

  if (*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxJacobian_state) == 0)
  {
     CCTK_WARN(0, "GRHydro requires storage for Jacobians. Please tell your Coordinates implementation to store the Jacobians!");
  }
  
  if (*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian_state) == 0)
  {
     CCTK_WARN(0, "GRHydro requires storage for inverse Jacobians. Please tell your Coordinates implementation to store the inverse Jacobians!");
  }
  
}



