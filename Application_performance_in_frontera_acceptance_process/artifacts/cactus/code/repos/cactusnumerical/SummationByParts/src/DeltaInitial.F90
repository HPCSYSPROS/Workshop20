! Check if there are enough grid cells and ghost zones to compute the finite
! differences.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine SBP_DeltaInitial(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

!$OMP PARALLEL WORKSHARE
  sbp_dx = CCTK_DELTA_SPACE(1)*h_scaling(1)
  sbp_dy = CCTK_DELTA_SPACE(2)*h_scaling(2)
  sbp_dz = CCTK_DELTA_SPACE(3)*h_scaling(3)
!$OMP END PARALLEL WORKSHARE

end subroutine SBP_DeltaInitial
