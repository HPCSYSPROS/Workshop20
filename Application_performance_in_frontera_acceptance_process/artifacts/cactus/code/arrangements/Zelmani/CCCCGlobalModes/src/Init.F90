#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"



subroutine CCCCGlobalModes_Init(CCTK_ARGUMENTS)
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(do_shibata.ne.0) then
    Ixx = 0.0d0;
    Iyy = 0.0d0;
    Ixy = 0.0d0;
    eta_plus = 0.0d0;
    eta_cross = 0.0d0;
  endif

  if(do_saijo.ne.0) then
    di_re = 0.0d0;
    di_im = 0.0d0;
    quad_re = 0.0d0;
    quad_im = 0.0d0;
    sextu_re = 0.0d0;
    sextu_im = 0.0d0;
    total_mass = 0.0d0;
  endif

  if(do_CoM.ne.0) then
    Mx = 0.0d0;
    My = 0.0d0;
    Mz = 0.0d0;
    Mr = 0.0d0;
  endif


  total_mass = 0.0d0

end subroutine CCCCGlobalModes_Init
