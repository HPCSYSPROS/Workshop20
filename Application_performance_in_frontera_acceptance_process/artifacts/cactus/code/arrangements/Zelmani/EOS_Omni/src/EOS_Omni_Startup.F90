#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EOS_Omni_Startup(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none
  character*256 :: infoline

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  if(poly_gamma_initial .gt. 0d0) then
    poly_gamma_ini = poly_gamma_initial
  else
    poly_gamma_ini = poly_gamma
  end if

  poly_k_cgs = poly_k * rho_gf**poly_gamma_ini / press_gf

  gl_k_cgs   = gl_k * rho_gf**poly_gamma_ini / press_gf

  hybrid_k1_cgs = hybrid_k1 * rho_gf**poly_gamma_ini / press_gf

  hybrid_k2_cgs = hybrid_k1_cgs * &
       (hybrid_rho_nuc * inv_rho_gf)**(hybrid_gamma1-hybrid_gamma2)

  pwp_k2 = pwp_k1 * &
       pwp_rho_1**(pwp_gamma1-pwp_gamma2)

  pwp_k3 = pwp_k2 * &
       pwp_rho_2**(pwp_gamma2-pwp_gamma3)
  
  pwp_k4 = pwp_k3 * &
       pwp_rho_3**(pwp_gamma3-pwp_gamma4)

  pwp_a1 = 0.d0
  
  pwp_a2 = pwp_a1 + pwp_k1 * (pwp_rho_1)**(pwp_gamma1-1.0d0) / (pwp_gamma1 - 1.0d0)  &
                  - pwp_k2 * (pwp_rho_1)**(pwp_gamma2-1.0d0) / (pwp_gamma2 - 1.0d0)
  
  pwp_a3 = pwp_a2 + pwp_k2 * (pwp_rho_2)**(pwp_gamma2-1.0d0) / (pwp_gamma2 - 1.0d0)  &
                  - pwp_k3 * (pwp_rho_2)**(pwp_gamma3-1.0d0) / (pwp_gamma3 - 1.0d0)

  pwp_a4 = pwp_a3 + pwp_k3 * (pwp_rho_3)**(pwp_gamma3-1.0d0) / (pwp_gamma3 - 1.0d0)  &
                  - pwp_k4 * (pwp_rho_3)**(pwp_gamma4-1.0d0) / (pwp_gamma4 - 1.0d0)

end subroutine EOS_Omni_Startup

subroutine EOS_Omni_Get_Energy_Shift(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  call nuc_eos_c_get_energy_shift(energy_shift,eos_tempmin,eos_tempmax,&
       eos_yemin,eos_yemax)

end subroutine EOS_Omni_Get_Energy_Shift
