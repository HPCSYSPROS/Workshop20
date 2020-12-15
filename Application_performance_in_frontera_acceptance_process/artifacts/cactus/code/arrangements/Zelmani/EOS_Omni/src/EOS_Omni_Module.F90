module EOS_Omni_Module

  implicit none

  ! conversion factors between cgs and M_Sun = c = G = 1
  ! see EOS_Omni/doc/units.py
  real*8,parameter :: rho_gf = 1.61887093132742d-18
  real*8,parameter :: press_gf = 1.80123683248503d-39
  real*8,parameter :: eps_gf = 1.11265005605362d-21
  real*8,parameter :: time_gf = 2.03040204956746d05
  real*8,parameter :: mass_gf =  5.02916918125126d-34
  real*8,parameter :: length_gf = 6.77269222552442d-06
  
  ! Inverses of the numbers above, calculated manually instead of by
  ! the compiler
  
  real*8,parameter :: inv_rho_gf = 6.17714470405638d17
  real*8,parameter :: inv_press_gf = 5.55174079257738d38
  real*8,parameter :: inv_eps_gf = 8.98755178736818d20
  real*8,parameter :: inv_time_gf = 4.92513293223396d-6
  real*8,parameter :: inv_mass_gf = 1.98840000000000d33
  real*8,parameter :: inv_length_gf = 1.47651770773117d05

  real*8,parameter :: clite = 2.99792458d10
  real*8,parameter :: cliteinv2 = 1.11265005605362d-21

  ! These values are set by EOS_Omni_Startup
  real*8 :: hybrid_k2 = 0.0d0
  real*8 :: poly_k_cgs = 0.0d0
  real*8 :: gl_k_cgs   = 0.0d0
 
  real*8 :: hybrid_k1_cgs = 0.0d0
  real*8 :: hybrid_k2_cgs = 0.0d0
  
  real*8 :: pwp_k2 = 0.0d0
  real*8 :: pwp_k3 = 0.0d0
  real*8 :: pwp_k4 = 0.0d0
  
  real*8 :: pwp_a1 = 0.0d0
  real*8 :: pwp_a2 = 0.0d0
  real*8 :: pwp_a3 = 0.0d0
  real*8 :: pwp_a4 = 0.0d0

  ! stuff for the cold, tabulated EOS with a gamma law
  ! set by the reader routine
  integer :: coldeos_nrho = 0
  real*8 :: coldeos_gammath = 0.0d0
  real*8 :: coldeos_rhomin = 0.0d0
  real*8 :: coldeos_rhomax = 0.0d0
  real*8 :: coldeos_low_kappa = 0.0d0
  real*8 :: coldeos_low_gamma = 0.0d0
  real*8 :: coldeos_kappa = 0.0d0
  real*8 :: coldeos_thfac = 1.0d0
  real*8 :: coldeos_dlrho = 1.0d0
  real*8 :: coldeos_dlrhoi = 1.0d0
  real*8, allocatable :: coldeos_logrho(:)
  real*8, allocatable :: coldeos_eps(:)
  real*8, allocatable :: coldeos_gamma(:)
  real*8, allocatable :: coldeos_cs2(:)


  ! stuff for the barotropic tabulated EOS with a gamma law
  ! set by the reader routine
  integer :: barotropiceos_nrho = 0
  real*8 :: barotropiceos_rhomin = 0.0d0
  real*8 :: barotropiceos_rhomax = 0.0d0
  real*8 :: barotropiceos_low_kappa = 0.0d0
  real*8 :: barotropiceos_low_gamma = 0.0d0
  real*8 :: barotropiceos_kappa = 0.0d0
  real*8 :: barotropiceos_thfac = 1.0d0
  real*8 :: barotropiceos_dlrho = 1.0d0
  real*8 :: barotropiceos_dlrhoi = 1.0d0
  real*8 :: barotropiceos_energyshift = 0.0d0
  real*8, allocatable :: barotropiceos_logrho(:)
  real*8, allocatable :: barotropiceos_logpress(:)
  real*8, allocatable :: barotropiceos_logeps(:)
  real*8, allocatable :: barotropiceos_temp(:)
  real*8, allocatable :: barotropiceos_ye(:)

  ! actual poly_gamma_ini used when using different EOS for initial data and run
  real*8 :: poly_gamma_ini

  ! energy shift and other vars if using nuc_eos
  real*8 :: energy_shift
  real*8 :: eos_tempmin,eos_tempmax
  real*8 :: eos_yemin,eos_yemax

end module EOS_Omni_Module
