#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! eoskey:
! 1 --- polytropic EOS
! 2 --- gamma-law EOS
! 3 --- hybrid EOS
! 4 --- finite-T microphysical NSE EOS

subroutine EOS_Omni_EOS_Press(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: press(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           press(i) = press_gf * poly_k_cgs * &
                (rho(i)*inv_rho_gf)**poly_gamma
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           press(i) = (gl_gamma - 1.0d0) * rho(i) * eps(i)
        enddo

     case (3)
        ! hybrid EOS
        do i=1,npoints
           if(rho(i).gt.hybrid_rho_nuc) then
              hybrid_local_gamma = hybrid_gamma2
              hybrid_local_k_cgs = hybrid_k2_cgs
           else
              hybrid_local_gamma = hybrid_gamma1
              hybrid_local_k_cgs = hybrid_k1_cgs
           endif
           hybrid_p_poly = press_gf * hybrid_local_k_cgs * &
                     (rho(i) * inv_rho_gf)**hybrid_local_gamma
           hybrid_p_th = - press_gf * hybrid_local_k_cgs * (hybrid_gamma_th - 1.d0) /      &
                (hybrid_local_gamma - 1.0d0) * (rho(i) * inv_rho_gf)**hybrid_local_gamma + &
                (hybrid_gamma_th - 1.0d0) * rho(i) * eps(i) -                              &
                (hybrid_gamma_th - 1.d0) * (hybrid_local_gamma - hybrid_gamma1) /          &
                (hybrid_gamma1 - 1.d0) / (hybrid_gamma2 - 1.d0) *                          &
                press_gf * hybrid_k1_cgs * inv_rho_gf**hybrid_gamma1 *                     &
                hybrid_rho_nuc**(hybrid_gamma1 - 1.d0) * rho(i)
           hybrid_p_th = max(zero, hybrid_p_th)
           press(i) = hybrid_p_poly + hybrid_p_th
        enddo

     case (4)
        ! nuc eos
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps(npoints,rho,temp,ye,&
                eps,press,keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_press(npoints,rho,temp,ye,eps,press,&
                rf_precision,keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not suppported!")
           STOP
        endif
     case (5)
        ! cold tabular EOS with gamma law
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              press(i) = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = press(i) / (coldeos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho = log10(rho(i))
              ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th
        enddo

     case (6)
        ! barotropic tabular EOS with gamma law
        do i=1,npoints
           if(rho(i).lt.barotropiceos_rhomin) then
              press(i) = barotropiceos_low_kappa * rho(i)**barotropiceos_low_gamma
              eps(i) = press(i) / (barotropiceos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.barotropiceos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho = log10(rho(i))
              ir = 2 + &
                   INT( (xrho - barotropiceos_logrho(1) - 1.0d-10) &
                   * barotropiceos_dlrhoi )
           endif
           ir = max(2, min(ir,barotropiceos_nrho))

           xprs = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

           xenr = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)

           xtemp = (barotropiceos_temp(ir) - barotropiceos_temp(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_temp(ir-1)

           xye = (barotropiceos_ye(ir) - barotropiceos_ye(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_ye(ir-1)

           press_cold = 10.0**xprs
           eps_cold = 10.0**xenr - barotropiceos_energyshift

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif
           press_th = coldeos_thfac*(barotropiceos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_Press

subroutine EOS_Omni_EOS_PressOMP(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
#ifdef _OPENMP
  use OMP_LIB
#endif

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: press(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma,  hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma
  integer :: num_threads, my_thread_num, slice_len, slice_min, slice_max
  CCTK_INT :: my_anyerr

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           !$OMP PARALLEL DO
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
           !$OMP END PARALLEL DO
        endif
        !$OMP PARALLEL DO
        do i=1,npoints
           press(i) = press_gf * poly_k_cgs * &
                (rho(i)*inv_rho_gf)**poly_gamma
        enddo
        !$OMP END PARALLEL DO
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           !$OMP PARALLEL DO
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
           !$OMP END PARALLEL DO
        endif
        !$OMP PARALLEL DO
        do i=1,npoints
           press(i) = (gl_gamma - 1.0d0) * rho(i) * eps(i)
        enddo
        !$OMP END PARALLEL DO

     case (3)
        ! hybrid EOS
        !$OMP PARALLEL DO PRIVATE(hybrid_local_gamma,hybrid_local_k_cgs,&
        !$OMP                     hybrid_p_poly, hybrid_p_th)                          
        do i=1,npoints
           if(rho(i).gt.hybrid_rho_nuc) then
              hybrid_local_gamma = hybrid_gamma2
              hybrid_local_k_cgs = hybrid_k2_cgs
           else
              hybrid_local_gamma = hybrid_gamma1
              hybrid_local_k_cgs = hybrid_k1_cgs
           endif
           hybrid_p_poly = press_gf * hybrid_local_k_cgs * &
                     (rho(i) * inv_rho_gf)**hybrid_local_gamma
           hybrid_p_th = - press_gf * hybrid_local_k_cgs * (hybrid_gamma_th - 1.d0) /      &
                (hybrid_local_gamma - 1.0d0) * (rho(i) * inv_rho_gf)**hybrid_local_gamma + &
                (hybrid_gamma_th - 1.0d0) * rho(i) * eps(i) -                              &
                (hybrid_gamma_th - 1.d0) * (hybrid_local_gamma - hybrid_gamma1) /          &
                (hybrid_gamma1 - 1.d0) / (hybrid_gamma2 - 1.d0) *                          &
                press_gf * hybrid_k1_cgs * inv_rho_gf**hybrid_gamma1 *                     &
                hybrid_rho_nuc**(hybrid_gamma1 - 1.d0) * rho(i)
           hybrid_p_th = max(zero, hybrid_p_th)
           press(i) = hybrid_p_poly + hybrid_p_th
        enddo
        !$OMP END PARALLEL DO

     case (4)
        ! nuc eos
        my_anyerr = 0
        !$OMP PARALLEL REDUCTION(+: my_anyerr)
#ifdef _OPENMP
        num_threads = omp_get_num_threads()
        my_thread_num = omp_get_thread_num()
#else
        num_threads = 1
        my_thread_num = 0
#endif
        slice_len = (npoints + num_threads-1)/num_threads
        slice_min = 1 + my_thread_num * slice_len
        slice_max = min(slice_min + slice_len - 1, npoints)
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps(slice_max-slice_min+1,&
                rho(slice_min:slice_max),temp(slice_min:slice_max),&
                ye(slice_min:slice_max),eps(slice_min:slice_max),&
                press(slice_min:slice_max),keyerr(slice_min:slice_max),&
                my_anyerr)
        else
           call nuc_eos_m_kt0_press(slice_max-slice_min+1,&
                rho(slice_min:slice_max),temp(slice_min:slice_max),&
                ye(slice_min:slice_max),eps(slice_min:slice_max),&
                press(slice_min:slice_max),&
                rf_precision,keyerr(slice_min:slice_max),anyerr)
        endif
        !$OMP END PARALLEL
        ! return 0/1 for false/true rather than zero/nonzero just in case a
        ! caller relied on this
        if(my_anyerr.ne.0) then
           anyerr = 1
        end if
     case (5)
        ! cold tabular EOS with gamma law
        !$OMP PARALLEL DO PRIVATE(xrho,ir,gamma,eps_cold,eps_th, &
        !$OMP                     press_cold,press_th)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              press(i) = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = press(i) / (coldeos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              !$OMP CRITICAL
              anyerr = 1
              !$OMP END CRITICAL
           else
              xrho = log10(rho(i))
              ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th
        enddo
        !$OMP END PARALLEL DO

     case (6)
        ! barotropic tabular EOS with gamma law
        !$OMP PARALLEL DO PRIVATE(xrho,xprs,xenr,xye,ir, &
        !$OMP press_cold,eps_cold,eps_th,press_th)
        do i=1,npoints
           if(rho(i).lt.barotropiceos_rhomin) then
              press(i) = barotropiceos_low_kappa * rho(i)**barotropiceos_low_gamma
              eps(i) = press(i) / (barotropiceos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.barotropiceos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho = log10(rho(i))
              ir = 2 + &
                   INT( (xrho - barotropiceos_logrho(1) - 1.0d-10) &
                   * barotropiceos_dlrhoi )
           endif
           ir = max(2, min(ir,barotropiceos_nrho))

           xprs = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

           xenr = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)

           xtemp = (barotropiceos_temp(ir) - barotropiceos_temp(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_temp(ir-1)

           xye = (barotropiceos_ye(ir) - barotropiceos_ye(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_ye(ir-1)

           press_cold = 10.0**xprs
           eps_cold = 10.0**xenr - barotropiceos_energyshift

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif
           press_th = coldeos_thfac*(barotropiceos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th

        enddo
        !$OMP END PARALLEL DO


     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_PressOMP


subroutine EOS_Omni_EOS_DPressByDRho(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,dpdrhoe,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpdrhoe(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_dp_poly, hybrid_dp_th1, hybrid_dp_th2
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for cold tabulated EOS + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma, cs2, h

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdrhoe(i) = press_gf * poly_k_cgs *  &
                poly_gamma * inv_rho_gf *        & 
                (rho(i)*inv_rho_gf) ** (poly_gamma - 1.d0)
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdrhoe(i) = (gl_gamma-1.0d0) * &
                eps(i)
        enddo
     case (3)
        ! hybrid EOS
        do i=1,npoints
           if(rho(i).gt.hybrid_rho_nuc) then
              hybrid_local_gamma = hybrid_gamma2
              hybrid_local_k_cgs = hybrid_k2_cgs
           else
              hybrid_local_gamma = hybrid_gamma1
              hybrid_local_k_cgs = hybrid_k1_cgs
           endif
           hybrid_dp_poly = hybrid_local_gamma * press_gf *                 &
                hybrid_local_k_cgs * rho(i)**(hybrid_local_gamma - 1.0d0) * &
                inv_rho_gf**hybrid_local_gamma

           hybrid_dp_th1 = - hybrid_local_gamma * press_gf * hybrid_local_k_cgs * &
                (hybrid_gamma_th - 1.d0) / (hybrid_local_gamma - 1.d0) *      &
                rho(i)**(hybrid_local_gamma - 1.d0) * inv_rho_gf**hybrid_local_gamma

           hybrid_dp_th2  = (hybrid_gamma_th - 1.d0) * eps(i)                        &
                - (hybrid_gamma_th - 1.d0) * (hybrid_local_gamma - hybrid_gamma1) /  &
                (hybrid_gamma1 - 1.d0) / (hybrid_gamma2 - 1.d0) *                    &
                press_gf * hybrid_k1_cgs * inv_rho_gf**hybrid_gamma1 *               &
                hybrid_rho_nuc**(hybrid_gamma1 - 1.d0)

           dpdrhoe(i) = hybrid_dp_poly + max(0.0d0,hybrid_dp_th1 + hybrid_dp_th2)

        enddo

     case (4)
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,xdpderho,dpdrhoe,xmunu,&
                keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,xdpderho,dpdrhoe,xmunu,rf_precision,&
                keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not supported!")
           STOP
        endif
        
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! by definition dPdrho is at constant internal energy
        ! and entropy (the latter, because T=0)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              dpdrhoe(i) = coldeos_low_kappa * coldeos_low_gamma * &
                   rho(i)**(coldeos_low_gamma-1.0d0)
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho = log10(rho(i))
              ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)
           
           ! this is the cold speed of sound squared
           cs2 = (coldeos_cs2(ir) - coldeos_cs2(ir-1)) / &
                (coldeos_cs2(ir) - coldeos_cs2(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_cs2(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th
           
           ! this is the cold enthalpy, because it belongs to the
           ! cold speed of sound squared
           h = 1.0d0 + eps_cold + press_cold / rho(i)
           dpdrhoe(i) = cs2*h + press_th / rho(i)

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
   end select
     
end subroutine EOS_Omni_EOS_DPressByDRho

subroutine EOS_Omni_EOS_DPressByDEps(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,dpdepsrho,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpdepsrho(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law                                                              
  integer :: ir
  real*8 :: eps_cold, eps_th

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdepsrho(i) = 0.0d0
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdepsrho(i) = (gl_gamma - 1.0d0) * &
                rho(i)
        enddo
     case (3)
        ! hybrid EOS
        do i=1,npoints
           dpdepsrho(i) = (hybrid_gamma_th - 1.0d0) * rho(i)
        enddo
     case (4)
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,dpdepsrho,xdpdrhoe,xmunu,&
                keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,dpdepsrho,xdpdrhoe,xmunu,rf_precision,&
                keyerr,anyerr)
        else 
           call CCTK_ERROR("This keytemp is not supported!")
           STOP
        endif
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! only the gamma law has non-zero dPdeps
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              dpdepsrho(i) = 0.0d0
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho = log10(rho(i))
              ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           dpdepsrho(i) = coldeos_thfac * (coldeos_gammath - 1.0d0)*rho(i)

        enddo
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
   end select
     
end subroutine EOS_Omni_EOS_DPressByDEps

subroutine EOS_Omni_EOS_cs2(eoskey,keytemp,rf_precision,npoints,&
                            rho,eps,temp,ye,cs2,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: cs2(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: xpress,xdpdrhoe,xdpderho
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt
  ! temporary vars for cold tabulated EOS + gamma law                                                   
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma, cs2_cold, cs2_th, h, h_cold

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           xpress = press_gf*poly_k_cgs * &
                (rho(i)*inv_rho_gf)**(poly_gamma)
           cs2(i) = poly_gamma * xpress / rho(i) / &
                (1 + eps(i) + xpress/rho(i))
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           xpress = (gl_gamma-1.0d0)*rho(i)*eps(i)
           xdpdrhoe = (gl_gamma-1.0d0)*eps(i)
           xdpderho = (gl_gamma-1.0d0)*rho(i)
           cs2(i) = (xdpdrhoe + xpress * xdpderho / (rho(i)**2)) / &
                (1.0d0 + eps(i) + xpress/rho(i))
        enddo
     case(3)
        ! hybrid EOS
        do i=1,npoints
           if(rho(i).gt.hybrid_rho_nuc) then
              hybrid_local_gamma = hybrid_gamma2
              hybrid_local_k_cgs = hybrid_k2_cgs
                 else
              hybrid_local_gamma = hybrid_gamma1
              hybrid_local_k_cgs = hybrid_k1_cgs
           endif
           ! first calculate the pressure
           hybrid_p_poly = press_gf * hybrid_local_k_cgs * &
                     (rho(i) * inv_rho_gf)**hybrid_local_gamma
           hybrid_p_th = - press_gf * hybrid_local_k_cgs * (hybrid_gamma_th - 1.d0) /      &
                (hybrid_local_gamma - 1.0d0) * (rho(i) * inv_rho_gf)**hybrid_local_gamma + &
                (hybrid_gamma_th - 1.0d0) * rho(i) * eps(i) -                              &
                (hybrid_gamma_th - 1.d0) * (hybrid_local_gamma - hybrid_gamma1) /          &
                (hybrid_gamma1 - 1.d0) / (hybrid_gamma2 - 1.d0) *                          &
                press_gf * hybrid_k1_cgs * inv_rho_gf**hybrid_gamma1 *                     &
                hybrid_rho_nuc**(hybrid_gamma1 - 1.d0) * rho(i)
           hybrid_p_th = max(zero, hybrid_p_th)
           xpress = hybrid_p_poly + hybrid_p_th 
           cs2(i) = (hybrid_local_gamma * hybrid_p_poly + hybrid_gamma_th * hybrid_p_th) / &
                rho(i) / (1.0d0 + eps(i) + xpress/rho(i))
        enddo
     case(4)
        ! nuc_eos
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps_cs2(npoints,rho,temp,ye,&
                eps,xprs,cs2,keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_press_cs2(npoints,rho,temp,ye,&
                eps,xprs,cs2,rf_precision,keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not supported!")
           STOP
        endif
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! by definition dPdrho is at constant internal energy
        ! and entropy (the latter, because T=0)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              xprs = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = xprs / (coldeos_low_gamma - 1.0d0) / rho(i)
              cs2(i) = coldeos_low_kappa * coldeos_low_gamma * &
                   rho(i)**(coldeos_low_gamma-1.0d0) / &
                   (1.0d0 + eps(i) + xprs)
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
                 else
              xrho = log10(rho(i))
              ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

           ! this is the cold speed of sound squared                                                    
           cs2_cold = (coldeos_cs2(ir) - coldeos_cs2(ir-1)) / &
                (coldeos_cs2(ir) - coldeos_cs2(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_cs2(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th

           xdpdrhoe = coldeos_thfac*(coldeos_gammath - 1.0d0)*eps_th
           xdpderho = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)
           cs2_th = (xdpdrhoe + press_th * xdpderho / (rho(i)**2))

           h = 1.0d0 + eps(i) + (press_cold+press_th) / rho(i)
           h_cold = 1.0d0 + eps_cold + press_cold / rho(i)
           
           cs2(i) = (cs2_cold * h_cold + cs2_th) / h
        enddo
           
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
   end select
     
end subroutine EOS_Omni_EOS_cs2


subroutine EOS_Omni_EOS_eps_from_press(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,xeps,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints),press(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: xeps(npoints)

  ! local vars
  integer        :: i
  character(256) :: warnstring
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe

  if(keytemp.eq.1) then
     anyerr = 1
     keyerr(:) = -1
  else
     anyerr    = 0
     keyerr(:) = 0
  endif

  select case (eoskey)
     case (1)
        ! polytropic EOS
        do i=1,npoints
           xeps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
        enddo
     case (2)
        ! gamma-law EOS
        do i=1,npoints
           xeps(i) = press(i) / (gl_gamma - 1.0d0) / rho(i)
        enddo
     case (3)
        ! hybrid EOS
        write(warnstring,*) "EOS_Omni_EpsFromPress call not supported for hybrid EOS"
        call CCTK_ERROR(warnstring)
        STOP
     case (4)
        ! nuc EOS
        write(warnstring,*) "EOS_Omni_EpsFromPress call not supported for nuc_eos yet"
        call CCTK_ERROR(warnstring)
        STOP
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_eps_from_press


subroutine EOS_Omni_EOS_RhoFromPressEpsTempEnt(eoskey,&
                           ikeytemp,rf_precision,&
                           npoints,rho,eps,temp,entropy,ye,press,keyerr,anyerr)

  ! This routine returns the density and spec. internal energy based on inputs of
  ! pressure, temperature, and Y_e.
  ! The current implementation is robust, but very slow; it should be used only
  ! for initial data where speed does not matter.

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,ikeytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: ye(npoints),press(npoints)
  CCTK_REAL, intent(inout) :: rho(npoints),temp(npoints),eps(npoints)
  CCTK_REAL, intent(inout) :: entropy(npoints)


  ! local vars
  integer        :: i
  character(256) :: warnstring
  integer :: keytemp
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for cold tabulated EOS
  integer :: ir
  real*8  :: gamma
  
  ! helper vars
  logical :: cont
  integer :: counter
  real*8  :: rho_guess, rho_guess2
  real*8  :: press_guess,press_guess2, mydpdrho
  real*8  :: fac
  keytemp = ikeytemp

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        do i=1,npoints
           rho(i) = (press(i) / poly_k)**(1.0/poly_gamma)
           eps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
        enddo
     case (2)
        ! gamma-law EOS
        do i=1,npoints
           rho(i) = press(i) / ((gl_gamma - 1.0d0)*eps(i))
        enddo
     case (3)
        ! hybrid EOS
        write(warnstring,*) "EOS_Omni_RestMassDensityFromPressEpsTemp not supported for hybrid EOS"
        call CCTK_ERROR(warnstring)
        STOP
     case (4)
        ! nuc EOS
        if(ikeytemp.eq.2) then
           call CCTK_ERROR("This function does not work yet when coming in with entropy")
           STOP
        else if(ikeytemp.eq.1) then
           keytemp = 1
        else
           call CCTK_ERROR("This function does not work yet when coming in with this keytemp")
           STOP
        endif
        
        call CCTK_ERROR("This routine does not work with the new nuc_eos_cxx [yet]")
        STOP

        keytemp = 1
        do i=1,npoints
           fac = 1.0d0
           counter = 0
           cont = .true.
           xprs = press(i) 
           press_guess = xprs
           rho_guess = rho(i) 
           xprs = press(i) 
           xtemp = temp(i)
           xye = ye(i)
           do while(cont)
              counter = counter + 1
              rho_guess2 = rho_guess * 1.0001d0
              call nuc_eos_m_kt1_short(1,rho_guess2,xtemp,xye,xenr,&
                   press_guess2,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu,&
                   rf_precision,keyerr(i),anyerr)
              call nuc_eos_m_kt1_short(1,rho_guess2,xtemp,xye,xenr,&
                   press_guess2,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu,&
                   rf_precision,keyerr(i),anyerr)
              mydpdrho = (press_guess2-press_guess)/(rho_guess2-rho_guess)
              if (mydpdrho.lt.0.0d0) then
                 !$OMP CRITICAL
                 write(warnstring,"(A25,1P10E15.6)") "Issue with table, dpdrho.lt.0",&
                      rho_guess,xtemp,xye
                 call CCTK_ERROR(warnstring)
                 !$OMP END CRITICAL
                 STOP
              endif
              
              if (abs(1.0d0-press_guess/xprs).lt.rf_precision) then
                 cont = .false.
                 rho(i) = rho_guess
                 eps(i) = xenr
              else
                 if (fac*(xprs-press_guess)/mydpdrho/rho_guess.gt.0.1d0) then
                    rho_guess = 0.99d0*rho_guess
                 else
                    rho_guess = rho_guess + fac*(xprs-press_guess)/mydpdrho
                 endif
              endif

              if (counter.gt.100) then
                 fac = 0.01d0
              endif

              if (rho_guess.lt.1.0d3.or.counter.gt.100000) then
                 keyerr(i) = 473
                 anyerr = 1
                 return
              endif
           enddo
        enddo
        
     case (5)
        ! cold tabulated EOS
        ! deal only with case in which thermal pressure is zero
        if(keytemp.ne.1) then
           call CCTK_ERROR("finding rho(press) for tabulated cold EOS only possible if keytemp=1")
           STOP
        endif

        do i=1,npoints
           fac = 1.0d0
           counter = 0
           cont = .true.
           xprs = press(i) 
           press_guess = xprs
           rho_guess = rho(i) 
           xprs = press(i) 
           do while(cont)
              counter = counter + 1
              rho_guess2 = rho_guess * 1.0001d0

              ! press 2
              if(rho_guess2.lt.coldeos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess2.gt.coldeos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho = log10(rho_guess2)
                 ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
              endif
              ir = max(2, min(ir,coldeos_nrho))
           
              gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

              press_guess2 = coldeos_kappa * rho_guess2**gamma
              ! end press 2
              ! press 1
              if(rho_guess.lt.coldeos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess.gt.coldeos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho = log10(rho_guess)
                 ir = 2 + INT( (xrho - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
              endif
              ir = max(2, min(ir,coldeos_nrho))
           
              gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

              press_guess = coldeos_kappa * rho_guess**gamma
              ! end press 1

              ! derivative
              mydpdrho = (press_guess2-press_guess)/(rho_guess2-rho_guess)
              if (mydpdrho.lt.0.0d0) then
                 !$OMP CRITICAL
                 write(warnstring,"(A25,1P10E15.6)") "Issue with table, dpdrho.lt.0",&
                      rho_guess,rho_guess2
                 call CCTK_ERROR(warnstring)
                 !$OMP END CRITICAL
                 STOP
              endif
              
              if (abs(1.0d0-press_guess/xprs).lt.rf_precision) then
                 cont = .false.
                 rho(i) = rho_guess
                 eps(i) = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)
              else
                 if (fac*(xprs-press_guess)/mydpdrho/rho_guess.gt.0.1d0) then
                    rho_guess = 0.99d0*rho_guess
                 else
                    rho_guess = rho_guess + fac*(xprs-press_guess)/mydpdrho
                 endif
              endif

              if (counter.gt.100) then
                 fac = 0.01d0
              endif

              if (rho_guess.lt.1.0d-15.or.counter.gt.100000) then
                 keyerr(i) = 473
                 anyerr = 1
                 return
              endif
           enddo

        enddo

     case (6)
        ! barotropic tabulated EOS
        ! deal only with case in which thermal pressure is zero
        if(keytemp.ne.1) then
           call CCTK_ERROR("finding rho(press) for tabulated barotropic EOS only possible if keytemp=1")
           STOP
        endif

        do i=1,npoints
           fac = 1.0d0
           counter = 0
           cont = .true.
           xprs = press(i) 
           press_guess = xprs
           rho_guess = rho(i) 
           do while(cont)
              counter = counter + 1
              rho_guess2 = rho_guess * 1.0001d0

              ! press 2
              if(rho_guess2.lt.barotropiceos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess2.gt.barotropiceos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho = log10(rho_guess2)
                 ir = 2 + INT( (xrho - barotropiceos_logrho(1) - 1.0d-10) * barotropiceos_dlrhoi )
              endif
              ir = max(2, min(ir,barotropiceos_nrho))
           
              xprs = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                   (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                   (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)
              
              press_guess2 = 10.0d0**xprs
              ! end press 2
              ! press 1
              if(rho_guess.lt.barotropiceos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess.gt.barotropiceos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho = log10(rho_guess)
                 ir = 2 + INT( (xrho - barotropiceos_logrho(1) - 1.0d-10) * barotropiceos_dlrhoi )
              endif
              ir = max(2, min(ir,barotropiceos_nrho))
           
              xprs = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

              press_guess = 10.0d0**xprs
              ! end press 1

 !             write(6,"(i7,1P10E15.6)") counter, rho_guess, press_guess, press_guess2, press(i)

              ! derivative
              mydpdrho = (press_guess2-press_guess)/(rho_guess2-rho_guess)
              if (mydpdrho.lt.0.0d0) then
                 !$OMP CRITICAL
                 write(warnstring,"(A25,1P10E15.6)") "Issue with table, dpdrho.lt.0",&
                      rho_guess,rho_guess2
                 call CCTK_ERROR(warnstring)
                 !$OMP END CRITICAL
                 STOP
              endif
              
              if (abs(1.0d0-press_guess/press(i)).lt.rf_precision) then
                 cont = .false.
                 rho(i) = rho_guess
                 xenr = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                      (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                      (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)
                 eps(i) = 10.0d0**xenr - barotropiceos_energyshift
              else
                 if (fac*(press(i)-press_guess)/mydpdrho/rho_guess.gt.0.1d0) then
                    rho_guess = 0.99d0*rho_guess
                 else
                    rho_guess = rho_guess + fac*(press(i)-press_guess)/mydpdrho
                 endif
              endif

              if (counter.gt.100) then
                 fac = 0.01d0
              endif

              if (rho_guess.lt.1.0d-15.or.counter.gt.100000) then
!                 !$OMP CRITICAL
!                 write(warnstring,"(A25,1P10E15.6)") "rho(p) issue", rho_guess,press(i),press_guess
!                 call CCTK_ERROR(warnstring)
!                 STOP
!                 !$OMP END CRITICAL
                 keyerr(i) = 473
                 anyerr = 1
                 return
              endif
           enddo

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_RhoFromPressEpsTempEnt


subroutine EOS_Omni_EOS_PressEpsTempYe_from_Rho(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints)
  CCTK_REAL, intent(out)   :: eps(npoints), ye(npoints) 
  CCTK_REAL, intent(out)   :: temp(npoints),press(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr
  real*8           :: xprs
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (6)
        ! barotropic tabular EOS with gamma law
        do i=1,npoints
           if(rho(i).lt.barotropiceos_rhomin) then
              press(i) = barotropiceos_low_kappa * rho(i)**barotropiceos_low_gamma
              eps(i) = press(i) / (barotropiceos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.barotropiceos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho = log10(rho(i))
              ir = 2 + &
                   INT( (xrho - barotropiceos_logrho(1) - 1.0d-10) &
                   * barotropiceos_dlrhoi )
           endif
           ir = max(2, min(ir,barotropiceos_nrho))

           xprs = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

           xenr = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)

           xtemp = (barotropiceos_temp(ir) - barotropiceos_temp(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_temp(ir-1)

           xye = (barotropiceos_ye(ir) - barotropiceos_ye(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho - barotropiceos_logrho(ir-1)) + barotropiceos_ye(ir-1)

           press_cold = 10.0**xprs
           eps_cold = 10.0**xenr - barotropiceos_energyshift

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif
           press_th = coldeos_thfac*(barotropiceos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th

           temp(i) = xtemp
           ye(i) = xye

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented for EOS_Omni_EOS_PressEpsTempYe_from_Rho!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

   end subroutine EOS_Omni_EOS_PressEpsTempYe_from_Rho
