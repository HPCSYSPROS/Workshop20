#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! eoskey:
! 1 --- polytropic EOS
! 2 --- gamma-law EOS
! 3 --- hybrid EOS
! 4 --- finite-T microphysical NSE EOS


subroutine EOS_Omni_EOS_short(eoskey,keytemp,rf_precision,npoints,&
     rho,eps,temp,ye,press,entropy,cs2,dedt,dpderho,dpdrhoe,munu,&
     keyerr,anyerr)

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
  CCTK_REAL, intent(inout) :: entropy(npoints)
  CCTK_REAL, intent(out)   :: cs2(npoints)
  CCTK_REAL, intent(out)   :: dedt(npoints)
  CCTK_REAL, intent(out)   :: dpderho(npoints)
  CCTK_REAL, intent(out)   :: dpdrhoe(npoints)
  CCTK_REAL, intent(out)   :: munu(npoints)

  ! local vars
#if 0
  integer          :: i
#endif
  character(256)   :: warnstring
#if 0
  ! temporary vars for nuc_eos
  real*8           :: xrho,xye,xtemp,xenr,xent
  real*8           :: xprs,xmunu,xcs2
  real*8           :: xdedt,xdpderho,xdpdrhoe
#endif

  if(eoskey.ne.4) then
     ! CRITICAL is not really required but might prevent multiple threads from
     ! outputting error messages at the same time and cluttering the log file
     !$OMP CRITICAL
     write(warnstring,"(A8,i5)") "eoskey: ", eoskey
     call CCTK_WARN(1,warnstring)
     call CCTK_ERROR("EOS_Omni_EOS_short currently does not work for this eoskey")
     !$OMP END CRITICAL
     STOP
  endif

  anyerr    = 0
  keyerr(:) = 0

  if(keytemp.eq.1) then
     call nuc_eos_m_kt1_short(npoints,rho,temp,ye,eps,press,&
          entropy,cs2,dedt,dpderho,dpdrhoe,munu,keyerr,anyerr)
  else if(keytemp.eq.0) then
     call nuc_eos_m_kt0_short(npoints,rho,temp,ye,eps,press,&
          entropy,cs2,dedt,dpderho,dpdrhoe,munu,rf_precision,&
          keyerr,anyerr)
  else if (keytemp.eq.2) then
     call nuc_eos_m_kt2_short(npoints,rho,temp,ye,eps,press,&
          entropy,cs2,dedt,dpderho,dpdrhoe,munu,rf_precision,&
          keyerr,anyerr)
  else
     !$OMP CRITICAL
     call CCTK_ERROR("This keytemp is not supported")
     !$OMP END CRITICAL
     STOP
  endif

end subroutine EOS_Omni_EOS_short

subroutine EOS_Omni_EOS_full(eoskey,keytemp,rf_precision,npoints,&
     rho,eps,temp,ye,press,entropy,cs2,dedt,dpderho,dpdrhoe,&
     xa,xh,xn,xp,abar,zbar,mue,mun,mup,muhat,keyerr,anyerr)

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
  CCTK_REAL, intent(inout) :: entropy(npoints)
  CCTK_REAL, intent(out)   :: cs2(npoints)
  CCTK_REAL, intent(out)   :: dedt(npoints)
  CCTK_REAL, intent(out)   :: dpderho(npoints)
  CCTK_REAL, intent(out)   :: dpdrhoe(npoints)
  CCTK_REAL, intent(out)   :: xa(npoints)
  CCTK_REAL, intent(out)   :: xh(npoints)
  CCTK_REAL, intent(out)   :: xn(npoints)
  CCTK_REAL, intent(out)   :: xp(npoints)
  CCTK_REAL, intent(out)   :: abar(npoints)
  CCTK_REAL, intent(out)   :: zbar(npoints)
  CCTK_REAL, intent(out)   :: mue(npoints)
  CCTK_REAL, intent(out)   :: mun(npoints)
  CCTK_REAL, intent(out)   :: mup(npoints)
  CCTK_REAL, intent(out)   :: muhat(npoints)

  ! local vars
  character(256)   :: warnstring

  if(eoskey.ne.4) then
     !$OMP CRITICAL
     write(warnstring,"(A8,i5)") "eoskey: ", eoskey
     call CCTK_WARN(1,warnstring)
     call CCTK_ERROR("EOS_Omni_EOS_full currently does not work for this eoskey")
     !$OMP END CRITICAL
     STOP
  endif

  anyerr    = 0
  keyerr(:) = 0

  if(keytemp.eq.1) then
     call nuc_eos_m_kt1_full(npoints,rho,temp,ye,eps,press,&
          entropy,cs2,dedt,dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,&
          mue,mun,mup,muhat,keyerr,anyerr)
  else if(keytemp.eq.0) then
     call nuc_eos_m_kt0_full(npoints,rho,temp,ye,eps,press,&
          entropy,cs2,dedt,dpderho,dpdrhoe,xa,xh,xn,xp,abar,zbar,&
          mue,mun,mup,muhat,rf_precision,&
          keyerr,anyerr)
  else
     !$OMP CRITICAL
     call CCTK_ERROR("This keytemp is not supported")
     !$OMP END CRITICAL
     STOP
  endif

end subroutine EOS_Omni_EOS_full


subroutine EOS_Omni_EOS_dpderho_dpdrhoe(eoskey,keytemp,rf_precision,npoints,&
     rho,eps,temp,ye,dpderho,dpdrhoe,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpderho(npoints)
  CCTK_REAL, intent(out)   :: dpdrhoe(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma
  real*8           :: hybrid_local_k_cgs
  real*8           :: hybrid_dp_poly,hybrid_dp_th1,hybrid_dp_th2

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
           dpderho(i) = 0.0d0
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
           dpderho(i) = (gl_gamma - 1.0d0) * &
                rho(i)
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

           dpdrhoe(i) = hybrid_dp_poly + hybrid_dp_th1 + hybrid_dp_th2
           dpderho(i) = (hybrid_gamma_th - 1.0d0) * rho(i)
        enddo
    case (4)
       if(keytemp.ne.0) then
          call CCTK_ERROR("Keytemp other than 0 not supported for dpdrhoe, dpderho")
          STOP
       else
          call nuc_eos_m_kt0_dpdrhoe_dpderho(npoints,&
               rho,temp,ye,eps,dpdrhoe,dpderho,rf_precision,&
               keyerr,anyerr)
       endif
   case DEFAULT
      write(warnstring,*) "eoskey ",eoskey," not implemented!"
      call CCTK_ERROR(warnstring)
      STOP
   end select

 end subroutine EOS_Omni_EOS_dpderho_dpdrhoe


 subroutine EOS_Omni_EOS_DEpsByDRho_DEpsByDPress(eoskey,keytemp,rf_precision,&
      npoints,rho,eps,temp,ye,depsdrho,depsdpress,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: depsdrho(npoints)
  CCTK_REAL, intent(out)   :: depsdpress(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma
  real*8           :: hybrid_local_k_cgs
  real*8           :: hybrid_dp_poly,hybrid_dp_th1,hybrid_dp_th2
  ! temporary vars for nuc_eos
  real*8           :: xdpderho,xdpdrhoe


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
           depsdpress(i) = 1.0d0/(poly_gamma - 1.0d0)/rho(i)
           depsdrho(i) = depsdpress(i) * press_gf * poly_k_cgs *  &
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
           depsdpress(i) = 1.0/( (gl_gamma - 1.0d0) * &
                rho(i))
           depsdrho(i) = -eps(i)/rho(i)
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

           xdpdrhoe = hybrid_dp_poly + hybrid_dp_th1 + hybrid_dp_th2
           xdpderho = (hybrid_gamma_th - 1.0d0) * rho(i)
           depsdpress(i) = 1.0 / xdpderho
           depsdrho(i) = - xdpdrhoe * depsdpress(i)
        enddo
    case (4)
      write(warnstring,*) "depsdrho and depsdpress not implemented yet for hot nuclear EOS"
      call CCTK_ERROR(warnstring)
      STOP
   case DEFAULT
      write(warnstring,*) "eoskey ",eoskey," not implemented!"
      call CCTK_ERROR(warnstring)
      STOP
   end select

 end subroutine EOS_Omni_EOS_DEpsByDRho_DEpsByDPress


subroutine EOS_Omni_EOS_Press_cs2(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,cs2,keyerr,anyerr)

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
  CCTK_REAL, intent(out)   :: cs2(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th, pwp_a
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8           :: xrho, xpress
  real*8           :: xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law
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
           press(i) = press_gf * poly_k_cgs * &
                (rho(i)*inv_rho_gf)**poly_gamma
           cs2(i) = poly_gamma * press(i) / rho(i) / &
                (1 + eps(i) + press(i)/rho(i))
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
           xdpdrhoe = (gl_gamma-1.0d0)*eps(i)
           xdpderho = (gl_gamma-1.0d0)*rho(i)
           cs2(i) = (xdpdrhoe + press(i) * xdpderho / (rho(i)**2)) / &
                (1.0d0 + eps(i) + press(i)/rho(i))
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

           cs2(i) = (hybrid_local_gamma * hybrid_p_poly + hybrid_gamma_th * hybrid_p_th) / &
                rho(i) / (1.0d0 + eps(i) + press(i)/rho(i))
        enddo

     case (4)
        ! nuc eos
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps_cs2(npoints,rho,temp,ye,&
                eps,press,cs2,keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_press_cs2(npoints,rho,temp,ye,eps,press,&
                cs2,rf_precision,keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not suppported!")
           STOP
        endif
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! by definition dPdrho is at constant internal energy
        ! and entropy (the latter, because T=0)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              press(i) = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = press(i) / (coldeos_low_gamma - 1.0) / rho(i) 
              cs2(i) = coldeos_low_kappa * coldeos_low_gamma * &
                   rho(i)**(coldeos_low_gamma-1.0d0) / &
                   (1.0d0 + eps(i) + press(i))
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
           press(i) = press_cold + press_th

           xdpdrhoe = coldeos_thfac*(coldeos_gammath - 1.0d0)*eps_th
           xdpderho = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)
           cs2_th = (xdpdrhoe + press_th * xdpderho / (rho(i)**2))
           
           h = 1.0d0 + eps(i) + (press_cold+press_th) / rho(i)
           h_cold = 1.0d0 + eps_cold + press_cold / rho(i)
           
           cs2(i) = (cs2_cold * h_cold + cs2_th) / h

        enddo

     case (6)
        ! barotropic tabular EOS with gamma law
        write(warnstring,*) "eoskey ",eoskey," for barotropic EOS not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     
     case (7)
        ! piecwise polytrope hybrid EOS
        do i=1,npoints
           if(rho(i).lt.pwp_rho_1) then
              hybrid_local_gamma = pwp_gamma1
              hybrid_local_k_cgs = pwp_k1
              pwp_a = pwp_a1
           else if (rho(i).lt.pwp_rho_2) then
              hybrid_local_gamma = pwp_gamma2
              hybrid_local_k_cgs = pwp_k2
              pwp_a = pwp_a2
           else if (rho(i).lt.pwp_rho_3) then
              hybrid_local_gamma = pwp_gamma3
              hybrid_local_k_cgs = pwp_k3
              pwp_a = pwp_a3
           else
              hybrid_local_gamma = pwp_gamma4
              hybrid_local_k_cgs = pwp_k4
              pwp_a = pwp_a4
           endif
           hybrid_p_poly = hybrid_local_k_cgs * &
                     (rho(i))**hybrid_local_gamma
           eps_cold = pwp_a + hybrid_local_k_cgs * &
                     (rho(i))**(hybrid_local_gamma-1.0d0) / (hybrid_local_gamma-1.0d0)
           eps_th = eps(i) - eps_cold

           hybrid_p_th = (pwp_gamma_th - 1.0d0) * rho(i) * eps_th
           hybrid_p_th = max(zero, hybrid_p_th)
           xpress = hybrid_p_poly + hybrid_p_th
                     
           xdpdrhoe = hybrid_local_gamma * hybrid_local_k_cgs * (rho(i))**(hybrid_local_gamma-1.0d0) &
                      + (pwp_gamma_th-1.0d0) * max(0.0d0, eps_th)

           xdpderho = (pwp_gamma_th - 1.0d0) * rho(i)
           
           cs2(i) = (xdpdrhoe + xpress*xdpderho/(rho(i)*rho(i)))/(1.0d0 + eps(i) + xpress/rho(i))
           press(i) = xpress
       enddo
     
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_Press_cs2

subroutine EOS_Omni_EOS_eps_temp_cs2_from_press(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,cs2,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  ! note that keytemp is ignored in this routine, since we are using
  ! pressure, density, and Y_e as independent variables
  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(in)    :: press(npoints)
  CCTK_REAL, intent(out)   :: cs2(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th, pwp_a
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8           :: xrho, xpress
  real*8           :: xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma, cs2_cold, cs2_th, h, h_cold

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        do i=1,npoints
           eps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
           cs2(i) = poly_gamma * press(i) / rho(i) / &
                (1 + eps(i) + press(i)/rho(i))
        enddo
           
     case (2)
        ! gamma-law EOS
        do i=1,npoints
           eps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
           xdpdrhoe = (gl_gamma-1.0d0)*eps(i)
           xdpderho = (gl_gamma-1.0d0)*rho(i)
           cs2(i) = (xdpdrhoe + press(i) * xdpderho / (rho(i)**2)) / &
                (1.0d0 + eps(i) + press(i)/rho(i))
        enddo

     case (3)
        ! hybrid EOS
        call CCTK_ERROR("EOS(rho,press,Y_e) not supported for hybrid EOS!")
        STOP

     case (4)
        ! nuc eos
        call CCTK_ERROR("EOS(rho,press,Y_e) not yet supported for nuc_eos!")
        STOP

     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! by definition dPdrho is at constant internal energy
        ! and entropy (the latter, because T=0)
        call CCTK_ERROR("EOS(rho,press,Y_e) not yet supported for cold EOS + thermal gamma law!")
        STOP

     case (6)
        ! barotropic tabular EOS with gamma law
        call CCTK_ERROR("EOS(rho,press,Y_e) not yet supported for barotropic EOS!")
        STOP
     
     case (7)
        ! piecewise polytrope hybrid EOS
        call CCTK_ERROR("EOS(rho,press,Y_e) not yet supported for piecewise polytrope + hybrid EOS!")
        STOP

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_eps_temp_cs2_from_press

