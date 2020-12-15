#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! eoskey:
! 1 --- polytropic EOS
! 2 --- gamma-law EOS
! 3 --- hybrid EOS
! 4 --- finite-T microphysical NSE EOS

subroutine EOS_Omni_EOS_Press_f_hrho_v2_rhoW(eoskey,keytemp,rf_precision,npoints,&
                              hrho,v2,rhoW,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: hrho(npoints),v2(npoints)
  CCTK_REAL, intent(in)    :: rhoW(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: press(npoints)

  ! local vars
  integer          :: i
  real*8           :: gtmp1,gtmp2
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th
  real*8,parameter :: zero = 0.0d0

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        call CCTK_ERROR("Polytropic EOS not implemented for press_f_hro_v2_rhoW")
        STOP
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           call CCTK_ERROR("keytemp.eq.1 not implemented for press_f_hro_v2_rhoW")
           STOP
        endif
        do i=1,npoints
           gtmp1 = 1.0d0 - v2(i)
           gtmp2 = rhoW(i)*sqrt(gtmp1)
           press(i) = (gl_gamma - 1.0d0)/(gl_gamma) * &
                (hrho(i)*gtmp1 - gtmp2)
           eps(i) = press(i) / (gl_gamma - 1.0d0) / gtmp2
        enddo
     case (3)
        ! hybrid EOS
        call CCTK_ERROR("Hybrid EOS not implemented for press_f_hro_v2_rhoW")
        STOP

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_Press_f_hrho_v2_rhoW


subroutine EOS_Omni_EOS_dpdhrho_f_hrho_v2_rhoW(eoskey,keytemp,rf_precision,npoints,&
                              hrho,v2,rhoW,eps,temp,ye,dpdhrho,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: hrho(npoints),v2(npoints)
  CCTK_REAL, intent(in)    :: rhoW(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpdhrho(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th
  real*8,parameter :: zero = 0.0d0

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        call CCTK_ERROR("Polytropic EOS not implemented for press_f_hro_v2_rhoW")
        STOP
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           call CCTK_ERROR("keytemp.eq.1 not implemented for press_f_hro_v2_rhoW")
           STOP
        endif
        do i=1,npoints
           dpdhrho(i) = (gl_gamma - 1.0d0) * (1.0d0 - v2(i)) / (gl_gamma)
        enddo
     case (3)
        ! hybrid EOS
        call CCTK_ERROR("Hybrid EOS not implemented for press_f_hro_v2_rhoW")
        STOP

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_dpdhrho_f_hrho_v2_rhoW


subroutine EOS_Omni_EOS_dpdv2_f_hrho_v2_rhoW(eoskey,keytemp,rf_precision,npoints,&
                              hrho,v2,rhoW,eps,temp,ye,dpdv2,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: hrho(npoints),v2(npoints)
  CCTK_REAL, intent(in)    :: rhoW(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpdv2(npoints)

  ! local vars
  integer          :: i
  character(256)   :: warnstring
  real*8           :: hybrid_local_gamma, hybrid_local_k_cgs, &
                      hybrid_p_poly, hybrid_p_th
  real*8,parameter :: zero = 0.0d0

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        call CCTK_ERROR("Polytropic EOS not implemented for press_f_hro_v2_rhoW")
        STOP
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           call CCTK_ERROR("keytemp.eq.1 not implemented for press_f_hro_v2_rhoW")
           STOP
        endif
        do i=1,npoints
           dpdv2(i) = (gl_gamma - 1.0d0) * (0.5d0 * rhoW(i) &
                / sqrt(1.0d0 - v2(i)) - hrho(i)) / gl_gamma
        enddo
     case (3)
        ! hybrid EOS
        call CCTK_ERROR("Hybrid EOS not implemented for press_f_hro_v2_rhoW")
        STOP

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_dpdv2_f_hrho_v2_rhoW


