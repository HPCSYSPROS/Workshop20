/*@@
@file      GRHydro_EOSResetHydro.F90
@date      Sat Jan 26 01:36:57 2002
@author     Ian Hawke
@desc 
   This routine will reset the specific internal energy using the polytropic
   EOS that will be used at evolution. This is wanted if the EoS changes
   between setting up the initial data and evolving
  @enddesc 
  @@*/
  
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"
 
#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

 /*@@
   @routine    GRHydro_EOSResetHydro
   @date       Sat Jan 26 01:38:12 2002
   @author     Ian Hawke
   @desc 
   see above
   @enddesc 
   @calls     
   @calledby   
   @history 
   
   @endhistory 

@@*/


subroutine GRHydro_EoSChangeGamma(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  
  CCTK_REAL :: local_gamma

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

!!$  Set up the fluid constants
! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       1.0d0,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
  
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       1.0d0,1.0d0,xtemp,xye,xpress,xeps,keyerr,anyerr)

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if
  
  local_Gamma = 1.0d0 + xpress/xeps
  press = poly_k * &
       rho**local_Gamma 
  eps = press / (rho * (local_Gamma - 1.d0))

!!$  Change the pressure and specific internal energy

!!$  Get the conserved variables. Hardwired polytrope EoS!!!
!!$  Note that this call also sets pressure and eps
    
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        call prim2conpolytype(GRHydro_polytrope_handle,g11(i,j,k),g12(i,j,k),&
             g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
             sdetg(i,j,k), dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

      end do
    end do
  end do


end subroutine GRHydro_EoSChangeGamma

 /*@@
   @routine    GRHydro_EoSChangeK
   @date       Mon Oct 20 12:56:14 2003
   @author     Ian Hawke
   @desc 
   Reset the hydro variables when K is changed.
   Unlike the routine above, this actually gives a solution to
   the constraints.

   Only two cases are given as the general case is transcendental.
   We find this by holding rho * enthalpy fixed and assuming a
   polytropic EOS.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_EoSChangeK(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  
  CCTK_REAL :: local_gamma, local_k
  
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: Q

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

!!$  Set up the fluid constants
! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       1.0d0,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
  
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       1.0d0,1.0d0,xtemp,xye,xpress,xeps,keyerr,anyerr)

  local_Gamma = 1.0d0 + xpress/xeps
  local_K = xpress

  call CCTK_INFO("Adjusting EOS via K change!")

  if (abs(local_Gamma - 2.d0) < 1.d-10) then

    rho = -0.5d0/local_k+sqrt(0.25d0/local_k**2+(rho+initial_k*rho**2)/local_k)

  else if (abs(local_Gamma - 3.d0) < 1.d-10) then

     call CCTK_ERROR("This part of the code is not tested!")
     STOP

!!$ This code is probably just wrong. We have never used it anyway.
    Q = -9.d0 * local_k**2 * rho * (2.d0 + 3.d0 * initial_k * rho**2) + &
         sqrt(local_k**3 * (32.d0 + 81.d0 * local_k * rho**2 * &
         (2.d0 + 3.d0 * initial_k * rho**2)**2))
    
    Q = Q**(1.d0/3.d0)
    
    rho = (2**(7.d0/3.d0) * local_k - 2**(2.d0/3.d0) * Q**2) / &
         (6.d0 * local_k * Q)

  else
    call CCTK_ERROR("EoSChangeK only knows how to do Gamma=2 or 3!")
    STOP
  end if
 
  press = local_k * rho**local_gamma
  eps = local_k / (local_gamma - 1.d0) * rho**(local_gamma-1.0d0)

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        call prim2conpolytype(GRHydro_polytrope_handle,g11(i,j,k),g12(i,j,k),&
             g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
             sdetg(i,j,k), dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

      end do
    end do
  end do



end subroutine GRHydro_EoSChangeK



 /*@@
   @routine    GRHydro_EoSChangeGammaK_Shibata
   @date       Jan. 2005
   @author     Christian D. Ott
   @desc 
   Reset the hydro variables when K and Gamma are changed.

   This is according to Shibata in astro-ph/0412243 (PRD71 024014) in
   which he switches K and Gamma after initial data setup,
   but keeps the internal energy constant.

   Note: this works only with one refinement level

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_EoSChangeGammaK_Shibata(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  
  CCTK_REAL :: local_Gamma, local_k, eos_k_initial

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars

!!$  Set up the fluid constants
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       1.0d0,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
  
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       1.0d0,1.0d0,xtemp,xye,xpress,xeps,keyerr,anyerr)

  ! save memory when MP is not used
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  local_Gamma = 1.0d0 + xpress/xeps
  local_K    = xpress

  eos_k_initial = initial_k

  press = (local_Gamma - 1.d0) / (initial_Gamma - 1.0d0 ) * eos_k_initial * &
               rho ** initial_Gamma

  eps = eos_k_initial * & 
     rho ** initial_Gamma / &
     (rho * (initial_Gamma - 1.0d0))

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        call prim2con(GRHydro_eos_handle,g11(i,j,k),g12(i,j,k),&
             g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
             sdetg(i,j,k), dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

      end do
    end do
  end do


end subroutine GRHydro_EoSChangeGammaK_Shibata
