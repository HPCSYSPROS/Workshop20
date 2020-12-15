 /*@@
   @file      GRHydro_EigenproblemM.F90
   @date      August 30, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke, Pedro Montero, Joachim Frieben
   @desc
   Computes the spectral decomposition of a given state.
   Implements the analytical scheme devised by J. M. Ibanez
   et al., "Godunov Methods: Theory and Applications", New
   York, 2001, 485-503. The optimized method for computing
   the Roe flux in the special relativistic case is due to
   M. A. Aloy et al., Comput. Phys. Commun. 120 (1999)
   115-121, and has been extended to the general relativistic
   case as employed in this subroutine by J. Frieben, J. M.
   Ibanez, and J. Pons (in preparation).
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

module GRHydro_EigenproblemM
  implicit none


 /*@@
   @routine    eigenvaluesM
   @date       Aug 30, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc
   Computes the eigenvalues of the Jacobian matrix evaluated
   at the given state.
   @enddesc
   @calls
   @calledby
   @history
   Culled from the routines in GR3D, author Mark Miller.
   @endhistory

@@*/

CONTAINS

subroutine eigenvaluesM(handle,rho,velx,vely,velz,eps,press,w_lorentz,&
     Bvcx,Bvcy,Bvcz,lam,gxx,gxy,gxz,gyy,gyz,gzz,u,alp,beta)
  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, intent(in) :: rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL, intent(in) :: Bvcx,Bvcy,Bvcz
  CCTK_REAL, intent(out) :: lam(5)
  CCTK_REAL, intent(in) :: gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL, intent(in) :: alp,beta,u

  CCTK_REAL cs2,one,two,U2
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_INT, intent(in) :: handle
  CCTK_REAL dpdrho,dpdeps,press

  CCTK_REAL Bvcxlow,Bvcylow,Bvczlow,Bvc2,rhohstar,va2
  CCTK_REAL Bdotv,b2

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars            

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities


  call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdeps,keyerr,anyerr)

  call EOS_Omni_DPressByDRho(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdrho,keyerr,anyerr)

  cs2 = (dpdrho + press * dpdeps / (rho**2))/ &
       (1.0d0 + eps + press/rho)

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

!!$ Lower the B-field, and square of the magnitude
  Bvcxlow = gxx*Bvcx + gxy*Bvcy + gxz*Bvcz
  Bvcylow = gxy*Bvcx + gyy*Bvcy + gyz*Bvcz
  Bvczlow = gxz*Bvcx + gyz*Bvcy + gzz*Bvcz
  Bvc2 = Bvcxlow*Bvcx + Bvcylow*Bvcy + Bvczlow*Bvcz

  Bdotv = Bvcxlow*velx + Bvcylow*vely + Bvczlow*velz
  w = w_lorentz

  b2=Bvc2/w**2+Bdotv**2

 
!!$ rhohstar is the term that appears in Tmunu as well = rho*enth + b^2
  rhohstar = rho*(1.0+eps)+press+b2

!!$ Alfven velocity squared
  va2 = b2/(rhohstar)

!!$ The following combination always comes up in the wavespeed calculation:
!!$ U2 = v_a^2 + c_s^2(1-v_a^2)
!!$ In the unmagnetized case, it reduces to cs2
  U2 = va2+cs2*(1.d0-va2)

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-U2) + sqrt(U2*(one-v2)*&
       (u*(one-v2*U2) - velx**2*(one-U2))))/(one-v2*U2)
  lamm_nobeta = (velx*(one-U2) - sqrt(U2*(one-v2)*&
       (u*(one-v2*U2) - velx**2*(one-U2))))/(one-v2*U2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp

end subroutine eigenvaluesM

subroutine eigenvaluesM_hot(handle,ii,jj,kk,rho,velx,vely,velz,eps,press,w_lorentz,&
                             Bvcx,Bvcy,Bvcz,temp,ye,lam,gxx,gxy,gxz,gyy,gyz,gzz,u,alp,beta)
  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, intent(in) :: rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL, intent(in) :: Bvcx,Bvcy,Bvcz
  CCTK_REAL, intent(out) :: lam(5)
  CCTK_REAL, intent(in) :: gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL, intent(in) :: alp,beta,u
  CCTK_REAL, intent(in) :: temp,ye
  CCTK_REAL, intent(in) :: press

  CCTK_REAL cs2,one,two,U2
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_INT, intent(in) :: handle, ii,jj,kk

  CCTK_REAL Bvcxlow,Bvcylow,Bvczlow,Bvc2,rhohstar,va2
  CCTK_REAL Bdotv,b2

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars            

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities

  call EOS_Omni_cs2(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,temp,ye,cs2,keyerr,anyerr)

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

!!$ Lower the B-field, and square of the magnitude
  Bvcxlow = gxx*Bvcx + gxy*Bvcy + gxz*Bvcz
  Bvcylow = gxy*Bvcx + gyy*Bvcy + gyz*Bvcz
  Bvczlow = gxz*Bvcx + gyz*Bvcy + gzz*Bvcz
  Bvc2 = Bvcxlow*Bvcx + Bvcylow*Bvcy + Bvczlow*Bvcz

  Bdotv = Bvcxlow*velx + Bvcylow*vely + Bvczlow*velz
  w = w_lorentz

  b2=Bvc2/w**2+Bdotv**2

 
!!$ rhohstar is the term that appears in Tmunu as well = rho*enth + b^2
  rhohstar = rho*(1.0+eps)+press+b2

!!$ Alfven velocity squared
  va2 = b2/(rhohstar)

!!$ The following combination always comes up in the wavespeed calculation:
!!$ U2 = v_a^2 + c_s^2(1-v_a^2)
!!$ In the unmagnetized case, it reduces to cs2
  U2 = va2+cs2*(1.d0-va2)

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-U2) + sqrt(U2*(one-v2)*&
       (u*(one-v2*U2) - velx**2*(one-U2))))/(one-v2*U2)
  lamm_nobeta = (velx*(one-U2) - sqrt(U2*(one-v2)*&
       (u*(one-v2*U2) - velx**2*(one-U2))))/(one-v2*U2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp

end subroutine eigenvaluesM_hot
 
!!$ WE need to implement eigenproblem and eigenproblem_leftright!!!!


 /*@@
   @routine    eigenvalues_generalM
   @date       Aug 30, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc
   Computes the eigenvalues of the Jacobian matrix evaluated
   at the given state.
   @enddesc
   @calls
   @calledby
   @history
   Culled from the routines in GR3D, author Mark Miller.
   @endhistory

@@*/

subroutine eigenvalues_generalM(&
     velx,vely,velz,cs2,va2,&
     lam,&
     gxx,gxy,gxz,gyy,gyz,gzz,&
     u,alp,beta)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL velx,vely,velz
  CCTK_REAL lam(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL alp,beta,u,U2

  CCTK_REAL cs2,va2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = one / sqrt(one - v2)

  U2 = va2+cs2*(1-va2)

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-U2) + sqrt(U2*(one-v2)*&
       (u*(one-v2*U2) - velx**2*(one-U2))))/(one-v2*U2)
  lamm_nobeta = (velx*(one-U2) - sqrt(U2*(one-v2)*&
       (u*(one-v2*U2) - velx**2*(one-U2))))/(one-v2*U2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp

end subroutine eigenvalues_generalM


!!$ We'll need to implement eigenproblem_general

end module GRHydro_EigenproblemM

