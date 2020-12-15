 /*@@
   @file      GRHydro_Eigenproblem.F90
   @date      Sat Jan 26 01:25:44 2002
   @author    Ian Hawke, Pedro Montero, Joachim Frieben
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

module GRHydro_Eigenproblem
  implicit none


 /*@@
   @routine    eigenvalues
   @date       Sat Jan 26 01:26:20 2002
   @author     Ian Hawke
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

subroutine eigenvalues(handle,rho,velx,vely,velz,eps, &
     w_lorentz,lam,gxx,gxy,gxz,gyy,gyz,gzz,u,alp,beta)
  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL lam(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL alp,beta,u

  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_INT handle
  CCTK_REAL dpdrho,dpdeps,press
  character*256 :: warnline

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0;xtemp=0.0d0;xye=0.0d0
! end EOS Omni vars

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities

!  call EOS_Omni_cs2(handle,keytemp,GRHydro_eos_rf_prec,n,&
!       rho,eps,xtemp,xye,cs2,keyerr,anyerr)
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,press,keyerr,anyerr)

  call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdeps,keyerr,anyerr)

  call EOS_Omni_DPressByDRho(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdrho,keyerr,anyerr)

  cs2 = (dpdrho + press * dpdeps / (rho**2))/ &
       (1.0d0 + eps + press/rho)

  if(cs2.lt.0.0d0) then
     if (abs(cs2) .gt. 1.0d-4) then
       !$OMP CRITICAL
        write(warnline,'(a60,6g16.7)') 'abs(cs2), rho, dpdrho, press*dpdeps/rho**2, eps, press/rho: ', abs(cs2), rho, dpdrho, press * dpdeps / (rho**2), eps, press/rho
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(1,"cs2 < 0! Check speed of sound calculation!")
        !$OMP END CRITICAL
        cs2 = 0.0d0
     else
        cs2 = 0.0d0
     endif
  endif


  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = w_lorentz

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp


end subroutine eigenvalues


subroutine eigenvalues_hot(handle,keytemp,ii,jj,kk,rho,velx,vely,velz,eps, &
     temp,ye,w_lorentz,lam,gxx,gxy,gxz,gyy,gyz,gzz,u,alp,beta)
  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL lam(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL alp,beta,u
  CCTK_REAL temp,ye

  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_INT handle,ii,jj,kk

  character(len=512) :: warnline

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xeps
  n=1;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xeps=0.0d0
! end EOS Omni vars

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities
  call EOS_Omni_cs2(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,temp,ye,cs2,keyerr,anyerr)
  if(anyerr.ne.0) then
     !$OMP CRITICAL
     call CCTK_WARN(1,"EOS ERROR in eigenvalues_hot")
     write(warnline,"(A10,i5,A10,i5)") "keyerr: ", keyerr, "keytemp: ",keytemp
     call CCTK_WARN(1,warnline)
     write(warnline,"(1P10E15.6)") rho,eps,temp,ye,cs2
     call CCTK_WARN(0,warnline)
     !$OMP END CRITICAL
  endif

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = w_lorentz

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp


end subroutine eigenvalues_hot

 /*@@
   @routine    eigenproblem
   @date       Sat Jan 26 01:27:59 2002
   @author     Ian Hawke, Pedro Montero, Joachim Frieben
   @desc
   Despite the name this routine currently actually returns the
   Roe flux given the input Roe average state.
   @enddesc
   @calls
   @calledby
   @history
   Culled and altered from the routines in GR3D, author Mark Miller.
   @endhistory

@@*/

subroutine eigenproblem(handle,rho,velx,vely,velz,eps,&
     w_lorentz,gxx,gxy,gxz,gyy,gyz,gzz,u,&
     alp,beta,qdiff1,qdiff2,qdiff3,qdiff4,qdiff5,&
     roeflux1,roeflux2,roeflux3,roeflux4,roeflux5)

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL lam(5),p(5,5),q(5,5),dw(5),rflux(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz,u
  CCTK_REAL alp,beta,roeflux1,roeflux2,roeflux3,roeflux4,roeflux5

  CCTK_REAL du(5),aa(5,5),qdiff1,qdiff2,qdiff3,qdiff4,qdiff5

  integer i,j
  CCTK_REAL paug(5,6),tmp1,tmp2,sump,summ,f_du(5)
  integer ii,jj,kk
  CCTK_REAL leivec1(5),leivec2(5),leivec3(5),leivecp(5),leivecm(5)
  CCTK_REAL reivec1(5),reivec2(5),reivec3(5),reivecp(5),reivecm(5)
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL press,dpdrho,dpdeps,enthalpy,kappa
  CCTK_REAL axp,axm,vxp,vxm,cxx,cxy,cxz,gam,xsi,dlt,vxa,vxb
  CCTK_INT handle

!!$  Warning, warning. Nasty hack follows

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1)
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars


  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,press,keyerr,anyerr)

  call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdeps,keyerr,anyerr)

  call EOS_Omni_DPressByDRho(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdrho,keyerr,anyerr)

  cs2 = (dpdrho + press * dpdeps / (rho**2))/ &
       (1.0d0 + eps + press/rho)
!  call EOS_Omni_cs2(handle,keytemp,GRHydro_eos_rf_prec,n,&
!       rho,eps,xtemp,xye,cs2,keyerr,anyerr)

!  if (cs2<0) cs2=0 ! this does not modify the roe crashing problem with shocktube
  enthalpy = one + eps + press / rho

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = w_lorentz

!!$Calculate eigenvalues and put them in conventional order

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp

  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

!!$  lam(1) = lamm
!!$  lam(2) = lam1
!!$  lam(3) = lam2
!!$  lam(4) = lam3
!!$  lam(5) = lamp

  lam(1) = lamm
  lam(5) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(2) = lamp

!!$Compute some auxiliary quantities

  axp = (u - velx*velx)/(u - velx*lamp_nobeta)
  axm = (u - velx*velx)/(u - velx*lamm_nobeta)
  vxp = (velx - lamp_nobeta)/(u - velx * lamp_nobeta)
  vxm = (velx - lamm_nobeta)/(u - velx * lamm_nobeta)

!!$Calculate associated right eigenvectors

  kappa = dpdeps / (dpdeps - rho * cs2)

  reivec1(1) = kappa / (enthalpy * w)
  reivec1(2) = vlowx
  reivec1(3) = vlowy
  reivec1(4) = vlowz
  reivec1(5) = one - reivec1(1)

  reivec2(1) = w * vlowy
  reivec2(2) = enthalpy * (gxy + two * w * w * vlowx * vlowy)
  reivec2(3) = enthalpy * (gyy + two * w * w * vlowy * vlowy)
  reivec2(4) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec2(5) = vlowy * w * (two * w * enthalpy - one)

  reivec3(1) = w * vlowz
  reivec3(2) = enthalpy * (gxz + two * w * w * vlowx * vlowz)
  reivec3(3) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec3(4) = enthalpy * (gzz + two * w * w * vlowz * vlowz)
  reivec3(5) = vlowz * w * (two * w * enthalpy - one)

  reivecp(1) = one
  reivecp(2) = enthalpy * w * (vlowx - vxp)
  reivecp(3) = enthalpy * w * vlowy
  reivecp(4) = enthalpy * w * vlowz
  reivecp(5) = enthalpy * w * axp - one

  reivecm(1) = one
  reivecm(2) = enthalpy * w * (vlowx - vxm)
  reivecm(3) = enthalpy * w * vlowy
  reivecm(4) = enthalpy * w * vlowz
  reivecm(5) = enthalpy * w * axm - one

!!$Calculate associated left eigenvectors if requested

  if (ANALYTICAL) then

    cxx = gyy * gzz - gyz * gyz
    cxy = gxz * gyz - gxy * gzz
    cxz = gxy * gyz - gxz * gyy
    gam = gxx * cxx + gxy * cxy + gxz * cxz
    xsi = cxx - gam * velx * velx
    dlt = enthalpy**3 * w * (kappa - one) * (vxm - vxp) * xsi

    tmp1 = w / (kappa - one)

    leivec1(1) = tmp1 * (enthalpy - w)
    leivec1(2) = tmp1 * w * velx
    leivec1(3) = tmp1 * w * vely
    leivec1(4) = tmp1 * w * velz
    leivec1(5) =-tmp1 * w

    tmp1 = one / (xsi * enthalpy)

    leivec2(1) = (gyz * vlowz - gzz * vlowy) * tmp1
    leivec2(2) = (gzz * vlowy - gyz * vlowz) * tmp1 * velx
    leivec2(3) = (gzz * (one - velx * vlowx) + gxz * vlowz * velx) * tmp1
    leivec2(4) = (gyz * (velx * vlowx - one) - gxz * velx * vlowy) * tmp1
    leivec2(5) = (gyz * vlowz - gzz * vlowy) * tmp1

    leivec3(1) = (gyz * vlowy - gyy * vlowz) * tmp1
    leivec3(2) = (gyy * vlowz - gyz * vlowy) * tmp1 * velx
    leivec3(3) = (gyz * (velx * vlowx - one) - gxy * velx * vlowz) * tmp1
    leivec3(4) = (gyy * (one - velx * vlowx) + gxy * velx * vlowy) * tmp1
    leivec3(5) = (gyz * vlowy - gyy * vlowz) * tmp1

    tmp1 = enthalpy * enthalpy / dlt
    tmp2 = w * w * xsi

    leivecp(1) = - (enthalpy * w * vxm * xsi + (one - kappa) * (vxm * &
      (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxm) * tmp1
    leivecp(2) = - (cxx * (one - kappa * axm) + (two * kappa - one) * vxm * &
      (tmp2 * velx - cxx * velx)) * tmp1
    leivecp(3) = - (cxy * (one - kappa * axm) + (two * kappa - one) * vxm * &
      (tmp2 * vely - cxy * velx)) * tmp1
    leivecp(4) = - (cxz * (one - kappa * axm) + (two * kappa - one) * vxm * &
      (tmp2 * velz - cxz * velx)) * tmp1
    leivecp(5) = - ((one - kappa) * (vxm * (tmp2 - cxx) - gam * velx) - &
      kappa * tmp2 * vxm) * tmp1

    leivecm(1) = (enthalpy * w * vxp * xsi + (one - kappa) * (vxp * &
      (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxp) * tmp1
    leivecm(2) = (cxx * (one - kappa * axp) + (two * kappa - one) * vxp * &
      (tmp2 * velx - cxx * velx)) * tmp1
    leivecm(3) = (cxy * (one - kappa * axp) + (two * kappa - one) * vxp * &
      (tmp2 * vely - cxy * velx)) * tmp1
    leivecm(4) = (cxz * (one - kappa * axp) + (two * kappa - one) * vxp * &
      (tmp2 * velz - cxz * velx)) * tmp1
    leivecm(5) = ((one - kappa) * (vxp * (tmp2 - cxx) - gam * velx) - &
      kappa * tmp2 * vxp) * tmp1

  endif

  du(1) = qdiff1
  du(2) = qdiff2
  du(3) = qdiff3
  du(4) = qdiff4
  du(5) = qdiff5

  if (ANALYTICAL) then

    if (FAST) then

      sump = 0.0d0
      summ = 0.0d0

      do i=1,5
         sump = sump + (abs(lamp) - abs(lam1)) * leivecp(i) * du(i)
         summ = summ + (abs(lamm) - abs(lam1)) * leivecm(i) * du(i)
      enddo

      vxa = sump + summ
      vxb =-(sump * vxp + summ * vxm)

      rflux(1) = abs(lam1) * du(1) + vxa
      rflux(2) = abs(lam1) * du(2) + enthalpy * w * (vlowx * vxa + vxb)
      rflux(3) = abs(lam1) * du(3) + enthalpy * w * (vlowy * vxa)
      rflux(4) = abs(lam1) * du(4) + enthalpy * w * (vlowz * vxa)
      rflux(5) = abs(lam1) * du(5) + enthalpy * w * (velx  * vxb + vxa) - vxa


    else

!!$Form Jacobian matrix in characteristic form from right eigenvectors.
!!$Invert to get the characteristic jumps given the conserved variable
!!$jumps

!!$  p(:,1) = reivecm(:)
!!$  p(:,2) = reivec1(:)
!!$  p(:,3) = reivec2(:)
!!$  p(:,4) = reivec3(:)
!!$  p(:,5) = reivecp(:)
!!$
!!$  write(*,*) p
!!$  stop

      p(:,1) = reivecm(:)
      p(:,2) = reivecp(:)
      p(:,3) = reivec2(:)
      p(:,4) = reivec3(:)
      p(:,5) = reivec1(:)

      q(1,:) = leivecm(:)
      q(2,:) = leivecp(:)
      q(3,:) = leivec2(:)
      q(4,:) = leivec3(:)
      q(5,:) = leivec1(:)

      do i=1,5
        dw(i) = 0.0d0
        do j=1,5
          dw(i) = dw(i) + q(i,j) * du(j)
        enddo
      enddo

!!$Calculate the Roe flux from the standard formula

      do i = 1, 5
        rflux(i) = 0.d0
        do j = 1, 5
          rflux(i) = rflux(i) + p(i,j) * abs(lam(j)) * dw(j)
        end do
      end do
    endif

  else

    p(:,1) = reivecm(:)
    p(:,2) = reivecp(:)
    p(:,3) = reivec2(:)
    p(:,4) = reivec3(:)
    p(:,5) = reivec1(:)

    do i=1,5
      dw(i) = du(i)
      do j=1,5
        aa(i,j) = p(i,j)
      end do
    enddo

    do ii=1,5
      paug(ii,1) = p(ii,1)
      paug(ii,2) = p(ii,2)
      paug(ii,3) = p(ii,3)
      paug(ii,4) = p(ii,4)
      paug(ii,5) = p(ii,5)
    enddo

    paug(1, 6) = du(1)
    paug(2, 6) = du(2)
    paug(3, 6) = du(3)
    paug(4, 6) = du(4)
    paug(5, 6) = du(5)

    do ii=1,5
      tmp1 = paug(ii,ii)
      do jj=ii,6
        paug(ii,jj) = paug(ii,jj)/tmp1
      enddo

      do jj=ii+1,5
        tmp1 = - (paug(jj,ii))
        do kk=ii,6
          paug(jj,kk) = paug(jj,kk) + tmp1*paug(ii,kk)
        enddo
      enddo

    enddo

    f_du(5) = paug(5,6)

    do ii=4,1,-1
      f_du(ii) = paug(ii,6)
      do jj=ii+1,5
        f_du(ii) = f_du(ii) - paug(ii,jj)*f_du(jj)
      enddo
    enddo

    dw(1) = f_du(1)
    dw(2) = f_du(2)
    dw(3) = f_du(3)
    dw(4) = f_du(4)
    dw(5) = f_du(5)

!!$  dw(1) = f_du(1)
!!$  dw(2) = f_du(5)
!!$  dw(3) = f_du(3)
!!$  dw(4) = f_du(4)
!!$  dw(5) = f_du(2)

!!$Calculate the Roe flux from the standard formula

    do i = 1, 5
      rflux(i) = 0.d0
      do j = 1, 5
        rflux(i) = rflux(i) + p(i,j) * abs(lam(j)) * dw(j)
      end do
    end do

  endif

  roeflux1 = rflux(1)
  roeflux2 = rflux(2)
  roeflux3 = rflux(3)
  roeflux4 = rflux(4)
  roeflux5 = rflux(5)

  if(roeflux1.ne.roeflux1) then
     call CCTK_WARN(0, "Found NaNs in roeflux1, aborting")
  endif

end subroutine eigenproblem

 /*@@
   @routine    eigenproblem_leftright
   @date       Sat Jan 26 01:27:59 2002
   @author     Ian Hawke, Pedro Montero, Joachim Frieben
   @desc
   Returns the left and right eigenvectors.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine eigenproblem_leftright(handle,rho,velx,vely,velz,eps,&
     w_lorentz,gxx,gxy,gxz,gyy,gyz,gzz,u,&
     alp,beta,lambda,levec,revec)

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL lambda(5),levec(5,5),revec(5,5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz,u
  CCTK_REAL alp,beta

  CCTK_REAL tmp1,tmp2
  CCTK_REAL leivec1(5),leivec2(5),leivec3(5),leivecp(5),leivecm(5)
  CCTK_REAL reivec1(5),reivec2(5),reivec3(5),reivecp(5),reivecm(5)
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL press,dpdeps,enthalpy,kappa
  CCTK_REAL axp,axm,vxp,vxm,cxx,cxy,cxz,gam,xsi,dlt
  CCTK_INT handle

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1)
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,press,keyerr,anyerr)

  call EOS_Omni_cs2(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,cs2,keyerr,anyerr)

  call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
       rho,eps,xtemp,xye,dpdeps,keyerr,anyerr)

!  if (cs2<0) cs2=0 ! this does not modify the roe crashing problem with shocktube
  enthalpy = one + eps + press / rho

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = w_lorentz

!!$Calculate eigenvalues and put them in conventional order

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp

  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

!!$  lam(1) = lamm
!!$  lam(2) = lam1
!!$  lam(3) = lam2
!!$  lam(4) = lam3
!!$  lam(5) = lamp

  lambda(1) = lamm
  lambda(2) = lam1
  lambda(3) = lam2
  lambda(4) = lam3
  lambda(5) = lamp

!!$Compute some auxiliary quantities

  axp = (u - velx*velx)/(u - velx*lamp_nobeta)
  axm = (u - velx*velx)/(u - velx*lamm_nobeta)
  vxp = (velx - lamp_nobeta)/(u - velx * lamp_nobeta)
  vxm = (velx - lamm_nobeta)/(u - velx * lamm_nobeta)

!!$Calculate associated right eigenvectors

  kappa = dpdeps / (dpdeps - rho * cs2)

  reivec1(1) = kappa / (enthalpy * w)
  reivec1(2) = vlowx
  reivec1(3) = vlowy
  reivec1(4) = vlowz
  reivec1(5) = one - reivec1(1)

  reivec2(1) = w * vlowy
  reivec2(2) = enthalpy * (gxy + two * w * w * vlowx * vlowy)
  reivec2(3) = enthalpy * (gyy + two * w * w * vlowy * vlowy)
  reivec2(4) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec2(5) = vlowy * w * (two * w * enthalpy - one)

  reivec3(1) = w * vlowz
  reivec3(2) = enthalpy * (gxz + two * w * w * vlowx * vlowz)
  reivec3(3) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec3(4) = enthalpy * (gzz + two * w * w * vlowz * vlowz)
  reivec3(5) = vlowz * w * (two * w * enthalpy - one)

  reivecp(1) = one
  reivecp(2) = enthalpy * w * (vlowx - vxp)
  reivecp(3) = enthalpy * w * vlowy
  reivecp(4) = enthalpy * w * vlowz
  reivecp(5) = enthalpy * w * axp - one

  reivecm(1) = one
  reivecm(2) = enthalpy * w * (vlowx - vxm)
  reivecm(3) = enthalpy * w * vlowy
  reivecm(4) = enthalpy * w * vlowz
  reivecm(5) = enthalpy * w * axm - one

  revec(1,:) = reivecm
  revec(2,:) = reivec1
  revec(3,:) = reivec2
  revec(4,:) = reivec3
  revec(5,:) = reivecp

!!$Calculate associated left eigenvectors if requested

  cxx = gyy * gzz - gyz * gyz
  cxy = gxz * gyz - gxy * gzz
  cxz = gxy * gyz - gxz * gyy
  gam = gxx * cxx + gxy * cxy + gxz * cxz
  xsi = cxx - gam * velx * velx
  dlt = enthalpy**3 * w * (kappa - one) * (vxm - vxp) * xsi
  
  tmp1 = w / (kappa - one)
  
  leivec1(1) = tmp1 * (enthalpy - w)
  leivec1(2) = tmp1 * w * velx
  leivec1(3) = tmp1 * w * vely
  leivec1(4) = tmp1 * w * velz
  leivec1(5) =-tmp1 * w
  
  tmp1 = one / (xsi * enthalpy)
  
  leivec2(1) = (gyz * vlowz - gzz * vlowy) * tmp1
  leivec2(2) = (gzz * vlowy - gyz * vlowz) * tmp1 * velx
  leivec2(3) = (gzz * (one - velx * vlowx) + gxz * vlowz * velx) * tmp1
  leivec2(4) = (gyz * (velx * vlowx - one) - gxz * velx * vlowy) * tmp1
  leivec2(5) = (gyz * vlowz - gzz * vlowy) * tmp1
  
  leivec3(1) = (gyz * vlowy - gyy * vlowz) * tmp1
  leivec3(2) = (gyy * vlowz - gyz * vlowy) * tmp1 * velx
  leivec3(3) = (gyz * (velx * vlowx - one) - gxy * velx * vlowz) * tmp1
  leivec3(4) = (gyy * (one - velx * vlowx) + gxy * velx * vlowy) * tmp1
  leivec3(5) = (gyz * vlowy - gyy * vlowz) * tmp1
  
  tmp1 = enthalpy * enthalpy / dlt
  tmp2 = w * w * xsi
  
  leivecp(1) = - (enthalpy * w * vxm * xsi + (one - kappa) * (vxm * &
       (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxm) * tmp1
  leivecp(2) = - (cxx * (one - kappa * axm) + (two * kappa - one) * vxm * &
       (tmp2 * velx - cxx * velx)) * tmp1
  leivecp(3) = - (cxy * (one - kappa * axm) + (two * kappa - one) * vxm * &
       (tmp2 * vely - cxy * velx)) * tmp1
  leivecp(4) = - (cxz * (one - kappa * axm) + (two * kappa - one) * vxm * &
       (tmp2 * velz - cxz * velx)) * tmp1
  leivecp(5) = - ((one - kappa) * (vxm * (tmp2 - cxx) - gam * velx) - &
       kappa * tmp2 * vxm) * tmp1
  
  leivecm(1) = (enthalpy * w * vxp * xsi + (one - kappa) * (vxp * &
       (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxp) * tmp1
  leivecm(2) = (cxx * (one - kappa * axp) + (two * kappa - one) * vxp * &
       (tmp2 * velx - cxx * velx)) * tmp1
  leivecm(3) = (cxy * (one - kappa * axp) + (two * kappa - one) * vxp * &
       (tmp2 * vely - cxy * velx)) * tmp1
  leivecm(4) = (cxz * (one - kappa * axp) + (two * kappa - one) * vxp * &
       (tmp2 * velz - cxz * velx)) * tmp1
  leivecm(5) = ((one - kappa) * (vxp * (tmp2 - cxx) - gam * velx) - &
       kappa * tmp2 * vxp) * tmp1

  levec(1,:) = leivecm
  levec(2,:) = leivec1
  levec(3,:) = leivec2
  levec(4,:) = leivec3
  levec(5,:) = leivecp

end subroutine eigenproblem_leftright

 /*@@
   @routine    eigenvalues
   @date       Sat Jan 26 01:26:20 2002
   @author     Ian Hawke
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

subroutine eigenvalues_general(&
     velx,vely,velz,cs2,&
     lam,&
     gxx,gxy,gxz,gyy,gyz,gzz,&
     u,alp,beta)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL velx,vely,velz
  CCTK_REAL lam(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL alp,beta,u

  CCTK_REAL cs2,one,two
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

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp

end subroutine eigenvalues_general

 /*@@
   @routine    eigenproblem_general
   @date       Sat Jan 26 01:27:59 2002
   @author     Ian Hawke, Pedro Montero, Joachim Frieben
   @desc
   Despite the name this routine currently actually returns the
   Roe flux given the input Roe average state.
   @enddesc
   @calls
   @calledby
   @history
   Culled and altered from the routines in GR3D, author Mark Miller.
   @endhistory

@@*/

subroutine eigenproblem_general(&
     rho,velx,vely,velz,eps,&
     press,cs2,dpdeps,&
     gxx,gxy,gxz,gyy,gyz,gzz,&
     u,alp,beta,&
     qdiff1,qdiff2,qdiff3,qdiff4,qdiff5,&
     roeflux1,roeflux2,roeflux3,roeflux4,roeflux5)

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL rho,velx,vely,velz,eps
  CCTK_REAL lam(5),p(5,5),q(5,5),dw(5),rflux(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz,u
  CCTK_REAL alp,beta,roeflux1,roeflux2,roeflux3,roeflux4,roeflux5

  CCTK_REAL du(5),aa(5,5),qdiff1,qdiff2,qdiff3,qdiff4,qdiff5

  integer i,j
  CCTK_REAL paug(5,6),tmp1,tmp2,sump,summ,f_du(5)
  integer ii,jj,kk
  CCTK_REAL leivec1(5),leivec2(5),leivec3(5),leivecp(5),leivecm(5)
  CCTK_REAL reivec1(5),reivec2(5),reivec3(5),reivecp(5),reivecm(5)
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL press,dpdeps,enthalpy,kappa
  CCTK_REAL axp,axm,vxp,vxm,cxx,cxy,cxz,gam,xsi,dlt,vxa,vxb

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities
  enthalpy = one + eps + press / rho

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = one / sqrt(one - v2)

!!$Calculate eigenvalues and put them in conventional order

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp

  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

!!$  lam(1) = lamm
!!$  lam(2) = lam1
!!$  lam(3) = lam2
!!$  lam(4) = lam3
!!$  lam(5) = lamp

  lam(1) = lamm
  lam(5) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(2) = lamp

!!$Compute some auxiliary quantities

  axp = (u - velx*velx)/(u - velx*lamp_nobeta)
  axm = (u - velx*velx)/(u - velx*lamm_nobeta)
  vxp = (velx - lamp_nobeta)/(u - velx * lamp_nobeta)
  vxm = (velx - lamm_nobeta)/(u - velx * lamm_nobeta)

!!$Calculate associated right eigenvectors

  kappa = dpdeps / (dpdeps - rho * cs2)

  reivec1(1) = kappa / (enthalpy * w)
  reivec1(2) = vlowx
  reivec1(3) = vlowy
  reivec1(4) = vlowz
  reivec1(5) = one - reivec1(1)

  reivec2(1) = w * vlowy
  reivec2(2) = enthalpy * (gxy + two * w * w * vlowx * vlowy)
  reivec2(3) = enthalpy * (gyy + two * w * w * vlowy * vlowy)
  reivec2(4) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec2(5) = vlowy * w * (two * w * enthalpy - one)

  reivec3(1) = w * vlowz
  reivec3(2) = enthalpy * (gxz + two * w * w * vlowx * vlowz)
  reivec3(3) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec3(4) = enthalpy * (gzz + two * w * w * vlowz * vlowz)
  reivec3(5) = vlowz * w * (two * w * enthalpy - one)

  reivecp(1) = one
  reivecp(2) = enthalpy * w * (vlowx - vxp)
  reivecp(3) = enthalpy * w * vlowy
  reivecp(4) = enthalpy * w * vlowz
  reivecp(5) = enthalpy * w * axp - one

  reivecm(1) = one
  reivecm(2) = enthalpy * w * (vlowx - vxm)
  reivecm(3) = enthalpy * w * vlowy
  reivecm(4) = enthalpy * w * vlowz
  reivecm(5) = enthalpy * w * axm - one

!!$Calculate associated left eigenvectors if requested

  if (ANALYTICAL) then

    cxx = gyy * gzz - gyz * gyz
    cxy = gxz * gyz - gxy * gzz
    cxz = gxy * gyz - gxz * gyy
    gam = gxx * cxx + gxy * cxy + gxz * cxz
    xsi = cxx - gam * velx * velx
    dlt = enthalpy**3 * w * (kappa - one) * (vxm - vxp) * xsi

    tmp1 = w / (kappa - one)

    leivec1(1) = tmp1 * (enthalpy - w)
    leivec1(2) = tmp1 * w * velx
    leivec1(3) = tmp1 * w * vely
    leivec1(4) = tmp1 * w * velz
    leivec1(5) =-tmp1 * w

    tmp1 = one / (xsi * enthalpy)

    leivec2(1) = (gyz * vlowz - gzz * vlowy) * tmp1
    leivec2(2) = (gzz * vlowy - gyz * vlowz) * tmp1 * velx
    leivec2(3) = (gzz * (one - velx * vlowx) + gxz * vlowz * velx) * tmp1
    leivec2(4) = (gyz * (velx * vlowx - one) - gxz * velx * vlowy) * tmp1
    leivec2(5) = (gyz * vlowz - gzz * vlowy) * tmp1

    leivec3(1) = (gyz * vlowy - gyy * vlowz) * tmp1
    leivec3(2) = (gyy * vlowz - gyz * vlowy) * tmp1 * velx
    leivec3(3) = (gyz * (velx * vlowx - one) - gxy * velx * vlowz) * tmp1
    leivec3(4) = (gyy * (one - velx * vlowx) + gxy * velx * vlowy) * tmp1
    leivec3(5) = (gyz * vlowy - gyy * vlowz) * tmp1

    tmp1 = enthalpy * enthalpy / dlt
    tmp2 = w * w * xsi

    leivecp(1) = - (enthalpy * w * vxm * xsi + (one - kappa) * (vxm * &
      (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxm) * tmp1
    leivecp(2) = - (cxx * (one - kappa * axm) + (two * kappa - one) * vxm * &
      (tmp2 * velx - cxx * velx)) * tmp1
    leivecp(3) = - (cxy * (one - kappa * axm) + (two * kappa - one) * vxm * &
      (tmp2 * vely - cxy * velx)) * tmp1
    leivecp(4) = - (cxz * (one - kappa * axm) + (two * kappa - one) * vxm * &
      (tmp2 * velz - cxz * velx)) * tmp1
    leivecp(5) = - ((one - kappa) * (vxm * (tmp2 - cxx) - gam * velx) - &
      kappa * tmp2 * vxm) * tmp1

    leivecm(1) = (enthalpy * w * vxp * xsi + (one - kappa) * (vxp * &
      (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxp) * tmp1
    leivecm(2) = (cxx * (one - kappa * axp) + (two * kappa - one) * vxp * &
      (tmp2 * velx - cxx * velx)) * tmp1
    leivecm(3) = (cxy * (one - kappa * axp) + (two * kappa - one) * vxp * &
      (tmp2 * vely - cxy * velx)) * tmp1
    leivecm(4) = (cxz * (one - kappa * axp) + (two * kappa - one) * vxp * &
      (tmp2 * velz - cxz * velx)) * tmp1
    leivecm(5) = ((one - kappa) * (vxp * (tmp2 - cxx) - gam * velx) - &
      kappa * tmp2 * vxp) * tmp1

  endif

  du(1) = qdiff1
  du(2) = qdiff2
  du(3) = qdiff3
  du(4) = qdiff4
  du(5) = qdiff5

  if (ANALYTICAL) then

    if (FAST) then

      sump = 0.0d0
      summ = 0.0d0

      do i=1,5
         sump = sump + (abs(lamp) - abs(lam1)) * leivecp(i) * du(i)
         summ = summ + (abs(lamm) - abs(lam1)) * leivecm(i) * du(i)
      enddo

      vxa = sump + summ
      vxb =-(sump * vxp + summ * vxm)

      rflux(1) = abs(lam1) * du(1) + vxa
      rflux(2) = abs(lam1) * du(2) + enthalpy * w * (vlowx * vxa + vxb)
      rflux(3) = abs(lam1) * du(3) + enthalpy * w * (vlowy * vxa)
      rflux(4) = abs(lam1) * du(4) + enthalpy * w * (vlowz * vxa)
      rflux(5) = abs(lam1) * du(5) + enthalpy * w * (velx  * vxb + vxa) - vxa

    else

!!$Form Jacobian matrix in characteristic form from right eigenvectors.
!!$Invert to get the characteristic jumps given the conserved variable
!!$jumps

!!$  p(:,1) = reivecm(:)
!!$  p(:,2) = reivec1(:)
!!$  p(:,3) = reivec2(:)
!!$  p(:,4) = reivec3(:)
!!$  p(:,5) = reivecp(:)
!!$
!!$  write(*,*) p
!!$  stop

      p(:,1) = reivecm(:)
      p(:,2) = reivecp(:)
      p(:,3) = reivec2(:)
      p(:,4) = reivec3(:)
      p(:,5) = reivec1(:)

      q(1,:) = leivecm(:)
      q(2,:) = leivecp(:)
      q(3,:) = leivec2(:)
      q(4,:) = leivec3(:)
      q(5,:) = leivec1(:)

      do i=1,5
        dw(i) = 0.0d0
        do j=1,5
          dw(i) = dw(i) + q(i,j) * du(j)
        enddo
      enddo

!!$Calculate the Roe flux from the standard formula

      do i = 1, 5
        rflux(i) = 0.d0
        do j = 1, 5
          rflux(i) = rflux(i) + p(i,j) * abs(lam(j)) * dw(j)
        end do
      end do
    endif

  else

    p(:,1) = reivecm(:)
    p(:,2) = reivecp(:)
    p(:,3) = reivec2(:)
    p(:,4) = reivec3(:)
    p(:,5) = reivec1(:)

    do i=1,5
      dw(i) = du(i)
      do j=1,5
        aa(i,j) = p(i,j)
      end do
    enddo

    do ii=1,5
      paug(ii,1) = p(ii,1)
      paug(ii,2) = p(ii,2)
      paug(ii,3) = p(ii,3)
      paug(ii,4) = p(ii,4)
      paug(ii,5) = p(ii,5)
    enddo

    paug(1, 6) = du(1)
    paug(2, 6) = du(2)
    paug(3, 6) = du(3)
    paug(4, 6) = du(4)
    paug(5, 6) = du(5)

    do ii=1,5
      tmp1 = paug(ii,ii)
      do jj=ii,6
        paug(ii,jj) = paug(ii,jj)/tmp1
      enddo

      do jj=ii+1,5
        tmp1 = - (paug(jj,ii))
        do kk=ii,6
          paug(jj,kk) = paug(jj,kk) + tmp1*paug(ii,kk)
        enddo
      enddo

    enddo

    f_du(5) = paug(5,6)

    do ii=4,1,-1
      f_du(ii) = paug(ii,6)
      do jj=ii+1,5
        f_du(ii) = f_du(ii) - paug(ii,jj)*f_du(jj)
      enddo
    enddo

    dw(1) = f_du(1)
    dw(2) = f_du(2)
    dw(3) = f_du(3)
    dw(4) = f_du(4)
    dw(5) = f_du(5)
!!$
!!$  dw(1) = f_du(1)
!!$  dw(2) = f_du(5)
!!$  dw(3) = f_du(3)
!!$  dw(4) = f_du(4)
!!$  dw(5) = f_du(2)

!!$Calculate the Roe flux from the standard formula

    do i = 1, 5
      rflux(i) = 0.d0
      do j = 1, 5
        rflux(i) = rflux(i) + p(i,j) * abs(lam(j)) * dw(j)
      end do
    end do

  endif

  roeflux1 = rflux(1)
  roeflux2 = rflux(2)
  roeflux3 = rflux(3)
  roeflux4 = rflux(4)
  roeflux5 = rflux(5)

end subroutine eigenproblem_general

end module GRHydro_Eigenproblem

