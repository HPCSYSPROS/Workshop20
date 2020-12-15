/*@@
@file      GRHydro_FluxM.F90
@date      August 30, 2010
@author    Joshua Faber, Scott Noble, Bruno Mundim, Pedro Montero, Ian Hawke
@desc 
The routine to calculate the numerical flux function given a 
  specific state
  @enddesc 
  @@*/
  
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
 
subroutine num_x_fluxM(dens,sx,sy,sz,tau,Bx,By,Bz,&
     densf,sxf,syf,szf,tauf,Bxf,Byf,Bzf,vxt,vyt,vzt,pressstar,&
     bsubx,bsuby,bsubz,ab0,w,det,alp,beta)

  implicit none

  CCTK_REAL :: dens,sx,sy,sz,tau,Bx,By,Bz
  CCTK_REAL :: densf,sxf,syf,szf,tauf,Bxf,Byf,Bzf
  CCTK_REAL :: vxt,vyt,vzt,bsubx,bsuby,bsubz,ab0,w
  CCTK_REAL :: det,alp,beta,pressstar
  CCTK_REAL :: velm
  CCTK_REAL :: sdet,psipstar

  sdet=sqrt(det)
  psipstar=pressstar*sdet

!!$ We actually need all three values of vtilde = v^i - beta^i/alp, as well as
  velm=vxt+beta/alp

!!$ GRHydro splits off alpha for later in the calculation, so we have psi^6 * [Anton eq.42]
!!$ In the notation of Anton et al.:  [psi^6 D] vtilde^i
  densf = dens * vxt

  sxf = sx*vxt+psipstar-bsubx*Bx/w

  syf = sy*vxt-bsuby*Bx/w

  szf = sz*vxt-bsubz*Bx/w

!!$ [psi^6 tau] vtilde^i +p* v^i - alp b^0 B^i/w
  tauf = tau*vxt + psipstar*velm - ab0*Bx/w

!!$ [psi^6 (B^k vtilde^i - B^i vtilde^k)]
  bxf = 0.0
  byf = By * vxt - Bx*vyt
  bzf = Bz * vxt - Bx*vzt

end subroutine num_x_fluxM

