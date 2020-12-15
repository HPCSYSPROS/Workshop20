/*@@
@file      GRHydro_Flux.F90
@date      Sat Jan 26 01:36:57 2002
@author    Pedro Montero, Ian Hawke
@desc 
The routine to calculate the numerical flux function given a 
  specific state
  @enddesc 
  @@*/
  
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
 
 /*@@
   @routine    num_flux
   @date       Sat Jan 26 01:38:12 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   A general routine that works out the correct direction and returns
   the appropriate flux function
   @enddesc 
   @calls     
   @calledby   
   @history 
   Called from GR3D, original author Mark Miller.
   @endhistory 

@@*/
 
subroutine num_flux(xoffset,yoffset,zoffset,dens,sx,sy,sz,tau,&
     densf,sxf,syf,szf,tauf,velx,vely,velz,press,det,alp,beta)

  implicit none
  
  CCTK_REAL :: dens,sx,sy,sz,tau,velx,vely,velz,det
  CCTK_REAL :: densf,sxf,syf,szf,tauf
  CCTK_REAL :: alp,beta,press 
  integer :: xoffset,yoffset,zoffset

  densf = dens * ( (velx*xoffset + vely*yoffset + &
       velz*zoffset)  - beta / alp )
               
  sxf = sx * ( (velx*xoffset + vely*yoffset +&
       velz*zoffset)  - beta/alp) + &
       sqrt(det)*press*xoffset
               
  syf = sy * ( (velx*xoffset + vely*yoffset +&
       velz*zoffset)  - beta/alp) + &
       sqrt(det)*press*yoffset
              
  szf = sz * ( (velx*xoffset + vely*yoffset +&
       velz*zoffset)  - beta/alp ) + &
       sqrt(det)*press*zoffset
               
  tauf = tau * ( (velx*xoffset + vely*yoffset +&
       velz*zoffset)  - beta/alp ) + &
       sqrt(det)*press*(velx*xoffset +&
       vely*yoffset + velz*zoffset) 
 
end subroutine num_flux

 /*@@
   @routine    num_x_flux
   @date       Sat Jan 26 01:38:55 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   The numerical flux just along the x direction. Used for
   fluxes along the y,z direction by permuting the appropriate 
   arguments in the subroutine call.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Culled and altered from GR3D, original author Mark Miller.
   @endhistory 

@@*/


subroutine num_x_flux(dens,sx,sy,sz,tau,&
     densf,sxf,syf,szf,tauf,vel,press,det,alp,beta)

  implicit none

  CCTK_REAL :: dens,sx,sy,sz,tau,vel,det
  CCTK_REAL :: densf,sxf,syf,szf,tauf
  CCTK_REAL :: alp,beta,press 
  CCTK_REAL :: velmbetainvalp,sqrtdetpress

  velmbetainvalp = vel - beta / alp
  sqrtdetpress = sqrt(det)*press
  densf = dens * ( velmbetainvalp )

  sxf = sx * velmbetainvalp + sqrtdetpress

  syf = sy * velmbetainvalp

  szf = sz * velmbetainvalp

  tauf = tau * velmbetainvalp + sqrtdetpress*vel
 
end subroutine num_x_flux












