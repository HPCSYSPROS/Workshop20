/*@@
@file      UtilsM.F
@date      Aug 30, 2010
@author    Joshua Faber, Scott Noble, Bruno Mundim
@desc 
Utility functions for other thorns. 
@enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


subroutine calc_vlow_blow(gxx,gxy,gxz,gyy,gyz,gzz, &
     velx,vely,velz,Bvecx,Bvecy,Bvecz, &
     velxlow,velylow,velzlow,Bvecxlow,Bvecylow,Bveczlow, &
     Bdotv,b2,v2,w,bxlow,bylow,bzlow)
  
!!$ Calculates v_i (see Anton Eq. 5) and B_i (Bvecxlow)- undensitized!
!!$ calculates B^i v_i [Anton eq. 44] and b^2 [LHS of Anton eq. 46] 
!!$ Calculates w (Lorentz factor) as (1-v^i v_i)^{-1/2}
!!$ Calculates b_i (bxlow)

  CCTK_REAL :: gxx,gxy,gxz,gyy,gyz,gzz
  CCTK_REAL :: velx,vely,velz,Bvecx,Bvecy,Bvecz
  CCTK_REAL :: velxlow,velylow,velzlow
  CCTK_REAL :: Bvecxlow,Bvecylow,Bveczlow
  CCTK_REAL :: Bdotv,v2,w,b2,bxlow,bylow,bzlow

!!$ vel_i  = g_ij v^j
!!$ B_i = g_ij B^i

  velxlow = gxx*velx + gxy*vely + gxz*velz
  velylow = gxy*velx + gyy*vely + gyz*velz
  velzlow = gxz*velx + gyz*vely + gzz*velz
  Bvecxlow = gxx*Bvecx + gxy*Bvecy + gxz*Bvecz
  Bvecylow = gxy*Bvecx + gyy*Bvecy + gyz*Bvecz
  Bveczlow = gxz*Bvecx + gyz*Bvecy + gzz*Bvecz
 
!!$  B^i v_i (= b^0/u^0)
  Bdotv = velxlow*Bvecx+velylow*Bvecy+velzlow*Bvecz
  
!!$v^2 = v_i v^i; w=(1-v^2)^{-1/2}

  v2 = velxlow*velx + velylow*vely + velzlow*velz
  w = 1.d0/sqrt(1.d0-v2)

!!$b^2 = B^i B_i / w^2 + (b^0/u^0)^2

  b2=(Bvecx*Bvecxlow+Bvecy*Bvecylow+Bvecz*Bveczlow)/w**2+Bdotv**2

!!$ b_i = B_i/w +w*(B dot v)*v_i
  bxlow = Bvecxlow/w+w*Bdotv*velxlow
  bylow = Bvecylow/w+w*Bdotv*velylow
  bzlow = Bveczlow/w+w*Bdotv*velzlow

end subroutine calc_vlow_blow
