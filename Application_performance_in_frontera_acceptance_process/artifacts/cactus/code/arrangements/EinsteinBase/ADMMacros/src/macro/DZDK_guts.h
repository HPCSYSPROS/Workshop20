/*@@
  @header   DZDG_guts.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the first derivatives of the 
  extrinsic with respect to z
  @enddesc
@@*/

#ifndef DZDK_GUTS
#define DZDK_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"
  
      if (local_spatial_order.eq.2) then   
        DZDK_DZDKXX = ADM_DZ_2(kxx,i,j,k)
        DZDK_DZDKXY = ADM_DZ_2(kxy,i,j,k)
        DZDK_DZDKXZ = ADM_DZ_2(kxz,i,j,k)
        DZDK_DZDKYY = ADM_DZ_2(kyy,i,j,k)
        DZDK_DZDKYZ = ADM_DZ_2(kyz,i,j,k)
        DZDK_DZDKZZ = ADM_DZ_2(kzz,i,j,k)
      else
        DZDK_DZDKXX = ADM_DZ_4(kxx,i,j,k)
        DZDK_DZDKXY = ADM_DZ_4(kxy,i,j,k)
        DZDK_DZDKXZ = ADM_DZ_4(kxz,i,j,k)
        DZDK_DZDKYY = ADM_DZ_4(kyy,i,j,k)
        DZDK_DZDKYZ = ADM_DZ_4(kyz,i,j,k)
        DZDK_DZDKZZ = ADM_DZ_4(kzz,i,j,k)
      end if
#endif

#ifdef CCODE

      DZDK_OO2DZ = 1/(2*cctkGH->cctk_delta_space[2]);
    
      DZDK_DZDKXX = DZDK_OO2DZ*(DZDK_KXX_KP - DZDK_KXX_KM);
      DZDK_DZDKXY = DZDK_OO2DZ*(DZDK_KXY_KP - DZDK_KXY_KM);
      DZDK_DZDKXZ = DZDK_OO2DZ*(DZDK_KXZ_KP - DZDK_KXZ_KM);
      DZDK_DZDKYY = DZDK_OO2DZ*(DZDK_KYY_KP - DZDK_KYY_KM);
      DZDK_DZDKYZ = DZDK_OO2DZ*(DZDK_KYZ_KP - DZDK_KYZ_KM);
      DZDK_DZDKZZ = DZDK_OO2DZ*(DZDK_KZZ_KP - DZDK_KZZ_KM);

#endif

#endif
