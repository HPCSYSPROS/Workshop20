/*@@
  @header   DYDG_guts.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the first derivatives of the 
  extrinsic curvature with respect to y
  @enddesc
@@*/

#ifndef DYDK_GUTS
#define DYDK_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"
  
      if (local_spatial_order.eq.2) then 
        DYDK_DYDKXX = ADM_DY_2(kxx,i,j,k)
        DYDK_DYDKXY = ADM_DY_2(kxy,i,j,k)
        DYDK_DYDKXZ = ADM_DY_2(kxz,i,j,k)
        DYDK_DYDKYY = ADM_DY_2(kyy,i,j,k)
        DYDK_DYDKYZ = ADM_DY_2(kyz,i,j,k)
        DYDK_DYDKZZ = ADM_DY_2(kzz,i,j,k)
      else
        DYDK_DYDKXX = ADM_DY_4(kxx,i,j,k)
        DYDK_DYDKXY = ADM_DY_4(kxy,i,j,k)
        DYDK_DYDKXZ = ADM_DY_4(kxz,i,j,k)
        DYDK_DYDKYY = ADM_DY_4(kyy,i,j,k)
        DYDK_DYDKYZ = ADM_DY_4(kyz,i,j,k)
        DYDK_DYDKZZ = ADM_DY_4(kzz,i,j,k)
      end if
#endif


#ifdef CCODE

      DYDK_OO2DY = 1/(2*cctkGH->cctk_delta_space[1]);
    
      DYDK_DYDKXX = DYDK_OO2DY*(DYDK_KXX_JP - DYDK_KXX_JM);
      DYDK_DYDKXY = DYDK_OO2DY*(DYDK_KXY_JP - DYDK_KXY_JM);
      DYDK_DYDKXZ = DYDK_OO2DY*(DYDK_KXZ_JP - DYDK_KXZ_JM);
      DYDK_DYDKYY = DYDK_OO2DY*(DYDK_KYY_JP - DYDK_KYY_JM);
      DYDK_DYDKYZ = DYDK_OO2DY*(DYDK_KYZ_JP - DYDK_KYZ_JM);
      DYDK_DYDKZZ = DYDK_OO2DY*(DYDK_KZZ_JP - DYDK_KZZ_JM);

#endif

#endif
