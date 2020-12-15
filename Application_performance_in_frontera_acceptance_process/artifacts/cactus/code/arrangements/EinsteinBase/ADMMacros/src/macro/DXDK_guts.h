/*@@
  @header   DXDK_guts.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the first derivatives of the 
  extrinsic curvature with respect to x
  @enddesc
@@*/

#ifndef DXDK_GUTS
#define DXDK_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"

      if (local_spatial_order.eq.2) then
        DXDK_DXDKXX = ADM_DX_2(kxx,i,j,k)
        DXDK_DXDKXY = ADM_DX_2(kxy,i,j,k)
        DXDK_DXDKXZ = ADM_DX_2(kxz,i,j,k)
        DXDK_DXDKYY = ADM_DX_2(kyy,i,j,k)
        DXDK_DXDKYZ = ADM_DX_2(kyz,i,j,k)
        DXDK_DXDKZZ = ADM_DX_2(kzz,i,j,k)
      else
        DXDK_DXDKXX = ADM_DX_4(kxx,i,j,k)
        DXDK_DXDKXY = ADM_DX_4(kxy,i,j,k)
        DXDK_DXDKXZ = ADM_DX_4(kxz,i,j,k)
        DXDK_DXDKYY = ADM_DX_4(kyy,i,j,k)
        DXDK_DXDKYZ = ADM_DX_4(kyz,i,j,k)
        DXDK_DXDKZZ = ADM_DX_4(kzz,i,j,k)
      end if

#endif

#ifdef CCODE

      DXDK_OO2DX = 1/(2*cctkGH->cctk_delta_space[0]);
    
      DXDK_DXDKXX = DXDK_OO2DX*(DXDK_KXX_IP - DXDK_KXX_IM);
      DXDK_DXDKXY = DXDK_OO2DX*(DXDK_KXY_IP - DXDK_KXY_IM);
      DXDK_DXDKXZ = DXDK_OO2DX*(DXDK_KXZ_IP - DXDK_KXZ_IM);
      DXDK_DXDKYY = DXDK_OO2DX*(DXDK_KYY_IP - DXDK_KYY_IM);
      DXDK_DXDKYZ = DXDK_OO2DX*(DXDK_KYZ_IP - DXDK_KYZ_IM);
      DXDK_DXDKZZ = DXDK_OO2DX*(DXDK_KZZ_IP - DXDK_KZZ_IM);

#endif

#endif
