/*@@
  @header   DXDCG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the first derivatives of the 
  conformal metric with respect to x

  The macro is defined in terms of standard variables in

  Requires: lower conformal metric at i+1,i-1 ; dx0

  Provides: derivative of lower conformal metric wrt x

  @seefile DXDCG_declare.h
  @enddesc
@@*/

#ifndef DXDCG_GUTS
#define DXDCG_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"

      if (local_spatial_order.eq.2) then
        DXDCG_DXDCGXX = ADM_DX_2(gxx,i,j,k)
        DXDCG_DXDCGXY = ADM_DX_2(gxy,i,j,k)
        DXDCG_DXDCGXZ = ADM_DX_2(gxz,i,j,k)
        DXDCG_DXDCGYY = ADM_DX_2(gyy,i,j,k)
        DXDCG_DXDCGYZ = ADM_DX_2(gyz,i,j,k)
        DXDCG_DXDCGZZ = ADM_DX_2(gzz,i,j,k)
      else
        DXDCG_DXDCGXX = ADM_DX_4(gxx,i,j,k)
        DXDCG_DXDCGXY = ADM_DX_4(gxy,i,j,k)
        DXDCG_DXDCGXZ = ADM_DX_4(gxz,i,j,k)
        DXDCG_DXDCGYY = ADM_DX_4(gyy,i,j,k)
        DXDCG_DXDCGYZ = ADM_DX_4(gyz,i,j,k)
        DXDCG_DXDCGZZ = ADM_DX_4(gzz,i,j,k)
      end if
#endif

#ifdef CCODE
 
      DXDCG_OO2DX = 1/(2*cctkGH->cctk_delta_space[0]);
    
      DXDCG_DXDCGXX = DXDCG_OO2DX*(DXDCG_GXX_IP - DXDCG_GXX_IM);
      DXDCG_DXDCGXY = DXDCG_OO2DX*(DXDCG_GXY_IP - DXDCG_GXY_IM);
      DXDCG_DXDCGXZ = DXDCG_OO2DX*(DXDCG_GXZ_IP - DXDCG_GXZ_IM);
      DXDCG_DXDCGYY = DXDCG_OO2DX*(DXDCG_GYY_IP - DXDCG_GYY_IM);
      DXDCG_DXDCGYZ = DXDCG_OO2DX*(DXDCG_GYZ_IP - DXDCG_GYZ_IM);
      DXDCG_DXDCGZZ = DXDCG_OO2DX*(DXDCG_GZZ_IP - DXDCG_GZZ_IM);

#endif

#endif
