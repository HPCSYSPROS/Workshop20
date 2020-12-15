/*@@
  @header   DA_guts.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate all first spatial derivative of lapse
  @enddesc
@@*/

#ifndef DA_GUTS
#define DA_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"

      if (local_spatial_order.eq.2) then
        DA_DXDA = ADM_DX_2(alp,i,j,k)
        DA_DYDA = ADM_DY_2(alp,i,j,k)
        DA_DZDA = ADM_DZ_2(alp,i,j,k)
      else
        DA_DXDA = ADM_DX_4(alp,i,j,k)
        DA_DYDA = ADM_DY_4(alp,i,j,k)
        DA_DZDA = ADM_DZ_4(alp,i,j,k)
      end if

#endif

#ifdef CCODE

      DA_OO2DX = 1/(2*DA_DX);
      DA_OO2DY = 1/(2*DA_DY);
      DA_OO2DZ = 1/(2*DA_DZ);

      DA_DXDA = DA_OO2DX*(DA_A_IP - DA_A_IM);
      DA_DYDA = DA_OO2DY*(DA_A_JP - DA_A_JM);
      DA_DZDA = DA_OO2DZ*(DA_A_KP - DA_A_KM);

#endif

#endif

