/*@@
  @header   DDA_guts.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate all second spatial derivative of lapse
  @enddesc
@@*/

#ifndef DDA_GUTS
#define DDA_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"

      if (local_spatial_order.eq.2) then
        DDA_DXXDA = ADM_DXX_2(alp,i,j,k)
        DDA_DXYDA = ADM_DXY_2(alp,i,j,k)
        DDA_DXZDA = ADM_DXZ_2(alp,i,j,k)
        DDA_DYYDA = ADM_DYY_2(alp,i,j,k)
        DDA_DYZDA = ADM_DYZ_2(alp,i,j,k)
        DDA_DZZDA = ADM_DZZ_2(alp,i,j,k)
      else
        DDA_DXXDA = ADM_DXX_4(alp,i,j,k)
        DDA_DXYDA = ADM_DXY_4(alp,i,j,k)
        DDA_DXZDA = ADM_DXZ_4(alp,i,j,k)
        DDA_DYYDA = ADM_DYY_4(alp,i,j,k)
        DDA_DYZDA = ADM_DYZ_4(alp,i,j,k)
        DDA_DZZDA = ADM_DZZ_4(alp,i,j,k)
      end if

#endif

#ifdef CCODE

      DDA_OODX2 = 1/(DDA_DX*DDA_DX);
      DDA_OODY2 = 1/(DDA_DY*DDA_DY);
      DDA_OODZ2 = 1/(DDA_DZ*DDA_DZ);
      DDA_OO4DXDY = 1/(4*DDA_DX*DDA_DY);
      DDA_OO4DXDZ = 1/(4*DDA_DX*DDA_DZ);
      DDA_OO4DYDZ = 1/(4*DDA_DY*DDA_DZ);

      DDA_DXXDA = DDA_OODX2*(DDA_A_IP  - 2*DDA_A + DDA_A_IM);
      DDA_DYYDA = DDA_OODY2*(DDA_A_JP  - 2*DDA_A + DDA_A_JM);
      DDA_DZZDA = DDA_OODZ2*(DDA_A_KP  - 2*DDA_A + DDA_A_KM);

      DDA_DXYDA = DDA_OO4DXDY*(DDA_A_IPJP-DDA_A_IPJM-DDA_A_IMJP+DDA_A_IMJM);
      DDA_DXZDA = DDA_OO4DXDZ*(DDA_A_IPKP-DDA_A_IPKM-DDA_A_IMKP+DDA_A_IMKM);
      DDA_DYZDA = DDA_OO4DYDZ*(DDA_A_JPKP-DDA_A_JPKM-DDA_A_JMKP+DDA_A_JMKM);

#endif

#endif

