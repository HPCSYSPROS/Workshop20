/*@@
  @header   DZDCG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the first derivatives of the 
  conformal metric with respect to z

  The macro is defined in terms of standard variables in

  Requires: lower conformal metric at k+1,k-1 ; dz0

  Provides: derivative of lower conformal metric wrt z

  @seefile DZDCG_declare.h
  @enddesc
@@*/

#ifndef DZDCG_GUTS
#define DZDCG_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"
      if (local_spatial_order.eq.2) then     
        DZDCG_DZDCGXX = ADM_DZ_2(gxx,i,j,k)
        DZDCG_DZDCGXY = ADM_DZ_2(gxy,i,j,k)
        DZDCG_DZDCGXZ = ADM_DZ_2(gxz,i,j,k)
        DZDCG_DZDCGYY = ADM_DZ_2(gyy,i,j,k)
        DZDCG_DZDCGYZ = ADM_DZ_2(gyz,i,j,k)
        DZDCG_DZDCGZZ = ADM_DZ_2(gzz,i,j,k)
      else
        DZDCG_DZDCGXX = ADM_DZ_4(gxx,i,j,k)
        DZDCG_DZDCGXY = ADM_DZ_4(gxy,i,j,k)
        DZDCG_DZDCGXZ = ADM_DZ_4(gxz,i,j,k)
        DZDCG_DZDCGYY = ADM_DZ_4(gyy,i,j,k)
        DZDCG_DZDCGYZ = ADM_DZ_4(gyz,i,j,k)
        DZDCG_DZDCGZZ = ADM_DZ_4(gzz,i,j,k)
      end if
#endif

#ifdef CCODE

      DZDCG_OO2DZ = 1/(2*cctkGH->cctk_delta_space[2]);
    
      DZDCG_DZDCGXX = DZDCG_OO2DZ*(DZDCG_GXX_KP - DZDCG_GXX_KM);
      DZDCG_DZDCGXY = DZDCG_OO2DZ*(DZDCG_GXY_KP - DZDCG_GXY_KM);
      DZDCG_DZDCGXZ = DZDCG_OO2DZ*(DZDCG_GXZ_KP - DZDCG_GXZ_KM);
      DZDCG_DZDCGYY = DZDCG_OO2DZ*(DZDCG_GYY_KP - DZDCG_GYY_KM);
      DZDCG_DZDCGYZ = DZDCG_OO2DZ*(DZDCG_GYZ_KP - DZDCG_GYZ_KM);
      DZDCG_DZDCGZZ = DZDCG_OO2DZ*(DZDCG_GZZ_KP - DZDCG_GZZ_KM);

#endif

#endif
