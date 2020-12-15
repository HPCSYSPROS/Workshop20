/*@@
  @header   DYYDG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the (first and) second derivatives of the 
  physical metric with respect to y

  The macro is defined in terms of standard variables in
  @seefile DYYDG_declare.h

  The macro uses @seefile DXDG_guts.h and @seefile DXDG_declare.h
  @enddesc
@@*/

#ifndef DYYDG_GUTS
#define DYYDG_GUTS

#include "DYDG_guts.h"

#ifdef FCODE 

#include "ADM_Derivative.h"

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DYYDG_FAC  = 0

      /* Now calculate the second deriatives */
      if (local_spatial_order.eq.2) then 
        DYYDG_DYYDGXX = 2*DYDCG_DYDCGXX*DYDG_FAC + DYYDG_FAC*DYDG_GXX \
                      + DYDG_PSI4*ADM_DYY_2(gxx,i,j,k)

        DYYDG_DYYDGXY = 2*DYDCG_DYDCGXY*DYDG_FAC + DYYDG_FAC*DYDG_GXY \
                      + DYDG_PSI4*ADM_DYY_2(gxy,i,j,k)

        DYYDG_DYYDGXZ = 2*DYDCG_DYDCGXZ*DYDG_FAC + DYYDG_FAC*DYDG_GXZ \
                      + DYDG_PSI4*ADM_DYY_2(gxz,i,j,k)

        DYYDG_DYYDGYY = 2*DYDCG_DYDCGYY*DYDG_FAC + DYYDG_FAC*DYDG_GYY \
                      + DYDG_PSI4*ADM_DYY_2(gyy,i,j,k)

        DYYDG_DYYDGYZ = 2*DYDCG_DYDCGYZ*DYDG_FAC + DYYDG_FAC*DYDG_GYZ \
                      + DYDG_PSI4*ADM_DYY_2(gyz,i,j,k)

        DYYDG_DYYDGZZ = 2*DYDCG_DYDCGZZ*DYDG_FAC + DYYDG_FAC*DYDG_GZZ \
                      + DYDG_PSI4*ADM_DYY_2(gzz,i,j,k)
      else
        DYYDG_DYYDGXX = 2*DYDCG_DYDCGXX*DYDG_FAC + DYYDG_FAC*DYDG_GXX \
                      + DYDG_PSI4*ADM_DYY_4(gxx,i,j,k)

        DYYDG_DYYDGXY = 2*DYDCG_DYDCGXY*DYDG_FAC + DYYDG_FAC*DYDG_GXY \
                      + DYDG_PSI4*ADM_DYY_4(gxy,i,j,k)

        DYYDG_DYYDGXZ = 2*DYDCG_DYDCGXZ*DYDG_FAC + DYYDG_FAC*DYDG_GXZ \
                      + DYDG_PSI4*ADM_DYY_4(gxz,i,j,k)

        DYYDG_DYYDGYY = 2*DYDCG_DYDCGYY*DYDG_FAC + DYYDG_FAC*DYDG_GYY \
                      + DYDG_PSI4*ADM_DYY_4(gyy,i,j,k)

        DYYDG_DYYDGYZ = 2*DYDCG_DYDCGYZ*DYDG_FAC + DYYDG_FAC*DYDG_GYZ \
                      + DYDG_PSI4*ADM_DYY_4(gyz,i,j,k)

        DYYDG_DYYDGZZ = 2*DYDCG_DYDCGZZ*DYDG_FAC + DYYDG_FAC*DYDG_GZZ \
                      + DYDG_PSI4*ADM_DYY_4(gzz,i,j,k)
      end if

#endif

#ifdef CCODE

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DYYDG_FAC   = 0;

      /* Now calculate the second deriatives */
      DYYDG_DYYDGXX = 2*DYDCG_DYDCGXX*DYDG_FAC+DYYDG_FAC*DYDG_GXX+DYDG_PSI4
      *DYYDG_OODY2*(DYDCG_GXX_JP-2*DYDG_GXX+DYDCG_GXX_JM);

      DYYDG_DYYDGXY = 2*DYDCG_DYDCGXY*DYDG_FAC+DYYDG_FAC*DYDG_GXY+DYDG_PSI4
      *DYYDG_OODY2*(DYDCG_GXY_JP-2*DYDG_GXY+DYDCG_GXY_JM);

      DYYDG_DYYDGXZ = 2*DYDCG_DYDCGXZ*DYDG_FAC+DYYDG_FAC*DYDG_GXZ+DYDG_PSI4
      *DYYDG_OODY2*(DYDCG_GXZ_JP-2*DYDG_GXZ+DYDCG_GXZ_JM);

      DYYDG_DYYDGYY = 2*DYDCG_DYDCGYY*DYDG_FAC+DYYDG_FAC*DYDG_GYY+DYDG_PSI4
      *DYYDG_OODY2*(DYDCG_GYY_JP-2*DYDG_GYY+DYDCG_GYY_JM);

      DYYDG_DYYDGYZ = 2*DYDCG_DYDCGYZ*DYDG_FAC+DYYDG_FAC*DYDG_GYZ+DYDG_PSI4
      *DYYDG_OODY2*(DYDCG_GYZ_JP-2*DYDG_GYZ+DYDCG_GYZ_JM);

      DYYDG_DYYDGZZ = 2*DYDCG_DYDCGZZ*DYDG_FAC+DYYDG_FAC*DYDG_GZZ+DYDG_PSI4
      *DYYDG_OODY2*(DYDCG_GZZ_JP-2*DYDG_GZZ+DYDCG_GZZ_JM);

#endif

#endif
