/*@@
  @header   DYZDG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the (first and) second derivatives of the 
  physical metric with respect to y,z

  The macro is defined in terms of standard variables in
  @seefile DYZDG_declare.h

  The macro uses @seefile DYDG_guts.h , @seefile DZDG_guts.h and 
  @seefile DYDG_declare.h , @seefile DZDG_declare.h
  @enddesc
@@*/

#ifndef DYZDG_GUTS
#define DYZDG_GUTS

#include "DYDG_guts.h"
#include "DZDG_guts.h"

#ifdef FCODE 

#include "ADM_Derivative.h"

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DYZDG_FAC  = 0
      
      /* Now calculate the second deriatives */
      if (local_spatial_order.eq.2) then
        DYZDG_DYZDGXX = DZDCG_DZDCGXX*DYDG_FAC + DYDCG_DYDCGXX*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GXX + DYDG_PSI4*ADM_DYZ_2(gxx,i,j,k)

        DYZDG_DYZDGXY = DZDCG_DZDCGXY*DYDG_FAC + DYDCG_DYDCGXY*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GXY + DYDG_PSI4*ADM_DYZ_2(gxy,i,j,k)

        DYZDG_DYZDGXZ = DZDCG_DZDCGXZ*DYDG_FAC + DYDCG_DYDCGXZ*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GXZ + DYDG_PSI4*ADM_DYZ_2(gxz,i,j,k)

        DYZDG_DYZDGYY = DZDCG_DZDCGYY*DYDG_FAC + DYDCG_DYDCGYY*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GYY + DYDG_PSI4*ADM_DYZ_2(gyy,i,j,k)

        DYZDG_DYZDGYZ = DZDCG_DZDCGYZ*DYDG_FAC + DYDCG_DYDCGYZ*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GYZ + DYDG_PSI4*ADM_DYZ_2(gyz,i,j,k)

        DYZDG_DYZDGZZ = DZDCG_DZDCGZZ*DYDG_FAC + DYDCG_DYDCGZZ*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GZZ + DYDG_PSI4*ADM_DYZ_2(gzz,i,j,k)
      else 
        DYZDG_DYZDGXX = DZDCG_DZDCGXX*DYDG_FAC + DYDCG_DYDCGXX*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GXX + DYDG_PSI4*ADM_DYZ_4(gxx,i,j,k)

        DYZDG_DYZDGXY = DZDCG_DZDCGXY*DYDG_FAC + DYDCG_DYDCGXY*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GXY + DYDG_PSI4*ADM_DYZ_4(gxy,i,j,k)

        DYZDG_DYZDGXZ = DZDCG_DZDCGXZ*DYDG_FAC + DYDCG_DYDCGXZ*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GXZ + DYDG_PSI4*ADM_DYZ_4(gxz,i,j,k)

        DYZDG_DYZDGYY = DZDCG_DZDCGYY*DYDG_FAC + DYDCG_DYDCGYY*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GYY + DYDG_PSI4*ADM_DYZ_4(gyy,i,j,k)

        DYZDG_DYZDGYZ = DZDCG_DZDCGYZ*DYDG_FAC + DYDCG_DYDCGYZ*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GYZ + DYDG_PSI4*ADM_DYZ_4(gyz,i,j,k)

        DYZDG_DYZDGZZ = DZDCG_DZDCGZZ*DYDG_FAC + DYDCG_DYDCGZZ*DZDG_FAC \
                      + DYZDG_FAC*DYDG_GZZ + DYDG_PSI4*ADM_DYZ_4(gzz,i,j,k)
      end if
#endif

#ifdef CCODE

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DYZDG_FAC   = 0;

      /* Now calculate the second deriatives */
      DYZDG_DYZDGXX = DZDCG_DZDCGXX*DYDG_FAC+DYDCG_DYDCGXX*DZDG_FAC+DYZDG_FAC*DYDG_GXX
         +DYDG_PSI4*DYZDG_OO4DYDZ*(DYZDG_GXX_JPKP-DYZDG_GXX_JPKM-DYZDG_GXX_JMKP
         +DYZDG_GXX_JMKM);

      DYZDG_DYZDGXY = DZDCG_DZDCGXY*DYDG_FAC+DYDCG_DYDCGXY*DZDG_FAC+DYZDG_FAC*DYDG_GXY
         +DYDG_PSI4*DYZDG_OO4DYDZ*(DYZDG_GXY_JPKP-DYZDG_GXY_JPKM-DYZDG_GXY_JMKP
         +DYZDG_GXY_JMKM);

      DYZDG_DYZDGXZ = DZDCG_DZDCGXZ*DYDG_FAC+DYDCG_DYDCGXZ*DZDG_FAC+DYZDG_FAC*DYDG_GXZ
         +DYDG_PSI4*DYZDG_OO4DYDZ*(DYZDG_GXZ_JPKP-DYZDG_GXZ_JPKM-DYZDG_GXZ_JMKP
         +DYZDG_GXZ_JMKM);

      DYZDG_DYZDGYY = DZDCG_DZDCGYY*DYDG_FAC+DYDCG_DYDCGYY*DZDG_FAC+DYZDG_FAC*DYDG_GYY
         +DYDG_PSI4*DYZDG_OO4DYDZ*(DYZDG_GYY_JPKP-DYZDG_GYY_JPKM-DYZDG_GYY_JMKP
         +DYZDG_GYY_JMKM);

      DYZDG_DYZDGYZ = DZDCG_DZDCGYZ*DYDG_FAC+DYDCG_DYDCGYZ*DZDG_FAC+DYZDG_FAC*DYDG_GYZ
         +DYDG_PSI4*DYZDG_OO4DYDZ*(DYZDG_GYZ_JPKP-DYZDG_GYZ_JPKM-DYZDG_GYZ_JMKP
         +DYZDG_GYZ_JMKM);

      DYZDG_DYZDGZZ = DZDCG_DZDCGZZ*DYDG_FAC+DYDCG_DYDCGZZ*DZDG_FAC+DYZDG_FAC*DYDG_GZZ
         +DYDG_PSI4*DYZDG_OO4DYDZ*(DYZDG_GZZ_JPKP-DYZDG_GZZ_JPKM-DYZDG_GZZ_JMKP
         +DYZDG_GZZ_JMKM);

#endif

#endif

