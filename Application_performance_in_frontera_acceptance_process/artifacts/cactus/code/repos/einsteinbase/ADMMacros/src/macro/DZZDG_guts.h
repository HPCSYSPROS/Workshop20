/*@@
  @header   DZZDG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the (first and) second derivatives of the 
  physical metric with respect to z

  The macro is defined in terms of standard variables in
  @seefile DZZDG_declare.h

  The macro uses @seefile DZDG_guts.h and @seefile DZDG_declare.h
  @enddesc
@@*/

#ifndef DZZDG_GUTS
#define DZZDG_GUTS

#include "DZDG_guts.h"

#ifdef FCODE 

#include "ADM_Derivative.h"

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DZZDG_FAC  = 0

      if (local_spatial_order.eq.2) then
        DZZDG_DZZDGXX = 2*DZDCG_DZDCGXX*DZDG_FAC + DZZDG_FAC*DZDG_GXX \
        + DZDG_PSI4*ADM_DZZ_2(gxx,i,j,k)

        DZZDG_DZZDGXY = 2*DZDCG_DZDCGXY*DZDG_FAC + DZZDG_FAC*DZDG_GXY \
        + DZDG_PSI4*ADM_DZZ_2(gxy,i,j,k)

        DZZDG_DZZDGXZ = 2*DZDCG_DZDCGXZ*DZDG_FAC + DZZDG_FAC*DZDG_GXZ \
        + DZDG_PSI4*ADM_DZZ_2(gxz,i,j,k)

        DZZDG_DZZDGYY = 2*DZDCG_DZDCGYY*DZDG_FAC + DZZDG_FAC*DZDG_GYY \
        + DZDG_PSI4*ADM_DZZ_2(gyy,i,j,k)

        DZZDG_DZZDGYZ = 2*DZDCG_DZDCGYZ*DZDG_FAC + DZZDG_FAC*DZDG_GYZ \
        + DZDG_PSI4*ADM_DZZ_2(gyz,i,j,k)

        DZZDG_DZZDGZZ = 2*DZDCG_DZDCGZZ*DZDG_FAC + DZZDG_FAC*DZDG_GZZ \
        + DZDG_PSI4*ADM_DZZ_2(gzz,i,j,k)
      else
        DZZDG_DZZDGXX = 2*DZDCG_DZDCGXX*DZDG_FAC + DZZDG_FAC*DZDG_GXX \
        + DZDG_PSI4*ADM_DZZ_4(gxx,i,j,k)

        DZZDG_DZZDGXY = 2*DZDCG_DZDCGXY*DZDG_FAC + DZZDG_FAC*DZDG_GXY \
        + DZDG_PSI4*ADM_DZZ_4(gxy,i,j,k)

        DZZDG_DZZDGXZ = 2*DZDCG_DZDCGXZ*DZDG_FAC + DZZDG_FAC*DZDG_GXZ \
        + DZDG_PSI4*ADM_DZZ_4(gxz,i,j,k)

        DZZDG_DZZDGYY = 2*DZDCG_DZDCGYY*DZDG_FAC + DZZDG_FAC*DZDG_GYY \
        + DZDG_PSI4*ADM_DZZ_4(gyy,i,j,k)

        DZZDG_DZZDGYZ = 2*DZDCG_DZDCGYZ*DZDG_FAC + DZZDG_FAC*DZDG_GYZ \
        + DZDG_PSI4*ADM_DZZ_4(gyz,i,j,k)

        DZZDG_DZZDGZZ = 2*DZDCG_DZDCGZZ*DZDG_FAC + DZZDG_FAC*DZDG_GZZ \
        + DZDG_PSI4*ADM_DZZ_4(gzz,i,j,k)
      end if
#endif

#ifdef CCODE

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DZZDG_FAC   = 0;

      /* Now calculate the second deriatives */
      DZZDG_DZZDGXX = 2*DZDCG_DZDCGXX*DZDG_FAC+DZZDG_FAC*DZDG_GXX+DZDG_PSI4
                *DZZDG_OODZ2*(DZDCG_GXX_KP-2*DZDG_GXX+DZDCG_GXX_KM);

      DZZDG_DZZDGXY = 2*DZDCG_DZDCGXY*DZDG_FAC+DZZDG_FAC*DZDG_GXY+DZDG_PSI4
                *DZZDG_OODZ2*(DZDCG_GXY_KP-2*DZDG_GXY+DZDCG_GXY_KM);

      DZZDG_DZZDGXZ = 2*DZDCG_DZDCGXZ*DZDG_FAC+DZZDG_FAC*DZDG_GXZ+DZDG_PSI4
                *DZZDG_OODZ2*(DZDCG_GXZ_KP-2*DZDG_GXZ+DZDCG_GXZ_KM);

      DZZDG_DZZDGYY = 2*DZDCG_DZDCGYY*DZDG_FAC+DZZDG_FAC*DZDG_GYY+DZDG_PSI4
                *DZZDG_OODZ2*(DZDCG_GYY_KP-2*DZDG_GYY+DZDCG_GYY_KM);

      DZZDG_DZZDGYZ = 2*DZDCG_DZDCGYZ*DZDG_FAC+DZZDG_FAC*DZDG_GYZ+DZDG_PSI4
                *DZZDG_OODZ2*(DZDCG_GYZ_KP-2*DZDG_GYZ+DZDCG_GYZ_KM);

      DZZDG_DZZDGZZ = 2*DZDCG_DZDCGZZ*DZDG_FAC+DZZDG_FAC*DZDG_GZZ+DZDG_PSI4
                *DZZDG_OODZ2*(DZDCG_GZZ_KP-2*DZDG_GZZ+DZDCG_GZZ_KM);

#endif

#endif
