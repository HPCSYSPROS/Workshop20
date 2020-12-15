/*@@
  @header   DXZDG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the (first and) second derivatives of the 
  physical metric with respect to x,z

  The macro is defined in terms of standard variables in
  @seefile DXZDG_declare.h

  The macro uses @seefile DXDG_guts.h , @seefile DZDG_guts.h and 
  @seefile DXDG_declare.h , @seefile DZDG_declare.h
  @enddesc
@@*/

#ifndef DXZDG_GUTS
#define DXZDG_GUTS

#include "DXDG_guts.h"
#include "DZDG_guts.h"

#ifdef FCODE 

#include "ADM_Derivative.h"

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DXZDG_FAC  = 0

      /* Now calculate the second deriatives */
      if (local_spatial_order.eq.2) then 
        DXZDG_DXZDGXX = DZDCG_DZDCGXX*DXDG_FAC + DXDCG_DXDCGXX*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GXX + DXDG_PSI4*ADM_DXZ_2(gxx,i,j,k)

        DXZDG_DXZDGXY = DZDCG_DZDCGXY*DXDG_FAC + DXDCG_DXDCGXY*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GXY + DXDG_PSI4*ADM_DXZ_2(gxy,i,j,k)

        DXZDG_DXZDGXZ = DZDCG_DZDCGXZ*DXDG_FAC + DXDCG_DXDCGXZ*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GXZ + DXDG_PSI4*ADM_DXZ_2(gxz,i,j,k)

        DXZDG_DXZDGYY = DZDCG_DZDCGYY*DXDG_FAC + DXDCG_DXDCGYY*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GYY + DXDG_PSI4*ADM_DXZ_2(gyy,i,j,k)

        DXZDG_DXZDGYZ = DZDCG_DZDCGYZ*DXDG_FAC + DXDCG_DXDCGYZ*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GYZ + DXDG_PSI4*ADM_DXZ_2(gyz,i,j,k)

        DXZDG_DXZDGZZ = DZDCG_DZDCGZZ*DXDG_FAC + DXDCG_DXDCGZZ*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GZZ + DXDG_PSI4*ADM_DXZ_2(gzz,i,j,k)

      else

        DXZDG_DXZDGXX = DZDCG_DZDCGXX*DXDG_FAC + DXDCG_DXDCGXX*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GXX + DXDG_PSI4*ADM_DXZ_4(gxx,i,j,k)

        DXZDG_DXZDGXY = DZDCG_DZDCGXY*DXDG_FAC + DXDCG_DXDCGXY*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GXY + DXDG_PSI4*ADM_DXZ_4(gxy,i,j,k)

        DXZDG_DXZDGXZ = DZDCG_DZDCGXZ*DXDG_FAC + DXDCG_DXDCGXZ*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GXZ + DXDG_PSI4*ADM_DXZ_4(gxz,i,j,k)

        DXZDG_DXZDGYY = DZDCG_DZDCGYY*DXDG_FAC + DXDCG_DXDCGYY*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GYY + DXDG_PSI4*ADM_DXZ_4(gyy,i,j,k)

        DXZDG_DXZDGYZ = DZDCG_DZDCGYZ*DXDG_FAC + DXDCG_DXDCGYZ*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GYZ + DXDG_PSI4*ADM_DXZ_4(gyz,i,j,k)

        DXZDG_DXZDGZZ = DZDCG_DZDCGZZ*DXDG_FAC + DXDCG_DXDCGZZ*DZDG_FAC \
                      + DXZDG_FAC*DXDG_GZZ + DXDG_PSI4*ADM_DXZ_4(gzz,i,j,k)
      end if

#endif

#ifdef CCODE

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DXZDG_FAC   = 0;

      /* Now calculate the second deriatives */
      DXZDG_DXZDGXX = DZDCG_DZDCGXX*DXDG_FAC+DXDCG_DXDCGXX*DZDG_FAC+DXZDG_FAC*DXDG_GXX+
        DXDG_PSI4*DXZDG_OO4DXDZ*(DXZDG_GXX_IPKP-DXZDG_GXX_IPKM-DXZDG_GXX_IMKP+
        DXZDG_GXX_IMKM);

      DXZDG_DXZDGXY = DZDCG_DZDCGXY*DXDG_FAC+DXDCG_DXDCGXY*DZDG_FAC+DXZDG_FAC*DXDG_GXY+
        DXDG_PSI4*DXZDG_OO4DXDZ*(DXZDG_GXY_IPKP-DXZDG_GXY_IPKM-DXZDG_GXY_IMKP+
        DXZDG_GXY_IMKM);

      DXZDG_DXZDGXZ = DZDCG_DZDCGXZ*DXDG_FAC+DXDCG_DXDCGXZ*DZDG_FAC+DXZDG_FAC*DXDG_GXZ+
        DXDG_PSI4*DXZDG_OO4DXDZ*(DXZDG_GXZ_IPKP-DXZDG_GXZ_IPKM-DXZDG_GXZ_IMKP+
        DXZDG_GXZ_IMKM);

      DXZDG_DXZDGYY = DZDCG_DZDCGYY*DXDG_FAC+DXDCG_DXDCGYY*DZDG_FAC+DXZDG_FAC*DXDG_GYY+
        DXDG_PSI4*DXZDG_OO4DXDZ*(DXZDG_GYY_IPKP-DXZDG_GYY_IPKM-DXZDG_GYY_IMKP+
        DXZDG_GYY_IMKM);

      DXZDG_DXZDGYZ = DZDCG_DZDCGYZ*DXDG_FAC+DXDCG_DXDCGYZ*DZDG_FAC+DXZDG_FAC*DXDG_GYZ+
        DXDG_PSI4*DXZDG_OO4DXDZ*(DXZDG_GYZ_IPKP-DXZDG_GYZ_IPKM-DXZDG_GYZ_IMKP+
        DXZDG_GYZ_IMKM);

      DXZDG_DXZDGZZ = DZDCG_DZDCGZZ*DXDG_FAC+DXDCG_DXDCGZZ*DZDG_FAC+DXZDG_FAC*DXDG_GZZ+
        DXDG_PSI4*DXZDG_OO4DXDZ*(DXZDG_GZZ_IPKP-DXZDG_GZZ_IPKM-DXZDG_GZZ_IMKP+
        DXZDG_GZZ_IMKM);

#endif

#endif
