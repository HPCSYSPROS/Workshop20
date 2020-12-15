/*@@
  @header   DXXDG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the (first and) second derivatives of the 
  physical metric with respect to x

  The macro is defined in terms of standard variables in
  @seefile DXXDG_declare.h

  The macro uses @seefile DXDG_guts.h and @seefile DXDG_declare.h
  @enddesc
@@*/

#ifndef DXXDG_GUTS
#define DXXDG_GUTS

#include "DXDG_guts.h"

#ifdef FCODE 

#include "ADM_Derivative.h"

      /* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
      DXXDG_FAC  = 0
      
      /* Now calculate the second derivatives */

      if (local_spatial_order.eq.2) then
        DXXDG_DXXDGXX = 2*DXDCG_DXDCGXX*DXDG_FAC + DXXDG_FAC*DXDG_GXX \
                      + DXDG_PSI4*ADM_DXX_2(gxx,i,j,k)

        DXXDG_DXXDGXY = 2*DXDCG_DXDCGXY*DXDG_FAC + DXXDG_FAC*DXDG_GXY \
                      + DXDG_PSI4*ADM_DXX_2(gxy,i,j,k)

        DXXDG_DXXDGXZ = 2*DXDCG_DXDCGXZ*DXDG_FAC + DXXDG_FAC*DXDG_GXZ \
                      + DXDG_PSI4*ADM_DXX_2(gxz,i,j,k)

        DXXDG_DXXDGYY = 2*DXDCG_DXDCGYY*DXDG_FAC + DXXDG_FAC*DXDG_GYY \
                      + DXDG_PSI4*ADM_DXX_2(gyy,i,j,k)

        DXXDG_DXXDGYZ = 2*DXDCG_DXDCGYZ*DXDG_FAC + DXXDG_FAC*DXDG_GYZ \
                      + DXDG_PSI4*ADM_DXX_2(gyz,i,j,k)

        DXXDG_DXXDGZZ = 2*DXDCG_DXDCGZZ*DXDG_FAC + DXXDG_FAC*DXDG_GZZ \
                      + DXDG_PSI4*ADM_DXX_2(gzz,i,j,k)
      else
        DXXDG_DXXDGXX = 2*DXDCG_DXDCGXX*DXDG_FAC + DXXDG_FAC*DXDG_GXX \
                      + DXDG_PSI4*ADM_DXX_4(gxx,i,j,k)

        DXXDG_DXXDGXY = 2*DXDCG_DXDCGXY*DXDG_FAC + DXXDG_FAC*DXDG_GXY \
                      + DXDG_PSI4*ADM_DXX_4(gxy,i,j,k)

        DXXDG_DXXDGXZ = 2*DXDCG_DXDCGXZ*DXDG_FAC + DXXDG_FAC*DXDG_GXZ \
                      + DXDG_PSI4*ADM_DXX_4(gxz,i,j,k)

        DXXDG_DXXDGYY = 2*DXDCG_DXDCGYY*DXDG_FAC + DXXDG_FAC*DXDG_GYY \
                      + DXDG_PSI4*ADM_DXX_4(gyy,i,j,k)

        DXXDG_DXXDGYZ = 2*DXDCG_DXDCGYZ*DXDG_FAC + DXXDG_FAC*DXDG_GYZ \
                      + DXDG_PSI4*ADM_DXX_4(gyz,i,j,k)

        DXXDG_DXXDGZZ = 2*DXDCG_DXDCGZZ*DXDG_FAC + DXXDG_FAC*DXDG_GZZ \
                      + DXDG_PSI4*ADM_DXX_4(gzz,i,j,k)
      end if

#endif

#ifdef CCODE

/* Factor involving 2nd derivative of conformal factor, nowadays zero */ 
DXXDG_FAC   = 0;

/* Now calculate the second derivatives */
      DXXDG_DXXDGXX = 2*DXDCG_DXDCGXX*DXDG_FAC+DXXDG_FAC*DXDG_GXX+DXDG_PSI4
        *DXXDG_OODX2*(DXDCG_GXX_IP-2*DXDG_GXX+DXDCG_GXX_IM);

      DXXDG_DXXDGXY = 2*DXDCG_DXDCGXY*DXDG_FAC+DXXDG_FAC*DXDG_GXY+DXDG_PSI4
        *DXXDG_OODX2*(DXDCG_GXY_IP-2*DXDG_GXY+DXDCG_GXY_IM);

      DXXDG_DXXDGXZ = 2*DXDCG_DXDCGXZ*DXDG_FAC+DXXDG_FAC*DXDG_GXZ+DXDG_PSI4
        *DXXDG_OODX2*(DXDCG_GXZ_IP-2*DXDG_GXZ+DXDCG_GXZ_IM);

      DXXDG_DXXDGYY = 2*DXDCG_DXDCGYY*DXDG_FAC+DXXDG_FAC*DXDG_GYY+DXDG_PSI4
        *DXXDG_OODX2*(DXDCG_GYY_IP-2*DXDG_GYY+DXDCG_GYY_IM);

      DXXDG_DXXDGYZ = 2*DXDCG_DXDCGYZ*DXDG_FAC+DXXDG_FAC*DXDG_GYZ+DXDG_PSI4
        *DXXDG_OODX2*(DXDCG_GYZ_IP-2*DXDG_GYZ+DXDCG_GYZ_IM);

      DXXDG_DXXDGZZ = 2*DXDCG_DXDCGZZ*DXDG_FAC+DXXDG_FAC*DXDG_GZZ+DXDG_PSI4
        *DXXDG_OODX2*(DXDCG_GZZ_IP-2*DXDG_GZZ+DXDCG_GZZ_IM);
 

#endif

#endif
