/*@@
  @header   LIEG_guts.h
  @date     July 2000
  @author   Gabrielle Allen, Miguel Alcubierre
  @desc

  Macro to calculate the Lie derivative of the lower 
  conformal metric.  Notice that the advection term
  can be calculated in several different ways.

  IMPORTANT:  THE C VERSION ONLY HAS CENTERED DIFFERENCES!

  @enddesc
@@*/

#ifndef LIEG_GUTS
#define LIEG_GUTS

#include "DB_guts.h"
#include "DG_guts.h"

#ifdef FCODE

#include "ADM_Derivative.h"

/*  Advection. */

      LIEG_LGXX = 0.0D0
      LIEG_LGYY = 0.0D0
      LIEG_LGZZ = 0.0D0
      LIEG_LGXY = 0.0D0
      LIEG_LGXZ = 0.0D0
      LIEG_LGYZ = 0.0D0

      /* Do the x-direction */

      if (admmacros_advectionx .eq. 0) then /* CENTER */

         LIEG_LGXX = LIEG_LGXX + DXDCG_DXDCGXX*LIEG_BX 
         LIEG_LGYY = LIEG_LGYY + DXDCG_DXDCGYY*LIEG_BX 
         LIEG_LGZZ = LIEG_LGZZ + DXDCG_DXDCGZZ*LIEG_BX 
         LIEG_LGXY = LIEG_LGXY + DXDCG_DXDCGXY*LIEG_BX
         LIEG_LGXZ = LIEG_LGXZ + DXDCG_DXDCGXZ*LIEG_BX 
         LIEG_LGYZ = LIEG_LGYZ + DXDCG_DXDCGYZ*LIEG_BX 

      else if (admmacros_advectionx .eq. 1) then /* UPWIND1 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DX_P1(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DX_P1(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DX_P1(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DX_P1(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DX_P1(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DX_P1(gyz,i,j,k)

      else if (admmacros_advectionx .eq. -1) then /* UPWIND1 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DX_M1(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DX_M1(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DX_M1(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DX_M1(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DX_M1(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DX_M1(gyz,i,j,k)

      else if (admmacros_advectionx .eq. 2) then /* UPWIND2 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DX_P2(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DX_P2(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DX_P2(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DX_P2(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DX_P2(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DX_P2(gyz,i,j,k)

      else if (admmacros_advectionx .eq. -2) then /* UPWIND2 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DX_M2(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DX_M2(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DX_M2(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DX_M2(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DX_M2(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DX_M2(gyz,i,j,k)

      end if

        
      /* Do the y-direction */

      if (admmacros_advectiony .eq. 0) then /* CENTER */

         LIEG_LGXX = LIEG_LGXX + DYDCG_DYDCGXX*LIEG_BY
         LIEG_LGYY = LIEG_LGYY + DYDCG_DYDCGYY*LIEG_BY 
         LIEG_LGZZ = LIEG_LGZZ + DYDCG_DYDCGZZ*LIEG_BY 
         LIEG_LGXY = LIEG_LGXY + DYDCG_DYDCGXY*LIEG_BY
         LIEG_LGXZ = LIEG_LGXZ + DYDCG_DYDCGXZ*LIEG_BY 
         LIEG_LGYZ = LIEG_LGYZ + DYDCG_DYDCGYZ*LIEG_BY 

      else if (admmacros_advectiony .eq. 1) then /* UPWIND1 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DY_P1(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DY_P1(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DY_P1(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DY_P1(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DY_P1(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DY_P1(gyz,i,j,k)

      else if (admmacros_advectiony .eq. -1) then /* UPWIND1 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DY_M1(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DY_M1(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DY_M1(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DY_M1(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DY_M1(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DY_M1(gyz,i,j,k)

      else if (admmacros_advectiony .eq. 2) then /* UPWIND2 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DY_P2(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DY_P2(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DY_P2(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DY_P2(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DY_P2(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DY_P2(gyz,i,j,k)

      else if (admmacros_advectiony .eq. -2) then /* UPWIND2 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DY_M2(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DY_M2(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DY_M2(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DY_M2(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DY_M2(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DY_M2(gyz,i,j,k)

      end if


      /* Do the z-direction */

      if (admmacros_advectionz .eq. 0) then /* CENTER */

         LIEG_LGXX = LIEG_LGXX + DZDCG_DZDCGXX*LIEG_BZ 
         LIEG_LGYY = LIEG_LGYY + DZDCG_DZDCGYY*LIEG_BZ
         LIEG_LGZZ = LIEG_LGZZ + DZDCG_DZDCGZZ*LIEG_BZ 
         LIEG_LGXY = LIEG_LGXY + DZDCG_DZDCGXY*LIEG_BZ
         LIEG_LGXZ = LIEG_LGXZ + DZDCG_DZDCGXZ*LIEG_BZ 
         LIEG_LGYZ = LIEG_LGYZ + DZDCG_DZDCGYZ*LIEG_BZ 

      else if ((admmacros_advectionz .eq. 1)) then /* UPWIND1 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DZ_P1(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DZ_P1(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DZ_P1(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DZ_P1(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DZ_P1(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DZ_P1(gyz,i,j,k)

      else if ((admmacros_advectionz .eq. -1)) then /* UPWIND1 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DZ_M1(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DZ_M1(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DZ_M1(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DZ_M1(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DZ_M1(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DZ_M1(gyz,i,j,k)

      else if (admmacros_advectionz .eq. 2) then /* UPWIND2 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DZ_P2(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DZ_P2(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DZ_P2(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DZ_P2(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DZ_P2(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DZ_P2(gyz,i,j,k)

      else if (admmacros_advectionz .eq. -2) then /* UPWIND2 */

         LIEG_LGXX = LIEG_LGXX + ADM_ADV_DZ_M2(gxx,i,j,k)
         LIEG_LGYY = LIEG_LGYY + ADM_ADV_DZ_M2(gyy,i,j,k)
         LIEG_LGZZ = LIEG_LGZZ + ADM_ADV_DZ_M2(gzz,i,j,k)
         LIEG_LGXY = LIEG_LGXY + ADM_ADV_DZ_M2(gxy,i,j,k)
         LIEG_LGXZ = LIEG_LGXZ + ADM_ADV_DZ_M2(gxz,i,j,k)
         LIEG_LGYZ = LIEG_LGYZ + ADM_ADV_DZ_M2(gyz,i,j,k)

      end if 

/*  Extra terms in the Lie derivative.  */

      LIEG_LGXX = LIEG_LGXX + 2.0D0*(DXDB_DXDBX*LIEG_GXX \
        + DXDB_DXDBY*LIEG_GXY + DXDB_DXDBZ*LIEG_GXZ)

      LIEG_LGYY = LIEG_LGYY + 2.0D0*(DYDB_DYDBX*LIEG_GXY \
        + DYDB_DYDBY*LIEG_GYY + DYDB_DYDBZ*LIEG_GYZ)

      LIEG_LGZZ = LIEG_LGZZ + 2.0D0*(DZDB_DZDBX*LIEG_GXZ \
        + DZDB_DZDBY*LIEG_GYZ + DZDB_DZDBZ*LIEG_GZZ)

      LIEG_LGXY = LIEG_LGXY + DYDB_DYDBX*LIEG_GXX + DXDB_DXDBY*LIEG_GYY \
        + (DXDB_DXDBX + DYDB_DYDBY)*LIEG_GXY \
        + DYDB_DYDBZ*LIEG_GXZ + DXDB_DXDBZ*LIEG_GYZ

      LIEG_LGXZ = LIEG_LGXZ + DZDB_DZDBX*LIEG_GXX + DXDB_DXDBZ*LIEG_GZZ \
        + (DXDB_DXDBX + DZDB_DZDBZ)*LIEG_GXZ \
        + DZDB_DZDBY*LIEG_GXY + DXDB_DXDBY*LIEG_GYZ

      LIEG_LGYZ = LIEG_LGYZ + DZDB_DZDBY*LIEG_GYY + DYDB_DYDBZ*LIEG_GZZ \
        + (DYDB_DYDBY + DZDB_DZDBZ)*LIEG_GYZ \
        + DZDB_DZDBX*LIEG_GXY + DYDB_DYDBX*LIEG_GXZ

#endif

#ifdef CCODE

/* Advection */

      LIEG_LGXX = DXDCG_DXDCGXX*LIEG_BX + DYDCG_DYDCGXX*LIEG_BY
         + DZDCG_DZDCGXX*LIEG_BZ;

      LIEG_LGYY = DXDCG_DXDCGYY*LIEG_BX + DYDCG_DYDCGYY*LIEG_BY
         + DZDCG_DZDCGYY*LIEG_BZ;

      LIEG_LGZZ = DXDCG_DXDCGZZ*LIEG_BX + DYDCG_DYDCGZZ*LIEG_BY
         + DZDCG_DZDCGZZ*LIEG_BZ;

      LIEG_LGXY = DXDCG_DXDCGXY*LIEG_BX + DYDCG_DYDCGXY*LIEG_BY
         + DZDCG_DZDCGXY*LIEG_BZ;

      LIEG_LGXZ = DXDCG_DXDCGXZ*LIEG_BX + DYDCG_DYDCGXZ*LIEG_BY
         + DZDCG_DZDCGXZ*LIEG_BZ;

      LIEG_LGYZ = DXDCG_DXDCGYZ*LIEG_BX + DYDCG_DYDCGYZ*LIEG_BY
         + DZDCG_DZDCGYZ*LIEG_BZ;

/* Extra terms in the Lie derivative */

      LIEG_LGXX = LIEG_LGXX + 2*(DXDB_DXDBX*LIEG_GXX
         + DXDB_DXDBY*LIEG_GXY + DXDB_DXDBZ*LIEG_GXZ);

      LIEG_LGYY = LIEG_LGYY + 2*(DYDB_DYDBX*LIEG_GXY
         + DYDB_DYDBY*LIEG_GYY + DYDB_DYDBZ*LIEG_GYZ);

      LIEG_LGZZ = LIEG_LGZZ + 2*(DZDB_DZDBX*LIEG_GXZ
         + DZDB_DZDBY*LIEG_GYZ + DZDB_DZDBZ*LIEG_GZZ);

      LIEG_LGXY = LIEG_LGXY + DYDB_DYDBX*LIEG_GXX + DXDB_DXDBY*LIEG_GYY
         + (DXDB_DXDBX + DYDB_DYDBY)*LIEG_GXY
         + DYDB_DYDBZ*LIEG_GXZ + DXDB_DXDBZ*LIEG_GYZ;

      LIEG_LGXZ = LIEG_LGXZ + DZDB_DZDBX*LIEG_GXX + DXDB_DXDBZ*LIEG_GZZ
         + (DXDB_DXDBX + DZDB_DZDBZ)*LIEG_GXZ
         + DZDB_DZDBY*LIEG_GXY + DXDB_DXDBY*LIEG_GYZ;

      LIEG_LGYZ = LIEG_LGYZ + DZDB_DZDBY*LIEG_GYY + DYDB_DYDBZ*LIEG_GZZ
         + (DYDB_DYDBY + DZDB_DZDBZ)*LIEG_GYZ
         + DZDB_DZDBX*LIEG_GXY + DYDB_DYDBX*LIEG_GXZ;

#endif

#endif
