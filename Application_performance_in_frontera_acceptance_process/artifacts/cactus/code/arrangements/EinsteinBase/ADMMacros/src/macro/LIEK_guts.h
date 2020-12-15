/*@@
  @header   LIEK_guts.h
  @date     Jul 2000
  @author   Gabrielle Allen + Miguel Alcubierre
  @desc

  Macro to calculate the Lie derivative of the lower 
  extrinsic curvature along the shift vector.   Notice
  that the advection term can be calculated in several
  different ways.

  IMPORTANT:  THE C VERSION ONLY HAS CENTERED DIFFERENCES!

  @enddesc
@@*/

#ifndef LIEK_GUTS
#define LIEK_GUTS

#include "DB_guts.h"
#include "DK_guts.h"

#ifdef FCODE

#include "ADM_Derivative.h"

/*     Advection. */

      LIEK_LKXX = 0.0D0
      LIEK_LKYY = 0.0D0
      LIEK_LKZZ = 0.0D0
      LIEK_LKXY = 0.0D0
      LIEK_LKXZ = 0.0D0
      LIEK_LKYZ = 0.0D0

      if (admmacros_advectionx .eq. 0) then /* CENTER X */
         
         LIEK_LKXX = LIEK_LKXX + DXDK_DXDKXX*LIEK_BX
         LIEK_LKYY = LIEK_LKYY + DXDK_DXDKYY*LIEK_BX
         LIEK_LKZZ = LIEK_LKZZ + DXDK_DXDKZZ*LIEK_BX
         LIEK_LKXY = LIEK_LKXY + DXDK_DXDKXY*LIEK_BX
         LIEK_LKXZ = LIEK_LKXZ + DXDK_DXDKXZ*LIEK_BX
         LIEK_LKYZ = LIEK_LKYZ + DXDK_DXDKYZ*LIEK_BX

      else if (admmacros_advectionx .eq. 1) then /* UPWIND1 X */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DX_P1(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DX_P1(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DX_P1(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DX_P1(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DX_P1(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DX_P1(kyz,i,j,k)

      else if (admmacros_advectionx .eq. -1) then /* UPWIND1 X */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DX_M1(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DX_M1(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DX_M1(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DX_M1(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DX_M1(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DX_M1(kyz,i,j,k)

      else if (admmacros_advectionx .eq. 2) then /* UPWIND2 X */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DX_P2(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DX_P2(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DX_P2(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DX_P2(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DX_P2(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DX_P2(kyz,i,j,k)

      else if (admmacros_advectionx .eq. -2) then /* UPWIND2 X */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DX_M2(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DX_M2(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DX_M2(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DX_M2(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DX_M2(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DX_M2(kyz,i,j,k)

      end if


      if (admmacros_advectiony .eq. 0) then /* CENTER Y */

         LIEK_LKXX = LIEK_LKXX + DYDK_DYDKXX*LIEK_BY
         LIEK_LKYY = LIEK_LKYY + DYDK_DYDKYY*LIEK_BY
         LIEK_LKZZ = LIEK_LKZZ + DYDK_DYDKZZ*LIEK_BY
         LIEK_LKXY = LIEK_LKXY + DYDK_DYDKXY*LIEK_BY
         LIEK_LKXZ = LIEK_LKXZ + DYDK_DYDKXZ*LIEK_BY
         LIEK_LKYZ = LIEK_LKYZ + DYDK_DYDKYZ*LIEK_BY

      else if (admmacros_advectiony .eq. 1) then /* UPWIND1 Y */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_P1(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_P1(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_P1(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_P1(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_P1(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_P1(kyz,i,j,k)

      else if (admmacros_advectiony .eq. -1) then /* UPWIND1 Y */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_M1(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_M1(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_M1(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_M1(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_M1(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_M1(kyz,i,j,k)

      else if (admmacros_advectiony .eq. 2) then /* UPWIND2 Y */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_P2(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_P2(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_P2(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_P2(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_P2(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_P2(kyz,i,j,k)

      else if (admmacros_advectiony .eq. -2) then /* UPWIND2 Y */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_M2(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_M2(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_M2(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_M2(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_M2(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_M2(kyz,i,j,k)

      end if


      if (admmacros_advectionz .eq. 0) then /* CENTER Z */

         LIEK_LKXX = LIEK_LKXX + DZDK_DZDKXX*LIEK_BZ
         LIEK_LKYY = LIEK_LKYY + DZDK_DZDKYY*LIEK_BZ
         LIEK_LKZZ = LIEK_LKZZ + DZDK_DZDKZZ*LIEK_BZ
         LIEK_LKXY = LIEK_LKXY + DZDK_DZDKXY*LIEK_BZ
         LIEK_LKXZ = LIEK_LKXZ + DZDK_DZDKXZ*LIEK_BZ
         LIEK_LKYZ = LIEK_LKYZ + DZDK_DZDKYZ*LIEK_BZ

      else if (admmacros_advectionz .eq. 1) then /* UPWIND1 Z */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_P1(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_P1(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_P1(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_P1(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_P1(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_P1(kyz,i,j,k)

      else if (admmacros_advectionz .eq. -1) then /* UPWIND1 Z */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_M1(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_M1(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_M1(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_M1(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_M1(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_M1(kyz,i,j,k)

      else if (admmacros_advectionz .eq. 2) then /* UPWIND2 Z */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_P2(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_P2(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_P2(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_P2(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_P2(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_P2(kyz,i,j,k)

      else if (admmacros_advectionz .eq. -2) then /* UPWIND2 Z */

         LIEK_LKXX = LIEK_LKXX + ADM_ADV_DY_M2(kxx,i,j,k)
         LIEK_LKYY = LIEK_LKYY + ADM_ADV_DY_M2(kyy,i,j,k)
         LIEK_LKZZ = LIEK_LKZZ + ADM_ADV_DY_M2(kzz,i,j,k)
         LIEK_LKXY = LIEK_LKXY + ADM_ADV_DY_M2(kxy,i,j,k)
         LIEK_LKXZ = LIEK_LKXZ + ADM_ADV_DY_M2(kxz,i,j,k)
         LIEK_LKYZ = LIEK_LKYZ + ADM_ADV_DY_M2(kyz,i,j,k)

      end if

      /*    Extra terms in the Lie derivative. */

      LIEK_LKXX = LIEK_LKXX + 2.0D0*(DXDB_DXDBX*LIEK_KXX \
      + DXDB_DXDBY*LIEK_KXY + DXDB_DXDBZ*LIEK_KXZ)

      LIEK_LKYY = LIEK_LKYY + 2.0D0*(DYDB_DYDBX*LIEK_KXY \
      + DYDB_DYDBY*LIEK_KYY + DYDB_DYDBZ*LIEK_KYZ)

      LIEK_LKZZ = LIEK_LKZZ + 2.0D0*(DZDB_DZDBX*LIEK_KXZ \
      + DZDB_DZDBY*LIEK_KYZ + DZDB_DZDBZ*LIEK_KZZ)

      LIEK_LKXY = LIEK_LKXY + DYDB_DYDBX*LIEK_KXX + DXDB_DXDBY*LIEK_KYY \
      + (DXDB_DXDBX + DYDB_DYDBY)*LIEK_KXY \
      + DYDB_DYDBZ*LIEK_KXZ + DXDB_DXDBZ*LIEK_KYZ

      LIEK_LKXZ = LIEK_LKXZ + DZDB_DZDBX*LIEK_KXX + DXDB_DXDBZ*LIEK_KZZ \
      + (DXDB_DXDBX + DZDB_DZDBZ)*LIEK_KXZ \
      + DZDB_DZDBY*LIEK_KXY + DXDB_DXDBY*LIEK_KYZ

      LIEK_LKYZ = LIEK_LKYZ + DZDB_DZDBY*LIEK_KYY + DYDB_DYDBZ*LIEK_KZZ \
      + (DYDB_DYDBY + DZDB_DZDBZ)*LIEK_KYZ \
      + DZDB_DZDBX*LIEK_KXY + DYDB_DYDBX*LIEK_KXZ


#endif

#ifdef CCODE

/* Advection */

      LIEK_LKXX = DXDK_DXDKXX*LIEK_BX + DYDK_DYDKXX*LIEK_BY
         + DZDK_DZDKXX*LIEK_BZ;

      LIEK_LKYY = DXDK_DXDKYY*LIEK_BX + DYDK_DYDKYY*LIEK_BY
         + DZDK_DZDKYY*LIEK_BZ;

      LIEK_LKZZ = DXDK_DXDKZZ*LIEK_BX + DYDK_DYDKZZ*LIEK_BY
         + DZDK_DZDKZZ*LIEK_BZ;

      LIEK_LKXY = DXDK_DXDKXY*LIEK_BX + DYDK_DYDKXY*LIEK_BY
         + DZDK_DZDKXY*LIEK_BZ;

      LIEK_LKXZ = DXDK_DXDKXZ*LIEK_BX + DYDK_DYDKXZ*LIEK_BY
         + DZDK_DZDKXZ*LIEK_BZ;

      LIEK_LKYZ = DXDK_DXDKYZ*LIEK_BX + DYDK_DYDKYZ*LIEK_BY
         + DZDK_DZDKYZ*LIEK_BZ;

/* Extra terms in the Lie derivative */

      LIEK_LKXX = LIEK_LKXX + 2*(DXDB_DXDBX*LIEK_KXX
         + DXDB_DXDBY*LIEK_KXY + DXDB_DXDBZ*LIEK_KXZ);

      LIEK_LKYY = LIEK_LKYY + 2*(DYDB_DYDBX*LIEK_KXY
         + DYDB_DYDBY*LIEK_KYY + DYDB_DYDBZ*LIEK_KYZ);

      LIEK_LKZZ = LIEK_LKZZ + 2*(DZDB_DZDBX*LIEK_KXZ
         + DZDB_DZDBY*LIEK_KYZ + DZDB_DZDBZ*LIEK_KZZ);

      LIEK_LKXY = LIEK_LKXY + DYDB_DYDBX*LIEK_KXX + DXDB_DXDBY*LIEK_KYY
         + (DXDB_DXDBX + DYDB_DYDBY)*LIEK_KXY
         + DYDB_DYDBZ*LIEK_KXZ + DXDB_DXDBZ*LIEK_KYZ;

      LIEK_LKXZ = LIEK_LKXZ + DZDB_DZDBX*LIEK_KXX + DXDB_DXDBZ*LIEK_KZZ
         + (DXDB_DXDBX + DZDB_DZDBZ)*LIEK_KXZ
         + DZDB_DZDBY*LIEK_KXY + DXDB_DXDBY*LIEK_KYZ;

      LIEK_LKYZ = LIEK_LKYZ + DZDB_DZDBY*LIEK_KYY + DYDB_DYDBZ*LIEK_KZZ
         + (DYDB_DYDBY + DZDB_DZDBZ)*LIEK_KYZ
         + DZDB_DZDBX*LIEK_KXY + DYDB_DYDBX*LIEK_KXZ;

#endif

#endif
