/*@@
  @header   CDXCDK_declare.h
  @date     Aug 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro to compute first conformal
  derivatives of the extrinsic curvature with respect to x
  @enddesc
@@*/

#ifndef CDXCDK_DECLARE
#define CDXCDK_DECLARE

#include "DXDK_declare.h"
#include "CHR2_declare.h"

#ifdef FCODE

/* Input variables */
#undef  CDXCDK_KXX 
#define CDXCDK_KXX kxx(i,j,k)
#undef  CDXCDK_KXY 
#define CDXCDK_KXY kxy(i,j,k)
#undef  CDXCDK_KXZ 
#define CDXCDK_KXZ kxz(i,j,k)
#undef  CDXCDK_KYY 
#define CDXCDK_KYY kyy(i,j,k)
#undef  CDXCDK_KYZ 
#define CDXCDK_KYZ kyz(i,j,k)
#undef  CDXCDK_KZZ 
#define CDXCDK_KZZ kzz(i,j,k)

/* Output variables */ 
#undef  CDXCDK_CDXCDKXX
#define CDXCDK_CDXCDKXX  cdxcdk_cdxcdkxx
#undef  CDXCDK_CDXCDKXY
#define CDXCDK_CDXCDKXY  cdxcdk_cdxcdkxy
#undef  CDXCDK_CDXCDKXZ
#define CDXCDK_CDXCDKXZ  cdxcdk_cdxcdkxz
#undef  CDXCDK_CDXCDKYY
#define CDXCDK_CDXCDKYY  cdxcdk_cdxcdkyy
#undef  CDXCDK_CDXCDKYZ
#define CDXCDK_CDXCDKYZ  cdxcdk_cdxcdkyz
#undef  CDXCDK_CDXCDKZZ
#define CDXCDK_CDXCDKZZ  cdxcdk_cdxcdkzz

/* Declare output variables */
      CCTK_REAL CDXCDK_CDXCDKXX
      CCTK_REAL CDXCDK_CDXCDKXY
      CCTK_REAL CDXCDK_CDXCDKXZ
      CCTK_REAL CDXCDK_CDXCDKYY
      CCTK_REAL CDXCDK_CDXCDKYZ
      CCTK_REAL CDXCDK_CDXCDKZZ

#endif


#ifdef CCODE

/* Input variables */
#undef  CDXCDK_KXX 
#define CDXCDK_KXX kxx[ijk]
#undef  CDXCDK_KXY 
#define CDXCDK_KXY kxy[ijk]
#undef  CDXCDK_KXZ 
#define CDXCDK_KXZ kxz[ijk]
#undef  CDXCDK_KYY 
#define CDXCDK_KYY kyy[ijk]
#undef  CDXCDK_KYZ 
#define CDXCDK_KYZ kyz[ijk]
#undef  CDXCDK_KZZ 
#define CDXCDK_KZZ kzz[ijk]

/* Output variables */ 
#undef  CDXCDK_CDXCDKXX
#define CDXCDK_CDXCDKXX  cdxcdk_cdxcdkxx
#undef  CDXCDK_CDXCDKXY
#define CDXCDK_CDXCDKXY  cdxcdk_cdxcdkxy
#undef  CDXCDK_CDXCDKXZ
#define CDXCDK_CDXCDKXZ  cdxcdk_cdxcdkxz
#undef  CDXCDK_CDXCDKYY
#define CDXCDK_CDXCDKYY  cdxcdk_cdxcdkyy
#undef  CDXCDK_CDXCDKYZ
#define CDXCDK_CDXCDKYZ  cdxcdk_cdxcdkyz
#undef  CDXCDK_CDXCDKZZ
#define CDXCDK_CDXCDKZZ  cdxcdk_cdxcdkzz

/* Declare output variables */
CCTK_REAL CDXCDK_CDXCDKXX;
CCTK_REAL CDXCDK_CDXCDKXY;
CCTK_REAL CDXCDK_CDXCDKXZ;
CCTK_REAL CDXCDK_CDXCDKYY;
CCTK_REAL CDXCDK_CDXCDKYZ;
CCTK_REAL CDXCDK_CDXCDKZZ;

#endif

#endif
