/*@@
  @header   DYDB_declare.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro @seefile DYDB_guts.h to compute first 
  derivatives of the shift with respect to y
  @enddesc
@@*/

#ifndef DYDB_DECLARE
#define DYDB_DECLARE

#ifdef FCODE

/* Input variables */
#undef  DYDB_BX
#define DYDB_BX betax(i,j,k)
#undef  DYDB_BY
#define DYDB_BY betay(i,j,k)
#undef  DYDB_BZ
#define DYDB_BZ betaz(i,j,k)
#undef  DYDB_BX_JP 
#define DYDB_BX_JP betax(i,j+1,k)
#undef  DYDB_BX_JPP 
#define DYDB_BX_JPP betax(i,j+2,k)
#undef  DYDB_BX_JM
#define DYDB_BX_JM betax(i,j-1,k)
#undef  DYDB_BX_JMM
#define DYDB_BX_JMM betax(i,j-2,k)
#undef  DYDB_BY_JP 
#define DYDB_BY_JP betay(i,j+1,k)
#undef  DYDB_BY_JPP 
#define DYDB_BY_JPP betay(i,j+2,k)
#undef  DYDB_BY_JM 
#define DYDB_BY_JM betay(i,j-1,k)
#undef  DYDB_BY_JMM 
#define DYDB_BY_JMM betay(i,j-2,k)
#undef  DYDB_BZ_JP 
#define DYDB_BZ_JP betaz(i,j+1,k)
#undef  DYDB_BZ_JPP 
#define DYDB_BZ_JPP betaz(i,j+2,k)
#undef  DYDB_BZ_JM
#define DYDB_BZ_JM betaz(i,j-1,k)
#undef  DYDB_BZ_JMM
#define DYDB_BZ_JMM betaz(i,j-2,k)

/* Output variables */ 
#undef  DYDB_DYDBX
#define DYDB_DYDBX  dydb_dydbx
#undef  DYDB_DYDBY
#define DYDB_DYDBY  dydb_dydby
#undef  DYDB_DYDBZ
#define DYDB_DYDBZ  dydb_dydbz

/* Internal variables */
#undef  DYDB_DY
#define DYDB_DY dy
#undef  DYDB_OO2DY   
#define DYDB_OO2DY dydb_oo2dy

/* Declare internal variables */
      CCTK_REAL DYDB_OO2DY

/* Declare output variables */
      CCTK_REAL DYDB_DYDBX
      CCTK_REAL DYDB_DYDBY
      CCTK_REAL DYDB_DYDBZ

#endif

#ifdef CCODE

/* Input variables */
#undef  DYDB_BX_JP 
#define DYDB_BX_JP betax[ dj+ijk]
#undef  DYDB_BX_JM
#define DYDB_BX_JM betax[-dj+ijk]
#undef  DYDB_BY_JP 
#define DYDB_BY_JP betay[ dj+ijk]
#undef  DYDB_BY_JM 
#define DYDB_BY_JM betay[-dj+ijk]
#undef  DYDB_BZ_JP 
#define DYDB_BZ_JP betaz[ dj+ijk]
#undef  DYDB_BZ_JM
#define DYDB_BZ_JM betaz[-dj+ijk]

/* Output variables */ 
#undef  DYDB_DYDBX
#define DYDB_DYDBX  delbeta21
#undef  DYDB_DYDBY
#define DYDB_DYDBY  delbeta22
#undef  DYDB_DYDBZ
#define DYDB_DYDBZ  delbeta23

/* Internal variables */
#undef  DYDB_OO2DY   
#define DYDB_OO2DY dydb_oo2dy

/* Declare internal variables */
CCTK_REAL DYDB_OO2DY;

/* Declare output variables */
CCTK_REAL DYDB_DYDBX;
CCTK_REAL DYDB_DYDBY;
CCTK_REAL DYDB_DYDBZ;

#endif

#endif
