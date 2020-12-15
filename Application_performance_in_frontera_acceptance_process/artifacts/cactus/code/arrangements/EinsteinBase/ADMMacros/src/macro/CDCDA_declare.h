/*@@
  @header   CDCDA_declare.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro to calculate all second covariant spatial derivative of lapse

  That is alpha_{;ij}
  @enddesc
@@*/

#ifndef CDCDA_DECLARE
#define CDCDA_DECLARE

#include "DA_declare.h"
#include "DDA_declare.h"
#include "CHR2_declare.h"

/* Output variables */ 
#undef  CDCDA_CDXXDA
#define CDCDA_CDXXDA cdcda_cdxxda
#undef  CDCDA_CDXYDA
#define CDCDA_CDXYDA cdcda_cdxyda
#undef  CDCDA_CDXZDA
#define CDCDA_CDXZDA cdcda_cdxzda
#undef  CDCDA_CDYYDA
#define CDCDA_CDYYDA cdcda_cdyyda
#undef  CDCDA_CDYZDA
#define CDCDA_CDYZDA cdcda_cdyzda
#undef  CDCDA_CDZZDA
#define CDCDA_CDZZDA cdcda_cdzzda

#ifdef FCODE

/* Declare output variables */
      CCTK_REAL CDCDA_CDXXDA
      CCTK_REAL CDCDA_CDXYDA
      CCTK_REAL CDCDA_CDXZDA
      CCTK_REAL CDCDA_CDYYDA
      CCTK_REAL CDCDA_CDYZDA
      CCTK_REAL CDCDA_CDZZDA

#endif

#ifdef CCODE

/* Declare output variables */
      CCTK_REAL CDCDA_CDXXDA;
      CCTK_REAL CDCDA_CDXYDA;
      CCTK_REAL CDCDA_CDXZDA;
      CCTK_REAL CDCDA_CDYYDA;
      CCTK_REAL CDCDA_CDYZDA;
      CCTK_REAL CDCDA_CDZZDA;

#endif

#endif
