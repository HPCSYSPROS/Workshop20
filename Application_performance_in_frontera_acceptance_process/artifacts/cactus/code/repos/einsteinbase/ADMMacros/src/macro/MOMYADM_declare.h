/*@@
  @header   MOMYADM_declare.h
  @date     Aug 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro to calculate vacuum part of 
  y-Momentum constraint
 
  @enddesc
@@*/

#ifndef MOMYADM_DECLARE
#define MOMYADM_DECLARE

#include "UPPERMET_declare.h"
#include "CDK_declare.h"

#ifdef FCODE

/* Output variables */ 
#undef  MOMYADM_MOMYADM
#define MOMYADM_MOMYADM momyadm_momyadm

/* Declare output variables */
      CCTK_REAL MOMYADM_MOMYADM

#endif


#ifdef CCODE

/* Output variables */
#undef  MOMYADM_MOMYADM
#define MOMYADM_MOMYADM momyadm_momyadm

/* Declare output variables */
      CCTK_REAL MOMYADM_MOMYADM;

#endif

#endif
