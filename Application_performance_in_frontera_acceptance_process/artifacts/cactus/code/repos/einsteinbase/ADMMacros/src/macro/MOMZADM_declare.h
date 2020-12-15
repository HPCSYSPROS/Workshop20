/*@@
  @header   MOMZADM_declare.h
  @date     Aug 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro to calculate vacuum part of 
  z-Momentum constraint
 
  @enddesc
@@*/

#ifndef MOMZADM_DECLARE
#define MOMZADM_DECLARE

#include "UPPERMET_declare.h"
#include "CDK_declare.h"

#ifdef FCODE

/* Output variables */ 
#undef  MOMZADM_MOMZADM
#define MOMZADM_MOMZADM momzadm_momzadm

/* Declare output variables */
      CCTK_REAL MOMZADM_MOMZADM

#endif


#ifdef CCODE

/* Output variables */
#undef  MOMZADM_MOMZADM
#define MOMZADM_MOMZADM momzadm_momzadm

/* Declare output variables */
      CCTK_REAL MOMZADM_MOMZADM;

#endif

#endif
