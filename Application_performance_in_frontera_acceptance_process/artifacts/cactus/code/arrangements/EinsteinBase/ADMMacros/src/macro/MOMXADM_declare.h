/*@@
  @header   MOMXADM_declare.h
  @date     Aug 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro to calculate vacuum part of 
  x-Momentum constraint
 
  @enddesc
@@*/

#ifndef MOMXADM_DECLARE
#define MOMXADM_DECLARE

#include "UPPERMET_declare.h"
#include "CDK_declare.h"

#ifdef FCODE

/* Output variables */ 
#undef  MOMXADM_MOMXADM
#define MOMXADM_MOMXADM momxadm_momxadm

/* Declare output variables */
      CCTK_REAL MOMXADM_MOMXADM

#endif


#ifdef CCODE

/* Output variables */
#undef  MOMXADM_MOMXADM
#define MOMXADM_MOMXADM momxadm_momxadm

/* Declare output variables */
      CCTK_REAL MOMXADM_MOMXADM;

#endif

#endif
