/*@@
  @header   HAMADM_declare.h
  @date     Aug 98
  @author   Gabrielle Allen
  @desc
  Declarations for macro to calculate vacuum part of Hamiltonian
  constraint
 
  @enddesc
@@*/

#ifndef HAMADM_DECLARE
#define HAMADM_DECLARE

#include "TRRICCI_declare.h"
#include "TRKK_declare.h"
#include "TRK_declare.h"

#ifdef FCODE

/* Output variables */ 
#undef  HAMADM_HAMADM
#define HAMADM_HAMADM hamadm_ham_adm

#undef  HAMADM_HAMADMABS
#define HAMADM_HAMADMABS hamadm_ham_adm_abs

/* Declare output variables */
      CCTK_REAL HAMADM_HAMADM
      CCTK_REAL HAMADM_HAMADMABS

#endif


#ifdef CCODE

/* Output variables */
#undef  HAMADM_HAMADM
#define HAMADM_HAMADM hamadm_ham_adm

#undef  HAMADM_HAMADMABS
#define HAMADM_HAMADMABS hamadm_ham_adm_abs

/* Declare output variables */
      CCTK_REAL HAMADM_HAMADM;
      CCTK_REAL HAMADM_HAMADMABS;

#endif

#endif
