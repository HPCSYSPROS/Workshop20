/*@@
  @header   HAMADM_guts.h
  @date     Aug 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the spacetime part of the 
  Hamiltonian Constraint. That is:

       R - K^i_j K^j_i + trK^2 

  @enddesc
@@*/

#ifndef HAMADM_GUTS
#define HAMADM_GUTS

#include "TRRICCI_guts.h"
#include "TRKK_guts.h"
#include "TRK_guts.h"

#ifdef FCODE 

      HAMADM_HAMADM = TRRICCI_TRRICCI - TRKK_TRKK + TRK_TRK**2
      HAMADM_HAMADMABS = abs(TRRICCI_TRRICCI) + abs(TRKK_TRKK) + TRK_TRK**2

#endif

#ifdef CCODE

      HAMADM_HAMADM = TRRICCI_TRRICCI - TRKK_TRKK + TRK_TRK*TRK_TRK;
      HAMADM_HAMADMABS = fabs(TRRICCI_TRRICCI) + fabs(TRKK_TRKK) +
                         TRK_TRK*TRK_TRK;

#endif

#endif

