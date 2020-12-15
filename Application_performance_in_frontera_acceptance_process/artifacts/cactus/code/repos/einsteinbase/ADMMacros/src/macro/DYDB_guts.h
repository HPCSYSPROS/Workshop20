/*@@
  @header   DYDB_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc

  Macro to calculate the first derivatives of the 
  shift with respect to y

  The macro is defined in terms of standard variables in
  @seefile DYDB_declare.h
  @enddesc
@@*/

#ifndef DYDB_GUTS
#define DYDB_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"

      if (local_spatial_order.eq.2) then
        DYDB_DYDBX = ADM_DY_2(betax,i,j,k)
        DYDB_DYDBY = ADM_DY_2(betay,i,j,k)
        DYDB_DYDBZ = ADM_DY_2(betaz,i,j,k)
      else
        DYDB_DYDBX = ADM_DY_4(betax,i,j,k)
        DYDB_DYDBY = ADM_DY_4(betay,i,j,k)
        DYDB_DYDBZ = ADM_DY_4(betaz,i,j,k)
      end if

#endif

#ifdef CCODE

      DYDB_OO2DY = 1/(2*cctkGH->cctk_delta_space[1]);
    
      DYDB_DYDBX = DYDB_OO2DY*(DYDB_BX_JP - DYDB_BX_JM);
      DYDB_DYDBY = DYDB_OO2DY*(DYDB_BY_JP - DYDB_BY_JM);
      DYDB_DYDBZ = DYDB_OO2DY*(DYDB_BZ_JP - DYDB_BZ_JM);

#endif

#endif

