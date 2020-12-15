/*@@
  @header   DZDB_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc

  Macro to calculate the first derivatives of the 
  shift with respect to z

  The macro is defined in terms of standard variables in
  @seefile DZDB_declare.h
  @enddesc
@@*/

#ifndef DZDB_GUTS
#define DZDB_GUTS

#ifdef FCODE 

#include "ADM_Derivative.h"
      if (local_spatial_order.eq.2) then
        DZDB_DZDBX = ADM_DZ_2(betax,i,j,k)
        DZDB_DZDBY = ADM_DZ_2(betay,i,j,k)
        DZDB_DZDBZ = ADM_DZ_2(betaz,i,j,k)
      else
        DZDB_DZDBX = ADM_DZ_4(betax,i,j,k)
        DZDB_DZDBY = ADM_DZ_4(betay,i,j,k)
        DZDB_DZDBZ = ADM_DZ_4(betaz,i,j,k)
      end if
#endif

#ifdef CCODE

      DZDB_OO2DZ = 1/(2*cctkGH->cctk_delta_space[2]);
    
      DZDB_DZDBX = DZDB_OO2DZ*(DZDB_BX_KP - DZDB_BX_KM);
      DZDB_DZDBY = DZDB_OO2DZ*(DZDB_BY_KP - DZDB_BY_KM);
      DZDB_DZDBZ = DZDB_OO2DZ*(DZDB_BZ_KP - DZDB_BZ_KM);

#endif

#endif

