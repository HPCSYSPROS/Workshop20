/*@@
  @header   DZDG_guts.h
  @date     Jun 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the first derivatives of the 
  physical metric with respect to z

  The macro is defined in terms of standard variables in
  @seefile DZDG_declare.h
  @enddesc
@@*/

#ifndef DZDG_GUTS
#define DZDG_GUTS

#include "DZDCG_guts.h"

#ifdef FCODE 

      DZDG_PSI4 = 1
      DZDG_FAC  = 0

      DZDG_DZDGXX = DZDCG_DZDCGXX*DZDG_PSI4 + DZDG_FAC*DZDG_GXX
      DZDG_DZDGXY = DZDCG_DZDCGXY*DZDG_PSI4 + DZDG_FAC*DZDG_GXY
      DZDG_DZDGXZ = DZDCG_DZDCGXZ*DZDG_PSI4 + DZDG_FAC*DZDG_GXZ
      DZDG_DZDGYY = DZDCG_DZDCGYY*DZDG_PSI4 + DZDG_FAC*DZDG_GYY
      DZDG_DZDGYZ = DZDCG_DZDCGYZ*DZDG_PSI4 + DZDG_FAC*DZDG_GYZ
      DZDG_DZDGZZ = DZDCG_DZDCGZZ*DZDG_PSI4 + DZDG_FAC*DZDG_GZZ

#endif

#ifdef CCODE

      DZDG_PSI4 = 1;

      DZDG_FAC  = 0;

      DZDG_DZDGXX = DZDCG_DZDCGXX*DZDG_PSI4 + DZDG_FAC*DZDG_GXX;
      DZDG_DZDGXY = DZDCG_DZDCGXY*DZDG_PSI4 + DZDG_FAC*DZDG_GXY;
      DZDG_DZDGXZ = DZDCG_DZDCGXZ*DZDG_PSI4 + DZDG_FAC*DZDG_GXZ;
      DZDG_DZDGYY = DZDCG_DZDCGYY*DZDG_PSI4 + DZDG_FAC*DZDG_GYY;
      DZDG_DZDGYZ = DZDCG_DZDCGYZ*DZDG_PSI4 + DZDG_FAC*DZDG_GYZ;
      DZDG_DZDGZZ = DZDCG_DZDCGZZ*DZDG_PSI4 + DZDG_FAC*DZDG_GZZ;

#endif

#endif
