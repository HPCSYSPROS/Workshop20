/*@@
  @header   DCGDT_guts.h
  @date     Jul 98
  @author   Gabrielle Allen
  @desc
  Macro to calculate the source term in the evolution equation for the
  conformal 3-metric. That is

  d g~_ij/dt =( - 2 alpha K_ij + L_beta g_ij )/Psi^4

  where g~ is the conformal metric
 
  @enddesc
@@*/

#ifndef DCGDT_GUTS
#define DCGDT_GUTS

#ifdef FCODE 

      DCGDT_DCGXXDT = - 2*DCGDT_A*DCGDT_KXX
      DCGDT_DCGXYDT = - 2*DCGDT_A*DCGDT_KXY
      DCGDT_DCGXZDT = - 2*DCGDT_A*DCGDT_KXZ
      DCGDT_DCGYYDT = - 2*DCGDT_A*DCGDT_KYY
      DCGDT_DCGYZDT = - 2*DCGDT_A*DCGDT_KYZ
      DCGDT_DCGZZDT = - 2*DCGDT_A*DCGDT_KZZ

      IF (shift_state .ne. 0) THEN

#include "LIEG_guts.h"

        DCGDT_DCGXXDT = DCGDT_DCGXXDT + LIEG_LGXX
        DCGDT_DCGXYDT = DCGDT_DCGXYDT + LIEG_LGXY
        DCGDT_DCGXZDT = DCGDT_DCGXZDT + LIEG_LGXZ
        DCGDT_DCGYYDT = DCGDT_DCGYYDT + LIEG_LGYY
        DCGDT_DCGYZDT = DCGDT_DCGYZDT + LIEG_LGYZ
        DCGDT_DCGZZDT = DCGDT_DCGZZDT + LIEG_LGZZ

      END IF

#endif

#ifdef CCODE

      DCGDT_DCGXXDT = - 2*DCGDT_A*DCGDT_KXX;
      DCGDT_DCGXYDT = - 2*DCGDT_A*DCGDT_KXY;
      DCGDT_DCGXZDT = - 2*DCGDT_A*DCGDT_KXZ;
      DCGDT_DCGYYDT = - 2*DCGDT_A*DCGDT_KYY;
      DCGDT_DCGYZDT = - 2*DCGDT_A*DCGDT_KYZ;
      DCGDT_DCGZZDT = - 2*DCGDT_A*DCGDT_KZZ;

      if (*shift_state != 0)
      {

#include "LIEG_guts.h"

        DCGDT_DCGXXDT = DCGDT_DCGXXDT + LIEG_LGXX;
        DCGDT_DCGXYDT = DCGDT_DCGXYDT + LIEG_LGXY;
        DCGDT_DCGXZDT = DCGDT_DCGXZDT + LIEG_LGXZ;
        DCGDT_DCGYYDT = DCGDT_DCGYYDT + LIEG_LGYY;
        DCGDT_DCGYZDT = DCGDT_DCGYZDT + LIEG_LGYZ;
        DCGDT_DCGZZDT = DCGDT_DCGZZDT + LIEG_LGZZ;

      }

#endif

#endif

