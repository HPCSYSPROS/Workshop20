/*@@
  @file      Initialise.c
  @date      May 12 2001
  @author    Gabrielle Allen
  @desc
  Initialise grid variables
  @enddesc
@@*/

#include <cctk.h>
#include <cctk_Arguments.h>

void Time_Initialise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  *courant_wave_speed = 0;
  *courant_min_time = 0;
  *courant_dt = 0;
}
