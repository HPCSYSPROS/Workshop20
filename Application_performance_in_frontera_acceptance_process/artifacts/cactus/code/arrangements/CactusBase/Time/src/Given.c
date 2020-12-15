/*@@
  @file      Given.c
  @date      September 4 1999
  @author    Gabrielle Allen
  @desc
  Standard specification of timestep
  @enddesc
@@*/

#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void Time_Given(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  cctkGH->cctk_delta_time = timestep;
}
