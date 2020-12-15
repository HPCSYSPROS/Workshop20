#include <stdlib.h>
#include <stdio.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <assert.h>


void CCCCGlobalModes_CheckStatus(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(compute_every < 1) {
    *do_stuff = 0;
    return;
  }

  *do_stuff = ( (((cctk_iteration) % compute_every) == 0) || cctk_iteration == 0) &&
    (cctk_time >= (start_time*2.03e2)) && (do_shibata || do_saijo || do_CoM);


  return;
}
