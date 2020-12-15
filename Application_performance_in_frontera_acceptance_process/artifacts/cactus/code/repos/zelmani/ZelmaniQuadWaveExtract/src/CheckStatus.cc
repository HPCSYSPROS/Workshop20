#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"
#include "assert.h"

#include "util_Table.h"


extern "C" { void ZelmaniQuadWaveExtract_CheckStatus(CCTK_ARGUMENTS);
}

void ZelmaniQuadWaveExtract_CheckStatus(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  assert(compute_every > 0);

  if(compute_every < 1) {
    *dostuff = 0;
    return;
  }

  *dostuff =
    (((cctk_iteration-1) % compute_every == 0 && cctk_time >= start_time*2.03e2)
     ||
     cctk_iteration-1 == 0);

  //  fprintf(stderr,"\n cctk_iteration: %d  %d\n",cctk_iteration, *dostuff);

  return;
}

