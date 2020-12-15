#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"
#include "ZelmaniQuadWaveExtract.hh"




extern "C" { void ZelmaniQuadWaveExtract_SetIntegrandPointer(CCTK_ARGUMENTS);
}

ZelmaniQuadWaveIntegrandPointer_t ZelmaniQuadWaveIntegrandPointer = NULL;


void ZelmaniQuadWaveExtract_SetIntegrandPointer(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  //set up function pointer

  if (CCTK_EQUALS(integrand_type,"legacy")) {
    ZelmaniQuadWaveIntegrandPointer = ZelmaniQuadWaveLegacy;
  } else if (CCTK_EQUALS(integrand_type,"legacy_radius_criterion")) {
    ZelmaniQuadWaveIntegrandPointer = ZelmaniQuadWaveLegacyRadiusCriterion;
  } else if (CCTK_EQUALS(integrand_type,"legacy_density_criterion")) {
    ZelmaniQuadWaveIntegrandPointer = ZelmaniQuadWaveLegacyDensityCriterion;
  } else if (CCTK_EQUALS(integrand_type,"edens_radius_criterion")) {
    ZelmaniQuadWaveIntegrandPointer = ZelmaniQuadWaveEdensRadiusCriterion;
  } else if (CCTK_EQUALS(integrand_type,"shibata")) {
    ZelmaniQuadWaveIntegrandPointer = ZelmaniQuadWaveShibata;
  } else {
    CCTK_WARN(0,"Integrand Type not defined");
  }




  return;
}

