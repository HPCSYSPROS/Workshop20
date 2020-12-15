#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdlib.h>
#include <math.h>

void GRHydro_CSPD(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL kappa;
  const CCTK_REAL drho = 1.;
  CCTK_REAL deps = atof("nan");
  CCTK_REAL xtemp = atof("nan");
  const CCTK_REAL xye = atof("nan");
  CCTK_INT keyerr = 0;
  CCTK_INT anyerr = 0;
  
  // HACK: get kappa out of EOS assuming P = kappa * rho^Gamma
  EOS_Omni_press(*GRHydro_eos_handle, 0, 0., 1, &drho, &deps, &xtemp, &xye,
                 &kappa, &keyerr, &anyerr);
  kappa *= pdeplete;

  // compute new rho to satisfy ADM Hamiltonina constraint for reduced kappa
  // this assumes Gamma == 2 and v^i == 0
  #pragma omp parallel
  CCTK_LOOP3_ALL(CSPD, cctkGH, i,j,k) {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    const CCTK_REAL E = rho[idx]*(1.+eps[idx]);
    rho[idx] = 0.5/kappa * (-1.+sqrt(1.+4.*kappa*E));
  } CCTK_ENDLOOP3_ALL(CSPD);

  #pragma omp parallel
  {
    CCTK_INT keyerr[cctk_lsh[0]];
    CCTK_INT anyerr;
    for(int k = 0 ; k < cctk_lsh[2] ; k++) {
      for(int j = 0 ; j < cctk_lsh[1] ; j++) {
        const int idx = CCTK_GFINDEX3D(cctkGH, 0,j,k);
        EOS_Omni_press(*GRHydro_eos_handle, 0, 0., cctk_lsh[0], &rho[idx],
                        &eps[idx], &temperature[idx], &Y_e[idx], &press[idx],
                        keyerr, &anyerr);
        EOS_Omni_EpsFromPress(*GRHydro_eos_handle, 0, 0., cctk_lsh[0], &rho[idx],
                              &eps[idx], &temperature[idx], &Y_e[idx],
                              &press[idx], &eps[idx], keyerr, &anyerr);
      }
    }
  }
}
