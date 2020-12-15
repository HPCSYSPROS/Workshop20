#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"

void ADMBase_LapseStatic(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  int npoints;

  if (CCTK_ActiveTimeLevelsGN(cctkGH, "ADMBase::lapse") > 1) {
    npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
    memcpy(alp, alp_p, npoints * sizeof *alp);
  }

  if (CCTK_ActiveTimeLevelsGN(cctkGH, "ADMBase::dtlapse") > 1) {
    npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
    memcpy(dtalp, dtalp_p, npoints * sizeof *dtalp);
  }
}

void ADMBase_ShiftStatic(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  int npoints;

  if (CCTK_ActiveTimeLevelsGN(cctkGH, "ADMBase::shift") > 1) {
    npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
    memcpy(betax, betax_p, npoints * sizeof *betax);
    memcpy(betay, betay_p, npoints * sizeof *betay);
    memcpy(betaz, betaz_p, npoints * sizeof *betaz);
  }

  if (CCTK_ActiveTimeLevelsGN(cctkGH, "ADMBase::dtshift") > 1) {
    npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
    memcpy(dtbetax, dtbetax_p, npoints * sizeof *dtbetax);
    memcpy(dtbetay, dtbetay_p, npoints * sizeof *dtbetay);
    memcpy(dtbetaz, dtbetaz_p, npoints * sizeof *dtbetaz);
  }
}

void ADMBase_Static(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  int npoints;

  if (CCTK_ActiveTimeLevelsGN(cctkGH, "ADMBase::metric") > 1) {
    npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
    memcpy(gxx, gxx_p, npoints * sizeof *gxx);
    memcpy(gxy, gxy_p, npoints * sizeof *gxy);
    memcpy(gxz, gxz_p, npoints * sizeof *gxz);
    memcpy(gyy, gyy_p, npoints * sizeof *gyy);
    memcpy(gyz, gyz_p, npoints * sizeof *gyz);
    memcpy(gzz, gzz_p, npoints * sizeof *gzz);
  }

  if (CCTK_ActiveTimeLevelsGN(cctkGH, "ADMBase::curv") > 1) {
    npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
    memcpy(kxx, kxx_p, npoints * sizeof *kxx);
    memcpy(kxy, kxy_p, npoints * sizeof *kxy);
    memcpy(kxz, kxz_p, npoints * sizeof *kxz);
    memcpy(kyy, kyy_p, npoints * sizeof *kyy);
    memcpy(kyz, kyz_p, npoints * sizeof *kyz);
    memcpy(kzz, kzz_p, npoints * sizeof *kzz);
  }
}
