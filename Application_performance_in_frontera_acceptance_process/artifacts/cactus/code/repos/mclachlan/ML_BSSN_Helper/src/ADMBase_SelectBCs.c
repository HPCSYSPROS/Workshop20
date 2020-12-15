#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <assert.h>

static void select_bcs(const cGH *cctkGH, const char *gn);

void ML_BSSN_ADMBase_SelectBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  select_bcs(cctkGH, "ADMBase::metric");
  select_bcs(cctkGH, "ADMBase::curv");
  select_bcs(cctkGH, "ADMBase::lapse");
  select_bcs(cctkGH, "ADMBase::shift");
  select_bcs(cctkGH, "ADMBase::dtlapse");
  select_bcs(cctkGH, "ADMBase::dtshift");
}

static void select_bcs(const cGH *cctkGH, const char *gn) {
  DECLARE_CCTK_PARAMETERS;

  int ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, boundary_width,
                                       -1, gn, "none");
  assert(!ierr);
}
