#include <cctk.h>
#include <cctk_Arguments.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

static void register_constrained(const char *gn);

void ML_BSSN_RegisterConstrained(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  register_constrained("ADMBase::metric");
  register_constrained("ADMBase::curv");
  register_constrained("ADMBase::lapse");
  register_constrained("ADMBase::shift");
  register_constrained("ADMBase::dtlapse");
  register_constrained("ADMBase::dtshift");
}

static void register_constrained(const char *gn) {
  assert(gn);

  int gi = CCTK_GroupIndex(gn);
  int ierr = MoLRegisterConstrainedGroup(gi);
  assert(!ierr);
}
