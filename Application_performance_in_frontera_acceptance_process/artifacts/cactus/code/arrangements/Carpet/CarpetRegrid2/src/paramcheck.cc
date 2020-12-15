#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

namespace CarpetRegrid2 {

using namespace Carpet;

extern "C" {
void CarpetRegrid2_ParamCheck(CCTK_ARGUMENTS);
}

void CarpetRegrid2_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  enum sym_t { sym_unknown, sym_90, sym_180, sym_parity };

  int num_params = 0;
  sym_t params = sym_unknown;
  char const *param = "";
  if (symmetry_rotating90) {
    ++num_params;
    params = sym_90;
    param = "symmetry_rotating90";
  }
  if (symmetry_rotating180) {
    ++num_params;
    params = sym_180;
    param = "symmetry_rotating180";
  }
  if (symmetry_parity) {
    ++num_params;
    params = sym_parity;
    param = "symmetry_parity";
  }

  int num_thorns = 0;
  sym_t thorns = sym_unknown;
  char const *thorn = "";
  if (CCTK_IsThornActive("RotatingSymmetry90")) {
    ++num_thorns;
    thorns = sym_90;
    thorn = "RotatingSymmetry90";
  }
  if (CCTK_IsThornActive("RotatingSymmetry90r")) {
    ++num_thorns;
    thorns = sym_90;
    thorn = "RotatingSymmetry90r";
  }
  if (CCTK_IsThornActive("RotatingSymmetry180")) {
    ++num_thorns;
    thorns = sym_180;
    thorn = "RotatingSymmetry180";
  }
  if (CCTK_IsThornActive("ParitySymmetry")) {
    ++num_thorns;
    thorns = sym_parity;
    thorn = "ParitySymmetry";
  }

  if (num_params > 1) {
    CCTK_PARAMWARN("Too many of the symmetry parameters at least two of "
                   "symmetry_rotating90, symmetry_rotating180, and "
                   "parity_symmetry are specified.  (At most one of these can "
                   "be specified.)");
  }

  if (num_thorns > 1) {
    CCTK_PARAMWARN("Too many of the symmetry thorns RotatingSymmetry90, "
                   "RotatingSymmetry90r, RotatingSymmetry180, and "
                   "ParitySymmetry are active.  (At most one of these can be "
                   "active.)");
  }

  if (params != sym_unknown and thorns != sym_unknown and params != thorns) {
    CCTK_VParamWarn(CCTK_THORNSTRING, "The symmetry parameters and the active "
                                      "thorns do not agree.  The symmetry "
                                      "parameter \"%s\" and the active thorn "
                                      "\"%s\" cannot be used together.",
                    param, thorn);
  }

  if ((num_centres >= 1 and num_levels_1 > maxreflevels) or
      (num_centres >= 2 and num_levels_2 > maxreflevels) or
      (num_centres >= 3 and num_levels_3 > maxreflevels) or
      (num_centres >= 4 and num_levels_4 > maxreflevels) or
      (num_centres >= 5 and num_levels_5 > maxreflevels) or
      (num_centres >= 6 and num_levels_6 > maxreflevels) or
      (num_centres >= 7 and num_levels_7 > maxreflevels) or
      (num_centres >= 8 and num_levels_8 > maxreflevels) or
      (num_centres >= 9 and num_levels_9 > maxreflevels) or
      (num_centres >= 10 and num_levels_10 > maxreflevels)) {
    CCTK_PARAMWARN("The number of requested refinement levels is larger than "
                   "the maximum number of levels specified by "
                   "Carpet::max_refinement_levels");
  }
}

} // namespace CarpetRegrid2
